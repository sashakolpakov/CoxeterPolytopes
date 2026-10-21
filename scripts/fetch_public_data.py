#!/usr/bin/env python3
"""Plan or restore checksum-verified bulk inputs from the public data release.

Uses only Python's standard library. Downloads and extraction are sequential;
an existing different file is never overwritten. Data integrity is not a
mathematical verification of the enclosed certificates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import urllib.request

MIB = 1024**2
ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(MIB), b""):
            h.update(block)
    return h.hexdigest()


def safe_target(root: Path, name: str) -> Path:
    if not isinstance(name, str) or "\\" in name or "\x00" in name:
        raise ValueError("Invalid manifest path")
    rel = PurePosixPath(name)
    if rel.as_posix() != name or rel.is_absolute() or not rel.parts or any(p in {"..", "."} for p in rel.parts):
        raise ValueError(f"Unsafe manifest path: {name}")
    target = root.joinpath(*rel.parts)
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes destination: {name}")
    for parent in [target, *target.parents]:
        if parent == root.parent:
            break
        if parent.is_symlink():
            raise ValueError(f"Refusing a symlink destination: {name}")
    return target


def existing_matches(path: Path, row: dict) -> bool:
    if not path.exists():
        return False
    if not path.is_file() or path.stat().st_size != row["bytes"] or digest(path) != row["sha256"]:
        raise ValueError(f"Existing file differs; preserve/resolve it before restoring: {path}")
    return True


def download(url: str, target: Path, part: dict) -> None:
    partial = target.with_suffix(target.suffix + ".download")
    if any(p.exists() or p.is_symlink() for p in (target, partial)):
        raise ValueError(f"Unfinished download exists; inspect it first: {partial}")
    request = urllib.request.Request(url, headers={"User-Agent": "CoxeterPolytopes-data/1"})
    h, count = hashlib.sha256(), 0
    owned = False
    try:
        with partial.open("xb") as out:
            owned = True
            with urllib.request.urlopen(request, timeout=60) as response:
                while block := response.read(MIB):
                    count += len(block)
                    if count > part["bytes"]:
                        raise ValueError("Download exceeds its manifest size")
                    out.write(block)
                    h.update(block)
                    if count % (64 * MIB) < MIB:
                        check_disk(target.parent)
        if count != part["bytes"] or h.hexdigest() != part["sha256"]:
            raise ValueError(f"Archive checksum/size mismatch: {part['name']}")
        os.link(partial, target)
        partial.unlink()
    except BaseException:
        if owned:
            partial.unlink(missing_ok=True)
        raise


def restore_part(archive: Path, part: dict, destination: Path, wanted: set[str]) -> int:
    rows = {row["path"]: row for row in part["files"]}
    if len(rows) != len(part["files"]):
        raise ValueError("Duplicate manifest paths")
    seen, restored = set(), 0
    with tarfile.open(archive, "r|gz") as tar:
        for member in tar:
            if not member.isfile() or member.name not in rows or member.name in seen:
                raise ValueError("Unexpected, linked, or duplicate archive member")
            seen.add(member.name)
            row = rows[member.name]
            if member.size != row["bytes"]:
                raise ValueError("Archive member size disagrees with manifest")
            if member.name not in wanted:
                continue
            target = safe_target(destination, member.name)
            if existing_matches(target, row):
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            temp = target.with_name(target.name + ".public-data.partial")
            if temp.exists() or temp.is_symlink():
                raise ValueError(f"Unfinished restore exists; inspect it first: {temp}")
            h, count = hashlib.sha256(), 0
            owned = False
            try:
                with temp.open("xb") as out:
                    owned = True
                    source = tar.extractfile(member)
                    while block := source.read(MIB):
                        out.write(block)
                        h.update(block)
                        count += len(block)
                        if count % (64 * MIB) < MIB:
                            check_disk(destination)
                if count != row["bytes"] or h.hexdigest() != row["sha256"]:
                    raise ValueError(f"Restored bytes disagree: {member.name}")
                # No existing destination can be replaced by a concurrent writer.
                os.link(temp, target)
                temp.unlink()
                restored += 1
            except BaseException:
                if owned:
                    temp.unlink(missing_ok=True)
                raise
    if seen != set(rows):
        raise ValueError("Archive omits manifest members")
    return restored


def check_disk(root: Path) -> None:
    if shutil.disk_usage(root).free < 20 * 1024**3:
        raise ValueError("Restoration stopped: less than 20 GiB disk free")


def validate_manifest(data: dict, destination: Path) -> list[dict]:
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported public-data manifest version")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", data.get("repository", "")) or not re.fullmatch(r"[A-Za-z0-9_.-]+", data.get("release_tag", "")):
        raise ValueError("Invalid release location")
    names, paths, rows = set(), set(), []
    def fields(item):
        if type(item.get("bytes")) is not int or item["bytes"] < 0 or not re.fullmatch(r"[0-9a-f]{64}", item.get("sha256", "")):
            raise ValueError("Invalid size or SHA-256 in manifest")
    for part in data["parts"]:
        name = part["name"]
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_-]+(?:\.[A-Za-z0-9_-]+)*\.tar\.gz", name) or name in names:
            raise ValueError("Unsafe or duplicate archive filename")
        names.add(name)
        fields(part)
        for row in part["files"]:
            fields(row)
            safe_target(destination, row["path"])
            if row["path"] in paths:
                raise ValueError("Manifest repeats a destination path")
            paths.add(row["path"])
            rows.append(row)
    if type(data.get("original_bytes")) is not int or data["original_bytes"] != sum(row["bytes"] for row in rows):
        raise ValueError("Manifest total disagrees with its file sizes")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true", help="Inspect only; no network or writes")
    mode.add_argument("--restore", action="store_true", help="Download and restore missing inputs")
    parser.add_argument("--path", action="append", default=[], help="Restore only this exact repository-relative path; repeatable")
    parser.add_argument("--manifest", type=Path, default=ROOT / "repro/public-data-manifest.json")
    parser.add_argument("--destination", type=Path, default=ROOT)
    parser.add_argument("--from-directory", type=Path, help="Use already downloaded archives instead of the network")
    parser.add_argument("--keep-cache", action="store_true")
    args = parser.parse_args()
    data = json.loads(args.manifest.read_text())
    destination = args.destination.resolve()
    all_rows = validate_manifest(data, destination)
    all_paths = {row["path"] for row in all_rows}
    if len(all_paths) != len(all_rows):
        raise ValueError("Manifest repeats a destination path")
    wanted = set(args.path) if args.path else all_paths
    if not wanted <= all_paths:
        raise ValueError(f"Unknown public input paths: {sorted(wanted-all_paths)}")
    for path in wanted:
        safe_target(destination, path)
    parts = [part for part in data["parts"] if any(row["path"] in wanted for row in part["files"])]
    missing = [row for row in all_rows if row["path"] in wanted and not existing_matches(safe_target(destination, row["path"]), row)]
    summary = {"files_selected": len(wanted), "files_missing": len(missing),
               "uncompressed_bytes_missing": sum(row["bytes"] for row in missing),
               "archive_parts_selected": len(parts),
               "archive_bytes_selected": sum(part["bytes"] for part in parts),
               "release": f"https://github.com/{data['repository']}/releases/tag/{data['release_tag']}",
               "mathematical_verification": False}
    print(json.dumps(summary, indent=2), flush=True)
    if args.plan or not missing:
        return
    if not destination.exists():
        destination.mkdir(parents=True)
    needed = sum(row["bytes"] for row in missing) + max(part["bytes"] for part in parts) + 2 * 1024**3
    if shutil.disk_usage(destination).free < max(20 * 1024**3, needed):
        raise ValueError("Not enough disk space for the sequential restoration; select fewer --path inputs")
    cache = safe_target(destination, ".publication/downloads")
    cache.mkdir(parents=True, exist_ok=True)
    missing_names = {row["path"] for row in missing}
    for part in parts:
        if not any(row["path"] in missing_names for row in part["files"]):
            continue
        if Path(part["name"]).name != part["name"]:
            raise ValueError("Unsafe archive filename")
        if args.from_directory:
            archive = args.from_directory / part["name"]
        else:
            archive = cache / part["name"]
            if not archive.exists():
                url = f"https://github.com/{data['repository']}/releases/download/{data['release_tag']}/{part['name']}"
                download(url, archive, part)
        if archive.stat().st_size != part["bytes"] or digest(archive) != part["sha256"]:
            raise ValueError(f"Archive checksum/size mismatch: {archive}")
        count = restore_part(archive, part, destination, wanted)
        print(json.dumps({"archive": part["name"], "restored_files": count}), flush=True)
        if not args.keep_cache and not args.from_directory:
            archive.unlink()


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, tarfile.TarError) as error:
        raise SystemExit(str(error))
