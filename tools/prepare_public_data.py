#!/usr/bin/env python3
"""Package large research files as streamed, verified release assets.

This is publication I/O, not a mathematical checker. Originals are never
modified. One archive is compressed and verified at a time, using 1 MiB
buffers. The small manifest is the public fetcher's only extraction allowlist.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

MIB = 1024**2
ROOT = Path(__file__).resolve().parents[1]
PUBLIC_DIRS = {
    "data", "docs", "formal", "manuscripts", "references", "repro",
    "results", "schemas", "scripts", "src", "tests", "tools",
}
EXCLUDED_DIRS = {
    ".git", ".lake", ".venv", "__pycache__", ".proof-batches",
    "node_modules", "_build", ".pytest_cache", ".mypy_cache",
}
EXCLUDED_SUFFIXES = {".pyc", ".olean", ".ilean", ".o", ".aux", ".log", ".toc"}


def public_file(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    if not rel.parts or rel.parts[0] not in PUBLIC_DIRS:
        return False
    if any(p in EXCLUDED_DIRS for p in rel.parts) or path.is_symlink():
        return False
    if path.suffix in EXCLUDED_SUFFIXES:
        return False
    if rel.parts[:2] == ("docs", "history"):
        return False
    # Third-party downloads remain local; their README files give provenance.
    if rel.parts[:2] == ("data", "external") and path.name != "README.md":
        return False
    if "AGENTS" in path.name or path.name.startswith(("GOAL", "RESEARCH_PLAN")):
        return False
    return path.is_file()


def discover() -> list[dict]:
    rows = []
    for name in sorted(PUBLIC_DIRS):
        for base, dirs, files in os.walk(ROOT / name, followlinks=False):
            dirs[:] = sorted(d for d in dirs if d not in EXCLUDED_DIRS)
            for filename in sorted(files):
                p = Path(base) / filename
                if (name in {"data", "results"} and public_file(p)
                        and p.stat().st_size >= MIB):
                    rows.append({"path": p.relative_to(ROOT).as_posix(),
                                 "bytes": p.stat().st_size})
    return sorted(rows, key=lambda r: r["path"])


def groups(rows: list[dict]) -> list[list[dict]]:
    out, current, size = [], [], 0
    for row in rows:
        if current and size + row["bytes"] > 512 * MIB:
            out.append(current)
            current, size = [], 0
        current.append(row)
        size += row["bytes"]
    if current:
        out.append(current)
    return out


def check_resources() -> None:
    if shutil.disk_usage(ROOT).free < 20 * 1024**3:
        raise RuntimeError("Publication packaging stopped: less than 20 GiB disk free")
    # Packaging uses small fixed buffers; proof/Lean admission gates are separate.
    if shutil.which("vm_stat"):
        vm = subprocess.run(["vm_stat"], check=True, text=True,
                            capture_output=True).stdout
        import re
        page = int(re.search(r"page size of (\d+) bytes", vm).group(1))
        count = sum(int(re.search(re.escape(k) + r":\s*(\d+)", vm).group(1))
                    for k in ("Pages free", "Pages inactive", "Pages speculative"))
        if page * count < 3 * 1024**3:
            raise RuntimeError("Publication packaging stopped: less than 3 GiB RAM available")


def digest_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for index, block in enumerate(iter(lambda: f.read(MIB), b"")):
            h.update(block)
            if index % 64 == 0:
                check_resources()
    return h.hexdigest()


class HashedReader:
    def __init__(self, source):
        self.source = source
        self.digest = hashlib.sha256()
        self.since_check = 0

    def read(self, size=-1):
        if size < 0 or size > MIB:
            raise ValueError("Archive reader requires bounded reads")
        block = self.source.read(size)
        self.digest.update(block)
        self.since_check += len(block)
        if self.since_check >= 64 * MIB:
            check_resources()
            self.since_check = 0
        return block


def verify_archive(path: Path, records: list[dict]) -> None:
    expected = {row["path"]: row for row in records}
    if len(expected) != len(records):
        raise RuntimeError("Duplicate archive inventory")
    seen = set()
    with tarfile.open(path, "r|gz") as tar:
        for member in tar:
            check_resources()
            if not member.isfile() or member.name not in expected or member.name in seen:
                raise RuntimeError("Unexpected or duplicate archive member")
            row = expected[member.name]
            if member.size != row["bytes"]:
                raise RuntimeError("Archive member size mismatch")
            h = hashlib.sha256()
            with tar.extractfile(member) as source:
                for index, block in enumerate(iter(lambda: source.read(MIB), b"")):
                    h.update(block)
                    if index % 64 == 0:
                        check_resources()
            if h.hexdigest() != row["sha256"]:
                raise RuntimeError("Archive round-trip mismatch")
            seen.add(member.name)
    if seen != set(expected):
        raise RuntimeError("Archive omitted an expected member")


def pack(group: list[dict], target: Path) -> dict:
    check_resources()
    records = []
    partial = target.with_suffix(target.suffix + ".partial")
    if any(p.exists() or p.is_symlink() for p in (target, partial)):
        raise RuntimeError(f"Refusing to overwrite archive: {target}")
    start = time.monotonic()
    with partial.open("xb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out,
                           compresslevel=1, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w|", format=tarfile.PAX_FORMAT,
                              copybufsize=MIB) as tar:
                for row in group:
                    check_resources()
                    p = ROOT / row["path"]
                    before = p.stat()
                    if before.st_size != row["bytes"]:
                        raise RuntimeError(f"Source changed: {p}")
                    info = tarfile.TarInfo(row["path"])
                    info.size, info.mode, info.mtime = row["bytes"], 0o644, 0
                    with p.open("rb") as source:
                        reader = HashedReader(source)
                        tar.addfile(info, reader)
                    after = p.stat()
                    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                        raise RuntimeError(f"Source changed while packaging: {p}")
                    records.append({**row, "sha256": reader.digest.hexdigest()})
    # Full streaming decode verifies every original byte before publication.
    verify_archive(partial, records)
    os.link(partial, target)  # atomic, and never replaces a competing writer
    partial.unlink()
    return {"name": target.name, "bytes": target.stat().st_size,
            "sha256": digest_file(target), "files": records,
            "roundtrip_verified": True,
            "seconds": round(time.monotonic() - start, 3)}


def validate_completed(state: dict, batches: list[list[dict]], work: Path) -> dict:
    parts = state["parts"]
    completed = {part["name"]: part for part in parts}
    expected = {f"research-data-{i:03d}.tar.gz": group
                for i, group in enumerate(batches, 1)}
    if len(completed) != len(parts) or not completed.keys() <= expected.keys():
        raise RuntimeError("Invalid or stale completed-part inventory")
    for name, part in completed.items():
        if [(r["path"], r["bytes"]) for r in part["files"]] != [
                (r["path"], r["bytes"]) for r in expected[name]]:
            raise RuntimeError("Bulk input inventory changed; use a new release directory")
        for row in part["files"]:
            source = ROOT / row["path"]
            if source.is_symlink() or not source.is_file() or source.stat().st_size != row["bytes"] or digest_file(source) != row["sha256"]:
                raise RuntimeError(f"Completed source changed: {row['path']}")
        archive = work / name
        if archive.is_symlink() or not archive.is_file() or archive.stat().st_size != part["bytes"] or digest_file(archive) != part["sha256"]:
            raise RuntimeError(f"Completed archive absent or changed: {name}")
        verify_archive(archive, part["files"])
    return completed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true")
    parser.add_argument("--max-parts", type=int)
    args = parser.parse_args()
    if args.max_parts is not None and args.max_parts <= 0:
        parser.error("--max-parts must be positive")
    rows = discover()
    batches = groups(rows)
    if args.plan:
        print(json.dumps({"large_files": len(rows),
                          "original_bytes": sum(r["bytes"] for r in rows),
                          "archive_parts": len(batches),
                          "max_raw_part_bytes": max((sum(r["bytes"] for r in g) for g in batches), default=0)}, indent=2))
        return
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise RuntimeError("Use tools/run_publication_io.py for monitored publication I/O")
    work = ROOT / ".publication" / "data"
    work.mkdir(parents=True, exist_ok=True)
    state_path = work / "state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"parts": []}
    completed = validate_completed(state, batches, work)
    new = 0
    for index, group in enumerate(batches, 1):
        name = f"research-data-{index:03d}.tar.gz"
        if name in completed:
            continue
        result = pack(group, work / name)
        state["parts"].append(result)
        temp = state_path.with_suffix(".tmp")
        temp.write_text(json.dumps(state, indent=2) + "\n")
        temp.replace(state_path)
        print(json.dumps({"part": name, "files": len(result["files"]),
                          "bytes": result["bytes"], "seconds": result["seconds"]}), flush=True)
        new += 1
        if args.max_parts is not None and new >= args.max_parts:
            break
    if len(state["parts"]) == len(batches):
        if discover() != rows:
            raise RuntimeError("Bulk input inventory changed during packaging")
        validate_completed(state, batches, work)
        manifest = {"schema_version": 1,
                    "repository": "sashakolpakov/CoxeterPolytopes",
                    "release_tag": "research-data-2026-10-21",
                    "scope": "Large scientific inputs and certificates; original files unchanged. This is data integrity, not mathematical verification.",
                    "original_bytes": sum(r["bytes"] for r in rows),
                    "parts": sorted(state["parts"], key=lambda p: p["name"])}
        dest = ROOT / "repro/public-data-manifest.json"
        if dest.exists():
            if json.loads(dest.read_text()) != manifest:
                raise RuntimeError("Refusing to overwrite a different public manifest")
        else:
            with dest.open("x") as output:
                output.write(json.dumps(manifest, indent=2) + "\n")
        print(json.dumps({"manifest": str(dest.relative_to(ROOT)),
                          "sha256": digest_file(dest), "complete": True}), flush=True)


if __name__ == "__main__":
    main()
