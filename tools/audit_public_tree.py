#!/usr/bin/env python3
"""Reject build debris, oversized blobs and secret markers in the Git index.

Reads staged Git blobs, not potentially different working files. Keeps the
six intentional paper PDFs, source-generated mathematical certificates and
verification records. This checks publication hygiene, not mathematics.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import PurePosixPath
import re
import subprocess

PAPERS = {
    "manuscripts/build/finite_volume57.pdf",
    "manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf",
    "manuscripts/build/compact_arithmetic_dimension_bound.pdf",
    "manuscripts/build/compact_farthest_vertex_bounds.pdf",
    "manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf",
    "manuscripts/build/compact_ra4_volume_constraints.pdf",
}
SUFFIXES = {".pyc", ".pyo", ".olean", ".ilean", ".o", ".bc", ".ir",
            ".trace", ".aux", ".toc", ".lof", ".lot", ".fls", ".fdb_latexmk",
            ".synctex", ".nav", ".snm", ".vrb", ".xdv", ".sqlite-shm", ".sqlite-wal",
            ".so", ".dylib", ".dll", ".a", ".exe", ".partial", ".download", ".tmp"}
CACHE_DIRS = {".lake", ".proof-batches", ".publication", ".venv", "venv",
              "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache",
              "_build", "node_modules"}
SECRET = re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH |DSA |ENCRYPTED )?PRIVATE KEY-----"
                    rb"|-----BEGIN PGP PRIVATE " rb"KEY BLOCK-----"
                    rb"|gh[pousr]_[A-Za-z0-9]{36,}"
                    rb"|github_pat_[A-Za-z0-9_]{70,}"
                    rb"|AKIA[A-Z0-9]{16}")


def debris(name: str) -> bool:
    p = PurePosixPath(name)
    return (bool(set(p.parts) & CACHE_DIRS) or p.suffix in SUFFIXES
            or name.endswith((".olean.hash", ".ilean.hash", ".setup.json", ".synctex.gz"))
            or name.startswith("build/")
            or (name.startswith("formal/") and p.suffix in {".c", ".hash"})
            or (name.startswith("manuscripts/build/") and p.suffix in {".log", ".out"}))


def main() -> None:
    raw = subprocess.check_output(["git", "ls-files", "--stage", "-z"])
    entries = []
    for record in raw.split(b"\0"):
        if not record:
            continue
        meta, path = record.split(b"\t", 1)
        mode, oid, stage = meta.split()
        if stage != b"0":
            raise SystemExit("Unmerged index entries; finish the merge before publication")
        entries.append((path.decode(), mode, oid))
    if not entries:
        raise SystemExit("Empty Git index; nothing has been audited")
    issues, total, largest, manifest = [], 0, (0, ""), None
    blob_list = hashlib.sha256()
    # Bound Git's mapped pack windows and delta cache after the first commit
    # has packed the repository. Default caches can exceed the I/O watchdog.
    child = subprocess.Popen(["git", "-c", "core.packedGitLimit=16m",
                              "-c", "core.packedGitWindowSize=1m",
                              "-c", "core.deltaBaseCacheLimit=8m",
                              "cat-file", "--batch"], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)
    try:
        for name, mode, oid in entries:
            if mode not in {b"100644", b"100755"}:
                issues.append([name, "symlink/submodule or unsupported mode"])
            if debris(name):
                issues.append([name, "compiler output, cache or temporary build file"])
            child.stdin.write(oid + b"\n")
            child.stdin.flush()
            header = child.stdout.readline().split()
            if len(header) != 3 or header[1] != b"blob":
                raise RuntimeError("Unexpected Git blob response")
            size = remaining = int(header[2])
            total += size
            largest = max(largest, (size, name))
            if size >= 100 * 1024**2:
                issues.append([name, "blob reaches GitHub's 100 MiB limit"])
            blob_list.update(name.encode() + b"\0" + oid + b"\n")
            tail, payload = b"", bytearray()
            while remaining:
                block = child.stdout.read(min(1024**2, remaining))
                if not block:
                    raise RuntimeError("Truncated Git blob response")
                remaining -= len(block)
                if SECRET.search(tail + block):
                    issues.append([name, "possible private key or credential; inspect locally"])
                tail = block[-256:]
                if name == "repro/public-data-manifest.json":
                    if size > 4 * 1024**2:
                        raise RuntimeError("Unexpectedly large publication manifest")
                    payload.extend(block)
            if child.stdout.read(1) != b"\n":
                raise RuntimeError("Invalid Git blob terminator")
            if payload:
                manifest = json.loads(payload)
    finally:
        child.stdin.close()
        child.stdout.close()
        child.wait()
    if child.returncode:
        raise RuntimeError("Git blob audit failed")
    names = {name for name, _, _ in entries}
    for name in sorted(PAPERS - names):
        issues.append([name, "required current paper PDF absent from index"])
    if manifest is None:
        issues.append(["repro/public-data-manifest.json", "bulk-data manifest absent from index"])
    else:
        if manifest.get("schema_version") != 1 or not manifest.get("parts"):
            raise RuntimeError("Missing or unsupported bulk-data manifest")
        archives, paths, original_bytes = set(), set(), 0
        def valid_fields(row):
            return (type(row.get("bytes")) is int and row["bytes"] >= 0
                    and isinstance(row.get("sha256"), str)
                    and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]))
        for part in manifest["parts"]:
            archive = part["name"]
            if not isinstance(archive, str) or not re.fullmatch(r"[A-Za-z0-9_-]+(?:\.[A-Za-z0-9_-]+)*\.tar\.gz", archive) or archive in archives or not valid_fields(part) or not part.get("files"):
                raise RuntimeError("Invalid/duplicate bulk archive record")
            archives.add(archive)
            for row in part["files"]:
                path = row["path"]
                if not isinstance(path, str) or "\\" in path or "\0" in path:
                    raise RuntimeError("Invalid bulk-data path")
                rel = PurePosixPath(path)
                if rel.is_absolute() or rel.as_posix() != path or not rel.parts or ".." in rel.parts or rel.parts[0] not in {"data", "results"} or path in paths or not valid_fields(row):
                    raise RuntimeError("Invalid/duplicate bulk-data file record")
                paths.add(path)
                original_bytes += row["bytes"]
                if path in names:
                    issues.append([path, "release-only bulk data accidentally staged"])
        if type(manifest.get("original_bytes")) is not int or manifest["original_bytes"] != original_bytes:
            raise RuntimeError("Bulk-data total disagrees with file records")
    print(json.dumps({"publication_hygiene_only": True, "index_files": len(entries),
                      "total_blob_bytes": total, "largest_blob_bytes": largest[0],
                      "largest_blob_path": largest[1],
                      "path_and_blob_list_sha256": blob_list.hexdigest(),
                      "issues": issues}, indent=2))
    if issues:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
