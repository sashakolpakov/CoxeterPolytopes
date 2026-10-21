#!/usr/bin/env python3
"""Reclaim APFS allocation without changing certificate or transcript bytes.

Only static regular JSON / JSONL files are eligible. Each file is handled alone,
with streaming hashes and a same-directory temporary copy. This never parses a
certificate, loads a transcript into memory, deletes history, or touches SQLite.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import uuid

GIB = 1024**3
MIB = 1024**2
COMPRESSED = 0x20
IMMUTABLE_OR_APPEND = 0x00000006 | 0x00060000
WORKSPACE = Path(__file__).resolve().parent.parent
SESSIONS = Path.home() / ".codex" / "sessions"


def identity(st: os.stat_result) -> tuple[int, ...]:
    return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb", buffering=0) as source:
        while block := source.read(MIB):
            h.update(block)
    return h.hexdigest()


def open_now(path: Path) -> bool:
    result = subprocess.run(
        ["/usr/sbin/lsof", "-n", "-P", "-t", "--", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
        timeout=30,
    )
    # lsof returns 1 when no process has the file open. Fail closed otherwise.
    return result.returncode != 1 or bool(result.stdout.strip())


def candidate_paths(root: Path, suffix: str, age: float, limit: int) -> list[Path]:
    # Metadata only: a few thousand paths, never file contents in memory.
    now = time.time()
    candidates: list[tuple[int, str]] = []
    for base, _dirs, files in os.walk(root):
        for name in files:
            if not name.endswith(suffix):
                continue
            path = Path(base) / name
            st = path.lstat()
            if (
                not path.is_symlink()
                and path.is_file()
                and st.st_nlink == 1
                and 2 * MIB <= st.st_size <= limit
                and now - st.st_mtime >= age
                and not st.st_flags & (COMPRESSED | IMMUTABLE_OR_APPEND)
            ):
                candidates.append((st.st_blocks * 512, str(path)))
    # Largest bounded files first: substantial recovery in short batches.
    candidates.sort(reverse=True)
    return [Path(path) for _, path in candidates]


def compress(path: Path, min_age: float) -> dict[str, object]:
    original = path.lstat()
    if (
        path.is_symlink()
        or original.st_nlink != 1
        or original.st_flags & (COMPRESSED | IMMUTABLE_OR_APPEND)
        or time.time() - original.st_mtime < min_age
        or open_now(path)
    ):
        return {"status": "skipped_changed_or_open"}
    original_attributes = set(os.listxattr(path))
    if "com.apple.ResourceFork" in original_attributes:
        # Preserve non-compression resource-fork contents without interpretation.
        return {"status": "skipped_existing_resource_fork"}
    temp = path.with_name("." + path.name + ".compress-" + uuid.uuid4().hex)
    try:
        subprocess.run(
            ["/usr/bin/ditto", "--hfsCompression", "--nocache", "--noclone", str(path), str(temp)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=True,
            timeout=300,
        )
        # On this APFS setup ditto creates decmpfs / resource-fork attributes but
        # leaves UF_COMPRESSED unset. Setting it makes normal reads decompress.
        os.chflags(temp, original.st_flags | COMPRESSED)
        source_hash = digest(path)
        copy_hash = digest(temp)
        after_hash = path.lstat()
        copied = temp.stat()
        if identity(original) != identity(after_hash):
            raise RuntimeError("source changed while copying; original preserved")
        if source_hash != copy_hash or original.st_size != copied.st_size:
            raise RuntimeError("decompressed byte comparison failed; original preserved")
        if (copied.st_mode, copied.st_uid, copied.st_gid) != (
            original.st_mode,
            original.st_uid,
            original.st_gid,
        ):
            raise RuntimeError("file permissions or ownership changed; original preserved")
        for attribute in original_attributes:
            if os.getxattr(path, attribute) != os.getxattr(temp, attribute):
                raise RuntimeError("original extended attribute changed; original preserved")
        before_bytes = original.st_blocks * 512
        after_bytes = copied.st_blocks * 512
        if after_bytes >= before_bytes:
            return {"status": "skipped_no_saving"}
        # Hash reads can change access time. Restore source times on replacement.
        os.utime(temp, ns=(original.st_atime_ns, original.st_mtime_ns))
        if identity(original) != identity(path.lstat()) or open_now(path):
            raise RuntimeError("source changed or was opened before replacement")
        os.replace(temp, path)
        final = path.stat()
        if final.st_size != original.st_size or final.st_mtime_ns != original.st_mtime_ns:
            raise RuntimeError("unexpected replacement metadata")
        return {
            "status": "compressed_verified",
            "sha256": source_hash,
            "logical_bytes": original.st_size,
            "allocated_before": before_bytes,
            "allocated_after": final.st_blocks * 512,
            "reclaimed_bytes": before_bytes - final.st_blocks * 512,
        }
    finally:
        if temp.exists():
            temp.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=("workspace", "sessions"), required=True)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--max-files", type=int, default=32)
    parser.add_argument("--max-input-gib", type=float, default=8)
    parser.add_argument("--max-file-mib", type=int, default=512)
    parser.add_argument("--target-free-gib", type=float, default=60)
    parser.add_argument("--reserve-gib", type=float, default=4)
    parser.add_argument("--log", type=Path)
    args = parser.parse_args()
    if sys.platform != "darwin" or not hasattr(os, "chflags"):
        parser.error("this exact-byte compression procedure requires macOS/APFS")
    if args.max_files < 1 or not 2 <= args.max_file_mib <= 512:
        parser.error("positive batch size and 2..512 MiB per-file limit required")
    def interrupted(signum: int, _frame: object) -> None:
        raise KeyboardInterrupt(f"cleanup interrupted by signal {signum}")
    for sig in (signal.SIGTERM, signal.SIGHUP):
        signal.signal(sig, interrupted)
    root = WORKSPACE if args.scope == "workspace" else SESSIONS
    suffix = ".json" if args.scope == "workspace" else ".jsonl"
    min_age = 600 if args.scope == "workspace" else 7 * 86400
    candidates = candidate_paths(root, suffix, min_age, args.max_file_mib * MIB)
    initial_free = shutil.disk_usage(root).free
    report: dict[str, object] = {
        "scope": args.scope,
        "apply": args.apply,
        "candidate_files": len(candidates),
        "initial_free_bytes": initial_free,
        "files": [],
    }
    records: list[dict[str, object]] = report["files"]  # type: ignore[assignment]
    input_bytes = 0
    started = time.monotonic()
    for path in candidates:
        free = shutil.disk_usage(root).free
        if free < max(4, args.reserve_gib) * GIB:
            report["stop_reason"] = "reserve_floor"
            break
        if free >= args.target_free_gib * GIB:
            report["stop_reason"] = "target_reached"
            break
        if len(records) >= args.max_files or input_bytes >= args.max_input_gib * GIB:
            report["stop_reason"] = "batch_limit"
            break
        st = path.stat()
        if input_bytes + st.st_size > args.max_input_gib * GIB:
            continue
        input_bytes += st.st_size
        if args.apply:
            try:
                record = compress(path, min_age)
            except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
                record = {"status": "error_original_preserved", "error": str(exc)}
        else:
            record = {"status": "would_compress", "logical_bytes": st.st_size}
        record["path"] = str(path.relative_to(root))
        records.append(record)
        if args.apply:
            # Keep output compact; no content or transcript excerpts are emitted.
            print(json.dumps({"completed_files": len(records), "last_status": record["status"],
                              "free_gib": round(shutil.disk_usage(root).free / GIB, 2)}), flush=True)
    report["final_free_bytes"] = shutil.disk_usage(root).free
    report["free_space_delta_bytes"] = report["final_free_bytes"] - initial_free  # type: ignore[operator]
    report["reclaimed_allocation_bytes"] = sum(int(r.get("reclaimed_bytes", 0)) for r in records)
    report["input_bytes"] = input_bytes
    report["elapsed_seconds"] = round(time.monotonic() - started, 2)
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        args.log.write_text(json.dumps(report, indent=2) + "\n")
    compact = {k: v for k, v in report.items() if k != "files"}
    compact["status_counts"] = {status: sum(r["status"] == status for r in records)
                               for status in sorted({str(r["status"]) for r in records})}
    print(json.dumps(compact, sort_keys=True))


if __name__ == "__main__":
    main()
