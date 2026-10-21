#!/usr/bin/env python3
"""Read-only scientific-payload inspection before public release upload.

Checks all manifest-bound original bytes and decodes gzip, Zstandard and
NumPy ZIP containers in bounded chunks. This is a credential-marker and
format/provenance inspection, not a guarantee about every kind of private
information and not a mathematical checker. Findings contain no matched values.
"""
from __future__ import annotations

import ast
from compression import zstd
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import struct
import zipfile

ROOT = Path(__file__).resolve().parents[1]
CHUNK = 1024**2
PATTERNS = {
    "private_key": re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH |DSA |ENCRYPTED )?PRIVATE KEY-----"),
    "github_token": re.compile(rb"gh[pousr]_[A-Za-z0-9]{36,}|github_pat_[A-Za-z0-9_]{70,}"),
    "aws_access_key": re.compile(rb"AKIA[A-Z0-9]{16}"),
    "credential_field": re.compile(
        rb'"(?:password|passwd|api_key|api_secret|access_token|refresh_token|private_key|client_secret)"\s*:\s*"[^"\s]{8,}"',
        re.IGNORECASE),
}


def scan(source):
    h, count, tail, head, findings = hashlib.sha256(), 0, b"", b"", set()
    while block := source.read(CHUNK):
        count += len(block)
        if count > 64 * 1024**3:
            raise ValueError("Single decoded stream exceeds inspection limit")
        h.update(block)
        if not head:
            head = block[:4096]
        for name, pattern in PATTERNS.items():
            if pattern.search(tail + block):
                findings.add(name)
        tail = block[-512:]
    return {"bytes": count, "sha256": h.hexdigest(),
            "marker_categories": sorted(findings)}, head


def npy_header(stream):
    magic = stream.read(8)
    if magic[:6] != b"\x93NUMPY" or magic[6] not in {1, 2, 3}:
        raise ValueError("Unexpected NumPy member format")
    n = 2 if magic[6] == 1 else 4
    size = struct.unpack("<H" if n == 2 else "<I", stream.read(n))[0]
    if size > 64 * 1024:
        raise ValueError("Oversized NumPy header")
    header = ast.literal_eval(stream.read(size).decode("utf-8"))
    return {k: header[k] for k in ("descr", "fortran_order", "shape")}


def main():
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise SystemExit("Use tools/run_publication_io.py for this read-only inspection")
    manifest_path = ROOT / "repro/public-data-manifest.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    rows = [row for part in manifest["parts"] for row in part["files"]]
    report = {"scope": "Read-only payload bytes, format and known credential markers; not a universal privacy guarantee or mathematical proof.",
              "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
              "files": [], "findings": []}
    for index, row in enumerate(rows, 1):
        name = row["path"]
        path = ROOT / name
        if path.is_symlink() or not path.is_file() or path.stat().st_size != row["bytes"]:
            raise ValueError("Original missing or changed: " + name)
        with path.open("rb") as source:
            checked, head = scan(source)
        if checked["sha256"] != row["sha256"]:
            raise ValueError("Original SHA differs: " + name)
        item = {"path": name, **checked}
        if path.suffix == ".json":
            if not head.lstrip().startswith((b"{", b"[")):
                raise ValueError("Unexpected JSON prefix: " + name)
            item["format"] = "JSON research record"
            item["initial_keys"] = list(dict.fromkeys(k.decode("utf-8", "replace") for k in re.findall(rb'"([^"\r\n]{1,100})"\s*:', head)))[:24]
        elif path.suffix in {".gz", ".zst"}:
            opener = gzip.open if path.suffix == ".gz" else zstd.open
            with opener(path, "rb") as source:
                decoded, decoded_head = scan(source)
            item["decoded"] = decoded
            item["format"] = "compressed research stream"
            item["decoded_prefix"] = decoded_head[:120].decode("utf-8", "replace")
        elif path.suffix == ".npz":
            item["format"] = "NumPy arrays"
            item["members"] = []
            with zipfile.ZipFile(path) as container:
                for info in container.infolist():
                    if not info.filename.endswith(".npy") or info.is_dir():
                        raise ValueError("Unexpected NumPy archive member: " + name)
                    with container.open(info) as source:
                        header = npy_header(source)
                    with container.open(info) as source:
                        decoded, _ = scan(source)
                    item["members"].append({"name": info.filename, "header": header, **decoded})
                    if "O" in str(header["descr"]):
                        report["findings"].append({"path": name, "category": "NumPy object dtype requires manual provenance review"})
        elif path.suffix == ".sqlite":
            if not head.startswith(b"SQLite format 3\x00"):
                raise ValueError("Unexpected SQLite magic: " + name)
            item["format"] = "SQLite research database"
            with sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True) as connection:
                item["schema"] = connection.execute("SELECT type,name,sql FROM sqlite_master ORDER BY type,name LIMIT 100").fetchall()
        elif path.suffix == ".pdf":
            if not head.startswith(b"%PDF-"):
                raise ValueError("Unexpected PDF magic: " + name)
            item["format"] = "reference PDF; separate provenance needed"
            report["findings"].append({"path": name, "category": "reference PDF needs separate provenance/redistribution review"})
        else:
            raise ValueError("Unreviewed format: " + name)
        inspected = [item, item.get("decoded", {})] + item.get("members", [])
        for stream in inspected:
            if stream.get("marker_categories"):
                report["findings"].append({"path": name, "category": "credential markers", "markers": stream["marker_categories"]})
        report["files"].append(item)
        if index % 25 == 0 or index == len(rows):
            print(json.dumps({"inspected_files": index, "total": len(rows), "findings": len(report["findings"])}), flush=True)
    report["checked_original_bytes"] = sum(r["bytes"] for r in report["files"])
    report["status"] = "NO_MARKER_FINDINGS" if not report["findings"] else "REQUIRES_MANUAL_FINDING_REVIEW"
    target = ROOT / ".publication/payload_contents_inspection_20261021.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "files": len(rows), "report": str(target.relative_to(ROOT))}), flush=True)


if __name__ == "__main__":
    main()
