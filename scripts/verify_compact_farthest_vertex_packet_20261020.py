#!/usr/bin/env python3
"""Read-only source-integrity replay; not mathematical or Lean verification."""

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="results/compact_farthest_vertex_dimension24_packet_20261020.json",
    )
    args = parser.parse_args()
    manifest_path = (ROOT / args.manifest).resolve()
    if not manifest_path.is_relative_to(ROOT):
        raise ValueError("Manifest must be inside this repository")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != "compact-farthest-vertex-source-packet-v1":
        raise ValueError("Unsupported manifest schema")
    entries = manifest["files"]
    if not entries:
        raise ValueError("Empty source packet")
    seen = set()
    total_bytes = 0
    roles = {}
    for entry in entries:
        relative = entry["path"]
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT) or relative in seen:
            raise ValueError("Unsafe or duplicate packet path: " + relative)
        seen.add(relative)
        expected_hash = entry["sha256"]
        if len(expected_hash) != 64 or any(
            c not in "0123456789abcdef" for c in expected_hash
        ):
            raise ValueError("Invalid hash: " + relative)
        digest = hashlib.sha256()
        size = 0
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(65536), b""):
                digest.update(block)
                size += len(block)
        if size != entry["bytes"] or digest.hexdigest() != expected_hash:
            raise ValueError("Packet bytes changed: " + relative)
        total_bytes += size
        role = entry["role"]
        roles[role] = roles.get(role, 0) + 1
    print(json.dumps({
        "status": "SOURCE_INTEGRITY_PASS",
        "files": len(entries),
        "bytes": total_bytes,
        "roles": roles,
        "mathematical_or_lean_verification": False,
        "scope": "All listed files match the packet; this does not prove the theorem.",
    }, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError) as error:
        print("SOURCE_INTEGRITY_FAIL: " + str(error), file=sys.stderr)
        sys.exit(1)
