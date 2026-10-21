"""Build or verify the exact H9/14 local-polytopal-link manifest."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.combinatorics.local_link import (
    build_local_link_manifest,
    verify_local_link_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
DEFAULT_OUTPUT = ROOT / "data" / "enumerations" / "h9_14_local_links.json"


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def generate(source_path: Path, output_path: Path) -> dict[str, object]:
    source_manifest = json.loads(source_path.read_text(encoding="utf-8"))
    records = source_manifest.get("records")
    if not isinstance(records, list):
        raise ValueError("support manifest records must be a list")
    graph6_values = [
        record["graph6"]
        for record in records
        if record.get("survives_second_layer") is True
    ]
    source = {
        "path": str(source_path.resolve().relative_to(ROOT)),
        "sha256": _file_sha256(source_path),
        "schema_version": source_manifest.get("schema_version"),
        "second_layer_survivors": len(graph6_values),
    }
    manifest = build_local_link_manifest(graph6_values, source)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def verify(output_path: Path) -> dict[str, object]:
    manifest = json.loads(output_path.read_text(encoding="utf-8"))
    return verify_local_link_manifest(manifest)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        summary = verify(args.output)
        mode = "verified"
    else:
        manifest = generate(args.source, args.output)
        summary = verify_local_link_manifest(manifest)
        mode = "generated and verified"
    print(json.dumps({"mode": mode, "summary": summary}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
