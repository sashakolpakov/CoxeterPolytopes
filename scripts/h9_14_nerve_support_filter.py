"""Build or independently verify the H9/14 nerve-support manifest."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.combinatorics.missing_edges import verify_manifest
from coxeter_polytopes.combinatorics.nerve_support import (
    build_support_manifest,
    verify_support_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data" / "enumerations" / "h9_14_missing_edges.json"
DEFAULT_OUTPUT = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def generate(source_path: Path, output_path: Path) -> dict[str, object]:
    source_manifest = json.loads(source_path.read_text(encoding="utf-8"))
    verify_manifest(source_manifest)
    graph6_values = [
        record["graph6"]
        for record in source_manifest["records"]
        if record["survives_graph_prefilter"]
    ]
    source = {
        "path": str(source_path),
        "sha256": _file_sha256(source_path),
        "schema_version": source_manifest["schema_version"],
        "graph_prefilter_survivors": len(graph6_values),
    }
    manifest = build_support_manifest(graph6_values, source)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def verify(path: Path) -> dict[str, object]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    return verify_support_manifest(manifest)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="verify the archived output without rebuilding it",
    )
    args = parser.parse_args(argv)

    if args.verify_only:
        summary = verify(args.output)
        mode = "verified"
    else:
        manifest = generate(args.source, args.output)
        summary = verify_support_manifest(manifest)
        mode = "generated and verified"
    print(
        json.dumps(
            {"mode": mode, "output": str(args.output), "counts": summary},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
