"""Build or verify forced-A4-empty minimal-transversal certificates."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.combinatorics.forced_a0_transversals import (
    build_forced_a0_manifest,
    verify_forced_a0_manifest,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
DEFAULT_OUTPUT = (
    ROOT / "data" / "enumerations" / "h9_14_forced_a0_transversals.json"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_graphs(source_path: Path) -> tuple[dict[str, object], tuple[str, ...]]:
    source_manifest = json.loads(source_path.read_text(encoding="utf-8"))
    verify_support_manifest(source_manifest)
    records = source_manifest.get("records")
    if not isinstance(records, list):
        raise ValueError("support manifest records must be a list")
    graph6_values = tuple(
        record["graph6"]
        for record in records
        if record.get("survives_second_layer") is True
    )
    return source_manifest, graph6_values


def generate(source_path: Path, output_path: Path) -> dict[str, object]:
    source_manifest, graph6_values = _source_graphs(source_path)
    resolved_source = source_path.resolve()
    try:
        portable_path = str(resolved_source.relative_to(ROOT))
    except ValueError:
        portable_path = str(resolved_source)
    source = {
        "path": portable_path,
        "sha256": _file_sha256(source_path),
        "schema_version": source_manifest.get("schema_version"),
        "second_layer_survivors": len(graph6_values),
    }
    manifest = build_forced_a0_manifest(graph6_values, source)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    verify_forced_a0_manifest(manifest, graph6_values)
    return manifest


def verify(source_path: Path, output_path: Path) -> dict[str, object]:
    _, graph6_values = _source_graphs(source_path)
    manifest = json.loads(output_path.read_text(encoding="utf-8"))
    if manifest.get("source", {}).get("sha256") != _file_sha256(source_path):
        raise ValueError("forced-a=0 source digest does not match")
    return verify_forced_a0_manifest(manifest, graph6_values)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        summary = verify(args.source, args.output)
        mode = "verified"
    else:
        manifest = generate(args.source, args.output)
        summary = manifest["summary"]
        mode = "generated and verified"
    print(json.dumps({"mode": mode, "summary": summary}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
