#!/usr/bin/env python3
"""Freeze the exhaustive order-nine triangle-free X-graph quotient for record 17."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record17_x_graphs.json"
ORDER = 9
EXPECTED_COUNT = 1_897


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _stream_digest(words: Sequence[str]) -> str:
    return sha256(("".join(word + "\n" for word in words)).encode()).hexdigest()


def _nauty_version(executable: str) -> str:
    resolved = str(Path(executable).resolve())
    match = re.search(r"/nauty/([^/]+)/", resolved)
    return match.group(1) if match else "unknown"


def build() -> dict[str, object]:
    executable = shutil.which("geng")
    if executable is None:
        fallback = Path("/opt/homebrew/bin/geng")
        if fallback.is_file():
            executable = str(fallback)
        else:
            raise RuntimeError("nauty geng is required to rebuild the archive")
    completed = subprocess.run(
        [executable, "-q", "-t", str(ORDER)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    words = completed.stdout.splitlines()
    if len(words) != EXPECTED_COUNT:
        raise AssertionError(f"geng returned {len(words)} graphs, not 1897")
    document: dict[str, object] = {
        "schema_version": 1,
        "kind": "h9_14_record17_triangle_free_x_graph_stream",
        "generator": {
            "program": "nauty geng",
            "version": _nauty_version(executable),
            "command": "geng -q -t 9",
            "canonicalization": "one graph6 word per unlabeled graph",
            "trusted_scope": (
                "exhaustive isomorphism quotient only; decoding and every "
                "retained graph property are replayed independently"
            ),
        },
        "order": ORDER,
        "graph_count": len(words),
        "raw_stream_sha256": _stream_digest(words),
        "graph6": words,
    }
    document["certificate_sha256"] = _canonical_digest(document)
    return document


def verify(document: Mapping[str, object]) -> None:
    if document.get("schema_version") != 1 or document.get("kind") != (
        "h9_14_record17_triangle_free_x_graph_stream"
    ):
        raise ValueError("wrong record-17 X-graph archive metadata")
    words = document.get("graph6")
    if not isinstance(words, list) or not all(isinstance(word, str) for word in words):
        raise ValueError("record-17 archive has malformed graph6 words")
    if (
        document.get("order") != ORDER
        or document.get("graph_count") != EXPECTED_COUNT
        or len(words) != EXPECTED_COUNT
        or len(set(words)) != EXPECTED_COUNT
    ):
        raise ValueError("record-17 archive graph count changed")
    if document.get("raw_stream_sha256") != _stream_digest(words):
        raise ValueError("record-17 raw graph stream digest changed")
    if document.get("certificate_sha256") != _canonical_digest(document):
        raise ValueError("record-17 X-graph archive digest changed")


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    document = (
        json.loads(OUTPUT.read_text(encoding="utf-8"))
        if args.verify_only
        else build()
    )
    if not args.verify_only:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    verify(document)
    print(
        json.dumps(
            {
                "certificate_sha256": document["certificate_sha256"],
                "graph_count": document["graph_count"],
                "raw_stream_sha256": document["raw_stream_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
