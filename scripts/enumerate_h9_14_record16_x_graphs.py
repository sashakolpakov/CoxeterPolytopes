#!/usr/bin/env python3
"""Freeze the tiny triangle-free graph streams used by record 16."""

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
OUTPUT = ROOT / "data/enumerations/h9_14_record16_x_graphs.json"
ORDERS = (6, 7, 8)
EXPECTED_COUNTS = {6: 38, 7: 107, 8: 410}


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
        raise RuntimeError("nauty geng is required to rebuild the X-graph archive")
    streams: dict[str, object] = {}
    for order in ORDERS:
        command = [executable, "-q", "-t", str(order)]
        completed = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        words = completed.stdout.splitlines()
        if len(words) != EXPECTED_COUNTS[order]:
            raise AssertionError(
                f"geng order {order} returned {len(words)} graphs"
            )
        streams[str(order)] = {
            "command": f"geng -q -t {order}",
            "graph_count": len(words),
            "raw_stream_sha256": _stream_digest(words),
            "graph6": words,
        }
    document: dict[str, object] = {
        "schema_version": 1,
        "kind": "h9_14_record16_triangle_free_x_graph_streams",
        "generator": {
            "program": "nauty geng",
            "version": _nauty_version(executable),
            "canonicalization": "one graph6 word per unlabeled graph",
            "trusted_scope": (
                "exhaustive generation only; every retained graph property "
                "is independently replayed"
            ),
        },
        "streams": streams,
    }
    document["certificate_sha256"] = _canonical_digest(document)
    return document


def verify(document: Mapping[str, object]) -> None:
    if document.get("schema_version") != 1 or document.get("kind") != (
        "h9_14_record16_triangle_free_x_graph_streams"
    ):
        raise ValueError("wrong record-16 X-graph archive metadata")
    streams = document.get("streams")
    if not isinstance(streams, Mapping):
        raise ValueError("record-16 X-graph archive lacks streams")
    for order, expected in EXPECTED_COUNTS.items():
        stream = streams.get(str(order))
        if not isinstance(stream, Mapping):
            raise ValueError(f"record-16 archive lacks order {order}")
        words = stream.get("graph6")
        if not isinstance(words, list) or not all(
            isinstance(word, str) for word in words
        ):
            raise ValueError("record-16 archive has malformed graph6 words")
        if len(words) != expected or stream.get("graph_count") != expected:
            raise ValueError("record-16 archive graph count changed")
        if len(set(words)) != expected:
            raise ValueError("record-16 archive contains duplicate graph6 words")
        if stream.get("raw_stream_sha256") != _stream_digest(words):
            raise ValueError("record-16 raw graph stream digest changed")
    if document.get("certificate_sha256") != _canonical_digest(document):
        raise ValueError("record-16 X-graph archive digest changed")


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        document = build()
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    verify(document)
    print(
        json.dumps(
            {
                "certificate_sha256": document["certificate_sha256"],
                "counts": {
                    key: value["graph_count"]
                    for key, value in document["streams"].items()
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
