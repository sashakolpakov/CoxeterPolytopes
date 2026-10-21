#!/usr/bin/env python3
"""Run the exact global Vinberg diagnostic on the 1,076-wall silver prefix.

The input is the frozen ``x=5+3sqrt(2)`` manifest.  This script consumes its
literal production-order wall list, rebuilds the exact Gram matrix, and runs
the sparse rank-eight edge, finite-endpoint, connected-affine, and full
parabolic-neighbourhood audit.  Failure is diagnostic only: later walls may
close or cut off a currently open edge.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Sequence

from coxeter_polytopes.verify.silver_vinberg_volume import (
    audit_silver_vinberg_wall_gram,
)
from scripts.diagnose_silver_k3_726_global_vinberg import (
    ROOT,
    _audit_record,
    _canonical_digest,
    _decode_roots,
    _note,
    build_wall_gram,
)


DEFAULT_INPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_5_plus_3sqrt2.json"
)
DEFAULT_OUTPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_1076_global_diagnostic.json"
)
_AFFINE_PROGRESS = re.compile(r"affine growth order (\d+); next=(\d+), affine=(\d+)")


def _histogram_sha256(histogram: dict[str, int]) -> str:
    return sha256(
        json.dumps(histogram, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build_diagnostic(input_path: Path = DEFAULT_INPUT) -> dict[str, object]:
    document = json.loads(input_path.read_text(encoding="utf-8"))
    input_certificate = document.get("certificate_sha256")
    if not isinstance(input_certificate, str):
        raise ValueError("input certificate has no certificate_sha256")
    roots = _decode_roots(document)
    if len(roots) != 1076:
        raise ValueError("the x=5+3sqrt(2) manifest does not contain 1,076 roots")
    _note(f"decoded {len(roots)} literal production-order roots")
    wall_gram = build_wall_gram(roots)
    pairing_histogram = dict(
        sorted(
            Counter(
                str(wall_gram[left][right])
                for left in range(len(wall_gram))
                for right in range(left)
            ).items()
        )
    )
    _note(
        "wall Gram complete; distinct off-diagonal values="
        f"{len(pairing_histogram)}"
    )

    affine_growth: list[int] = []

    def progress(message: str) -> None:
        _note(message)
        match = _AFFINE_PROGRESS.fullmatch(message)
        if match is not None:
            order, following, affine = map(int, match.groups())
            if order != len(affine_growth) + 1 or affine != 0:
                raise AssertionError("the affine-growth transcript is inconsistent")
            affine_growth.append(following)

    audit = audit_silver_vinberg_wall_gram(wall_gram, 9, progress=progress)
    record = _audit_record(
        audit,
        input_path=input_path,
        input_certificate=input_certificate,
    )
    record["schema_version"] = "1.0.0"
    record["kind"] = "silver_k3_1076_global_vinberg_prefix_diagnostic"
    record["connected_affine_growth_next_by_order"] = affine_growth
    record["off_diagonal_pairing_histogram"] = pairing_histogram
    record["off_diagonal_pairing_histogram_sha256"] = _histogram_sha256(
        pairing_histogram
    )
    record.pop("diagnostic_sha256")
    record["diagnostic_sha256"] = _canonical_digest(record)
    return record


def _summary(record: dict[str, object]) -> dict[str, object]:
    return {
        key: record[key]
        for key in (
            "diagnostic_sha256",
            "input_certificate_sha256",
            "facet_count",
            "backtracking_positive_prefix_visits_by_order_1_through_8",
            "elliptic_rank_eight_edge_count",
            "elliptic_rank_nine_finite_vertex_count",
            "finite_endpoint_histogram",
            "connected_affine_growth_next_by_order",
            "connected_affine_subdiagram_count",
            "ideal_vertex_count",
            "finite_volume_open_edge_count",
            "invalid_finite_angle_pair_count",
            "off_diagonal_pairing_histogram_sha256",
            "compact_stopping_available",
            "finite_volume_stopping_available",
        )
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args(argv)
    record = build_diagnostic(args.input)
    if args.write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(record if args.full else _summary(record), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
