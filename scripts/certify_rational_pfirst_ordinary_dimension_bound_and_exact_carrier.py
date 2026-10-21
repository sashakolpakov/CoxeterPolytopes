#!/usr/bin/env python3
"""Replay the ordinary P-first n<=9 bound and the two carrier branches."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_rational_pfirst_anticanonical_cusp_saturation_closure import (
    payload as cusp_payload,
)
from scripts.certify_rational_pfirst_timelike_anticanonical_vertex_closure import (
    payload as timelike_payload,
)

DEPENDENCIES = {
    "timelike_report": {
        "path": "results/rational_pfirst_timelike_anticanonical_vertex_closure.md",
        "sha256": "6450ab0629c539fcb3fff00e6aa2bdd0ddf194f17bdaaa2c1f2265d2b4aebb29",
    },
    "timelike_hostile_audit": {
        "path": (
            "results/"
            "rational_pfirst_timelike_anticanonical_vertex_closure_hostile_audit.md"
        ),
        "sha256": "2d7eddae058f81f06fded910b4d8ed45219b30e011af5be6a232f20655df0499",
    },
    "cusp_report": {
        "path": "results/rational_pfirst_anticanonical_cusp_saturation_closure.md",
        "sha256": "9d4142c8c4647bfced451ac1ff29184a6a278d9a51cdb1a0638276f3eb813032",
    },
    "cusp_hostile_audit": {
        "path": (
            "results/"
            "rational_pfirst_anticanonical_cusp_saturation_closure_hostile_audit.md"
        ),
        "sha256": "1cc36249914c82ba42145e64cb0e1a43a888d4ac0c02d6b02084cade3ea498d8",
    },
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def rank_record(n: int) -> dict[str, object]:
    canonical = (-3,) + (1,) * n
    anticanonical = tuple(-entry for entry in canonical)
    square = 9 - n

    assert dot(canonical, canonical) == square
    assert dot(anticanonical, anticanonical) == square
    assert all(entry % 2 for entry in anticanonical)
    assert anticanonical[0] > 0

    if n <= 8:
        branch = "timelike_exact_carrier"
        finite_volume_possible = True
        exact_carrier = True
    elif n == 9:
        branch = "characteristic_cusp_exact_carrier"
        finite_volume_possible = True
        exact_carrier = True
    else:
        branch = "spacelike_wall_cone_infinite_volume"
        finite_volume_possible = False
        exact_carrier = False
        assert square < 0

    return {
        "hyperbolic_dimension": n,
        "rank": n + 1,
        "canonical": list(canonical),
        "anticanonical": list(anticanonical),
        "canonical_square": square,
        "ordinary_facet_species": [
            {
                "curve_square": -1,
                "canonical_pairing": -1,
                "anticanonical_wall_slack": 1,
            },
            {
                "curve_square": -2,
                "canonical_pairing": 0,
                "anticanonical_wall_slack": 0,
            },
        ],
        "escape_polynomial": {
            "expression": "A^2 + 2*t*(A.F) + t^2*(9-n)",
            "quadratic_coefficient": square,
            "strict_wall_slack_expression": "A.C + t*(F.C) > 0",
        },
        "branch": branch,
        "finite_volume_possible_in_ordinary_lane": finite_volume_possible,
        "exact_rational_carrier_proved": exact_carrier,
    }


def verify_positive_dependencies() -> dict[str, object]:
    for item in DEPENDENCIES.values():
        observed = sha256_file(ROOT / str(item["path"]))
        assert observed == item["sha256"], (
            item["path"],
            item["sha256"],
            observed,
        )

    timelike = timelike_payload()
    assert timelike["exceptional_counts_covered"] == list(range(3, 9))
    assert all(
        record["flex_cluster_carrier"]["full_nef_return"]
        for record in timelike["records"]
    )
    assert timelike["theorem_scope"][
        "classification_free_exact_carrier"
    ] == "PROVED"

    cusp = cusp_payload()
    assert cusp["saturated_root_squares"] == [-1] + [-2] * 9
    assert cusp["characteristic_ideal_vertex_count"] == 1
    assert cusp["explicit_realization"]["full_nef_chamber_returned"] is True
    assert cusp["scope"]["hamming_saturation_outcome"] == (
        "MAXIMAL_AND_REALIZABLE"
    )

    return {
        "dependency_hashes_verified": True,
        "timelike_exceptional_counts": list(range(3, 9)),
        "timelike_all_exact_full_nef_returns": True,
        "rank_ten_characteristic_cusp_exact_full_nef_return": True,
    }


def payload() -> dict[str, object]:
    positive = verify_positive_dependencies()
    records = [rank_record(n) for n in range(3, 16)]

    assert [row["hyperbolic_dimension"] for row in records[:7]] == list(
        range(3, 10)
    )
    assert all(
        row["exact_rational_carrier_proved"] for row in records[:7]
    )
    assert all(
        not row["finite_volume_possible_in_ordinary_lane"]
        for row in records[7:]
    )

    result: dict[str, object] = {
        "schema": (
            "rational_pfirst_ordinary_dimension_bound_and_exact_carrier/v1"
        ),
        "standard_marking": {
            "lattice": "I_(1,n)",
            "canonical": "K=(-3,1,...,1)",
            "canonical_square": "9-n",
            "hyperbolic_dimensions": "n>=3",
        },
        "wall_cone_lemma": {
            "spacelike_vector_in_homogeneous_wall_cone": (
                "forces a strict open ideal patch and infinite volume"
            ),
            "past_causal_vector_in_homogeneous_wall_cone": (
                "forces a spacelike positive combination"
            ),
            "future_null_vector_in_finite_volume_wall_cone": (
                "is an extreme ideal ray"
            ),
            "surface_or_reflection_maximality_used": False,
        },
        "ordinary_color_identity": {
            "minus_one_wall_F_pairing": 1,
            "minus_two_wall_F_pairing": 0,
            "F_square": "9-n",
        },
        "sharp_dimension_bound": {
            "hyperbolic_dimension_at_most": 9,
            "rank_at_most": 10,
            "n_at_least_10": "IMPOSSIBLE_BY_FINITE_VOLUME",
        },
        "records": records,
        "positive_branch_composition": positive,
        "dependencies": DEPENDENCIES,
        "theorem_scope": {
            "ordinary_minus_one_minus_two_standard_marked_lane": "CLOSED",
            "relative_reflection_maximality_needed_for_carrier": True,
            "dimension_bound_needs_maximality": False,
            "canonical_escape_is_automatic_after_finite_volume": True,
            "exact_carriers_for_n_3_through_9": True,
            "twice_primitive_minus_four_adapted_host_lane": "OPEN",
            "nonstandard_characteristic_orbits": "OPEN",
        },
    }
    canonical = json.dumps(
        result, sort_keys=True, separators=(",", ":")
    ).encode()
    result["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    result = payload()
    default = ROOT / (
        "data/enumerations/"
        "rational_pfirst_ordinary_dimension_bound_and_exact_carrier.json"
    )
    target = args.output or default
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        assert target.read_text(encoding="utf-8") == rendered
    elif args.output:
        target.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
