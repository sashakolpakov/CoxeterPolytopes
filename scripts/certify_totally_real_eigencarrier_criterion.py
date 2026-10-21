#!/usr/bin/env python3
"""Freeze the general totally-real eigencarrier criterion and its baselines.

The general Hodge, shell-finiteness, and arithmetic-uniformity implications
are mathematical theorems proved in the accompanying report.  This script
checks their exact numerical specializations against the already frozen
Bugaenko, silver-K3, and golden-OG10 certificates.  It deliberately does not
turn those cited theorems into numerical claims.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.trace_lattice import INV_SQRT5


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/candidates/totally_real_eigencarrier_criterion.json"

INPUTS = {
    "bugaenko_polytope": ROOT / "data/known_examples/bugaenko_compact_h8.json",
    "bugaenko_k3": ROOT / "data/trace_lattices/bugaenko_k3_eigen_nef.json",
    "silver_k3": ROOT / "data/candidates/silver_k3_carrier.json",
    "silver_stabilizer": (
        ROOT / "data/candidates/silver_k3_norm2_infinite_stabilizer.json"
    ),
    "og10_carrier": ROOT / "data/candidates/og10_golden_h10_carrier.json",
    "og10_walls": ROOT / "data/candidates/og10_golden_complete_wall_shells.json",
    "og10_stabilizer": (
        ROOT / "data/candidates/og10_golden_h9_infinite_stabilizer.json"
    ),
}


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(name: str) -> dict[str, Any]:
    return json.loads(INPUTS[name].read_text(encoding="utf-8"))


def trace_rank(field_degree: int, field_rank: int) -> int:
    if field_degree < 2 or field_rank < 2:
        raise ValueError("an eigencarrier needs degree >= 2 and field rank >= 2")
    return field_degree * field_rank


def carrier_dimension_ceiling(b2: int, field_degree: int) -> int:
    """Largest n allowed by rank for a full-Picard H^n eigencarrier."""

    if b2 < 3 or field_degree < 2:
        raise ValueError("invalid Beauville rank or totally-real degree")
    return (b2 - 2) // field_degree - 1


def period_dimension(b2: int, trace_lattice_rank: int) -> int:
    """Complex dimension of the generic full-Picard period domain."""

    result = b2 - trace_lattice_rank - 2
    if result < 0:
        raise ValueError("the trace lattice leaves no positive period plane")
    return result


def _golden_shells(beauville_square: int, *, even: bool) -> list[Golden]:
    """Enumerate the exact totally-positive golden shell in a proved box.

    For alpha=-phi/sqrt(5), Tr(alpha*(a+b*phi))=-(a+b).  Positivity gives
    -s*phi < b < s/phi for s=-beauville_square, so |a|,|b| <= 4*s is a
    deliberately redundant finite box.
    """

    if beauville_square >= 0:
        raise ValueError("a wall square must be negative")
    alpha = -PHI * INV_SQRT5
    bound = -4 * beauville_square
    result: list[Golden] = []
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            if even and (a % 2 or b % 2):
                continue
            value = Golden(a, b)
            if (
                value.sign(0) > 0
                and value.sign(1) > 0
                and (alpha * value).trace() == beauville_square
            ):
                result.append(value)
    return sorted(result, key=lambda value: (value.a, value.b))


def _golden_json(values: list[Golden]) -> list[dict[str, str]]:
    return [value.to_json() for value in values]


def build() -> dict[str, Any]:
    bugaenko_polytope = _read("bugaenko_polytope")
    bugaenko = _read("bugaenko_k3")
    silver = _read("silver_k3")
    silver_stabilizer = _read("silver_stabilizer")
    og10 = _read("og10_carrier")
    og10_walls = _read("og10_walls")
    og10_stabilizer = _read("og10_stabilizer")

    raw_minus_two = _golden_shells(-2, even=False)
    even_minus_two = _golden_shells(-2, even=True)
    even_minus_four = _golden_shells(-4, even=True)

    silver_model = silver["carrier_models"][0]
    checks = {
        "rank_formula_bugaenko": trace_rank(2, 9) == 18,
        "rank_formula_silver": trace_rank(2, 10) == 20,
        "rank_formula_og10": trace_rank(2, 11) == 22,
        "k3_quadratic_ceiling_H9": carrier_dimension_ceiling(22, 2) == 9,
        "og10_quadratic_ceiling_H10": carrier_dimension_ceiling(24, 2) == 10,
        "bugaenko_period_dimension_two": period_dimension(22, 18) == 2,
        "silver_period_dimension_zero": period_dimension(22, 20) == 0,
        "og10_period_dimension_zero": period_dimension(24, 22) == 0,
        "bugaenko_exact_k3_embedding": (
            bugaenko["trace_lattice_primitive_k3_embedding_verified"] is True
        ),
        "bugaenko_complete_compact_chamber": (
            bugaenko_polytope["dimension"] == 8
            and bugaenko_polytope["certificate_kind"] == "compact_reflective"
            and bugaenko["vinberg_completion_verified"] is True
            and bugaenko["compact_vertex_count"] == 41
            and bugaenko["eigen_nef_realization"]
            ["restricted_chamber_equals_bugaenko_chamber"]
            is True
        ),
        "golden_raw_minus_two_shell_exact": [str(value) for value in raw_minus_two]
        == ["1 + 1*phi", "2", "3 - 1*phi", "4 - 2*phi", "5 - 3*phi"],
        "golden_even_minus_two_shell_exact": [
            str(value) for value in even_minus_two
        ]
        == ["2", "4 - 2*phi"],
        "golden_even_minus_four_shell_exact": [
            str(value) for value in even_minus_four
        ]
        == [
            "2 + 2*phi",
            "4",
            "6 - 2*phi",
            "8 - 4*phi",
            "10 - 6*phi",
        ],
        "og10_three_visible_wall_orbits": (
            og10_walls["complete_restricted_search_norms"]
            == [
                {"a": "2", "b": "0"},
                {"a": "4", "b": "0"},
                {"a": "6", "b": "-2"},
            ]
        ),
        "silver_admissible_field_signatures": (
            silver_model["field_identity_inertia"] == [9, 1, 0]
            and silver_model["field_conjugate_inertia"] == [10, 0, 0]
            and silver_model["trace_inertia"] == [1, 19, 0]
        ),
        "silver_primitive_k3_embedding": (
            silver["k3_embedding"]["primitive_embedding_verified"] is True
        ),
        "silver_infinite_stabilizer_obstruction": (
            silver_stabilizer["verified"] is True
            and silver_stabilizer["arithmetic_conclusion"]
            ["complete_chamber_stabilizer"]
            == "infinite"
            and silver_stabilizer["arithmetic_conclusion"]
            ["complete_chamber_hyperbolic_volume"]
            == "infinite"
        ),
        "og10_admissible_field_signatures": (
            og10["field_signatures"]
            == {"identity": [10, 1, 0], "conjugate": [11, 0, 0]}
            and og10["trace_lattice"]["signature"] == [1, 21]
        ),
        "og10_primitive_carrier": (
            og10["verified"] is True
            and og10["checks"]["period_complement_positive_rank_2"] is True
        ),
        "og10_infinite_stabilizer_obstruction": (
            og10_stabilizer["verified"] is True
            and og10_stabilizer["group_theoretic_conclusion"]
            ["complete_chamber_stabilizer"]
            == "infinite"
            and og10_stabilizer["ambient_consequence"]
            ["golden_H10_chamber_volume"]
            == "infinite"
        ),
    }

    result: dict[str, Any] = {
        "schema_version": "1.0.0",
        "kind": "totally_real_eigencarrier_criterion_and_baseline_audit",
        "theorem_scope": {
            "weak_full_picard_criterion": (
                "For a fixed admissible trace datum and deformation type, a "
                "weak full-Picard eigencarrier exists iff its even trace "
                "lattice admits an orientation-compatible primitive embedding "
                "in the Beauville lattice."
            ),
            "finite_shell_reduction": (
                "For each allowed negative Beauville wall square, visibility "
                "forces h(delta,delta) into a finite set of totally positive "
                "elements cut out by one trace equation."
            ),
            "uniformity_obstruction": (
                "An admissible form with every non-distinguished real place "
                "definite is F-anisotropic.  For a complete invariant integral "
                "reflection arrangement, finite chamber volume is equivalent "
                "to compactness and to finiteness of its arithmetic chamber "
                "stabilizer."
            ),
            "not_automatic": [
                "primitive embedding of an abstract Coxeter Gram without a fixed integral scaling",
                "MBM/effective-wall membership",
                "equality with a nef, movable, or Kahler chamber",
                "Coxeter angles for arbitrary projected walls",
                "finite-sidedness, finite volume, or compactness",
            ],
        },
        "rank_ladder": [
            {
                "name": "Bugaenko golden H8",
                "field_degree": 2,
                "field_rank": 9,
                "trace_rank": 18,
                "carrier_type": "K3",
                "b2": 22,
                "period_dimension": 2,
                "chamber_status": "complete compact 11-facet eigen-nef chamber",
            },
            {
                "name": "silver H9",
                "field_degree": 2,
                "field_rank": 10,
                "trace_rank": 20,
                "carrier_type": "K3",
                "b2": 22,
                "period_dimension": 0,
                "chamber_status": "complete norm-two chamber has infinite volume",
            },
            {
                "name": "golden H10",
                "field_degree": 2,
                "field_rank": 11,
                "trace_rank": 22,
                "carrier_type": "OG10",
                "b2": 24,
                "period_dimension": 0,
                "chamber_status": "complete norm-two movable chamber has infinite volume",
            },
        ],
        "dimension_ceilings": {
            "K3_degree_2": 9,
            "OG10_degree_2": 10,
            "formula": "floor((b2-2)/[F:Q])-1",
        },
        "golden_shell_replay": {
            "trace_scalar": (-PHI * INV_SQRT5).to_json(),
            "raw_square_minus_2": _golden_json(raw_minus_two),
            "even_square_minus_2": _golden_json(even_minus_two),
            "even_square_minus_4": _golden_json(even_minus_four),
        },
        "baseline_outcome": {
            "positive": (
                "Bugaenko H8 is an exact compact K3 eigen-nef realization."
            ),
            "negative_design_filter": (
                "The silver K3 H9 and golden OG10 H10 intrinsic norm-two "
                "carriers are eliminated from finite volume by infinite "
                "complete-chamber stabilizers, not by finite Vinberg prefixes."
            ),
        },
        "external_theorem_inputs": [
            "surjectivity of the marked K3/hyperkahler period map and the projectivity criterion",
            "K3 Weyl transitivity or a deformation-type-specific Kahler/movable chamber theorem",
            "boundedness/classification of MBM wall squares for the fixed deformation type",
            "Borel-Harish-Chandra arithmeticity and the anisotropic compactness criterion",
            "Vinberg reflection-chamber theory",
        ],
        "inputs": {
            name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": _file_sha256(path),
            }
            for name, path in INPUTS.items()
        },
        "checks": checks,
        "verified": all(checks.values()),
    }
    result["certificate_sha256"] = _canonical_sha256(result)
    return result


def verify(archived: dict[str, Any]) -> dict[str, Any]:
    fresh = build()
    if archived != fresh:
        raise AssertionError("the archived eigencarrier criterion has drifted")
    if not fresh["verified"]:
        failed = [name for name, value in fresh["checks"].items() if not value]
        raise AssertionError(f"failed eigencarrier checks: {failed}")
    return {
        "certificate_sha256": fresh["certificate_sha256"],
        "rank_ladder": [row["trace_rank"] for row in fresh["rank_ladder"]],
        "dimension_ceilings": fresh["dimension_ceilings"],
        "baseline_statuses": [row["chamber_status"] for row in fresh["rank_ladder"]],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    if args.write:
        OUTPUT.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
