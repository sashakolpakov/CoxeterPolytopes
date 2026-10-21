#!/usr/bin/env python3
"""Exact finite checks for the canonical-orthogonal/K3-cover rank theorem.

The surface-classification, double-cover, and complex K3 monodromy inputs are
theoretical.  This script freezes the rank arithmetic, the two root-lift
identities, and the provenance of every calibrated carrier.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data/enumerations/canonical_orthogonal_k3_cover_calibration.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trace_rank(field_degree: int, hyperbolic_dimension: int) -> int:
    assert field_degree >= 1
    assert hyperbolic_dimension >= 2
    return field_degree * (hyperbolic_dimension + 1)


def carrier_ceiling(characteristic: int) -> int:
    return 20 if characteristic == 0 else 22


def dimension_ceiling(field_degree: int, characteristic: int) -> int:
    return carrier_ceiling(characteristic) // field_degree - 1


def coble_rank(boundary_components: int) -> int:
    assert boundary_components >= 1
    return 10 + boundary_components


def coble_dimension(boundary_components: int) -> int:
    return coble_rank(boundary_components) - 1


def coble_k3_boundary_ceiling(characteristic: int) -> int:
    """Boundary ceiling supplied by the K3 Picard-rank argument alone."""

    return carrier_ceiling(characteristic) - 10


def coble_k3_dimension_ceiling(characteristic: int) -> int:
    return coble_dimension(coble_k3_boundary_ceiling(characteristic))


def known_coble_boundary_ceiling(characteristic: int) -> int:
    """General AG bound for strong Coble surfaces (characteristic != 2).

    This is a frozen theoretical input, not proved by the finite arithmetic
    replay.  Kondo, Section 3, cites Dolgachev--Kondo, Corollary 9.1.5.
    """

    assert characteristic != 2
    return 10


def known_coble_dimension_ceiling(characteristic: int) -> int:
    return coble_dimension(known_coble_boundary_ceiling(characteristic))


def single_root_reflection(root_square: int = -2) -> sp.Matrix:
    """Reflection on the rank-one root line."""

    gram = sp.Matrix([[root_square]])
    root = sp.Matrix([1])
    return sp.eye(1) - sp.Rational(2) * root * (root.T * gram) / root_square


def split_root_restriction_identity() -> bool:
    """Check s_c s_c' restricts to s_(c+c') on the invariant subspace.

    Use a formal invariant vector x with x.c=x.c'=a.  The value of x^2 is
    irrelevant.  The computation is made in the span of x,c,c'.
    """

    a, z = sp.symbols("a z")
    gram = sp.Matrix([[z, a, a], [a, -2, 0], [a, 0, -2]])
    c = sp.Matrix([0, 1, 0])
    cp = sp.Matrix([0, 0, 1])
    x = sp.Matrix([1, 0, 0])
    v = c + cp

    def reflection(root: sp.Matrix) -> sp.Matrix:
        square = (root.T * gram * root)[0]
        return sp.eye(3) - sp.Rational(2) * root * (root.T * gram) / square

    paired = reflection(c) * reflection(cp)
    summed = reflection(v)
    return paired * x == summed * x and paired * v == summed * v


def certify_pullback_pairing() -> None:
    """Freeze the scaled extended-Picard pullback calculation.

    If h=B/2 downstairs and R is the ramification curve upstairs, then
    pi^*B=2R.  For x.h=b and x^2=a, the source form multiplied by two is
    exactly the target form on pi^*x,R.
    """

    a, b = sp.symbols("a b")
    scaled_source = 2 * sp.Matrix([[a, b], [b, -1]])
    target = sp.Matrix([[2 * a, 2 * b], [2 * b, -2]])
    assert scaled_source == target


def certify_sources(data: dict) -> None:
    for relative, expected in data["source_pins"].items():
        assert sha256(ROOT / relative) == expected, relative


def certify_frozen_inputs() -> None:
    bugaenko = json.loads(
        (ROOT / "data/trace_lattices/bugaenko_k3_eigen_nef.json").read_text()
    )
    factor = bugaenko["eigen_nef_realization"]["field_reflection_factorization"]
    assert all(factor.values())
    assert bugaenko["trace_lattice"] == "<2> + <-2> + E8(-1) + E8(-1)"

    type_ii = json.loads(
        (ROOT / "data/enumerations/coble_type_ii_b3_plane_inverse.json").read_text()
    )
    assert type_ii["blowup_surface"]["strong_Coble"]
    assert type_ii["blowup_surface"]["Picard_rank"] == 11
    assert type_ii["same_polytope_return"]["exact_same_marked_H10_30"]

    type_i = json.loads(
        (
            ROOT
            / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json"
        ).read_text()
    )
    invariants = type_i["full_NS_lattice"][
        "doubled_extended_Picard_nikulin_invariants"
    ]
    assert invariants == {"a": 10, "delta": 1, "r": 12}
    assert "H11" in type_i["claim"]["theorem"]

    type_vi = json.loads(
        (
            ROOT
            / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
        ).read_text()
    )
    assert type_vi["full_NS_lattice"]["rank"] == 15
    assert type_vi["full_NS_lattice"]["hyperbolic_dimension"] == 14

    type_mii = json.loads(
        (ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json").read_text()
    )
    stratum = type_mii["AG_native_stratum"]
    assert stratum["Picard_rank"] == 18
    assert stratum["boundary_component_count"] == 8
    assert stratum["hyperbolic_dimension"] == 17
    assert "supersingular" in stratum["covering_K3"]

    enriques = json.loads(
        (
            ROOT
            / "data/known_examples/enriques_finite_automorphism_nef_chambers.json"
        ).read_text()
    )
    assert enriques["lattice"]["rank"] == 10
    assert len(enriques["types"]) == 7
    assert enriques["global_checks"]["facet_count_multiset"] == [
        12,
        12,
        20,
        20,
        20,
        20,
        20,
    ]


def certify() -> None:
    data = json.loads(LEDGER.read_text())
    certify_sources(data)
    certify_frozen_inputs()

    assert single_root_reflection() == sp.Matrix([[-1]])
    assert split_root_restriction_identity()
    certify_pullback_pairing()

    for row in data["rows"]:
        assert row["trace_rank"] == trace_rank(row["d"], row["n"])
        assert row["trace_rank"] <= carrier_ceiling(row["characteristic"])
        assert row["canonical_orthogonality"].startswith("automatic")
        assert row["algebraic_correspondences"]
        if "boundary_components" in row:
            boundary_count = row["boundary_components"]
            assert row["trace_rank"] == coble_rank(boundary_count)
            assert row["n"] == coble_dimension(boundary_count)

    assert dimension_ceiling(1, 0) == 19
    assert dimension_ceiling(2, 0) == 9
    assert dimension_ceiling(1, 3) == 21
    assert coble_k3_boundary_ceiling(0) == 10
    assert coble_k3_dimension_ceiling(0) == 19
    assert coble_k3_boundary_ceiling(3) == 12
    assert coble_k3_dimension_ceiling(3) == 21
    assert known_coble_boundary_ceiling(0) == 10
    assert known_coble_boundary_ceiling(3) == 10
    assert known_coble_dimension_ceiling(0) == 19
    assert known_coble_dimension_ceiling(3) == 19

    summary = data["summary"]
    assert summary["all_rows_pass_canonical_orthogonal_rank_hypothesis"]
    assert summary["all_rows_have_algebraic_correspondence_lifts"]
    assert not summary["universal_raw_coxeter_to_carrier_compiler"]

    coble = data["strong_coble_cover_theorem"]
    assert coble["canonical_cover_cohomological_boundary_ceiling"] == {
        "complex": 10,
        "positive_characteristic_not_2": 12,
    }
    assert coble["canonical_cover_cohomological_dimension_ceiling"] == {
        "complex": 19,
        "positive_characteristic_not_2": 21,
    }
    assert coble["best_known_classwide_boundary_ceiling"] == 10
    assert coble["best_known_classwide_dimension_ceiling"] == 19
    assert not coble["uses_finite_automorphism_classification"]

    sources = data["theoretical_sources"]
    assert sources["kondo_coble_2021"]["arxiv"] == "2102.00591"
    assert sources["huybrechts_k3_lectures"]["doi"] == (
        "10.1017/CBO9781316594193"
    )
    assert not sources["surface_classification_input"]["machine_verified"]


def main() -> int:
    certify()
    print("PASS: canonical-orthogonal and K3-cover rank calibration is exact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
