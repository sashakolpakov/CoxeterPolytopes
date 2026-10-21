#!/usr/bin/env python3
"""Exact source-stabilizer no-go for Shimada's complete ``X_1`` chamber.

Every subgroup of the genuine 720-element ample stabilizer whose common fixed
lattice has rank at least eleven is exhausted.  The unique high class is the
symplectic involution class.  Its complete 1,094-wall restriction has strict
facet witnesses but fails the Coxeter pair-angle gate, so no volume claim is
made.  The full abstract weighted-wall graph group is intentionally outside
this certificate's denominator.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import sympy as sp

from scripts.certify_shimada_aut3singk3_x2_rational_involution_h10_37 import (
    _common_fixed_rank,
    _fold,
    _matrix,
    _rank,
)
from scripts.ingest_shimada_aut3singk3_x1_1094_stock import (
    SCHEMA as STOCK_SCHEMA,
    validate_stock,
)


ROOT = Path(__file__).resolve().parents[1]
STOCK = ROOT / "data/enumerations/shimada_aut3singk3_x1_1094_wall_stock.json"
OUTPUT = ROOT / "data/enumerations/shimada_aut3singk3_x1_symplectic_h11_no_go.json"
SCHEMA = "coxeter.shimada.aut3singk3.x1.symplectic_h11_source_group_no_go.v1"
STOCK_SHA256 = "8fab28bffb6c617b09d81fd9735581abd39a28f89cd6d760bedaefbbdf8b330e"
IDENTITY_INDEX = 7
SYMPLECTIC_REPRESENTATIVE = 3
EXPECTED_FIXED_RANK_HISTOGRAM = {
    2: 144,
    4: 324,
    6: 90,
    8: 80,
    10: 36,
    12: 45,
    20: 1,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _public_fold(fold: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in fold.items() if not key.startswith("_")}


def _high_subgroup_denominator(group: Sequence[sp.Matrix]) -> dict[str, Any]:
    identity = sp.eye(20)
    fixed_ranks = [20 - _rank(matrix.T - identity) for matrix in group]
    histogram = Counter(fixed_ranks)
    if dict(sorted(histogram.items())) != EXPECTED_FIXED_RANK_HISTOGRAM:
        raise ArithmeticError("the X_1 fixed-rank histogram changed")

    arrays = np.asarray(
        [[[int(value) for value in row] for row in matrix.tolist()] for matrix in group],
        dtype=np.int64,
    )
    maximum_entry = int(np.max(np.abs(arrays)))
    product_bound = 20 * maximum_entry**2
    if product_bound >= 2**63:
        raise ArithmeticError("int64 product bound failed")
    index_of = {
        tuple(int(value) for value in arrays[index].ravel()): index
        for index in range(len(group))
    }
    if len(index_of) != 720:
        raise ArithmeticError("duplicate group matrix")

    high = [
        index
        for index, rank in enumerate(fixed_ranks)
        if rank >= 11 and index != IDENTITY_INDEX
    ]
    identity_array = np.eye(20, dtype=np.int64)
    if not (
        len(high) == 45
        and all(fixed_ranks[index] == 12 for index in high)
        and all(np.array_equal(arrays[index] @ arrays[index], identity_array) for index in high)
    ):
        raise ArithmeticError("the high symplectic involution denominator changed")

    representative = arrays[SYMPLECTIC_REPRESENTATIVE]
    centralizer = [
        index
        for index, matrix in enumerate(arrays)
        if np.array_equal(matrix @ representative, representative @ matrix)
    ]
    if len(centralizer) != 16 or 720 // len(centralizer) != 45:
        raise ArithmeticError("the 45 high involutions are not one conjugacy class")

    high_set = set(high)
    v4_groups: set[frozenset[int]] = set()
    for position, left in enumerate(high):
        for right in high[position + 1 :]:
            product = arrays[left] @ arrays[right]
            product_index = index_of.get(tuple(int(value) for value in product.ravel()))
            if product_index is None:
                raise ArithmeticError("the published stabilizer list is not closed on a high pair")
            if product_index not in high_set:
                continue
            if not np.array_equal(product, arrays[right] @ arrays[left]):
                raise ArithmeticError("a high involutive product failed to commute")
            subgroup = frozenset((IDENTITY_INDEX, left, right, product_index))
            if len(subgroup) == 4:
                v4_groups.add(subgroup)

    v4_ranks = [
        _common_fixed_rank([group[index] for index in sorted(subgroup - {IDENTITY_INDEX})])
        for subgroup in v4_groups
    ]
    if len(v4_groups) != 30 or Counter(v4_ranks) != {8: 30}:
        raise ArithmeticError("the high-element V4 closure changed")

    return {
        "fixed_rank_histogram": {str(key): value for key, value in sorted(histogram.items())},
        "high_nonidentity_element_count": len(high),
        "unique_high_element_class": {
            "surface_type": "symplectic involution",
            "cardinality": 45,
            "fixed_rank": 12,
            "hyperbolic_dimension": 11,
            "representative_zero_based": SYMPLECTIC_REPRESENTATIVE,
            "centralizer_order": len(centralizer),
            "conjugacy_class_size": 45,
        },
        "rational_involution_class": {
            "cardinality": 36,
            "fixed_rank": 10,
            "hyperbolic_dimension": 9,
            "excluded_by_dimension_cut": True,
        },
        "high_closed_V4_count": len(v4_groups),
        "high_closed_V4_common_fixed_rank_histogram": {"8": 30},
        "only_C2_can_be_a_proper_fixed_rank_at_least_11_subgroup": True,
        "proof_of_exhaustion": (
            "If Fix(H) has rank at least 11, each element of H has individual fixed "
            "rank at least 11.  The exact 720-matrix census leaves only 45 rank-12 "
            "involutions, forming one class by the order-16 centralizer.  If H has two "
            "distinct nonidentity elements, their product is another high element; the "
            "complete high-pair closure gives one of 30 V4 groups, all of common fixed "
            "rank 8.  Thus every proper source-group section above H9 is conjugate to "
            "the displayed symplectic C2 section."
        ),
        "exact_int64_audit": {
            "maximum_matrix_entry_absolute": maximum_entry,
            "single_product_entry_absolute_bound": product_bound,
            "below_signed_int64_limit": True,
        },
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(STOCK) != STOCK_SHA256:
        raise ArithmeticError("the frozen X_1 stock changed")
    stock = json.loads(STOCK.read_text(encoding="utf-8"))
    validate_stock(stock)
    if stock["schema"] != STOCK_SCHEMA:
        raise ArithmeticError("wrong stock schema")

    gram = _matrix(stock["surface"]["published_neron_severi_gram"])
    ample = sp.Matrix(1, 20, stock["surface"]["published_ample_vector_a1"])
    walls = [
        sp.Matrix(1, 20, wall)
        for orbit in stock["complete_induced_chamber"]["wall_orbits"]
        for wall in orbit["primitive_outward_wall_vectors_in_published_NS_basis"]
    ]
    group = [
        _matrix(matrix)
        for matrix in stock["ample_stabilizer"][
            "all_integral_20_by_20_matrices_in_supplement_order"
        ]
    ]
    denominator = _high_subgroup_denominator(group)
    fold = _fold(
        SYMPLECTIC_REPRESENTATIVE,
        "symplectic",
        gram,
        ample,
        walls,
        group,
    )
    if not (
        fold["wall_orbit_count"] == 594
        and fold["positive_ray_deduplicated_restriction_count"] == 580
        and fold["restriction_square_sign_histogram"]
        == {"negative": 162, "positive": 390, "zero": 28}
        and fold["strict_facet_count"] == 162
        and fold["facet_gram_rank"] == 12
        and fold["positive_wall_gram_inertia"] == [11, 1, 150]
        and fold["pair_gate"]["pair_count"] == 13_041
        and fold["pair_gate"]["coxeter_edge_counts"]["forbidden"] == 915
        and not fold["pair_gate"]["all_pairs_coxeter"]
    ):
        raise ArithmeticError("the X_1 fixed-section screen changed")

    checks = {
        "complete_1094_wall_source_replayed": len(walls) == 1094,
        "all_720_genuine_stabilizer_matrices_replayed": len(group) == 720,
        "proper_source_group_high_subgroups_exhausted": denominator[
            "only_C2_can_be_a_proper_fixed_rank_at_least_11_subgroup"
        ],
        "all_162_negative_restrictions_have_strict_exact_facet_witnesses": len(
            fold["strict_facet_witnesses"]
        )
        == 162,
        "full_pair_angle_gate_exhausted": fold["pair_gate"]["pair_count"] == 13_041,
        "unique_high_section_is_non_Coxeter": not fold["pair_gate"]["all_pairs_coxeter"],
        "no_volume_inference_made": True,
    }
    if not all(checks.values()):
        raise ArithmeticError("an X_1 no-go check failed")
    return {
        "schema": SCHEMA,
        "id": "SHIMADA-AUT3SINGK3-X1-SYMP-H11-NO-GO",
        "source_stock": str(STOCK.relative_to(ROOT)),
        "source_stock_sha256": STOCK_SHA256,
        "denominator": {
            "included": "all subgroups of the genuine Aut(X_1,a_1)=PGL_2(F_9) stabilizer whose fixed rank is at least 11",
            "excluded": "the full abstract oriented weighted-wall graph group, if larger",
            "reason_for_scope": "the source theorem supplies exactly the 720 genuine matrices; the dense full weighted-graph computation was deferred under the Record-31 memory priority",
        },
        "high_subgroup_census": denominator,
        "unique_high_fixed_section": _public_fold(fold),
        "conclusion": {
            "Coxeter_polytope_above_H9_in_source_group_denominator": False,
            "reason": "the unique H11 section has 915 forbidden pairs among 162 independently strict facets",
            "compact_example_excluded": False,
            "finite_volume_noncompact_example_excluded_outside_this_denominator": False,
            "volume_status_of_folded_cone": "not tested and irrelevant after the Coxeter-angle failure",
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def validate_certificate(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("all_checks_pass") is not True:
        raise ArithmeticError("wrong or failed X_1 certificate")
    if payload["unique_high_fixed_section"]["pair_gate"]["coxeter_edge_counts"]["forbidden"] != 915:
        raise ArithmeticError("the archived forbidden-pair count changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise ArithmeticError(f"archive mismatch: {OUTPUT}")
        print(f"archive matches: {OUTPUT}")
        return
    OUTPUT.write_text(encoded, encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "fixed_dimension": 11,
        "strict_facets": 162,
        "forbidden_pairs": 915,
        "coxeter_survivor": False,
    }, indent=2))


if __name__ == "__main__":
    main()
