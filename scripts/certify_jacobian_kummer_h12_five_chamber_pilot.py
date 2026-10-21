#!/usr/bin/env python3
"""Certify five explicit five-chamber JK-H12 path gluings.

This is a pilot, not a census of five-vertex colored trees.  Starting from
the exact 172-wall rank-13 seed, it chooses one properly colored P5 for each
number B=0,...,4 of edges colored by W77--W172.  Every five-copy boundary is
rebuilt from reflection matrices; no facet formula is assumed.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base


BASE_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_small_subtree_census.py"
BASE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_small_subtree_census.json"
FOUR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py"
FOUR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_four_vertex_subtree_census.json"
PATH_CUBE = ROOT / "data/enumerations/jacobian_kummer_h12_path_cube_hostile_audit.json"
LINEAGE = ROOT / "data/enumerations/positive_example_lineage_header_controls.json"
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_five_chamber_pilot.json"

SCHEMA = "coxeter.jacobian_kummer.h12_five_chamber_pilot.v1"
STATUS_DATE = "2026-08-08"
RANK = 13
SEED_FACETS = 172
STOCK_START = 60
STOCK_END = 172
LOW_END = 76
TARGET_FACETS = (752, 761, 770, 779, 788)
LABEL_CODES = {"m2": 2, "m3": 3, "m4": 4, "m6": 6, "infinity": 0}

PINNED_SHA256 = {
    "small_tree_builder": "30ccb14a9240c32dfd8baeb70f1a3482a3941f5a6e4a9523094f4ebc4630d125",
    "small_tree_archive": "796575ebbaaf8ec0ba3054911397d83834d8939e1eb4c20cd24aae1b56aee56d",
    "four_tree_builder": "474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c",
    "four_tree_archive": "2591d5f078421b5b4902045d37fa57e33548e2450f090f5ea8513c564c1a4b24",
    "path_cube_hostile_audit": "5c81e4c311e96a64a86ec377d10aa40b54dcd94f685432fee11d9bc9f8ea5663",
    "project_lineage_snapshot": "58e00b546f50a9f3a3bcd4a21e16414aad6bb0d65f471475b6bfc199a57421aa",
}

# Absolute zero-based seed-wall indices.  Thus 60 is W61 and 76 is W77.
# All four colors in every path are distinct, a stronger condition than
# proper edge coloring for a path.
PILOTS: tuple[tuple[int, tuple[int, int, int, int]], ...] = (
    (0, (60, 61, 62, 63)),
    (1, (60, 61, 62, 76)),
    (2, (60, 61, 76, 77)),
    (3, (60, 76, 77, 78)),
    (4, (76, 77, 78, 79)),
)

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(payload: Any) -> str:
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(rendered.encode()).hexdigest()


def _content(vector: Iterable[int]) -> int:
    values = tuple(abs(int(value)) for value in vector if value)
    return math.gcd(*values) if values else 0


def _primitive_key_and_sign(vector: Sequence[int]) -> tuple[Vector, int]:
    entries = tuple(int(value) for value in vector)
    content = _content(entries)
    if not content:
        raise ArithmeticError("zero vector has no primitive hyperplane")
    primitive = tuple(value // content for value in entries)
    first = next(value for value in primitive if value)
    if first > 0:
        return primitive, 1
    return tuple(-value for value in primitive), -1


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((index for index in range(work.rows) if work[index, index]), None)
        if pivot is None:
            pair = next(
                (
                    (left, right)
                    for left in range(work.rows)
                    for right in range(left + 1, work.rows)
                    if work[left, right]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            order = list(pair) + [
                index for index in range(work.rows) if index not in pair
            ]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [index for index in range(work.rows) if index != pivot]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        positive += 1 if diagonal > 0 else 0
        negative += 1 if diagonal < 0 else 0
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = sp.simplify(work[1:, 1:] - column * column.T / diagonal)
    return positive, negative, 0


def _pinned_inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    paths = {
        "small_tree_builder": BASE_SCRIPT,
        "small_tree_archive": BASE_ARCHIVE,
        "four_tree_builder": FOUR_SCRIPT,
        "four_tree_archive": FOUR_ARCHIVE,
        "path_cube_hostile_audit": PATH_CUBE,
        "project_lineage_snapshot": LINEAGE,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a pinned five-chamber prerequisite changed: {hashes}")

    small = json.loads(BASE_ARCHIVE.read_text(encoding="utf-8"))
    four = json.loads(FOUR_ARCHIVE.read_text(encoding="utf-8"))
    path_cube = json.loads(PATH_CUBE.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not (
        small["all_checks_pass"]
        and four["all_checks_pass"]
        and path_cube["all_checks_pass"]
        and lineage["all_checks_pass"]
    ):
        raise ArithmeticError("a pinned project predecessor is not green")
    return small, four, {"path_cube": path_cube, "lineage": lineage}


def _stock_gate(gram: Matrix, roots: tuple[Vector, ...]) -> dict[str, Any]:
    squares = tuple(base._pair(gram, root, root) for root in roots)
    stock = tuple(range(STOCK_START, STOCK_END))
    doubling_count = 0
    for left in stock:
        labels = []
        for right in range(SEED_FACETS):
            if right == left:
                continue
            labels.append(
                base._relation_label(
                    base._pair(gram, roots[left], roots[right]),
                    squares[left],
                    squares[right],
                )
            )
        if "m3" in labels:
            raise ArithmeticError(f"W{left + 1} is not a doubling wall")
        doubling_count += 1
    disjoint_pairs = 0
    for offset, left in enumerate(stock):
        for right in stock[offset + 1 :]:
            label = base._relation_label(
                base._pair(gram, roots[left], roots[right]),
                squares[left],
                squares[right],
            )
            if not (label == "parallel" or label.startswith("ultraparallel:")):
                raise ArithmeticError("two stock walls meet in hyperbolic space")
            disjoint_pairs += 1
    return {
        "wall_range_one_based": "W61--W172",
        "wall_count": len(stock),
        "A_wall_count_W61_through_W76": LOW_END - STOCK_START,
        "B_wall_count_W77_through_W172": STOCK_END - LOW_END,
        "doubling_wall_count": doubling_count,
        "pairwise_disjoint_pair_count": disjoint_pairs,
        "all_doubling": doubling_count == 112,
        "all_pairwise_disjoint_in_Allcock_sense": disjoint_pairs == math.comb(112, 2),
    }


def _path_transformations(
    gram: Matrix, roots: tuple[Vector, ...], colors: tuple[int, int, int, int]
) -> tuple[tuple[Matrix, ...], dict[int, Matrix]]:
    reflections = {
        color: base._reflection_matrix(gram, roots[color]) for color in colors
    }
    gram_matrix = sp.Matrix(gram)
    for reflection in reflections.values():
        matrix = sp.Matrix(reflection)
        if matrix.T * gram_matrix * matrix != gram_matrix:
            raise ArithmeticError("a selected reflection is not a lattice isometry")
    transformations = [base._identity(RANK)]
    for color in colors:
        transformations.append(base._matmul(transformations[-1], reflections[color]))
    if len(set(transformations)) != 5:
        raise ArithmeticError("the path word does not give five distinct chambers")
    for transformation in transformations:
        matrix = sp.Matrix(transformation)
        if matrix.T * gram_matrix * matrix != gram_matrix:
            raise ArithmeticError("a path transformation is not a lattice isometry")
    return tuple(transformations), reflections


def _five_path_boundary(
    gram: Matrix,
    roots: tuple[Vector, ...],
    seed_interior: Vector,
    colors: tuple[int, int, int, int],
    core_smith: list[int],
) -> dict[str, Any]:
    if not (
        all(STOCK_START <= color < STOCK_END for color in colors)
        and len(set(colors)) == 4
    ):
        raise ArithmeticError("the pilot is not a properly colored P5")
    transformations, reflections = _path_transformations(gram, roots, colors)
    transform_index = {transformation: index for index, transformation in enumerate(transformations)}

    adjacency: set[tuple[int, int, int]] = set()
    internal: set[Vector] = set()
    for chamber, transformation in enumerate(transformations):
        for color in colors:
            neighbor = base._matmul(transformation, reflections[color])
            if neighbor not in transform_index:
                continue
            other = transform_index[neighbor]
            if chamber < other:
                adjacency.add((chamber, other, color))
                left_key = base._primitive_key(
                    base._matvec(transformation, roots[color])
                )
                right_key = base._primitive_key(base._matvec(neighbor, roots[color]))
                if left_key != right_key:
                    raise ArithmeticError("the two sides of an internal mirror disagree")
                internal.add(left_key)
    expected_adjacency = {
        (index, index + 1, color) for index, color in enumerate(colors)
    }
    if adjacency != expected_adjacency or len(internal) != 4:
        raise ArithmeticError("the selected transformations do not induce the marked P5")

    union_interior = tuple(
        sum(
            base._matvec(transformation, seed_interior)[coordinate]
            for transformation in transformations
        )
        for coordinate in range(RANK)
    )
    occurrences: defaultdict[Vector, set[int]] = defaultdict(set)
    for transformation in transformations:
        for seed_root in roots:
            image = base._matvec(transformation, seed_root)
            key, sign = _primitive_key_and_sign(image)
            occurrences[key].add(sign)
    canceled = {key for key, signs in occurrences.items() if signs == {-1, 1}}
    if canceled != internal:
        raise ArithmeticError("orientation cancellation and P5 adjacency disagree")
    if any(len(signs) != 1 for key, signs in occurrences.items() if key not in internal):
        raise ArithmeticError("an exterior mirror has inconsistent local orientations")

    boundary: list[Vector] = []
    for key in sorted(set(occurrences) - internal):
        pairing = base._pair(gram, union_interior, key)
        if pairing == 0:
            raise ArithmeticError("the summed five-chamber interior lies on a wall")
        boundary.append(key if pairing > 0 else tuple(-value for value in key))
    boundary_roots = tuple(boundary)
    boundary_keys = {base._primitive_key(root) for root in boundary_roots}
    core_keys = {base._primitive_key(root) for root in roots[:60]}
    if not core_keys <= boundary_keys or core_smith != [1] * RANK:
        raise ArithmeticError("the fixed Smith-1^13 Picard core left the boundary")

    covectors = tuple(base._covector(gram, root) for root in boundary_roots)
    squares = tuple(
        sum(root[index] * covector[index] for index in range(RANK))
        for root, covector in zip(boundary_roots, covectors, strict=True)
    )
    primitive_integral = 0
    for root, square, covector in zip(
        boundary_roots, squares, covectors, strict=True
    ):
        if not (
            _content(root) == 1
            and square < 0
            and all((2 * value) % square == 0 for value in covector)
        ):
            raise ArithmeticError("a boundary normal is not primitive integral reflective")
        primitive_integral += 1
    interior_pairings = tuple(
        base._pair(gram, union_interior, root) for root in boundary_roots
    )
    interior_square = base._pair(gram, union_interior, union_interior)
    if not (interior_square > 0 and all(pairing > 0 for pairing in interior_pairings)):
        raise ArithmeticError("the five-chamber boundary lacks a strict timelike interior")

    ordinary_histogram: Counter[str] = Counter()
    ordinary_word = bytearray()
    for left in range(len(boundary_roots)):
        for right in range(left):
            label = base._relation_label(
                sum(
                    boundary_roots[left][index] * covectors[right][index]
                    for index in range(RANK)
                ),
                squares[left],
                squares[right],
            )
            ordinary = (
                "infinity"
                if label == "parallel" or label.startswith("ultraparallel:")
                else label
            )
            if ordinary not in LABEL_CODES:
                raise ArithmeticError("a boundary pair is not an ordinary Coxeter relation")
            ordinary_histogram[ordinary] += 1
            ordinary_word.append(LABEL_CODES[ordinary])
    facet_count = len(boundary_roots)
    if sum(ordinary_histogram.values()) != facet_count * (facet_count - 1) // 2:
        raise ArithmeticError("the ordinary Coxeter matrix is incomplete")
    return {
        "tree_shape": "P5",
        "edge_colors_one_based": [color + 1 for color in colors],
        "proper_edge_coloring": True,
        "all_four_edge_colors_distinct": True,
        "chamber_count": len(transformations),
        "marked_chamber_adjacency": [
            [left + 1, right + 1, color + 1]
            for left, right, color in sorted(adjacency)
        ],
        "internal_mirror_count": len(internal),
        "orientation_cancellation_internal_mirror_count": len(canceled),
        "facet_count": facet_count,
        "strict_interior_square": interior_square,
        "minimum_strict_interior_facet_pairing": min(interior_pairings),
        "primitive_integral_reflection_root_count": primitive_integral,
        "full_linear_rank": RANK,
        "fixed_boundary_core_one_based": "W1--W60",
        "fixed_boundary_core_Smith_diagonal": core_smith,
        "marked_boundary_root_sha256": base._canonical_sha256(
            [list(root) for root in boundary_roots]
        ),
        "marked_ordinary_Coxeter_matrix_sha256": hashlib.sha256(
            len(boundary_roots).to_bytes(2, "big") + bytes(ordinary_word)
        ).hexdigest(),
        "ordinary_Coxeter_relation_histogram": {
            label: ordinary_histogram[label]
            for label in ("m2", "m3", "m4", "m6", "infinity")
        },
        "maximum_boundary_coordinate_absolute": max(
            abs(value) for root in boundary_roots for value in root
        ),
    }


def build() -> dict[str, Any]:
    small, four, project = _pinned_inputs()
    carrier, gram, roots, interior, bounded_corpus = base._load_inputs()
    if not (
        len(gram) == RANK
        and len(roots) == SEED_FACETS
        and _exact_inertia(sp.Matrix(gram)) == (1, 12, 0)
        and int(sp.Matrix(gram).det()) == 256
    ):
        raise ArithmeticError("the raw rank-13 hyperbolic seed changed")
    stock = _stock_gate(gram, roots)
    if not (stock["all_doubling"] and stock["all_pairwise_disjoint_in_Allcock_sense"]):
        raise ArithmeticError("the Allcock stock gate failed")

    core_matrix = sp.Matrix(roots[:60])
    smith = smith_normal_form(core_matrix, domain=ZZ)
    core_smith = [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]
    if core_smith != [1] * RANK:
        raise ArithmeticError("the seed Picard core is not primitive of index one")

    rows = []
    for b_count, colors in PILOTS:
        if sum(color >= LOW_END for color in colors) != b_count:
            raise ArithmeticError("a pilot has the wrong B-edge count")
        row = _five_path_boundary(gram, roots, interior, colors, core_smith)
        expected = 172 + 145 * (row["chamber_count"] - 1) + 9 * b_count
        if row["facet_count"] != expected:
            raise ArithmeticError(
                f"five-chamber formula correction at B={b_count}: "
                f"observed {row['facet_count']}, predicted {expected}"
            )
        rows.append(
            {
                "pilot_number": len(rows) + 1,
                "B_edge_occurrence_count": b_count,
                "expected_facet_count_172_plus_145_times_4_plus_9B": expected,
                **row,
            }
        )
    facet_counts = [row["facet_count"] for row in rows]
    if facet_counts != list(TARGET_FACETS) or len(set(facet_counts)) != 5:
        raise ArithmeticError("the five pilot facet counts changed")

    prior_project_counts = {
        172,
        *(
            int(row["facet_count"])
            for row in small["complete_one_edge_census"]["type_representatives"]
        ),
        *(
            int(value)
            for value in small["complete_three_vertex_coloured_path_census"][
                "diagram_type_facet_histogram"
            ]
        ),
        *(int(value) for value in four["theorem"]["type_facet_histogram"]),
    }
    lineage_rows = project["lineage"]["ledger"]["rows"]
    snapshot_h12_counts = {
        int(row["facets"]) for row in lineage_rows if int(row["dimension"]) == 12
    }
    corpus_counts = {
        int(value) for value in bounded_corpus["ordinary_H12_facet_histogram"]
    } | {
        int(value) for value in bounded_corpus["index_two_H12_facet_histogram"]
    }
    targets = set(TARGET_FACETS)
    if not (
        targets.isdisjoint(prior_project_counts)
        and targets.isdisjoint(snapshot_h12_counts)
        and targets.isdisjoint(corpus_counts)
        and project["path_cube"]["universal_family"]["facet_formula"]
        == "f(N)=886*N+55 for every integer N>=1"
    ):
        raise ArithmeticError("a five-chamber facet header failed bounded lineage")

    checks = {
        "raw_172_wall_rank13_seed_loaded_without_four_tree_boundary_data": True,
        "seed_lattice_has_signature_1_12_and_determinant_256": True,
        "all_112_stock_walls_are_doubling_and_pairwise_disjoint": True,
        "one_properly_colored_P5_constructed_for_each_B_count_0_through_4": len(rows) == 5,
        "all_five_marked_chamber_sets_are_embedded_P5_trees": all(
            row["chamber_count"] == 5
            and row["internal_mirror_count"] == 4
            and row["orientation_cancellation_internal_mirror_count"] == 4
            for row in rows
        ),
        "direct_facet_counts_are_752_761_770_779_788": facet_counts == list(TARGET_FACETS),
        "pilot_formula_172_plus_145_times_V_minus_1_plus_9B_is_exact": True,
        "all_boundary_pairs_are_ordinary_Coxeter_relations": True,
        "all_boundary_normals_are_primitive_integral_reflection_roots": all(
            row["primitive_integral_reflection_root_count"] == row["facet_count"]
            for row in rows
        ),
        "all_boundaries_have_strict_timelike_interiors": all(
            row["strict_interior_square"] > 0
            and row["minimum_strict_interior_facet_pairing"] > 0
            for row in rows
        ),
        "all_boundaries_contain_fixed_Smith_1_power_13_core": all(
            row["fixed_boundary_core_Smith_diagonal"] == [1] * RANK for row in rows
        ),
        "five_ordinary_Coxeter_matrices_are_pairwise_nonisomorphic_by_facet_count": len(
            set(facet_counts)
        )
        == 5,
        "Allcock_index_and_volume_multiplier_are_five": True,
        "all_five_pilots_are_finite_volume_noncompact": True,
        "target_headers_absent_from_pinned_project_predecessors_and_CoxIter_lanes": True,
        "worldwide_priority_explicitly_open": True,
        "claim_is_pilot_not_full_five_vertex_census": True,
        "carrier_claim_is_weak_not_full_nef_or_surface_action": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"five-chamber pilot failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": {
            **PINNED_SHA256,
            "published_Dolgachev_Shimada_data": base.PINNED_SHA256[
                "published_Dolgachev_Shimada_data"
            ],
            "H12_fold_orbits": base.PINNED_SHA256["H12_fold_orbits"],
            "full_Picard_K3_carrier": base.PINNED_SHA256[
                "full_Picard_K3_carrier"
            ],
            "bounded_CoxIter_corpus_audit": base.PINNED_SHA256[
                "bounded_CoxIter_corpus_audit"
            ],
        },
        "result": {
            "decision": "PASS_FREEZE",
            "hyperbolic_dimension": 12,
            "explicit_five_chamber_polytope_count": 5,
            "facet_counts": facet_counts,
            "pairwise_nonisomorphic_ordinary_Coxeter_matrix_count": 5,
            "reflection_index_and_volume_multiplier": 5,
            "finite_volume": True,
            "compact": False,
            "worldwide_priority": "open",
        },
        "seed_and_Allcock_stock": {
            "seed_id": "JK-H12/172",
            "seed_facet_count": SEED_FACETS,
            "seed_lattice": "U(2)+D5(-1)^2+<-4>",
            "seed_lattice_signature": [1, 12, 0],
            "seed_lattice_determinant": 256,
            "seed_finite_volume_noncompact": True,
            "stock": stock,
            "fixed_core_Smith_diagonal": core_smith,
        },
        "pilot_theorem": {
            "scope": "one explicit P5 for each B=0,...,4; not a five-vertex census",
            "verified_formula": "f(P5)=172+145*(5-1)+9*B=752+9*B on the five listed pilots",
            "B_to_facet_count": {str(index): count for index, count in enumerate(facet_counts)},
            "ordinary_nonisomorphism": "the five ordinary Coxeter matrices have distinct orders 752,761,770,779,788",
            "Allcock_conclusion": "each P5 is a five-copy fundamental polytope of reflection index and volume multiplier five",
            "finite_volume_noncompact_reason": "finite index in the nonuniform finite-covolume seed reflection group",
        },
        "pilots": rows,
        "lineage": {
            "pinned_project_predecessor_facet_counts": sorted(prior_project_counts),
            "pinned_project_ledger_H12_facet_counts": sorted(snapshot_h12_counts),
            "target_counts_absent_from_both_project_sets": True,
            "CoxIter_commit": bounded_corpus["commit"],
            "ordinary_H12_file_count": bounded_corpus["ordinary_file_count"],
            "index_two_file_count": bounded_corpus["index_two_file_count"],
            "pinned_CoxIter_H12_facet_counts": sorted(corpus_counts),
            "target_headers_absent_from_pinned_CoxIter_lanes": True,
            "evidence_grade": "B1: bounded project/header absence only",
            "worldwide_priority": "open",
        },
        "carrier_scope": {
            "same_full_Picard_K3_family": True,
            "K3_moduli_dimension": 7,
            "fixed_index_one_divisor_lattice_core": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_realization_claimed": False,
            "surface_automorphism_realization_claimed": False,
        },
        "claim_boundary": {
            "full_five_vertex_coloring_or_shape_census_claimed": False,
            "general_facet_formula_for_every_five_vertex_tree_claimed": False,
            "five_explicit_directly_certified_examples_claimed": True,
            "worldwide_novelty_claimed": False,
        },
        "exact_arithmetic": {
            "backend": "Python arbitrary-precision integers and exact SymPy rationals",
            "maximum_boundary_coordinate_absolute": max(
                row["maximum_boundary_coordinate_absolute"] for row in rows
            ),
            "all_products_exact_without_fixed_width_overflow": True,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("five-chamber pilot archive differs from exact replay")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("facet_counts=752,761,770,779,788")
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
