#!/usr/bin/env python3
"""Certify the complete five-vertex K1,4 JK-H12 subtree census.

The edge colours are four-element subsets of the frozen 112-wall Allcock
stock W61--W172.  We quotient the complete set by the exact 384-element seed
action, rebuild every five-copy boundary over the rank-13 Picard lattice, and
compare the resulting ordinary Coxeter matrices by exact invariants.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any, Iterable, Sequence

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base
import scripts.certify_jacobian_kummer_h12_four_vertex_subtree_census as four


BASE_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_small_subtree_census.py"
BASE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_small_subtree_census.json"
FOUR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py"
FOUR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_four_vertex_subtree_census.json"
PILOT_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_chamber_pilot.py"
PILOT_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_chamber_pilot.json"
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"

SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_star_census.v1"
STATUS_DATE = "2026-08-08"
RANK = 13
STOCK_START = 60
STOCK_END = 172
STOCK_SIZE = 112
LOW_END = 76
RAW_COLORING_COUNT = math.comb(STOCK_SIZE, 4)
INT64_MAX = np.iinfo(np.int64).max
LABEL_NAMES = ("m2", "m3", "m4", "m6", "infinity")

PINNED_SHA256 = {
    "small_tree_builder": "30ccb14a9240c32dfd8baeb70f1a3482a3941f5a6e4a9523094f4ebc4630d125",
    "small_tree_archive": "796575ebbaaf8ec0ba3054911397d83834d8939e1eb4c20cd24aae1b56aee56d",
    "four_tree_builder": "474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c",
    "four_tree_archive": "2591d5f078421b5b4902045d37fa57e33548e2450f090f5ea8513c564c1a4b24",
    "five_chamber_pilot_builder": "d1e678dc2eea7c67998fff1d5fcef34a8ad45281fac87424f74f3876ea4c9f6d",
    "five_chamber_pilot_archive": "809302f4aa10fcf723c68061c34a4d857c9d5c64a4c03e04227c4c54c917ebec",
}

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
Permutation = tuple[int, ...]
ColorFour = tuple[int, int, int, int]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(rendered.encode()).hexdigest()


def _content(vector: Iterable[int]) -> int:
    values = tuple(abs(int(value)) for value in vector if value)
    return math.gcd(*values) if values else 0


def _primitive_key(vector: Sequence[int]) -> Vector:
    content = _content(vector)
    if not content:
        raise ArithmeticError("zero vector has no primitive hyperplane")
    primitive = tuple(int(value) // content for value in vector)
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _combination_rank(colors: Sequence[int]) -> int:
    """Colexicographic rank of a four-subset of range(112)."""

    if len(colors) != 4 or tuple(colors) != tuple(sorted(set(colors))):
        raise ArithmeticError("invalid four-subset")
    return sum(math.comb(int(value), index + 1) for index, value in enumerate(colors))


def _restricted_actions() -> tuple[Permutation, ...]:
    actions = base._folded_chamber_actions()
    restricted = tuple(
        tuple(action[index] - STOCK_START for index in range(STOCK_START, STOCK_END))
        for action in actions
    )
    if not (
        len(restricted) == 384
        and len(set(restricted)) == 384
        and all(sorted(action) == list(range(STOCK_SIZE)) for action in restricted)
    ):
        raise ArithmeticError("the restricted 112-wall action changed")
    return restricted


def _cycle_lengths(permutation: Permutation) -> tuple[int, ...]:
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        if current != start:
            raise ArithmeticError("a restricted action is not a permutation")
        lengths.append(length)
    return tuple(sorted(lengths))


def _fixed_four_subsets(permutation: Permutation) -> int:
    coefficients = [1, 0, 0, 0, 0]
    for length in _cycle_lengths(permutation):
        for degree in range(4, length - 1, -1):
            coefficients[degree] += coefficients[degree - length]
    return coefficients[4]


def _star_color_orbits(actions: Sequence[Permutation]) -> tuple[dict[str, Any], ...]:
    seen = bytearray(RAW_COLORING_COUNT)
    rows: list[dict[str, Any]] = []
    for colors in combinations(range(STOCK_SIZE), 4):
        rank = _combination_rank(colors)
        if seen[rank]:
            continue
        orbit = {
            tuple(sorted(action[color] for color in colors)) for action in actions
        }
        for image in orbit:
            seen[_combination_rank(image)] = 1
        rows.append(
            {
                "colors_zero_based_stock": colors,
                "unrooted_coloring_orbit_size_under_seed_action": len(orbit),
            }
        )
    burnside_sum = sum(_fixed_four_subsets(action) for action in actions)
    if not (
        len(seen) == RAW_COLORING_COUNT
        and sum(seen) == RAW_COLORING_COUNT
        and sum(row["unrooted_coloring_orbit_size_under_seed_action"] for row in rows)
        == RAW_COLORING_COUNT
        and burnside_sum % len(actions) == 0
        and burnside_sum // len(actions) == len(rows)
    ):
        raise ArithmeticError("the complete K1,4 color-orbit census failed")
    return tuple(rows)


class ExactStarBoundaryEngine:
    def __init__(self, gram: Matrix, roots: tuple[Vector, ...], interior: Vector) -> None:
        self.gram_tuple = gram
        self.gram = np.array(gram, dtype=np.int64)
        self.roots_tuple = roots
        self.roots = np.array(roots, dtype=np.int64)
        self.interior = np.array(interior, dtype=np.int64)
        self.identity = np.eye(RANK, dtype=np.int64)
        self.reflections = {
            index: np.array(base._reflection_matrix(gram, roots[index]), dtype=np.int64)
            for index in range(STOCK_START, STOCK_END)
        }
        gram_maximum = int(np.max(np.abs(self.gram)))
        reflection_maximum = max(
            int(np.max(np.abs(reflection))) for reflection in self.reflections.values()
        )
        root_maximum = int(np.max(np.abs(self.roots)))
        interior_maximum = int(np.max(np.abs(self.interior)))
        image_bound = RANK * reflection_maximum * root_maximum
        transformed_interior_bound = RANK * reflection_maximum * interior_maximum
        summed_interior_bound = interior_maximum + 4 * transformed_interior_bound
        self.int64_preflight_bounds = {
            "reflection_entry_absolute": reflection_maximum,
            "two_reflection_product_entry_absolute": RANK * reflection_maximum**2,
            "transformed_seed_root_entry_absolute": image_bound,
            "sum_of_five_transformed_interiors_entry_absolute": summed_interior_bound,
            "orientation_pairing_absolute": (
                RANK**2 * gram_maximum * image_bound * summed_interior_bound
            ),
            "boundary_pairing_absolute": RANK**2 * gram_maximum * image_bound**2,
            "summed_interior_square_absolute": (
                RANK**2 * gram_maximum * summed_interior_bound**2
            ),
        }
        if max(self.int64_preflight_bounds.values()) >= INT64_MAX:
            raise ArithmeticError("the five-chamber star int64 preflight failed")
        core_smith = smith_normal_form(sp.Matrix(roots[:60]), domain=ZZ)
        self.core_smith = [
            abs(int(core_smith[index, index]))
            for index in range(min(core_smith.rows, core_smith.cols))
            if core_smith[index, index]
        ]
        self.core_keys = {_primitive_key(root) for root in roots[:60]}
        if self.core_smith != [1] * RANK:
            raise ArithmeticError("the fixed W1--W60 Smith core changed")

    def boundary(self, stock_colors: ColorFour) -> tuple[np.ndarray, dict[str, Any]]:
        colors = tuple(color + STOCK_START for color in stock_colors)
        transformations = (self.identity,) + tuple(
            self.reflections[color] for color in colors
        )
        transform_keys = {four._matrix_key(matrix) for matrix in transformations}
        if len(transform_keys) != 5:
            raise ArithmeticError("a star has fewer than five chambers")

        internal: set[Vector] = set()
        adjacency: set[tuple[int, int, int]] = set()
        for chamber, transformation in enumerate(transformations):
            for color in colors:
                neighbor_key = four._matrix_key(transformation @ self.reflections[color])
                if neighbor_key not in transform_keys:
                    continue
                other = next(
                    index
                    for index, candidate in enumerate(transformations)
                    if four._matrix_key(candidate) == neighbor_key
                )
                if chamber < other:
                    adjacency.add((chamber, other, color))
                internal.add(_primitive_key(transformation @ self.roots[color]))
        expected_adjacency = {(0, index + 1, color) for index, color in enumerate(colors)}
        if adjacency != expected_adjacency or len(internal) != 4:
            raise ArithmeticError("the transformations do not induce the marked K1,4")

        union_interior = sum(
            (transformation @ self.interior for transformation in transformations),
            start=np.zeros(RANK, dtype=np.int64),
        )
        boundary: dict[Vector, Vector] = {}
        for transformation in transformations:
            images = self.roots @ transformation.T
            pairings = images @ self.gram @ union_interior
            for image, pairing in zip(images, pairings, strict=True):
                key = _primitive_key(image)
                if key in internal:
                    continue
                if pairing == 0:
                    raise ArithmeticError("the summed star interior lies on a boundary")
                oriented = image if pairing > 0 else -image
                boundary.setdefault(key, tuple(int(value) for value in oriented))

        ordered_keys = tuple(sorted(boundary))
        if not self.core_keys <= set(ordered_keys):
            raise ArithmeticError("the fixed Picard core left a five-chamber boundary")
        roots = np.array([boundary[key] for key in ordered_keys], dtype=np.int64)
        covectors = roots @ self.gram
        squares = np.einsum("ij,ij->i", roots, covectors)
        if not np.all(squares < 0):
            raise ArithmeticError("a boundary root is not spacelike")
        primitive_integral = 0
        for root, square, covector in zip(roots, squares, covectors, strict=True):
            if not (
                _content(root) == 1
                and all((2 * int(value)) % int(square) == 0 for value in covector)
            ):
                raise ArithmeticError("a boundary reflection is nonintegral")
            primitive_integral += 1
        interior_pairings = roots @ self.gram @ union_interior
        interior_square = int(union_interior @ self.gram @ union_interior)
        if not (interior_square > 0 and np.all(interior_pairings > 0)):
            raise ArithmeticError("a five-chamber star strict interior failed")
        return roots, {
            "facet_count": len(roots),
            "chamber_count": 5,
            "internal_mirror_count": len(internal),
            "marked_chamber_adjacency": [
                [left + 1, right + 1, color + 1]
                for left, right, color in sorted(adjacency)
            ],
            "strict_interior_square": interior_square,
            "minimum_strict_interior_facet_pairing": int(np.min(interior_pairings)),
            "primitive_integral_reflection_root_count": primitive_integral,
            "fixed_boundary_core_Smith_diagonal": self.core_smith,
            "full_linear_rank": RANK,
            "marked_boundary_root_sha256": _json_sha256(roots.tolist()),
            "maximum_boundary_coordinate_absolute": int(np.max(np.abs(roots))),
        }

    def ordinary_data(self, roots: np.ndarray) -> tuple[np.ndarray, bytes, dict[str, Any]]:
        return four.ExactBoundaryEngine.ordinary_coxeter_data(self, roots)


def _pinned_predecessors() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    paths = {
        "small_tree_builder": BASE_SCRIPT,
        "small_tree_archive": BASE_ARCHIVE,
        "four_tree_builder": FOUR_SCRIPT,
        "four_tree_archive": FOUR_ARCHIVE,
        "five_chamber_pilot_builder": PILOT_SCRIPT,
        "five_chamber_pilot_archive": PILOT_ARCHIVE,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a pinned five-vertex prerequisite changed: {hashes}")
    small = json.loads(BASE_ARCHIVE.read_text(encoding="utf-8"))
    four_archive = json.loads(FOUR_ARCHIVE.read_text(encoding="utf-8"))
    pilot = json.loads(PILOT_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        small["all_checks_pass"]
        and four_archive["all_checks_pass"]
        and pilot["all_checks_pass"]
        and small["selected_stock"]["wall_count"] == STOCK_SIZE
        and small["selected_stock"]["all_doubling"]
        and small["selected_stock"]["all_pairwise_disjoint_in_Allcock_sense"]
        and small["selected_stock"]["folded_chamber_action_order"] == 384
        and pilot["result"]["reflection_index_and_volume_multiplier"] == 5
        and pilot["result"]["finite_volume"]
        and not pilot["result"]["compact"]
    ):
        raise ArithmeticError("a pinned five-vertex predecessor is not green")
    return small, four_archive, pilot


def _root_gram_action_gate(
    gram: Matrix, roots: tuple[Vector, ...], actions: Sequence[Permutation]
) -> str:
    root_array = np.array(roots, dtype=np.int64)
    gram_array = np.array(gram, dtype=np.int64)
    root_gram = root_array @ gram_array @ root_array.T
    for restricted in actions:
        full = np.arange(172)
        full[STOCK_START:STOCK_END] = np.array(restricted) + STOCK_START
        # The restricted permutation alone does not specify W1--W60, so this
        # gate checks exactly the stock-stock block used by the coloring quotient.
        stock_image = full[STOCK_START:STOCK_END]
        if not np.array_equal(
            root_gram[STOCK_START:STOCK_END, STOCK_START:STOCK_END],
            root_gram[np.ix_(stock_image, stock_image)],
        ):
            raise ArithmeticError("a restricted seed action changes the stock Gram")
    return _json_sha256(root_gram[STOCK_START:STOCK_END, STOCK_START:STOCK_END].tolist())


def _public_color_row(row: dict[str, Any], orbit_number: int) -> dict[str, Any]:
    stock_colors = tuple(int(value) for value in row["colors_zero_based_stock"])
    absolute = tuple(value + STOCK_START for value in stock_colors)
    b_count = sum(value >= LOW_END for value in absolute)
    return {
        "orbit_number": orbit_number,
        "tree_shape": "K1,4",
        "edge_colors_one_based": [value + 1 for value in absolute],
        "edge_color_class_word": "A" * (4 - b_count) + "B" * b_count,
        "B_edge_occurrence_count": b_count,
        "unrooted_coloring_orbit_size_under_seed_action": row[
            "unrooted_coloring_orbit_size_under_seed_action"
        ],
    }


def build(*, verbose: bool = False) -> dict[str, Any]:
    small, four_archive, pilot = _pinned_predecessors()
    carrier, gram, roots, interior, bounded_corpus = base._load_inputs()
    actions = _restricted_actions()
    stock_gram_sha256 = _root_gram_action_gate(gram, roots, actions)
    color_rows = _star_color_orbits(actions)
    if len(color_rows) != 16_872:
        raise ArithmeticError("the complete K1,4 seed-orbit count changed")

    cycle_types = Counter(
        tuple(sorted(Counter(_cycle_lengths(action)).items())) for action in actions
    )
    burnside_sum = sum(_fixed_four_subsets(action) for action in actions)
    if burnside_sum != 6_478_848:
        raise ArithmeticError("the independent K1,4 Burnside sum changed")

    engine = ExactStarBoundaryEngine(gram, roots, interior)
    degree_groups: defaultdict[tuple[int, bytes], list[int]] = defaultdict(list)
    public_rows: list[dict[str, Any]] = []
    relation_histograms_by_facet: defaultdict[
        int, Counter[tuple[tuple[str, int], ...]]
    ] = defaultdict(Counter)
    coordinate_maximum = 0
    pairing_maximum = 0
    started = time.perf_counter()
    for index, color_row in enumerate(color_rows):
        stock_colors = tuple(int(value) for value in color_row["colors_zero_based_stock"])
        boundary, boundary_data = engine.boundary(stock_colors)
        labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        del labels
        public = _public_color_row(color_row, index + 1)
        b_count = public["B_edge_occurrence_count"]
        expected_facets = 172 + 145 * 4 + 9 * b_count
        if boundary_data["facet_count"] != expected_facets:
            raise ArithmeticError("the K1,4 facet formula changed")
        degree_groups[(len(boundary), degree_key)].append(index)
        relation_histograms_by_facet[len(boundary)][
            tuple(sorted(ordinary_data["ordinary_Coxeter_relation_histogram"].items()))
        ] += 1
        coordinate_maximum = max(
            coordinate_maximum, boundary_data["maximum_boundary_coordinate_absolute"]
        )
        pairing_maximum = max(
            pairing_maximum, ordinary_data["maximum_boundary_pairing_absolute"]
        )
        public_rows.append(
            {
                **public,
                **boundary_data,
                **ordinary_data,
                "proper_edge_coloring": True,
                "all_four_incident_edge_colors_distinct": True,
                "expected_facet_count_752_plus_9B": expected_facets,
                "fixed_boundary_core_one_based": "W1--W60",
                "_degree_key": degree_key,
            }
        )
        if verbose and ((index + 1) % 500 == 0 or index + 1 == len(color_rows)):
            elapsed = time.perf_counter() - started
            print(
                f"boundary_progress={index + 1}/{len(color_rows)} "
                f"elapsed_seconds={elapsed:.1f}",
                flush=True,
            )

    collision_groups = tuple(
        indices for indices in degree_groups.values() if len(indices) > 1
    )
    collision_indices = {index for group in collision_groups for index in group}
    refinement_keys: dict[int, bytes] = {}
    refinement_discrete = 0
    for replay_number, index in enumerate(sorted(collision_indices), start=1):
        color_row = color_rows[index]
        stock_colors = tuple(int(value) for value in color_row["colors_zero_based_stock"])
        boundary, _boundary_data = engine.boundary(stock_colors)
        labels, degree_key, _ordinary_data = engine.ordinary_data(boundary)
        if degree_key != public_rows[index]["_degree_key"]:
            raise ArithmeticError("a residual K1,4 degree profile changed on replay")
        refinement_key, refinement_data = four._canonical_color_refinement(labels)
        refinement_keys[index] = refinement_key
        refinement_discrete += int(
            refinement_data["color_refinement_discrete_canonical_form"]
        )
        public_rows[index].update(refinement_data)
        if verbose and replay_number % 100 == 0:
            print(
                f"refinement_progress={replay_number}/{len(collision_indices)}",
                flush=True,
            )

    unresolved_groups = []
    for group in collision_groups:
        buckets: defaultdict[bytes, list[int]] = defaultdict(list)
        for index in group:
            buckets[refinement_keys[index]].append(index)
        unresolved_groups.extend(bucket for bucket in buckets.values() if len(bucket) > 1)
    if unresolved_groups:
        sizes = Counter(len(group) for group in unresolved_groups)
        raise ArithmeticError(
            f"canonical color refinement left unresolved K1,4 groups: {sizes}"
        )

    exact_type_keys = []
    for index, row in enumerate(public_rows):
        degree_key = row.pop("_degree_key")
        refinement_key = refinement_keys.get(index, b"")
        exact_key = (
            row["facet_count"].to_bytes(2, "big")
            + degree_key
            + b"\0refinement\0"
            + refinement_key
        )
        row["ordinary_Coxeter_type_invariant_sha256"] = hashlib.sha256(
            exact_key
        ).hexdigest()
        row.setdefault("color_refinement_applied", False)
        exact_type_keys.append(exact_key)
    if len(set(exact_type_keys)) != len(public_rows):
        raise ArithmeticError("two K1,4 ordinary Coxeter invariants collided")

    type_facet_histogram = Counter(row["facet_count"] for row in public_rows)
    raw_facet_histogram = Counter()
    orbit_size_histogram = Counter()
    for row in public_rows:
        orbit_size = row["unrooted_coloring_orbit_size_under_seed_action"]
        raw_facet_histogram[row["facet_count"]] += orbit_size
        orbit_size_histogram[orbit_size] += 1
    expected_raw_facet_histogram = Counter(
        {
            752 + 9 * b_count: math.comb(16, 4 - b_count) * math.comb(96, b_count)
            for b_count in range(5)
        }
    )
    if not (
        raw_facet_histogram == expected_raw_facet_histogram
        and sum(raw_facet_histogram.values()) == RAW_COLORING_COUNT
        and set(type_facet_histogram) == {752, 761, 770, 779, 788}
    ):
        raise ArithmeticError("the K1,4 facet histograms changed")

    relation_inventory = {
        str(facets): {
            "distinct_relation_histogram_count": len(histograms),
            "histogram_multiplicity_sha256": _json_sha256(
                [
                    {"histogram": dict(histogram), "type_count": multiplicity}
                    for histogram, multiplicity in sorted(histograms.items())
                ]
            ),
        }
        for facets, histograms in sorted(relation_histograms_by_facet.items())
    }
    checks = {
        "frozen_112_wall_factory_and_predecessors_replayed": True,
        "exact_384_action_preserves_stock_Gram": True,
        "all_6210820_unordered_distinct_color_four_subsets_covered": True,
        "Burnside_sum_6478848_independently_gives_16872_orbits": True,
        "all_16872_orbit_representative_boundaries_rebuilt": len(public_rows) == 16_872,
        "every_chamber_set_is_an_embedded_K14_with_four_internal_mirrors": all(
            row["chamber_count"] == 5 and row["internal_mirror_count"] == 4
            for row in public_rows
        ),
        "every_boundary_has_a_strict_timelike_interior": all(
            row["strict_interior_square"] > 0
            and row["minimum_strict_interior_facet_pairing"] > 0
            for row in public_rows
        ),
        "every_boundary_root_is_primitive_integral_reflective": all(
            row["primitive_integral_reflection_root_count"] == row["facet_count"]
            for row in public_rows
        ),
        "every_boundary_contains_the_fixed_Smith_1_power_13_core": all(
            row["fixed_boundary_core_Smith_diagonal"] == [1] * RANK
            for row in public_rows
        ),
        "every_boundary_pair_is_an_ordinary_Coxeter_relation": True,
        "facet_formula_is_752_plus_9_times_B_edge_occurrences": all(
            row["facet_count"] == row["expected_facet_count_752_plus_9B"]
            for row in public_rows
        ),
        "exact_degree_profiles_plus_canonical_refinement_separate_every_type": True,
        "all_16872_ordinary_Coxeter_matrices_are_pairwise_nonisomorphic": len(
            set(exact_type_keys)
        )
        == 16_872,
        "Allcock_index_and_volume_multiplier_are_five": True,
        "all_descendants_are_finite_volume_noncompact": True,
        "same_full_Picard_K3_weak_carrier_retained": True,
        "worldwide_priority_explicitly_open": True,
        "claim_is_complete_K14_shape_not_all_five_vertex_shapes": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"five-vertex star census failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": {
            **PINNED_SHA256,
            "full_Picard_K3_carrier": base.PINNED_SHA256["full_Picard_K3_carrier"],
            "bounded_CoxIter_corpus_audit": base.PINNED_SHA256[
                "bounded_CoxIter_corpus_audit"
            ],
        },
        "seed": {
            "id": "JK-H12/172",
            "hyperbolic_dimension": 12,
            "facet_count": 172,
            "finite_volume": True,
            "compact": False,
            "selected_pairwise_disjoint_doubling_stock": "W61--W172",
            "stock_size": STOCK_SIZE,
            "folded_seed_action_order": len(actions),
            "stock_Gram_sha256": stock_gram_sha256,
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "source_carrier_theorem": carrier["theorem"],
        },
        "coloring_denominator": {
            "shape": "K1,4",
            "rule": "unordered four-subset of distinct incident edge colors",
            "raw_unrooted_proper_colorings": RAW_COLORING_COUNT,
            "seed_action_order": len(actions),
            "Burnside_fixed_subset_sum": burnside_sum,
            "Burnside_orbit_count": len(color_rows),
            "explicit_orbit_cover_count": len(color_rows),
            "explicit_orbit_size_sum": sum(
                row["unrooted_coloring_orbit_size_under_seed_action"]
                for row in color_rows
            ),
            "seed_orbit_size_histogram": {
                str(size): count for size, count in sorted(orbit_size_histogram.items())
            },
            "seed_action_cycle_type_histogram": {
                str(cycle_type): multiplicity
                for cycle_type, multiplicity in sorted(cycle_types.items())
            },
        },
        "theorem": {
            "scope": "complete K1,4 five-vertex shape only",
            "ordinary_Coxeter_type_count": len(public_rows),
            "pairwise_nonisomorphic": True,
            "facet_formula": "f(K1,4)=172+145*4+9*B=752+9*B",
            "type_facet_histogram": {
                str(facets): count for facets, count in sorted(type_facet_histogram.items())
            },
            "raw_unrooted_coloring_facet_histogram": {
                str(facets): count for facets, count in sorted(raw_facet_histogram.items())
            },
            "ordinary_relation_inventory_by_facet_count": relation_inventory,
            "reflection_index_and_volume_multiple": 5,
            "finite_volume": True,
            "compact": False,
        },
        "ordinary_Coxeter_nonisomorphism_certificate": {
            "parallel_and_ultraparallel_both_labelled_infinity": True,
            "initial_exact_degree_profile_class_count": len(degree_groups),
            "degree_profile_collision_group_count": len(collision_groups),
            "degree_profile_collision_member_count": len(collision_indices),
            "maximum_initial_collision_group_size": max(
                (len(group) for group in collision_groups), default=1
            ),
            "color_refined_member_count": len(refinement_keys),
            "color_refined_discrete_canonical_form_count": refinement_discrete,
            "all_residual_refinement_invariants_distinct_within_initial_class": True,
            "all_combined_invariants_distinct": True,
            "method": (
                "exact multiset of five ordinary-label degrees at every facet; "
                "canonical one-dimensional color refinement of the complete ordinary "
                "Coxeter matrix is applied only inside exact degree-profile collisions"
            ),
        },
        "Allcock_and_carrier": {
            "finite_subtree_theorem": (
                "every representative is a five-copy K1,4 fundamental polytope of "
                "reflection index and volume multiplier five"
            ),
            "noncompactness": (
                "finite-index subgroups of the nonuniform seed lattice are nonuniform"
            ),
            "same_full_Picard_K3_family": True,
            "fixed_index_one_divisor_lattice_core": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_effectivity_or_surface_action_claim": False,
        },
        "arithmetic_runtime_bounds": {
            "integer_backend": "numpy int64 after explicit absolute-entry bounds",
            "int64_maximum": INT64_MAX,
            "preflight_absolute_entry_bounds": engine.int64_preflight_bounds,
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_boundary_pairing_absolute": pairing_maximum,
            "all_pre_multiplication_bounds_below_2_power_63": True,
        },
        "lineage_and_claim_boundary": {
            "five_facet_headers_already_present_in_five_chamber_P5_pilot": pilot[
                "result"
            ]["facet_counts"],
            "positive_increment": (
                "a complete 16872-type K1,4 shape census, not merely five headers"
            ),
            "P5_complete_census_claimed": False,
            "degree_three_fork_complete_census_claimed": False,
            "all_five_vertex_shapes_complete_census_claimed": False,
            "worldwide_priority": "open",
            "bounded_CoxIter_commit": bounded_corpus["commit"],
            "weak_carrier_not_full_nef": True,
        },
        "orbit_representatives": public_rows,
        "checks": checks,
    }


def _probe(limit: int) -> None:
    start = time.perf_counter()
    _carrier, gram, roots, interior, _bounded = base._load_inputs()
    actions = _restricted_actions()
    rows = _star_color_orbits(actions)
    print(f"orbit_count={len(rows)} orbit_seconds={time.perf_counter() - start:.3f}")
    engine = ExactStarBoundaryEngine(gram, roots, interior)
    facet_histogram: Counter[int] = Counter()
    degree_keys = set()
    boundary_seconds = time.perf_counter()
    boundary_only = 0.0
    ordinary_only = 0.0
    for row in rows[:limit]:
        tick = time.perf_counter()
        boundary, data = engine.boundary(row["colors_zero_based_stock"])
        boundary_only += time.perf_counter() - tick
        tick = time.perf_counter()
        _labels, degree_key, _ordinary = engine.ordinary_data(boundary)
        ordinary_only += time.perf_counter() - tick
        facet_histogram[data["facet_count"]] += 1
        degree_keys.add((data["facet_count"], degree_key))
    elapsed = time.perf_counter() - boundary_seconds
    print(f"probed={limit} seconds={elapsed:.3f} seconds_per={elapsed/max(limit,1):.6f}")
    print(f"boundary_seconds={boundary_only:.3f} ordinary_seconds={ordinary_only:.3f}")
    print(f"facet_histogram={dict(sorted(facet_histogram.items()))}")
    print(f"degree_key_count={len(degree_keys)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=int)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    if arguments.probe is not None:
        _probe(arguments.probe)
        return
    payload = build(verbose=True)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file():
            raise SystemExit(f"archive is missing: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from the exact K1,4 rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print(f"K1_4_types={payload['theorem']['ordinary_Coxeter_type_count']}")
    print(
        "type_facet_histogram="
        + json.dumps(payload["theorem"]["type_facet_histogram"], sort_keys=True)
    )
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
