#!/usr/bin/env python3
"""Certify the complete B<=2 five-vertex P5 layer in the JK-H12 factory.

The four edge colours of a path are taken from the frozen 112-wall Allcock
stock.  The exact seed action preserves its 16-wall A part and 96-wall B
part.  Consequently the ten binary edge-class words modulo path reversal
are disjoint invariant shards.  This production archive closes the seven
shards with at most two B occurrences; the same command can run any one of
the three remaining shards, or all ten together.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import product
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any, Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base
import scripts.certify_jacobian_kummer_h12_four_vertex_subtree_census as four
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star


BASE_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_small_subtree_census.py"
BASE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_small_subtree_census.json"
FOUR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py"
FOUR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_four_vertex_subtree_census.json"
PILOT_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_chamber_pilot.py"
PILOT_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_chamber_pilot.json"
DENOMINATOR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_tree_denominators.py"
DENOMINATOR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_tree_denominators.json"
STAR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py"
STAR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json"

SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_p5_b_le_2_census.v1"
STATUS_DATE = "2026-08-08"
RANK = 13
STOCK_START = 60
STOCK_END = 172
STOCK_SIZE = 112
A_SIZE = 16
B_SIZE = 96
INT64_MAX = np.iinfo(np.int64).max

PINNED_SHA256 = {
    "small_tree_builder": "30ccb14a9240c32dfd8baeb70f1a3482a3941f5a6e4a9523094f4ebc4630d125",
    "small_tree_archive": "796575ebbaaf8ec0ba3054911397d83834d8939e1eb4c20cd24aae1b56aee56d",
    "four_tree_builder": "474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c",
    "four_tree_archive": "2591d5f078421b5b4902045d37fa57e33548e2450f090f5ea8513c564c1a4b24",
    "five_chamber_pilot_builder": "d1e678dc2eea7c67998fff1d5fcef34a8ad45281fac87424f74f3876ea4c9f6d",
    "five_chamber_pilot_archive": "809302f4aa10fcf723c68061c34a4d857c9d5c64a4c03e04227c4c54c917ebec",
    "five_vertex_denominator_builder": "1227008ddd3e337afe6b76ec1481ebd235d7d8e29201e5b772fef2c906e0c5bb",
    "five_vertex_denominator_archive": "2306ebe0048b576a5ba1aa145d5af500b6efdb13f16ffacb8174c7984d1d5d10",
    "five_vertex_star_builder": "bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb",
    "five_vertex_star_archive": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
}

# Lexicographically chosen representatives of binary words modulo reversal.
ALL_PATTERNS = (
    "AAAA",
    "AAAB",
    "AABA",
    "AABB",
    "ABAB",
    "ABBA",
    "BAAB",
    "ABBB",
    "BABB",
    "BBBB",
)
PRODUCTION_PATTERNS = ALL_PATTERNS[:7]
REMAINING_PATTERNS = ALL_PATTERNS[7:]
EXPECTED_PATTERN_ORBITS = {
    "AAAA": 142,
    "AAAB": 1_007,
    "AABA": 1_084,
    "AABB": 5_819,
    "ABAB": 6_304,
    "ABBA": 3_256,
    "BAAB": 3_081,
    "ABBB": 36_247,
    "BABB": 36_648,
    "BBBB": 108_159,
}

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
Permutation = tuple[int, ...]
ColorWord = tuple[int, int, int, int]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(rendered.encode()).hexdigest()


def _pattern_ranges(pattern: str) -> tuple[range, ...]:
    if pattern not in ALL_PATTERNS:
        raise ArithmeticError(f"unknown canonical P5 class word {pattern}")
    return tuple(range(A_SIZE) if letter == "A" else range(A_SIZE, STOCK_SIZE) for letter in pattern)


def _is_palindromic(pattern: str) -> bool:
    return pattern == pattern[::-1]


def _local_code(word: ColorWord, pattern: str) -> int:
    code = 0
    for color, letter in zip(word, pattern, strict=True):
        size = A_SIZE if letter == "A" else B_SIZE
        local = color if letter == "A" else color - A_SIZE
        if not 0 <= local < size:
            raise ArithmeticError("a color left its binary class shard")
        code = code * size + local
    return code


def _raw_pattern_count(pattern: str) -> int:
    count = 1
    for index, letter in enumerate(pattern):
        size = A_SIZE if letter == "A" else B_SIZE
        count *= size - int(index > 0 and pattern[index - 1] == letter)
    if _is_palindromic(pattern):
        if count % 2:
            raise ArithmeticError("a palindromic proper-word count is odd")
        count //= 2
    return count


def _class_fixed_counts(action: Permutation) -> tuple[int, int, int, int, int, int]:
    fixed_a = sum(action[color] == color for color in range(A_SIZE))
    fixed_b = sum(action[color] == color for color in range(A_SIZE, STOCK_SIZE))
    fixed2_a = sum(action[action[color]] == color for color in range(A_SIZE))
    fixed2_b = sum(
        action[action[color]] == color for color in range(A_SIZE, STOCK_SIZE)
    )
    return (
        fixed_a,
        fixed_b,
        fixed2_a,
        fixed2_b,
        fixed2_a - fixed_a,
        fixed2_b - fixed_b,
    )


def _burnside_pattern_orbits(pattern: str, actions: Sequence[Permutation]) -> dict[str, int]:
    identity_sum = 0
    reversal_sum = 0
    palindrome = _is_palindromic(pattern)
    for action in actions:
        fixed_a, fixed_b, fixed2_a, fixed2_b, two_a, two_b = _class_fixed_counts(action)
        fixed = {"A": fixed_a, "B": fixed_b}
        fixed2 = {"A": fixed2_a, "B": fixed2_b}
        genuine_two = {"A": two_a, "B": two_b}
        contribution = 1
        for index, letter in enumerate(pattern):
            contribution *= fixed[letter] - int(
                index > 0 and pattern[index - 1] == letter
            )
        identity_sum += contribution
        if palindrome:
            # x2 lies in a genuine two-cycle.  x1 is fixed by g^2; it must
            # differ from x2 only when the outer and inner classes agree.
            reversal_sum += genuine_two[pattern[1]] * (
                fixed2[pattern[0]] - int(pattern[0] == pattern[1])
            )
    denominator = len(actions) * (2 if palindrome else 1)
    numerator = identity_sum + reversal_sum
    if numerator % denominator:
        raise ArithmeticError("a class-word Burnside numerator does not divide")
    return {
        "identity_fixed_sum": identity_sum,
        "reversal_fixed_sum": reversal_sum,
        "Burnside_denominator": denominator,
        "orbit_count": numerator // denominator,
    }


def _pattern_color_orbits(
    pattern: str, actions: Sequence[Permutation]
) -> tuple[dict[str, Any], ...]:
    ranges = _pattern_ranges(pattern)
    dense_size = math.prod(len(values) for values in ranges)
    seen = bytearray(dense_size)
    palindrome = _is_palindromic(pattern)
    rows: list[dict[str, Any]] = []
    for candidate in product(*ranges):
        word = tuple(int(value) for value in candidate)
        if any(word[index] == word[index + 1] for index in range(3)):
            continue
        if palindrome and word > word[::-1]:
            continue
        code = _local_code(word, pattern)
        if seen[code]:
            continue
        orbit: set[ColorWord] = set()
        for action in actions:
            image = tuple(action[color] for color in word)
            if palindrome:
                image = min(image, image[::-1])
            orbit.add(image)
        for image in orbit:
            seen[_local_code(image, pattern)] = 1
        rows.append(
            {
                "edge_class_word": pattern,
                "colors_zero_based_stock": word,
                "unrooted_coloring_orbit_size_under_seed_action": len(orbit),
            }
        )
    raw_count = _raw_pattern_count(pattern)
    burnside = _burnside_pattern_orbits(pattern, actions)
    if not (
        sum(seen) == raw_count
        and sum(row["unrooted_coloring_orbit_size_under_seed_action"] for row in rows)
        == raw_count
        and len(rows) == EXPECTED_PATTERN_ORBITS[pattern]
        and burnside["orbit_count"] == len(rows)
    ):
        raise ArithmeticError(f"the complete {pattern} P5 color quotient changed")
    return tuple(rows)


class ExactP5BoundaryEngine(star.ExactStarBoundaryEngine):
    """Exact five-step path boundary reconstruction over the rank-13 lattice."""

    def __init__(self, gram: Matrix, roots: tuple[Vector, ...], interior: Vector) -> None:
        super().__init__(gram, roots, interior)
        self.gram_maximum = int(np.max(np.abs(self.gram)))
        reflection_maximum = max(
            int(np.max(np.abs(reflection))) for reflection in self.reflections.values()
        )
        root_maximum = int(np.max(np.abs(self.roots)))
        interior_maximum = int(np.max(np.abs(self.interior)))
        four_transform_bound = RANK**3 * reflection_maximum**4
        transformed_root_bound = RANK * four_transform_bound * root_maximum
        transformed_interior_bound = RANK * four_transform_bound * interior_maximum
        path_bounds = {
            "four_reflection_transform_entry_absolute": four_transform_bound,
            "transformed_seed_root_entry_absolute": transformed_root_bound,
            "transformed_seed_interior_entry_absolute": transformed_interior_bound,
        }
        if max(path_bounds.values()) >= INT64_MAX:
            raise ArithmeticError("the P5 transform int64 preflight failed")
        self.int64_preflight_bounds = {
            **self.int64_preflight_bounds,
            **path_bounds,
        }

    def _transformations(self, stock_colors: ColorWord) -> tuple[np.ndarray, ...]:
        transformations = [self.identity]
        for color in stock_colors:
            transformations.append(
                transformations[-1] @ self.reflections[color + STOCK_START]
            )
        if len({_matrix_key(matrix) for matrix in transformations}) != 5:
            raise ArithmeticError("a proper P5 word did not give five chambers")
        return tuple(transformations)

    def boundary(self, stock_colors: ColorWord) -> tuple[np.ndarray, dict[str, Any]]:
        if not (
            len(stock_colors) == 4
            and all(0 <= color < STOCK_SIZE for color in stock_colors)
            and all(stock_colors[index] != stock_colors[index + 1] for index in range(3))
        ):
            raise ArithmeticError("invalid properly colored P5 word")
        absolute = tuple(color + STOCK_START for color in stock_colors)
        transformations = self._transformations(stock_colors)
        transform_index = {
            _matrix_key(matrix): index for index, matrix in enumerate(transformations)
        }
        adjacency: set[tuple[int, int, int]] = set()
        internal: set[Vector] = set()
        for chamber, transformation in enumerate(transformations):
            for color in set(absolute):
                neighbor = transformation @ self.reflections[color]
                other = transform_index.get(_matrix_key(neighbor))
                if other is None:
                    continue
                left_key = star._primitive_key(transformation @ self.roots[color])
                right_key = star._primitive_key(neighbor @ self.roots[color])
                if left_key != right_key:
                    raise ArithmeticError("the two sides of an internal P5 mirror differ")
                internal.add(left_key)
                if chamber < other:
                    adjacency.add((chamber, other, color))
        expected_adjacency = {
            (index, index + 1, color) for index, color in enumerate(absolute)
        }
        if adjacency != expected_adjacency or len(internal) != 4:
            raise ArithmeticError("the transformations do not induce the marked P5")

        union_interior = sum(
            (transformation @ self.interior for transformation in transformations),
            start=np.zeros(RANK, dtype=np.int64),
        )
        boundary: dict[Vector, Vector] = {}
        maximum_orientation_pairing = 0
        maximum_transformed_root = 0
        for transformation in transformations:
            images = self.roots @ transformation.T
            image_maximum = int(np.max(np.abs(images)))
            maximum_transformed_root = max(maximum_transformed_root, image_maximum)
            orientation_bound = (
                RANK**2
                * self.gram_maximum
                * image_maximum
                * int(np.max(np.abs(union_interior)))
            )
            if orientation_bound >= INT64_MAX:
                raise ArithmeticError("a P5 orientation pairing could overflow int64")
            pairings = images @ self.gram @ union_interior
            maximum_orientation_pairing = max(
                maximum_orientation_pairing, int(np.max(np.abs(pairings)))
            )
            for image, pairing in zip(images, pairings, strict=True):
                key = star._primitive_key(image)
                if key in internal:
                    continue
                if pairing == 0:
                    raise ArithmeticError("the summed P5 interior lies on a boundary")
                oriented = image if pairing > 0 else -image
                existing = boundary.setdefault(key, tuple(int(value) for value in oriented))
                if existing != tuple(int(value) for value in oriented):
                    raise ArithmeticError("a P5 exterior wall received opposite orientations")

        ordered_keys = tuple(sorted(boundary))
        if not self.core_keys <= set(ordered_keys):
            raise ArithmeticError("the fixed Picard core left a P5 boundary")
        roots = np.array([boundary[key] for key in ordered_keys], dtype=np.int64)
        coordinate_maximum = int(np.max(np.abs(roots)))
        pair_bound = RANK**2 * self.gram_maximum * coordinate_maximum**2
        interior_maximum = int(np.max(np.abs(union_interior)))
        interior_bound = RANK**2 * self.gram_maximum * interior_maximum**2
        if max(pair_bound, interior_bound) >= INT64_MAX:
            raise ArithmeticError("a P5 boundary Gram product could overflow int64")
        covectors = roots @ self.gram
        squares = np.einsum("ij,ij->i", roots, covectors)
        if not np.all(squares < 0):
            raise ArithmeticError("a P5 boundary root is not spacelike")
        primitive_integral = 0
        for root, square, covector in zip(roots, squares, covectors, strict=True):
            if not (
                star._content(root) == 1
                and all((2 * int(value)) % int(square) == 0 for value in covector)
            ):
                raise ArithmeticError("a P5 boundary reflection is nonintegral")
            primitive_integral += 1
        interior_pairings = roots @ self.gram @ union_interior
        interior_square = int(union_interior @ self.gram @ union_interior)
        if not (interior_square > 0 and np.all(interior_pairings > 0)):
            raise ArithmeticError("a P5 strict interior failed")
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
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_transformed_seed_root_coordinate_absolute": maximum_transformed_root,
            "maximum_orientation_pairing_absolute": maximum_orientation_pairing,
            "maximum_chamber_transform_coordinate_absolute": max(
                int(np.max(np.abs(matrix))) for matrix in transformations
            ),
        }


def _matrix_key(matrix: np.ndarray) -> bytes:
    return np.ascontiguousarray(matrix, dtype=np.int64).tobytes()


def _pinned_predecessors() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    paths = {
        "small_tree_builder": BASE_SCRIPT,
        "small_tree_archive": BASE_ARCHIVE,
        "four_tree_builder": FOUR_SCRIPT,
        "four_tree_archive": FOUR_ARCHIVE,
        "five_chamber_pilot_builder": PILOT_SCRIPT,
        "five_chamber_pilot_archive": PILOT_ARCHIVE,
        "five_vertex_denominator_builder": DENOMINATOR_SCRIPT,
        "five_vertex_denominator_archive": DENOMINATOR_ARCHIVE,
        "five_vertex_star_builder": STAR_SCRIPT,
        "five_vertex_star_archive": STAR_ARCHIVE,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a pinned P5 prerequisite changed: {hashes}")
    small = json.loads(BASE_ARCHIVE.read_text(encoding="utf-8"))
    pilot = json.loads(PILOT_ARCHIVE.read_text(encoding="utf-8"))
    denominator = json.loads(DENOMINATOR_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        small["all_checks_pass"]
        and pilot["all_checks_pass"]
        and denominator["all_checks_pass"]
        and denominator["tree_shapes"]["P5"]["seed_and_tree_automorphism_orbits"]
        == 201_747
        and pilot["result"]["finite_volume"]
        and not pilot["result"]["compact"]
        and pilot["result"]["reflection_index_and_volume_multiplier"] == 5
    ):
        raise ArithmeticError("a frozen P5 predecessor is not green")
    return small, pilot, denominator


def _canonical_orbit_representative(
    word: ColorWord, pattern: str, actions: Sequence[Permutation]
) -> ColorWord:
    actual_pattern = "".join("A" if color < A_SIZE else "B" for color in word)
    if actual_pattern == pattern[::-1] and actual_pattern != pattern:
        word = word[::-1]
        actual_pattern = pattern
    if actual_pattern != pattern:
        raise ArithmeticError("a pilot word does not belong to its canonical shard")
    images = []
    for action in actions:
        image = tuple(action[color] for color in word)
        if _is_palindromic(pattern):
            image = min(image, image[::-1])
        images.append(image)
    return min(images)


def _selected_patterns(patterns: Iterable[str]) -> tuple[str, ...]:
    chosen = tuple(patterns)
    if not chosen or len(set(chosen)) != len(chosen):
        raise ArithmeticError("the P5 pattern selection is empty or repeated")
    if any(pattern not in ALL_PATTERNS for pattern in chosen):
        raise ArithmeticError("the P5 pattern selection is not canonical")
    return tuple(pattern for pattern in ALL_PATTERNS if pattern in set(chosen))


def build(
    *, patterns: Iterable[str] = PRODUCTION_PATTERNS, verbose: bool = False
) -> dict[str, Any]:
    patterns = _selected_patterns(patterns)
    small, pilot, denominator = _pinned_predecessors()
    carrier, gram, roots, interior, bounded_corpus = base._load_inputs()
    actions = star._restricted_actions()
    if not all(
        all((color < A_SIZE) == (action[color] < A_SIZE) for color in range(STOCK_SIZE))
        for action in actions
    ):
        raise ArithmeticError("the seed action does not preserve the A/B split")
    stock_gram_sha256 = star._root_gram_action_gate(gram, roots, actions)

    color_rows: list[dict[str, Any]] = []
    pattern_denominator: dict[str, dict[str, Any]] = {}
    for pattern in patterns:
        tick = time.perf_counter()
        rows = _pattern_color_orbits(pattern, actions)
        burnside = _burnside_pattern_orbits(pattern, actions)
        pattern_denominator[pattern] = {
            "raw_unrooted_proper_colorings": _raw_pattern_count(pattern),
            **burnside,
            "explicit_orbit_cover_count": len(rows),
            "explicit_orbit_size_sum": sum(
                row["unrooted_coloring_orbit_size_under_seed_action"] for row in rows
            ),
        }
        color_rows.extend(rows)
        if verbose:
            print(
                f"orbit_shard={pattern} orbits={len(rows)} "
                f"seconds={time.perf_counter() - tick:.2f}",
                flush=True,
            )
    expected_type_count = sum(EXPECTED_PATTERN_ORBITS[pattern] for pattern in patterns)
    if len(color_rows) != expected_type_count:
        raise ArithmeticError("the selected P5 orbit denominator changed")

    engine = ExactP5BoundaryEngine(gram, roots, interior)
    degree_groups: defaultdict[tuple[int, bytes], list[int]] = defaultdict(list)
    public_rows: list[dict[str, Any]] = []
    relation_histograms_by_facet: defaultdict[
        int, Counter[tuple[tuple[str, int], ...]]
    ] = defaultdict(Counter)
    coordinate_maximum = 0
    transformed_root_maximum = 0
    transform_maximum = 0
    pairing_maximum = 0
    boundary_started = time.perf_counter()
    for index, color_row in enumerate(color_rows):
        stock_colors = tuple(int(value) for value in color_row["colors_zero_based_stock"])
        boundary, boundary_data = engine.boundary(stock_colors)
        labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        del labels
        pattern = str(color_row["edge_class_word"])
        b_count = pattern.count("B")
        expected_facets = 752 + 9 * b_count
        if boundary_data["facet_count"] != expected_facets:
            raise ArithmeticError("the P5 facet formula changed")
        degree_groups[(len(boundary), degree_key)].append(index)
        relation_histograms_by_facet[len(boundary)][
            tuple(sorted(ordinary_data["ordinary_Coxeter_relation_histogram"].items()))
        ] += 1
        coordinate_maximum = max(
            coordinate_maximum, boundary_data["maximum_boundary_coordinate_absolute"]
        )
        transformed_root_maximum = max(
            transformed_root_maximum,
            boundary_data["maximum_transformed_seed_root_coordinate_absolute"],
        )
        transform_maximum = max(
            transform_maximum,
            boundary_data["maximum_chamber_transform_coordinate_absolute"],
        )
        pairing_maximum = max(
            pairing_maximum, ordinary_data["maximum_boundary_pairing_absolute"]
        )
        repeated_positions = 4 - len(set(stock_colors))
        public_rows.append(
            {
                "orbit_number": index + 1,
                "tree_shape": "P5",
                "edge_class_word": pattern,
                "edge_colors_one_based": [
                    color + STOCK_START + 1 for color in stock_colors
                ],
                "B_edge_occurrence_count": b_count,
                "unrooted_coloring_orbit_size_under_seed_action": color_row[
                    "unrooted_coloring_orbit_size_under_seed_action"
                ],
                "proper_edge_coloring": True,
                "nonadjacent_color_repetition_count": repeated_positions,
                **boundary_data,
                **ordinary_data,
                "expected_facet_count_752_plus_9B": expected_facets,
                "fixed_boundary_core_one_based": "W1--W60",
                "_degree_key": degree_key,
            }
        )
        if verbose and ((index + 1) % 500 == 0 or index + 1 == len(color_rows)):
            print(
                f"boundary_progress={index + 1}/{len(color_rows)} "
                f"elapsed_seconds={time.perf_counter() - boundary_started:.1f}",
                flush=True,
            )

    collision_groups = tuple(
        indices for indices in degree_groups.values() if len(indices) > 1
    )
    collision_indices = {index for group in collision_groups for index in group}
    refinement_keys: dict[int, bytes] = {}
    refinement_discrete = 0
    for replay_number, index in enumerate(sorted(collision_indices), start=1):
        stock_colors = tuple(
            int(value) for value in color_rows[index]["colors_zero_based_stock"]
        )
        boundary, _boundary_data = engine.boundary(stock_colors)
        labels, degree_key, _ordinary_data = engine.ordinary_data(boundary)
        if degree_key != public_rows[index]["_degree_key"]:
            raise ArithmeticError("a residual P5 degree profile changed on replay")
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
        raise ArithmeticError(
            "canonical color refinement left unresolved P5 groups: "
            f"{Counter(len(group) for group in unresolved_groups)}"
        )

    exact_type_keys: list[bytes] = []
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
        raise ArithmeticError("two selected P5 ordinary Coxeter invariants collided")

    representative_lookup = {
        (str(color_row["edge_class_word"]), tuple(color_row["colors_zero_based_stock"])):
        index
        for index, color_row in enumerate(color_rows)
    }
    pilot_recovery = []
    for pilot_row in pilot["pilots"]:
        stock_word = tuple(
            int(color) - STOCK_START - 1 for color in pilot_row["edge_colors_one_based"]
        )
        class_word = "".join("A" if color < A_SIZE else "B" for color in stock_word)
        pattern = min(class_word, class_word[::-1])
        if pattern not in patterns:
            continue
        representative = _canonical_orbit_representative(stock_word, pattern, actions)
        index = representative_lookup[(pattern, representative)]
        pilot_recovery.append(
            {
                "B_edge_occurrence_count": pattern.count("B"),
                "pilot_edge_colors_one_based": pilot_row["edge_colors_one_based"],
                "class_word_shard": pattern,
                "recovered_orbit_number": index + 1,
                "recovered_facet_count": public_rows[index]["facet_count"],
                "same_seed_action_and_path_reversal_orbit": True,
            }
        )

    star_archive = json.loads(STAR_ARCHIVE.read_text(encoding="utf-8"))
    star_degree_hashes = {
        (row["facet_count"], row["ordinary_degree_profile_multiset_sha256"])
        for row in star_archive["orbit_representatives"]
    }
    star_type_hashes = {
        row["ordinary_Coxeter_type_invariant_sha256"]
        for row in star_archive["orbit_representatives"]
    }
    cross_star_degree_hash_matches = sum(
        (row["facet_count"], row["ordinary_degree_profile_multiset_sha256"])
        in star_degree_hashes
        for row in public_rows
    )
    cross_star_type_hash_matches = sum(
        row["ordinary_Coxeter_type_invariant_sha256"] in star_type_hashes
        for row in public_rows
    )

    pattern_type_histogram = Counter(row["edge_class_word"] for row in public_rows)
    facet_type_histogram = Counter(row["facet_count"] for row in public_rows)
    raw_facet_histogram = Counter()
    orbit_size_histogram = Counter()
    repetition_histogram = Counter()
    for row in public_rows:
        orbit_size = row["unrooted_coloring_orbit_size_under_seed_action"]
        raw_facet_histogram[row["facet_count"]] += orbit_size
        orbit_size_histogram[orbit_size] += 1
        repetition_histogram[row["nonadjacent_color_repetition_count"]] += 1
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
    production_scope = patterns == PRODUCTION_PATTERNS
    expected_raw_total = sum(_raw_pattern_count(pattern) for pattern in patterns)
    checks = {
        "frozen_112_wall_factory_and_predecessors_replayed": True,
        "exact_384_action_preserves_stock_Gram_and_AB_split": True,
        "every_selected_binary_reversal_shard_has_exact_explicit_and_Burnside_denominator": all(
            pattern_denominator[pattern]["orbit_count"]
            == EXPECTED_PATTERN_ORBITS[pattern]
            == pattern_denominator[pattern]["explicit_orbit_cover_count"]
            and pattern_denominator[pattern]["raw_unrooted_proper_colorings"]
            == pattern_denominator[pattern]["explicit_orbit_size_sum"]
            for pattern in patterns
        ),
        "all_selected_orbit_representative_boundaries_rebuilt": len(public_rows)
        == expected_type_count,
        "every_chamber_set_is_an_embedded_P5_with_four_internal_mirrors": all(
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
        "exact_degree_profiles_plus_canonical_refinement_separate_every_selected_type": True,
        "all_selected_ordinary_Coxeter_matrices_are_pairwise_nonisomorphic": len(
            set(exact_type_keys)
        )
        == expected_type_count,
        "Allcock_index_and_volume_multiplier_are_five": True,
        "all_descendants_are_finite_volume_noncompact": True,
        "same_full_Picard_K3_weak_carrier_retained": True,
        "production_scope_is_exactly_B_at_most_2": (not production_scope)
        or (patterns == PRODUCTION_PATTERNS and expected_type_count == 20_693),
        "full_P5_census_not_claimed_for_production_scope": (not production_scope)
        or expected_type_count < 201_747,
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"five-vertex P5 shard census failed: {checks}")

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
            "A_stock": "W61--W76",
            "A_stock_size": A_SIZE,
            "B_stock": "W77--W172",
            "B_stock_size": B_SIZE,
            "folded_seed_action_order": len(actions),
            "stock_Gram_sha256": stock_gram_sha256,
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "source_carrier_theorem": carrier["theorem"],
        },
        "coloring_denominator": {
            "shape": "P5",
            "selected_binary_class_word_shards": list(patterns),
            "selected_shard_count": len(patterns),
            "selected_raw_unrooted_proper_colorings": expected_raw_total,
            "selected_seed_and_reversal_orbit_count": expected_type_count,
            "per_pattern": pattern_denominator,
            "all_ten_pattern_orbit_count": sum(EXPECTED_PATTERN_ORBITS.values()),
            "all_shape_denominator_P5_orbit_count": denominator["tree_shapes"]["P5"][
                "seed_and_tree_automorphism_orbits"
            ],
            "remaining_disjoint_patterns": [
                pattern for pattern in ALL_PATTERNS if pattern not in patterns
            ],
            "remaining_disjoint_orbit_count": sum(
                EXPECTED_PATTERN_ORBITS[pattern]
                for pattern in ALL_PATTERNS
                if pattern not in patterns
            ),
            "seed_orbit_size_histogram": {
                str(size): count for size, count in sorted(orbit_size_histogram.items())
            },
        },
        "theorem": {
            "scope": (
                "complete P5 binary class-word layer with B-edge occurrence at most 2"
                if production_scope
                else "complete selected P5 binary class-word shards"
            ),
            "ordinary_Coxeter_type_count": len(public_rows),
            "pairwise_nonisomorphic_within_selected_layer": True,
            "facet_formula": "f(P5)=172+145*4+9*B=752+9*B",
            "pattern_type_histogram": dict(sorted(pattern_type_histogram.items())),
            "type_facet_histogram": {
                str(facets): count for facets, count in sorted(facet_type_histogram.items())
            },
            "raw_unrooted_coloring_facet_histogram": {
                str(facets): count for facets, count in sorted(raw_facet_histogram.items())
            },
            "nonadjacent_color_repetition_type_histogram": {
                str(count): multiplicity
                for count, multiplicity in sorted(repetition_histogram.items())
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
        },
        "pilot_recovery": {
            "eligible_B_0_1_2_pilot_count": len(pilot_recovery),
            "rows": pilot_recovery,
            "project_positive_increment_over_the_five_path_pilot_within_this_layer": len(
                public_rows
            )
            - len(pilot_recovery),
        },
        "cross_shape_novelty_screen": {
            "comparison_archive": str(STAR_ARCHIVE.relative_to(ROOT)),
            "complete_star_type_count": len(star_archive["orbit_representatives"]),
            "degree_profile_SHA256_match_count": cross_star_degree_hash_matches,
            "combined_type_invariant_SHA256_match_count": cross_star_type_hash_matches,
            "interpretation": (
                "a collision-resistant archival screen, not a pure exact proof of "
                "cross-shape nonisomorphism when the match count is zero"
            ),
        },
        "Allcock_and_carrier": {
            "finite_subtree_theorem": (
                "every representative is a five-copy P5 fundamental polytope of "
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
            "integer_backend": "numpy int64 after explicit global and dynamic bounds",
            "int64_maximum": INT64_MAX,
            "preflight_absolute_entry_bounds": engine.int64_preflight_bounds,
            "maximum_chamber_transform_coordinate_absolute": transform_maximum,
            "maximum_transformed_seed_root_coordinate_absolute": transformed_root_maximum,
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_boundary_pairing_absolute": pairing_maximum,
            "all_pre_multiplication_bounds_below_2_power_63": True,
        },
        "resumable_full_census_boundary": {
            "all_ten_binary_reversal_shards": list(ALL_PATTERNS),
            "production_completed_shards": list(patterns),
            "remaining_shards": [
                pattern for pattern in ALL_PATTERNS if pattern not in patterns
            ],
            "remaining_orbit_count": sum(
                EXPECTED_PATTERN_ORBITS[pattern]
                for pattern in ALL_PATTERNS
                if pattern not in patterns
            ),
            "single_shard_CLI": "--pattern ABBB (or BABB or BBBB) --output PATH",
            "full_CLI": "--all-patterns --output PATH",
            "shards_are_pairwise_disjoint_because_the_seed_action_preserves_A_and_B": True,
            "a_final_cross_shard_exact_nonisomorphism_merge_is_still_required": not (
                patterns == ALL_PATTERNS
            ),
        },
        "lineage_and_claim_boundary": {
            "positive_increment_over_the_five_P5_pilots_inside_completed_layer": len(
                public_rows
            )
            - len(pilot_recovery),
            "complete_B_at_most_2_layer_claimed": production_scope,
            "complete_P5_census_claimed": patterns == ALL_PATTERNS,
            "degree_three_fork_complete_census_claimed": False,
            "all_five_vertex_shapes_complete_census_claimed": False,
            "cross_star_pure_exact_nonisomorphism_claimed": False,
            "worldwide_priority": "open",
            "bounded_CoxIter_commit": bounded_corpus["commit"],
            "weak_carrier_not_full_nef": True,
        },
        "orbit_representatives": public_rows,
        "checks": checks,
    }


def _probe(patterns: tuple[str, ...], limit: int) -> None:
    _small, _pilot, _denominator = _pinned_predecessors()
    _carrier, gram, roots, interior, _bounded = base._load_inputs()
    actions = star._restricted_actions()
    rows = []
    started = time.perf_counter()
    for pattern in patterns:
        rows.extend(_pattern_color_orbits(pattern, actions))
    print(f"orbit_count={len(rows)} orbit_seconds={time.perf_counter() - started:.3f}")
    engine = ExactP5BoundaryEngine(gram, roots, interior)
    facet_histogram: Counter[int] = Counter()
    degree_keys = set()
    tick = time.perf_counter()
    for row in rows[:limit]:
        boundary, data = engine.boundary(row["colors_zero_based_stock"])
        _labels, degree_key, _ordinary = engine.ordinary_data(boundary)
        facet_histogram[data["facet_count"]] += 1
        degree_keys.add((data["facet_count"], degree_key))
    elapsed = time.perf_counter() - tick
    print(f"probed={min(limit, len(rows))} seconds={elapsed:.3f}")
    print(f"seconds_per={elapsed/max(min(limit, len(rows)), 1):.6f}")
    print(f"facet_histogram={dict(sorted(facet_histogram.items()))}")
    print(f"degree_key_count={len(degree_keys)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--pattern", action="append", choices=ALL_PATTERNS)
    selection.add_argument("--all-patterns", action="store_true")
    parser.add_argument("--probe", type=int)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    if arguments.all_patterns:
        patterns = ALL_PATTERNS
    elif arguments.pattern:
        patterns = _selected_patterns(arguments.pattern)
    else:
        patterns = PRODUCTION_PATTERNS
    if patterns != PRODUCTION_PATTERNS and arguments.output == OUTPUT and arguments.probe is None:
        raise SystemExit("nonproduction shard runs require an explicit --output path")
    if arguments.probe is not None:
        _probe(patterns, arguments.probe)
        return
    payload = build(patterns=patterns, verbose=True)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file():
            raise SystemExit(f"archive is missing: {arguments.output}")
        archived = json.loads(arguments.output.read_text(encoding="utf-8"))
        if archived != payload:
            raise SystemExit("archive differs from the exact P5 shard rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print(f"selected_patterns={','.join(patterns)}")
    print(f"P5_types={payload['theorem']['ordinary_Coxeter_type_count']}")
    print(
        "type_facet_histogram="
        + json.dumps(payload["theorem"]["type_facet_histogram"], sort_keys=True)
    )
    print(f"payload_sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
