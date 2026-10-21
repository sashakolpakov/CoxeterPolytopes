#!/usr/bin/env python3
"""Certify every four-vertex subtree in the JK-H12 112-wall factory.

There are two unlabelled four-vertex trees: the path P4 and the star K1,3.
Path edge colors (a,b,c) satisfy a != b and b != c, allow a == c, and are
identified with (c,b,a).  Star colors are unordered triples of distinct
generators.  This script quotients both exact coloring sets by the frozen
384-element seed action before reconstructing any reflected boundary.

Every orbit representative is rebuilt over the rank-13 Picard lattice.  The
ordinary Coxeter matrix (parallel and ultraparallel both labelled infinity)
is checked directly.  Exact vertex-degree profiles separate almost every
representative; canonical color refinement separates the complete residue.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base


ROOT = Path(__file__).resolve().parents[1]
BASE_ARCHIVE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_small_subtree_census.json"
)
OUTPUT = (
    ROOT / "data/enumerations/jacobian_kummer_h12_four_vertex_subtree_census.json"
)

SCHEMA = "coxeter.jacobian_kummer.h12_four_vertex_subtree_census.v1"
STATUS_DATE = "2026-08-08"
BASE_SCRIPT_SHA256 = (
    "30ccb14a9240c32dfd8baeb70f1a3482a3941f5a6e4a9523094f4ebc4630d125"
)
BASE_ARCHIVE_SHA256 = (
    "796575ebbaaf8ec0ba3054911397d83834d8939e1eb4c20cd24aae1b56aee56d"
)

RANK = 13
SEED_FACETS = 172
STOCK_START = 60  # zero based: W61
STOCK_END = 172
STOCK_SIZE = STOCK_END - STOCK_START
LOW_END = 76  # zero based: W76 is the last 317-type wall
INT64_MAX = np.iinfo(np.int64).max
LABEL_NAMES = ("m2", "m3", "m4", "m6", "infinity")

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
Permutation = tuple[int, ...]
ColorTriple = tuple[int, int, int]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_sha256(payload: Any) -> str:
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return _bytes_sha256(rendered.encode())


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


def _matrix_key(matrix: np.ndarray) -> bytes:
    return np.ascontiguousarray(matrix, dtype=np.int64).tobytes()


def _restricted_actions(actions: Sequence[Permutation]) -> tuple[Permutation, ...]:
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


def _path_color_orbits(
    actions: Sequence[Permutation],
) -> tuple[dict[str, Any], ...]:
    seen = bytearray(STOCK_SIZE**3)
    rows: list[dict[str, Any]] = []
    for left in range(STOCK_SIZE):
        for middle in range(STOCK_SIZE):
            if left == middle:
                continue
            for right in range(STOCK_SIZE):
                if middle == right:
                    continue
                code = (left * STOCK_SIZE + middle) * STOCK_SIZE + right
                if seen[code]:
                    continue
                orbit: set[int] = set()
                for action in actions:
                    image_left = action[left]
                    image_middle = action[middle]
                    image_right = action[right]
                    orbit.add(
                        (image_left * STOCK_SIZE + image_middle) * STOCK_SIZE
                        + image_right
                    )
                    orbit.add(
                        (image_right * STOCK_SIZE + image_middle) * STOCK_SIZE
                        + image_left
                    )
                for image in orbit:
                    seen[image] = 1
                palindromic = left == right
                unrooted_size = len(orbit) if palindromic else len(orbit) // 2
                rows.append(
                    {
                        "shape": "P4",
                        "colors_zero_based": (
                            left + STOCK_START,
                            middle + STOCK_START,
                            right + STOCK_START,
                        ),
                        "outer_colors_repeated": palindromic,
                        "ordered_sequence_orbit_size_under_seed_and_reversal": len(
                            orbit
                        ),
                        "unrooted_coloring_orbit_size_under_seed_action": unrooted_size,
                    }
                )
    if not (
        sum(seen) == 1_379_952
        and len(rows) == 1_964
        and sum(
            row["unrooted_coloring_orbit_size_under_seed_action"] for row in rows
        )
        == 696_192
        and sum(row["outer_colors_repeated"] for row in rows) == 47
    ):
        raise ArithmeticError("the complete P4 color-orbit census changed")
    return tuple(rows)


def _star_color_orbits(
    actions: Sequence[Permutation],
) -> tuple[dict[str, Any], ...]:
    seen: set[ColorTriple] = set()
    rows: list[dict[str, Any]] = []
    for left in range(STOCK_SIZE):
        for middle in range(left + 1, STOCK_SIZE):
            for right in range(middle + 1, STOCK_SIZE):
                colors = (left, middle, right)
                if colors in seen:
                    continue
                orbit = {
                    tuple(sorted((action[left], action[middle], action[right])))
                    for action in actions
                }
                seen.update(orbit)
                rows.append(
                    {
                        "shape": "K1,3",
                        "colors_zero_based": (
                            left + STOCK_START,
                            middle + STOCK_START,
                            right + STOCK_START,
                        ),
                        "unrooted_coloring_orbit_size_under_seed_action": len(orbit),
                    }
                )
    if not (
        len(seen) == 227_920
        and len(rows) == 680
        and sum(
            row["unrooted_coloring_orbit_size_under_seed_action"] for row in rows
        )
        == 227_920
    ):
        raise ArithmeticError("the complete K1,3 color-orbit census changed")
    return tuple(rows)


class ExactBoundaryEngine:
    """Fast integer boundary reconstruction with explicit overflow bounds."""

    def __init__(
        self,
        gram: Matrix,
        roots: tuple[Vector, ...],
        interior: Vector,
    ) -> None:
        self.gram_tuple = gram
        self.roots_tuple = roots
        self.interior_tuple = interior
        self.gram = np.array(gram, dtype=np.int64)
        self.roots = np.array(roots, dtype=np.int64)
        self.interior = np.array(interior, dtype=np.int64)
        self.identity = np.eye(RANK, dtype=np.int64)
        self.reflections = {
            index: np.array(
                base._reflection_matrix(gram, roots[index]), dtype=np.int64
            )
            for index in range(STOCK_START, STOCK_END)
        }
        self.gram_maximum = int(np.max(np.abs(self.gram)))
        reflection_maximum = max(
            int(np.max(np.abs(matrix))) for matrix in self.reflections.values()
        )
        root_maximum = int(np.max(np.abs(self.roots)))
        interior_maximum = int(np.max(np.abs(self.interior)))
        # If every entry of A_i has absolute value at most M, every entry of
        # A_1...A_k is bounded by RANK^(k-1) M^k.  The internal-mirror test
        # below can form a product of four reflections.  A chamber transform
        # has at most three reflections and is also applied to every root and
        # to the seed interior; four transformed interiors are then summed.
        # These checks execute before any of those products.
        transform_bound = RANK**2 * reflection_maximum**3
        four_reflection_product_bound = RANK**3 * reflection_maximum**4
        image_bound = RANK**3 * reflection_maximum**3 * root_maximum
        transformed_interior_bound = (
            RANK**3 * reflection_maximum**3 * interior_maximum
        )
        summed_interior_bound = 4 * transformed_interior_bound
        self.int64_preflight_bounds = {
            "reflection_entry_absolute": reflection_maximum,
            "seed_root_entry_absolute": root_maximum,
            "seed_interior_entry_absolute": interior_maximum,
            "three_reflection_transform_entry_absolute": transform_bound,
            "four_reflection_product_entry_absolute": four_reflection_product_bound,
            "transformed_seed_root_entry_absolute": image_bound,
            "transformed_seed_interior_entry_absolute": transformed_interior_bound,
            "sum_of_four_transformed_interiors_entry_absolute": summed_interior_bound,
        }
        if max(self.int64_preflight_bounds.values()) >= INT64_MAX:
            raise ArithmeticError("the four-chamber int64 transform bound failed")

        core = sp.Matrix(roots[:60])
        smith = smith_normal_form(core, domain=ZZ)
        self.core_smith = [
            abs(int(smith[index, index]))
            for index in range(min(smith.rows, smith.cols))
            if smith[index, index]
        ]
        self.core_keys = {_primitive_key(root) for root in roots[:60]}
        if self.core_smith != [1] * RANK:
            raise ArithmeticError("the fixed W1--W60 Smith core changed")

    def _transformations(
        self, shape: str, colors: ColorTriple
    ) -> tuple[np.ndarray, ...]:
        left, middle, right = colors
        if shape == "P4":
            first = self.reflections[left]
            second = first @ self.reflections[middle]
            third = second @ self.reflections[right]
            return self.identity, first, second, third
        if shape == "K1,3":
            return (
                self.identity,
                self.reflections[left],
                self.reflections[middle],
                self.reflections[right],
            )
        raise ArithmeticError(f"unknown four-vertex shape {shape}")

    def boundary(self, shape: str, colors: ColorTriple) -> tuple[np.ndarray, dict[str, Any]]:
        transformations = self._transformations(shape, colors)
        transform_keys = {_matrix_key(matrix) for matrix in transformations}
        internal: set[Vector] = set()
        for transformation in transformations:
            for color in set(colors):
                if _matrix_key(transformation @ self.reflections[color]) in transform_keys:
                    internal.add(
                        _primitive_key(transformation @ self.roots[color])
                    )
        if len(internal) != 3:
            raise ArithmeticError(
                f"{shape} {colors} has {len(internal)} rather than three internal mirrors"
            )

        union_interior = sum(
            (transformation @ self.interior for transformation in transformations),
            start=np.zeros(RANK, dtype=np.int64),
        )
        boundary: dict[Vector, Vector] = {}
        maximum_orient_pairing = 0
        for transformation in transformations:
            images = self.roots @ transformation.T
            image_maximum = int(np.max(np.abs(images)))
            interior_maximum = int(np.max(np.abs(union_interior)))
            orient_bound = (
                RANK**2
                * self.gram_maximum
                * image_maximum
                * interior_maximum
            )
            if orient_bound >= INT64_MAX:
                raise ArithmeticError("an orientation pairing could overflow int64")
            pairings = images @ self.gram @ union_interior
            maximum_orient_pairing = max(
                maximum_orient_pairing, int(np.max(np.abs(pairings)))
            )
            for image, pairing in zip(images, pairings, strict=True):
                key = _primitive_key(image)
                if key in internal:
                    continue
                if pairing == 0:
                    raise ArithmeticError("the summed interior lies on a boundary")
                oriented = image if pairing > 0 else -image
                boundary.setdefault(key, tuple(int(value) for value in oriented))

        ordered_keys = tuple(sorted(boundary))
        if not self.core_keys <= set(ordered_keys):
            raise ArithmeticError("the fixed Picard core left a four-chamber boundary")
        roots = np.array([boundary[key] for key in ordered_keys], dtype=np.int64)
        coordinate_maximum = int(np.max(np.abs(roots)))
        interior_maximum = int(np.max(np.abs(union_interior)))
        interior_square_bound = (
            RANK**2 * self.gram_maximum * interior_maximum**2
        )
        pair_bound = (
            RANK**2 * self.gram_maximum * coordinate_maximum**2
        )
        if max(pair_bound, interior_square_bound) >= INT64_MAX:
            raise ArithmeticError("a boundary/interior Gram entry could overflow int64")
        oriented_pairings = roots @ self.gram @ union_interior
        if not (
            np.all(oriented_pairings > 0)
            and int(union_interior @ self.gram @ union_interior) > 0
        ):
            raise ArithmeticError("a four-chamber strict interior failed")
        return roots, {
            "facet_count": len(roots),
            "internal_mirror_count": len(internal),
            "strict_interior_square": int(
                union_interior @ self.gram @ union_interior
            ),
            "minimum_strict_interior_facet_pairing": int(
                np.min(oriented_pairings)
            ),
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_orientation_pairing_absolute": maximum_orient_pairing,
            "marked_boundary_root_sha256": _json_sha256(roots.tolist()),
        }

    def ordinary_coxeter_data(
        self, roots: np.ndarray
    ) -> tuple[np.ndarray, bytes, dict[str, Any]]:
        pairing = roots @ self.gram @ roots.T
        maximum_pairing = int(np.max(np.abs(pairing)))
        if 4 * maximum_pairing**2 >= INT64_MAX:
            raise ArithmeticError("a normalized-square comparison could overflow int64")
        squares = np.diag(pairing)
        if not np.all(squares < 0):
            raise ArithmeticError("a boundary root is not spacelike")
        if 3 * int(np.max(np.abs(squares))) ** 2 >= INT64_MAX:
            raise ArithmeticError("a root-square product could overflow int64")
        square_pairing = pairing * pairing
        square_products = squares[:, None] * squares[None, :]
        count = len(roots)
        diagonal = np.eye(count, dtype=bool)
        masks = (
            (pairing == 0) & ~diagonal,
            (4 * square_pairing == square_products) & ~diagonal,
            (2 * square_pairing == square_products) & ~diagonal,
            (4 * square_pairing == 3 * square_products) & ~diagonal,
            (square_pairing >= square_products) & ~diagonal,
        )
        coverage = sum(mask.astype(np.int8) for mask in masks)
        if not (
            np.all(pairing[~diagonal] >= 0)
            and np.all(coverage[~diagonal] == 1)
            and np.all(coverage[diagonal] == 0)
        ):
            raise ArithmeticError("a boundary pair is not an ordinary Coxeter relation")
        labels = np.full((count, count), -1, dtype=np.int8)
        for label, mask in enumerate(masks):
            labels[mask] = label

        covectors = roots @ self.gram
        primitive_integral = 0
        for root, square, covector in zip(
            roots, squares, covectors, strict=True
        ):
            if not (
                _content(root) == 1
                and all((2 * int(value)) % int(square) == 0 for value in covector)
            ):
                raise ArithmeticError("a transported boundary reflection is nonintegral")
            primitive_integral += 1

        degree_columns = [mask.sum(axis=1).astype(np.uint16) for mask in masks]
        profiles = np.column_stack(degree_columns)
        order = np.lexsort(tuple(profiles[:, column] for column in reversed(range(5))))
        profiles = np.ascontiguousarray(profiles[order], dtype="<u2")
        profile_bytes = profiles.tobytes()
        profile_rows = [tuple(int(value) for value in row) for row in profiles]
        histogram = {
            name: int(mask.sum() // 2)
            for name, mask in zip(LABEL_NAMES, masks, strict=True)
        }
        if sum(histogram.values()) != count * (count - 1) // 2:
            raise ArithmeticError("the ordinary relation inventory is incomplete")
        return labels, profile_bytes, {
            "ordinary_Coxeter_relation_histogram": histogram,
            "ordinary_degree_profile_multiset_sha256": _bytes_sha256(profile_bytes),
            "ordinary_degree_profile_class_count": len(set(profile_rows)),
            "maximum_ordinary_degree_profile_multiplicity": max(
                Counter(profile_rows).values()
            ),
            "primitive_integral_reflection_root_count": primitive_integral,
            "maximum_boundary_pairing_absolute": maximum_pairing,
        }


def _canonical_refinement_colors(labels: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Return canonically numbered stable 1-WL colors and round sizes.

    Initial degree vectors and every later signature are lexicographically
    sorted before IDs are assigned.  Inductively, the numeric IDs therefore
    commute with every simultaneous row/column permutation, rather than being
    dependent on first vertex occurrence.
    """
    count = len(labels)
    degree_profiles = np.column_stack(
        [(labels == label).sum(axis=1) for label in range(5)]
    )
    unique, colors = np.unique(degree_profiles, axis=0, return_inverse=True)
    del unique
    colors = colors.astype(np.int32)
    round_color_counts = [int(colors.max()) + 1]
    vertex_indices = np.arange(count)
    while round_color_counts[-1] < count:
        color_count = round_color_counts[-1]
        signatures: list[tuple[int, ...]] = []
        for vertex in range(count):
            combined = labels[vertex].astype(np.int64) * color_count + colors
            neighbor_counts = np.bincount(
                combined[vertex_indices != vertex], minlength=5 * color_count
            )
            signatures.append(
                (int(colors[vertex]), *(int(value) for value in neighbor_counts))
            )
        palette = {
            signature: index for index, signature in enumerate(sorted(set(signatures)))
        }
        new_colors = np.array(
            [palette[signature] for signature in signatures], dtype=np.int32
        )
        round_color_counts.append(len(palette))
        if np.array_equal(colors, new_colors):
            break
        colors = new_colors
    return colors, round_color_counts


def _canonical_color_refinement(labels: np.ndarray) -> tuple[bytes, dict[str, Any]]:
    """Return an isomorphism invariant of an ordinary edge-labelled matrix.

    When refinement becomes discrete, the resulting adjacency word is an
    exact canonical form.  Otherwise the canonical stable quotient is still
    a rigorous invariant, sufficient here because all residual words differ.
    """

    count = len(labels)
    colors, round_color_counts = _canonical_refinement_colors(labels)

    discrete = round_color_counts[-1] == count
    if discrete:
        canonical_order = np.argsort(colors)
        canonical = labels[np.ix_(canonical_order, canonical_order)]
        upper = np.ascontiguousarray(
            canonical[np.triu_indices(count, 1)], dtype=np.int8
        ).tobytes()
        invariant = b"discrete\0" + count.to_bytes(2, "big") + upper
    else:
        upper_left, upper_right = np.triu_indices(count, 1)
        left_colors = colors[upper_left]
        right_colors = colors[upper_right]
        quotient_rows = np.column_stack(
            (
                labels[upper_left, upper_right],
                np.minimum(left_colors, right_colors),
                np.maximum(left_colors, right_colors),
            )
        )
        quotient_types, multiplicities = np.unique(
            quotient_rows, axis=0, return_counts=True
        )
        invariant_payload = {
            "round_color_counts": round_color_counts,
            "cell_size_histogram": dict(sorted(Counter(colors.tolist()).items())),
            "quotient_edge_inventory": [
                [*(int(value) for value in row), int(multiplicity)]
                for row, multiplicity in zip(
                    quotient_types, multiplicities, strict=True
                )
            ],
        }
        invariant = b"stable\0" + json.dumps(
            invariant_payload, sort_keys=True, separators=(",", ":")
        ).encode()
    return invariant, {
        "color_refinement_applied": True,
        "color_refinement_round_color_counts": round_color_counts,
        "color_refinement_final_cell_count": round_color_counts[-1],
        "color_refinement_maximum_final_cell_size": max(
            Counter(colors.tolist()).values()
        ),
        "color_refinement_discrete_canonical_form": discrete,
        "color_refinement_invariant_sha256": _bytes_sha256(invariant),
    }


def _audit_refinement_permutation_invariance(
    labels: np.ndarray, expected_invariant: bytes, expected_discrete: bool
) -> dict[str, Any]:
    """Check canonical color IDs and the final word under vertex relabelling."""

    count = len(labels)
    colors, rounds = _canonical_refinement_colors(labels)
    permutations = {
        "reverse": np.arange(count - 1, -1, -1),
        "even_vertices_then_odd_vertices": np.concatenate(
            (np.arange(0, count, 2), np.arange(1, count, 2))
        ),
    }
    if expected_discrete != (rounds[-1] == count):
        raise ArithmeticError("a selected refinement audit case changed type")
    for name, permutation in permutations.items():
        if len(np.unique(permutation)) != count:
            raise ArithmeticError(f"{name} is not a vertex permutation")
        permuted = labels[np.ix_(permutation, permutation)]
        permuted_colors, permuted_rounds = _canonical_refinement_colors(permuted)
        invariant, data = _canonical_color_refinement(permuted)
        if not (
            permuted_rounds == rounds
            and np.array_equal(permuted_colors, colors[permutation])
            and invariant == expected_invariant
            and data["color_refinement_discrete_canonical_form"]
            == expected_discrete
        ):
            raise ArithmeticError(
                f"canonical refinement failed the {name} permutation audit"
            )
    return {
        "facet_count": count,
        "discrete": expected_discrete,
        "tested_vertex_permutations": list(permutations),
        "canonical_numeric_color_IDs_equivariant": True,
        "exact_final_invariant_word_unchanged": True,
        "round_color_counts_unchanged": True,
    }


def _public_color_row(row: dict[str, Any], orbit_number: int) -> dict[str, Any]:
    colors = tuple(int(value) for value in row["colors_zero_based"])
    public = {
        key: value
        for key, value in row.items()
        if not key.startswith("_") and key != "colors_zero_based"
    }
    public.update(
        {
            "orbit_number": orbit_number,
            "edge_colors_one_based": [value + 1 for value in colors],
            "edge_color_class_word": "".join(
                "A" if value < LOW_END else "B" for value in colors
            ),
            "B_edge_occurrence_count": sum(value >= LOW_END for value in colors),
        }
    )
    return public


def build() -> dict[str, Any]:
    actual_base_script = _sha256(Path(base.__file__))
    actual_base_archive = _sha256(BASE_ARCHIVE)
    if not (
        actual_base_script == BASE_SCRIPT_SHA256
        and actual_base_archive == BASE_ARCHIVE_SHA256
    ):
        raise ArithmeticError(
            f"a frozen small-tree prerequisite changed: "
            f"{actual_base_script}, {actual_base_archive}"
        )
    previous = json.loads(BASE_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        previous["all_checks_pass"]
        and previous["selected_stock"]["wall_count"] == 112
        and previous["selected_stock"]["all_doubling"]
        and previous["selected_stock"]["all_pairwise_disjoint_in_Allcock_sense"]
        and previous["selected_stock"]["folded_chamber_action_order"] == 384
    ):
        raise ArithmeticError("the frozen 112-wall factory is not green")

    carrier, gram, roots, interior, bounded_corpus = base._load_inputs()
    actions = base._folded_chamber_actions()
    restricted = _restricted_actions(actions)
    path_rows = list(_path_color_orbits(restricted))
    star_rows = list(_star_color_orbits(restricted))
    color_rows = path_rows + star_rows
    if len(color_rows) != 2_644:
        raise ArithmeticError("the combined four-vertex color quotient changed")

    engine = ExactBoundaryEngine(gram, roots, interior)
    degree_groups: defaultdict[tuple[int, bytes], list[int]] = defaultdict(list)
    public_rows: list[dict[str, Any]] = []
    relation_histograms_by_facet: defaultdict[int, set[tuple[tuple[str, int], ...]]] = (
        defaultdict(set)
    )
    coordinate_maximum = 0
    pairing_maximum = 0
    for index, row in enumerate(color_rows):
        shape = row["shape"]
        colors = tuple(int(value) for value in row["colors_zero_based"])
        boundary, boundary_data = engine.boundary(shape, colors)
        labels, degree_key, coxeter_data = engine.ordinary_coxeter_data(boundary)
        del labels
        b_count = sum(color >= LOW_END for color in colors)
        expected_facets = 607 + 9 * b_count
        if boundary_data["facet_count"] != expected_facets:
            raise ArithmeticError("the four-vertex facet formula changed")
        degree_groups[(len(boundary), degree_key)].append(index)
        relation_histograms_by_facet[len(boundary)].add(
            tuple(sorted(coxeter_data["ordinary_Coxeter_relation_histogram"].items()))
        )
        coordinate_maximum = max(
            coordinate_maximum,
            boundary_data["maximum_boundary_coordinate_absolute"],
        )
        pairing_maximum = max(
            pairing_maximum, coxeter_data["maximum_boundary_pairing_absolute"]
        )
        public_rows.append(
            {
                **_public_color_row(row, index + 1),
                **boundary_data,
                **coxeter_data,
                "expected_facet_count_607_plus_9B": expected_facets,
                "fixed_boundary_core_one_based": "W1--W60",
                "fixed_boundary_core_Smith_diagonal": engine.core_smith,
                "full_linear_rank": RANK,
                "_degree_key": degree_key,
            }
        )

    collision_groups = tuple(
        indices for indices in degree_groups.values() if len(indices) > 1
    )
    collision_indices = {index for group in collision_groups for index in group}
    if not (
        len(degree_groups) == 2_535
        and len(collision_groups) == 102
        and len(collision_indices) == 211
        and max(map(len, collision_groups)) == 4
    ):
        raise ArithmeticError("the ordinary degree-profile residue changed")

    refinement_keys: dict[int, bytes] = {}
    refinement_discrete = 0
    refinement_audit_inputs: dict[bool, tuple[int, np.ndarray, bytes]] = {}
    for index in sorted(collision_indices):
        row = color_rows[index]
        colors = tuple(int(value) for value in row["colors_zero_based"])
        boundary, _boundary_data = engine.boundary(row["shape"], colors)
        labels, degree_key, _coxeter_data = engine.ordinary_coxeter_data(boundary)
        if degree_key != public_rows[index]["_degree_key"]:
            raise ArithmeticError("a residual degree profile changed on replay")
        refinement_key, refinement_data = _canonical_color_refinement(labels)
        refinement_keys[index] = refinement_key
        discrete = bool(
            refinement_data["color_refinement_discrete_canonical_form"]
        )
        refinement_discrete += int(discrete)
        refinement_audit_inputs.setdefault(
            discrete, (index, labels.copy(), refinement_key)
        )
        public_rows[index].update(refinement_data)

    if set(refinement_audit_inputs) != {False, True}:
        raise ArithmeticError("both refinement audit cases were not found")
    refinement_permutation_audit: dict[str, dict[str, Any]] = {}
    for discrete, case_name in ((True, "discrete"), (False, "nondiscrete")):
        index, labels, invariant = refinement_audit_inputs[discrete]
        audit = _audit_refinement_permutation_invariance(
            labels, invariant, discrete
        )
        audit["orbit_number"] = index + 1
        audit["shape"] = public_rows[index]["shape"]
        audit["invariant_sha256"] = _bytes_sha256(invariant)
        refinement_permutation_audit[case_name] = audit

    for group in collision_groups:
        if len({refinement_keys[index] for index in group}) != len(group):
            raise ArithmeticError("ordinary color refinement did not separate a residue")
    exact_type_keys = []
    for index, row in enumerate(public_rows):
        degree_key = row.pop("_degree_key")
        refinement = refinement_keys.get(index, b"")
        exact_key = (
            len(degree_key).to_bytes(8, "big")
            + degree_key
            + b"\0refinement\0"
            + refinement
        )
        row["ordinary_Coxeter_type_invariant_sha256"] = _bytes_sha256(exact_key)
        row.setdefault("color_refinement_applied", False)
        exact_type_keys.append(exact_key)
    if len(set(exact_type_keys)) != 2_644:
        raise ArithmeticError("two four-vertex ordinary Coxeter invariants collided")

    facet_type_histogram = Counter(row["facet_count"] for row in public_rows)
    shape_facet_histograms = {
        shape: Counter(
            row["facet_count"] for row in public_rows if row["shape"] == shape
        )
        for shape in ("P4", "K1,3")
    }
    raw_coloring_facet_histogram = Counter()
    for row in public_rows:
        raw_coloring_facet_histogram[row["facet_count"]] += row[
            "unrooted_coloring_orbit_size_under_seed_action"
        ]
    expected_facet_types = Counter({607: 23, 616: 170, 625: 844, 634: 1_607})
    expected_shape_types = {
        "P4": Counter({607: 17, 616: 127, 625: 626, 634: 1_194}),
        "K1,3": Counter({607: 6, 616: 43, 625: 218, 634: 413}),
    }
    expected_raw = Counter({607: 2_480, 616: 47_616, 625: 293_376, 634: 580_640})
    if not (
        facet_type_histogram == expected_facet_types
        and shape_facet_histograms == expected_shape_types
        and raw_coloring_facet_histogram == expected_raw
        and sum(raw_coloring_facet_histogram.values()) == 924_112
        and set(relation_histograms_by_facet) == {607, 616, 625, 634}
        and all(len(values) == 1 for values in relation_histograms_by_facet.values())
    ):
        raise ArithmeticError("the complete four-vertex facet census changed")

    corpus_counts = {
        int(value) for value in bounded_corpus["ordinary_H12_facet_histogram"]
    } | {
        int(value) for value in bounded_corpus["index_two_H12_facet_histogram"]
    }
    new_counts = {607, 616, 625, 634}
    if new_counts & corpus_counts:
        raise ArithmeticError("a four-vertex header appeared in the bounded corpus")

    checks = {
        "frozen_112_wall_factory_replayed": True,
        "exact_384_seed_action_restricted_to_stock": True,
        "all_1379952_admissible_ordered_P4_color_words_covered": True,
        "P4_reversal_quotient_has_696192_unrooted_colorings_and_1964_orbits": True,
        "all_227920_unordered_distinct_K13_color_triples_covered_in_680_orbits": True,
        "all_2644_orbit_representative_boundaries_rebuilt": len(public_rows) == 2_644,
        "every_boundary_has_three_internal_mirrors_and_strict_interior": all(
            row["internal_mirror_count"] == 3 and row["strict_interior_square"] > 0
            for row in public_rows
        ),
        "every_boundary_root_is_primitive_integral_reflective": all(
            row["primitive_integral_reflection_root_count"] == row["facet_count"]
            for row in public_rows
        ),
        "every_boundary_has_fixed_Smith_1_to_13_core": all(
            row["fixed_boundary_core_Smith_diagonal"] == [1] * RANK
            for row in public_rows
        ),
        "every_boundary_pair_is_an_ordinary_Coxeter_relation": True,
        "facet_formula_is_607_plus_9_times_B_edge_occurrences": all(
            row["facet_count"] == row["expected_facet_count_607_plus_9B"]
            for row in public_rows
        ),
        "2535_degree_profiles_leave_exactly_102_groups_and_211_members": True,
        "canonical_color_refinement_separates_every_residual_group": True,
        "canonical_refinement_IDs_and_words_pass_discrete_and_nondiscrete_permutations": all(
            row["canonical_numeric_color_IDs_equivariant"]
            and row["exact_final_invariant_word_unchanged"]
            for row in refinement_permutation_audit.values()
        ),
        "all_2644_ordinary_Coxeter_matrices_are_pairwise_nonisomorphic": len(
            set(exact_type_keys)
        )
        == 2_644,
        "Allcock_index_and_volume_multiplier_are_four": True,
        "all_descendants_are_finite_volume_noncompact": True,
        "new_607_616_625_634_headers_absent_from_pinned_CoxIter_lanes": True,
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"four-vertex census failed: {checks}")

    relation_histogram_table = {
        str(facets): dict(next(iter(values)))
        for facets, values in sorted(relation_histograms_by_facet.items())
    }
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pinned_inputs": {
            "small_subtree_builder_sha256": BASE_SCRIPT_SHA256,
            "small_subtree_archive_sha256": BASE_ARCHIVE_SHA256,
            "full_Picard_carrier_sha256": base.PINNED_SHA256[
                "full_Picard_K3_carrier"
            ],
            "bounded_CoxIter_audit_sha256": base.PINNED_SHA256[
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
            "stock_size": 112,
            "folded_seed_action_order": 384,
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "source_carrier_theorem": carrier["theorem"],
        },
        "coloring_denominator": {
            "P4": {
                "rule": "ordered (a,b,c) with a!=b and b!=c; a=c allowed; quotient by reversal",
                "admissible_ordered_words": 1_379_952,
                "unrooted_colored_paths_before_seed_action": 696_192,
                "seed_action_orbits_after_reversal": len(path_rows),
                "palindromic_outer_color_orbits": sum(
                    row["outer_colors_repeated"] for row in path_rows
                ),
                "unrooted_seed_orbit_size_histogram": {
                    str(key): value
                    for key, value in sorted(
                        Counter(
                            row["unrooted_coloring_orbit_size_under_seed_action"]
                            for row in path_rows
                        ).items()
                    )
                },
            },
            "K1,3": {
                "rule": "unordered triple of three distinct edge colors",
                "unrooted_colored_stars_before_seed_action": 227_920,
                "seed_action_orbits": len(star_rows),
                "seed_orbit_size_histogram": {
                    str(key): value
                    for key, value in sorted(
                        Counter(
                            row["unrooted_coloring_orbit_size_under_seed_action"]
                            for row in star_rows
                        ).items()
                    )
                },
            },
            "total_unrooted_colorings_before_seed_action": 924_112,
            "total_seed_and_tree_automorphism_orbits": len(public_rows),
        },
        "theorem": {
            "four_vertex_tree_shape_count": 2,
            "P4_ordinary_Coxeter_type_count": len(path_rows),
            "K13_ordinary_Coxeter_type_count": len(star_rows),
            "total_pairwise_nonisomorphic_ordinary_Coxeter_polytopes": len(
                public_rows
            ),
            "P4_and_K13_types_are_mutually_disjoint": True,
            "facet_formula": "f(T)=607+9*B(T), where B counts edge-color occurrences in W77--W172",
            "type_facet_histogram": {
                str(key): value for key, value in sorted(facet_type_histogram.items())
            },
            "shape_type_facet_histograms": {
                shape: {
                    str(key): value for key, value in sorted(histogram.items())
                }
                for shape, histogram in shape_facet_histograms.items()
            },
            "raw_unrooted_coloring_facet_histogram": {
                str(key): value
                for key, value in sorted(raw_coloring_facet_histogram.items())
            },
            "ordinary_relation_histogram_by_facet_count": relation_histogram_table,
            "reflection_index_and_volume_multiple": 4,
            "finite_volume": True,
            "compact": False,
            "new_to_project_facet_counts": [607, 616, 625, 634],
        },
        "ordinary_Coxeter_nonisomorphism_certificate": {
            "parallel_and_ultraparallel_both_labelled_infinity": True,
            "initial_exact_degree_profile_class_count": len(degree_groups),
            "degree_profile_collision_group_count": len(collision_groups),
            "degree_profile_collision_member_count": len(collision_indices),
            "maximum_initial_collision_group_size": max(
                map(len, collision_groups)
            ),
            "color_refined_member_count": len(refinement_keys),
            "color_refined_discrete_canonical_form_count": refinement_discrete,
            "all_residual_refinement_invariants_distinct_within_initial_class": True,
            "all_2644_combined_invariants_distinct": True,
            "permutation_invariance_audit": refinement_permutation_audit,
            "method": (
                "canonical one-dimensional color refinement of the complete "
                "ordinary Coxeter matrix; a discrete stable coloring gives a full "
                "canonical adjacency word, while a non-discrete stable quotient is "
                "used only as a separating isomorphism invariant"
            ),
        },
        "allcock_and_carrier": {
            "finite_subtree_theorem": (
                "every representative is a four-copy Coxeter fundamental polytope "
                "of reflection index and volume multiplier four"
            ),
            "noncompactness": (
                "finite-index subgroups of the nonuniform seed lattice are nonuniform"
            ),
            "direct_gates": [
                "exact internal-mirror deletion and primitive hyperplane deduplication",
                "strict timelike summed interior",
                "complete ordinary Coxeter pair matrix",
                "primitive integral reflection roots",
                "fixed W1--W60 boundary core with Smith diagonal 1^13",
            ],
            "same_full_Picard_K3_family": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_effectivity_or_surface_action_claim": False,
        },
        "arithmetic_runtime_bounds": {
            "integer_backend": "numpy int64 after explicit absolute-entry bounds",
            "int64_maximum": INT64_MAX,
            "preflight_absolute_entry_bounds": engine.int64_preflight_bounds,
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_boundary_pairing_absolute": pairing_maximum,
            "all_pre_multiplication_bounds_below_2^63": True,
        },
        "bounded_lineage": {
            "CoxIter_commit": bounded_corpus["commit"],
            "ordinary_H12_file_count": bounded_corpus["ordinary_file_count"],
            "index_two_file_count": bounded_corpus["index_two_file_count"],
            "headers_607_616_625_634_absent": True,
            "evidence_grade": "B1: bounded same-header absence only",
            "worldwide_priority": "open",
            "warning": (
                "Allcock already proved exponential H12 abundance; these exact "
                "surface-derived diagrams have not received a publication-complete audit"
            ),
        },
        "orbit_representatives": public_rows,
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.exists():
            raise SystemExit(f"archive is missing: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from the exact four-vertex rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("P4_types=1964")
    print("K1_3_types=680")
    print("total_types=2644")
    print("new_facet_counts=607,616,625,634")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
