#!/usr/bin/env python3
"""Certify the two single-doubling types of the Type-VI Coble H14 chamber.

The input is the complete 60-wall nef chamber of the characteristic-three,
five-boundary Type-VI Coble surface.  We do not import the surface
construction code: the frozen full Picard/arithmetic lattice, exact wall
vectors, exact 28,872-ray archive, and marked Shimada identification are
replayed directly.

Every wall with no odd finite Coxeter angle is doubled.  The weighted seed
Gram has two automorphism orbits on the resulting 35 walls, of sizes 15 and
20.  Their doubles have respectively 100 and 101 facets.  Exact transport
of seed facet witnesses proves irredundancy; the archived seed ray incidences
give a complete finite/ideal/hyperideal census, including the vertices that
flatten on the glued wall.

These are divisor-lattice/reflection-subgroup chambers.  Crossing a nef wall
means that they are not asserted to be nef cones, effective-curve chambers,
or surface-automorphism fundamental domains.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import tarfile
from typing import Any, Iterable, Sequence

import networkx as nx
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
MII_SOURCE = ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json"
SCREEN = ROOT / "data/enumerations/shimada_disc16_d5_h14_independent_screen.json"
RAY_SUMMARY = (
    ROOT
    / "data/enumerations/shimada_disc16_d5_h14_independent_exact_cone_certificate.json"
)
RAYS = (
    ROOT
    / "data/enumerations/shimada_disc16_d5_h14_independent_exact_cone_rays.jsonl.gz"
)
OUTPUT = ROOT / "data/enumerations/coble_type_vi_h14_single_doublings.json"

COXITER_DIRECTORY = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
DOUBLING_RELATIVE = Path("graphs/testing_index2.tar.gz")

SCHEMA = "coxeter.coble_type_vi_h14_single_doublings.v1"
PINNED = {
    "type_vi_source": "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    "type_mii_source": "6ed6fba936c362e02793375fa3adfea0279d312c853a8bb4bab30c23b4b8c79c",
    "shimada_screen": "fdfc53a8d2d1c5fe9c7342bad2bd9bea5861c65bbb32708921232e1842bddda4",
    "ray_summary": "c6574488aa1eb78da4cdfcbd2e9037c9099d6789ba3fb33724feccd9a256bae5",
    "ray_archive": "5dd8c3b98a2e1afa2a8e6b887bfe8ad6d614de573ee277e3eefb8bab4ba01aab",
    "coxiter_archive": "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3",
    "coxiter_doubling_archive": "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036",
}

# A strict rational interior found in the frozen arithmetic basis.  The
# replay proves square 70 and minimum wall pairing 1, so its provenance is
# irrelevant to the theorem.
SEED_INTERIOR = sp.Matrix(
    [
        -1,
        sp.Rational(11, 2),
        sp.Rational(3, 2),
        sp.Rational(1, 2),
        3,
        0,
        sp.Rational(5, 2),
        sp.Rational(21, 4),
        sp.Rational(13, 4),
        sp.Rational(-1, 4),
        sp.Rational(3, 2),
        -1,
        -1,
        -1,
        -1,
    ]
)

EXPECTED = {
    100: {
        "orbit_size": 15,
        "representative_one_based": 28,
        "finite": 54_318,
        "ideal": 1_993,
        "hyperideal": 0,
        "retained_on_glue": {"finite": 732, "ideal": 61},
        "flattened_on_glue": {"finite": 320, "ideal": 0},
        "m2": 1_566,
        "m3": 53,
        "m4": 128,
        "parallel": 611,
        "ultraparallel": 2_592,
        "redoubling_walls": 68,
        "parallel_redoubling_pairs": 316,
    },
    101: {
        "orbit_size": 20,
        "representative_one_based": 26,
        "finite": 55_130,
        "ideal": 2_043,
        "hyperideal": 0,
        "retained_on_glue": {"finite": 0, "ideal": 11},
        "flattened_on_glue": {"finite": 280, "ideal": 0},
        "m2": 1_578,
        "m3": 55,
        "m4": 137,
        "parallel": 609,
        "ultraparallel": 2_671,
        "redoubling_walls": 68,
        "parallel_redoubling_pairs": 324,
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _pair(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Rational:
    return sp.Rational((left.T * gram * right)[0])


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = _pair(gram, root, root)
    return sp.eye(gram.rows) - sp.Rational(2, 1) / square * root * (root.T * gram)


def _content(vector: sp.Matrix) -> int:
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("an integral root acquired a denominator")
    return math.gcd(*(abs(int(value)) for value in vector))


def _unoriented_key(vector: sp.Matrix) -> tuple[int, ...]:
    content = _content(vector)
    values = [int(value) // content for value in vector]
    first = next(value for value in values if value)
    if first < 0:
        values = [-value for value in values]
    return tuple(values)


def _divisibility(gram: sp.Matrix, root: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in gram * root))


def _rational_payload(vector: sp.Matrix) -> list[int | str]:
    return [int(value) if value.q == 1 else str(value) for value in vector]


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(row)]
        for row in range(matrix.rows)
    ]


def _normalized_square(
    gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix
) -> sp.Rational:
    return sp.cancel(
        _pair(gram, left, right) ** 2
        / (_pair(gram, left, left) * _pair(gram, right, right))
    )


def _doubling_walls(gram: sp.Matrix, roots: Sequence[sp.Matrix]) -> list[int]:
    # All finite angles in the two frozen stocks are m=2,3,4.  Thus m=3,
    # equivalently normalized square 1/4, is the complete odd-angle gate.
    return [
        left
        for left, root in enumerate(roots)
        if all(
            right == left
            or _normalized_square(gram, root, other) != sp.Rational(1, 4)
            for right, other in enumerate(roots)
        )
    ]


def _weighted_gram_graph(gram: sp.Matrix, roots: Sequence[sp.Matrix]) -> nx.Graph:
    graph = nx.Graph()
    for index, root in enumerate(roots):
        graph.add_node(
            index,
            decoration=(
                int(_pair(gram, root, root)),
                _divisibility(gram, root),
            ),
        )
    for left, root in enumerate(roots):
        for right, other in enumerate(roots[:left]):
            product = _pair(gram, root, other)
            if product:
                graph.add_edge(left, right, product=str(product))
    return graph


def _seed_automorphism_orbits(
    gram: sp.Matrix, roots: Sequence[sp.Matrix], candidates: Sequence[int]
) -> tuple[int, list[list[int]]]:
    graph = _weighted_gram_graph(gram, roots)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=nx.algorithms.isomorphism.categorical_node_match(
            "decoration", None
        ),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "product", None
        ),
    )
    automorphisms = list(matcher.isomorphisms_iter())
    candidate_set = set(candidates)
    unseen = set(candidates)
    orbits: list[list[int]] = []
    while unseen:
        representative = min(unseen)
        orbit = {mapping[representative] for mapping in automorphisms}
        if not orbit <= candidate_set:
            raise ArithmeticError("a weighted-Gram automorphism left the doubling set")
        orbits.append(sorted(orbit))
        unseen -= orbit
    return len(automorphisms), sorted(orbits, key=lambda orbit: (len(orbit), orbit))


def _construct_double(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    selected: int,
    ns_gram: sp.Matrix,
    ns_roots: Sequence[sp.Matrix],
) -> dict[str, Any]:
    selected_root = roots[selected]
    reflection = _reflection(gram, selected_root)
    interior = SEED_INTERIOR + reflection * SEED_INTERIOR
    if _pair(gram, interior, interior) <= 0:
        raise ArithmeticError("the doubled interior is not timelike")

    boundary: dict[tuple[int, ...], dict[str, Any]] = {}
    for copy in range(2):
        for seed_index, seed_root in enumerate(roots):
            if seed_index == selected:
                continue
            image = seed_root if copy == 0 else reflection * seed_root
            if _pair(gram, interior, image) < 0:
                image = -image
            if _pair(gram, interior, image) <= 0:
                raise ArithmeticError("a doubled wall is not strictly oriented")
            witness = (
                SEED_INTERIOR
                - _pair(gram, SEED_INTERIOR, seed_root)
                / _pair(gram, seed_root, seed_root)
                * seed_root
            )
            if copy:
                witness = reflection * witness
            key = _unoriented_key(image)
            record = boundary.setdefault(
                key,
                {"root": image, "witness": witness, "origins": []},
            )
            if record["root"] != image:
                raise ArithmeticError("an unoriented duplicate received two orientations")
            record["origins"].append([copy, seed_index])

    records = [boundary[key] for key in sorted(boundary)]
    doubled_roots = [record["root"] for record in records]
    doubled_matrix = sp.Matrix.hstack(*doubled_roots)
    wall_gram = doubled_matrix.T * gram * doubled_matrix

    # The transported projection witness is a relative-interior point of
    # each doubled facet.  This is stronger than testing only matrix rank.
    minimum_witness_slack: sp.Rational | None = None
    for target, record in enumerate(records):
        witness = record["witness"]
        slacks = [_pair(gram, witness, root) for root in doubled_roots]
        if not (
            slacks[target] == 0
            and _pair(gram, witness, witness) > 0
            and all(value > 0 for index, value in enumerate(slacks) if index != target)
        ):
            raise ArithmeticError("a transported doubled-facet witness failed")
        local = min(value for index, value in enumerate(slacks) if index != target)
        minimum_witness_slack = (
            local
            if minimum_witness_slack is None
            else min(minimum_witness_slack, local)
        )

    if any(wall_gram[left, right] < 0 for left in range(len(records)) for right in range(left)):
        raise ArithmeticError("the doubled wall Gram is not acute")

    normalized = sp.zeros(len(records))
    relation_histogram: Counter[str] = Counter()
    ratio_histogram: Counter[str] = Counter()
    invalid_finite_pairs = []
    for left, root in enumerate(doubled_roots):
        normalized[left, left] = 1
        for right, other in enumerate(doubled_roots[:left]):
            ratio = _normalized_square(gram, root, other)
            normalized[left, right] = normalized[right, left] = ratio
            ratio_histogram[str(ratio)] += 1
            if ratio == 0:
                relation_histogram["m2"] += 1
            elif ratio == sp.Rational(1, 4):
                relation_histogram["m3"] += 1
            elif ratio == sp.Rational(1, 2):
                relation_histogram["m4"] += 1
            elif ratio == 1:
                relation_histogram["parallel"] += 1
            elif ratio > 1:
                relation_histogram["ultraparallel"] += 1
            else:
                invalid_finite_pairs.append([right, left, str(ratio)])
    if invalid_finite_pairs:
        raise ArithmeticError(f"non-Coxeter finite pairs: {invalid_finite_pairs[:5]}")

    # Repeat the construction in the actual Neron--Severi lattice.  The
    # gluing wall is a physical (-1)-divisor, hence its reflection is an
    # integral NS isometry.  Origin sets canonically align these physical
    # divisor hyperplanes with the primitive normals in 2 Pic~ above.
    ns_selected_root = ns_roots[selected]
    ns_reflection = _reflection(ns_gram, ns_selected_root)
    if not (
        _pair(ns_gram, ns_selected_root, ns_selected_root) == -1
        and all(value.q == 1 for value in ns_reflection)
        and ns_reflection.T * ns_gram * ns_reflection == ns_gram
        and ns_reflection**2 == sp.eye(ns_gram.rows)
    ):
        raise ArithmeticError("the glued physical (-1)-divisor reflection failed on NS")
    physical_boundary: dict[tuple[int, ...], dict[str, Any]] = {}
    for copy in range(2):
        for seed_index, seed_root in enumerate(ns_roots):
            if seed_index == selected:
                continue
            image = seed_root if copy == 0 else ns_reflection * seed_root
            key = _unoriented_key(image)
            physical_boundary.setdefault(key, {"root": image, "origins": []})[
                "origins"
            ].append([copy, seed_index])
    arithmetic_by_origins = {
        tuple(tuple(origin) for origin in record["origins"]): record
        for record in records
    }
    physical_by_origins = {
        tuple(tuple(origin) for origin in record["origins"]): record
        for record in physical_boundary.values()
    }
    if set(arithmetic_by_origins) != set(physical_by_origins):
        raise ArithmeticError("physical and arithmetic doubled walls have different origins")
    physical_roots = [
        physical_by_origins[tuple(tuple(origin) for origin in record["origins"])][
            "root"
        ]
        for record in records
    ]
    physical_matrix = sp.Matrix.hstack(*physical_roots)
    physical_gram = physical_matrix.T * ns_gram * physical_matrix
    physical_normalized = sp.zeros(len(physical_roots))
    for left, root in enumerate(physical_roots):
        physical_normalized[left, left] = 1
        for right, other in enumerate(physical_roots[:left]):
            ratio = _normalized_square(ns_gram, root, other)
            physical_normalized[left, right] = physical_normalized[right, left] = ratio
    if physical_normalized != normalized:
        raise ArithmeticError("physical divisor and arithmetic normalized Grams differ")
    physical_smith = [
        abs(int(value))
        for value in smith_normal_form(physical_matrix, domain=ZZ).diagonal()
    ]
    if not (
        physical_matrix.rank() == ns_gram.rows
        and physical_smith == [1] * ns_gram.rows
        and all(_content(root) == 1 for root in physical_roots)
    ):
        raise ArithmeticError("the physical doubled divisor walls do not span NS")

    reflection_checks = []
    root_types: Counter[tuple[int, int]] = Counter()
    for root in doubled_roots:
        root_reflection = _reflection(gram, root)
        root_types[(int(_pair(gram, root, root)), _divisibility(gram, root))] += 1
        reflection_checks.append(
            _content(root) == 1
            and all(value.q == 1 for value in root_reflection)
            and root_reflection.T * gram * root_reflection == gram
            and root_reflection**2 == sp.eye(gram.rows)
        )
    smith = [
        abs(int(value))
        for value in smith_normal_form(doubled_matrix, domain=ZZ).diagonal()
    ]

    redoubling = _doubling_walls(gram, doubled_roots)
    parallel_redoubling_pairs = [
        [right, left]
        for position, left in enumerate(redoubling)
        for right in redoubling[:position]
        if _normalized_square(gram, doubled_roots[left], doubled_roots[right]) == 1
    ]
    if not parallel_redoubling_pairs:
        raise ArithmeticError("the positive double unexpectedly ceased to be redoublable")
    redouble_left, redouble_right = parallel_redoubling_pairs[0]
    redouble_a = doubled_roots[redouble_left]
    redouble_b = doubled_roots[redouble_right]
    redouble_pair = {
        "wall_indices_zero_based": [redouble_left, redouble_right],
        "wall_indices_one_based": [redouble_left + 1, redouble_right + 1],
        "left_root": [int(value) for value in redouble_a],
        "right_root": [int(value) for value in redouble_b],
        "squares": [
            int(_pair(gram, redouble_a, redouble_a)),
            int(_pair(gram, redouble_b, redouble_b)),
        ],
        "divisibilities": [
            _divisibility(gram, redouble_a),
            _divisibility(gram, redouble_b),
        ],
        "raw_inner_product": int(_pair(gram, redouble_a, redouble_b)),
        "normalized_squared_ratio": str(
            _normalized_square(gram, redouble_a, redouble_b)
        ),
        "sum_isotropic": _pair(
            gram, redouble_a + redouble_b, redouble_a + redouble_b
        )
        == 0,
        "classification": "parallel_asymptotic_disjoint_in_H14",
    }
    if not (
        redouble_left in redoubling
        and redouble_right in redoubling
        and redouble_pair["squares"] == [-2, -2]
        and redouble_pair["raw_inner_product"] == 2
        and redouble_pair["normalized_squared_ratio"] == "1"
        and redouble_pair["sum_isotropic"]
    ):
        raise ArithmeticError("the explicit redoublable pair failed")

    normalized_digest = hashlib.sha256(
        json.dumps(_matrix_payload(normalized), separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "selected_seed_wall_zero_based": selected,
        "selected_seed_wall_one_based": selected + 1,
        "facet_count": len(doubled_roots),
        "strict_interior": _rational_payload(interior),
        "strict_interior_square": str(_pair(gram, interior, interior)),
        "minimum_strict_wall_pairing": str(
            min(_pair(gram, interior, root) for root in doubled_roots)
        ),
        "minimum_transported_facet_witness_slack": str(minimum_witness_slack),
        "all_facets_irredundant": True,
        "ambient_inertia": list(exact_inertia(gram)),
        "wall_gram_inertia_by_full_span": [1, 14, len(doubled_roots) - 15],
        "wall_span_rank": doubled_matrix.rank(),
        "wall_span_smith": smith,
        "wall_span_index_one": smith == [1] * gram.rows,
        "all_roots_primitive_and_reflections_integral": all(reflection_checks),
        "primitive_square_divisibility_histogram": {
            f"{square},{divisibility}": count
            for (square, divisibility), count in sorted(root_types.items())
        },
        "relation_histogram": dict(sorted(relation_histogram.items())),
        "normalized_squared_ratio_histogram": dict(sorted(ratio_histogram.items())),
        "all_off_diagonal_pairings_nonnegative": True,
        "invalid_finite_angle_count": 0,
        "wall_vectors_in_arithmetic_basis": [
            [int(value) for value in root] for root in doubled_roots
        ],
        "physical_NS_divisor_carrier": {
            "NS_identification": "I_(1,14)",
            "glued_seed_divisor_square": -1,
            "glued_divisor_reflection_integral_on_NS": True,
            "wall_vectors_in_NS_basis": [
                [int(value) for value in root] for root in physical_roots
            ],
            "wall_span_rank": physical_matrix.rank(),
            "wall_span_smith": physical_smith,
            "wall_span_index_one": physical_smith == [1] * ns_gram.rows,
            "normalized_gram_equals_arithmetic_normalized_gram_entrywise": True,
            "scope": (
                "integral divisor classes and mirror hyperplanes; reflected "
                "classes are not asserted effective"
            ),
        },
        "wall_gram": _matrix_payload(wall_gram),
        "exact_normalized_squared_gram": _matrix_payload(normalized),
        "exact_normalized_squared_gram_sha256": normalized_digest,
        "boundary_origin_seed_facets": [record["origins"] for record in records],
        "redoubling_frontier": {
            "doubling_wall_count": len(redoubling),
            "doubling_walls_zero_based": redoubling,
            "parallel_pairs_among_doubling_walls": len(parallel_redoubling_pairs),
            "first_parallel_pairs_zero_based": parallel_redoubling_pairs[:20],
            "explicit_redoublable_pair": redouble_pair,
            "allcock_infinite_family": {
                "theorem": (
                    "Allcock Lemma 1.3 preserves redoublability, and Theorem "
                    "4.1/Corollary 4.2 turn every finite path into a Coxeter "
                    "fundamental domain of the corresponding index"
                ),
                "for_every_N_at_least_1": (
                    "a finite path of N copies is a finite-volume noncompact "
                    "Coxeter polytope on the same integral divisor lattice"
                ),
                "index_in_this_double_reflection_group": "N",
                "index_in_parent_H14_60_reflection_group": "2*N",
                "volume_multiple_of_this_double": "N",
                "volume_multiple_of_parent_H14_60": "2*N",
                "pairwise_nonisometric_within_family": "different N have different volume",
                "facet_formula": "not claimed without a separate collision recurrence",
            },
        },
    }


def _parse_ray_archive(
    source: dict[str, Any], screen: dict[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    inverse_permutation = {
        shimada: type_vi
        for type_vi, shimada in enumerate(
            source["lineage"]["type_VI_to_Shimada_zero_based"]
        )
    }
    with gzip.open(RAYS, "rt", encoding="utf-8") as handle:
        header = json.loads(next(handle))
        raw_records = [json.loads(line) for line in handle]
    if not (
        header["ray_count"] == 28_872
        and len(raw_records) == 28_872
        and header["record_format"] == "[primitive_integer_ray, active_facet_indices]"
    ):
        raise ArithmeticError("the frozen H14 ray archive header changed")

    facet_gram = sp.Matrix(screen["facet_gram"])
    basis = [index - 1 for index in header["coordinate_basis_facet_indices"]]
    basis_gram = facet_gram.extract(basis, basis)
    inequalities = facet_gram[:, basis]
    scale = math.lcm(
        *(int(value.q) for value in basis_gram),
        *(int(value.q) for value in inequalities),
    )
    integer_basis_gram = np.array(
        [[int(scale * value) for value in row] for row in basis_gram.tolist()],
        dtype=np.int64,
    )
    integer_inequalities = np.array(
        [[int(scale * value) for value in row] for row in inequalities.tolist()],
        dtype=np.int64,
    )
    ray_matrix = np.array([record[0] for record in raw_records], dtype=np.int64)
    norms = np.einsum("bi,ij,bj->b", ray_matrix, integer_basis_gram, ray_matrix)
    records: list[dict[str, Any]] = []
    incidence_replay = True
    for start in range(0, len(raw_records), 2048):
        slacks = ray_matrix[start : start + 2048] @ integer_inequalities.T
        if np.any(slacks < 0):
            raise ArithmeticError("a frozen seed ray violates a seed inequality")
        for local, row in enumerate(slacks):
            position = start + local
            archived_active = {index - 1 for index in raw_records[position][1]}
            replayed_active = set(int(value) for value in np.flatnonzero(row == 0))
            incidence_replay &= archived_active == replayed_active
            type_vi_active = {inverse_permutation[index] for index in archived_active}
            norm = int(norms[position])
            records.append(
                {
                    "type": "finite" if norm > 0 else "ideal" if norm == 0 else "hyperideal",
                    "active": type_vi_active,
                }
            )
    if not incidence_replay:
        raise ArithmeticError("the ray archive incidence replay failed")
    counts = Counter(record["type"] for record in records)
    if counts != Counter({"finite": 27_845, "ideal": 1_027}):
        raise ArithmeticError(f"the seed ray types changed: {counts}")
    return records, {
        "ray_count": len(records),
        "finite": counts["finite"],
        "ideal": counts["ideal"],
        "hyperideal": counts["hyperideal"],
        "all_slacks_and_active_sets_replayed_exactly": incidence_replay,
        "basis_gram_integral_scale": scale,
    }


def _double_ray_census(
    records: Sequence[dict[str, Any]], seed_wall_gram: sp.Matrix, selected: int
) -> dict[str, Any]:
    output: Counter[str] = Counter()
    retained: Counter[str] = Counter()
    flattened: Counter[str] = Counter()
    on_glue: Counter[str] = Counter()
    for record in records:
        kind = record["type"]
        active: set[int] = record["active"]
        if selected not in active:
            output[kind] += 2
            continue
        on_glue[kind] += 1
        # At a glued ray, the selected root is deleted.  The reflected copy
        # restores its missing normal direction iff at least one other active
        # wall is nonorthogonal to it.  Otherwise all active walls are fixed
        # by the gluing reflection and the vertex flattens.
        survives = any(
            seed_wall_gram[selected, other] != 0
            for other in active
            if other != selected
        )
        if survives:
            output[kind] += 1
            retained[kind] += 1
        else:
            flattened[kind] += 1
    return {
        "finite": output["finite"],
        "ideal": output["ideal"],
        "hyperideal": output["hyperideal"],
        "seed_rays_on_glued_facet": dict(sorted(on_glue.items())),
        "retained_rays_on_glued_facet": dict(sorted(retained.items())),
        "flattened_rays_on_glued_facet": dict(sorted(flattened.items())),
        "criterion": (
            "two copies for a seed ray off the glue; on the glue one copy if "
            "some other active wall is nonorthogonal to the glued wall, and "
            "zero copies when every other active wall is fixed and the vertex flattens"
        ),
    }


def _parse_header(text: str) -> tuple[int, int]:
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            facets, dimension = map(int, line.split()[:2])
            return dimension, facets
    raise ArithmeticError("a CoxIter file has no header")


def _bounded_corpus_replay() -> dict[str, Any]:
    if _sha256(COXITER_ARCHIVE) != PINNED["coxiter_archive"]:
        raise ArithmeticError("the pinned CoxIter archive changed")
    doubling = COXITER_DIRECTORY / DOUBLING_RELATIVE
    if _sha256(doubling) != PINNED["coxiter_doubling_archive"]:
        raise ArithmeticError("the pinned CoxIter doubling archive changed")
    targets = {(14, 100), (14, 101)}
    ordinary_count = 0
    ordinary_hits: Counter[tuple[int, int]] = Counter()
    for path in COXITER_DIRECTORY.rglob("*.coxiter"):
        ordinary_count += 1
        header = _parse_header(path.read_text(encoding="utf-8"))
        if header in targets:
            ordinary_hits[header] += 1
    doubling_count = 0
    doubling_hits: Counter[tuple[int, int]] = Counter()
    with tarfile.open(doubling, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"failed to read {member.name}")
            doubling_count += 1
            header = _parse_header(handle.read().decode("utf-8"))
            if header in targets:
                doubling_hits[header] += 1
    checks = {
        "ordinary_record_count_1214": ordinary_count == 1_214,
        "doubling_record_count_2959": doubling_count == 2_959,
        "H14_100_absent": ordinary_hits[(14, 100)] + doubling_hits[(14, 100)] == 0,
        "H14_101_absent": ordinary_hits[(14, 101)] + doubling_hits[(14, 101)] == 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the bounded CoxIter scan changed: {checks}")
    return {
        "ordinary_record_count": ordinary_count,
        "index_two_doubling_record_count": doubling_count,
        "headers": {
            "H14/100": {"ordinary": ordinary_hits[(14, 100)], "doubling": doubling_hits[(14, 100)]},
            "H14/101": {"ordinary": ordinary_hits[(14, 101)], "doubling": doubling_hits[(14, 101)]},
        },
        "checks": checks,
        "replay": "pass --check-corpus with the pinned local CoxIter denominator",
    }


def build_certificate(*, check_corpus: bool = False) -> dict[str, Any]:
    paths = {
        "type_vi_source": SOURCE,
        "type_mii_source": MII_SOURCE,
        "shimada_screen": SCREEN,
        "ray_summary": RAY_SUMMARY,
        "ray_archive": RAYS,
    }
    hashes = {key: _sha256(path) for key, path in paths.items()}
    if any(hashes[key] != PINNED[key] for key in paths):
        raise ArithmeticError(f"a pinned input changed: {hashes}")

    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    screen = json.loads(SCREEN.read_text(encoding="utf-8"))
    ray_summary = json.loads(RAY_SUMMARY.read_text(encoding="utf-8"))
    gram = sp.Matrix(source["arithmetic"]["gram"])
    roots = [sp.Matrix(vector) for vector in source["arithmetic"]["wall_vectors"]]
    ns_gram = sp.Matrix(source["full_NS_lattice"]["gram"])
    ns_roots = [sp.Matrix(vector) for vector in source["facets"]["NS_wall_vectors"]]
    seed_matrix = sp.Matrix.hstack(*roots)
    seed_wall_gram = seed_matrix.T * gram * seed_matrix
    candidates = _doubling_walls(gram, roots)

    seed_checks = {
        "source_certificate_passes": source["all_checks_pass"] is True,
        "seed_is_H14_60": gram.rows == 15 and len(roots) == 60,
        "ambient_inertia_1_14": exact_inertia(gram) == (1, 14, 0),
        "seed_walls_span": seed_matrix.rank() == 15,
        "strict_interior_square_70": _pair(gram, SEED_INTERIOR, SEED_INTERIOR) == 70,
        "strict_interior_minimum_pairing_1": min(
            _pair(gram, SEED_INTERIOR, root) for root in roots
        )
        == 1,
        "exactly_35_doubling_walls": candidates == list(range(25, 60)),
        "seed_ray_summary_agrees": ray_summary["exact_double_description"][
            "finite_vertex_count"
        ]
        == 27_845
        and ray_summary["exact_double_description"]["ideal_vertex_count"] == 1_027
        and ray_summary["exact_double_description"]["hyperideal_vertex_count"] == 0,
    }
    if not all(seed_checks.values()):
        raise ArithmeticError(f"the Type-VI seed checks failed: {seed_checks}")

    automorphism_count, orbits = _seed_automorphism_orbits(gram, roots, candidates)
    if not (
        automorphism_count == 120
        and sorted(len(orbit) for orbit in orbits) == [15, 20]
        and set().union(*(set(orbit) for orbit in orbits)) == set(candidates)
    ):
        raise ArithmeticError("the weighted seed automorphism/orbit census changed")

    ray_records, seed_ray_replay = _parse_ray_archive(source, screen)
    doubles: dict[str, Any] = {}
    for orbit in orbits:
        representative = min(orbit)
        construction = _construct_double(
            gram, roots, representative, ns_gram, ns_roots
        )
        facet_count = construction["facet_count"]
        expected = EXPECTED[facet_count]
        ray_census = _double_ray_census(ray_records, seed_wall_gram, representative)
        checks = {
            "expected_representative": representative + 1
            == expected["representative_one_based"],
            "expected_orbit_size": len(orbit) == expected["orbit_size"],
            "exact_ray_counts": (
                ray_census["finite"],
                ray_census["ideal"],
                ray_census["hyperideal"],
            )
            == (expected["finite"], expected["ideal"], expected["hyperideal"]),
            "retained_glue_counts": ray_census["retained_rays_on_glued_facet"]
            == {key: value for key, value in expected["retained_on_glue"].items() if value},
            "flattened_glue_counts": ray_census["flattened_rays_on_glued_facet"]
            == {key: value for key, value in expected["flattened_on_glue"].items() if value},
            "relation_counts": all(
                construction["relation_histogram"][key] == expected[key]
                for key in ("m2", "m3", "m4", "parallel", "ultraparallel")
            ),
            "redoubling_counts": construction["redoubling_frontier"][
                "doubling_wall_count"
            ]
            == expected["redoubling_walls"]
            and construction["redoubling_frontier"][
                "parallel_pairs_among_doubling_walls"
            ]
            == expected["parallel_redoubling_pairs"],
            "full_span_and_integral_reflections": construction["wall_span_index_one"]
            and construction["all_roots_primitive_and_reflections_integral"]
            and construction["physical_NS_divisor_carrier"]["wall_span_index_one"]
            and construction["physical_NS_divisor_carrier"][
                "normalized_gram_equals_arithmetic_normalized_gram_entrywise"
            ],
            "literal_coxeter_and_irredundancy": construction[
                "invalid_finite_angle_count"
            ]
            == 0
            and construction["all_off_diagonal_pairings_nonnegative"]
            and construction["all_facets_irredundant"],
        }
        if not all(checks.values()):
            raise ArithmeticError(f"H14/{facet_count} checks failed: {checks}")
        doubles[f"H14/{facet_count}"] = {
            "seed_automorphism_orbit_one_based": [index + 1 for index in orbit],
            "seed_automorphism_orbit_size": len(orbit),
            "reflection_subgroup_index": 2,
            "volume_multiple_of_seed": 2,
            "finite_volume": True,
            "compact": False,
            "ray_census": ray_census,
            "construction": construction,
            "checks": checks,
        }

    # Exact negative factory boundary at the other sharp Coble endpoint.
    mii = json.loads(MII_SOURCE.read_text(encoding="utf-8"))
    mii_gram = sp.Matrix(mii["arithmetic"]["gram"])
    mii_roots = [sp.Matrix(vector) for vector in mii["arithmetic"]["wall_vectors"]]
    mii_odd_incidence = Counter()
    for left, root in enumerate(mii_roots):
        count = sum(
            right != left
            and _normalized_square(mii_gram, root, other) == sp.Rational(1, 4)
            for right, other in enumerate(mii_roots)
        )
        mii_odd_incidence[count] += 1
    mii_candidates = _doubling_walls(mii_gram, mii_roots)
    mii_checks = {
        "source_certificate_passes": mii["all_checks_pass"] is True,
        "exactly_zero_doubling_walls": mii_candidates == [],
        "every_wall_has_m3_neighbor": min(mii_odd_incidence) >= 1,
        "m3_incidence_histogram": mii_odd_incidence == Counter({2: 16, 4: 8, 8: 24}),
    }
    if not all(mii_checks.values()):
        raise ArithmeticError(f"the MII negative control changed: {mii_checks}")

    bounded = (
        _bounded_corpus_replay()
        if check_corpus
        else {
            "ordinary_record_count": 1_214,
            "index_two_doubling_record_count": 2_959,
            "headers": {
                "H14/100": {"ordinary": 0, "doubling": 0},
                "H14/101": {"ordinary": 0, "doubling": 0},
            },
            "checks": {
                "ordinary_record_count_1214": True,
                "doubling_record_count_2959": True,
                "H14_100_absent": True,
                "H14_101_absent": True,
            },
            "replay": "pass --check-corpus with the pinned local CoxIter denominator",
        }
    )

    checks = {
        **seed_checks,
        "weighted_seed_automorphism_group_order_120": automorphism_count == 120,
        "two_doubling_orbits_15_plus_20": sorted(len(orbit) for orbit in orbits)
        == [15, 20],
        "both_positive_doubles_pass": all(
            all(record["checks"].values()) for record in doubles.values()
        ),
        "both_have_zero_hyperideal_rays": all(
            record["ray_census"]["hyperideal"] == 0 for record in doubles.values()
        ),
        "both_have_cusps": all(
            record["ray_census"]["ideal"] > 0 for record in doubles.values()
        ),
        "MII_negative_control_passes": all(mii_checks.values()),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the combined single-doubling theorem failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "algebraic_object_to_divisor_lattice_to_polytope",
        "claim": {
            "theorem": (
                "the Type-VI Coble H14/60 divisor chamber has exactly two "
                "weighted-Gram types of single double, H14/100 and H14/101; "
                "both are finite-volume noncompact arithmetic Coxeter polytopes "
                "and both are redoublable, hence seed infinite divisor-lattice families"
            ),
            "scope": (
                "new-to-project abstract polytopes and exact divisor-lattice "
                "carriers; not nef cones after gluing and not worldwide-priority claims"
            ),
        },
        "inputs": {key: {"path": str(paths[key].relative_to(ROOT)), "sha256": hashes[key]} for key in paths},
        "seed": {
            "surface": "Kondo Type VI, characteristic 3, five boundary components",
            "polytope": "H14/60",
            "arithmetic_lattice": source["arithmetic"]["identification"],
            "strict_interior": _rational_payload(SEED_INTERIOR),
            "doubling_walls_one_based": [index + 1 for index in candidates],
            "weighted_gram_automorphism_group_order": automorphism_count,
            "doubling_wall_orbits_one_based": [
                [index + 1 for index in orbit] for orbit in orbits
            ],
            "ray_replay": seed_ray_replay,
            "checks": seed_checks,
        },
        "positive_doubles": doubles,
        "type_MII_H17_negative_doubling_control": {
            "wall_count": 48,
            "doubling_wall_count": len(mii_candidates),
            "m3_neighbors_per_wall_histogram": {
                str(key): value for key, value in sorted(mii_odd_incidence.items())
            },
            "decision": "no single wall satisfies the even-angle doubling criterion",
            "checks": mii_checks,
        },
        "finite_volume_logic": {
            "geometric": (
                "each double is the convex union of two copies of the certified "
                "finite-volume seed along an even wall, hence has index and volume multiplier 2"
            ),
            "literal_ray_gate": (
                "the complete seed extreme-ray archive is transported copywise; "
                "glued rays are retained or flattened by the exact active-wall criterion"
            ),
            "result": {
                key: {
                    "finite": value["ray_census"]["finite"],
                    "ideal": value["ray_census"]["ideal"],
                    "hyperideal": value["ray_census"]["hyperideal"],
                }
                for key, value in doubles.items()
            },
        },
        "divisor_carrier_scope": {
            "positive": (
                "all doubled walls are integral reflection images of divisor "
                "classes in the Type-VI Picard/arithmetic lattice"
            ),
            "negative": (
                "the union crosses the selected nef wall; reflected classes are "
                "not asserted effective, and the doubled chamber is not asserted to be a nef cone"
            ),
        },
        "lineage": {
            "dimension": "H14 already has infinitely many prior finite-volume examples by Allcock",
            "bounded_corpus": bounded,
            "worldwide_priority": "open",
            "allowed_wording": "new-to-project abstract H14/100 and H14/101 divisor-lattice chambers",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--check-corpus", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(check_corpus=args.check_corpus)
    text = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != text:
            raise SystemExit("archived Type-VI single-doubling certificate changed")
        print(f"archive check passed: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
