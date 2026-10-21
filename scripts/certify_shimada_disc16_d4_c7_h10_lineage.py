#!/usr/bin/env python3
"""Certify the exact lattice and reflection lineage of ``D0004_C007``.

The archived Shimada discriminant-16 fixed section has 27 Coxeter facets in
H10.  This replay keeps two integral lattices separate:

* the saturated divisor-fixed lattice ``L`` has determinant 512, and the 27
  primitive mirrors generate ``L``;
* the same reflections preserve a unique nontrivial integral overlattice
  ``M`` of index 4 and determinant 32.

An explicit isometry identifies ``M`` with
``U + D6(-1) + A1(-1)^3`` and with Scharlau--Kirschmer List10/Gram10 entry
20.  A fresh exact Vinberg enumeration recovers its 14-wall fundamental
reflection chamber.  Solomon's formula then gives the exact index of the
27-wall subgroup.

The CoxIter comparison is deliberately bounded and is not a worldwide
priority claim for the literal 27-wall chamber.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import re
import tarfile
from typing import Any, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from coxeter_polytopes.search.disc16_d4c007_maximal_vinberg import (
    divisibility,
    height_zero_simple_roots,
    lattice_gram,
    pair,
    run_search,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _factorization,
    _relation_matrix,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    _connected_elliptic_masks,
    _parse_coxiter,
    _parse_coxiter_text,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CASE = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_atlas_archived_high_sections"
    / "D0004_C007.json"
)
CONE = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_exact_cone_certificate.json"
)
VINBERG = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_vinberg_face_certificate.json"
)
WEIGHTED_INVENTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_weighted_wall_symmetry_inventory.json"
)
WEIGHTED_INVENTORY_SHA256 = (
    "e2b1da1b0fbb8bac659dcb77bd59a640df17d95996e066c5065808ded9160fb8"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_lineage.json"
)

COXITER_CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_RELATIVE = Path("graphs/testing_index2.tar.gz")
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)

SK_DIRECTORY = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
SK_ARCHIVE = Path("/private/tmp/Scharlau-Kirschmer-Lorentz.zip")
SK_ARCHIVE_SHA256 = (
    "74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5"
)
SK_LIST10_SHA256 = (
    "9e627c4e1303015fd96acdcf6ccbeb264ebadba49abfcd761705f54d211b2e21"
)
SK_GRAM10_SHA256 = (
    "b8186e1fa2771c517a6040998f3b234568c3edde3ec8b7bb8412865d7710af18"
)
SK_GENUS = "II_{10,1}([2^{5}]_{1})"


# Columns give a basis of M in the archived basis of L.  Thus M=B Z^11 inside
# L tensor Q, while L=Z^11 is the index-four sublattice; det(B)=1/4 records
# [M:L]=4.
OVERLATTICE_BASIS = sp.eye(11)
OVERLATTICE_BASIS[4, 5] = sp.Rational(1, 2)
OVERLATTICE_BASIS[5, 5] = sp.Rational(1, 2)
OVERLATTICE_BASIS[8, 10] = sp.Rational(1, 2)
OVERLATTICE_BASIS[10, 10] = sp.Rational(1, 2)


# Q^T G_M Q is the standard U + D6(-1) + A1(-1)^3 Gram matrix.
OVERLATTICE_TO_STANDARD = sp.Matrix(
    [
        [10, 100, -20, 0, 10, 31, -11, -51, -6, 0, 0],
        [13, 130, -27, 0, 13, 40, -14, -66, -8, 0, 0],
        [16, 160, -34, 1, 16, 49, -17, -81, -10, 0, 0],
        [12, 120, -24, -4, 16, 36, -12, -60, -6, 1, 1],
        [7, 70, -15, 2, 5, 22, -7, -36, -5, 0, 0],
        [8, 80, -16, -3, 12, 23, -8, -40, -4, 0, 0],
        [12, 120, -24, -4, 16, 36, -12, -60, -6, 1, 0],
        [2, 20, -4, -1, 3, 6, -2, -10, -1, 0, 0],
        [-8, -79, 16, 2, -10, -24, 8, 40, 4, 0, 0],
        [-6, -60, 12, 2, -8, -18, 6, 30, 3, 0, 0],
        [5, 50, -10, -2, 7, 15, -5, -25, -2, 0, 0],
    ]
)


# S^T Gram10[20] S is the negative of the standard Lorentz Gram.  The sign
# difference is exactly the Scharlau--Kirschmer signature convention (10,1).
STANDARD_TO_SK_SIGNED = sp.Matrix(
    [
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0],
        [0, 0, 0, 0, -1, 1, -1, 0, 0, 0, 0],
        [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _rational(value: sp.Rational | int) -> int | str:
    rational = sp.Rational(value)
    return int(rational) if rational.q == 1 else str(rational)


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [[_rational(value) for value in row] for row in matrix.tolist()]


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(matrix.rows, matrix.cols))
        if smith[index, index]
    ]


def _root_type(gram: sp.Matrix, root: sp.Matrix) -> tuple[int, int]:
    square = int((root.T * gram * root)[0])
    root_divisibility = math.gcd(*(abs(int(value)) for value in gram * root))
    return square, root_divisibility


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.eye(gram.rows) - sp.Rational(2, square) * root * (root.T * gram)


def _integral_reflection(gram: sp.Matrix, root: sp.Matrix) -> bool:
    reflection = _reflection(gram, root)
    return (
        all(value.is_Integer is True for value in reflection)
        and reflection.T * gram * reflection == gram
        and reflection**2 == sp.eye(gram.rows)
        and reflection * root == -root
    )


def _standard_gram() -> sp.Matrix:
    gram = sp.zeros(11)
    gram[0, 1] = gram[1, 0] = 1
    for index in range(6):
        gram[2 + index, 2 + index] = -2
    for left, right in ((0, 1), (1, 2), (2, 3), (3, 4), (3, 5)):
        gram[2 + left, 2 + right] = gram[2 + right, 2 + left] = 1
    for index in range(8, 11):
        gram[index, index] = -2
    return gram


def _parse_sk_list10(path: Path) -> tuple[int, dict[str, Any]]:
    pattern = re.compile(
        r'<(\d+),\s*"([^"]+)",\s*\[\s*([^]]+)\s*\],\s*(\d+)>,'
    )
    entries: list[dict[str, Any]] = []
    for match in pattern.finditer(path.read_text(encoding="utf-8")):
        entries.append(
            {
                "determinant": int(match.group(1)),
                "genus_symbol": match.group(2),
                "finite_face_counts_facets_to_vertices": [
                    int(value.strip()) for value in match.group(3).split(",")
                ],
                "cusp_count": int(match.group(4)),
            }
        )
    positions = [
        index for index, entry in enumerate(entries) if entry["genus_symbol"] == SK_GENUS
    ]
    if positions != [19]:
        raise ArithmeticError("Scharlau--Kirschmer List10 entry 20 moved")
    return positions[0] + 1, entries[positions[0]]


def _corpus_audit(
    corpus: Path, candidate_relations: Sequence[Sequence[str]]
) -> dict[str, Any]:
    ordinary_count = 0
    ordinary_h10: Counter[int] = Counter()
    ordinary_candidates: list[dict[str, Any]] = []
    ordinary_matches: list[dict[str, Any]] = []
    for path in sorted(corpus.rglob("*.coxiter")):
        facets, dimension, relations = _parse_coxiter(path)
        ordinary_count += 1
        if dimension == 10:
            ordinary_h10[facets] += 1
        if (facets, dimension) != (27, 10):
            continue
        relative = path.relative_to(corpus).as_posix()
        ordinary_candidates.append({"path": relative})
        if relations == candidate_relations:
            ordinary_matches.append(
                {"path": relative, "entrywise_relation_matrix_match": True}
            )

    doubling_path = corpus / DOUBLING_ARCHIVE_RELATIVE
    doubling_count = 0
    doubling_h10: Counter[int] = Counter()
    doubling_candidates: list[dict[str, Any]] = []
    doubling_matches: list[dict[str, Any]] = []
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"failed to read {member.name}")
            facets, dimension, relations = _parse_coxiter_text(
                handle.read().decode("utf-8")
            )
            doubling_count += 1
            if dimension == 10:
                doubling_h10[facets] += 1
            if (facets, dimension) != (27, 10):
                continue
            doubling_candidates.append({"path": member.name})
            if relations == candidate_relations:
                doubling_matches.append(
                    {"path": member.name, "entrywise_relation_matrix_match": True}
                )
    return {
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "ordinary_file_count": ordinary_count,
        "ordinary_H10_facet_histogram": {
            str(facets): count for facets, count in sorted(ordinary_h10.items())
        },
        "ordinary_H10_27_candidates": ordinary_candidates,
        "ordinary_exact_weighted_diagram_matches": ordinary_matches,
        "doubling_archive_relative_path": DOUBLING_ARCHIVE_RELATIVE.as_posix(),
        "doubling_archive_sha256": _sha256(doubling_path),
        "doubling_file_count": doubling_count,
        "doubling_H10_facet_histogram": {
            str(facets): count for facets, count in sorted(doubling_h10.items())
        },
        "doubling_H10_27_candidates": doubling_candidates,
        "doubling_exact_weighted_diagram_matches": doubling_matches,
    }


def _discriminant_overlattice_audit(
    gram: sp.Matrix, reflections: Sequence[sp.Matrix]
) -> dict[str, Any]:
    diagonal, left, _ = smith_normal_decomp(
        DomainMatrix.from_Matrix(gram).convert_to(ZZ)
    )
    diagonal = diagonal.to_Matrix()
    left = left.to_Matrix()
    moduli = [abs(int(diagonal[index, index])) for index in range(gram.rows)]
    elements = list(product(*(range(modulus) for modulus in moduli)))
    index_of = {element: index for index, element in enumerate(elements)}
    left_inverse = left.inv()
    gram_inverse = gram.inv()
    discriminant_gram = left_inverse.T * gram_inverse * left_inverse
    denominator = math.lcm(*(int(value.q) for value in discriminant_gram))
    scaled_discriminant_gram = [
        [
            int(denominator * discriminant_gram[row, column])
            for column in range(gram.cols)
        ]
        for row in range(gram.rows)
    ]
    linear_forms = [
        tuple(
            sum(
                element[row] * scaled_discriminant_gram[row][column]
                for row in range(gram.rows)
            )
            for column in range(gram.cols)
        )
        for element in elements
    ]

    def add(left_index: int, right_index: int) -> int:
        return index_of[
            tuple(
                (elements[left_index][index] + elements[right_index][index])
                % moduli[index]
                for index in range(gram.rows)
            )
        ]

    def bilinear_numerator(left_index: int, right_index: int) -> int:
        return sum(
            linear_forms[left_index][index] * elements[right_index][index]
            for index in range(gram.rows)
        )

    def integral_pair(left_index: int, right_index: int) -> bool:
        return bilinear_numerator(left_index, right_index) % denominator == 0

    actions: list[list[int]] = []
    for reflection in reflections:
        action = left * gram * reflection * gram_inverse * left_inverse
        if any(value.is_Integer is not True for value in action):
            raise ArithmeticError("reflection did not act integrally on A_L")
        action_integer = [
            [int(action[row, column]) for column in range(gram.cols)]
            for row in range(gram.rows)
        ]
        permutation: list[int] = []
        for element in elements:
            image = tuple(
                sum(
                    action_integer[row][column] * element[column]
                    for column in range(gram.cols)
                )
                % moduli[row]
                for row in range(gram.rows)
            )
            permutation.append(index_of[image])
        actions.append(permutation)

    def orbit(seed: int) -> frozenset[int]:
        result = {seed}
        queue = [seed]
        while queue:
            element = queue.pop()
            for action in actions:
                image = action[element]
                if image not in result:
                    result.add(image)
                    queue.append(image)
        return frozenset(result)

    def additive_closure(generators: Sequence[int]) -> frozenset[int]:
        subgroup = {0}
        for generator in generators:
            old = list(subgroup)
            cyclic: list[int] = []
            current = generator
            while current:
                cyclic.append(current)
                current = add(current, generator)
            subgroup = {
                add(element, multiple)
                for element in old
                for multiple in (0, *cyclic)
            }
        return frozenset(subgroup)

    self_integral = [
        index for index in range(len(elements)) if integral_pair(index, index)
    ]
    even_norm = [
        index
        for index in self_integral
        if (bilinear_numerator(index, index) // denominator) % 2 == 0
    ]
    distinct_orbits: dict[frozenset[int], list[int]] = {}
    invariant_integral_closures: dict[frozenset[int], list[int]] = {}
    for seed in self_integral:
        seed_orbit = orbit(seed)
        distinct_orbits.setdefault(seed_orbit, []).append(seed)
        if all(
            integral_pair(left_index, right_index)
            for left_index in seed_orbit
            for right_index in seed_orbit
        ):
            closure = additive_closure(tuple(seed_orbit))
            invariant_integral_closures.setdefault(closure, []).append(seed)

    nontrivial = [
        closure for closure in invariant_integral_closures if len(closure) > 1
    ]
    if len(nontrivial) != 1:
        raise ArithmeticError("invariant integral overlattice is not unique")
    maximal = nontrivial[0]

    basis_generator_cosets: list[int] = []
    for column in (5, 10):
        dual_coordinate = left * gram * OVERLATTICE_BASIS[:, column]
        if any(value.is_Integer is not True for value in dual_coordinate):
            raise ArithmeticError("overlattice basis column is not in L dual")
        element = tuple(
            int(dual_coordinate[index]) % moduli[index]
            for index in range(gram.rows)
        )
        basis_generator_cosets.append(index_of[element])
    basis_closure = additive_closure(basis_generator_cosets)
    if basis_closure != maximal:
        raise ArithmeticError("explicit basis does not generate the unique closure")

    maximal_even = all(
        (bilinear_numerator(index, index) // denominator) % 2 == 0
        for index in maximal
    )
    if not maximal_even:
        raise ArithmeticError("unique invariant integral overlattice is not even")

    return {
        "discriminant_group_smith_moduli": moduli,
        "discriminant_group_order": len(elements),
        "discriminant_bilinear_common_denominator": denominator,
        "self_integral_coset_count": len(self_integral),
        "even_norm_coset_count": len(even_norm),
        "self_integral_reflection_orbit_size_histogram": {
            str(size): count
            for size, count in sorted(
                Counter(len(seed_orbit) for seed_orbit in distinct_orbits).items()
            )
        },
        "invariant_integral_subgroup_count": len(
            invariant_integral_closures
        ),
        "invariant_integral_subgroup_orders": sorted(
            len(closure) for closure in invariant_integral_closures
        ),
        "unique_nontrivial_subgroup_cosets_in_smith_coordinates": [
            list(elements[index]) for index in sorted(maximal)
        ],
        "nonzero_seeds_whose_orbit_closure_is_integral": [
            list(elements[index])
            for index in sorted(invariant_integral_closures[maximal])
        ],
        "explicit_basis_generator_cosets_in_smith_coordinates": [
            list(elements[index]) for index in basis_generator_cosets
        ],
        "unique_nontrivial_subgroup_order": len(maximal),
        "unique_nontrivial_subgroup_even": maximal_even,
        "exhaustion_argument": (
            "For every self-integral discriminant coset x, close its orbit "
            "under all 27 reflections additively. An invariant integral "
            "overlattice containing x must contain this closure. Only x=0 "
            "or the three nonzero elements of the displayed order-four "
            "subgroup have integral closure; hence the trivial and displayed "
            "subgroups are the only invariant integral subgroups."
        ),
        "only_trivial_and_unique_order_four_subgroups_survive": (
            len(invariant_integral_closures) == 2
            and sorted(len(closure) for closure in invariant_integral_closures)
            == [1, 4]
        ),
    }


def build_certificate(
    coxiter_corpus: Path = COXITER_CORPUS,
    coxiter_archive: Path = COXITER_ARCHIVE,
    sk_directory: Path = SK_DIRECTORY,
    sk_archive: Path = SK_ARCHIVE,
) -> dict[str, Any]:
    case = json.loads(CASE.read_text(encoding="utf-8"))
    cone = json.loads(CONE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    weighted_inventory = json.loads(
        WEIGHTED_INVENTORY.read_text(encoding="utf-8")
    )
    weighted_chamber = next(
        record
        for record in weighted_inventory["all_chambers"]
        if record["chamber_index"] == 4
    )
    if not (
        case["id"] == "D0004_C007"
        and case["case"]["subgroup_order"] == 4
        and case["case"]["fixed_rank"] == 11
        and case["exact_reduction"]["hyperbolic_facet_count"] == 27
        and case["coxeter_gate"]["all_pairs_coxeter"]
        and cone["decision"]
        == {
            "classification": "finite_volume_noncompact",
            "compact": False,
            "criterion": (
                "complete pointed projective cone has no spacelike/"
                "hyperideal extreme ray"
            ),
            "finite_volume": True,
        }
        and cone["exact_double_description"]["finite_vertex_count"] == 895
        and cone["exact_double_description"]["ideal_vertex_count"] == 67
        and cone["exact_double_description"]["hyperideal_vertex_count"] == 0
        and vinberg["decision"]["finite_volume"]
        and not vinberg["decision"]["compact"]
        and vinberg["literal_endpoint_audit"]["edge_count"] == 6_221
        and vinberg["literal_endpoint_audit"]["bad_edge_count"] == 0
        and vinberg["literal_endpoint_audit"]
        ["bad_parabolic_neighborhood_count"]
        == 0
        and cone["source"]["case_sha256"] == _sha256(CASE)
        and _sha256(WEIGHTED_INVENTORY) == WEIGHTED_INVENTORY_SHA256
        and weighted_chamber
        == {
            "archived_stabilizer_order": 48,
            "chamber_index": 4,
            "groups_equal": True,
            "wall_count": 158,
            "weighted_wall_group_order": 48,
        }
    ):
        raise ArithmeticError("an upstream D0004_C007 gate failed")

    fixed_gram = _matrix(case["saturated_fixed_lattice"]["gram"])
    fixed_smith = _smith_diagonal(fixed_gram)
    facet_indices = [
        index - 1
        for index in case["exact_reduction"]["hyperbolic_facet_orbit_indices"]
    ]
    fixed_roots: list[sp.Matrix] = []
    fixed_reflections: list[sp.Matrix] = []
    fixed_types: Counter[tuple[int, int]] = Counter()
    mirror_records: list[dict[str, Any]] = []
    for facet_position, orbit_index in enumerate(facet_indices, start=1):
        source = case["folded_walls"]["orbit_records"][orbit_index]
        root = sp.Matrix(source["primitive_fixed_normal_coordinates"])
        reflection = _reflection(fixed_gram, root)
        root_type = _root_type(fixed_gram, root)
        if not _integral_reflection(fixed_gram, root):
            raise ArithmeticError("a fixed-lattice mirror is not reflective")
        fixed_roots.append(root)
        fixed_reflections.append(reflection)
        fixed_types[root_type] += 1
        mirror_records.append(
            {
                "facet_position": facet_position,
                "wall_orbit_index": orbit_index + 1,
                "source_wall_indices": source["source_wall_indices"],
                "primitive_fixed_coordinates": [int(value) for value in root],
                "fixed_lattice_square": root_type[0],
                "fixed_lattice_divisibility": root_type[1],
                "reflection_integral_on_fixed_lattice": True,
            }
        )
    fixed_root_matrix = sp.Matrix.hstack(*fixed_roots)
    root_span_smith = _smith_diagonal(fixed_root_matrix)

    discriminant_audit = _discriminant_overlattice_audit(
        fixed_gram, fixed_reflections
    )
    basis = OVERLATTICE_BASIS
    stable_gram = basis.T * fixed_gram * basis
    stable_smith = _smith_diagonal(stable_gram)
    standard_gram = _standard_gram()
    delta_witness = standard_gram.inv()[:, 8]
    delta_witness_square = (delta_witness.T * standard_gram * delta_witness)[0]
    if sp.Matrix(lattice_gram()) != standard_gram:
        raise ArithmeticError("maximal Vinberg module changed lattice Gram")
    isometry = OVERLATTICE_TO_STANDARD

    stable_types: Counter[tuple[int, int]] = Counter()
    stable_contents: Counter[int] = Counter()
    for record, root, reflection in zip(
        mirror_records, fixed_roots, fixed_reflections, strict=True
    ):
        stable_coordinates = basis.inv() * root
        if any(value.is_Integer is not True for value in stable_coordinates):
            raise ArithmeticError("a fixed root left the stable lattice")
        content = math.gcd(*(abs(int(value)) for value in stable_coordinates))
        primitive_stable = stable_coordinates / content
        stable_type = _root_type(stable_gram, primitive_stable)
        stable_reflection = basis.inv() * reflection * basis
        target_coordinates = isometry.inv() * primitive_stable
        if not (
            all(value.is_Integer is True for value in stable_reflection)
            and _integral_reflection(stable_gram, primitive_stable)
            and all(value.is_Integer is True for value in target_coordinates)
        ):
            raise ArithmeticError("a mirror did not preserve the stable lattice")
        stable_types[stable_type] += 1
        stable_contents[content] += 1
        record.update(
            {
                "stable_coordinates_before_primitive_normalization": [
                    int(value) for value in stable_coordinates
                ],
                "stable_content": content,
                "primitive_stable_coordinates": [
                    int(value) for value in primitive_stable
                ],
                "primitive_standard_coordinates": [
                    int(value) for value in target_coordinates
                ],
                "stable_lattice_square": stable_type[0],
                "stable_lattice_divisibility": stable_type[1],
                "reflection_integral_on_stable_lattice": True,
            }
        )

    sk_list_path = sk_directory / "List10"
    sk_gram_path = sk_directory / "Gram10"
    sk_position, sk_entry = _parse_sk_list10(sk_list_path)
    sk_grams = ast.literal_eval(sk_gram_path.read_text(encoding="utf-8"))
    sk_gram = sp.Matrix(11, 11, sk_grams[sk_position - 1])
    sk_isometry = STANDARD_TO_SK_SIGNED

    maximal_search = run_search(4, attempt_stopping=True)
    maximal_roots = [sp.Matrix(root) for root in maximal_search.roots]
    maximal_root_types = Counter(
        (pair(root, root), divisibility(root)) for root in maximal_search.roots
    )
    maximal_reflections_integral = all(
        _integral_reflection(standard_gram, root) for root in maximal_roots
    )
    maximal_positive_gram = [
        [str(-pair(left, right)) for right in maximal_search.roots]
        for left in maximal_search.roots
    ]
    maximal_relations = _relation_matrix(maximal_positive_gram)
    maximal_relation_histogram = Counter(
        maximal_relations[left][right]
        for left in range(len(maximal_relations))
        for right in range(left + 1, len(maximal_relations))
    )

    candidate_relations = _relation_matrix(case["coxeter_gate"]["positive_wall_gram"])
    candidate_connected, candidate_component_search = _connected_elliptic_masks(
        case["coxeter_gate"]["positive_wall_gram"]
    )
    candidate_euler = _euler_from_connected_components(
        candidate_relations, candidate_connected
    )
    maximal_connected = _connected_masks_direct(maximal_relations)
    maximal_euler = _euler_from_connected_components(
        maximal_relations, maximal_connected
    )
    candidate_chi = Fraction(candidate_euler["rational_euler_characteristic"])
    maximal_chi = Fraction(maximal_euler["rational_euler_characteristic"])
    subgroup_index = candidate_chi / maximal_chi
    volume_coefficient = -candidate_chi * Fraction(32, 945)
    candidate_full_face_vector = [
        cone["exact_double_description"]["finite_vertex_count"]
        + cone["exact_double_description"]["ideal_vertex_count"],
        *reversed(candidate_euler["finite_special_subset_count_by_rank"][1:10]),
        1,
    ]
    maximal_full_face_vector = [
        maximal_search.stopping.finite_vertex_count
        + maximal_search.stopping.cusp_count,
        *reversed(maximal_euler["finite_special_subset_count_by_rank"][1:10]),
        1,
    ]

    corpus_audit = _corpus_audit(coxiter_corpus, candidate_relations)
    expected_ordinary_h10 = {"11": 14, "12": 8, "13": 4, "14": 1}
    expected_doubling_h10 = {
        "12": 3,
        "13": 4,
        "14": 3,
        "15": 1,
        "16": 3,
        "19": 1,
        "20": 3,
        "28": 2,
        "44": 2,
    }
    expected_candidate_special = [
        1,
        27,
        264,
        1_435,
        5_034,
        12_061,
        19_934,
        22_234,
        15_784,
        6_221,
        895,
    ]
    expected_maximal_special = [
        1,
        14,
        89,
        340,
        868,
        1_554,
        1_986,
        1_800,
        1_114,
        424,
        70,
    ]
    expected_candidate_full_faces = [
        962,
        6_221,
        15_784,
        22_234,
        19_934,
        12_061,
        5_034,
        1_435,
        264,
        27,
        1,
    ]
    expected_maximal_full_faces = [
        75,
        424,
        1_114,
        1_800,
        1_986,
        1_554,
        868,
        340,
        89,
        14,
        1,
    ]
    checks = {
        "case_and_both_volume_certificates_replay": True,
        "archived_stabilizer_is_the_full_weighted_wall_group": (
            weighted_chamber["groups_equal"]
            and weighted_chamber["archived_stabilizer_order"] == 48
            and weighted_chamber["weighted_wall_group_order"] == 48
        ),
        "fixed_lattice_invariants_are_expected": (
            fixed_gram.det() == 512
            and exact_inertia(fixed_gram) == (1, 10, 0)
            and fixed_smith == [1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4]
        ),
        "all_27_fixed_mirrors_are_integral_reflections": all(
            record["reflection_integral_on_fixed_lattice"]
            for record in mirror_records
        ),
        "twenty_seven_roots_generate_the_fixed_lattice": (
            fixed_root_matrix.rank() == 11 and root_span_smith == [1] * 11
        ),
        "fixed_root_types_are_expected": fixed_types
        == Counter({(-2, 2): 10, (-2, 1): 6, (-4, 2): 8, (-8, 4): 3}),
        "discriminant_exhaustion_is_complete": (
            discriminant_audit["discriminant_group_order"] == 512
            and discriminant_audit["self_integral_coset_count"] == 256
            and discriminant_audit["only_trivial_and_unique_order_four_subgroups_survive"]
        ),
        "explicit_overlattice_basis_has_index_four": basis.det()
        == sp.Rational(1, 4),
        "stable_lattice_invariants_are_expected": (
            stable_gram.det() == 32
            and exact_inertia(stable_gram) == (1, 10, 0)
            and all(int(stable_gram[index, index]) % 2 == 0 for index in range(11))
            and stable_smith == [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
            and delta_witness_square == sp.Rational(-1, 2)
        ),
        "all_27_reflections_preserve_the_stable_lattice": all(
            record["reflection_integral_on_stable_lattice"]
            for record in mirror_records
        ),
        "three_fixed_roots_become_twice_primitive": stable_contents
        == Counter({1: 24, 2: 3}),
        "stable_root_types_are_expected": stable_types
        == Counter({(-2, 2): 10, (-2, 1): 9, (-4, 2): 8}),
        "explicit_standard_isometry_replays": (
            isometry.det() == 1
            and isometry.T * stable_gram * isometry == standard_gram
        ),
        "scharlau_kirschmer_sources_are_pinned": (
            _sha256(sk_archive) == SK_ARCHIVE_SHA256
            and _sha256(sk_list_path) == SK_LIST10_SHA256
            and _sha256(sk_gram_path) == SK_GRAM10_SHA256
        ),
        "scharlau_kirschmer_entry20_is_explicitly_isometric": (
            sk_position == 20
            and sk_isometry.det() == -1
            and sk_isometry.T * sk_gram * sk_isometry == -standard_gram
        ),
        "maximal_vinberg_height_zero_system_has_type_B6_B3_A1": (
            len(height_zero_simple_roots()) == 10
            and all(
                pair(left, right) >= 0
                for position, left in enumerate(height_zero_simple_roots())
                for right in height_zero_simple_roots()[position + 1 :]
            )
        ),
        "maximal_vinberg_search_has_14_integral_roots": (
            len(maximal_search.roots) == 14 and maximal_reflections_integral
        ),
        "maximal_root_types_are_expected": maximal_root_types
        == Counter({(-2, 1): 8, (-4, 2): 4, (-2, 2): 2}),
        "maximal_relation_histogram_is_expected": maximal_relation_histogram
        == Counter({"2": 76, "3": 9, "4": 4, "0": 2}),
        "maximal_cone_is_finite_volume_noncompact": (
            maximal_search.stopping.fired
            and maximal_search.stopping.finite_vertex_count == 70
            and maximal_search.stopping.cusp_count == 5
            and maximal_search.stopping.hyperideal_ray_count == 0
            and maximal_search.stopping.opposite_sheet_ray_count == 0
        ),
        "scharlau_kirschmer_closes_maximal_root_completeness": sk_entry
        == {
            "determinant": 32,
            "genus_symbol": SK_GENUS,
            "finite_face_counts_facets_to_vertices": expected_maximal_special[1:],
            "cusp_count": 5,
        },
        "candidate_special_counts_match_literal_vinberg": (
            candidate_euler["finite_special_subset_count_by_rank"]
            == expected_candidate_special
            and expected_candidate_special[9]
            == vinberg["literal_endpoint_audit"]["edge_count"]
            and expected_candidate_special[10]
            == vinberg["literal_endpoint_audit"]["finite_vertex_count"]
        ),
        "maximal_special_counts_match_scharlau_kirschmer": (
            maximal_euler["finite_special_subset_count_by_rank"]
            == expected_maximal_special
            and expected_maximal_special[1:]
            == sk_entry["finite_face_counts_facets_to_vertices"]
        ),
        "candidate_euler_is_expected": candidate_chi
        == Fraction(-527, 1_658_880),
        "maximal_euler_is_expected": maximal_chi
        == Fraction(-527, 3_344_302_080),
        "index_is_expected": subgroup_index == 2_016,
        "volume_is_expected": volume_coefficient == Fraction(527, 48_988_800),
        "candidate_full_face_vector_reconstructs": candidate_full_face_vector
        == expected_candidate_full_faces,
        "maximal_full_face_vector_reconstructs": maximal_full_face_vector
        == expected_maximal_full_faces,
        "coxiter_sources_are_pinned": (
            _sha256(coxiter_archive) == COXITER_ARCHIVE_SHA256
            and corpus_audit["doubling_archive_sha256"]
            == DOUBLING_ARCHIVE_SHA256
        ),
        "ordinary_corpus_scope_is_expected": (
            corpus_audit["ordinary_file_count"] == 1_214
            and corpus_audit["ordinary_H10_facet_histogram"]
            == expected_ordinary_h10
        ),
        "doubling_corpus_scope_is_expected": (
            corpus_audit["doubling_file_count"] == 2_959
            and corpus_audit["doubling_H10_facet_histogram"]
            == expected_doubling_h10
        ),
        "no_bounded_corpus_H10_27_record": (
            corpus_audit["ordinary_H10_27_candidates"] == []
            and corpus_audit["ordinary_exact_weighted_diagram_matches"] == []
            and corpus_audit["doubling_H10_27_candidates"] == []
            and corpus_audit["doubling_exact_weighted_diagram_matches"] == []
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"D0004_C007 lineage checks failed: {checks}")

    return {
        "schema": "shimada_disc16_d4_c7_h10_27_lineage_v1",
        "candidate": {
            "id": "D16-D4-C7-H10-27",
            "source_chamber_index": 4,
            "archived_subgroup_class_index": 7,
            "hyperbolic_dimension": 10,
            "facet_count": 27,
            "finite_vertex_count": 895,
            "ideal_vertex_count": 67,
            "compact": False,
            "finite_volume": True,
            "source_chamber_archived_group_equals_full_weighted_wall_group": True,
        },
        "sources": {
            "archived_fixed_section_case": str(CASE.relative_to(ROOT)),
            "archived_fixed_section_case_sha256": _sha256(CASE),
            "exact_cone_certificate": str(CONE.relative_to(ROOT)),
            "exact_cone_certificate_sha256": _sha256(CONE),
            "literal_vinberg_certificate": str(VINBERG.relative_to(ROOT)),
            "literal_vinberg_certificate_sha256": _sha256(VINBERG),
            "full_weighted_wall_inventory": str(
                WEIGHTED_INVENTORY.relative_to(ROOT)
            ),
            "full_weighted_wall_inventory_sha256": _sha256(
                WEIGHTED_INVENTORY
            ),
            "vinberg_primary_source": (
                "https://doi.org/10.1070/SM1972v016n01ABEH001346"
            ),
            "scharlau_kirschmer_source": (
                "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html"
            ),
            "scharlau_kirschmer_archive_sha256": _sha256(sk_archive),
            "coxiter_repository_commit": COXITER_COMMIT,
            "coxiter_archive_sha256": _sha256(coxiter_archive),
        },
        "saturated_divisor_fixed_lattice": {
            "name": "L",
            "rank": 11,
            "inertia": list(exact_inertia(fixed_gram)),
            "gram_in_archived_basis": _matrix_payload(fixed_gram),
            "determinant": int(fixed_gram.det()),
            "smith_diagonal": fixed_smith,
            "discriminant_group": "(Z/2)^5 + (Z/4)^2",
            "even": True,
            "not_the_same_integral_lattice_as_the_stable_overlattice": True,
        },
        "primitive_candidate_mirrors": {
            "facet_count": 27,
            "fixed_lattice_square_divisibility_histogram": {
                f"{square},{root_divisibility}": count
                for (square, root_divisibility), count in sorted(fixed_types.items())
            },
            "fixed_root_matrix_rank": fixed_root_matrix.rank(),
            "fixed_root_span_smith_diagonal": root_span_smith,
            "roots_generate_entire_fixed_lattice": True,
            "records": mirror_records,
        },
        "unique_reflection_stable_integral_overlattice": {
            "name": "M",
            "scope_of_uniqueness": (
                "integral overlattices of L invariant under all 27 candidate "
                "reflections"
            ),
            "discriminant_exhaustion": discriminant_audit,
            "index_over_L": 4,
            "basis_columns_in_L_tensor_Q": _matrix_payload(basis),
            "basis_description_one_based": (
                "replace e6 by (e5+e6)/2 and e11 by (e9+e11)/2"
            ),
            "gram": _matrix_payload(stable_gram),
            "determinant": int(stable_gram.det()),
            "smith_diagonal": stable_smith,
            "even": True,
            "two_elementary_invariants": {"r": 11, "a": 5, "delta": 1},
            "delta_one_witness_in_standard_dual": [
                _rational(value) for value in delta_witness
            ],
            "delta_one_witness_square": str(delta_witness_square),
            "candidate_root_content_histogram": {
                str(content): count for content, count in sorted(stable_contents.items())
            },
            "primitive_candidate_root_square_divisibility_histogram": {
                f"{square},{root_divisibility}": count
                for (square, root_divisibility), count in sorted(stable_types.items())
            },
            "all_27_reflections_integral": True,
            "standard_lattice": "U + D6(-1) + A1(-1)^3",
            "standard_gram": _matrix_payload(standard_gram),
            "explicit_isometry_convention": "Q^T * G_M * Q = G_standard",
            "explicit_isometry_matrix": _matrix_payload(isometry),
            "explicit_isometry_determinant": int(isometry.det()),
        },
        "known_full_reflection_chamber": {
            "lattice": "U + D6(-1) + A1(-1)^3",
            "scharlau_kirschmer_list10_entry": sk_position,
            "scharlau_kirschmer_record": sk_entry,
            "scharlau_kirschmer_gram": _matrix_payload(sk_gram),
            "signed_isometry_convention": (
                "S^T * Gram10[20] * S = -G_standard"
            ),
            "signed_isometry_matrix": _matrix_payload(sk_isometry),
            "signed_isometry_determinant": int(sk_isometry.det()),
            "all_root_vinberg_control_vector": [1, 1, *([0] * 9)],
            "height_zero_type": "B6 + B3 + A1",
            "height_zero_roots": [
                list(root) for root in height_zero_simple_roots()
            ],
            "searched_height_levels": [
                {
                    "control": level.control,
                    "norm_abs": level.norm_abs,
                    "height": str(Fraction(level.control**2, level.norm_abs)),
                    "candidate_count": level.candidate_count,
                    "accepted_roots": [list(root) for root in level.accepted_roots],
                }
                for level in maximal_search.levels
            ],
            "explicit_roots_in_standard_coordinates": [
                list(root) for root in maximal_search.roots
            ],
            "facet_count": len(maximal_search.roots),
            "primitive_square_divisibility_histogram": {
                f"{square},{root_divisibility}": count
                for (square, root_divisibility), count in sorted(
                    maximal_root_types.items()
                )
            },
            "all_reflections_integral": maximal_reflections_integral,
            "primitive_root_norm_exhaustion": (
                "M is even with discriminant exponent two. For a primitive "
                "reflection root r, div(r) divides two and -r^2 divides "
                "2*div(r); hence the only possible negative squares are "
                "-2 and -4, both enumerated."
            ),
            "normalized_relation_matrix": maximal_relations,
            "relation_histogram": {
                relation: count
                for relation, count in sorted(maximal_relation_histogram.items())
            },
            "finite_vertex_count": maximal_search.stopping.finite_vertex_count,
            "ideal_vertex_count": maximal_search.stopping.cusp_count,
            "hyperideal_ray_count": maximal_search.stopping.hyperideal_ray_count,
            "finite_volume": True,
            "compact": False,
            "completeness_argument": (
                "The standard Vinberg enumeration produces 14 fundamental "
                "roots. The pinned Scharlau--Kirschmer record for the "
                "explicitly isometric class-number-one lattice states that "
                "the fundamental domain of W(M) has exactly 14 facets; its "
                "entire finite-face vector and five cusps also agree. Hence "
                "no fundamental root remains unenumerated."
            ),
            "full_integral_root_reflection_group": True,
        },
        "solomon_index": {
            "method": "complete Solomon sum over every finite special subgroup",
            "candidate_connected_elliptic_components": [
                [index + 1 for index in range(27) if mask & (1 << index)]
                for mask in sorted(candidate_connected)
            ],
            "candidate_component_search_by_order": candidate_component_search,
            "candidate_27_wall_group": candidate_euler,
            "full_14_wall_group": maximal_euler,
            "inclusion": (
                "W_27 < W(M), because every candidate mirror is an integral "
                "root reflection of M"
            ),
            "index_formula": "chi(W_27) / chi(W(M))",
            "index": subgroup_index.numerator,
            "index_prime_factorization": _factorization(subgroup_index.numerator),
            "hyperbolic_covolume_formula": (
                "vol = -(32*pi^5/945)*chi in H10"
            ),
            "candidate_hyperbolic_covolume": "pi^5 * 527/48988800",
            "full_group_hyperbolic_covolume": (
                "pi^5 * 527/98761420800"
            ),
            "candidate_full_face_vector": candidate_full_face_vector,
            "maximal_full_face_vector": maximal_full_face_vector,
        },
        "bounded_weighted_diagram_audit": {
            "scope": (
                "pinned CoxIter ordinary corpus and its separate reflection-"
                "doubling archive only"
            ),
            "candidate_relation_matrix": candidate_relations,
            "coxiter_snapshot": corpus_audit,
            "literal_27_facet_record_present": False,
            "does_not_prove_worldwide_absence": True,
        },
        "lineage_decision": {
            "new_positive_algebraic_extraction": True,
            "new_fixed_lattice_claimed": False,
            "new_reflection_stable_lattice": False,
            "new_maximal_reflection_group": False,
            "new_arithmetic_commensurability_class": False,
            "literal_27_wall_subgroup_chamber_priority": (
                "open after bounded corpus exclusion"
            ),
            "divisor_fixed_section_realization_priority": "open",
            "worldwide_priority_claimed": False,
            "novelty_boundary": (
                "the positive divisor-lattice extraction and possibly the "
                "literal index-2016 27-wall subgroup chamber"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coxiter-corpus", type=Path, default=COXITER_CORPUS)
    parser.add_argument("--coxiter-archive", type=Path, default=COXITER_ARCHIVE)
    parser.add_argument("--sk-directory", type=Path, default=SK_DIRECTORY)
    parser.add_argument("--sk-archive", type=Path, default=SK_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(
        arguments.coxiter_corpus,
        arguments.coxiter_archive,
        arguments.sk_directory,
        arguments.sk_archive,
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
