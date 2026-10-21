#!/usr/bin/env python3
"""Certify the literature lineage of the Borcherds H13/50 and H14/98 sections.

The companion construction and Vinberg certificates establish two literal
finite-volume Coxeter polytopes.  This supplement answers the separate
lineage question.

* ``SS2-2A-H13-50`` lives on the known even 2-elementary lattice
  ``(14,10,0)``.  Alexeev's 16-wall full reflection chamber is recovered
  directly as a facet of his published ``(15,11,1)`` matrix.  Thus the
  50-wall group is finite index in a known full reflection group; the index
  and worldwide priority of the literal 50-wall chamber remain open.
* ``SS2-2C-H14-98`` is exactly Alexeev's published ``(15,9,1)`` 98-wall
  chamber after doubled-dual short/long-root exchange.  An explicit
  98-letter permutation checks every Gram entry.  The fixed divisor lattice
  has a unique reflection-stable integral index-two overlattice, the odd
  doubled-dual lattice in Scharlau--Kirschmer List14 entry 7.

The bounded CoxIter scan is deliberately not promoted to a worldwide
novelty statement.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, deque
import hashlib
import json
from pathlib import Path
import re
import tarfile
from typing import Any, Iterable
import zipfile

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import DiagramEnumerator
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_construction.json"
)
VINBERG = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_vinberg.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_lineage.json"
)
DEFAULT_ALEXEEV_DATA = Path(
    "/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt"
)
DEFAULT_ALEXEEV_TEX = Path("/private/tmp/alexeev_2209_source/g0.tex")
DEFAULT_SK_ARCHIVE = Path("/private/tmp/Scharlau-Kirschmer-Lorentz.zip")
DEFAULT_COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")

SCHEMA = "coxeter.borcherds_char2_h13_h14_lineage.v1"
CONSTRUCTION_SHA256 = (
    "a9087200b360c4bed75cd93c987bd7e7d493ff1ce852bb9c5ea4dc3490f1fe11"
)
VINBERG_SHA256 = (
    "032c8343a6dd7f4904b9daecb074b8b0bbbf1a94f9431252884484c98d2e8456"
)
ALEXEEV_DATA_SHA256 = (
    "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
)
ALEXEEV_TEX_SHA256 = (
    "006e2a1eac2aa3c4916ef66aceb27e8e80a35145162e6785a4874843f249f1df"
)
SK_ARCHIVE_SHA256 = (
    "74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5"
)
LIST13_SHA256 = "4a57d3082cb50ab3f9ee2e294358c84e32c3e729c275fb7c846c2d8183dc40d8"
GRAM13_SHA256 = "b0ab3e916036ab463e5f1acbc44ff6eeb140e4871ca799cd43e8d248312f77f8"
LIST14_SHA256 = "be1bb032fe4ecc66045df3251ab10a61ecd5c654f0727f61ae983f8fb5eb9231"
GRAM14_SHA256 = "3ea06c7b4d8628f6a3816cc788af7e4fd8875c5b0c073ae8b9ca06f4062ac736"
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
COXITER_DOUBLING_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)

H13_FULL_FACE_VECTOR = [
    16, 120, 560, 1820, 4364, 7965, 11232,
    12272, 10300, 6488, 2920, 842, 111,
]
H14_FULL_FACE_VECTOR = [
    98, 1666, 14616, 83559, 342888, 1056573, 2490674,
    4501714, 6177500, 6299076, 4598006, 2249968, 646786, 75898,
]
SK_LIST13_ENTRY10 = (
    16,
    "II_{13,1}(2^{-4})",
    H13_FULL_FACE_VECTOR,
    6,
)
SK_LIST14_ENTRY7 = (
    64,
    "I_{14,1}(2^{6})",
    H14_FULL_FACE_VECTOR,
    2114,
)

# Candidate-wall index -> Alexeev-wall index.  The corresponding root has
# the opposite length: candidate -4 maps to Alexeev -2 and conversely.
H14_ALEXEEV_PERMUTATION = [
    4, 0, 14, 29, 11, 10, 18, 8, 13, 2, 1, 24, 31, 6, 5, 77,
    26, 19, 28, 16, 23, 27, 22, 32, 36, 21, 56, 59, 12, 15, 67, 3,
    17, 54, 68, 78, 39, 87, 88, 7, 81, 71, 52, 38, 85, 74, 35, 64,
    90, 20, 70, 96, 57, 66, 79, 58, 69, 33, 76, 51, 47, 62, 86, 94,
    55, 89, 42, 41, 73, 92, 80, 75, 25, 34, 37, 97, 9, 30, 63, 91,
    72, 84, 65, 50, 53, 49, 45, 44, 83, 46, 48, 61, 60, 93, 40, 43,
    82, 95,
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _matrix_digest(matrix: sp.Matrix) -> str:
    raw = json.dumps(
        _matrix_payload(matrix), separators=(",", ":"), sort_keys=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def _smith_absolute(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(smith[index, index])) for index in range(matrix.rows)]


def _mod_one(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - sp.floor(value)


def _mod_two(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - 2 * sp.floor(value / 2)


def _even_discriminant_q_histogram(matrix: sp.Matrix) -> dict[str, int]:
    """Return q(x)=x^2 mod 2Z for an even lattice's discriminant group."""

    inverse = matrix.inv()
    generators = [
        tuple(_mod_one(inverse[row, column]) for row in range(matrix.rows))
        for column in range(matrix.cols)
    ]
    zero = tuple(sp.Rational(0) for _ in range(matrix.rows))
    classes = {zero}
    queue = deque([zero])
    while queue:
        current = queue.popleft()
        for generator in generators:
            candidate = tuple(
                _mod_one(left + right)
                for left, right in zip(current, generator, strict=True)
            )
            if candidate not in classes:
                classes.add(candidate)
                queue.append(candidate)
    histogram = Counter(
        str(_mod_two((sp.Matrix(vector).T * matrix * sp.Matrix(vector))[0]))
        for vector in classes
    )
    return {key: histogram[key] for key in sorted(histogram, key=sp.Rational)}


def _edge_histogram(negative_gram: sp.Matrix) -> dict[str, int]:
    labels: Counter[str] = Counter()
    for left in range(negative_gram.rows):
        for right in range(left + 1, negative_gram.cols):
            product = sp.Rational(negative_gram[left, right])
            ratio = sp.factor(
                product**2
                / (negative_gram[left, left] * negative_gram[right, right])
            )
            label = (
                "2" if product == 0
                else "3" if ratio == sp.Rational(1, 4)
                else "4" if ratio == sp.Rational(1, 2)
                else "6" if ratio == sp.Rational(3, 4)
                else "infinity" if ratio == 1
                else "ultraparallel" if ratio > 1
                else "forbidden"
            )
            labels[label] += 1
    return {key: labels[key] for key in sorted(labels)}


def _h13_full_chamber(gram_15_11_1: sp.Matrix) -> dict[str, Any]:
    section_root = 21  # zero based; member of Alexeev's size-36 root orbit
    projected_squares = {
        index: sp.factor(
            gram_15_11_1[index, index]
            + gram_15_11_1[section_root, index] ** 2 / 2
        )
        for index in range(gram_15_11_1.rows)
        if index != section_root
    }
    negative_projected = [
        index for index, square in projected_squares.items() if square < 0
    ]
    orthogonal = [
        index
        for index in range(gram_15_11_1.rows)
        if index != section_root and gram_15_11_1[section_root, index] == 0
    ]
    full_gram = gram_15_11_1.extract(orthogonal, orthogonal)
    enumerator = DiagramEnumerator(-full_gram, dimension=13)
    face_vector = [
        len(enumerator.finite_subdiagrams(rank)) for rank in range(1, 14)
    ]
    affine = enumerator.connected_affine_subdiagrams()
    neighborhoods: set[int] = set()
    for component_mask in affine:
        component = [
            index for index in range(16) if component_mask & (1 << index)
        ]
        zero_neighbors = [
            index
            for index in range(16)
            if not component_mask & (1 << index)
            and all(full_gram[index, node] == 0 for node in component)
        ]
        neighborhood = component_mask | sum(1 << index for index in zero_neighbors)
        nodes = [index for index in range(16) if neighborhood & (1 << index)]
        if full_gram.extract(nodes, nodes).rank() == 12:
            neighborhoods.add(neighborhood)
    checks = {
        "section_root_is_minus_two": gram_15_11_1[section_root, section_root] == -2,
        "exactly_16_negative_projected_walls": (
            negative_projected == orthogonal and len(orthogonal) == 16
        ),
        "root_norm_histogram_12_long_4_short": Counter(
            int(full_gram[index, index]) for index in range(16)
        ) == Counter({-4: 12, -2: 4}),
        "edge_histogram_expected": _edge_histogram(full_gram)
        == {"2": 104, "3": 14, "4": 2},
        "finite_face_vector_matches_SK": face_vector == H13_FULL_FACE_VECTOR,
        "six_maximal_parabolic_neighborhoods": len(neighborhoods) == 6,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Alexeev H13 facet replay changed: {checks}")
    return {
        "parent_case": "(15,11,1)",
        "parent_wall_count": 66,
        "section_root_zero_based": section_root,
        "section_root_orbit_size": 36,
        "orthogonal_parent_wall_indices_zero_based": orthogonal,
        "facet_gram": _matrix_payload(full_gram),
        "facet_gram_sha256": _matrix_digest(full_gram),
        "facet_count": 16,
        "root_norm_histogram": {"-4": 12, "-2": 4},
        "coxeter_edge_histogram": _edge_histogram(full_gram),
        "finite_special_subset_counts_by_rank": face_vector,
        "connected_affine_component_count": len(affine),
        "maximal_parabolic_neighborhood_count": len(neighborhoods),
        "checks": checks,
    }


def _reflection_matrix_mod_two(gram: sp.Matrix, root: sp.Matrix) -> list[list[int]]:
    square = int((root.T * gram * root)[0])
    action = sp.eye(gram.rows) - sp.Rational(2, square) * root * (root.T * gram)
    if any(not value.is_Integer for value in action):
        raise ArithmeticError("a production wall ceased to be an integral reflection")
    return [
        [int(action[row, column]) & 1 for column in range(gram.cols)]
        for row in range(gram.rows)
    ]


def _matvec_mod_two(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(value * vector[column] for column, value in enumerate(row)) & 1
        for row in matrix
    )


def _h14_stable_overlattice(case: dict[str, Any]) -> dict[str, Any]:
    gram = sp.Matrix(case["fixed_lattice"]["gram"])
    roots = [
        sp.Matrix(record["primitive_fixed_coordinates"])
        for record in case["restricted_chamber"]["primitive_facet_lattice_records"]
    ]
    actions = [_reflection_matrix_mod_two(gram, root) for root in roots]
    integral_classes: list[tuple[int, ...]] = []
    stable_classes: list[tuple[int, ...]] = []
    for bits in range(1, 1 << gram.rows):
        vector = tuple((bits >> index) & 1 for index in range(gram.rows))
        column = sp.Matrix(vector)
        if any(int(value) & 1 for value in gram * column):
            continue
        if int((column.T * gram * column)[0]) % 4:
            continue
        integral_classes.append(vector)
        if all(_matvec_mod_two(action, vector) == vector for action in actions):
            stable_classes.append(vector)

    expected = (0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if stable_classes != [expected]:
        raise ArithmeticError(f"stable index-two closure changed: {stable_classes}")
    vector = sp.Matrix(expected)
    pivot = expected.index(1)
    basis = sp.Matrix.hstack(
        vector / 2,
        *(sp.eye(gram.rows)[:, index] for index in range(gram.rows) if index != pivot),
    )
    closure_gram = basis.T * gram * basis
    checks = {
        "127_integral_nonzero_order_two_classes": len(integral_classes) == 127,
        "unique_W98_stable_class": len(stable_classes) == 1,
        "generator_square_minus_one": (vector.T * gram * vector)[0] / 4 == -1,
        "basis_has_index_two": abs(basis.det()) == sp.Rational(1, 2),
        "closure_integral_odd": (
            all(value.is_Integer for value in closure_gram)
            and any(int(closure_gram[index, index]) & 1 for index in range(15))
        ),
        "closure_signature_1_14": exact_inertia(closure_gram) == (1, 14, 0),
        "closure_determinant_64": int(closure_gram.det()) == 64,
        "closure_smith_1_power_9_2_power_6": (
            _smith_absolute(closure_gram) == [1] * 9 + [2] * 6
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"H14 stable-overlattice audit failed: {checks}")
    return {
        "integral_nonzero_order_two_overlattice_class_count": len(integral_classes),
        "W98_stable_class_count": len(stable_classes),
        "unique_stable_class_vector_mod_2": list(expected),
        "adjoined_vector": "(e_2+e_3)/2 in the archived fixed basis",
        "adjoined_vector_square": -1,
        "overlattice_index": 2,
        "basis_columns_in_fixed_Q_space": _matrix_payload(basis),
        "gram": _matrix_payload(closure_gram),
        "determinant": 64,
        "smith_diagonal_absolute": [1] * 9 + [2] * 6,
        "signature_positive_negative": [1, 14],
        "parity": "odd",
        "checks": checks,
    }


def _verify_h14_doubled_dual(
    candidate_gram: sp.Matrix, alexeev_gram: sp.Matrix
) -> dict[str, Any]:
    permutation = H14_ALEXEEV_PERMUTATION
    if sorted(permutation) != list(range(98)):
        raise ArithmeticError("the frozen H14 permutation is not bijective")
    transformed = sp.zeros(98)
    for left in range(98):
        alex_left = permutation[left]
        left_scale = sp.Rational(1, 1) if alexeev_gram[alex_left, alex_left] == -2 else sp.Rational(1, 2)
        for right in range(98):
            alex_right = permutation[right]
            right_scale = sp.Rational(1, 1) if alexeev_gram[alex_right, alex_right] == -2 else sp.Rational(1, 2)
            # S^dagger=S^*(2): retain a -2 root, halve a div-2 -4 root,
            # and multiply the bilinear form by two.
            transformed[left, right] = (
                2
                * left_scale
                * right_scale
                * alexeev_gram[alex_left, alex_right]
            )
    direct_counts = []
    enumerator = DiagramEnumerator(-candidate_gram, dimension=14)
    for rank in (1, 2, 3):
        direct_counts.append(len(enumerator.finite_subdiagrams(rank)))
    checks = {
        "candidate_diagonal_84_long_14_short": Counter(candidate_gram.diagonal())
        == Counter({-4: 84, -2: 14}),
        "alexeev_diagonal_84_short_14_long": Counter(alexeev_gram.diagonal())
        == Counter({-2: 84, -4: 14}),
        "every_root_length_is_exchanged": all(
            candidate_gram[index, index]
            * alexeev_gram[permutation[index], permutation[index]]
            == 8
            for index in range(98)
        ),
        "exact_doubled_dual_weighted_gram_equality": transformed == candidate_gram,
        "direct_finite_ranks_1_2_3": direct_counts == H14_FULL_FACE_VECTOR[:3],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"H14 doubled-dual match failed: {checks}")
    return {
        "source_case": "(15,9,1)",
        "source_wall_count": 98,
        "candidate_to_alexeev_permutation_zero_based": permutation,
        "root_exchange": {
            "candidate_minus4_to_alexeev_minus2": 84,
            "candidate_minus2_to_alexeev_minus4": 14,
        },
        "formula": (
            "S^dagger=S^*(2): keep each Alexeev -2 root, replace each "
            "divisibility-2 -4 root by half of itself, and double the form"
        ),
        "source_facet_gram_sha256": _matrix_digest(alexeev_gram),
        "transformed_facet_gram_sha256": _matrix_digest(transformed),
        "candidate_facet_gram_sha256": _matrix_digest(candidate_gram),
        "direct_finite_special_counts_ranks_1_2_3": direct_counts,
        "checks": checks,
    }


def _header(lines: Iterable[str]) -> tuple[int, int]:
    content = next(
        line.split("#", 1)[0].strip()
        for line in lines
        if line.split("#", 1)[0].strip()
    )
    facets, dimension = map(int, content.split()[:2])
    return dimension, facets


def _bounded_coxiter_scan(root: Path) -> dict[str, Any]:
    targets = {(13, 50), (14, 98)}
    ordinary: Counter[tuple[int, int]] = Counter()
    ordinary_matches = {target: [] for target in targets}
    for path in sorted(root.rglob("*.coxiter")):
        key = _header(path.read_text(encoding="utf-8", errors="replace").splitlines())
        ordinary[key] += 1
        if key in targets:
            ordinary_matches[key].append(path.relative_to(root).as_posix())
    doubling_path = root / "graphs/testing_index2.tar.gz"
    doubling: Counter[tuple[int, int]] = Counter()
    doubling_matches = {target: [] for target in targets}
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not (member.isfile() and member.name.endswith(".coxiter")):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"could not read {member.name}")
            key = _header(raw.decode("utf-8", errors="replace") for raw in handle)
            doubling[key] += 1
            if key in targets:
                doubling_matches[key].append(member.name)
    rows = {}
    for dimension, facets in sorted(targets):
        key = (dimension, facets)
        rows[f"H{dimension}_facets_{facets}"] = {
            "ordinary_header_match_count": ordinary[key],
            "ordinary_header_matches": ordinary_matches[key],
            "doubling_header_match_count": doubling[key],
            "doubling_header_matches": doubling_matches[key],
            "absent_from_both_by_dimension_and_facet_count": not (
                ordinary[key] or doubling[key]
            ),
        }
    archive_path = root.parent.with_suffix(".zip")
    checks = {
        "ordinary_snapshot_1214_files": sum(ordinary.values()) == 1214,
        "doubling_snapshot_2959_files": sum(doubling.values()) == 2959,
        "doubling_archive_hash_pinned": _sha256(doubling_path) == COXITER_DOUBLING_SHA256,
        "ordinary_archive_hash_pinned": (
            archive_path.exists() and _sha256(archive_path) == COXITER_ARCHIVE_SHA256
        ),
        "both_literal_sizes_absent": all(
            record["absent_from_both_by_dimension_and_facet_count"]
            for record in rows.values()
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"bounded CoxIter scan changed: {checks}")
    return {
        "scope": (
            "pinned 1,214-file ordinary CoxIter snapshot and separate "
            "2,959-file index-two reflection-doubling snapshot"
        ),
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "ordinary_archive_path": str(archive_path),
        "ordinary_archive_sha256": COXITER_ARCHIVE_SHA256,
        "doubling_archive_path": str(doubling_path),
        "doubling_archive_sha256": COXITER_DOUBLING_SHA256,
        "targets": rows,
        "warning": (
            "facet-count absence is only a bounded negative check; H14 is "
            "nevertheless an exact match to Alexeev outside this corpus"
        ),
        "checks": checks,
    }


def _source_materials(
    alexeev_data_path: Path,
    alexeev_tex_path: Path,
    sk_archive_path: Path,
) -> tuple[dict[str, Any], sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    if _sha256(alexeev_data_path) != ALEXEEV_DATA_SHA256:
        raise ArithmeticError("Alexeev all_cases source hash changed")
    if _sha256(alexeev_tex_path) != ALEXEEV_TEX_SHA256:
        raise ArithmeticError("Alexeev TeX source hash changed")
    source = alexeev_data_path.read_text(encoding="utf-8")
    tex = alexeev_tex_path.read_text(encoding="utf-8")
    required_tex_fragments = [
        "There are $66$ roots",
        "$(14,10,0)$ and $(14,10,1)$",
        "$\\Gamma_r$ has $98$ roots",
        "\\PGL(2,7)",
        "the doubled dual $S^\\dag:=S^*(2)$",
    ]
    if not all(fragment in tex for fragment in required_tex_fragments):
        raise ArithmeticError("Alexeev's pinned lineage statements changed")
    gram_15_11_1 = _extract_matrix(source, (15, 11, 1))
    gram_15_9_1 = _extract_matrix(source, (15, 9, 1))
    alexeev_matrix_inventory = []
    assignment_pattern = re.compile(
        r"^coxeter_matrix\[\((\d+),(\d+),(\d+)\)\] = matrix\((\[.*\])\)$"
    )
    for line in source.splitlines():
        match = assignment_pattern.fullmatch(line)
        if match is None:
            continue
        key = tuple(int(match.group(index)) for index in range(1, 4))
        rows = ast.literal_eval(match.group(4))
        alexeev_matrix_inventory.append({
            "lattice_invariants_r_a_delta": list(key),
            "hyperbolic_dimension": key[0] - 1,
            "facet_count": len(rows),
        })

    archive_raw = sk_archive_path.read_bytes()
    if _sha256_bytes(archive_raw) != SK_ARCHIVE_SHA256:
        raise ArithmeticError("Scharlau--Kirschmer archive hash changed")
    with zipfile.ZipFile(sk_archive_path) as archive:
        list13_raw = archive.read("List13")
        gram13_raw = archive.read("Gram13")
        list14_raw = archive.read("List14")
        gram14_raw = archive.read("Gram14")
    hashes = {
        "List13": _sha256_bytes(list13_raw),
        "Gram13": _sha256_bytes(gram13_raw),
        "List14": _sha256_bytes(list14_raw),
        "Gram14": _sha256_bytes(gram14_raw),
    }
    if hashes != {
        "List13": LIST13_SHA256,
        "Gram13": GRAM13_SHA256,
        "List14": LIST14_SHA256,
        "Gram14": GRAM14_SHA256,
    }:
        raise ArithmeticError(f"Scharlau--Kirschmer member hashes changed: {hashes}")
    list13 = ast.literal_eval(
        list13_raw.decode("ascii").replace("<", "(").replace(">", ")")
    )
    list14 = ast.literal_eval(
        list14_raw.decode("ascii").replace("<", "(").replace(">", ")")
    )
    if tuple(list13[9]) != SK_LIST13_ENTRY10:
        raise ArithmeticError("Scharlau--Kirschmer List13 entry 10 changed")
    if tuple(list14[6]) != SK_LIST14_ENTRY7:
        raise ArithmeticError("Scharlau--Kirschmer List14 entry 7 changed")
    sk13_inventory = [
        {
            "entry_one_based": index,
            "genus": entry[1],
            "facet_count": entry[2][0],
        }
        for index, entry in enumerate(list13, 1)
    ]
    source_inventory_checks = {
        "alexeev_18_explicit_matrix_assignments": len(alexeev_matrix_inventory) == 18,
        "alexeev_no_H13_50_assignment": not any(
            row["hyperbolic_dimension"] == 13 and row["facet_count"] == 50
            for row in alexeev_matrix_inventory
        ),
        "alexeev_only_50_wall_assignment_is_H19": [
            row
            for row in alexeev_matrix_inventory
            if row["facet_count"] == 50
        ] == [{
            "lattice_invariants_r_a_delta": [20, 2, 1],
            "hyperbolic_dimension": 19,
            "facet_count": 50,
        }],
        "SK_List13_has_17_entries": len(sk13_inventory) == 17,
        "SK_List13_has_no_50_wall_full_chamber": not any(
            row["facet_count"] == 50 for row in sk13_inventory
        ),
    }
    if not all(source_inventory_checks.values()):
        raise ArithmeticError(
            f"bounded H13 primary-source inventory changed: {source_inventory_checks}"
        )
    grams13 = ast.literal_eval(gram13_raw.decode("ascii"))
    grams14 = ast.literal_eval(gram14_raw.decode("ascii"))
    sk13 = sp.Matrix(14, 14, grams13[9])
    sk14 = sp.Matrix(15, 15, grams14[6])
    record = {
        "alexeev": {
            "paper": "Valery Alexeev, Reflective hyperbolic 2-elementary lattices, arXiv:2209.09110v4",
            "arxiv_url": "https://arxiv.org/abs/2209.09110",
            "all_cases_path": str(alexeev_data_path),
            "all_cases_sha256": ALEXEEV_DATA_SHA256,
            "tex_path": str(alexeev_tex_path),
            "tex_sha256": ALEXEEV_TEX_SHA256,
            "matrix_cases_reparsed": ["(15,11,1)", "(15,9,1)"],
        },
        "scharlau_kirschmer": {
            "source_page": "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html",
            "archive_path": str(sk_archive_path),
            "archive_sha256": SK_ARCHIVE_SHA256,
            "member_sha256": hashes,
            "list13_entry_one_based": 10,
            "list13_entry": list(SK_LIST13_ENTRY10),
            "list14_entry_one_based": 7,
            "list14_entry": list(SK_LIST14_ENTRY7),
        },
        "bounded_H13_literal_primary_source_comparison": {
            "scope": (
                "all 18 explicit Coxeter-matrix assignments in Alexeev's "
                "pinned ancillary file and all 17 full H13 chambers in "
                "Scharlau--Kirschmer List13"
            ),
            "alexeev_explicit_matrix_inventory": alexeev_matrix_inventory,
            "scharlau_kirschmer_List13_inventory": sk13_inventory,
            "literal_H13_50_found": False,
            "conclusion": (
                "the known (14,10,0) 16-wall full chamber is present, but "
                "the literal H13/50 subgroup chamber is absent from these "
                "two pinned primary-source tables"
            ),
            "warning": (
                "this is not an audit of all papers or unpublished subgroup chambers"
            ),
            "checks": source_inventory_checks,
        },
    }
    return record, gram_15_11_1, gram_15_9_1, sk13, sk14


def build_certificate(
    construction_path: Path = CONSTRUCTION,
    vinberg_path: Path = VINBERG,
    alexeev_data_path: Path = DEFAULT_ALEXEEV_DATA,
    alexeev_tex_path: Path = DEFAULT_ALEXEEV_TEX,
    sk_archive_path: Path = DEFAULT_SK_ARCHIVE,
    coxiter_root: Path = DEFAULT_COXITER_ROOT,
) -> dict[str, Any]:
    if _sha256(construction_path) != CONSTRUCTION_SHA256:
        raise ArithmeticError("pinned H13/H14 construction certificate changed")
    if _sha256(vinberg_path) != VINBERG_SHA256:
        raise ArithmeticError("pinned H13/H14 Vinberg certificate changed")
    construction = json.loads(construction_path.read_text(encoding="utf-8"))
    vinberg = json.loads(vinberg_path.read_text(encoding="utf-8"))
    construction_cases = {case["id"]: case for case in construction["cases"]}
    vinberg_cases = {case["id"]: case for case in vinberg["cases"]}
    if not construction["all_checks_pass"] or not vinberg["all_checks_pass"]:
        raise ArithmeticError("a production certificate reports a failed check")
    sources, parent_h13, alexeev_h14, sk13, sk14 = _source_materials(
        alexeev_data_path, alexeev_tex_path, sk_archive_path
    )

    h13_case = construction_cases["SS2-2A-H13-50"]
    h13_vinberg = vinberg_cases["SS2-2A-H13-50"]
    h13_gram = sp.Matrix(h13_case["fixed_lattice"]["gram"])
    h13_dual = 2 * h13_gram.inv()
    h13_full = _h13_full_chamber(parent_h13)
    h13_q = _even_discriminant_q_histogram(h13_dual)
    sk13_q = _even_discriminant_q_histogram(-sk13)
    h13_checks = {
        "production_and_vinberg_pass": (
            h13_case["all_checks_pass"] and h13_vinberg["all_checks_pass"]
        ),
        "finite_volume_noncompact": (
            h13_vinberg["finite_volume"] and not h13_vinberg["compact"]
        ),
        "fixed_lattice_is_even_14_10_0": (
            exact_inertia(h13_gram) == (1, 13, 0)
            and int(h13_gram.det()) == -1024
            and _smith_absolute(h13_gram) == [1] * 4 + [2] * 10
            and all(not int(h13_gram[index, index]) & 1 for index in range(14))
            and all((2 * value).is_Integer for value in h13_gram.inv())
            and all(value.is_Integer for value in h13_gram.inv().diagonal())
        ),
        "doubled_dual_integral_even_det16": (
            all(value.is_Integer for value in h13_dual)
            and all(not int(h13_dual[index, index]) & 1 for index in range(14))
            and int(h13_dual.det()) == -16
            and _smith_absolute(h13_dual) == [1] * 10 + [2] * 4
        ),
        "doubled_dual_discriminant_form_matches_SK_entry10": (
            h13_q == sk13_q == {"0": 6, "1": 10}
        ),
        "full_16_wall_chamber_replayed": all(h13_full["checks"].values()),
        "full_face_vector_and_cusps_match_SK_entry10": (
            h13_full["finite_special_subset_counts_by_rank"]
            == list(SK_LIST13_ENTRY10[2])
            and h13_full["maximal_parabolic_neighborhood_count"]
            == SK_LIST13_ENTRY10[3]
        ),
    }
    if not all(h13_checks.values()):
        raise ArithmeticError(f"H13 lineage checks failed: {h13_checks}")

    h14_case = construction_cases["SS2-2C-H14-98"]
    h14_vinberg = vinberg_cases["SS2-2C-H14-98"]
    h14_candidate_gram = sp.Matrix(h14_case["restricted_chamber"]["facet_gram"])
    h14_closure = _h14_stable_overlattice(h14_case)
    h14_match = _verify_h14_doubled_dual(h14_candidate_gram, alexeev_h14)
    closure_gram = sp.Matrix(h14_closure["gram"])
    h14_checks = {
        "production_and_vinberg_pass": (
            h14_case["all_checks_pass"] and h14_vinberg["all_checks_pass"]
        ),
        "finite_volume_noncompact": (
            h14_vinberg["finite_volume"] and not h14_vinberg["compact"]
        ),
        "unique_reflection_stable_index_two_closure": all(
            h14_closure["checks"].values()
        ),
        "closure_matches_SK_entry7_basic_genus_data": (
            exact_inertia(closure_gram) == (1, 14, 0)
            and int(closure_gram.det()) == 64
            and _smith_absolute(closure_gram) == [1] * 9 + [2] * 6
            and any(int(closure_gram[index, index]) & 1 for index in range(15))
            and int(sk14.det()) == -64
            and _smith_absolute(sk14) == [1] * 9 + [2] * 6
        ),
        "exact_Alexeev_doubled_dual_98_wall_match": all(h14_match["checks"].values()),
        "full_finite_vector_equals_SK_entry7": (
            h14_vinberg["finite_special_subset_counts_by_rank"]
            == list(SK_LIST14_ENTRY7[2])
        ),
        "full_cusp_count_equals_SK_entry7": (
            h14_vinberg["ideal_vertex_count"] == SK_LIST14_ENTRY7[3]
        ),
        "rank_two_count_is_1666_not_1656": (
            h14_vinberg["finite_special_subset_counts_by_rank"][1] == 1666
            and sum(
                h14_case["restricted_chamber"]["coxeter_edge_counts"][label]
                for label in ("2", "3", "4", "6")
            ) == 1666
        ),
    }
    if not all(h14_checks.values()):
        raise ArithmeticError(f"H14 lineage checks failed: {h14_checks}")

    coxiter = _bounded_coxiter_scan(coxiter_root)
    checks = {
        "production_certificate_pinned": True,
        "vinberg_certificate_pinned": True,
        "primary_sources_pinned_and_reparsed": True,
        "H13_known_lattice_and_full_chamber": all(h13_checks.values()),
        "H14_exact_Alexeev_rediscovery": all(h14_checks.values()),
        "bounded_CoxIter_scan_replayed": all(coxiter["checks"].values()),
    }
    return {
        "schema": SCHEMA,
        "id": "SS2-2A-H13-50+SS2-2C-H14-98-LINEAGE",
        "production_certificates": {
            "construction": {
                "path": str(construction_path.relative_to(ROOT)),
                "sha256": CONSTRUCTION_SHA256,
            },
            "vinberg": {
                "path": str(vinberg_path.relative_to(ROOT)),
                "sha256": VINBERG_SHA256,
                "method": vinberg["method"],
            },
        },
        "primary_sources": sources,
        "H13_50": {
            "id": "SS2-2A-H13-50",
            "literal_polytope": {
                "dimension": 13,
                "facet_count": 50,
                "finite_volume_noncompact": True,
                "finite_special_subset_counts_by_rank": h13_vinberg[
                    "finite_special_subset_counts_by_rank"
                ],
                "ideal_vertex_count": h13_vinberg["ideal_vertex_count"],
            },
            "fixed_lattice": {
                "nikulin_invariants": [14, 10, 0],
                "classification_method": (
                    "Nikulin uniqueness for indefinite even 2-elementary lattices"
                ),
                "determinant": -1024,
                "smith_diagonal_absolute": [1] * 4 + [2] * 10,
                "doubled_dual_gram": _matrix_payload(h13_dual),
                "doubled_dual_determinant": -16,
                "doubled_dual_smith_diagonal_absolute": [1] * 10 + [2] * 4,
                "doubled_dual_discriminant_q_histogram": h13_q,
                "SK_doubled_dual_entry": "List13/Gram13 entry 10",
            },
            "known_full_reflection_chamber": h13_full,
            "group_relation": {
                "statement": (
                    "the 50-wall reflection group is finite index in the full "
                    "reflection group of (14,10,0); its chamber is tiled by "
                    "copies of the known 16-wall chamber"
                ),
                "finite_index_proof": (
                    "both chambers have finite hyperbolic covolume on the same lattice"
                ),
                "index_value": None,
                "why_Solomon_does_not_supply_it": (
                    "both Coxeter Euler characteristics vanish in odd hyperbolic dimension"
                ),
            },
            "priority": {
                "dimension_previously_populated": True,
                "lattice_and_full_chamber_known": True,
                "literal_50_wall_diagram_absent_from_pinned_CoxIter_snapshots": True,
                "worldwide_literal_chamber_priority": "open",
                "claim_permitted": (
                    "new exact algebro-geometric construction candidate for a "
                    "50-wall finite-index subgroup chamber, not a new dimension or lattice"
                ),
            },
            "checks": h13_checks,
        },
        "H14_98": {
            "id": "SS2-2C-H14-98",
            "literal_polytope": {
                "dimension": 14,
                "facet_count": 98,
                "finite_volume_noncompact": True,
                "finite_special_subset_counts_by_rank": h14_vinberg[
                    "finite_special_subset_counts_by_rank"
                ],
                "ideal_vertex_count": h14_vinberg["ideal_vertex_count"],
            },
            "reflection_stable_overlattice": h14_closure,
            "exact_Alexeev_match": h14_match,
            "priority": {
                "dimension_previously_populated": True,
                "abstract_polytope": (
                    "known: exact doubled-dual realization of Alexeev (15,9,1)"
                ),
                "abstract_chamber_reflection_index": 1,
                "new_hyperbolic_polytope_claim": False,
                "possible_remaining_novelty": (
                    "the particular fixed-sublattice/Frobenius derivation only"
                ),
                "geometric_warning": (
                    "the 2c Frobenius wall symmetry is not in the geometric K3 subgroup"
                ),
            },
            "checks": h14_checks,
        },
        "bounded_CoxIter_comparison": coxiter,
        "decision": {
            "H13": (
                "positive finite-volume H13/50 example; known (14,10,0) "
                "arithmetic stock and known 16-wall full chamber; literal "
                "50-wall subgroup chamber priority and index remain open"
            ),
            "H14": (
                "positive but not new: exact Alexeev (15,9,1) 98-wall chamber "
                "under coodd doubled-dual root exchange"
            ),
            "new_abstract_hyperbolic_examples_certified_by_this_lineage_audit": 0,
            "literal_priority_candidates_still_alive": ["SS2-2A-H13-50"],
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--construction", type=Path, default=CONSTRUCTION)
    parser.add_argument("--vinberg", type=Path, default=VINBERG)
    parser.add_argument("--alexeev-data", type=Path, default=DEFAULT_ALEXEEV_DATA)
    parser.add_argument("--alexeev-tex", type=Path, default=DEFAULT_ALEXEEV_TEX)
    parser.add_argument("--sk-archive", type=Path, default=DEFAULT_SK_ARCHIVE)
    parser.add_argument("--coxiter-root", type=Path, default=DEFAULT_COXITER_ROOT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(
        args.construction,
        args.vinberg,
        args.alexeev_data,
        args.alexeev_tex,
        args.sk_archive,
        args.coxiter_root,
    )
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("rebuilt H13/H14 lineage differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "id": payload["id"],
        "H13": payload["decision"]["H13"],
        "H14": payload["decision"]["H14"],
        "all_checks_pass": payload["all_checks_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
