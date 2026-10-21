#!/usr/bin/env python3
"""Exact certificate for the compact complete K3 root return theorem.

The eight cases are the compact rank-three/rank-four 2-reflective lattices
classified by Nikulin and Vinberg.  The lattice and complete (-2)-curve
Gram matrices are reproduced in Roulleau, arXiv:1909.01909v5, Sections 3--4
and the appendix.

This script checks only the exact finite matrix consequences used by the
project report.  The exhaustion of all lattices is a cited classification
theorem, not a computer-search claim made by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import tarfile
from math import gcd
from typing import Any

from sympy import Matrix, Poly, ilcm, symbols


SOURCE_TAR_SHA256 = (
    "ec64645999be6cc0e0ad3ef03fbf046a5e641f2c180e370ef2c91d052ec805a7"
)
SOURCE_TEX_SHA256 = (
    "effc17df3a96e4ead5fa8cee247e0fdffe178b970e6c76586592577ce9b62589"
)
SOURCE_TEX_MEMBER = "GeomericModelsK3FiniteAuteGroup.tex"


CASES: dict[str, dict[str, Any]] = {
    "S1": {
        "ns_gram": [[6, 0, 0], [0, -2, 0], [0, 0, -2]],
        "wall_gram": [
            [-2, 0, 4, 6, 4, 0],
            [0, -2, 0, 4, 6, 4],
            [4, 0, -2, 0, 4, 6],
            [6, 4, 0, -2, 0, 4],
            [4, 6, 4, 0, -2, 0],
            [0, 4, 6, 4, 0, -2],
        ],
        "wall_basis": [0, 1, 2],
        "expected_vertices": 6,
        "expected_weyl_vector": True,
        "geometric_model": "double plane with three 6-tangent conics",
    },
    "S2": {
        "ns_gram": [[36, 0, 0], [0, -2, 1], [0, 1, -2]],
        "wall_gram": [
            [-2, 1, 7, 10, 7, 1],
            [1, -2, 1, 7, 10, 7],
            [7, 1, -2, 1, 7, 10],
            [10, 7, 1, -2, 1, 7],
            [7, 10, 7, 1, -2, 1],
            [1, 7, 10, 7, 1, -2],
        ],
        "wall_basis": [0, 1, 2],
        "expected_vertices": 6,
        "expected_weyl_vector": True,
        "geometric_model": "quartic with three split quadric sections",
    },
    "S3": {
        "ns_gram": [[12, 0, 0], [0, -2, 1], [0, 1, -2]],
        "wall_gram": [
            [-2, 1, 4, 1],
            [1, -2, 1, 4],
            [4, 1, -2, 1],
            [1, 4, 1, -2],
        ],
        "wall_basis": [0, 1, 2],
        "expected_vertices": 4,
        "expected_weyl_vector": True,
        "geometric_model": "quartic with two split hyperplane sections",
    },
    "S4": {
        "ns_gram": [[2, 1, 2], [1, -2, 1], [2, 1, -2]],
        "wall_gram": [
            [-2, 1, 3, 1],
            [1, -2, 1, 6],
            [3, 1, -2, 1],
            [1, 6, 1, -2],
        ],
        "wall_basis": [0, 1, 2],
        "expected_vertices": 4,
        "expected_weyl_vector": False,
        "geometric_model": (
            "double plane with one tritangent line and one 6-tangent conic"
        ),
    },
    "S5": {
        "ns_gram": [[4, 0, 0], [0, -2, 1], [0, 1, -2]],
        "wall_gram": [
            [-2, 1, 3, 0],
            [1, -2, 0, 3],
            [3, 0, -2, 1],
            [0, 3, 1, -2],
        ],
        "wall_basis": [0, 1, 2],
        "expected_vertices": 4,
        "expected_weyl_vector": True,
        "geometric_model": "double plane with two tritangent lines",
    },
    "S6": {
        "ns_gram": [[2, 2, 3], [2, -2, 1], [3, 1, -2]],
        "wall_gram": [
            [-2, 6, 1, 5, 5, 1],
            [6, -2, 5, 1, 1, 5],
            [1, 5, -2, 11, 0, 9],
            [5, 1, 11, -2, 9, 0],
            [5, 1, 0, 9, -2, 11],
            [1, 5, 9, 0, 11, -2],
        ],
        "wall_basis": [0, 2, 4],
        "expected_vertices": 6,
        "expected_weyl_vector": False,
        "geometric_model": (
            "double plane with a 6-tangent conic and two cuspidal cubics"
        ),
    },
    "L24": {
        "ns_gram": [
            [2, 1, 1, 1],
            [1, -2, 0, 0],
            [1, 0, -2, 0],
            [1, 0, 0, -2],
        ],
        "wall_gram": [
            [-2, 0, 0, 3, 1, 1],
            [0, -2, 0, 1, 3, 1],
            [0, 0, -2, 1, 1, 3],
            [3, 1, 1, -2, 0, 0],
            [1, 3, 1, 0, -2, 0],
            [1, 1, 3, 0, 0, -2],
        ],
        "wall_basis": [0, 1, 2, 3],
        "expected_vertices": 8,
        "expected_weyl_vector": True,
        "geometric_model": "double plane with three tritangent lines",
    },
    "L27": {
        "ns_gram": [
            [12, 2, 0, 0],
            [2, -2, 1, 0],
            [0, 1, -2, 1],
            [0, 0, 1, -2],
        ],
        "wall_gram": [
            [-2, 3, 1, 1, 1, 1, 1, 1],
            [3, -2, 1, 1, 1, 1, 1, 1],
            [1, 1, -2, 6, 4, 0, 4, 0],
            [1, 1, 6, -2, 0, 4, 0, 4],
            [1, 1, 4, 0, -2, 6, 4, 0],
            [1, 1, 0, 4, 6, -2, 0, 4],
            [1, 1, 4, 0, 4, 0, -2, 6],
            [1, 1, 0, 4, 0, 4, 6, -2],
        ],
        "wall_basis": [1, 3, 4, 6],
        "expected_vertices": 12,
        "expected_weyl_vector": False,
        "geometric_model": (
            "double plane with one tritangent line and three 6-tangent conics"
        ),
    },
}


# For each row, the columns of this unimodular matrix are the selected
# ``wall_basis`` roots written in the displayed ``ns_gram`` basis.  The
# rank-three and L24 matrices come directly from the curve formulae in
# Roulleau, Sections 3--4.  For L27, Roulleau explicitly says that
# A2,A4,A5,A7 generate NS; the matrix below is an exact integral congruence
# between his displayed Vinberg Gram and the corresponding curve-basis Gram.
# Freezing these changes of basis is important: equality of determinants alone
# would not prove that the listed walls span the displayed NS lattice.
NS_TO_WALL_BASIS: dict[str, list[list[int]]] = {
    "S1": [[0, 0, 1], [1, 0, -2], [0, 1, 0]],
    "S2": [[0, 0, 1], [1, 0, -5], [0, 1, -3]],
    "S3": [[0, 0, 1], [1, 0, -3], [0, 1, -2]],
    "S4": [[0, 0, 1], [1, 0, -1], [0, 1, 0]],
    "S5": [[0, 0, 1], [1, 0, -2], [0, 1, -1]],
    "S6": [[0, 0, 1], [1, 0, -1], [0, 1, 1]],
    "L24": [[0, 0, 0, 1], [1, 0, 0, -1], [0, 1, 0, 0], [0, 0, 1, 0]],
    "L27": [[-3, -1, -2, -1], [3, 1, 3, 0], [-3, 0, 0, -2], [2, 2, 2, 1]],
}


def _signature(matrix: Matrix) -> tuple[int, int, int]:
    variable = symbols("lambda")
    zero = matrix.cols - matrix.rank()
    polynomial = matrix.charpoly(variable).as_expr()
    if zero:
        polynomial = polynomial / variable**zero
    nonzero_polynomial = Poly(polynomial, variable)
    positive = negative = 0
    for (left, right), multiplicity in nonzero_polynomial.intervals():
        if right <= 0:
            negative += int(multiplicity)
        elif left >= 0:
            positive += int(multiplicity)
        else:  # pragma: no cover - zero factors were removed above
            raise AssertionError(f"root-isolation interval meets zero: {(left, right)}")
    assert positive + negative + zero == matrix.cols
    return positive, negative, zero


def _primitive_tuple(vector: Matrix) -> tuple[int, ...]:
    denominator = 1
    for entry in vector:
        denominator = int(ilcm(denominator, int(entry.q)))
    values = [int(entry * denominator) for entry in vector]
    divisor = 0
    for value in values:
        divisor = gcd(divisor, abs(value))
    values = [value // divisor for value in values]
    first = next(value for value in values if value)
    if first < 0:
        values = [-value for value in values]
    return tuple(values)


def _canonical_gram(matrix: Matrix) -> tuple[int, ...]:
    """Canonical flattened Gram under simultaneous row/column permutation."""

    rows = matrix.tolist()
    size = matrix.rows
    return min(
        tuple(
            rows[permutation[i]][permutation[j]]
            for i in range(size)
            for j in range(size)
        )
        for permutation in itertools.permutations(range(size))
    )


def _extreme_rays(wall_gram: Matrix, basis: list[int]) -> list[dict[str, Any]]:
    rank = len(basis)
    dimension = rank - 1
    ambient_gram = wall_gram.extract(basis, basis)
    coordinates = Matrix.hstack(
        *[
            ambient_gram.inv() * wall_gram.extract(basis, [index])
            for index in range(wall_gram.rows)
        ]
    )
    assert coordinates.T * ambient_gram * coordinates == wall_gram
    inequalities = coordinates.T * ambient_gram
    interior = coordinates * Matrix.ones(wall_gram.rows, 1)
    assert all(value > 0 for value in inequalities * interior)
    assert (interior.T * ambient_gram * interior)[0] > 0

    rays: dict[tuple[int, ...], dict[str, Any]] = {}
    for active in itertools.combinations(range(wall_gram.rows), dimension):
        nullspace = inequalities.extract(active, range(rank)).nullspace()
        if len(nullspace) != 1:
            continue
        vector = nullspace[0]
        values = inequalities * vector
        if all(value >= 0 for value in values):
            pass
        elif all(value <= 0 for value in values):
            vector = -vector
            values = -values
        else:
            continue
        key = _primitive_tuple(vector)
        zero_walls = [index for index, value in enumerate(values) if value == 0]
        norm = (vector.T * ambient_gram * vector)[0]
        assert norm > 0
        time_pairing = (vector.T * ambient_gram * interior)[0]
        assert time_pairing > 0
        rays[key] = {
            "active_walls": zero_walls,
            "positive_norm": str(norm),
            "positive_time_pairing": str(time_pairing),
        }
    return [rays[key] for key in sorted(rays)]


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check_source_archive(tar_path: Path, tex_path: Path) -> dict[str, Any]:
    """Verify the pinned primary-source archive and its classification text."""

    tar_bytes = tar_path.read_bytes()
    tex_bytes = tex_path.read_bytes()
    assert _sha256(tar_bytes) == SOURCE_TAR_SHA256
    assert _sha256(tex_bytes) == SOURCE_TEX_SHA256

    with tarfile.open(tar_path, mode="r:*") as archive:
        assert SOURCE_TEX_MEMBER in archive.getnames()
        source = archive.extractfile(SOURCE_TEX_MEMBER)
        assert source is not None
        archived_tex = source.read()

    assert archived_tex == tex_bytes
    assert b"there are $6$ types of lattices $S_{1},\\dots,S_{6}$" in tex_bytes
    assert b"there are only two lattices of rank $4$ with compact fundamental" in tex_bytes

    return {
        "archive_path": str(tar_path),
        "archive_sha256": _sha256(tar_bytes),
        "main_tex_path": str(tex_path),
        "main_tex_sha256": _sha256(tex_bytes),
        "archive_member_matches_main_tex": True,
        "rank_three_six_type_statement_present": True,
        "rank_four_two_type_statement_present": True,
        "all_checks_pass": True,
    }


def certify() -> dict[str, Any]:
    reports: list[dict[str, Any]] = []
    canonical_wall_grams: set[tuple[int, ...]] = set()
    for name, case in CASES.items():
        ns_gram = Matrix(case["ns_gram"])
        wall_gram = Matrix(case["wall_gram"])
        rank = ns_gram.rows
        dimension = rank - 1
        facets = wall_gram.rows

        assert ns_gram == ns_gram.T
        assert wall_gram == wall_gram.T
        assert _signature(ns_gram) == (1, rank - 1, 0)
        assert _signature(wall_gram) == (1, rank - 1, facets - rank)
        assert wall_gram.rank() == rank
        canonical_wall_gram = _canonical_gram(wall_gram)
        assert canonical_wall_gram not in canonical_wall_grams
        canonical_wall_grams.add(canonical_wall_gram)
        assert all(wall_gram[i, i] == -2 for i in range(facets))
        assert all(
            wall_gram[i, j] >= 0
            for i in range(facets)
            for j in range(facets)
            if i != j
        )

        basis = case["wall_basis"]
        root_basis_gram = wall_gram.extract(basis, basis)
        ns_to_wall_basis = Matrix(NS_TO_WALL_BASIS[name])
        assert ns_to_wall_basis.shape == (rank, rank)
        assert abs(ns_to_wall_basis.det()) == 1
        assert ns_to_wall_basis.T * ns_gram * ns_to_wall_basis == root_basis_gram
        assert _signature(root_basis_gram) == (1, rank - 1, 0)

        wall_coordinates_in_root_basis = Matrix.hstack(
            *[
                root_basis_gram.inv() * wall_gram.extract(basis, [index])
                for index in range(facets)
            ]
        )
        assert all(entry.q == 1 for entry in wall_coordinates_in_root_basis)
        assert wall_coordinates_in_root_basis.extract(
            range(rank), basis
        ) == Matrix.eye(rank)
        wall_coordinates_in_ns_basis = (
            ns_to_wall_basis * wall_coordinates_in_root_basis
        )
        assert (
            wall_coordinates_in_ns_basis.T
            * ns_gram
            * wall_coordinates_in_ns_basis
            == wall_gram
        )
        assert all(
            gcd(*[abs(int(value)) for value in wall_coordinates_in_ns_basis[:, j]])
            == 1
            for j in range(facets)
        )

        rays = _extreme_rays(wall_gram, basis)
        assert len(rays) == case["expected_vertices"]
        assert all(int(len(ray["active_walls"])) == dimension for ray in rays)
        assert all(not ray["positive_norm"].startswith("-") for ray in rays)
        active_facets = {
            wall for ray in rays for wall in ray["active_walls"]
        }
        assert active_facets == set(range(facets))

        finite_pair_counts = {"m2": 0, "m3": 0, "dotted": 0}
        finite_pairs: set[tuple[int, int]] = set()
        for i, j in itertools.combinations(range(facets), 2):
            value = int(wall_gram[i, j])
            if value == 0:
                finite_pair_counts["m2"] += 1
                finite_pairs.add((i, j))
            elif value == 1:
                finite_pair_counts["m3"] += 1
                finite_pairs.add((i, j))
            else:
                assert value >= 3  # no parallel pair; compactness is literal
                finite_pair_counts["dotted"] += 1

        pair_vertex_incidence: dict[tuple[int, int], int] = {}
        for ray in rays:
            for pair in itertools.combinations(ray["active_walls"], 2):
                normalized_pair = tuple(sorted(pair))
                pair_vertex_incidence[normalized_pair] = (
                    pair_vertex_incidence.get(normalized_pair, 0) + 1
                )
        assert set(pair_vertex_incidence) == finite_pairs
        expected_pair_vertex_incidence = 1 if dimension == 2 else 2
        assert set(pair_vertex_incidence.values()) == {
            expected_pair_vertex_incidence
        }

        augmented = wall_gram.row_join(Matrix.ones(facets, 1))
        has_weyl_vector = augmented.rank() == wall_gram.rank()
        assert has_weyl_vector is case["expected_weyl_vector"]
        weyl_vector_square: str | None = None
        if has_weyl_vector:
            wall_pairing_map = (
                wall_coordinates_in_root_basis.T * root_basis_gram
            )
            weyl_vector = wall_pairing_map.gauss_jordan_solve(
                Matrix.ones(facets, 1)
            )[0]
            assert not any(entry.free_symbols for entry in weyl_vector)
            assert wall_pairing_map * weyl_vector == Matrix.ones(facets, 1)
            square = (weyl_vector.T * root_basis_gram * weyl_vector)[0]
            assert square > 0
            weyl_vector_square = str(square)

        if dimension == 2:
            f_vector = [len(rays), facets]
        else:
            edges = len(finite_pairs)
            assert len(rays) - edges + facets == 2
            f_vector = [len(rays), edges, facets]

        reports.append(
            {
                "key": name,
                "rank": rank,
                "hyperbolic_dimension": dimension,
                "facets": facets,
                "f_vector": f_vector,
                "ns_determinant": int(ns_gram.det()),
                "wall_gram_rank": int(wall_gram.rank()),
                "wall_gram_inertia": list(_signature(wall_gram)),
                "integral_ns_change_of_basis_determinant": int(
                    ns_to_wall_basis.det()
                ),
                "walls_span_ns_integrally": True,
                "finite_pair_counts": finite_pair_counts,
                "all_extreme_rays_finite": True,
                "rational_lattice_weyl_vector": has_weyl_vector,
                "rational_weyl_vector_square": weyl_vector_square,
                "geometric_model": case["geometric_model"],
            }
        )

    assert sum(report["hyperbolic_dimension"] == 2 for report in reports) == 6
    assert sum(report["hyperbolic_dimension"] == 3 for report in reports) == 2
    assert len(canonical_wall_grams) == len(reports) == 8
    assert sum(report["rational_lattice_weyl_vector"] for report in reports) == 5
    assert {
        report["key"]
        for report in reports
        if not report["rational_lattice_weyl_vector"]
    } == {"S4", "S6", "L27"}

    return {
        "theorem": (
            "The compact complete norm-two root chambers of even integral "
            "hyperbolic lattices of rank at least three have exactly eight "
            "marked lattice/chamber types: six in H2 and two in H3."
        ),
        "classification_input": (
            "Nikulin rank-three and Vinberg rank-four exhaustion; matrix data "
            "from Roulleau arXiv:1909.01909v5"
        ),
        "classification_is_cited_not_recomputed": True,
        "source_pins": {
            "arxiv_source_tar_sha256": SOURCE_TAR_SHA256,
            "main_tex_sha256": SOURCE_TEX_SHA256,
        },
        "case_count": len(reports),
        "dimension_histogram": {"2": 6, "3": 2},
        "pairwise_distinct_unlabelled_wall_grams": True,
        "weyl_vector_count": 5,
        "new_beyond_weyl_vector_keys": ["S4", "S6", "L27"],
        "cases": reports,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-source",
        action="store_true",
        help="verify the pinned Roulleau source tarball and main TeX file",
    )
    parser.add_argument(
        "--source-tar",
        type=Path,
        default=Path("/private/tmp/roulleau1909.tar"),
    )
    parser.add_argument(
        "--source-tex",
        type=Path,
        default=Path(
            "/private/tmp/roulleau1909v5/GeomericModelsK3FiniteAuteGroup.tex"
        ),
    )
    args = parser.parse_args()

    report = certify()
    if args.check_source:
        report["source_archive_check"] = check_source_archive(
            args.source_tar, args.source_tex
        )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
