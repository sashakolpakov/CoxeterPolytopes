"""Match the Shimada-derived H14 chamber to Alexeev's (15,7,1) chamber.

The comparison is deliberately stronger than matching dimensions, facet
counts, or lattice invariants.  We reconstruct the integral Gram matrix of
the 60 outward primitive mirror normals from the fixed-lattice certificate,
parse the 60 by 60 matrix distributed with Alexeev's arXiv source, and find a
colour-preserving isomorphism of their weighted complete graphs.  The
resulting permutation is then checked entry by entry.

Alexeev source attachment (arXiv:2209.09110v4):
``all_cases_in_2el_paper.txt``.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import networkx as nx
import sympy as sp

from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import (
    OUTPUT as FIXED_LATTICE_CERTIFICATE,
)
from scripts.independent_screen_shimada_disc16_d5_h14 import ROOT


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_alexeev_exact_match.json"
)
DEFAULT_ALEXEEV_MATRIX_SOURCE = Path(
    "/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt"
)
MATRIX_KEY = (15, 7, 1)
EXPECTED_MATRIX_SHA256 = (
    "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
)
EXPECTED_SOURCE_TAR_SHA256 = (
    "10e4002f6518f8bf6d669946775bbe42731f69277a52847c5602fa60ea33c2a6"
)


def _matrix(payload: list[list[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(entry) for entry in row] for row in payload])


def _parse_alexeev_matrix(path: Path) -> sp.Matrix:
    key_text = ",".join(str(value) for value in MATRIX_KEY)
    prefix = f"coxeter_matrix[({key_text})] = matrix("
    matching = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(prefix)
    ]
    if len(matching) != 1:
        raise ArithmeticError(
            f"expected one Alexeev matrix for {MATRIX_KEY}, found {len(matching)}"
        )
    line = matching[0]
    if not line.endswith(")"):
        raise ArithmeticError("Alexeev matrix line has an unexpected format")
    payload = ast.literal_eval(line[len(prefix) : -1])
    matrix = sp.Matrix(payload)
    if matrix.shape != (60, 60) or matrix != matrix.T:
        raise ArithmeticError("Alexeev's (15,7,1) matrix is not symmetric 60 by 60")
    return matrix


def _shimada_oriented_roots_and_gram(
    certificate: dict[str, Any],
) -> tuple[sp.Matrix, sp.Matrix]:
    fixed_gram = _matrix(certificate["saturated_fixed_lattice"]["gram"])
    roots: list[sp.Matrix] = []
    for record in certificate["mirrors"]["records"]:
        root = sp.Matrix(record["primitive_fixed_coordinates"])
        scale = sp.Rational(record["presented_equals_scale_times_primitive"])
        # The presented normal has the chamber orientation.  Retain that
        # orientation while replacing it by a primitive integral vector.
        if scale < 0:
            root = -root
        roots.append(root)
    root_rows = sp.Matrix.vstack(*(root.T for root in roots))
    gram = root_rows * fixed_gram * root_rows.T
    if gram.shape != (60, 60) or gram != gram.T:
        raise ArithmeticError("the Shimada mirror Gram matrix is malformed")
    return root_rows, gram


def _node_signature(matrix: sp.Matrix, index: int) -> tuple[Any, ...]:
    interactions = Counter(
        (int(matrix[other, other]), int(matrix[index, other]))
        for other in range(matrix.rows)
        if other != index
    )
    return (
        int(matrix[index, index]),
        tuple(sorted((square, product, count) for (square, product), count in interactions.items())),
    )


def _weighted_graph(matrix: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    for index in range(matrix.rows):
        graph.add_node(index, signature=_node_signature(matrix, index))
    # Include zero products too.  This is a complete edge-coloured graph, so
    # an isomorphism certifies every matrix entry rather than merely the
    # incidence graph of nonorthogonal mirrors.
    for left in range(matrix.rows):
        for right in range(left + 1, matrix.cols):
            graph.add_edge(left, right, product=int(matrix[left, right]))
    return graph


def _histogram(values: list[int]) -> dict[str, int]:
    return {str(value): count for value, count in sorted(Counter(values).items())}


def build_certificate(
    fixed_lattice_path: Path = FIXED_LATTICE_CERTIFICATE,
    alexeev_matrix_path: Path = DEFAULT_ALEXEEV_MATRIX_SOURCE,
) -> dict[str, Any]:
    fixed_certificate = json.loads(fixed_lattice_path.read_text(encoding="utf-8"))
    alexeev_bytes = alexeev_matrix_path.read_bytes()
    alexeev_sha256 = hashlib.sha256(alexeev_bytes).hexdigest()
    if alexeev_sha256 != EXPECTED_MATRIX_SHA256:
        raise ArithmeticError(
            "Alexeev matrix attachment hash changed: "
            f"expected {EXPECTED_MATRIX_SHA256}, got {alexeev_sha256}"
        )

    shimada_root_rows, shimada_gram = _shimada_oriented_roots_and_gram(
        fixed_certificate
    )
    alexeev_gram = _parse_alexeev_matrix(alexeev_matrix_path)
    if any(shimada_gram[left, right] < 0 for left in range(60) for right in range(60) if left != right):
        raise ArithmeticError("the Shimada chamber orientations give a negative off-diagonal product")
    if any(alexeev_gram[left, right] < 0 for left in range(60) for right in range(60) if left != right):
        raise ArithmeticError("Alexeev's chamber matrix has a negative off-diagonal product")

    shimada_graph = _weighted_graph(shimada_gram)
    alexeev_graph = _weighted_graph(alexeev_gram)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        shimada_graph,
        alexeev_graph,
        node_match=nx.algorithms.isomorphism.categorical_node_match(
            "signature", None
        ),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("product", None),
    )
    if not matcher.is_isomorphic():
        raise ArithmeticError(
            "the Shimada and Alexeev weighted mirror diagrams are not isomorphic"
        )
    mapping = {int(left): int(right) for left, right in matcher.mapping.items()}
    if set(mapping) != set(range(60)) or set(mapping.values()) != set(range(60)):
        raise ArithmeticError("graph matcher did not return a permutation")
    if any(
        shimada_gram[left, right]
        != alexeev_gram[mapping[left], mapping[right]]
        for left in range(60)
        for right in range(60)
    ):
        raise ArithmeticError("the proposed permutation fails entrywise verification")

    main_indices = [
        index
        for index, record in enumerate(fixed_certificate["mirrors"]["records"])
        if record["primitive_square"] == -2
        and record["primitive_divisibility"] == 1
    ]
    main_basis = shimada_root_rows[main_indices, :]
    main_gram = shimada_gram.extract(main_indices, main_indices)
    main_graph = nx.Graph()
    main_graph.add_nodes_from(range(15))
    main_graph.add_edges_from(
        (left, right)
        for left in range(15)
        for right in range(left + 1, 15)
        if main_gram[left, right] == 1
    )
    k5_subdivision = nx.Graph()
    k5_subdivision.add_nodes_from(range(5))
    next_node = 5
    for left in range(5):
        for right in range(left + 1, 5):
            k5_subdivision.add_edge(left, next_node)
            k5_subdivision.add_edge(next_node, right)
            next_node += 1
    main_is_k5_subdivision = nx.is_isomorphic(main_graph, k5_subdivision)
    if not (
        len(main_indices) == 15
        and abs(int(main_basis.det())) == 1
        and main_gram.det() == 128
        and main_is_k5_subdivision
    ):
        raise ArithmeticError("the 15 divisibility-one roots are not Alexeev's main basis")

    shimada_diagonal = [int(shimada_gram[index, index]) for index in range(60)]
    alexeev_diagonal = [int(alexeev_gram[index, index]) for index in range(60)]
    shimada_off_diagonal = [
        int(shimada_gram[left, right])
        for left in range(60)
        for right in range(left + 1, 60)
    ]
    alexeev_off_diagonal = [
        int(alexeev_gram[left, right])
        for left in range(60)
        for right in range(left + 1, 60)
    ]
    permutation_one_based = [mapping[index] + 1 for index in range(60)]

    return {
        "schema": "shimada_disc16_d5_h14_alexeev_exact_match_v1",
        "sources": {
            "shimada_fixed_lattice_certificate": str(
                fixed_lattice_path.relative_to(ROOT)
            ),
            "alexeev": {
                "author": "Valery Alexeev",
                "title": (
                    "Reflective hyperbolic 2-elementary lattices, K3 surfaces "
                    "and hyperkahler varieties"
                ),
                "arxiv": "2209.09110v4",
                "arxiv_url": "https://arxiv.org/abs/2209.09110",
                "matrix_attachment": "all_cases_in_2el_paper.txt",
                "matrix_attachment_sha256": alexeev_sha256,
                "source_tar_sha256": EXPECTED_SOURCE_TAR_SHA256,
                "matrix_key": list(MATRIX_KEY),
            },
        },
        "intrinsic_lattice": {
            "two_elementary_invariants": {"r": 15, "a": 7, "delta": 1},
            "signature": [1, 14],
            "determinant": 128,
            "same_lattice_by_indefinite_two_elementary_uniqueness": True,
        },
        "weighted_gram_matrix_comparison": {
            "size": 60,
            "shimada_diagonal_histogram": _histogram(shimada_diagonal),
            "alexeev_diagonal_histogram": _histogram(alexeev_diagonal),
            "shimada_off_diagonal_histogram": _histogram(shimada_off_diagonal),
            "alexeev_off_diagonal_histogram": _histogram(alexeev_off_diagonal),
            "complete_edge_coloured_graph_isomorphic": True,
            "entrywise_equal_after_permutation": True,
            "shimada_facet_to_alexeev_root_permutation_one_based": permutation_one_based,
        },
        "alexeev_main_root_basis_comparison": {
            "shimada_divisibility_one_minus_two_facet_indices_one_based": [
                index + 1 for index in main_indices
            ],
            "mapped_alexeev_root_indices_one_based": [
                mapping[index] + 1 for index in main_indices
            ],
            "count": 15,
            "coordinate_basis_determinant": int(main_basis.det()),
            "gram_determinant": int(main_gram.det()),
            "subdiagram": "K5^(2)",
            "subdiagram_is_edge_subdivision_of_K5": main_is_k5_subdivision,
            "note": (
                "These are exactly the 15 (-2)-roots of divisibility one. "
                "They form an integral basis, as do Alexeev's 15 main roots."
            ),
        },
        "independent_finite_volume_cross_check": {
            "shimada_exact_ideal_extreme_rays": 1027,
            "alexeev_maximal_parabolic_subdiagrams": 1027,
            "counts_equal": True,
            "note": (
                "The equality is corroborating evidence; the entrywise Gram-matrix "
                "match is the decisive chamber identification."
            ),
        },
        "decision": {
            "same_intrinsic_lattice": True,
            "same_weighted_coxeter_diagram": True,
            "same_hyperbolic_chamber_up_to_isometry": True,
            "new_hyperbolic_example": False,
            "novelty_status": (
                "independent reconstruction of Alexeev's published (15,7,1) "
                "60-wall finite-volume noncompact chamber"
            ),
            "surface_lineage_novelty": (
                "not established by this comparison; the Shimada-derived route "
                "may be a new realization, but it is not a new chamber"
            ),
        },
        "checks": {
            "primary_attachment_hash_exact": True,
            "both_matrices_symmetric_60_by_60": True,
            "both_off_diagonals_nonnegative": True,
            "weighted_graph_isomorphism_found": True,
            "all_3600_entries_verified_after_permutation": True,
            "divisibility_one_main_roots_form_unimodular_coordinate_basis": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-lattice", type=Path, default=FIXED_LATTICE_CERTIFICATE)
    parser.add_argument(
        "--alexeev-matrices", type=Path, default=DEFAULT_ALEXEEV_MATRIX_SOURCE
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(arguments.fixed_lattice, arguments.alexeev_matrices)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(json.dumps(result["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
