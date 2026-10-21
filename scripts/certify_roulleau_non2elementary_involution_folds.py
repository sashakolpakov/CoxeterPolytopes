#!/usr/bin/env python3
"""Certify two non-2-elementary Roulleau divisor-chamber folds.

The source chambers are the complete ``(-2)``-curve configurations on the
K3 Neron--Severi lattices ``U+A2+E8`` and ``U+E8+A3`` in Roulleau's atlas.
In each graph the unique nontrivial symmetry exchanges two terminal curves.
This file checks the graph transcription, reconstructs the integral NS
action, takes its saturated fixed lattice, folds every wall orbit, and applies
Vinberg's exact edge-endpoint criterion to the resulting chamber.

This is deliberately a divisor-lattice construction.  No assertion that the
graph isometry is induced by an automorphism of a surface is needed or made.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path
import tarfile
from typing import Any, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.certify_alexeev_remaining_stable_overlattices import _q_histogram
from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _parse_coxiter,
    _permutation as _relation_permutation,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import _parse_coxiter_text
from scripts.search_k3_high_rank_involution_folds import (
    _coxeter_audit,
    _cycles,
    _geometry_certificate,
    _matrix_payload,
    _weighted_graph_automorphisms,
)
from scripts.verify_k3_quotient_nef_slices import _orbit_sum_matrix
from coxeter_polytopes.search.u_d10_vinberg import (
    lattice_gram as u_d10_lattice_gram,
    pair as u_d10_pair,
    run_search as run_u_d10_vinberg,
)
from coxeter_polytopes.search.u_e8_a1_vinberg import (
    lattice_gram as u_e8_a1_lattice_gram,
    pair as u_e8_a1_pair,
    run_search as run_u_e8_a1_vinberg,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_non2elementary_involution_folds.json"
SOURCE_INVENTORY = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
COXITER_CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
DOUBLING_ARCHIVE = COXITER_CORPUS / "graphs/testing_index2.tar.gz"
SK_ROOT = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"


CASES: tuple[dict[str, Any], ...] = (
    {
        "id": "R12_UA2E8_C2_FIXED",
        "source_stock": "R12_UA2E8",
        "lattice": "U+A2(-1)+E8(-1)",
        "rho": 12,
        "wall_count": 13,
        "section_lines": [9034, 9125],
        "table_lines": [10928, 10928],
        "surface_automorphism_group": "C2",
        "surface_NS_action_tag": "n",
        "chain": [1, 2, 4, 5, 6, 7, 8, 9, 10, 11],
        "extra_edges": [(4, 3), (11, 12), (12, 13), (13, 11)],
        "swapped_labels": [12, 13],
        "dependent_coordinates": [2, 4, 3, 6, 5, 4, 3, 2, 1, 0, -1, -1],
        "expected": {
            "dimension": 10,
            "facets": 12,
            "edges": 99,
            "finite_vertices": 18,
            "cusps": 1,
            "parabolic_audits": 2,
            "fixed_determinant": 2,
            "smith": [1] * 10 + [2],
            "root_types": {"-2,1": 11, "-2,2": 1},
            "relation_histogram": {"0": 1, "2": 55, "3": 10},
            "finite_special_counts": [1, 12, 65, 210, 450, 672, 714, 540, 285, 99, 18],
            "euler": "-1/183936614400",
            "sk_list": "List10",
            "sk_entry": 3,
            "sk_symbol": "II_{10,1}([2^{1}]_{1})",
            "known_coxiter_matches": [
                "graphs/10-thierryP_1-10.coxiter",
                "graphs/Tumarkin pyramids n+2/10-tum04_12_01-k=2-l=3.coxiter",
            ],
        },
    },
    {
        "id": "R13_UE8A3_C2_FIXED",
        "source_stock": "R13_UE8A3",
        "lattice": "U+E8(-1)+A3(-1)",
        "rho": 13,
        "wall_count": 14,
        "section_lines": [9491, 9583],
        "table_lines": [10945, 10945],
        "surface_automorphism_group": "C2",
        "surface_NS_action_tag": "n",
        "chain": [1, 2, 3, 5, 6, 7, 8, 9, 10, 11],
        "extra_edges": [(3, 4), (11, 12), (12, 13), (13, 14), (14, 11)],
        "swapped_labels": [12, 14],
        "dependent_coordinates": [2, 4, 6, 3, 5, 4, 3, 2, 1, 0, -1, -1, -1],
        "expected": {
            "dimension": 11,
            "facets": 13,
            "edges": 162,
            "finite_vertices": 26,
            "cusps": 2,
            "parabolic_audits": 3,
            "fixed_determinant": -4,
            "smith": [1] * 10 + [2, 2],
            "root_types": {"-4,2": 1, "-2,1": 11, "-2,2": 1},
            "relation_histogram": {"2": 66, "3": 10, "4": 2},
            "finite_special_counts": [1, 13, 78, 285, 705, 1242, 1596, 1506, 1035, 504, 162, 26],
            "euler": "0",
            "sk_list": "List11",
            "sk_entry": 6,
            "sk_symbol": "II_{11,1}([2^{2}]_{2})",
            "known_coxiter_matches": [
                "graphs/group units qf/11-1,1,1,1,1,1,1,1,1,1,1,1.coxiter",
                "graphs/Tumarkin pyramids n+2/11-tum04_12_02.coxiter",
            ],
        },
    },
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _curve_gram(case: dict[str, Any]) -> sp.Matrix:
    size = int(case["wall_count"])
    gram = -2 * sp.eye(size)
    chain = [int(label) for label in case["chain"]]
    edges = list(zip(chain, chain[1:])) + list(case["extra_edges"])
    for left, right in edges:
        gram[left - 1, right - 1] = gram[right - 1, left - 1] = 1
    return gram


def _smith(matrix: sp.Matrix) -> list[int]:
    return [abs(int(value)) for value in smith_normal_form(matrix, domain=ZZ).diagonal()]


def _fixed_basis(action: sp.Matrix) -> sp.Matrix:
    equations = DomainMatrix.from_Matrix(action - sp.eye(action.rows)).convert_to(ZZ)
    diagonal, _, right = smith_normal_decomp(equations)
    normal = diagonal.to_Matrix()
    rank = sum(normal[index, index] != 0 for index in range(action.rows))
    basis = right.to_Matrix()[:, rank:]
    if abs(int(right.to_Matrix().det())) != 1:
        raise ArithmeticError("Smith right transformation is not unimodular")
    if (action - sp.eye(action.rows)) * basis != sp.zeros(action.rows, basis.cols):
        raise ArithmeticError("Smith fixed basis is not in the kernel")
    return basis


def _divisibility(gram: sp.Matrix, root: sp.Matrix) -> int:
    result = 0
    for value in gram * root:
        result = gcd(result, abs(int(value)))
    return result


def _integral_reflection_audit(
    fixed_gram: sp.Matrix, roots: sp.Matrix, cycles: Sequence[Sequence[int]]
) -> dict[str, Any]:
    records = []
    for column, cycle in enumerate(cycles):
        root = roots[:, column]
        norm = int((root.T * fixed_gram * root)[0])
        divisibility = _divisibility(fixed_gram, root)
        reflection = sp.eye(fixed_gram.rows) - sp.Rational(2, norm) * root * (
            root.T * fixed_gram
        )
        integral = all(sp.Rational(value).q == 1 for value in reflection)
        records.append(
            {
                "wall_orbit": list(cycle),
                "fixed_lattice_coordinates": [int(value) for value in root],
                "norm": norm,
                "divisibility": divisibility,
                "reflection_is_integral": integral,
                "reflection_is_involutive_isometry": (
                    reflection * reflection == sp.eye(fixed_gram.rows)
                    and reflection.T * fixed_gram * reflection == fixed_gram
                ),
            }
        )
    return {
        "all_folded_reflections_are_integral": all(
            record["reflection_is_integral"]
            and record["reflection_is_involutive_isometry"]
            for record in records
        ),
        "checks": records,
    }


def _scan_coxiter(
    relation_matrix: list[list[str]], dimension: int, facets: int
) -> dict[str, Any]:
    ordinary_candidates = 0
    matches: list[dict[str, Any]] = []
    for path in sorted(COXITER_CORPUS.rglob("*.coxiter")):
        header = next(
            content
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if (content := line.split("#", 1)[0].strip())
        )
        header_facets, header_dimension = map(int, header.split()[:2])
        if (header_dimension, header_facets) != (dimension, facets):
            continue
        ordinary_candidates += 1
        parsed_facets, parsed_dimension, target = _parse_coxiter(path)
        if (parsed_dimension, parsed_facets) != (dimension, facets):
            raise ArithmeticError("CoxIter header changed during parse")
        permutation = _relation_permutation(relation_matrix, target)
        if permutation is not None:
            matches.append(
                {
                    "lane": "ordinary",
                    "path": str(path.relative_to(COXITER_CORPUS)),
                    "sha256": _sha256(path),
                    "permutation_zero_based": permutation,
                }
            )

    doubling_header_candidates = 0
    doubling_matches: list[dict[str, Any]] = []
    with tarfile.open(DOUBLING_ARCHIVE, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"cannot read {member.name}")
            raw = handle.read()
            lines = [
                content
                for line in raw.decode("utf-8", errors="replace").splitlines()
                if (content := line.split("#", 1)[0].strip())
            ]
            if not lines:
                continue
            header_facets, header_dimension = map(int, lines[0].split()[:2])
            if (header_dimension, header_facets) != (dimension, facets):
                continue
            doubling_header_candidates += 1
            parsed_facets, parsed_dimension, target = _parse_coxiter_text(
                raw.decode("utf-8", errors="replace")
            )
            if (parsed_dimension, parsed_facets) != (dimension, facets):
                raise ArithmeticError("doubling header changed during parse")
            permutation = _relation_permutation(relation_matrix, target)
            if permutation is not None:
                doubling_matches.append(
                    {
                        "lane": "index_two_doubling",
                        "path": member.name,
                        "sha256": hashlib.sha256(raw).hexdigest(),
                        "permutation_zero_based": permutation,
                    }
                )
    return {
        "coxiter_commit": COXITER_COMMIT,
        "ordinary_header_candidate_count": ordinary_candidates,
        "ordinary_exact_matches": matches,
        "doubling_archive_sha256": _sha256(DOUBLING_ARCHIVE),
        "doubling_header_candidate_count": doubling_header_candidates,
        "doubling_exact_matches": doubling_matches,
    }


def _exact_gram_permutation(
    source: Sequence[Sequence[int]], target: Sequence[Sequence[int]]
) -> list[int] | None:
    def graph(matrix: Sequence[Sequence[int]]) -> nx.Graph:
        result = nx.Graph()
        for index in range(len(matrix)):
            result.add_node(index, norm=int(matrix[index][index]))
        for left in range(len(matrix)):
            for right in range(left + 1, len(matrix)):
                result.add_edge(left, right, product=int(matrix[left][right]))
        return result

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph(source),
        graph(target),
        node_match=nx.algorithms.isomorphism.categorical_node_match("norm", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("product", None),
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        return None
    permutation = [mapping[index] for index in range(len(source))]
    if any(
        source[left][right] != target[permutation[left]][permutation[right]]
        for left in range(len(source))
        for right in range(len(source))
    ):
        raise ArithmeticError("Gram matcher did not preserve every entry")
    return permutation


def _analyze_case(case: dict[str, Any], stock: dict[str, Any]) -> dict[str, Any]:
    expected = case["expected"]
    full_gram = _curve_gram(case)
    size = full_gram.rows
    basis_indices = list(range(size - 1))
    basis_gram = full_gram.extract(basis_indices, basis_indices)
    wall_coordinates = [
        basis_gram.inv() * full_gram.extract(basis_indices, [wall])
        for wall in range(size)
    ]
    if [int(value) for value in wall_coordinates[-1]] != case["dependent_coordinates"]:
        raise ArithmeticError(f"dependent wall transcription failed for {case['id']}")
    if any(sp.Rational(value).q != 1 for wall in wall_coordinates for value in wall):
        raise ArithmeticError("a source wall is not integral in the claimed NS basis")

    group, search_nodes = _weighted_graph_automorphisms(full_gram)
    nonidentity = [element for element in group if element != tuple(range(size))]
    if len(group) != 2 or len(nonidentity) != 1:
        raise ArithmeticError("source graph does not have the expected unique involution")
    permutation = nonidentity[0]
    expected_swap = {int(label) - 1 for label in case["swapped_labels"]}
    moved = {index for index, image in enumerate(permutation) if index != image}
    if moved != expected_swap:
        raise ArithmeticError("unique graph involution swaps the wrong walls")

    action = sp.Matrix.hstack(
        *(wall_coordinates[permutation[index]] for index in basis_indices)
    )
    fixed_basis = _fixed_basis(action)
    fixed_gram = fixed_basis.T * basis_gram * fixed_basis
    cycles = _cycles(permutation)
    orbit_sums = _orbit_sum_matrix(size, cycles)
    orbit_gram = orbit_sums.T * full_gram * orbit_sums
    full_wall_matrix = sp.Matrix.hstack(*wall_coordinates)
    orbit_vectors = full_wall_matrix * orbit_sums
    fixed_roots = sp.Matrix.hstack(
        *(fixed_basis.gauss_jordan_solve(orbit_vectors[:, column])[0]
          for column in range(orbit_vectors.cols))
    )
    if any(sp.Rational(value).q != 1 for value in fixed_roots):
        raise ArithmeticError("a folded wall is not integral in the fixed lattice")
    if fixed_roots.T * fixed_gram * fixed_roots != orbit_gram:
        raise ArithmeticError("fixed-lattice and orbit-sum Grams disagree")

    reflection_audit = _integral_reflection_audit(fixed_gram, fixed_roots, cycles)
    facet_indices = list(range(len(cycles)))
    coxeter = _coxeter_audit(orbit_gram, facet_indices, cycles)
    action_checks = {
        "action_is_integral": all(sp.Rational(value).q == 1 for value in action),
        "action_is_unimodular": abs(int(action.det())) == 1,
        "action_has_order_two": action * action == sp.eye(action.rows),
        "action_preserves_NS_gram": action.T * basis_gram * action == basis_gram,
        "action_permutes_every_wall": all(
            action * wall_coordinates[index] == wall_coordinates[permutation[index]]
            for index in range(size)
        ),
    }
    ns_action = {
        "basis_wall_labels": [index + 1 for index in basis_indices],
        "basis_gram_determinant": int(basis_gram.det()),
        "action_matrix_columns": _matrix_payload(action),
        "determinant": int(action.det()),
        "fixed_rank": fixed_basis.cols,
        "fixed_basis_is_saturated_by_unimodular_Smith_change": True,
        "checks": action_checks,
        "all_integral_lattice_checks_pass": all(action_checks.values()),
        "surface_realization_status": (
            "proved_for_the_general_K3: Roulleau_records_Aut=C2_and_tag_n_"
            "(hyperelliptic_involution_nontrivial_on_NS); exhaustive_graph_"
            "Aut=C2_forces_its_wall_action_to_equal_this_unique_involution"
        ),
        "surface_realization_source": {
            "general_surface_automorphism_group": stock[
                "general_surface_automorphism_group"
            ],
            "roulleau_NS_action_tag": stock["roulleau_ns_action_tag"],
            "tag_definition_tex_lines": [10985, 10990],
        },
    }
    geometry = _geometry_certificate(
        int(case["rho"]), cycles, orbit_gram, coxeter, reflection_audit, ns_action
    )

    positive_gram = [[str(-orbit_gram[i, j]) for j in range(orbit_gram.cols)]
                     for i in range(orbit_gram.rows)]
    relations = _relation_matrix(positive_gram)
    connected = _connected_masks_direct(relations)
    solomon = _euler_from_connected_components(relations, connected)
    relation_histogram = Counter(
        relations[left][right]
        for left in range(len(relations))
        for right in range(left + 1, len(relations))
    )
    root_types = Counter(
        (int((fixed_roots[:, column].T * fixed_gram * fixed_roots[:, column])[0]),
         _divisibility(fixed_gram, fixed_roots[:, column]))
        for column in range(fixed_roots.cols)
    )
    coxiter = _scan_coxiter(relations, int(expected["dimension"]), int(expected["facets"]))
    ordinary_paths = sorted(record["path"] for record in coxiter["ordinary_exact_matches"])

    if case["id"] == "R12_UA2E8_C2_FIXED":
        search = run_u_e8_a1_vinberg(40)
        intrinsic_pair = u_e8_a1_pair
        intrinsic_gram = sp.Matrix(u_e8_a1_lattice_gram())
        intrinsic_name = "U+E8(-1)+A1(-1)"
        expected_stopping_height = sp.Rational(1, 2)
    elif case["id"] == "R13_UE8A3_C2_FIXED":
        search = run_u_d10_vinberg(40)
        intrinsic_pair = u_d10_pair
        intrinsic_gram = sp.Matrix(u_d10_lattice_gram())
        intrinsic_name = "U+D10(-1)"
        expected_stopping_height = sp.Rational(8)
    else:
        raise AssertionError("unconfigured intrinsic Vinberg case")
    vinberg_positive = [
        [-intrinsic_pair(left, right) for right in search.roots]
        for left in search.roots
    ]
    candidate_positive = [
        [-int(orbit_gram[left, right]) for right in range(orbit_gram.cols)]
        for left in range(orbit_gram.rows)
    ]
    permutation = _exact_gram_permutation(candidate_positive, vinberg_positive)
    if permutation is None:
        raise ArithmeticError(
            f"folded walls do not match intrinsic {intrinsic_name} Vinberg roots"
        )
    ordered_vinberg_roots = sp.Matrix.hstack(
        *(sp.Matrix(search.roots[permutation[index]]) for index in range(len(permutation)))
    )
    independent: list[int] = []
    for column in range(fixed_roots.cols):
        if fixed_roots[:, [*independent, column]].rank() > len(independent):
            independent.append(column)
    if len(independent) != fixed_gram.rows:
        raise ArithmeticError("folded roots do not span the fixed lattice")
    isometry = (
        ordered_vinberg_roots[:, independent]
        * fixed_roots[:, independent].inv()
    )
    intrinsic_vinberg_check = (
        search.stopping.fired
        and len(search.roots) == int(expected["facets"])
        and search.levels[-1].height == expected_stopping_height
        and search.stopping.finite_vertex_count == int(expected["finite_vertices"])
        and search.stopping.cusp_count == int(expected["cusps"])
        and all(sp.Rational(value).q == 1 for value in isometry)
        and abs(int(isometry.det())) == 1
        and isometry.T * intrinsic_gram * isometry == fixed_gram
        and isometry * fixed_roots == ordered_vinberg_roots
    )
    full_reflection_lineage = {
        "status": "exact_intrinsic_Vinberg_full_chamber_conjugacy",
        "index": 1,
        "intrinsic_lattice": intrinsic_name,
        "vinberg_stopping_height": str(search.levels[-1].height),
        "vinberg_root_count": len(search.roots),
        "finite_vertex_count": search.stopping.finite_vertex_count,
        "cusp_count": search.stopping.cusp_count,
        "candidate_to_vinberg_facet_permutation_zero_based": permutation,
        "fixed_lattice_to_intrinsic_isometry": _matrix_payload(isometry),
        "isometry_determinant": int(isometry.det()),
        "wall_by_wall_signed_Gram_and_vector_conjugacy": intrinsic_vinberg_check,
        "levels": [
            {
                "height": str(level.height),
                "keys": [list(key) for key in level.keys],
                "candidate_count": level.candidate_count,
                "accepted_roots": [list(root) for root in level.accepted_roots],
            }
            for level in search.levels
        ],
    }

    list_path = SK_ROOT / expected["sk_list"]
    list_line = list_path.read_text().splitlines()[int(expected["sk_entry"])]
    checks = {
        "source_full_wall_gram_rank_and_signature": (
            full_gram.rank() == int(case["rho"])
            and geometry["orbit_gram_inertia"] == [1, int(expected["dimension"]), 1]
        ),
        "unique_source_graph_involution": len(group) == 2,
        "integral_NS_action": ns_action["all_integral_lattice_checks_pass"],
        "fixed_lattice_invariants": (
            int(fixed_gram.det()) == int(expected["fixed_determinant"])
            and _smith(fixed_gram) == expected["smith"]
            and all(int(fixed_gram[index, index]) % 2 == 0 for index in range(fixed_gram.rows))
        ),
        "all_folded_reflections_integral": reflection_audit["all_folded_reflections_are_integral"],
        "exact_coxeter_angles": coxeter["all_pairs_are_exact_coxeter_relations"],
        "vinberg_finite_volume_noncompact": (
            geometry["finite_volume_certified"] and geometry["compact"] is False
        ),
        "vinberg_exact_counts": (
            geometry["vinberg_volume_criterion"]["elliptic_edge_count"] == expected["edges"]
            and geometry["finite_vertex_count"] == expected["finite_vertices"]
            and geometry["ideal_vertex_count"] == expected["cusps"]
            and geometry["vinberg_volume_criterion"]["parabolic_neighborhood_audit_count"]
            == expected["parabolic_audits"]
        ),
        "root_type_histogram": {
            f"{norm},{divisibility}": count
            for (norm, divisibility), count in sorted(root_types.items())
        } == expected["root_types"],
        "relation_histogram": dict(sorted(relation_histogram.items()))
        == expected["relation_histogram"],
        "solomon_characteristic": (
            solomon["finite_special_subset_count_by_rank"]
            == expected["finite_special_counts"]
            and solomon["rational_euler_characteristic"] == expected["euler"]
        ),
        "SK_list_row": expected["sk_symbol"] in list_line,
        "ordinary_CoxIter_matches": ordinary_paths == sorted(expected["known_coxiter_matches"]),
        "no_index_two_doubling_match": not coxiter["doubling_exact_matches"],
        "surface_involution_realizes_unique_graph_action": (
            stock["general_surface_automorphism_group"]
            == case["surface_automorphism_group"]
            and stock["roulleau_ns_action_tag"] == case["surface_NS_action_tag"]
            and case["surface_NS_action_tag"] == "n"
            and len(group) == 2
        ),
        "intrinsic_full_reflection_chamber_lineage": intrinsic_vinberg_check,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"case checks failed for {case['id']}: {checks}")

    return {
        "id": case["id"],
        "source": {
            "stock_key": case["source_stock"],
            "lattice": case["lattice"],
            "picard_rank": case["rho"],
            "complete_minus_two_curve_count": size,
            "tex_section_lines": case["section_lines"],
            "curve_gram": _matrix_payload(full_gram),
            "basis_wall_labels": [index + 1 for index in basis_indices],
            "dependent_last_wall_coordinates": case["dependent_coordinates"],
            "completeness": "Roulleau_published_complete_curve_configuration",
            "surface_automorphism_group": case["surface_automorphism_group"],
            "surface_NS_action_tag": case["surface_NS_action_tag"],
            "surface_action_tag_definition_tex_lines": [10985, 10990],
        },
        "graph_automorphism": {
            "group_order": len(group),
            "exhaustive_backtracking_node_count": search_nodes,
            "permutation_images_one_based": [value + 1 for value in permutation],
            "cycles": cycles,
        },
        "neron_severi_action": ns_action,
        "fixed_lattice": {
            "saturated_basis_columns": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "rank": fixed_gram.rows,
            "determinant": int(fixed_gram.det()),
            "smith_invariants": _smith(fixed_gram),
            "even": all(int(fixed_gram[index, index]) % 2 == 0 for index in range(fixed_gram.rows)),
            "discriminant_q_histogram_mod_2": _q_histogram(fixed_gram),
            "classification": {
                "method": "indefinite_even_2_elementary_lattice_classification",
                "scharlau_kirschmer_list": expected["sk_list"],
                "entry_one_based": expected["sk_entry"],
                "genus_symbol": expected["sk_symbol"],
                "class_number": 1,
            },
        },
        "folded_walls": {
            "wall_orbits": cycles,
            "root_coordinate_columns": _matrix_payload(fixed_roots),
            "negative_gram": _matrix_payload(orbit_gram),
            "positive_wall_gram": positive_gram,
            "normalized_relation_matrix": relations,
            "root_type_histogram": {
                f"{norm},{divisibility}": count
                for (norm, divisibility), count in sorted(root_types.items())
            },
            "reflection_audit": reflection_audit,
        },
        "geometry": geometry,
        "solomon": solomon,
        "lineage": {
            "scharlau_kirschmer_list_row": list_line.strip(),
            "coxiter": coxiter,
            "full_reflection_chamber": full_reflection_lineage,
            "abstract_chamber_status": "known_exact_CoxIter_relation_graph_match",
            "maximal_reflection_index_status": "exact_index_one",
            "algebro_geometric_status": (
                "exact_complete_K3_divisor_chamber_fixed_section_realization; "
                "literature_priority_of_this_realization_not_yet_closed"
            ),
            "surface_automorphism_required_for_polytope_proof": False,
            "surface_automorphism_realization": (
                "proved_by_Roulleau_C2_nontrivial_NS_action_plus_unique_C2_"
                "weighted_wall_graph_action"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate() -> dict[str, Any]:
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    if not inventory["all_checks_pass"]:
        raise ArithmeticError("Roulleau source inventory is not certified")
    stock_by_key = {stock["key"]: stock for stock in inventory["stocks"]}
    tex_lines = SOURCE_TEX.read_text().splitlines()
    for case in CASES:
        stock = stock_by_key[case["source_stock"]]
        if (
            stock["wall_count"] != case["wall_count"]
            or stock["picard_rank"] != case["rho"]
            or stock["source"]["section_lines"] != case["section_lines"]
        ):
            raise ArithmeticError(f"inventory metadata changed for {case['id']}")
        start, end = case["section_lines"]
        excerpt = "\n".join(tex_lines[start - 1:end])
        if (
            stock["tex_lattice_token"] not in excerpt
            or f"contains ${case['wall_count']}$" not in excerpt
        ):
            raise ArithmeticError(f"source wall-count witness missing for {case['id']}")
    records = [_analyze_case(case, stock_by_key[case["source_stock"]]) for case in CASES]
    return {
        "schema": "coxeter.roulleau_non2elementary_involution_folds.v1",
        "artifact_id": "roulleau_non2elementary_involution_folds",
        "source": {
            "inventory": str(SOURCE_INVENTORY.relative_to(ROOT)),
            "inventory_sha256": _sha256(SOURCE_INVENTORY),
            "roulleau_tex_sha256": _sha256(SOURCE_TEX),
            "scharlau_kirschmer": {
                name: _sha256(SK_ROOT / name)
                for name in ("List10", "List11", "Gram10", "Gram11")
            },
        },
        "cases": records,
        "summary": {
            "case_count": len(records),
            "finite_volume_noncompact_count": sum(
                record["geometry"]["finite_volume_certified"]
                and record["geometry"]["compact"] is False
                for record in records
            ),
            "hyperbolic_dimensions": [
                record["geometry"]["hyperbolic_dimension"] for record in records
            ],
            "facet_counts": [record["geometry"]["facet_count"] for record in records],
            "new_abstract_chamber_count": 0,
            "new_algebro_geometric_realization_priority": "open_literature_check",
        },
        "all_checks_pass": all(record["all_checks_pass"] for record in records),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        if certificate != json.loads(OUTPUT.read_text()):
            raise SystemExit("recomputed certificate differs from archive")
    elif args.output:
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
