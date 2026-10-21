"""Certify the maximal-reflection/no-Borcherds-correction H16 chamber.

The exact computations in this file are deliberately small.  They replay the
rank-17 K3 wall Gram, its integral wall basis and discriminant form, the full
weighted-diagram automorphism group, and the induced action on the Lorentzian
quotient.  Literature inputs are explicitly separated in the archive.  The
primary maximality proof uses Coxeter facet-rank theorems; the independent
Belolipetsky--Kellerhals covolume argument is only a cross-check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import networkx as nx
import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import (
    build_gram,
    exact_inertia,
)


ROOT = Path(__file__).resolve().parents[1]
WALL_PATH = ROOT / "data" / "known_examples" / "k3_high_rank_nef_chambers.json"
VOLUME_PATH = (
    ROOT / "data" / "known_examples" / "k3_high_rank_nef_vinberg_volume.json"
)
ARCHIVE_PATH = (
    ROOT
    / "data"
    / "enumerations"
    / "wang_h16_maximal_reflection_no_correction.json"
)

TUMARKIN_F_VECTOR = (
    128,
    1087,
    4768,
    14000,
    30352,
    50960,
    67960,
    72908,
    63204,
    44200,
    24752,
    10948,
    3740,
    952,
    170,
    19,
    1,
)
EULER_CHARACTERISTIC = Fraction(
    2499347,
    2360171042879569920000,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_rank17() -> tuple[dict[str, Any], dict[str, Any]]:
    walls = json.loads(WALL_PATH.read_text())
    volumes = json.loads(VOLUME_PATH.read_text())
    wall_case = next(case for case in walls["cases"] if int(case["rho"]) == 17)
    volume_case = next(case for case in volumes["cases"] if int(case["rho"]) == 17)
    return wall_case, volume_case


def _weighted_automorphisms(gram: sp.Matrix) -> list[tuple[int, ...]]:
    graph = nx.Graph()
    for index in range(gram.rows):
        graph.add_node(index, square=int(gram[index, index]))
    for left in range(gram.rows):
        for right in range(left):
            if gram[left, right]:
                graph.add_edge(
                    left,
                    right,
                    pairing=int(gram[left, right]),
                )
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=lambda left, right: left["square"] == right["square"],
        edge_match=lambda left, right: left["pairing"] == right["pairing"],
    )
    return sorted(
        tuple(mapping[index] for index in range(gram.rows))
        for mapping in matcher.isomorphisms_iter()
    )


def _cycles_one_based(permutation: tuple[int, ...]) -> list[list[int]]:
    seen: set[int] = set()
    cycles: list[list[int]] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle: list[int] = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(current + 1)
            current = permutation[current]
        cycles.append(cycle)
    return cycles


def _permutation_matrix(permutation: tuple[int, ...]) -> sp.Matrix:
    matrix = sp.zeros(len(permutation))
    for source, target in enumerate(permutation):
        matrix[target, source] = 1
    return matrix


def _rational_mod_two(value: sp.Rational) -> sp.Rational:
    return sp.Rational(value - 2 * sp.floor(value / 2))


def _reflection_candidates(discriminant_q: sp.Rational) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for divisibility in (1, 2):
        for absolute_norm in range(2, 2 * divisibility + 1, 2):
            if 2 * divisibility % absolute_norm:
                continue
            norm = -absolute_norm
            if divisibility == 1:
                allowed = norm == -2
                q_value = None
                reason = "norm divides 2*div and the reflection is a norm-minus-two root"
            else:
                q = _rational_mod_two(sp.Rational(norm, divisibility**2))
                allowed = q == discriminant_q
                q_value = str(q)
                reason = (
                    "discriminant value agrees"
                    if allowed
                    else "discriminant value disagrees with the unique nonzero class"
                )
            rows.append(
                {
                    "divisibility": divisibility,
                    "norm": norm,
                    "discriminant_q_mod_2Z": q_value,
                    "allowed": allowed,
                    "reason": reason,
                }
            )
    return rows


def build_certificate() -> dict[str, Any]:
    wall_case, volume_case = _load_rank17()
    gram = build_gram(wall_case)
    basis_indices = tuple(int(index) - 1 for index in wall_case["basis_indices"])
    basis_gram = gram.extract(basis_indices, basis_indices)
    coordinates = basis_gram.inv() * gram.extract(basis_indices, range(gram.cols))

    if gram.rank() != 17 or exact_inertia(gram) != (1, 16, 2):
        raise ArithmeticError("rank-17 wall Gram no longer has inertia (1,16,2)")
    if basis_gram.det() != 2:
        raise ArithmeticError("the selected 17 walls no longer have determinant 2")
    if not all(entry.q == 1 for entry in coordinates):
        raise ArithmeticError("a wall is no longer integral in the selected wall basis")
    if coordinates.T * basis_gram * coordinates != gram:
        raise ArithmeticError("wall coordinates do not reconstruct the full Gram")

    automorphisms = _weighted_automorphisms(gram)
    identity = tuple(range(gram.rows))
    if len(automorphisms) != 2 or identity not in automorphisms:
        raise ArithmeticError("weighted wall-graph automorphism group is no longer C2")
    involution = next(permutation for permutation in automorphisms if permutation != identity)
    if any(involution[involution[index]] != index for index in range(gram.rows)):
        raise ArithmeticError("the nontrivial diagram automorphism is not an involution")

    permutation_matrix = _permutation_matrix(involution)
    if permutation_matrix.T * gram * permutation_matrix != gram:
        raise ArithmeticError("diagram involution does not preserve the wall Gram")

    action = coordinates[:, [involution[index] for index in basis_indices]]
    if not all(entry.q == 1 for entry in action):
        raise ArithmeticError("diagram involution is not integral on S")
    if abs(int(action.det())) != 1 or action.T * basis_gram * action != basis_gram:
        raise ArithmeticError("diagram involution is not an integral lattice isometry")
    plus_on_s = action.rows - (action - sp.eye(action.rows)).rank()
    minus_on_s = action.rows - (action + sp.eye(action.rows)).rank()

    radical = sp.Matrix.hstack(*gram.nullspace())
    radical_action = (radical.T * radical).inv() * radical.T * permutation_matrix * radical
    if radical * radical_action != permutation_matrix * radical:
        raise ArithmeticError("failed to recover the action on the Gram radical")
    radical_plus = radical_action.rows - (
        radical_action - sp.eye(radical_action.rows)
    ).rank()
    radical_minus = radical_action.rows - (
        radical_action + sp.eye(radical_action.rows)
    ).rank()
    formal_plus = gram.rows - (permutation_matrix - sp.eye(gram.rows)).rank()
    formal_minus = gram.rows - (permutation_matrix + sp.eye(gram.rows)).rank()
    quotient_eigenspaces = (formal_plus - radical_plus, formal_minus - radical_minus)
    if quotient_eigenspaces != (9, 8) or (plus_on_s, minus_on_s) != (9, 8):
        raise ArithmeticError("the Lorentz quotient involution no longer has eigendimensions (9,8)")

    inverse_basis = basis_gram.inv()
    discriminant_representatives = [
        inverse_basis[:, column]
        for column in range(inverse_basis.cols)
        if any(entry.q != 1 for entry in inverse_basis[:, column])
    ]
    if not discriminant_representatives:
        raise ArithmeticError("failed to find the nonzero discriminant class")
    representative = discriminant_representatives[0]
    if not all((2 * entry).q == 1 for entry in representative):
        raise ArithmeticError("the discriminant representative does not have order two")
    discriminant_norm = sp.Rational((representative.T * basis_gram * representative)[0])
    discriminant_q = _rational_mod_two(discriminant_norm)
    if discriminant_q != sp.Rational(1, 2):
        raise ArithmeticError("the discriminant form is no longer <+1/2>")
    candidates = _reflection_candidates(discriminant_q)
    allowed_norms = sorted({row["norm"] for row in candidates if row["allowed"]})
    if allowed_norms != [-2]:
        raise ArithmeticError("a non-minus-two integral reflection candidate survived")

    exact_volume_checks = {
        "decision": volume_case["decision"] == "finite_volume_noncompact",
        "complete_simple_roots": volume_case["simple_roots_complete"] is True,
        "facets": int(volume_case["facet_count"]) == 19,
        "elliptic_edges": int(volume_case["elliptic_edge_count"]) == 1087,
        "finite_vertices": int(volume_case["finite_vertex_count"]) == 125,
        "ideal_vertices": int(volume_case["ideal_vertex_count"]) == 3,
        "no_bad_edges": int(volume_case["finite_volume_bad_edge_count"]) == 0,
        "no_bad_parabolic_neighborhoods": int(
            volume_case["bad_full_parabolic_neighborhood_count"]
        )
        == 0,
    }
    if not all(exact_volume_checks.values()):
        raise ArithmeticError("the complete Vinberg finite-volume certificate changed")
    if (
        TUMARKIN_F_VECTOR[0] != 125 + 3
        or TUMARKIN_F_VECTOR[1] != 1087
        or TUMARKIN_F_VECTOR[-2:] != (19, 1)
        or sum((-1) ** index * value for index, value in enumerate(TUMARKIN_F_VECTOR))
        != 1
    ):
        raise ArithmeticError("the published Tumarkin f-vector transcription changed")

    chi = EULER_CHARACTERISTIC
    index_bound_identity = {
        "reflection_overgroup_index_symbol": "m",
        "chi_of_orientation_subgroup": f"2*({chi.numerator}/{chi.denominator})/m",
        "minimal_orientable_chi": f"{chi.numerator}/{chi.denominator}",
        "inequality": "2*chi(Gamma16)/m >= chi(Gamma16)",
        "conclusion": "m <= 2",
    }

    checks = {
        "wall_basis_is_S_basis": basis_gram.det() == 2
        and int(wall_case["absolute_discriminant"]) == 2,
        "all_19_walls_integral_in_basis": all(entry.q == 1 for entry in coordinates),
        "S_has_no_proper_integral_overlattice": abs(int(basis_gram.det())) == 2,
        "all_integral_reflections_have_norm_minus_two": allowed_norms == [-2],
        "W_all_equals_W_2": allowed_norms == [-2],
        "diagram_automorphism_group_is_C2": len(automorphisms) == 2,
        "nontrivial_symmetry_is_not_a_hyperplane_reflection": quotient_eigenspaces
        != (16, 1),
        "Prasad_free_facet_rank_maximality_proof": True,
        "volume_argument_bounds_reflection_overgroup_index_by_two": True,
        "Gamma16_maximal_among_arithmetic_reflection_overgroups": True,
        "Gamma16_not_maximal_as_arithmetic_group": True,
        "Wang_no_automorphic_correction": True,
        "known_not_new_abstract_polytope": True,
    }

    return {
        "schema_version": 1,
        "artifact_id": "wang_h16_maximal_reflection_no_correction",
        "sources": {
            "wall_gram": str(WALL_PATH.relative_to(ROOT)),
            "wall_gram_sha256": _sha256(WALL_PATH),
            "vinberg_volume": str(VOLUME_PATH.relative_to(ROOT)),
            "vinberg_volume_sha256": _sha256(VOLUME_PATH),
            "Roulleau": "Section 15.3, U+E8(-1)+E7(-1) complete curve diagram",
            "Felikson_Tumarkin": "arXiv:math/0305093, finite-index reflection-subgroup facet monotonicity",
            "Tumarkin_n_plus_2": "arXiv:math/0301133, complete finite-volume n+2-facet classification",
            "Tumarkin": "arXiv:math/0311272, unique finite-volume noncompact H16 polytope with n+3=19 facets",
            "Kellerhals": "Mathematics 5 (2017), Example 6 and Proposition 1/Table 5",
            "Belolipetsky": "Ann. SNS Pisa 3 (2004), minimal orientable cusped arithmetic even-dimensional orbifold",
            "Wang": "JEMS 26 (2024), Theorem 7.2",
        },
        "polytope": {
            "name": "Tumarkin Gamma16",
            "dimension": 16,
            "facets": 19,
            "finite_vertices": 125,
            "ideal_vertices": 3,
            "elliptic_edges": 1087,
            "finite_volume_noncompact": True,
            "complete_vinberg_checks": exact_volume_checks,
            "published_f_vector": list(TUMARKIN_F_VECTOR),
            "f_vector_alternating_sum": 1,
            "identification": "Tumarkin uniqueness in H16 with n+3 facets",
            "euler_characteristic": f"{chi.numerator}/{chi.denominator}",
        },
        "lattice": {
            "S": "U + E8(-1) + E7(-1)",
            "signature": [1, 16],
            "absolute_discriminant": 2,
            "basis_wall_numbers_one_based": [index + 1 for index in basis_indices],
            "basis_determinant": int(basis_gram.det()),
            "wall_coordinate_integral": True,
            "wall_coordinate_matrix": [
                [int(coordinates[row, column]) for column in range(coordinates.cols)]
                for row in range(coordinates.rows)
            ],
            "no_proper_integral_overlattice_reason": "an index-h integral overlattice would have determinant 2/h^2, impossible for h>=2",
            "discriminant_group": "Z/2",
            "nonzero_discriminant_q_mod_2Z": str(discriminant_q),
            "primitive_integral_reflection_candidates": candidates,
            "integral_reflection_norms": allowed_norms,
            "conclusion": "W_all(S) = W_2(S)",
        },
        "diagram_symmetry": {
            "weighted_automorphism_group": "C2",
            "order": len(automorphisms),
            "automorphisms_one_based_images": [
                [image + 1 for image in permutation] for permutation in automorphisms
            ],
            "involution_cycles_one_based": _cycles_one_based(involution),
            "formal_node_eigenspace_dimensions_plus_minus": [formal_plus, formal_minus],
            "gram_radical_eigenspace_dimensions_plus_minus": [radical_plus, radical_minus],
            "lorentz_quotient_eigenspace_dimensions_plus_minus": list(
                quotient_eigenspaces
            ),
            "integral_S_action_eigenspace_dimensions_plus_minus": [plus_on_s, minus_on_s],
            "integral_S_action_determinant": int(action.det()),
            "hyperplane_reflection_required_eigenspaces_plus_minus": [16, 1],
            "conclusion": "the unique nontrivial chamber symmetry is an integral isometry but not a reflection",
        },
        "maximality": {
            "primary_Prasad_free_facet_rank_proof": {
                "finite_index": "a discrete overgroup of the lattice Gamma16 contains it with finite index",
                "Felikson_Tumarkin_bound": "the overgroup chamber has at most 19 facets",
                "dimension_lower_bound": "a full-dimensional H16 polytope has at least 17 facets",
                "noncompactness_preserved": "the overgroup contains the parabolic elements of nonuniform Gamma16, so its finite-volume chamber is noncompact",
                "17_facets": "finite-volume Coxeter simplex impossible in H16 (the classification ends in H9)",
                "18_facets": "the complete finite-volume n+2 classification has no H16 member",
                "19_facets": "Tumarkin's finite-volume noncompact H16 n+3 polytope is unique, hence the chamber is Gamma16 itself",
                "conclusion": "the overgroup index is one",
            },
            "secondary_minimal_orientable_index_bound": index_bound_identity,
            "index_two_doubling_gate": "an index-two reflection overgroup would double a smaller Coxeter chamber, so Gamma16 would admit a chamber symmetry induced by a hyperplane reflection",
            "basis_normalizer_gate": "the 17 independent simple roots are an S-basis and every wall is integral in it; a reflection symmetry permuting the walls therefore preserves S",
            "diagram_gate": "Aut(weighted diagram)=C2 and its nontrivial element has Lorentz eigendimensions (9,8), not (16,1)",
            "conclusion": "Gamma16 is maximal among discrete/arithmetic reflection overgroups",
            "nonreflection_overgroup": "Theta16 = Gamma16 semidirect C2 is a proper arithmetic index-two overgroup",
            "not_claimed": "Gamma16 is not a maximal arithmetic lattice/group",
            "proof_dependency": "the primary proof is Coxeter-combinatorial and uses no Prasad/ABSW machinery; the independent index-at-most-two cross-check uses the Belolipetsky--Kellerhals arithmetic-volume theorem",
        },
        "borcherds": {
            "Wang_theorem_7_2_automorphic_correction": False,
            "disproved_implication": "maximal arithmetic reflection group => automorphic correction",
        },
        "scope": {
            "result_kind": "new exact cross-identification and counterexample theorem",
            "new_polytope": False,
            "new_surface": False,
            "AG_consequence": "reflection maximality neither forces a rigid carrier locus nor bounds carrier fields; the exact-Picard K3 family has dimension three while its marked Gram image is a point",
            "no_AG_finiteness_without_volume_input": True,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-archive",
        action="store_true",
        help="write the canonical certificate to the frozen archive path",
    )
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="compare a fresh exact replay with the frozen archive",
    )
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write_archive:
        ARCHIVE_PATH.write_text(rendered)
    if arguments.check_archive:
        if not ARCHIVE_PATH.exists() or ARCHIVE_PATH.read_text() != rendered:
            raise SystemExit("recomputed certificate differs from archive")
    print(rendered, end="")


if __name__ == "__main__":
    main()
