#!/usr/bin/env python3
"""Identify and reject the rank-16 K3 graph-involution H14 fold.

The proposed action fixes walls 1,...,15 and swaps (16,18),(17,19).
This script checks three exact facts:

* it is the unique nonidentity weighted-graph symmetry fixing 15 walls;
* its action on the 2-elementary discriminant group of
  ``U + E8(-1) + D6(-1)`` is nontrivial; and
* its folded Coxeter diagram is exactly Vinberg's published H14 chamber,
  transcribed as ``14-vinb85.coxiter`` in Guglielmetti's Appendix A.

The geometric nonrealization also uses Kondō's source theorem and Torelli,
as explained in the accompanying report.  The source statements themselves
are not replaced by computation here.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

from scripts.certify_k3_high_rank_all_involution_folds import (
    _canonical_sha256,
    _cycles,
    _orbit_sum_gram,
    exact_weighted_gram_automorphisms,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
CERTIFICATE_PATH = ROOT / "data/candidates/k3_rho16_h14_fold_identification.json"

PROPOSED_CYCLES = [[index] for index in range(1, 16)] + [[16, 18], [17, 19]]

# Guglielmetti, CoxIter, Appendix A, Figure A.1.  The third entry 0 denotes
# the bold (parallel-wall) edge in CoxIter syntax.
VINBERG_14_EDGES = (
    (1, 2, 3),
    (2, 3, 3),
    (3, 4, 3),
    (4, 5, 3),
    (5, 6, 3),
    (6, 7, 3),
    (7, 8, 3),
    (8, 9, 3),
    (9, 10, 3),
    (10, 11, 3),
    (11, 12, 3),
    (12, 13, 3),
    (13, 14, 4),
    (14, 17, 0),
    (17, 1, 4),
    (11, 16, 3),
    (3, 15, 3),
)

VINBERG_14_F_VECTOR = (
    94,
    704,
    2695,
    6825,
    12579,
    17633,
    19215,
    16425,
    11009,
    5733,
    2275,
    665,
    135,
    17,
    1,
)

# Our folded wall label -> the label in 14-vinb85.coxiter.
OUR_TO_VINBERG = {
    1: 1,
    2: 2,
    3: 3,
    4: 15,
    5: 4,
    6: 5,
    7: 6,
    8: 7,
    9: 8,
    10: 9,
    11: 10,
    12: 11,
    13: 16,
    14: 12,
    15: 13,
    16: 14,
    17: 17,
}


def _load_rank16() -> dict[str, Any]:
    artifact = json.loads(SOURCE_PATH.read_text())
    return next(case for case in artifact["cases"] if int(case["rho"]) == 16)


def _permutation_from_cycles(size: int, cycles: Sequence[Sequence[int]]) -> tuple[int, ...]:
    answer = list(range(size))
    seen: set[int] = set()
    for cycle in cycles:
        values = [int(label) - 1 for label in cycle]
        if not values or seen.intersection(values):
            raise ValueError("cycles are empty or overlap")
        seen.update(values)
        for index, source in enumerate(values):
            answer[source] = values[(index + 1) % len(values)]
    if seen != set(range(size)):
        raise ValueError("cycles do not cover all walls")
    return tuple(answer)


def _fractional_part(value: sp.Rational) -> sp.Rational:
    return sp.Rational(value - sp.floor(value))


def _coset_key(vector: sp.Matrix) -> tuple[sp.Rational, ...]:
    return tuple(_fractional_part(sp.Rational(value)) for value in vector)


def _discriminant_classes(inverse_gram: sp.Matrix) -> tuple[tuple[sp.Rational, ...], ...]:
    zero = tuple(sp.Rational(0) for _ in range(inverse_gram.rows))
    classes = {zero}
    for column in range(inverse_gram.cols):
        generator = inverse_gram[:, column]
        previous = tuple(classes)
        classes.update(_coset_key(sp.Matrix(value) + generator) for value in previous)
    return tuple(sorted(classes))


def _wall_coordinates(wall_gram: sp.Matrix, basis: Sequence[int]) -> tuple[sp.Matrix, ...]:
    basis_gram = wall_gram.extract(basis, basis)
    inverse = basis_gram.inv()
    answer = tuple(
        inverse * wall_gram.extract(basis, [wall])
        for wall in range(wall_gram.rows)
    )
    if any(any(sp.Rational(value).q != 1 for value in vector) for vector in answer):
        raise AssertionError("a published wall is nonintegral in the archived basis")
    return answer


def _isometry_matrix(
    permutation: Sequence[int], wall_coordinates: Sequence[sp.Matrix], basis: Sequence[int]
) -> sp.Matrix:
    matrix = sp.zeros(len(basis))
    for column, wall in enumerate(basis):
        matrix[:, column] = wall_coordinates[permutation[wall]]
    return matrix


def _fold_edges(orbit_gram: sp.Matrix) -> tuple[tuple[int, int, int], ...]:
    answer: list[tuple[int, int, int]] = []
    for left in range(orbit_gram.rows):
        for right in range(left + 1, orbit_gram.cols):
            label = _coxeter_label(
                sp.Rational(orbit_gram[left, left]),
                sp.Rational(orbit_gram[right, right]),
                sp.Rational(orbit_gram[left, right]),
            )
            if label == "2":
                continue
            if label == "infinity":
                encoded = 0
            elif label in {"3", "4", "6"}:
                encoded = int(label)
            else:
                raise AssertionError(f"unexpected folded edge {label}")
            answer.append((left + 1, right + 1, encoded))
    return tuple(answer)


def _unoriented_edges(edges: Sequence[tuple[int, int, int]]) -> set[tuple[int, int, int]]:
    return {(min(left, right), max(left, right), label) for left, right, label in edges}


@lru_cache(maxsize=1)
def build_certificate() -> dict[str, Any]:
    case = _load_rank16()
    wall_gram = build_gram(case)
    basis = tuple(int(label) - 1 for label in case["basis_indices"])
    basis_gram = wall_gram.extract(basis, basis)
    wall_coordinates = _wall_coordinates(wall_gram, basis)
    automorphisms = exact_weighted_gram_automorphisms(wall_gram)
    proposed = _permutation_from_cycles(wall_gram.rows, PROPOSED_CYCLES)
    fixed_counts = {
        action: sum(action[index] == index for index in range(len(action)))
        for action in automorphisms
    }

    isometry = _isometry_matrix(proposed, wall_coordinates, basis)
    if isometry.T * basis_gram * isometry != basis_gram:
        raise AssertionError("the proposed graph action is not a lattice isometry")
    discriminant_classes = _discriminant_classes(basis_gram.inv())
    discriminant_action = tuple(
        discriminant_classes.index(_coset_key(isometry * sp.Matrix(value)))
        for value in discriminant_classes
    )

    orbit_gram = _orbit_sum_gram(wall_gram, PROPOSED_CYCLES)
    our_edges = _fold_edges(orbit_gram)
    transported_edges = tuple(
        (OUR_TO_VINBERG[left], OUR_TO_VINBERG[right], label)
        for left, right, label in our_edges
    )
    diagram_matches = _unoriented_edges(transported_edges) == _unoriented_edges(
        VINBERG_14_EDGES
    )

    fixed_wall_gram = wall_gram.extract(range(15), range(15))
    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "artifact_id": "k3_rho16_h14_fold_identification",
        "certificate_kind": "exact_nonrealization_and_published_diagram_identification",
        "source": {
            "path": str(SOURCE_PATH.relative_to(ROOT)),
            "sha256": hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
            "ambient_lattice": case["lattice"],
            "ambient_rho": 16,
            "wall_count": wall_gram.rows,
        },
        "proposed_action": {
            "cycles": PROPOSED_CYCLES,
            "permutation_images_zero_based": list(proposed),
            "weighted_graph_automorphism_group": "C2 x C2",
            "weighted_graph_automorphism_count": len(automorphisms),
            "all_automorphism_cycles": [_cycles(action) for action in automorphisms],
            "fixed_wall_count": fixed_counts[proposed],
            "unique_nonidentity_action_fixing_15_walls": sum(
                action != tuple(range(wall_gram.rows)) and count == 15
                for action, count in fixed_counts.items()
            )
            == 1,
        },
        "lattice_action": {
            "basis_wall_labels": [label + 1 for label in basis],
            "basis_determinant": int(basis_gram.det()),
            "isometry_matrix_columns": [
                [int(isometry[row, column]) for row in range(isometry.rows)]
                for column in range(isometry.cols)
            ],
            "isometry_order_two": isometry * isometry == sp.eye(isometry.rows),
            "discriminant_group": "(Z/2)^2",
            "discriminant_classes_half_coordinate_supports": [
                [index + 1 for index, value in enumerate(coset) if value == Fraction(1, 2)]
                for coset in discriminant_classes
            ],
            "discriminant_action": list(discriminant_action),
            "discriminant_action_is_nontrivial": discriminant_action
            != tuple(range(len(discriminant_classes))),
            "fixed_wall_basis_determinant": int(fixed_wall_gram.det()),
            "fixed_wall_basis_inertia": list(exact_inertia(fixed_wall_gram)),
        },
        "folded_diagram": {
            "hyperbolic_dimension": orbit_gram.rank() - 1,
            "facet_count": orbit_gram.rows,
            "orbit_gram_rank": orbit_gram.rank(),
            "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
            "edges_coxiter_encoding": [list(edge) for edge in our_edges],
            "edge_label_histogram": {
                "3": sum(label == 3 for _, _, label in our_edges),
                "4": sum(label == 4 for _, _, label in our_edges),
                "infinity": sum(label == 0 for _, _, label in our_edges),
            },
        },
        "published_identification": {
            "name": "14-vinb85.coxiter",
            "attribution": (
                "Vinberg's reflection chamber for the units of "
                "-x0^2+x1^2+...+x14^2; reproduced by Guglielmetti"
            ),
            "primary_locator": (
                "R. Guglielmetti, CoxIter, LMS J. Comput. Math. 18 (2015), "
                "Appendix A, Figure A.1, doi:10.1112/S1461157015000273"
            ),
            "published_edges_coxiter_encoding": [list(edge) for edge in VINBERG_14_EDGES],
            "our_label_to_published_label": {
                str(label): image for label, image in OUR_TO_VINBERG.items()
            },
            "exact_edge_labeled_graph_match": diagram_matches,
            "published_f_vector": list(VINBERG_14_F_VECTOR),
            "published_total_vertex_count": VINBERG_14_F_VECTOR[0],
            "published_edge_count": VINBERG_14_F_VECTOR[1],
            "published_ideal_vertex_count": 5,
            "published_finite_vertex_count": VINBERG_14_F_VECTOR[0] - 5,
            "published_finite_volume_noncompact": True,
            "published_arithmetic": True,
        },
        "algebraic_realization": {
            "kondo_identification": (
                "This is the unique graph reflection fixing 15 curve walls, hence "
                "Kondo's middle-horizontal reflection tau_2 in Section 4.4."
            ),
            "symplectic_involution": False,
            "symplectic_obstruction": (
                "Kondo observes that tau_2 preserves 15 rational curves; a "
                "symplectic K3 involution has only eight fixed points."
            ),
            "nonsymplectic_involution": False,
            "nonsymplectic_torelli_obstruction": (
                "For an involution on a K3 with exactly this NS, a nonsymplectic "
                "action must be -identity on T. It is therefore trivial on the "
                "2-elementary discriminant, whereas tau_2 swaps two nonzero "
                "discriminant classes. The actions cannot glue in the K3 lattice."
            ),
            "realized_K3_involution": False,
            "special_period_scope": (
                "This excludes an order-two surface automorphism for every K3 with "
                "exactly this Neron-Severi lattice. It does not classify possible "
                "higher-order special-period automorphisms inducing the same order-two "
                "action on NS; such an action would not realize the proposed involution quotient."
            ),
        },
        "conclusion": {
            "new_abstract_hyperbolic_coxeter_polytope": False,
            "new_algebro_geometric_quotient_example": False,
            "status": (
                "Exact rediscovery of Vinberg's published H14 chamber through a "
                "non-realizable symmetry of the K3 curve graph."
            ),
        },
        "all_checks_pass": (
            proposed in automorphisms
            and fixed_counts[proposed] == 15
            and discriminant_action == (0, 2, 1, 3)
            and orbit_gram.rank() == 15
            and exact_inertia(orbit_gram) == (1, 14, 2)
            and diagram_matches
            and VINBERG_14_F_VECTOR[0] - 5 == 89
            and VINBERG_14_F_VECTOR[1] == 704
        ),
    }
    payload["certificate_payload_sha256"] = _canonical_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, Any]) -> dict[str, bool]:
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256", None)
    fresh = build_certificate()
    checks = {
        "payload_digest_matches": digest == _canonical_sha256(body),
        "fresh_replay_matches": payload == fresh,
        "all_exact_checks_pass": payload.get("all_checks_pass") is True,
        "not_realized": payload.get("algebraic_realization", {}).get(
            "realized_K3_involution"
        )
        is False,
        "not_new": payload.get("conclusion", {}).get(
            "new_abstract_hyperbolic_coxeter_polytope"
        )
        is False,
    }
    if not all(checks.values()):
        raise AssertionError(f"rho16 H14 fold identification replay failed: {checks}")
    return checks


if __name__ == "__main__":
    print(json.dumps(build_certificate(), indent=2, sort_keys=True))
