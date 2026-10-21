#!/usr/bin/env python3
"""Exhaust every weighted-wall-graph involution in the rank-18/19 K3 data.

The search starts from the complete published (-2)-curve Gram matrices in
``k3_high_rank_nef_chambers.json``.  It enumerates the full automorphism group
of each *weighted complete Gram graph* by exact backtracking, selects every
nonidentity involution, forms its orbit-sum Gram matrix, and applies the full
fold and Vinberg--CoxIter audits from ``verify_k3_quotient_nef_slices``.

A weighted graph symmetry is not silently promoted to a surface
automorphism.  The final realization labels use the separate published facts
already recorded in ``k3_quotient_nef_slices.json``: the rank-18 symmetry is
induced by Kondō's involution, while the full rank-19 S3 graph action is
induced by surface automorphisms.  Novelty is then tested modulo that realized
action, not modulo arbitrary relabeling alone.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia
from scripts.verify_k3_quotient_nef_slices import analyze_case


ROOT = Path(__file__).resolve().parents[1]
HIGH_RANK_PATH = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
PUBLISHED_FOLDS_PATH = ROOT / "data/known_examples/k3_quotient_nef_slices.json"
CERTIFICATE_PATH = (
    ROOT / "data/candidates/k3_high_rank_all_involution_folds.json"
)

Permutation = tuple[int, ...]


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _canonical_sha256(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def exact_weighted_gram_automorphisms(gram: sp.Matrix) -> tuple[Permutation, ...]:
    """Enumerate every permutation ``p`` with ``G[p(i),p(j)]=G[i,j]``.

    The per-vertex multiset of Gram entries is only a pruning invariant.  The
    backtracking checks every relation to an already assigned vertex; hence a
    completed bijection preserves the entire matrix and the enumeration is
    exhaustive.
    """

    if gram.rows != gram.cols or gram != gram.T:
        raise ValueError("the wall Gram matrix must be square and symmetric")
    size = gram.rows
    relations = tuple(
        tuple(int(gram[left, right]) for right in range(size))
        for left in range(size)
    )
    signatures = tuple(
        tuple(sorted(relations[index][other] for other in range(size)))
        for index in range(size)
    )
    candidates = tuple(
        tuple(target for target in range(size) if signatures[target] == signatures[source])
        for source in range(size)
    )
    order = tuple(sorted(range(size), key=lambda index: (len(candidates[index]), index)))
    image = [-1] * size
    used: set[int] = set()
    answers: list[Permutation] = []

    def visit(depth: int) -> None:
        if depth == size:
            answer = tuple(image)
            if any(
                relations[left][right]
                != relations[answer[left]][answer[right]]
                for left in range(size)
                for right in range(size)
            ):
                raise AssertionError("backtracking admitted a false automorphism")
            answers.append(answer)
            return
        source = order[depth]
        for target in candidates[source]:
            if target in used:
                continue
            if any(
                image[other] >= 0
                and relations[source][other]
                != relations[target][image[other]]
                for other in range(size)
            ):
                continue
            image[source] = target
            used.add(target)
            visit(depth + 1)
            used.remove(target)
            image[source] = -1

    visit(0)
    return tuple(sorted(answers))


def _identity(size: int) -> Permutation:
    return tuple(range(size))


def _compose(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left o right``."""

    return tuple(left[right[index]] for index in range(len(left)))


def _inverse(permutation: Permutation) -> Permutation:
    answer = [-1] * len(permutation)
    for source, target in enumerate(permutation):
        answer[target] = source
    return tuple(answer)


def _conjugate(move: Permutation, action: Permutation) -> Permutation:
    return _compose(_compose(move, action), _inverse(move))


def _order(permutation: Permutation) -> int:
    current = _identity(len(permutation))
    for exponent in range(1, len(permutation) + 2):
        current = _compose(permutation, current)
        if current == _identity(len(permutation)):
            return exponent
    raise AssertionError("permutation order exceeded the symmetric-group bound")


def _permutation_from_cycles(size: int, cycles: Iterable[Iterable[int]]) -> Permutation:
    image = list(range(size))
    seen: set[int] = set()
    for raw_cycle in cycles:
        cycle = [int(label) - 1 for label in raw_cycle]
        if not cycle or any(index < 0 or index >= size for index in cycle):
            raise ValueError(f"invalid cycle {raw_cycle}")
        if seen.intersection(cycle):
            raise ValueError(f"overlapping cycles at {raw_cycle}")
        seen.update(cycle)
        for position, source in enumerate(cycle):
            image[source] = cycle[(position + 1) % len(cycle)]
    if seen != set(range(size)):
        raise ValueError("the cycle list does not cover every wall")
    return tuple(image)


def _cycles(permutation: Permutation) -> list[list[int]]:
    seen: set[int] = set()
    answer: list[list[int]] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle: list[int] = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(current + 1)
            current = permutation[current]
        answer.append(cycle)
    return answer


def _orbit_sum_gram(gram: sp.Matrix, cycles: Sequence[Sequence[int]]) -> sp.Matrix:
    incidence = sp.zeros(gram.rows, len(cycles))
    for orbit_index, cycle in enumerate(cycles):
        for label in cycle:
            incidence[int(label) - 1, orbit_index] = 1
    return incidence.T * gram * incidence


def _orbit_transport(
    published_cycles: Sequence[Sequence[int]],
    current_cycles: Sequence[Sequence[int]],
    conjugator: Permutation,
) -> list[int]:
    current_lookup = {
        frozenset(int(label) - 1 for label in cycle): index
        for index, cycle in enumerate(current_cycles)
    }
    return [
        current_lookup[
            frozenset(conjugator[int(label) - 1] for label in published_orbit)
        ]
        for published_orbit in published_cycles
    ]


def _source_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (HIGH_RANK_PATH, PUBLISHED_FOLDS_PATH)
    }


def _summarize_fold(
    *,
    rho: int,
    index: int,
    permutation: Permutation,
    automorphisms: tuple[Permutation, ...],
    published_permutation: Permutation,
    published_cycles: list[list[int]],
    gram: sp.Matrix,
    high_rank_case: dict[str, Any],
    published_case: dict[str, Any],
) -> dict[str, Any]:
    cycles = _cycles(permutation)
    orbit_gram = _orbit_sum_gram(gram, cycles)
    conjugators = tuple(
        move
        for move in automorphisms
        if _conjugate(move, published_permutation) == permutation
    )
    if not conjugators:
        raise AssertionError("an involution left the published conjugacy class")
    conjugator = conjugators[0]
    transport = _orbit_transport(published_cycles, cycles, conjugator)
    published_orbit_gram = _orbit_sum_gram(gram, published_cycles)
    transported_gram_matches = all(
        published_orbit_gram[left, right]
        == orbit_gram[transport[left], transport[right]]
        for left in range(len(transport))
        for right in range(len(transport))
    )
    if not transported_gram_matches:
        raise AssertionError("conjugate actions produced nonisometric orbit Grams")

    # This is deliberately a direct audit for every action, even though
    # conjugacy would suffice to transport the representative certificate.
    direct_case = dict(published_case)
    direct_case["involution_cycles"] = cycles
    audit = analyze_case(direct_case, high_rank_case)

    two_orbit_rows: list[dict[str, Any]] = []
    for orbit_index, cycle in enumerate(cycles):
        if len(cycle) != 2:
            continue
        left, right = cycle[0] - 1, cycle[1] - 1
        intersection = int(gram[left, right])
        normal_square = int(orbit_gram[orbit_index, orbit_index])
        two_orbit_rows.append(
            {
                "orbit": cycle,
                "curve_intersection": intersection,
                "orbit_normal_square": normal_square,
                "role": (
                    "disjoint_reflective_facet"
                    if normal_square < 0 and intersection == 0
                    else "redundant_timelike_inequality"
                    if normal_square > 0
                    else "other"
                ),
            }
        )

    exact_published_action = permutation == published_permutation
    dimension = int(audit["hyperbolic_dimension"])
    gate_checks = {
        "weighted_gram_automorphism": audit["checks"]["published_wall_gram_preserved"],
        "nonidentity_involution": permutation != _identity(len(permutation))
        and audit["checks"]["action_is_an_involution"],
        "fixed_space_is_lorentzian": audit["checks"]["lorentzian_fixed_space"],
        "hyperbolic_dimension_gt_9": dimension > 9,
        "all_facet_pairs_are_exact_coxeter_relations": (
            audit["checks"]["no_ultraparallel_edges"]
        ),
        "all_folded_facet_reflections_integral": audit["checks"][
            "folded_reflections_are_integral"
        ],
        "all_retained_walls_irredundant": audit["checks"]["every_wall_is_a_facet"],
        "vinberg_coxiter_finite_volume": audit["checks"][
            "vinberg_finite_volume_criterion"
        ]
        and audit["checks"]["literal_coxiter_prop_1_10_finite_covolume"],
        "direct_full_audit_passes": audit["all_checks_pass"],
        "exactly_conjugate_to_published_action": bool(conjugators),
        "orbit_gram_transport_matches": transported_gram_matches,
    }
    geometric_checks = [
        value for key, value in gate_checks.items() if key != "hyperbolic_dimension_gt_9"
    ]
    target_passes = all(gate_checks.values())
    if rho == 18:
        realization = {
            "status": "published_surface_automorphism",
            "reason": (
                "This is exactly Kondō's realized rank-18 graph involution; "
                "the graph census alone is not being used as realization evidence."
            ),
        }
    else:
        realization = {
            "status": "published_surface_automorphism",
            "reason": (
                "Kondō's full S3 graph action is realized on the surface; this "
                "transposition is either the archived one or its conjugate by "
                "that realized S3 action."
            ),
        }

    return {
        "fold_id": f"rho{rho}_involution_{index}",
        "ambient_rho": rho,
        "permutation_images_zero_based": list(permutation),
        "cycles": cycles,
        "equals_archived_published_action": exact_published_action,
        "published_action_conjugator_count": len(conjugators),
        "chosen_conjugator_images_zero_based": list(conjugator),
        "published_orbit_to_this_orbit_zero_based": transport,
        "fixed_rank": int(audit["fixed_rank"]),
        "hyperbolic_dimension": dimension,
        "orbit_gram": _matrix_payload(orbit_gram),
        "orbit_gram_sha256": _canonical_sha256(_matrix_payload(orbit_gram)),
        "orbit_gram_rank": int(orbit_gram.rank()),
        "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
        "facet_orbits": audit["facet_orbits"],
        "redundant_timelike_orbits": audit["redundant_timelike_orbits"],
        "two_curve_orbit_audit": two_orbit_rows,
        "facet_gram_rank": int(audit["facet_gram_rank"]),
        "facet_gram_inertia": audit["facet_gram_inertia"],
        "coxeter_edges": audit["coxeter_edges"],
        "coxeter_edge_counts": audit["coxeter_edge_counts"],
        "fold_reflection_checks": audit["fold_reflection_checks"],
        "facet_support_checks": audit["facet_support_checks"],
        "finite_vertex_count": int(audit["finite_vertex_count"]),
        "ideal_vertex_count": int(audit["ideal_vertex_count"]),
        "vinberg_volume_criterion": audit["vinberg_volume_criterion"],
        "literal_coxiter_prop_1_10": audit["literal_coxiter_prop_1_10"],
        "finite_volume": audit["finite_volume"],
        "compact": audit["compact"],
        "gate_checks": gate_checks,
        "all_geometric_gates_except_target_dimension_pass": all(geometric_checks),
        "passes_dimension_gt_9_target": target_passes,
        "algebraic_realization": realization,
        "novelty": (
            "archived_published_fold"
            if exact_published_action
            else "realized_S3_conjugate_of_archived_published_fold"
        ),
        "certified_new_fold": False,
    }


@lru_cache(maxsize=1)
def build_certificate() -> dict[str, Any]:
    high_rank = _load(HIGH_RANK_PATH)
    published = _load(PUBLISHED_FOLDS_PATH)
    high_rank_by_rho = {
        int(case["rho"]): case
        for case in high_rank["cases"]
        if int(case["rho"]) in (18, 19)
    }
    published_by_rho = {
        int(case["ambient_rho"]): case for case in published["cases"]
    }

    case_payloads: list[dict[str, Any]] = []
    all_folds: list[dict[str, Any]] = []
    for rho in (18, 19):
        high_rank_case = high_rank_by_rho[rho]
        published_case = published_by_rho[rho]
        gram = build_gram(high_rank_case)
        automorphisms = exact_weighted_gram_automorphisms(gram)
        identity = _identity(gram.rows)
        automorphism_set = set(automorphisms)
        group_checks = {
            "identity_present": identity in automorphism_set,
            "all_actions_preserve_gram": all(
                all(
                    gram[left, right] == gram[action[left], action[right]]
                    for left in range(gram.rows)
                    for right in range(gram.rows)
                )
                for action in automorphisms
            ),
            "closed_under_composition": all(
                _compose(left, right) in automorphism_set
                for left in automorphisms
                for right in automorphisms
            ),
            "closed_under_inverse": all(
                _inverse(action) in automorphism_set for action in automorphisms
            ),
        }
        involutions = tuple(
            action
            for action in automorphisms
            if action != identity and _compose(action, action) == identity
        )
        published_cycles = [
            [int(label) for label in cycle]
            for cycle in published_case["involution_cycles"]
        ]
        published_permutation = _permutation_from_cycles(
            gram.rows, published_cycles
        )
        if published_permutation not in involutions:
            raise AssertionError("the archived action is absent from the census")

        involution_classes = {
            tuple(sorted(_conjugate(move, action) for move in automorphisms))
            for action in involutions
        }
        fold_payloads = [
            _summarize_fold(
                rho=rho,
                index=index,
                permutation=action,
                automorphisms=automorphisms,
                published_permutation=published_permutation,
                published_cycles=published_cycles,
                gram=gram,
                high_rank_case=high_rank_case,
                published_case=published_case,
            )
            for index, action in enumerate(involutions, start=1)
        ]
        all_folds.extend(fold_payloads)
        order_histogram = Counter(_order(action) for action in automorphisms)
        expected_group = "C2" if rho == 18 else "S3"
        expected_order_histogram = {1: 1, 2: 1} if rho == 18 else {1: 1, 2: 3, 3: 2}
        case_checks = {
            **group_checks,
            "expected_group_order": len(automorphisms) == (2 if rho == 18 else 6),
            "expected_element_order_histogram": dict(order_histogram)
            == expected_order_histogram,
            "expected_nonidentity_involution_count": len(involutions)
            == (1 if rho == 18 else 3),
            "single_involution_conjugacy_class": len(involution_classes) == 1,
            "published_action_in_census": published_permutation in involutions,
            "all_direct_fold_audits_pass": all(
                fold["all_geometric_gates_except_target_dimension_pass"]
                for fold in fold_payloads
            ),
            "no_new_fold_modulo_published_realized_action": not any(
                fold["certified_new_fold"] for fold in fold_payloads
            ),
        }
        case_payloads.append(
            {
                "ambient_rho": rho,
                "wall_count": gram.rows,
                "weighted_graph_automorphism_group": expected_group,
                "automorphism_count": len(automorphisms),
                "element_order_histogram": {
                    str(order): count for order, count in sorted(order_histogram.items())
                },
                "automorphism_images_zero_based": [list(action) for action in automorphisms],
                "nonidentity_involution_count": len(involutions),
                "nonidentity_involution_conjugacy_class_count": len(involution_classes),
                "archived_published_cycles": published_cycles,
                "folds": fold_payloads,
                "checks": case_checks,
                "all_checks_pass": all(case_checks.values()),
            }
        )

    target_folds = [
        fold["fold_id"] for fold in all_folds if fold["passes_dimension_gt_9_target"]
    ]
    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "artifact_id": "k3_high_rank_all_involution_folds",
        "certificate_kind": "exhaustive_exact_weighted_graph_involution_fold_census",
        "source_hashes": _source_hashes(),
        "search_scope": {
            "ambient_picard_ranks": [18, 19],
            "actions": "all nonidentity involutions of the complete weighted wall Gram graph",
            "target_hyperbolic_dimension": ">9",
            "volume_test": "Vinberg/CoxIter Proposition 1.10 literal endpoint audit on complete folded walls",
        },
        "cases": case_payloads,
        "totals": {
            "weighted_graph_automorphisms": sum(
                case["automorphism_count"] for case in case_payloads
            ),
            "nonidentity_involutions": len(all_folds),
            "direct_finite_volume_fold_certificates": sum(
                bool(fold["finite_volume"]) for fold in all_folds
            ),
            "dimension_gt_9_target_actions": target_folds,
            "dimension_gt_9_geometric_conjugacy_classes": 1,
            "certified_new_folds": 0,
        },
        "realization_boundary": (
            "Graph automorphism was only the combinatorial gate. Realization is "
            "asserted here solely because Kondō separately realizes the rank-18 "
            "involution and the full rank-19 S3 action by K3 automorphisms."
        ),
        "novelty_conclusion": (
            "The exhaustive rank-18/19 weighted-wall-graph involution census "
            "contains exactly the two published fold classes: the H9 rank-18 "
            "fold and the H10 rank-19 transposition fold. The other two H10 "
            "actions are realized S3-conjugates, so no additional Coxeter "
            "polytope is obtained."
        ),
        "new_hyperbolic_example": False,
        "all_checks_pass": all(case["all_checks_pass"] for case in case_payloads),
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
        "four_nonidentity_involutions": payload.get("totals", {}).get(
            "nonidentity_involutions"
        )
        == 4,
        "no_new_fold": payload.get("totals", {}).get("certified_new_folds") == 0,
    }
    if not all(checks.values()):
        raise AssertionError(f"K3 involution-fold replay failed: {checks}")
    return checks


if __name__ == "__main__":
    print(json.dumps(build_certificate(), indent=2, sort_keys=True))
