#!/usr/bin/env python3
"""Complete exactly-two-cross-edge scout for Record 30 ``a=2`` orbit 10.

This is deliberately a *minimal shell* computation.  It first reconstructs
every sparse Lanner placement on the four higher minimal nonfaces, with all
other meeting pairs right-angled.  It then adds exactly two finite non-right
edges, labelled 3, 4, or 5, which connect the three components of the
minimal-nonface hypergraph.  Every actual rank-nine vertex is retested.

For each local survivor, the script searches the eight largest principal
blocks which contain no dotted entry.  The block is scaled by two, making its
entries integral in Z[sqrt(2),sqrt(5)], and its determinant is reduced through
an explicitly checked finite-field specialization.  A nonzero residue is an
exact rank-at-least-eleven certificate, not a floating-point heuristic.

Nothing here treats assignments with three or more cross-component edges.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import argparse
import json
from pathlib import Path
from typing import Mapping

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)


N = 14
DIMENSION = 9
DOTTED_PAIRS = ((0, 11), (1, 12), (2, 13))
HIGHER_NONFACES = (
    (0, 3, 4, 5),
    (6, 7, 8, 11),
    (3, 4, 5, 9, 10),
    (6, 7, 8, 9, 10),
)
MINIMAL_NONFACES = DOTTED_PAIRS + HIGHER_NONFACES

MAIN_COMPONENT = (0, 3, 4, 5, 6, 7, 8, 9, 10, 11)
PAIR_COMPONENT_ONE = (1, 12)
PAIR_COMPONENT_TWO = (2, 13)
COMPONENTS = (MAIN_COMPONENT, PAIR_COMPONENT_ONE, PAIR_COMPONENT_TWO)

FINITE_LABELS = (3, 4, 5)
LabelState = tuple[tuple[Pair, int], ...]
CrossChoice = tuple[Pair, int]

# 1009 is prime; 439^2 = 2 and 244^2 = 5 modulo 1009.  The scaled Gram
# entries therefore have a ring specialization to F_1009.
MODULUS = 1009
SQRT2_MOD = 439
SQRT5_MOD = 244

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a2_orbit10_two_cross.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "df4216dc9d8141d2a3d9466fae0fe87f580a80058151e01dc5cb8386f679a94a"
)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return _canonical_sha256(payload)


def _actual_vertices() -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(face) for face in MINIMAL_NONFACES)
    return tuple(
        chosen
        for chosen in combinations(range(N), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=1_048_576)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def _restriction_finite(
    restriction: tuple[int, ...], labels: dict[Pair, int]
) -> bool:
    return _finite_code(
        len(restriction),
        tuple(
            labels[pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


def _next_support(
    supports: tuple[tuple[int, ...], ...], remaining: set[int], union: set[int]
) -> int:
    if not union:
        return max(remaining, key=lambda index: (len(supports[index]), -index))
    return max(
        remaining,
        key=lambda index: (
            len(set(supports[index]) & union),
            len(supports[index]),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def sparse_states(
    actual_vertices: tuple[tuple[int, ...], ...]
) -> tuple[tuple[LabelState, ...], tuple[dict[str, object], ...]]:
    """Reconstruct the complete sparse terminal list and join transcript."""

    supports = HIGHER_NONFACES
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages: list[dict[str, object]] = []
    while remaining:
        support_index = _next_support(supports, remaining, union)
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(facet for facet in vertex if facet in union)
                    for vertex in actual_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        options = _placement_options(support)
        passing: set[LabelState] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    _restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "option_count": len(options),
                "raw_extension_count": len(options)
                * (stages[-1]["passing_state_count"] if stages else 1),
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "restriction_count": len(restrictions),
            }
        )
        remaining.remove(support_index)
    return states, tuple(stages)


def _component_cross_edges() -> tuple[tuple[Pair, ...], ...]:
    return tuple(
        tuple(pair(left, right) for left in COMPONENTS[i] for right in COMPONENTS[j])
        for i, j in ((0, 1), (0, 2), (1, 2))
    )


def _two_edge_placements() -> tuple[tuple[Pair, Pair], ...]:
    """Every literal two-edge tree on the three component vertices."""

    categories = _component_cross_edges()
    return tuple(
        (first, second)
        for left_category, right_category in combinations(range(3), 2)
        for first in categories[left_category]
        for second in categories[right_category]
    )


def _base_labels(state: LabelState) -> dict[Pair, int]:
    dotted = {pair(*edge) for edge in DOTTED_PAIRS}
    labels = {
        pair(left, right): 2
        for left, right in combinations(range(N), 2)
        if pair(left, right) not in dotted
    }
    labels.update(dict(state))
    return labels


def _scaled_gram_value(label: int) -> int:
    if label == 2:
        return 0
    if label == 3:
        return -2
    if label == 4:
        return -2 * SQRT2_MOD
    if label == 5:
        return -(1 + SQRT5_MOD)
    raise ValueError(f"unsupported finite label {label}")


def _det_mod(matrix: list[list[int]]) -> int:
    p = MODULUS
    work = [[value % p for value in row] for row in matrix]
    determinant = 1
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        value = work[column][column]
        determinant = determinant * value % p
        inverse = pow(value, p - 2, p)
        for row in range(column + 1, len(work)):
            multiplier = work[row][column] * inverse % p
            if not multiplier:
                continue
            for entry in range(column, len(work)):
                work[row][entry] = (
                    work[row][entry] - multiplier * work[column][entry]
                ) % p
    return determinant % p


def _fixed_minor_witness(labels: dict[Pair, int]) -> dict[str, object] | None:
    for excluded in product(*DOTTED_PAIRS):
        order = tuple(vertex for vertex in range(N) if vertex not in excluded)
        if len(order) != 11:
            raise AssertionError("a dotted-free principal order changed")
        matrix = [
            [
                4
                if left == right
                else _scaled_gram_value(labels[pair(left, right)])
                for right in order
            ]
            for left in order
        ]
        residue = _det_mod(matrix)
        if residue:
            return {
                "excluded_facets": list(excluded),
                "principal_order": list(order),
                "scaled_determinant_mod_1009": residue,
            }
    return None


def build() -> dict[str, object]:
    if pow(SQRT2_MOD, 2, MODULUS) != 2 or pow(SQRT5_MOD, 2, MODULUS) != 5:
        raise AssertionError("finite-field radical specialization changed")
    # Trial division is ample for this small fixed modulus and keeps primality
    # part of the executable certificate.
    if any(MODULUS % divisor == 0 for divisor in range(2, int(MODULUS**0.5) + 1)):
        raise AssertionError("the Gram witness modulus is not prime")

    vertices = _actual_vertices()
    states, stages = sparse_states(vertices)
    placements = _two_edge_placements()
    if len(vertices) != 128 or len(states) != 648 or len(placements) != 560:
        raise AssertionError("the orbit-10 shell cardinalities changed")

    vertex_sets = tuple(frozenset(vertex) for vertex in vertices)
    affected = {
        edge: tuple(
            index
            for index, vertex in enumerate(vertex_sets)
            if set(edge) <= vertex
        )
        for category in _component_cross_edges()
        for edge in category
    }

    individual_trial_count = 0
    individually_passing_count = 0
    pair_residual_trial_count = 0
    local_survivors: list[tuple[int, Pair, int, Pair, int]] = []
    witness_histogram: Counter[tuple[tuple[int, ...], int]] = Counter()
    no_fixed_minor_witness: list[tuple[int, Pair, int, Pair, int]] = []
    per_state_local_histogram: Counter[int] = Counter()
    per_label_pair_histogram: Counter[tuple[int, int]] = Counter()

    for state_index, state in enumerate(states):
        labels = _base_labels(state)
        individually_valid: dict[CrossChoice, bool] = {}
        for category in _component_cross_edges():
            for edge in category:
                for label in FINITE_LABELS:
                    individual_trial_count += 1
                    labels[edge] = label
                    valid = all(
                        _restriction_finite(vertices[index], labels)
                        for index in affected[edge]
                    )
                    labels[edge] = 2
                    individually_valid[(edge, label)] = valid
                    individually_passing_count += valid

        state_survivors = 0
        for first_edge, second_edge in placements:
            common_vertices = tuple(
                sorted(set(affected[first_edge]) & set(affected[second_edge]))
            )
            for first_label in FINITE_LABELS:
                if not individually_valid[(first_edge, first_label)]:
                    continue
                for second_label in FINITE_LABELS:
                    if not individually_valid[(second_edge, second_label)]:
                        continue
                    pair_residual_trial_count += 1
                    labels[first_edge] = first_label
                    labels[second_edge] = second_label
                    valid = all(
                        _restriction_finite(vertices[index], labels)
                        for index in common_vertices
                    )
                    if valid:
                        survivor = (
                            state_index,
                            first_edge,
                            first_label,
                            second_edge,
                            second_label,
                        )
                        local_survivors.append(survivor)
                        state_survivors += 1
                        per_label_pair_histogram[
                            tuple(sorted((first_label, second_label)))
                        ] += 1
                        witness = _fixed_minor_witness(labels)
                        if witness is None:
                            no_fixed_minor_witness.append(survivor)
                        else:
                            witness_histogram[
                                (
                                    tuple(witness["excluded_facets"]),
                                    int(witness["scaled_determinant_mod_1009"]),
                                )
                            ] += 1
                    labels[first_edge] = 2
                    labels[second_edge] = 2
        per_state_local_histogram[state_survivors] += 1

    # The digest is set-valued evidence, independent of either enumerator's
    # traversal order.
    survivor_payload = sorted(
        [index, list(first), first_label, list(second), second_label]
        for index, first, first_label, second, second_label in local_survivors
    )
    unwitnessed_payload = sorted(
        [index, list(first), first_label, list(second), second_label]
        for index, first, first_label, second, second_label in no_fixed_minor_witness
    )
    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record30_a2_orbit10_exact_two_cross_no_go",
        "claim_status": "frozen_exact_theorem_with_independent_replay",
        "claim": (
            "Every compact-Coxeter labeling in the exactly-two-cross-component "
            "finite-nonright-edge shell of the Record 30 a=2 orbit-10 nerve "
            "has Gram rank at least 11 and therefore cannot occur in H^9."
        ),
        "scope": {
            "record_index": 30,
            "addition_schedule": 2,
            "positive_A4_orbit_index": 10,
            "exactly_two_cross_component_nonright_meeting_edges": True,
            "cross_edge_labels": list(FINITE_LABELS),
            "three_or_more_cross_edges_excluded": False,
            "whole_orbit10_nerve_excluded": False,
            "whole_record30_a2_schedule_excluded": False,
            "finite_volume_noncompact_statement": False,
        },
        "nerve": {
            "minimal_nonfaces": [list(face) for face in MINIMAL_NONFACES],
            "components": [list(component) for component in COMPONENTS],
            "actual_rank_nine_vertex_count": len(vertices),
            "actual_vertices_sha256": _canonical_sha256(vertices),
        },
        "sparse_gate": {
            "stages": list(stages),
            "terminal_state_count": len(states),
            "terminal_states_sha256": _canonical_sha256(states),
        },
        "two_cross_shell": {
            "component_pair_edge_counts": [
                len(category) for category in _component_cross_edges()
            ],
            "literal_edge_pair_placement_count": len(placements),
            "literal_labelled_trial_count": len(states)
            * len(placements)
            * len(FINITE_LABELS) ** 2,
            "single_edge_monotonicity_trial_count": individual_trial_count,
            "single_edge_passing_count": individually_passing_count,
            "two_edge_residual_trial_count": pair_residual_trial_count,
            "local_survivor_count": len(local_survivors),
            "local_survivor_sha256": _canonical_sha256(survivor_payload),
            "local_survivor_histogram_per_sparse_state": {
                str(count): multiplicity
                for count, multiplicity in sorted(per_state_local_histogram.items())
            },
            "local_survivor_label_pair_histogram": {
                f"{labels[0]},{labels[1]}": multiplicity
                for labels, multiplicity in sorted(per_label_pair_histogram.items())
            },
        },
        "fixed_minor_gate": {
            "normalization": (
                "two times the twice-normalized Gram: diagonal 4, m3=-2, "
                "m4=-2*sqrt(2), m5=-(1+sqrt(5))"
            ),
            "coefficient_ring": "Z[sqrt(2),sqrt(5)]",
            "modulus": MODULUS,
            "sqrt2_residue": SQRT2_MOD,
            "sqrt5_residue": SQRT5_MOD,
            "candidate_dotted_free_principal_order_count": 8,
            "witnessed_rank_at_least_eleven_count": (
                len(local_survivors) - len(no_fixed_minor_witness)
            ),
            "witness_histogram": [
                {
                    "excluded_facets": list(excluded),
                    "scaled_determinant_mod_1009": residue,
                    "multiplicity": multiplicity,
                }
                for (excluded, residue), multiplicity in sorted(
                    witness_histogram.items()
                )
            ],
            "no_fixed_minor_witness_count": len(no_fixed_minor_witness),
            "no_fixed_minor_witness_sha256": _canonical_sha256(unwitnessed_payload),
            "no_fixed_minor_witnesses": unwitnessed_payload,
        },
        "conclusion": {
            "all_local_survivors_excluded_by_exact_rank_gate": not no_fixed_minor_witness,
            "new_hyperbolic_polytope": False,
        },
        "independent_replay": {
            "script": "scripts/independent_h9_14_record30_a2_orbit10_two_cross.py",
            "imports_production_enumerator": False,
            "imports_coxeter_label_module": False,
            "opposite_support_join_order": [3, 1, 2, 0],
            "second_modulus": 1129,
            "second_sqrt2_residue": 406,
            "second_sqrt5_residue": 474,
        },
    }
    document["certificate_sha256"] = _certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record30_a2_orbit10_exact_two_cross_no_go":
        raise ValueError("unexpected Record 30 a2 orbit-10 two-cross kind")
    rebuilt = build()
    if document != rebuilt:
        raise ValueError("the Record 30 orbit-10 two-cross certificate did not replay")
    digest = _certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("the Record 30 orbit-10 two-cross digest is invalid")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("the Record 30 orbit-10 two-cross certificate is not frozen")
    return {
        "terminal_sparse_state_count": rebuilt["sparse_gate"]["terminal_state_count"],
        "literal_labelled_trial_count": rebuilt["two_cross_shell"]["literal_labelled_trial_count"],
        "local_survivor_count": rebuilt["two_cross_shell"]["local_survivor_count"],
        "rank_at_least_eleven_count": rebuilt["fixed_minor_gate"]["witnessed_rank_at_least_eleven_count"],
        "unresolved_local_survivor_count": rebuilt["fixed_minor_gate"]["no_fixed_minor_witness_count"],
        "certificate_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(verify(json.loads(ARCHIVE.read_text(encoding="utf-8"))), indent=2, sort_keys=True))
    else:
        print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
