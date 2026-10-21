#!/usr/bin/env python3
"""Exact, deliberately non-classificatory certificate for known ``P_{6,10}``.

This reconstructs the sparse unit Gram matrix over ``Q(sqrt(2),sqrt(5))``
and its mixed golden crystallographic scaling.  Compactness is a pinned
source/CoxIter provenance claim in the artifact, not a conclusion drawn here
from the recorded f-vector.  In particular, this script does not certify the
disputed uniqueness assertion of arXiv:2608.02894v1.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "known_examples" / "p6_10_exact_certificate.json"

ORDINARY_EDGES = {
    (1, 2, 5), (2, 3, 3), (3, 4, 3), (4, 5, 3), (5, 6, 3),
    (9, 10, 5), (2, 7, 4), (2, 8, 4), (5, 9, 3), (6, 10, 5),
}
DOTTED_PAIRS = {(6, 7), (7, 8), (8, 9)}
MAPPING = {
    "1": "s8", "2": "s4", "3": "s6", "4": "s5", "5": "s3",
    "6": "s1", "7": "t0_s2", "8": "s2", "9": "t0_s1", "10": "s7",
}
HASHES = {
    "graphs/testing_index2.tar.gz": "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036",
    "testing_index2/6-per09-tumarkin-03_1.coxiter": "2d32aa1bceb38410bf2029b9acc54b5b85671cd0c3d4b72ed8cd6841fb35923c",
    "testing_index2/6-per09-tumarkin-03_0.coxiter": "f41d2cbcc728f4e063f098cd4c1ced0221204d155d0f1987edcfff2699148750",
    "graphs/cocompact/6-Bugaenko.coxiter": "0ed96024d9da1861c4dc689ba65df1449c8b4e906e444af922ffef6fb50e3b47",
    "graphs/cocompact/6-per09-tumarkin-03.coxiter": "d7a0c92eaa30333c0482eda887ba4944983cbb71ba2765f8689e6a4501d96337",
    "testing/tests.txt": "9151721b589aba2a531b322706b2e42aede5f3793dbd8980b8bc6e6bca830529",
    "CoxIter/index2.cpp": "84317a17bd52c71bda9b9ab0441383a2d8e4745dc7370bebf6c2e579196b995b",
}


def load_artifact() -> dict[str, Any]:
    return json.loads(DATA_PATH.read_text(encoding="utf-8"))


def _qfield(raw: list[int]) -> sp.Expr:
    """Decode coordinates on 1,sqrt(2),sqrt(5),sqrt(10).

    A fifth coordinate, when present, is the common denominator.
    """

    if len(raw) not in (4, 5):
        raise ValueError(f"bad multiquadratic coordinate {raw}")
    denominator = raw[4] if len(raw) == 5 else 1
    if not isinstance(denominator, int) or denominator == 0:
        raise ValueError(f"bad denominator in {raw}")
    q2, q5 = sp.sqrt(2), sp.sqrt(5)
    return sp.expand(
        (raw[0] + raw[1] * q2 + raw[2] * q5 + raw[3] * q2 * q5)
        / denominator
    )


def build_unit_gram(payload: dict[str, Any]) -> sp.Matrix:
    datum = payload["unit_gram"]
    if set(map(tuple, datum["ordinary_edges_i_j_m"])) != ORDINARY_EDGES:
        raise ValueError("ordinary edge stock differs from the pinned P_{6,10} diagram")

    gram = sp.eye(10) * _qfield(datum["diagonal"])
    seen: set[tuple[int, int]] = set()
    entries = datum["ordinary_gram_entry_by_m"]
    exact_cosines = {
        3: -sp.Rational(1, 2),
        4: -sp.sqrt(2) / 2,
        5: -(1 + sp.sqrt(5)) / 4,
    }
    for raw_i, raw_j, raw_m in datum["ordinary_edges_i_j_m"]:
        i, j, m = int(raw_i), int(raw_j), int(raw_m)
        pair = tuple(sorted((i, j)))
        if not (1 <= pair[0] < pair[1] <= 10) or pair in seen:
            raise ValueError(f"invalid or repeated ordinary edge {pair}")
        seen.add(pair)
        value = _qfield(entries[str(m)])
        if sp.simplify(value - exact_cosines[m]) != 0:
            raise ValueError(f"wrong exact m={m} Gram entry")
        gram[i - 1, j - 1] = gram[j - 1, i - 1] = value

    dotted_seen: set[tuple[int, int]] = set()
    for edge in datum["dotted_edges"]:
        i, j = map(int, edge["facets"])
        pair = tuple(sorted((i, j)))
        if pair in seen or pair in dotted_seen:
            raise ValueError(f"invalid or repeated dotted edge {pair}")
        dotted_seen.add(pair)
        weight = _qfield(edge["positive_weight"])
        if weight.is_positive is not True or sp.simplify(weight - 1).is_positive is not True:
            raise ValueError(f"dotted weight is not exactly greater than one: {pair}")
        gram[i - 1, j - 1] = gram[j - 1, i - 1] = -weight
    if dotted_seen != DOTTED_PAIRS:
        raise ValueError("dotted edge stock differs from the pinned P_{6,10} diagram")
    return gram


def _exact_nonsingular_inertia(matrix: sp.Matrix) -> tuple[int, int, list[sp.Expr]]:
    """Inertia from exact scalar Schur pivots; no floating point is used."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    pivots: list[sp.Expr] = []
    while work.rows:
        pivot = sp.radsimp(sp.simplify(work[0, 0]))
        if pivot == 0:
            raise ValueError("chosen principal block has a zero congruence pivot")
        if pivot.is_positive is True:
            positive += 1
        elif pivot.is_negative is True:
            negative += 1
        else:
            raise ValueError(f"SymPy could not decide exact algebraic sign of {pivot}")
        pivots.append(pivot)
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = (work[1:, 1:] - column * column.T / pivot).applyfunc(
            lambda value: sp.radsimp(sp.simplify(value))
        )
    return positive, negative, pivots


def _zphi_coefficients(value: sp.Expr) -> tuple[int, int]:
    """Return exact integral ``a,b`` with value=a+b*phi, or reject it."""

    q2, q5 = sp.sqrt(2), sp.sqrt(5)
    value = sp.radsimp(sp.expand(value))
    if value.has(q2, sp.sqrt(10)):
        raise ValueError(f"value did not descend to Q(sqrt(5)): {value}")
    conjugate = value.xreplace({q5: -q5})
    b = sp.simplify((value - conjugate) / q5)
    a = sp.simplify(value - b * (1 + q5) / 2)
    if a.is_Integer is not True or b.is_Integer is not True:
        raise ValueError(f"value is not in Z[phi]: {value} = {a} + {b} phi")
    return int(a), int(b)


def analyze() -> dict[str, Any]:
    payload = load_artifact()
    gram = build_unit_gram(payload)
    unit = payload["unit_gram"]

    rank = gram.rank()
    kernel = gram.nullspace()
    principal = gram[:7, :7]
    positive, negative, pivots = _exact_nonsingular_inertia(principal)
    schur = (
        gram[7:, 7:] - gram[:7, 7:].T * principal.inv() * gram[:7, 7:]
    ).applyfunc(lambda value: sp.radsimp(sp.simplify(value)))

    dotted = {
        tuple(edge["facets"]): _qfield(edge["positive_weight"])
        for edge in unit["dotted_edges"]
    }
    w67, w78 = dotted[(6, 7)], dotted[(7, 8)]

    golden = payload["golden_crystallographic_scaling"]
    d = list(map(int, golden["d"]))
    diagonal_scale = sp.diag(*[sp.sqrt(value) for value in d])
    mixed = (2 * diagonal_scale * gram * diagonal_scale).applyfunc(
        lambda value: sp.radsimp(sp.expand(value))
    )
    mixed_coefficients = [
        [_zphi_coefficients(mixed[i, j]) for j in range(10)] for i in range(10)
    ]
    reflection_coefficients = [
        [_zphi_coefficients(sp.simplify(mixed[i, j] / d[i])) for j in range(10)]
        for i in range(10)
    ]

    expected_ordinary = {
        int(label): tuple(value)
        for label, value in golden["ordinary_entry_by_m_as_a_plus_b_phi"].items()
    }
    expected_dotted = {
        tuple(edge["facets"]): tuple(edge["entry"])
        for edge in golden["dotted_entries_as_a_plus_b_phi"]
    }
    edge_scaling_ok = all(
        mixed_coefficients[i - 1][j - 1] == expected_ordinary[m]
        for i, j, m in ORDINARY_EDGES
    ) and all(
        mixed_coefficients[i - 1][j - 1] == value
        for (i, j), value in expected_dotted.items()
    )

    geometry = payload["geometry"]
    lineage = payload["index_two_lineage"]
    f_vector = list(map(int, geometry["f_vector_vertices_to_polytope"]))
    chi = sp.Rational(geometry["coxeter_orbifold_euler_characteristic"])
    parent_chi = sp.Rational(lineage["parent_coxeter_orbifold_euler_characteristic"])
    volume = sp.sympify(geometry["volume"], locals={"pi": sp.pi})

    checks = {
        "schema_and_id": payload["schema_version"] == "0.1.0"
        and payload["id"] == "bugaenko_compact_h6_index_two_double_p6_10",
        "field_basis": unit["field"] == "Q(sqrt(2),sqrt(5))"
        and unit["coordinate_basis"] == ["1", "sqrt(2)", "sqrt(5)", "sqrt(10)"],
        "rank": rank == unit["expected_rank"] == 7,
        "kernel": len(kernel) == unit["expected_kernel_dimension"] == 3
        and all(
            all(sp.radsimp(sp.simplify(value)) == 0 for value in gram * vector)
            for vector in kernel
        ),
        "inertia": [positive, negative, 10 - rank]
        == unit["expected_inertia_positive_negative_zero"]
        == [6, 1, 3],
        "zero_schur_complement": schur == sp.zeros(3),
        "dotted_weight_relation": sp.simplify(w78 - (w67**2 - 1)) == 0,
        "mixed_diagonal": [mixed[i, i] for i in range(10)] == golden["diagonal"],
        "mixed_edge_values": edge_scaling_ok,
        "mixed_entries_integral": len(mixed_coefficients) == 10,
        "reflection_coefficients_integral": len(reflection_coefficients) == 10,
        "f_vector": len(f_vector) == 7 and f_vector[5:] == [10, 1]
        and sum((-1) ** i * value for i, value in enumerate(f_vector)) == 1,
        "index_two_euler_and_volume": chi == 2 * parent_chi
        and sp.simplify(volume + chi * sp.Rational(8, 15) * sp.pi**3) == 0,
        "lineage_mapping": lineage["child_facet_to_index2_archive_node"] == MAPPING
        and len(set(MAPPING.values())) == 10,
        "lineage_hashes": lineage["sha256"] == HASHES,
        "lineage_operation": lineage["reflection_group_index"] == 2
        and lineage["doubling_wall_only_finite_nonorthogonal_label"] == 4
        and lineage["bugaenko_to_tumarkin_parent_node_permutation"]
        == [1, 2, 3, 4, 5, 6, 7, 9, 8],
        "source_pin": payload["source_pin"]["commit"]
        == "7a6dcc1a7918985c19fd159c8278042e8f51495b",
        "compactness_scope_honest": payload["compactness_provenance"]["status"]
        == "accepted_known_compact_example"
        and payload["compactness_provenance"]["local_certificate_scope"]
        == "Gram, algebraic scaling, and lineage consistency only",
        "uniqueness_unaccepted": payload["classification_boundary"][
            "arXiv_2608_02894v1_uniqueness"
        ] == "unaccepted_open_due_to_invalid_facet_filter"
        and payload["classification_boundary"]["nominal_combinatorial_types"] == 387
        and payload["classification_boundary"]["types_searched_after_disputed_filters"] == 11
        and payload["classification_boundary"]["worldwide_new_abstract_polytope_claim"] is False,
    }
    return {
        "id": payload["id"],
        "rank": rank,
        "kernel_dimension": len(kernel),
        "inertia_positive_negative_zero": [positive, negative, 10 - rank],
        "principal_block_determinant": str(sp.radsimp(principal.det())),
        "congruence_pivots": [str(pivot) for pivot in pivots],
        "mixed_diagonal": [int(mixed[i, i]) for i in range(10)],
        "compactness_status": payload["compactness_provenance"]["status"],
        "uniqueness_status": payload["classification_boundary"][
            "arXiv_2608_02894v1_uniqueness"
        ],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
