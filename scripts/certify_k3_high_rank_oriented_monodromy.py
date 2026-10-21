#!/usr/bin/env python3
"""Certify the oriented embedding and arithmetic period-group layer of K3RootFV.

The input is the frozen 19-row high-rank complete-root classification.  This
certificate does *not* present any of the infinite orthogonal groups.  It
checks the finite discriminant modules, counts their orthogonal groups, and
records the exact lattice-theoretic consequences of three source theorems:

* primitive embedding uniqueness in the K3 lattice;
* the discriminant-kernel description of the pointwise lattice-polarized K3
  arithmetic change-of-marking group; and
* surjectivity O(L) -> O(q_L) when the indefinite even lattice L splits U.

The last hypothesis is available on both sides of every row: N visibly splits
U, while Roulleau's irreducibility proof supplies U in every complement T.
The positive reflection in that U is also checked explicitly; it is what
allows every full-orthogonal statement to be sharpened to the relevant
orientation-preserving statement.
"""

from __future__ import annotations

import argparse
import functools
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
INVENTORY_PATH = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)
HIGH_CHAMBERS_PATH = (
    ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
)
HIGH_VINBERG_PATH = (
    ROOT / "data/known_examples/k3_high_rank_nef_vinberg_volume.json"
)
OUTPUT_PATH = (
    ROOT / "data/enumerations/k3_high_rank_oriented_monodromy.json"
)
LITERAL_VINBERG_PATH = (
    ROOT / "data/enumerations/k3_rank11_14_literal_vinberg_volume.json"
)
ROULLEAU_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
DOLGACHEV_SOURCE_GZ = Path("/private/tmp/dolgachev_mirror_9502005.tar")

PINNED_SOURCE_HASHES = {
    "roulleau_tex": (
        "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
    ),
    "dolgachev_arxiv_source_gzip": (
        "4bc625d0854b5c1c304d098e7247fd54597e0c781ea48c4320aa0d9f8b091baf"
    ),
}

# The root factors following the displayed U in each N.  E8 is included even
# though it contributes no discriminant group, because the reconstructed Gram
# matrix is used to check rank, determinant, and Smith form.
CASE_FACTORS: dict[str, tuple[str, ...]] = {
    "R11_UE8A1": ("E8", "A1"),
    "R11_UD8A1": ("D8", "A1"),
    "R11_UD4D4A1": ("D4", "D4", "A1"),
    "R11_UD4A1x5": ("D4", "A1", "A1", "A1", "A1", "A1"),
    "R12_UA2E8": ("A2", "E8"),
    "R12_UE8A1x2": ("E8", "A1", "A1"),
    "R12_UD8A1x2": ("D8", "A1", "A1"),
    "R12_UD4x2A1x2": ("D4", "D4", "A1", "A1"),
    "R13_UE8A3": ("E8", "A3"),
    "R13_UE8A1x3": ("E8", "A1", "A1", "A1"),
    "R13_UD8A1x3": ("D8", "A1", "A1", "A1"),
    "R14_UE8D4": ("E8", "D4"),
    "R14_UD8D4": ("D8", "D4"),
    "R14_UE8A1x4": ("E8", "A1", "A1", "A1", "A1"),
    "R15_UE8D4A1": ("E8", "D4", "A1"),
    "R16_UE8D6": ("E8", "D6"),
    "R17_UE8E7": ("E8", "E7"),
    "R18_UE8E8": ("E8", "E8"),
    "R19_UE8E8A1": ("E8", "E8", "A1"),
}

# For the 2-elementary cases, 2*q takes values in Z/4.  These are its values
# on the natural discriminant coordinates of the negative ADE factor.  A
# D_{2m} coordinate uses the two spinor classes; their sum is the vector class.
Q4_FACTOR_TABLE = {
    "A1": (0, 3),
    "D4": (0, 2, 2, 2),
    "D6": (0, 1, 1, 2),
    "D8": (0, 0, 0, 2),
    "E7": (0, 1),
}

# Frozen regression orders.  The certifier recomputes every entry from the
# standard ADE finite-quadratic-form tables encoded above.  These constants are
# a regression boundary, not an independent inverse-Gram derivation of q_N.
EXPECTED_DISCRIMINANT_ORTHOGONAL_ORDERS = {
    "R11_UE8A1": 1,
    "R11_UD8A1": 2,
    "R11_UD4D4A1": 72,
    "R11_UD4A1x5": 40320,
    "R12_UA2E8": 2,
    "R12_UE8A1x2": 2,
    "R12_UD8A1x2": 12,
    "R12_UD4x2A1x2": 1440,
    "R13_UE8A3": 2,
    "R13_UE8A1x3": 6,
    "R13_UD8A1x3": 120,
    "R14_UE8D4": 6,
    "R14_UD8D4": 120,
    "R14_UE8A1x4": 24,
    "R15_UE8D4A1": 6,
    "R16_UE8D6": 2,
    "R17_UE8E7": 1,
    "R18_UE8E8": 1,
    "R19_UE8E8A1": 1,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _compact_lattice_name(name: str) -> str:
    compact = name.replace(" ", "")
    return compact.replace("E8(-1)^2", "E8(-1)+E8(-1)")


def _cartan_a(n: int) -> sp.Matrix:
    gram = 2 * sp.eye(n)
    for i in range(n - 1):
        gram[i, i + 1] = gram[i + 1, i] = -1
    return gram


def _cartan_d(n: int) -> sp.Matrix:
    if n < 4:
        raise ValueError(n)
    gram = 2 * sp.eye(n)
    edges = [(0, 2), (1, 2)] + [(i, i + 1) for i in range(2, n - 1)]
    for i, j in edges:
        gram[i, j] = gram[j, i] = -1
    return gram


def _cartan_e(n: int) -> sp.Matrix:
    if n not in {6, 7, 8}:
        raise ValueError(n)
    gram = 2 * sp.eye(n)
    edges = [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)]
    for i, j in edges:
        gram[i, j] = gram[j, i] = -1
    return gram


def _root_gram(name: str) -> sp.Matrix:
    family, rank = name[0], int(name[1:])
    if family == "A":
        return _cartan_a(rank)
    if family == "D":
        return _cartan_d(rank)
    if family == "E":
        return _cartan_e(rank)
    raise ValueError(name)


def _block_diag(blocks: list[sp.Matrix]) -> sp.Matrix:
    return sp.diag(*blocks)


def _ns_gram(factors: tuple[str, ...]) -> sp.Matrix:
    u = sp.Matrix([[0, 1], [1, 0]])
    return _block_diag([u] + [-_root_gram(factor) for factor in factors])


def _smith_nonunit_invariants(gram: sp.Matrix) -> list[int]:
    smith = smith_normal_form(gram, domain=ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(smith.rows)
        if abs(int(smith[i, i])) > 1
    ]


def _factor_discriminant_label(factor: str) -> str | None:
    labels = {
        "A1": "(Z/2,q=-1/2)",
        "A2": "(Z/3,q=-2/3)",
        "A3": "(Z/4,q=-3/4)",
        "D4": "q_D4(-1) on (Z/2)^2",
        "D6": "q_D6(-1) on (Z/2)^2",
        "D8": "q_D8(-1) on (Z/2)^2",
        "E7": "(Z/2,q=-3/2)",
        "E8": None,
    }
    return labels[factor]


def _q4_form(factors: tuple[str, ...]) -> tuple[int, Callable[[int], int]]:
    pieces: list[tuple[int, tuple[int, ...]]] = []
    offset = 0
    for factor in factors:
        if factor == "E8":
            continue
        table = Q4_FACTOR_TABLE.get(factor)
        if table is None:
            raise ValueError(f"not a 2-elementary factor: {factor}")
        width = 1 if len(table) == 2 else 2
        pieces.append((offset, table))
        offset += width

    def q4(vector: int) -> int:
        value = 0
        for start, table in pieces:
            width = 1 if len(table) == 2 else 2
            value += table[(vector >> start) & ((1 << width) - 1)]
        return value % 4

    return offset, q4


def _insert_binary_basis(basis: tuple[int, ...], vector: int) -> tuple[int, ...]:
    rows = list(basis)
    value = vector
    while value:
        pivot = value.bit_length() - 1
        matching = next(
            (i for i, row in enumerate(rows) if row.bit_length() - 1 == pivot),
            None,
        )
        if matching is None:
            rows.append(value)
            break
        value ^= rows[matching]
    rows.sort(reverse=True)
    return tuple(rows)


def _binary_basis_accepts(basis: tuple[int, ...], vector: int) -> bool:
    value = vector
    for row in basis:
        value = min(value, value ^ row)
    return value != 0


@functools.lru_cache(maxsize=None)
def _two_elementary_orthogonal_order(factors: tuple[str, ...]) -> int:
    dimension, q4 = _q4_form(factors)
    if dimension == 0:
        return 1
    standard = [1 << i for i in range(dimension)]
    standard_q = [q4(vector) for vector in standard]

    def b2(left: int, right: int) -> int:
        difference = (q4(left ^ right) - q4(left) - q4(right)) % 4
        assert difference in {0, 2}
        return difference // 2

    target_b = [
        [b2(standard[i], standard[j]) for j in range(dimension)]
        for i in range(dimension)
    ]
    candidates = [
        [
            vector
            for vector in range(1, 1 << dimension)
            if q4(vector) == standard_q[i]
        ]
        for i in range(dimension)
    ]

    count = 0

    def extend(
        index: int,
        images: tuple[int, ...],
        row_basis: tuple[int, ...],
    ) -> None:
        nonlocal count
        if index == dimension:
            count += 1
            return
        for vector in candidates[index]:
            if not _binary_basis_accepts(row_basis, vector):
                continue
            if any(
                b2(vector, images[j]) != target_b[index][j]
                for j in range(index)
            ):
                continue
            extend(
                index + 1,
                images + (vector,),
                _insert_binary_basis(row_basis, vector),
            )

    extend(0, (), ())
    return count


def _discriminant_orthogonal_order(factors: tuple[str, ...]) -> int:
    nontrivial = tuple(factor for factor in factors if factor != "E8")
    if not nontrivial:
        return 1
    if nontrivial == ("A2",):
        # Both units of Z/3 square to one, so both preserve q=-2/3.
        return 2
    if nontrivial == ("A3",):
        # The two units of Z/4 are +/-1 and preserve q=-3/4.
        return 2
    return _two_elementary_orthogonal_order(nontrivial)


def _u_positive_reflection_certificate() -> dict[str, Any]:
    gram = sp.Matrix([[0, 1], [1, 0]])
    h = sp.Matrix([1, 1])
    square = int((h.T * gram * h)[0])
    reflection = sp.eye(2) - h * (h.T * gram)
    assert square == 2
    assert reflection.T * gram * reflection == gram
    assert reflection * reflection == sp.eye(2)
    assert int(reflection.det()) == -1
    return {
        "U_gram": [[0, 1], [1, 0]],
        "positive_vector_h": [1, 1],
        "h_square": square,
        "reflection_matrix": [
            [int(reflection[i, j]) for j in range(2)] for i in range(2)
        ],
        "reflection_determinant": -1,
        "integral_isometry": True,
        "discriminant_action": "trivial because U is unimodular",
        "orientation_effect": (
            "reverses one positive direction; it exchanges the two positive-"
            "cone/type-IV orientation components in the relevant factor"
        ),
    }


def check_pinned_sources() -> dict[str, Any]:
    if not ROULLEAU_TEX.exists():
        raise FileNotFoundError(ROULLEAU_TEX)
    if not DOLGACHEV_SOURCE_GZ.exists():
        raise FileNotFoundError(DOLGACHEV_SOURCE_GZ)
    assert _sha256(ROULLEAU_TEX) == PINNED_SOURCE_HASHES["roulleau_tex"]
    assert (
        _sha256(DOLGACHEV_SOURCE_GZ)
        == PINNED_SOURCE_HASHES["dolgachev_arxiv_source_gzip"]
    )

    roulleau_lines = ROULLEAU_TEX.read_text().splitlines()
    uniqueness = "\n".join(roulleau_lines[1140:1170])
    split_u_general = "\n".join(roulleau_lines[1180:1200])
    split_u_rank19 = "\n".join(roulleau_lines[1208:1229])
    assert "unique up to automorphisms" in uniqueness
    assert "contains a copy of $U$" in split_u_general
    assert "orthogonal complement contains" in split_u_rank19
    assert "a copy of~$U$" in split_u_rank19

    dolgachev = gzip.decompress(DOLGACHEV_SOURCE_GZ.read_bytes()).decode(
        errors="replace"
    )
    assert "Gamma_M = O(N)^*." in dolgachev
    assert "sigma (m) = m" in dolgachev
    assert "D_M$ is the union of two copies" in dolgachev
    return {
        "checked": True,
        "all_checks_pass": True,
        "sha256": {
            "roulleau_tex": _sha256(ROULLEAU_TEX),
            "dolgachev_arxiv_source_gzip": _sha256(DOLGACHEV_SOURCE_GZ),
        },
        "roulleau_witnesses": {
            "full_O_embedding_uniqueness_lines": [1142, 1168],
            "rank_11_18_complement_splits_U_lines": [1184, 1195],
            "rank_19_complement_splits_U_lines": [1221, 1227],
        },
        "dolgachev_witness": {
            "arxiv": "alg-geom/9502005",
            "result": "Proposition 3.3",
            "pointwise_polarization_change_of_marking_group": "discriminant kernel",
        },
        "nikulin_witness": {
            "status": "bibliographic_primary_web_input_not_locally_hash_pinned",
            "result": "Theorem 1.14.2",
            "primary_web": "https://www.mathnet.ru/eng/im1677",
            "doi": "10.1070/IM1980v014n01ABEH001060",
        },
    }


def _discriminant_label(factors: tuple[str, ...]) -> str:
    labels = [
        label
        for factor in factors
        if (label := _factor_discriminant_label(factor)) is not None
    ]
    return "0 (unimodular)" if not labels else " orthogonal_sum ".join(labels)


def build_artifact(source_check: dict[str, Any] | None = None) -> dict[str, Any]:
    inventory = _json(INVENTORY_PATH)
    high_chambers = _json(HIGH_CHAMBERS_PATH)
    high_vinberg = _json(HIGH_VINBERG_PATH)
    literal_vinberg = _json(LITERAL_VINBERG_PATH)
    assert inventory["all_checks_pass"] is True
    assert inventory["denominator"]["rank_at_least_11_stock_count"] == 19
    assert len(inventory["stocks"]) == 19
    assert {stock["key"] for stock in inventory["stocks"]} == set(CASE_FACTORS)
    assert high_chambers["artifact_status"] == "source_audited_gram_verified"
    assert {int(case["rho"]) for case in high_chambers["cases"]} == set(
        range(15, 20)
    )
    assert {int(case["rho"]) for case in high_vinberg["cases"]} == set(
        range(15, 20)
    )
    assert all(case["all_checks_pass"] for case in high_vinberg["cases"])
    assert all(
        case["decision"] == "finite_volume_noncompact"
        for case in high_vinberg["cases"]
    )
    assert literal_vinberg["all_checks_pass"] is True
    assert literal_vinberg["checks"]["all_literal_vinberg_audits_pass"] is True
    assert literal_vinberg["checks"]["all_finite_volume_noncompact"] is True
    assert literal_vinberg["totals"]["case_count"] == 14
    assert {case["key"] for case in literal_vinberg["cases"]} == {
        stock["key"]
        for stock in inventory["stocks"]
        if int(stock["picard_rank"]) <= 14
    }
    source_by_key = {stock["key"]: stock for stock in inventory["stocks"]}
    literal_by_key = {case["key"]: case for case in literal_vinberg["cases"]}
    high_chamber_by_rank = {
        int(case["rho"]): case for case in high_chambers["cases"]
    }
    high_vinberg_by_rank = {
        int(case["rho"]): case for case in high_vinberg["cases"]
    }
    reflection = _u_positive_reflection_certificate()

    rows: list[dict[str, Any]] = []
    for key, factors in CASE_FACTORS.items():
        base = source_by_key[key]
        gram = _ns_gram(factors)
        rho = int(base["picard_rank"])
        complement_rank = 22 - rho
        discriminant_invariants = _smith_nonunit_invariants(gram)
        discriminant_order = abs(int(gram.det()))
        discriminant_length = len(discriminant_invariants)
        orthogonal_order = _discriminant_orthogonal_order(factors)

        assert gram.rows == rho
        assert discriminant_order == base["absolute_discriminant"]
        assert orthogonal_order == EXPECTED_DISCRIMINANT_ORTHOGONAL_ORDERS[key]
        assert discriminant_length + 2 <= complement_rank
        if rho <= 14:
            literal_case = literal_by_key[key]
            assert int(literal_case["rho"]) == rho
            assert literal_case["facet_count"] == base["wall_count"]
            assert _compact_lattice_name(literal_case["lattice"]) == _compact_lattice_name(
                base["lattice"]
            )
            assert literal_case["all_checks_pass"] is True
            assert literal_case["decision"] == "finite_volume_noncompact"
        else:
            high_case = high_chamber_by_rank[rho]
            volume_case = high_vinberg_by_rank[rho]
            assert int(high_case["expected_rank"]) == rho
            assert high_case["wall_count"] == base["wall_count"]
            assert high_case["absolute_discriminant"] == base["absolute_discriminant"]
            assert _compact_lattice_name(high_case["lattice"]) == _compact_lattice_name(
                base["lattice"]
            )
            assert volume_case["facet_count"] == base["wall_count"]
        # The stronger U-splitting bound is exactly what Roulleau uses for
        # ranks <=18; rank 19 is its separately checked exception.
        if rho <= 18:
            assert discriminant_length <= 18 - rho
            complement_u_witness = "Roulleau_corollary_embedding_leaves_U"
        else:
            assert key == "R19_UE8E8A1"
            complement_u_witness = "Roulleau_rank19_direct_complement_contains_U"

        rows.append(
            {
                "key": key,
                "picard_rank": rho,
                "hyperbolic_dimension": rho - 1,
                "lattice": base["lattice"],
                "facet_count": base["wall_count"],
                "N_signature": [1, rho - 1],
                "T_rank": complement_rank,
                "T_signature": [2, 20 - rho],
                "N_splits_U": True,
                "T_splits_U": True,
                "T_splits_U_witness": complement_u_witness,
                "discriminant_group_invariants": discriminant_invariants,
                "absolute_discriminant": discriminant_order,
                "discriminant_length": discriminant_length,
                "q_N": _discriminant_label(factors),
                "q_T": "-q_N under the primitive unimodular gluing anti-isometry",
                "q_form_provenance": (
                    "standard ADE discriminant-form table encoded in the "
                    "certifier; not independently derived from the inverse Gram"
                ),
                "discriminant_form_orthogonal_group_order": orthogonal_order,
                "primitive_embedding_orbits": {
                    "under_full_O_Lambda_K3": 1,
                    "under_O_plus_Lambda_K3": 1,
                    "O_plus_no_split_reason": (
                        "compose a wrong-orientation transporter with the "
                        "pointwise complement-U positive reflection"
                    ),
                },
                "positive_cone_components": {
                    "number": 2,
                    "orbits_under_O_N": 1,
                    "orbits_under_O_plus_N": 2,
                    "stabilizer_of_either_component": "O^+(N)",
                },
                "complete_root_chambers_in_one_positive_component": {
                    "raw_cardinality": "countably infinite W_2(N)-torsor",
                    "W_2_action": "simply_transitive",
                    "orbits_under_W_2_N": 1,
                    "orbits_under_O_plus_N": 1,
                    "W_2_discriminant_action": "trivial",
                    "arithmetic_chamber_stabilizer": "Gamma_C=Stab_{O^+(N)}(C)",
                },
                "discriminant_representations": {
                    "O_plus_N_to_O_qN": "surjective",
                    "O_plus_T_to_O_qT": "surjective",
                    "surjectivity_reason": (
                        "Nikulin Theorem 1.14.2 plus U splitting; the explicit "
                        "trivial-discriminant positive reflection corrects orientation"
                    ),
                    "stable_index_N": orthogonal_order,
                    "stable_index_T": orthogonal_order,
                },
                "oriented_period_change_of_marking_data": {
                    "pointwise_N_polarized_T_period_image": (
                        "tilde O^+(T)=ker(O^+(T)->O(q_T))"
                    ),
                    "pointwise_discriminant_image": "identity",
                    "unmarked_ambient_chamber_normalizer": "Gamma^+_{j,C}",
                    "unmarked_chamber_normalizer_T_period_image": "O^+(T)",
                    "unmarked_T_discriminant_image": "all O(q_T)",
                    "arithmetic_group_index": orthogonal_order,
                    "coarse_period_quotient_degree": (
                        "not asserted: central isometries such as -id can act "
                        "ineffectively on the projectivized period domain"
                    ),
                    "ambient_normalizer_exact_sequence": (
                        "1 -> tilde O^+(T) -> Gamma^+_{j,C} -> "
                        "Stab_{O^+(N)}(C) -> 1"
                    ),
                    "geometric_monodromy_identification": (
                        "conditional_only_on_the_standard_moduli_stack_or_deck_group_"
                        "identification; no independent geometric-monodromy packet "
                        "is certified here"
                    ),
                },
                "generic_surface_automorphism_group_source_label": base[
                    "general_surface_automorphism_group"
                ],
                "separation_warning": (
                    "the source automorphism group is comparison-only: it is "
                    "neither the weighted-diagram group nor the oriented arithmetic "
                    "period/change-of-marking group; exact geometric monodromy is "
                    "not asserted"
                ),
            }
        )

    order_histogram: dict[str, int] = {}
    for row in rows:
        order = str(row["discriminant_form_orthogonal_group_order"])
        order_histogram[order] = order_histogram.get(order, 0) + 1

    return {
        "schema": "coxeter.k3_high_rank.oriented_period_groups.v2",
        "artifact_id": "K3-HIGH-RANK-ORIENTED-PERIOD-GROUPS-19",
        "all_checks_pass": True,
        "scope": {
            "case_count": 19,
            "picard_ranks": [11, 19],
            "input_class": "K3RootFV[11,20] finite-volume rows",
            "claim_level": (
                "exact lattice-theoretic orientation, embedding-orbit, chamber-"
                "orbit, and oriented arithmetic period/change-of-marking image "
                "classification"
            ),
        },
        "source_artifacts": {
            "dependency_graph": (
                "acyclic upstream-only: this artifact does not read or hash any "
                "downstream synthesis"
            ),
            "roulleau_19_stock_inventory": str(INVENTORY_PATH.relative_to(ROOT)),
            "roulleau_19_stock_inventory_sha256": _sha256(INVENTORY_PATH),
            "rank15_19_chamber_manifest": str(HIGH_CHAMBERS_PATH.relative_to(ROOT)),
            "rank15_19_chamber_manifest_sha256": _sha256(HIGH_CHAMBERS_PATH),
            "rank15_19_literal_vinberg_manifest": str(
                HIGH_VINBERG_PATH.relative_to(ROOT)
            ),
            "rank15_19_literal_vinberg_manifest_sha256": _sha256(
                HIGH_VINBERG_PATH
            ),
            "rank11_14_literal_vinberg_manifest": str(
                LITERAL_VINBERG_PATH.relative_to(ROOT)
            ),
            "rank11_14_literal_vinberg_manifest_sha256": _sha256(
                LITERAL_VINBERG_PATH
            ),
            "nikulin_surjectivity_input": {
                "citation": (
                    "V. V. Nikulin, Integral symmetric bilinear forms and some "
                    "of their applications, Math. USSR-Izv. 14 (1980), "
                    "Theorem 1.14.2"
                ),
                "doi": "10.1070/IM1980v014n01ABEH001060",
                "provenance": "bibliographic_primary_web_input_not_locally_hash_pinned",
                "primary_web": "https://www.mathnet.ru/eng/im1677",
                "used_special_case": (
                    "an even indefinite lattice that splits U has surjective "
                    "O(L)->O(q_L)"
                ),
            },
            "source_check": source_check
            or {
                "checked": False,
                "all_checks_pass": False,
                "reason": "build called without external pinned-source replay",
            },
        },
        "uniform_theorem": {
            "primitive_embedding_orbits_full_O": 1,
            "primitive_embedding_orbits_O_plus": 1,
            "positive_cone_component_orbits_under_O": 1,
            "positive_cone_component_orbits_under_O_plus": 2,
            "chamber_orbits_in_fixed_component_under_W2": 1,
            "chamber_orbits_in_fixed_component_under_O_plus": 1,
            "pointwise_N_polarized_T_period_image": "tilde O^+(T)",
            "pointwise_discriminant_image": "identity",
            "unmarked_ambient_group": "Gamma^+_{j,C}",
            "unmarked_chamber_normalizer_T_period_image": "O^+(T)",
            "unmarked_T_discriminant_image": "all O(q_T)",
            "geometric_monodromy_status": (
                "conditional on the standard moduli-stack/deck-group "
                "identification; not independently certified"
            ),
            "proof_dependencies": [
                "full-O primitive embedding uniqueness",
                "U splits both N and T",
                "Nikulin discriminant-representation surjectivity for U-split lattices",
                "primitive unimodular gluing extension criterion",
                "W_2 acts simply transitively on complete root chambers",
                "every (-2)-reflection acts trivially on the discriminant group",
                "Dolgachev Proposition 3.3 for pointwise lattice polarization",
            ],
        },
        "marking_conventions": {
            "fixed_source_marking": (
                "the two positive-cone components and every individual chamber "
                "remain distinct labels; the chamber set is infinite"
            ),
            "oriented_integral_isometry": (
                "components are isomorphic via an isometry carrying the chosen "
                "component to the chosen component, and all chambers in one "
                "component form one W_2 orbit"
            ),
            "ambient_O_plus_embedding": (
                "primitive embeddings do not split into two orbits because the "
                "complement-U reflection fixes N pointwise and reverses ambient "
                "positive-three-space orientation"
            ),
            "type_IV_components": (
                "O^+(T) preserves each of the two period-domain components; O(T) "
                "is transitive on them. Fixing the N positive-cone component and "
                "the ambient K3 orientation selects the compatible period component"
            ),
        },
        "u_positive_reflection": reflection,
        "finite_discriminant_quotients": {
            "order_histogram": dict(
                sorted(order_histogram.items(), key=lambda item: int(item[0]))
            ),
            "maximum_order": max(EXPECTED_DISCRIMINANT_ORTHOGONAL_ORDERS.values()),
            "maximum_order_case": "R11_UD4A1x5",
            "interpretation": (
                "[O^+(T):tilde O^+(T)] and the corresponding Picard-side "
                "stable-chamber-stabilizer index"
            ),
        },
        "finite_quadratic_form_provenance": {
            "ADE_q_forms": (
                "table-derived from the standard ADE discriminant-form values "
                "encoded in Q4_FACTOR_TABLE and _factor_discriminant_label"
            ),
            "smith_invariants": (
                "independently recomputed from the integral U plus negative-ADE "
                "Gram matrices"
            ),
            "orthogonal_group_orders": (
                "algorithmically enumerated from the table-derived finite forms; "
                "EXPECTED_DISCRIMINANT_ORTHOGONAL_ORDERS are regression values, "
                "not a second inverse-Gram derivation"
            ),
        },
        "literal_vinberg_status": {
            "all_19_closed": True,
            "ranks_11_14_case_count": 14,
            "ranks_15_19_case_count": 5,
            "new_rank11_14_report": "results/k3_rank11_14_literal_vinberg_volume.md",
            "remaining_literal_V_gap_count": 0,
        },
        "rows": rows,
        "presentation_boundary": {
            "explicit_finite_presentations_supplied": False,
            "reason": (
                "the local pinned sources, bibliographic Nikulin input, and this "
                "table-driven finite discriminant computation do "
                "not derive generators-and-relations presentations for the 19 "
                "infinite groups tilde O^+(T), or an independent geometric-"
                "monodromy packet; naming the arithmetic T-image as an exact "
                "discriminant kernel is not a presentation"
            ),
            "not_conflated": [
                "weighted Coxeter-diagram automorphisms",
                "integral chamber stabilizers",
                "generic single-surface automorphisms",
                "oriented arithmetic period/change-of-marking images",
                "geometric family monodromy",
                "unmarked arithmetic period groups",
            ],
        },
        "strict_M3c_residual": {
            "status": "open",
            "reasons": [
                "the ambient class is not the bounded raw IC0 specification",
                "no explicit finite presentations of the 19 infinite arithmetic groups",
                "no independently certified geometric-monodromy packet beyond the conditional standard stack identification",
            ],
            "literal_V_status": "closed for all 19 and no longer a residual gap",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--check-sources", action="store_true")
    args = parser.parse_args()

    frozen = _json(OUTPUT_PATH) if OUTPUT_PATH.exists() else None
    if args.check_sources:
        source_check = check_pinned_sources()
    elif frozen is not None:
        source_check = frozen["source_artifacts"]["source_check"]
    else:
        source_check = None
    artifact = build_artifact(source_check)

    if args.check_archive:
        if frozen is None:
            raise FileNotFoundError(OUTPUT_PATH)
        if artifact != frozen:
            raise SystemExit("oriented-period-group archive differs from exact replay")
        print("archive_check=PASS")
    else:
        OUTPUT_PATH.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
        print(f"wrote={OUTPUT_PATH}")
    print("case_count=19")
    print("O_plus_embedding_orbits_per_case=1")
    print("pointwise_T_period_image=tilde_O_plus_T")
    print("unmarked_chamber_normalizer_T_period_image=O_plus_T")
    print(
        "maximum_discriminant_quotient_order="
        f"{artifact['finite_discriminant_quotients']['maximum_order']}"
    )


if __name__ == "__main__":
    main()
