#!/usr/bin/env python3
"""Certify the common period-moduli K3 carrier criterion for S in H10.

The lattice

    S = U + E8(-2) + <-4>

is the very-general Neron--Severi lattice on a nine-dimensional K3 period
locus.  Hence every genuine H10 Coxeter polyhedron whose primitive negative
facet reflections preserve S is a weak full-Picard divisor carrier after an
S-marking.  This verifier derives the *exact* primitive negative reflection
alphabet from the full discriminant form, not just from its exponent-four
upper bound.  It also pins the complete 13-wall reflection chamber and the
known 123- and 84-wall finite-index calibrations.

This is a period-level carrier theorem.  It supplies neither a common
equation family nor a full-nef or surface-automorphism realization.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/ue8m2a4_common_k3_moduli_carrier.json"
SOURCES = {
    "shared_k3_carrier": (
        ROOT / "data/enumerations/f3_v4_jk_h10_shared_k3_carrier.json"
    ),
    "full_reflection_lineage": (
        ROOT
        / "data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json"
    ),
    "jk_index": (
        ROOT / "data/enumerations/jacobian_kummer_h10_reflection_index.json"
    ),
    "f3_vinberg_index": (
        ROOT
        / "data/enumerations/ss_fermat_quartic_char3_v4_h10_vinberg_certificate.json"
    ),
}
EXPECTED_FILE_SHA256 = {
    "shared_k3_carrier": (
        "b84dcc92f82e10cdde6eed98ef2b599c2985f1ec0bbb91c156c166782fe17291"
    ),
    "full_reflection_lineage": (
        "3cb8b0e9120d13479f9b9e3a32cfb63922e03f41f2f0fc1dc501b0a9f6654c29"
    ),
    "jk_index": (
        "c8b8a6562ae3b99357d293d3e8c9e23d231236bdf2850659ab6f462c4d28e3d8"
    ),
    "f3_vinberg_index": (
        "684e03ca1c1196c657cc707b80f82ac19f84d9510a408b5e8d9b63271b7272be"
    ),
}
EXPECTED_PAYLOAD_SHA256 = {
    "shared_k3_carrier": (
        "39574ca56421155388ec734076cf728a49f81c94bdfcdbc7d472b104e5014f08"
    ),
    "full_reflection_lineage": (
        "40f44ad15d4e7b2c1d3014bd19d43a178042635e4ddf5a851001e08f48d630d9"
    ),
    "jk_index": (
        "46d4f3ce47b5e98c8853139a47d27d3d44d7a8516d02270cfc67ff0c33eaed7d"
    ),
    "f3_vinberg_index": (
        "5d104aac54ae18f6a54e64683db920227fa803df9c2997e88b9410fd8bf97f42"
    ),
}
SCHEMA = "coxeter.k3.ue8m2a4_common_moduli_carrier.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_payload_sha256(payload: dict[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _payload_sha256(payload: dict[str, Any]) -> str:
    unsealed = dict(payload)
    unsealed.pop("payload_sha256", None)
    return _canonical_payload_sha256(unsealed)


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[index, index]))
        for index in range(min(diagonal.shape))
        if diagonal[index, index]
    ]


def _matrix(payload: list[list[int | str]]) -> sp.Matrix:
    return sp.Matrix(
        [[sp.Rational(value) for value in row] for row in payload]
    )


def _mod_two(value: sp.Rational) -> sp.Rational:
    return sp.Rational(value) % 2


def _load_sources() -> dict[str, dict[str, Any]]:
    documents = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }
    for key, path in SOURCES.items():
        if _sha256(path) != EXPECTED_FILE_SHA256[key]:
            raise ArithmeticError(f"pinned {key} file changed")
        if _canonical_payload_sha256(documents[key]) != EXPECTED_PAYLOAD_SHA256[key]:
            raise ArithmeticError(f"pinned {key} payload changed")
    if documents["shared_k3_carrier"].get("all_checks_pass") is not True:
        raise ArithmeticError("shared K3 carrier source is not green")
    if documents["full_reflection_lineage"].get("all_checks_pass") is not True:
        raise ArithmeticError("full reflection lineage source is not green")
    if documents["jk_index"].get("all_checks_pass") is not True:
        raise ArithmeticError("JK index source is not green")
    if not all(documents["f3_vinberg_index"]["checks"].values()):
        raise ArithmeticError("F3 Vinberg/index source is not green")
    return documents


def _primitive_reflection_record(
    gram: sp.Matrix, vector: sp.Matrix
) -> dict[str, Any]:
    content = math.gcd(*(abs(int(value)) for value in vector))
    square = int((vector.T * gram * vector)[0])
    divisibility = math.gcd(*(abs(int(value)) for value in gram * vector))
    reflection = sp.simplify(
        sp.eye(gram.rows)
        - sp.Rational(2, square) * vector * (vector.T * gram)
    )
    dual_class = vector / divisibility
    exact_order = next(
        order
        for order in range(1, divisibility + 1)
        if all(value.is_Integer for value in order * dual_class)
    )
    return {
        "coordinates": [int(value) for value in vector],
        "square": square,
        "divisibility": divisibility,
        "primitive": content == 1,
        "dual_class_exact_order": exact_order,
        "dual_class_q_mod_2": str(
            _mod_two(sp.Rational(square, divisibility * divisibility))
        ),
        "reflection_integral": all(value.is_Integer for value in reflection),
        "reflection_preserves_S": reflection.T * gram * reflection == gram,
        "hyperbolic_mirror_normal": square < 0,
    }


def _exponent_four_upper_bound() -> list[dict[str, int]]:
    candidates: list[dict[str, int]] = []
    for divisibility in (1, 2, 4):
        for absolute_square in range(2, 2 * divisibility + 1, 2):
            if 2 * divisibility % absolute_square == 0:
                candidates.append(
                    {
                        "square": -absolute_square,
                        "divisibility": divisibility,
                    }
                )
    return candidates


def _discriminant_form_audit(gram: sp.Matrix) -> dict[str, Any]:
    """Enumerate A_S using the visible U + E8(-2) + <-4> splitting."""

    U = sp.Matrix([[0, 1], [1, 0]])
    e8_half = sp.simplify(gram[2:10, 2:10] / 2)
    block_checks = {
        "visible_U_block_unimodular": gram[:2, :2] == U,
        "visible_E8_half_integral_even_unimodular": (
            all(value.is_Integer for value in e8_half)
            and abs(int(e8_half.det())) == 1
            and all(e8_half[index, index] % 2 == 0 for index in range(8))
        ),
        "visible_last_block_minus_four": gram[10, 10] == -4,
        "orthogonal_block_splitting": (
            gram[:2, 2:] == sp.zeros(2, 9)
            and gram[2:10, :2] == sp.zeros(8, 2)
            and gram[2:10, 10] == sp.zeros(8, 1)
            and gram[10, :10] == sp.zeros(1, 10)
        ),
    }

    histogram: dict[tuple[int, str], int] = {}
    class_keys: set[tuple[sp.Rational, ...]] = set()
    representatives = 0
    all_are_dual = True
    for bits in itertools.product(range(2), repeat=8):
        for last in range(4):
            alpha = sp.Matrix(
                [
                    0,
                    0,
                    *(sp.Rational(bit, 2) for bit in bits),
                    sp.Rational(last, 4),
                ]
            )
            all_are_dual &= all(value.is_Integer for value in gram * alpha)
            class_keys.add(tuple(sp.Rational(value) % 1 for value in alpha))
            order = next(
                candidate
                for candidate in range(1, 5)
                if all(value.is_Integer for value in candidate * alpha)
            )
            q_value = str(_mod_two(sp.Rational((alpha.T * gram * alpha)[0])))
            histogram[(order, q_value)] = histogram.get((order, q_value), 0) + 1
            representatives += 1

    rendered_histogram = {
        f"order_{order}_q_{q_value}": count
        for (order, q_value), count in sorted(histogram.items())
    }
    available = {
        order: {
            sp.Rational(q_value)
            for (entry_order, q_value), count in histogram.items()
            if entry_order == order and count
        }
        for order in (1, 2, 4)
    }

    candidate_records: list[dict[str, Any]] = []
    for candidate in _exponent_four_upper_bound():
        square = candidate["square"]
        divisibility = candidate["divisibility"]
        self_pairing_compatible = square % divisibility == 0
        required_q = _mod_two(
            sp.Rational(square, divisibility * divisibility)
        )
        discriminant_compatible = required_q in available[divisibility]
        allowed = self_pairing_compatible and discriminant_compatible
        if not self_pairing_compatible:
            status = "excluded_by_d_divides_r_square_before_discriminant_form"
        elif not discriminant_compatible:
            status = "excluded_by_full_discriminant_form"
        else:
            status = "allowed_and_explicitly_realized"
        candidate_records.append(
            {
                **candidate,
                "self_pairing_divisibility_compatible": self_pairing_compatible,
                "required_discriminant_q_mod_2": str(required_q),
                "available_q_mod_2_at_exact_order": [
                    str(value) for value in sorted(available[divisibility])
                ],
                "discriminant_form_compatible": discriminant_compatible,
                "survives_all_arithmetic_filters": allowed,
                "status": status,
            }
        )

    return {
        "representative_model": (
            "alpha=(0,0,b_1/2,...,b_8/2,c/4), with b_i in {0,1} "
            "and c in {0,1,2,3}"
        ),
        "representative_count": representatives,
        "distinct_classes_mod_S": len(class_keys),
        "all_representatives_lie_in_S_dual": all_are_dual,
        "exact_order_q_mod_2_histogram": rendered_histogram,
        "block_checks": block_checks,
        "exponent_four_and_reflection_integrality_naive_upper_bound": (
            _exponent_four_upper_bound()
        ),
        "after_d_divides_self_pairing_filter": [
            {"square": row["square"], "divisibility": row["divisibility"]}
            for row in candidate_records
            if row["self_pairing_divisibility_compatible"]
        ],
        "candidate_filter": candidate_records,
        "proof": (
            "For primitive r with d=div_S(r), r/d has exact order d in A_S, "
            "d divides r^2 because r^2 is one of the pairings, and "
            "q(r/d)=r^2/d^2 mod 2Z. Reflection integrality is equivalent to "
            "|r^2| dividing 2d. The self-pairing divisibility first removes "
            "(-2,4); the enumerated exact-order q-values then remove (-2,2) "
            "and (-8,4). Explicit primitive vectors realize the three survivors."
        ),
    }


def build_certificate() -> dict[str, Any]:
    sources = _load_sources()
    shared = sources["shared_k3_carrier"]
    maximal = sources["full_reflection_lineage"]
    jk_index_source = sources["jk_index"]
    f3_index_source = sources["f3_vinberg_index"]

    lattice_source = shared["shared_Picard_lattice"]
    gram = _matrix(lattice_source["gram"])
    smith = _smith(gram)
    discriminant = _discriminant_form_audit(gram)

    witness_vectors = {
        "minus2_div1": sp.Matrix([1, -1, *([0] * 9)]),
        "minus4_div2": sp.Matrix([0, 0, 1, *([0] * 8)]),
        "minus4_div4": sp.Matrix([*([0] * 10), 1]),
    }
    witnesses = {
        key: _primitive_reflection_record(gram, vector)
        for key, vector in witness_vectors.items()
    }
    positive_boundary = _primitive_reflection_record(
        gram, sp.Matrix([1, 1, *([0] * 9)])
    )
    alphabet = [
        {"square": -2, "divisibility": 1},
        {"square": -4, "divisibility": 2},
        {"square": -4, "divisibility": 4},
    ]

    full = maximal["vinberg_maximal_integral_reflection_chamber"]
    base = jk_index_source["vinberg_f2_h10_13_wall_chamber"]
    chi_full = Fraction(base["rational_euler_characteristic"])
    f3_index = int(
        f3_index_source["exact_euler_volume_cross_check"]
        ["even_subgroup_index_in_f2_power_10"]
    )
    jk_index = int(jk_index_source["index"]["value"])
    chi_f3 = chi_full * f3_index
    chi_jk = Fraction(
        jk_index_source["subgroup_84_wall_chamber"]
        ["rational_euler_characteristic"]
    )
    chambers = {row["id"]: row for row in shared["chambers"]}

    calibrations = {
        "F3-V4-H10-123": {
            "facets": chambers["F3-V4-H10-123"]["facet_count"],
            "wall_types": chambers["F3-V4-H10-123"]
            ["primitive_mirror_type_histogram"],
            "wall_span_index": chambers["F3-V4-H10-123"]["wall_span"]
            ["index_in_full_Picard_lattice"],
            "finite_vertices": chambers["F3-V4-H10-123"]["literal_Vinberg"]
            ["finite_vertices"],
            "cusps": chambers["F3-V4-H10-123"]["literal_Vinberg"]["cusps"],
            "elliptic_rank9_edges": chambers["F3-V4-H10-123"]
            ["literal_Vinberg"]["elliptic_rank9_edges"],
            "classification": chambers["F3-V4-H10-123"]["literal_Vinberg"]
            ["classification"],
            "index_in_W_S": f3_index,
            "rational_euler_characteristic": str(chi_f3),
        },
        "JK-H10-0007": {
            "facets": chambers["JK-H10-0007"]["facet_count"],
            "wall_types": chambers["JK-H10-0007"]
            ["primitive_mirror_type_histogram"],
            "wall_span_index": chambers["JK-H10-0007"]["wall_span"]
            ["index_in_full_Picard_lattice"],
            "finite_vertices": chambers["JK-H10-0007"]["literal_Vinberg"]
            ["finite_vertices"],
            "cusps": chambers["JK-H10-0007"]["literal_Vinberg"]["cusps"],
            "elliptic_rank9_edges": chambers["JK-H10-0007"]
            ["literal_Vinberg"]["elliptic_rank9_edges"],
            "classification": chambers["JK-H10-0007"]["literal_Vinberg"]
            ["classification"],
            "index_in_W_S": jk_index,
            "rational_euler_characteristic": str(chi_jk),
        },
    }

    observed_filter = {
        (row["square"], row["divisibility"]): (
            row["self_pairing_divisibility_compatible"],
            row["discriminant_form_compatible"],
            row["survives_all_arithmetic_filters"],
        )
        for row in discriminant["candidate_filter"]
    }
    expected_filter = {
        (-2, 1): (True, True, True),
        (-2, 2): (True, False, False),
        (-4, 2): (True, True, True),
        (-2, 4): (False, False, False),
        (-4, 4): (True, True, True),
        (-8, 4): (True, False, False),
    }
    alphabet_keys = {
        f"{row['square']},{row['divisibility']}" for row in alphabet
    }
    checks = {
        "S_even_rank11_signature_1_10": (
            gram.rows == 11
            and exact_inertia(gram) == (1, 10, 0)
            and all(gram[index, index] % 2 == 0 for index in range(11))
        ),
        "S_determinant_1024": gram.det() == 1024,
        "S_discriminant_group_2x8_plus_4": smith
        == [1, 1] + [2] * 8 + [4],
        "full_discriminant_form_pinned": lattice_source["discriminant_form"]
        == "q_plus(F2^8) + <-1/4>",
        "all_1024_discriminant_classes_enumerated": (
            discriminant["representative_count"] == 1024
            and discriminant["distinct_classes_mod_S"] == 1024
            and discriminant["representative_count"] == abs(int(gram.det()))
            and discriminant["all_representatives_lie_in_S_dual"]
            and all(discriminant["block_checks"].values())
            and sum(discriminant["exact_order_q_mod_2_histogram"].values())
            == 1024
        ),
        "exact_order_q_histogram": discriminant[
            "exact_order_q_mod_2_histogram"
        ]
        == {
            "order_1_q_0": 1,
            "order_2_q_0": 255,
            "order_2_q_1": 256,
            "order_4_q_3/4": 240,
            "order_4_q_7/4": 272,
        },
        "self_pairing_then_full_form_filter_naive_upper_bound_exactly": (
            observed_filter == expected_filter
            and discriminant["after_d_divides_self_pairing_filter"]
            == [
                {"square": -2, "divisibility": 1},
                {"square": -2, "divisibility": 2},
                {"square": -4, "divisibility": 2},
                {"square": -4, "divisibility": 4},
                {"square": -8, "divisibility": 4},
            ]
        ),
        "all_three_surviving_types_have_explicit_witnesses": (
            {
                (row["square"], row["divisibility"])
                for row in witnesses.values()
            }
            == {(-2, 1), (-4, 2), (-4, 4)}
            and all(
                row["primitive"]
                and row["dual_class_exact_order"] == row["divisibility"]
                and row["reflection_integral"]
                and row["reflection_preserves_S"]
                and row["hyperbolic_mirror_normal"]
                for row in witnesses.values()
            )
        ),
        "positive_norm_integral_reflection_out_of_scope": (
            positive_boundary["square"] == 2
            and positive_boundary["divisibility"] == 1
            and positive_boundary["primitive"]
            and positive_boundary["reflection_integral"]
            and positive_boundary["reflection_preserves_S"]
            and not positive_boundary["hyperbolic_mirror_normal"]
        ),
        "period_locus_is_nine_dimensional_full_Picard": (
            lattice_source["period_consequence"]["K3_moduli_dimension"] == 9
            and lattice_source["period_consequence"]
            ["very_general_full_Picard_lattice"]
            == "NS(X)=S"
            and all(
                lattice_source["primitive_K3_embedding"]["checks"].values()
            )
        ),
        "full_reflection_chamber_is_complete_13_wall_nonuniform": (
            full["facet_count"] == 13
            and full["finite_vertices"] == 29
            and full["cusps"] == 3
            and full["hyperideal_rays"] == 0
            and full["vinberg_stopping_fired"] is True
            and full["primitive_root_type_histogram"]
            == {
                "square_-2_div_1": 1,
                "square_-4_div_2": 11,
                "square_-4_div_4": 1,
            }
        ),
        "calibrations_are_green_and_use_exact_alphabet": (
            all(all(row["checks"].values()) for row in chambers.values())
            and all(
                set(row["wall_types"]) <= alphabet_keys
                for row in calibrations.values()
            )
        ),
        "subgroup_indices_exact_65802240_and_16450560": [
            calibrations[key]["index_in_W_S"]
            for key in ("F3-V4-H10-123", "JK-H10-0007")
        ]
        == [65_802_240, 16_450_560],
        "subgroup_indices_equal_exact_Euler_ratios": (
            chi_f3 / chi_full == 65_802_240
            and chi_jk / chi_full == 16_450_560
            and chi_f3 == Fraction(-527, 2880)
            and chi_jk == Fraction(-527, 11520)
        ),
        "index_source_pins_same_full_chamber_lineage": (
            jk_index_source["sources"]["maximal_chamber_audit_sha256"]
            == EXPECTED_FILE_SHA256["full_reflection_lineage"]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError([key for key, value in checks.items() if not value])

    certificate: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": {
            "common_period_moduli_criterion": (
                "Every genuine finite-sided oriented H10 Coxeter polyhedron "
                "admitting an orientation-compatible realization by primitive "
                "negative vectors whose reflections preserve "
                "S=U+E8(-2)+<-4> is a weak period-level full-Picard divisor "
                "carrier, after choosing an S-marking of a very general K3 "
                "surface on the nine-dimensional S-polarized period locus."
            ),
            "converse_for_this_marked_period_locus": (
                "Every primitive negative Picard vector defining an integral "
                "hyperbolic mirror reflection has type (-2,div1), (-4,div2), "
                "or (-4,div4). Self-pairing divisibility excludes (-2,div4), "
                "and the full discriminant form excludes (-2,div2) and "
                "(-8,div4)."
            ),
            "finite_volume_subgroup_criterion": (
                "Let W(S) be generated by every primitive negative integral "
                "mirror reflection of S. Its complete 13-wall chamber has "
                "finite volume and three cusps. If a genuine Coxeter "
                "polyhedron P is the fundamental domain of its facet group "
                "Gamma<=W(S), then vol(P) is finite exactly when [W(S):Gamma] "
                "is finite, and every finite-volume such P is noncompact."
            ),
            "compact_no_go_inside_S": (
                "No compact H10 Coxeter polytope admits a primitive-negative "
                "integral-mirror realization in S."
            ),
        },
        "lattice": {
            "name": "S=U+E8(-2)+<-4>",
            "rank": 11,
            "signature": [1, 10, 0],
            "determinant": 1024,
            "smith": smith,
            "discriminant_group": "(Z/2)^8 + Z/4",
            "discriminant_form": lattice_source["discriminant_form"],
            "discriminant_exponent": 4,
            "full_discriminant_form_audit": discriminant,
            "primitive_negative_integral_reflection_alphabet": alphabet,
            "explicit_negative_type_witnesses": witnesses,
            "positive_norm_integral_reflection_boundary_witness": positive_boundary,
            "alphabet_scope": (
                "This is the exhaustive alphabet for primitive negative "
                "hyperbolic-mirror normals in S, not for positive-norm lattice "
                "reflections."
            ),
        },
        "period_carrier": {
            "primitive_K3_embedding": lattice_source["primitive_K3_embedding"]
            ["embedding_description"],
            "orthogonal_complement": lattice_source["primitive_K3_embedding"]
            ["orthogonal_complement"],
            "moduli_dimension": lattice_source["period_consequence"]
            ["K3_moduli_dimension"],
            "very_general_full_Picard": lattice_source["period_consequence"]
            ["very_general_full_Picard_lattice"],
            "level": (
                "period-level: the primitive embedding and K3 period "
                "surjectivity give the S-polarized locus; no explicit common "
                "Weierstrass, quartic, or other equation family is certified"
            ),
            "wall_scope": (
                "An S-marking makes every transported wall an algebraic "
                "divisor class in full NS. Equations, effective representatives, "
                "and irreducible curves for arbitrary walls are not supplied."
            ),
        },
        "calibrations": calibrations,
        "full_reflection_group": {
            "name": "W(S), the full integral reflection group; classical f_2^10 chamber",
            "fundamental_facets": full["facet_count"],
            "finite_vertices": full["finite_vertices"],
            "cusps": full["cusps"],
            "nonuniform": True,
            "primitive_wall_types": full["primitive_root_type_histogram"],
            "rational_euler_characteristic": str(chi_full),
            "proof_of_subgroup_criterion": (
                "For a genuine Coxeter polyhedron, the Poincare polyhedron "
                "theorem makes P a fundamental domain of its facet-reflection "
                "group Gamma. Primitive reduction puts Gamma in W(S). Modulo "
                "measure-zero boundaries, a Gamma domain is a disjoint union "
                "of [W(S):Gamma] W(S)-chambers, so its volume is finite iff "
                "the index is finite. A finite-index subgroup of nonuniform "
                "W(S) contains positive powers of parabolic elements and is "
                "nonuniform; hence its Coxeter domain is not compact."
            ),
        },
        "promotion_boundary": {
            "minus_two_effectivity": (
                "For each chamber separately, a Weyl re-marking can move a "
                "generic strict interior point into the ample chamber; K3 "
                "Riemann--Roch then makes each oriented (-2) facet effective."
            ),
            "reducibility": (
                "An effective (-2) divisor class can be reducible and need not "
                "itself be an irreducible nef-wall curve."
            ),
            "full_nef_criterion": (
                "After chamber-specific Weyl re-marking, full-nef equality "
                "requires every irredundant facet to be a primitive (-2)-root "
                "hyperplane and no other (-2)-root hyperplane to meet the "
                "chamber interior."
            ),
            "minus_four_obstruction": (
                "An irredundant primitive (-4,div2) or (-4,div4) facet cannot "
                "be a full K3 nef wall, although its reflection is integral."
            ),
            "volume_gate": (
                "The period carrier alone supplies no volume conclusion. Use "
                "a complete Vinberg endpoint/parabolic-neighborhood audit or "
                "an independent finite-index proof in W(S)."
            ),
        },
        "scope": {
            "new_abstract_polytopes": "none",
            "new_result": (
                "a bounded common period-moduli carrier criterion, exact wall "
                "alphabet, and finite-index/compactness boundary for S"
            ),
            "equation_family": "not asserted",
            "surface_automorphisms": "not asserted",
            "full_nef_for_calibrations": "not asserted",
            "simultaneous_ample_marking": (
                "not asserted; different chambers may require different Weyl markings"
            ),
            "worldwide_priority": "none",
        },
        "sources": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "file_sha256": EXPECTED_FILE_SHA256[key],
                "payload_sha256": EXPECTED_PAYLOAD_SHA256[key],
            }
            for key, path in SOURCES.items()
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    certificate["payload_sha256"] = _payload_sha256(certificate)
    return certificate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    certificate = build_certificate()
    if args.write:
        OUTPUT.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != certificate:
            raise ArithmeticError("archived common S-carrier certificate changed")
    print(json.dumps(certificate, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
