#!/usr/bin/env python3
"""Certify a common full-Picard K3 carrier of the known H10/34 chamber.

The A15-C2 and characteristic-two V4 fixed sections have the same normalized
Coxeter Gram matrix but different integral fixed lattices.  This certificate
closes their reflection-stable overlattice denominators.  The unique even
index-four A15 closure is the Picard lattice U+D4(-1)+A1(-1)^5.  Independently,
the doubled dual of the unique maximal (odd index-eight) V4 closure is exactly
the same marked lattice and sends all 34 oriented mirrors to the A15 mirrors.

An explicit isometry attaches the chamber to Roulleau's actual 90-wall K3
stock.  Five published curve reflections put the sum-interior in the nef cone;
adding Roulleau's ample class D18 gives an ample class still strictly inside
all 34 imported halfspaces.  This is a new algebraic carrier of a known
polytope, not a new abstract chamber or a full-nef/automorphism theorem.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _overlattice_basis,
    _q_histogram,
    _stable_subgroups,
)
from scripts.certify_rank12_five_chamber_simultaneous_carrier import (
    _binary_discriminant_isometry,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
A15_SOURCE = ROOT / "data/enumerations/alexeev_20_2_1_15_9_1_fixed_section_census.json"
A15_VINBERG = ROOT / "data/enumerations/alexeev_20_2_1_15_9_1_literal_vinberg_census.json"
SS2_SOURCE = ROOT / "data/enumerations/borcherds_char2_v4_h10_34_certificate.json"
K3_SOURCE = ROOT / "data/enumerations/roulleau_r11_ud4a1x5_90_denominator.json"
OUTPUT = ROOT / "data/enumerations/h10_34_ud4a1x5_k3_carrier.json"
SCHEMA = "coxeter.h10_34.ud4a1x5_common_k3_carrier.v1"

EXPECTED_SHA256 = {
    "a15_source": "0294897befa7357f485d45df17d5e3e203656a3605c9c3231ab362d3a30754a6",
    "a15_vinberg": "1ca27c1287a05824fcbede3d13a0351ea9f3bc7889463bb024229faaf847a283",
    "ss2_source": "af88ac64a7011a798f8d6f7dd35b75d79088a00da74f321e7b8044184d80f5dc",
    "k3_source": "795d094c61c07eac28a58081e6360ba0d889ff3d2d9de576cadbbcd673e78a39",
}

# Coordinates in the A15 maximal stable closure go to Roulleau's archived
# curve basis A1,...,A11.  The replay below checks determinant and congruence.
A15_M_TO_ROULLEAU = sp.Matrix(
    [
        [0, 0, 5, 0, 4, 4, 0, 0, 3, -10, -18],
        [2, 1, -2, 0, -2, 0, -2, 0, -2, 2, 4],
        [0, 0, -3, 0, -2, -2, 0, 0, -2, 4, 10],
        [2, 1, -1, 0, -1, 1, -2, 0, -2, 1, 0],
        [0, -1, 2, 0, 1, 1, 0, 0, 1, -2, -4],
        [0, 0, 1, 0, 1, 0, 1, 0, 1, -1, -3],
        [0, 0, 2, 0, 2, 1, 0, 0, 1, -2, -8],
        [-1, 0, 1, 0, 1, 1, 0, 0, 1, -1, -3],
        [-1, 0, 2, 0, 1, 1, 1, 0, 1, -3, -7],
        [0, 0, 2, 1, 2, 2, 0, 0, 2, -5, -10],
        [0, 0, 2, 1, 2, 2, 0, 1, 1, -4, -9],
    ]
)

# Coordinates in N=2*M_SS2^{-1} go to the A15 maximal stable closure.  More
# strongly than lattice isometry, this sends every primitive oriented mirror
# to the exact A15 mirror under the archived 34-facet permutation.
SS2_N_TO_A15_M = sp.Matrix(
    [
        [-50, 37, 66, 44, 22, -76, -197, 319, -254, 114, 188],
        [269, -202, -360, -241, -122, 410, 1057, -1722, 1362, -608, -1002],
        [376, -282, -502, -336, -170, 573, 1478, -2406, 1904, -851, -1402],
        [342, -257, -457, -306, -155, 522, 1347, -2192, 1735, -776, -1278],
        [78, -58, -103, -69, -35, 119, 308, -499, 396, -178, -293],
        [-140, 105, 186, 124, 62, -213, -551, 896, -710, 318, 524],
        [-100, 74, 132, 88, 44, -152, -394, 638, -507, 228, 376],
        [-188, 142, 252, 169, 86, -288, -744, 1210, -958, 428, 706],
        [-130, 98, 174, 116, 59, -199, -514, 836, -662, 296, 488],
        [-36, 27, 48, 32, 16, -55, -142, 231, -183, 82, 135],
        [84, -63, -112, -75, -38, 128, 330, -537, 425, -190, -313],
    ]
)

EXPECTED_WEYL_SEQUENCE_ZERO_BASED = [19, 3, 13, 15, 2]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _payload_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(normal[index, index]))
        for index in range(min(normal.rows, normal.cols))
    ]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def _histogram(counter: Counter[tuple[int, ...]]) -> dict[str, int]:
    return {
        ",".join(str(value) for value in key): count
        for key, count in sorted(counter.items())
    }


def _stable_census(
    identifier: str, gram: sp.Matrix, roots: list[sp.Matrix]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    reflections = [_reflection(gram, root) for root in roots]
    audit = _stable_subgroups(gram, reflections)
    groups = audit.pop("groups_internal")
    left_inverse = audit.pop("left_smith_inverse")
    rows: list[dict[str, Any]] = []
    internal: list[dict[str, Any]] = []
    for position, group in enumerate(groups):
        subgroup = audit["stable_subgroups"][position]
        basis = _overlattice_basis(
            gram, left_inverse, subgroup["cosets_in_smith_coordinates"]
        )
        closure_gram = sp.simplify(basis.T * gram * basis)
        closure_roots: list[sp.Matrix] = []
        root_types: Counter[tuple[int, int]] = Counter()
        actions: list[sp.Matrix] = []
        for root, reflection in zip(roots, reflections):
            raw = basis.inv() * root
            if any(not value.is_Integer for value in raw):
                raise ArithmeticError(f"{identifier} root left stable closure")
            primitive = (raw / _content(raw)).applyfunc(int)
            closure_roots.append(primitive)
            square = int((primitive.T * closure_gram * primitive)[0])
            divisibility = _content(closure_gram * primitive)
            root_types[(square, divisibility)] += 1
            action = sp.simplify(basis.inv() * reflection * basis)
            actions.append(action)
        root_hnf = hermite_normal_form(sp.Matrix.hstack(*closure_roots))

        scaled_duals: dict[str, Any] = {}
        for scale in (2, 4):
            dual = sp.simplify(scale * closure_gram.inv())
            integral = all(value.is_Integer for value in dual)
            record: dict[str, Any] = {
                "scale": scale,
                "operation": f"Gram -> {scale}*Gram^(-1)",
                "integral": integral,
            }
            if integral:
                dual = dual.applyfunc(int)
                dual_roots = []
                dual_checks = []
                dual_types: Counter[tuple[int, int]] = Counter()
                for primitive, action in zip(closure_roots, actions):
                    divisibility = _content(closure_gram * primitive)
                    covector = (closure_gram * primitive / divisibility).applyfunc(int)
                    dual_roots.append(covector)
                    dual_reflection = _reflection(dual, covector)
                    dual_checks.append(
                        dual_reflection == action.T
                        and all(value.is_Integer for value in dual_reflection)
                    )
                    dual_types[
                        (
                            int((covector.T * dual * covector)[0]),
                            _content(dual * covector),
                        )
                    ] += 1
                dual_hnf = hermite_normal_form(sp.Matrix.hstack(*dual_roots))
                record.update(
                    {
                        "gram": _payload(dual),
                        "signature": list(exact_inertia(dual)),
                        "determinant": int(dual.det()),
                        "smith_invariants": _smith(dual),
                        "parity": (
                            "even"
                            if all(int(dual[i, i]) % 2 == 0 for i in range(11))
                            else "odd"
                        ),
                        "content": math.gcd(*(abs(int(value)) for value in dual)),
                        "primitive_mirror_types_square_divisibility": _histogram(
                            dual_types
                        ),
                        "root_span_index": abs(int(dual_hnf.det())),
                        "all_contragredient_reflections_integral": all(dual_checks),
                    }
                )
            scaled_duals[str(scale)] = record

        checks = {
            "index_basis_law": abs(sp.Rational(1, basis.det())) == len(group),
            "determinant_index_square_law": (
                closure_gram.det() == gram.det() / len(group) ** 2
            ),
            "all_reflections_integral": all(
                value.is_Integer for action in actions for value in action
            ),
            "all_reflections_are_isometries": all(
                action.T * closure_gram * action == closure_gram
                for action in actions
            ),
            "all_34_roots_span_full_rank": root_hnf.rows == root_hnf.cols == 11,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"stable closure failed for {identifier}/{position}")
        row = {
            "stable_subgroup_position_zero_based": position,
            "index_over_fixed_lattice": len(group),
            "basis_columns_in_fixed_lattice": _payload(basis),
            "gram": _payload(closure_gram),
            "signature": list(exact_inertia(closure_gram)),
            "determinant": int(closure_gram.det()),
            "smith_invariants": _smith(closure_gram),
            "parity": (
                "even"
                if all(int(closure_gram[i, i]) % 2 == 0 for i in range(11))
                else "odd"
            ),
            "primitive_mirror_types_square_divisibility": _histogram(root_types),
            "root_span_index": abs(int(root_hnf.det())),
            "scaled_duals": scaled_duals,
            "checks": checks,
        }
        rows.append(row)
        internal.append(
            {
                "basis": basis,
                "gram": closure_gram,
                "roots": closure_roots,
                "actions": actions,
            }
        )
    return {
        "id": identifier,
        "discriminant_module_audit": audit,
        "stable_integral_subgroup_orders": audit["stable_integral_subgroup_orders"],
        "stable_even_subgroup_orders": audit["stable_even_subgroup_orders"],
        "maximal_positions_zero_based": audit[
            "maximal_stable_subgroup_positions_zero_based"
        ],
        "closures": rows,
    }, internal


def _k3_glue(picard: sp.Matrix) -> dict[str, Any]:
    """Give an exact primitive embedding in the even unimodular K3 lattice."""

    u = sp.Matrix([[0, 1], [1, 0]])
    u2 = sp.Matrix([[0, 2], [2, 0]])
    d4 = -sp.Matrix(
        [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]
    )
    complement = sp.diag(u, u2, d4, sp.diag(-2, -2, -2))
    anti = _binary_discriminant_isometry(picard, -complement)
    if anti is None:
        raise ArithmeticError("Picard/complement discriminant anti-isometry failed")

    mod_n, left_n, _ = _discriminant_data(picard)
    mod_t, left_t, _ = _discriminant_data(-complement)
    idx_n = [i for i, modulus in enumerate(mod_n) if modulus > 1]
    idx_t = [i for i, modulus in enumerate(mod_t) if modulus > 1]
    generators = [sp.eye(22).col(i) for i in range(22)]
    glue_vectors = []
    for source_index, image in zip(idx_n, anti["basis_images_as_columns"]):
        snf_n = sp.zeros(11, 1)
        snf_n[source_index] = 1
        dual_n = picard.inv() * left_n * snf_n
        snf_t = sp.zeros(11, 1)
        for bit, target_index in zip(image, idx_t):
            snf_t[target_index] = bit
        dual_t = (-complement).inv() * left_t * snf_t
        vector = dual_n.col_join(dual_t)
        glue_vectors.append(vector)
        generators.append(vector)
    denominator = math.lcm(
        *(int(value.q) for vector in generators for value in vector)
    )
    integer_generators = sp.Matrix.hstack(
        *(denominator * vector for vector in generators)
    ).applyfunc(int)
    glue_basis = hermite_normal_form(integer_generators) / denominator
    split = sp.diag(picard, complement)
    k3_gram = sp.simplify(glue_basis.T * split * glue_basis)
    n_embedding = sp.simplify(glue_basis.inv() * sp.eye(22)[:, :11])
    t_embedding = sp.simplify(glue_basis.inv() * sp.eye(22)[:, 11:])
    checks = {
        "all_128_discriminant_classes_checked": anti["all_elements_checked"] == 128,
        "anti_isometry_bijective": anti["bijective"],
        "seven_order_two_glue_generators": len(glue_vectors) == 7,
        "glue_index_128": abs(glue_basis.det()) == sp.Rational(1, 128),
        "even_unimodular_signature_3_19": (
            all(value.is_Integer for value in k3_gram)
            and all(int(k3_gram[i, i]) % 2 == 0 for i in range(22))
            and k3_gram.det() == -1
            and exact_inertia(k3_gram) == (3, 19, 0)
        ),
        "Picard_embedding_integral_isometric_primitive": (
            all(value.is_Integer for value in n_embedding)
            and n_embedding.T * k3_gram * n_embedding == picard
            and _smith(n_embedding) == [1] * 11
        ),
        "complement_embedding_integral_isometric_primitive": (
            all(value.is_Integer for value in t_embedding)
            and t_embedding.T * k3_gram * t_embedding == complement
            and _smith(t_embedding) == [1] * 11
        ),
        "orthogonal_images": (
            n_embedding.T * k3_gram * t_embedding == sp.zeros(11, 11)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"primitive K3 glue failed: {checks}")
    return {
        "complement_lattice": "U+U(2)+D4(-1)+A1(-1)^3",
        "complement_gram": _payload(complement),
        "discriminant_anti_isometry": anti,
        "glue_vectors": [
            [int(x) if x.q == 1 else str(x) for x in vector]
            for vector in glue_vectors
        ],
        "glued_even_unimodular_gram": _payload(k3_gram),
        "Picard_embedding_columns": _payload(n_embedding),
        "complement_embedding_columns": _payload(t_embedding),
        "classification": "the unique even unimodular lattice of signature (3,19)",
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "a15_source": A15_SOURCE,
        "a15_vinberg": A15_VINBERG,
        "ss2_source": SS2_SOURCE,
        "k3_source": K3_SOURCE,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")
    a15_source = json.loads(A15_SOURCE.read_text())
    a15_vinberg = json.loads(A15_VINBERG.read_text())
    ss2_source = json.loads(SS2_SOURCE.read_text())
    k3_source = json.loads(K3_SOURCE.read_text())

    a15 = next(row for row in a15_source["cases"] if row["id"] == "A15_9_1_C01")
    a15_v = next(row for row in a15_vinberg["cases"] if row["id"] == "A15_9_1_C01")
    a15_records = [
        row
        for row in a15["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    ]
    gram_a15_l = sp.Matrix(a15["saturated_fixed_lattice"]["gram"])
    roots_a15_l = [
        sp.Matrix(row["primitive_integral_same_direction_coordinates"])
        for row in a15_records
    ]
    gram_ss2_l = sp.Matrix(ss2_source["fixed_lattice"]["gram"])
    ss2_records = ss2_source["restricted_chamber"][
        "primitive_facet_lattice_records"
    ]
    roots_ss2_l = [sp.Matrix(row["primitive_fixed_coordinates"]) for row in ss2_records]

    a15_census, a15_internal = _stable_census(
        "A15-C2-H10-34", gram_a15_l, roots_a15_l
    )
    ss2_census, ss2_internal = _stable_census(
        "SS2-V4-H10-34", gram_ss2_l, roots_ss2_l
    )
    denominator_checks = {
        "A15_exact_stable_orders": (
            a15_census["stable_integral_subgroup_orders"] == [1, 4]
            and a15_census["stable_even_subgroup_orders"] == [1, 4]
            and a15_census["maximal_positions_zero_based"] == [1]
        ),
        "SS2_exact_stable_orders": (
            ss2_census["stable_integral_subgroup_orders"] == [1, 2, 4, 8]
            and ss2_census["stable_even_subgroup_orders"] == [1, 4]
            and ss2_census["maximal_positions_zero_based"] == [3]
        ),
    }
    if not all(denominator_checks.values()):
        raise ArithmeticError("stable denominator changed")

    a15_m_data = a15_internal[1]
    gram_a15_m = a15_m_data["gram"]
    basis_a15_m = a15_m_data["basis"]
    roots_a15_m = []
    for record, root_l in zip(a15_records, roots_a15_l):
        raw = basis_a15_m.inv() * root_l
        orientation = int(sp.Rational(record["presentation_scale_from_primitive"]))
        roots_a15_m.append(orientation * (raw / _content(raw)).applyfunc(int))

    ss2_m_data = ss2_internal[3]
    gram_ss2_m = ss2_m_data["gram"]
    basis_ss2_m = ss2_m_data["basis"]
    gram_ss2_n = sp.simplify(2 * gram_ss2_m.inv()).applyfunc(int)
    roots_ss2_m = []
    roots_ss2_n = []
    for root_l in roots_ss2_l:
        raw = basis_ss2_m.inv() * root_l
        root_m = (raw / _content(raw)).applyfunc(int)
        divisibility = _content(gram_ss2_m * root_m)
        root_n = (gram_ss2_m * root_m / divisibility).applyfunc(int)
        roots_ss2_m.append(root_m)
        roots_ss2_n.append(root_n)

    permutation = [
        value - 1
        for value in ss2_source["exact_deduplication"][
            "source_to_target_facet_permutation_one_based"
        ]
    ]
    common_bridge_checks = {
        "A15_closure_even_det128": (
            gram_a15_m.det() == 128
            and exact_inertia(gram_a15_m) == (1, 10, 0)
            and _smith(gram_a15_m) == [1] * 4 + [2] * 7
            and all(int(gram_a15_m[i, i]) % 2 == 0 for i in range(11))
        ),
        "SS2_maximal_closure_odd_det16": (
            gram_ss2_m.det() == 16
            and exact_inertia(gram_ss2_m) == (1, 10, 0)
            and _smith(gram_ss2_m) == [1] * 7 + [2] * 4
            and any(int(gram_ss2_m[i, i]) % 2 for i in range(11))
        ),
        "SS2_doubled_dual_even_det128": (
            gram_ss2_n.det() == 128
            and exact_inertia(gram_ss2_n) == (1, 10, 0)
            and _smith(gram_ss2_n) == [1] * 4 + [2] * 7
            and all(int(gram_ss2_n[i, i]) % 2 == 0 for i in range(11))
        ),
        "bridge_unimodular_isometry": (
            SS2_N_TO_A15_M.det() == 1
            and SS2_N_TO_A15_M.T * gram_a15_m * SS2_N_TO_A15_M
            == gram_ss2_n
        ),
        "all_34_oriented_roots_identified": all(
            SS2_N_TO_A15_M * roots_ss2_n[index]
            == roots_a15_m[permutation[index]]
            for index in range(34)
        ),
        "exact_discriminant_form": (
            _q_histogram(gram_a15_m)
            == _q_histogram(gram_ss2_n)
            == {"0": 36, "1": 28, "1/2": 28, "3/2": 36}
        ),
    }
    if not all(common_bridge_checks.values()):
        raise ArithmeticError(f"common bridge failed: {common_bridge_checks}")

    gram_k3 = sp.Matrix(k3_source["lattice"]["basis_gram_A1_through_A11"])
    walls_k3 = [
        sp.Matrix(row)
        for row in k3_source["walls"]["coordinate_rows_in_A1_through_A11"]
    ]
    target_checks = {
        "actual_Roulleau_full_Picard_stock": (
            k3_source["all_checks_pass"]
            and k3_source["lattice"]["name"] == "U+D4(-1)+A1(-1)^5"
            and k3_source["hyperbolic_geometry"]["status"]
            == "published_finite_volume_noncompact_full_nef_chamber"
            and len(walls_k3) == 90
        ),
        "A15_closure_to_K3_unimodular_isometry": (
            A15_M_TO_ROULLEAU.det() == -1
            and A15_M_TO_ROULLEAU.T * gram_k3 * A15_M_TO_ROULLEAU
            == gram_a15_m
        ),
    }
    if not all(target_checks.values()):
        raise ArithmeticError(f"Roulleau target failed: {target_checks}")

    roots_k3_unframed = [A15_M_TO_ROULLEAU * root for root in roots_a15_m]
    interior_unframed = sum(roots_k3_unframed, sp.zeros(11, 1))
    weyl = sp.eye(11)
    interior = interior_unframed
    sequence = []
    for _ in range(1000):
        pairings = [int((interior.T * gram_k3 * wall)[0]) for wall in walls_k3]
        negative = [index for index, value in enumerate(pairings) if value < 0]
        if not negative:
            break
        chosen = min(negative, key=lambda index: (pairings[index], index))
        reflection = _reflection(gram_k3, walls_k3[chosen])
        interior = reflection * interior
        weyl = reflection * weyl
        sequence.append(chosen)
    else:
        raise ArithmeticError("Weyl reduction did not terminate")
    roots_k3 = [weyl * root for root in roots_k3_unframed]
    d18 = sp.Matrix([1, 5, 2, 4, 0, 0, 0, 0, 0, 0, 0])
    ample = interior + d18
    nef_pairings = [int((interior.T * gram_k3 * wall)[0]) for wall in walls_k3]
    ample_nef_pairings = [int((ample.T * gram_k3 * wall)[0]) for wall in walls_k3]
    imported_pairings = [int((interior.T * gram_k3 * root)[0]) for root in roots_k3]
    ample_imported_pairings = [int((ample.T * gram_k3 * root)[0]) for root in roots_k3]
    root_matrix = sp.Matrix.hstack(*roots_k3).applyfunc(int)
    root_hnf = hermite_normal_form(root_matrix)
    target_types: Counter[tuple[int, int]] = Counter()
    for root in roots_k3:
        target_types[
            (int((root.T * gram_k3 * root)[0]), _content(gram_k3 * root))
        ] += 1
    ample_checks = {
        "deterministic_five_step_Weyl_reduction": sequence
        == EXPECTED_WEYL_SEQUENCE_ZERO_BASED,
        "Weyl_action_integral_isometry": (
            all(value.is_Integer for value in weyl)
            and abs(weyl.det()) == 1
            and weyl.T * gram_k3 * weyl == gram_k3
        ),
        "roots_generate_full_Picard_lattice": root_hnf == sp.eye(11),
        "target_root_types": target_types
        == Counter({(-2, 2): 16, (-4, 2): 10, (-2, 1): 8}),
        "sum_interior_big_nef": (
            interior == sum(roots_k3, sp.zeros(11, 1))
            and int((interior.T * gram_k3 * interior)[0]) == 1648
            and min(nef_pairings) == 0
            and sum(value == 0 for value in nef_pairings) == 5
            and min(imported_pairings) == 4
        ),
        "D18_published_ample_class": (
            int((d18.T * gram_k3 * d18)[0]) == 18
            and min(int((d18.T * gram_k3 * wall)[0]) for wall in walls_k3) == 1
        ),
        "explicit_ample_in_imported_chamber": (
            int((ample.T * gram_k3 * ample)[0]) == 2178
            and min(ample_nef_pairings) == 1
            and min(ample_imported_pairings) == 4
        ),
    }
    if not all(ample_checks.values()):
        raise ArithmeticError(f"ample framing failed: {ample_checks}")

    inverse_permutation = [0] * 34
    for ss2_index, a15_index in enumerate(permutation):
        inverse_permutation[a15_index] = ss2_index
    facet_records = []
    for a15_index, root in enumerate(roots_k3):
        ss2_index = inverse_permutation[a15_index]
        reflection = _reflection(gram_k3, root)
        checks = {
            "primitive": _content(root) == 1,
            "exact_A15_transport": (
                root == weyl * A15_M_TO_ROULLEAU * roots_a15_m[a15_index]
            ),
            "exact_SS2_transport": (
                root
                == weyl
                * A15_M_TO_ROULLEAU
                * SS2_N_TO_A15_M
                * roots_ss2_n[ss2_index]
            ),
            "reflection_integral_unimodular_involution_isometry": (
                all(value.is_Integer for value in reflection)
                and reflection.det() == -1
                and reflection**2 == sp.eye(11)
                and reflection.T * gram_k3 * reflection == gram_k3
            ),
            "ample_positive_orientation": int((ample.T * gram_k3 * root)[0]) > 0,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"facet {a15_index} failed: {checks}")
        facet_records.append(
            {
                "facet_zero_based_in_A15_order": a15_index,
                "A15_source_orbit_index_one_based": a15_records[a15_index][
                    "orbit_index"
                ],
                "SS2_facet_zero_based": ss2_index,
                "SS2_source_orbit_index_one_based": ss2_records[ss2_index][
                    "orbit_index_one_based"
                ],
                "primitive_oriented_K3_coordinates": [int(value) for value in root],
                "square": int((root.T * gram_k3 * root)[0]),
                "divisibility": _content(gram_k3 * root),
                "big_nef_interior_pairing": imported_pairings[a15_index],
                "ample_pairing": ample_imported_pairings[a15_index],
                "reflection_matrix": _payload(reflection),
                "checks": checks,
            }
        )

    k3_glue = _k3_glue(gram_k3)
    vinberg_checks = {
        "A15_literal_audit_green": all(a15_v["checks"].values()),
        "SS2_literal_audit_green": ss2_source["vinberg"]["all_checks_pass"],
        "same_counts": (
            a15_v["literal_endpoint_audit"]["finite_vertex_count"]
            == ss2_source["vinberg"]["finite_vertex_count"]
            == 1774
            and a15_v["literal_endpoint_audit"]["ideal_vertex_count"]
            == ss2_source["vinberg"]["ideal_vertex_count"]
            == 146
            and a15_v["literal_endpoint_audit"]["edge_count"]
            == ss2_source["vinberg"]["elliptic_edge_count"]
            == 12792
            and a15_v["diagram_enumeration"]["connected_affine_component_count"]
            == ss2_source["vinberg"]["connected_affine_component_count"]
            == 426
        ),
        "zero_failures": (
            a15_v["literal_endpoint_audit"]["bad_edge_count"] == 0
            and a15_v["literal_endpoint_audit"][
                "bad_parabolic_neighborhood_count"
            ]
            == 0
            and ss2_source["vinberg"]["bad_edge_count"] == 0
            and ss2_source["vinberg"]["bad_parabolic_neighborhood_count"] == 0
        ),
        "finite_volume_noncompact": (
            a15_v["decision"]["finite_volume"]
            and not a15_v["decision"]["compact"]
            and ss2_source["vinberg"]["finite_volume"]
            and not ss2_source["vinberg"]["compact"]
        ),
    }
    checks = {
        "upstream_archives": (
            all(all(case["checks"].values()) for case in a15_source["cases"])
            and ss2_source["all_checks_pass"]
            and k3_source["all_checks_pass"]
        ),
        "stable_denominators": all(denominator_checks.values()),
        "common_integral_bridge": all(common_bridge_checks.values()),
        "actual_K3_target": all(target_checks.values()),
        "primitive_K3_embedding": all(k3_glue["checks"].values()),
        "ample_interior": all(ample_checks.values()),
        "all_34_facets": all(all(row["checks"].values()) for row in facet_records),
        "literal_Vinberg": all(vinberg_checks.values()),
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "id": "H10-34-UD4A1X5-COMMON-K3-CARRIER",
        "upstream": {
            key: {"path": str(path.relative_to(ROOT)), "sha256": EXPECTED_SHA256[key]}
            for key, path in paths.items()
        },
        "theorem": (
            "The known finite-volume noncompact H10/34 Coxeter polytope shared "
            "by A15-C2 and SS2-V4 has a common full-Picard complex-K3 carrier. "
            "The unique even index-four A15 stable closure is Roulleau's "
            "U+D4(-1)+A1(-1)^5 Picard lattice; the doubled dual of the unique "
            "maximal odd index-eight SS2 closure is exactly the same marked "
            "lattice. All 34 oriented mirrors/reflections transport and generate "
            "Picard, and an explicit square-2178 ample divisor lies strictly "
            "inside their imported chamber."
        ),
        "stable_closure_and_scaled_dual_census": {
            "method": (
                "all totally integral discriminant-bilinear subgroups invariant "
                "under all 34 facet reflections, closed under sums; scales 2 and "
                "4 tested on every resulting closure"
            ),
            "A15": a15_census,
            "SS2": ss2_census,
            "checks": denominator_checks,
        },
        "common_lattice_bridge": {
            "A15_maximal_stable_closure_basis": _payload(basis_a15_m),
            "A15_maximal_stable_closure_gram": _payload(gram_a15_m),
            "SS2_maximal_stable_closure_basis": _payload(basis_ss2_m),
            "SS2_maximal_stable_closure_gram": _payload(gram_ss2_m),
            "SS2_doubled_dual_gram": _payload(gram_ss2_n),
            "SS2_doubled_dual_to_A15_closure": _payload(SS2_N_TO_A15_M),
            "SS2_to_A15_facet_permutation_zero_based": permutation,
            "discriminant_q_mod_2_histogram": _q_histogram(gram_a15_m),
            "checks": common_bridge_checks,
        },
        "actual_K3_carrier": {
            "surface_stock": "Roulleau U+D4(-1)+A1(-1)^5 K3 family",
            "Picard_gram_in_A1_through_A11": _payload(gram_k3),
            "published_complete_nef_wall_count": 90,
            "A15_closure_to_Roulleau_basis": _payload(A15_M_TO_ROULLEAU),
            "SS2_doubled_dual_to_Roulleau_basis": _payload(
                A15_M_TO_ROULLEAU * SS2_N_TO_A15_M
            ),
            "primitive_K3_lattice_embedding": k3_glue,
            "period_consequence": (
                "The archived Roulleau stock already realizes this as full NS; "
                "independently, the displayed primitive K3-lattice glue and period "
                "surjectivity give a 9-dimensional lattice-polarized K3 family."
            ),
            "checks": target_checks,
        },
        "ample_framing": {
            "deterministic_rule": "reflect in the most negative published nef wall, ties by wall index",
            "Weyl_sequence_zero_based": sequence,
            "Weyl_sequence_one_based": [value + 1 for value in sequence],
            "Weyl_matrix": _payload(weyl),
            "big_nef_sum_interior_coordinates": [int(value) for value in interior],
            "big_nef_sum_interior_square": int((interior.T * gram_k3 * interior)[0]),
            "big_nef_minimum_90_wall_pairing": min(nef_pairings),
            "big_nef_minimum_imported_pairing": min(imported_pairings),
            "published_D18_coordinates": [int(value) for value in d18],
            "ample_coordinates": [int(value) for value in ample],
            "ample_square": int((ample.T * gram_k3 * ample)[0]),
            "ample_minimum_90_wall_pairing": min(ample_nef_pairings),
            "ample_minimum_imported_pairing": min(ample_imported_pairings),
            "target_mirror_types_square_divisibility": _histogram(target_types),
            "root_matrix_HNF": _payload(root_hnf),
            "generation_index": abs(int(root_hnf.det())),
            "checks": ample_checks,
        },
        "facet_records": facet_records,
        "literal_Vinberg_audit": {
            "criterion": a15_v["method"]["criterion"],
            "finite_vertices": 1774,
            "ideal_vertices": 146,
            "elliptic_rank9_edges": 12792,
            "endpoint_histogram": {"finite-finite": 6046, "finite-ideal": 5648, "ideal-ideal": 1098},
            "connected_affine_components": 426,
            "bad_edges": 0,
            "bad_parabolic_neighborhoods": 0,
            "classification": "finite_volume_noncompact_Coxeter_polytope",
            "checks": vinberg_checks,
        },
        "lineage_and_scope": {
            "abstract_polytope": "known A15-C2-H10-34 = SS2-V4-H10-34",
            "new_positive_gain": (
                "a common full-Picard algebraic K3 realization, exact equation-free "
                "lattice marking, and an explicit ample interior polarization"
            ),
            "not_a_new_abstract_polytope": True,
            "not_the_complete_90_wall_nef_chamber": True,
            "ten_square_minus_four_mirrors_are_not_curve_walls": True,
            "twenty_four_effective_square_minus_two_classes_not_asserted_irreducible": True,
            "reflections_not_asserted_surface_automorphisms": True,
            "worldwide_priority_claim": "none",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if not payload["all_checks_pass"]:
        raise ArithmeticError("H10/34 K3 carrier certificate failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise ArithmeticError("archived H10/34 K3 carrier payload changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()
