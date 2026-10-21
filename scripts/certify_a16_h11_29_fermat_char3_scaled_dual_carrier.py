#!/usr/bin/env python3
"""Certify a characteristic-three Fermat scaled-dual carrier of H11/29.

The A16_6 H11/29 chamber preserves a unique maximal integral odd
reflection-stable overclosure M of determinant -4.  Its doubled dual
N=M^vee(2) is even of determinant -1024.  A fully explicit unimodular
isometry identifies N with the C2-fixed divisor lattice on the supersingular
Fermat quartic in characteristic three,

    U + E8(-2) + <-2>^2.

All 29 oriented halfspaces and reflections are transported.  This is an
additional algebraic divisor-lattice carrier of the already known chamber,
not a new abstract polytope, nef chamber, or automorphism-fixed section.
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

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_118_dual_h11_29_construction.json"
)
STABLE = (
    ROOT / "data/enumerations/alexeev_118_dual_h11_29_stable_overlattices.json"
)
VINBERG = ROOT / "data/enumerations/alexeev_118_dual_h11_29_vinberg.json"
FERMAT_FIXED = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_survivor_fixed_lattices.json"
)
FERMAT_LINEAGE = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_lineage.json"
)
PRIOR_CARRIER = (
    ROOT / "data/enumerations/alexeev_h11_29_ud8_divisor_overlattice_carrier.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/a16_h11_29_fermat_char3_scaled_dual_carrier.json"
)

SCHEMA = "coxeter.a16_h11_29.fermat_char3_scaled_dual_carrier.v1"
EXPECTED_SHA256 = {
    "construction": "53432322ecc6210e047c9b6e9ca101a2b2f86cd5cfef758532f0975bd53f9a99",
    "stable": "3083308ff2b269b8e50b5466b0433dde70c41ef8df089c62dff486b42d14b45b",
    "vinberg": "f711b1869f9106ae0c4c6d48673ff8b3de84495dc8ac72a03787e0434e9497d9",
    "fermat_fixed": "431ddb789fcb30e843efe870a67f9ce276ccd75c4fb65ea4f01f51afdc836124",
    "fermat_lineage": "e1b20c853d3f2f254900c0dd46cdf187093292fb4356e1d6e9de23ca839faf89",
    "prior_carrier": "5d4c2299ecf2b8b853063c142472c46c61f8554197168b06d144924cdc6d07fd",
}


# Columns split a primitive U from U+<-2>^10.  The first two columns are
# u=(3,4,2,1^8,0) and v=(4,3,2,1^8,0).  This is the exceptional isotropic
# orbit: its positive complement has only four norm-two vectors.
EXCEPTIONAL_U_SPLIT_BASIS = sp.Matrix(
    [
        [3, 4, 4, -2, -2, -2, -2, -2, -2, -2, -2, 0],
        [4, 3, 4, -2, -2, -2, -2, -2, -2, -2, -2, 0],
        [2, 2, 7, -4, -4, -4, -4, -4, -4, -4, -4, 0],
        [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    ]
)


# PARI/GP qfisom output.  If K is the positive complement extracted by the
# split above and T=E8(2)+<2>^2, then X^T K X=T.  The permanent certificate
# checks this exact identity and does not require PARI.
POSITIVE_COMPLEMENT_ISOMETRY = sp.Matrix(
    [
        [-3, 0, 0, 0, 0, 0, 3, 0, 0, 5],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, -1, 0, 1],
        [0, -1, 0, 0, 0, 0, 1, 0, 0, 1],
        [-1, 0, 0, 0, 0, 1, 0, 0, 0, 1],
        [-1, 0, 0, 0, 1, -1, 0, 0, 0, 1],
        [-1, 0, 0, 1, -1, 0, 1, 0, 0, 1],
        [-1, 0, 1, -1, 0, 0, 1, 0, 0, 1],
        [-1, 1, -1, 0, 0, 0, 1, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    ]
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(payload: list[list[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[index, index])) for index in range(matrix.rows)]


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


def _standard_target_gram() -> sp.Matrix:
    e8 = 2 * sp.eye(8)
    for left, right in (
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (4, 6),
        (6, 7),
    ):
        e8[left, right] = e8[right, left] = -1
    return sp.diag(
        sp.Matrix([[0, 1], [1, 0]]),
        -2 * e8,
        -2,
        -2,
    )


def build_certificate() -> dict[str, Any]:
    paths = {
        "construction": CONSTRUCTION,
        "stable": STABLE,
        "vinberg": VINBERG,
        "fermat_fixed": FERMAT_FIXED,
        "fermat_lineage": FERMAT_LINEAGE,
        "prior_carrier": PRIOR_CARRIER,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    stable = json.loads(STABLE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    fixed = json.loads(FERMAT_FIXED.read_text(encoding="utf-8"))
    fermat_lineage = json.loads(FERMAT_LINEAGE.read_text(encoding="utf-8"))
    prior_carrier = json.loads(PRIOR_CARRIER.read_text(encoding="utf-8"))

    case = next(
        row for row in construction["cases"] if row["id"] == "A16_6_1_C01"
    )
    stable_case = next(
        row for row in stable["cases"] if row["id"] == "A16_6_1_C01"
    )
    vinberg_case = next(
        row for row in vinberg["audits"] if row["id"] == "A16_6_1_C01"
    )
    fermat_case = next(
        row for row in fixed["cases"] if row["id"] == "char3-C2-H11"
    )

    gram_l = _matrix(case["saturated_fixed_lattice"]["gram"])
    closure = stable_case["unique_maximal_reflection_stable_overlattice"]
    closure_basis = _matrix(closure["basis_columns_in_original_fixed_lattice"])
    gram_m = sp.simplify(closure_basis.T * gram_l * closure_basis)
    gram_n = sp.simplify(2 * gram_m.inv())
    gram_fermat = sp.Matrix(fermat_case["gram"])
    standard = _standard_target_gram()

    # First identify N with U+<-2>^10 through the archived
    # Scharlau--Kirschmer normal form of M.
    sk_gram = sp.diag(
        sp.eye(10),
        sp.Matrix([[0, -2], [-2, -4]]),
    )
    p6 = _matrix(
        stable["maximal_reflection_lineage"]["P6_to_Gram11_entry5"]
    )
    dual_sk_gram = -2 * sk_gram.inv()
    reorder_u = sp.zeros(12)
    reorder_u[11, 0] = 1
    reorder_u[10, 1] = 1
    reorder_u[11, 1] = 1
    for index in range(10):
        reorder_u[index, index + 2] = 1
    a1_model = sp.diag(sp.Matrix([[0, 1], [1, 0]]), *([-2] * 10))
    a1_basis_in_n = p6.inv().T * reorder_u

    # Then use the exceptional U orbit to turn the A1^10 complement into
    # E8(2)+A1(2)^2.
    split = EXCEPTIONAL_U_SPLIT_BASIS
    positive_complement = -(split.T * a1_model * split)[2:, 2:]
    target_positive_complement = -standard[2:, 2:]
    complement_isometry = POSITIVE_COMPLEMENT_ISOMETRY
    standard_basis_in_n = (
        a1_basis_in_n * split * sp.diag(sp.eye(2), complement_isometry)
    )

    # The archived Fermat map has columns equal to the standard target basis
    # in the actual saturated C2-fixed divisor basis.
    standard_basis_in_fermat = _matrix(
        fermat_lineage["fixed_lattice_identification"]["explicit_isometry_P"]
    )
    transport = sp.simplify(
        standard_basis_in_fermat * standard_basis_in_n.inv()
    )

    prior_overclosure = prior_carrier["lattice_and_embedding"]["overclosure_N"]
    lattice_checks = {
        "all_six_upstream_certificates_green": (
            construction["all_checks_pass"]
            and stable["all_checks_pass"]
            and vinberg["all_checks_pass"]
            and fixed["all_checks_pass"]
            and fermat_lineage["all_checks_pass"]
            and prior_carrier["all_checks_pass"]
        ),
        "source_fixed_lattice_rank12_det_minus256_even": (
            exact_inertia(gram_l) == (1, 11, 0)
            and gram_l.det() == -256
            and _smith(gram_l) == [1] * 6 + [2] * 4 + [4] * 2
            and all(int(gram_l[index, index]) % 2 == 0 for index in range(12))
        ),
        "unique_stable_odd_index_eight_closure": (
            closure["index_over_original_lattice"] == 8
            and abs(closure_basis.det()) == sp.Rational(1, 8)
            and gram_m == _matrix(closure["gram"])
            and gram_m.det() == -4
            and exact_inertia(gram_m) == (1, 11, 0)
            and _smith(gram_m) == [1] * 10 + [2] * 2
            and any(int(gram_m[index, index]) % 2 for index in range(12))
            and all(closure["checks"].values())
        ),
        "doubled_dual_integral_even_det_minus1024": (
            all(value.is_Integer for value in gram_n)
            and gram_n.det() == -1024
            and exact_inertia(gram_n) == (1, 11, 0)
            and _smith(gram_n) == [1] * 2 + [2] * 10
            and all(int(gram_n[index, index]) % 2 == 0 for index in range(12))
        ),
        "archived_SK_normal_form_and_dual_A1_model": (
            abs(p6.det()) == 1
            and p6.T * (-gram_m) * p6 == sk_gram
            and reorder_u.T * dual_sk_gram * reorder_u == a1_model
            and abs(a1_basis_in_n.det()) == 1
            and a1_basis_in_n.T * gram_n * a1_basis_in_n == a1_model
        ),
        "exceptional_primitive_U_split": (
            abs(split.det()) == 1
            and split[:, :2].T * a1_model * split[:, :2]
            == sp.Matrix([[0, 1], [1, 0]])
            and split[:, :2].T * a1_model * split[:, 2:] == sp.zeros(2, 10)
            and split.T * a1_model * split
            == sp.diag(sp.Matrix([[0, 1], [1, 0]]), -positive_complement)
        ),
        "positive_complement_is_E8_2_plus_A1_2_squared": (
            abs(complement_isometry.det()) == 1
            and complement_isometry.T
            * positive_complement
            * complement_isometry
            == target_positive_complement
        ),
        "standard_basis_in_N_is_unimodular_isometry": (
            abs(standard_basis_in_n.det()) == 1
            and standard_basis_in_n.T * gram_n * standard_basis_in_n
            == standard
        ),
        "target_is_actual_saturated_Fermat_fixed_divisor_lattice": (
            fermat_case["rank"] == 12
            and fermat_case["signature"] == [1, 11]
            and fermat_case["determinant_absolute"] == 1024
            and fermat_case["smith_diagonal_absolute"] == [1] * 2 + [2] * 10
            and fermat_case["even"] is True
            and fermat_case["smith_kernel"]["saturated_smith_kernel"] is True
            and fermat_case["source_metadata"]["strict_facet_count_before_volume_audit"]
            == 206
            and fermat_lineage["fixed_lattice_identification"]["target_lattice"]
            == "U + E8(-2) + A1(-2)^2"
        ),
        "archived_standard_to_Fermat_isometry": (
            abs(standard_basis_in_fermat.det()) == 1
            and standard_basis_in_fermat.T
            * gram_fermat
            * standard_basis_in_fermat
            == standard
        ),
        "explicit_N_to_Fermat_unimodular_isometry": (
            all(value.is_Integer for value in transport)
            and transport.det() == 1
            and transport.T * gram_fermat * transport == gram_n
        ),
        "arithmetically_different_from_prior_complex_K3_carrier": (
            prior_overclosure["isometry_class"] == "U+D8(-1)+A1(-1)^2"
            and prior_overclosure["determinant"] == -16
            and gram_n.det() == -1024
        ),
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError(
            "A16 Fermat scaled-dual bridge failed: "
            f"{[key for key, value in lattice_checks.items() if not value]}"
        )

    active_records = [
        row
        for row in case["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    ]
    if len(active_records) != 29:
        raise ArithmeticError("active A16_6 facet list changed")

    inverse_closure_basis = closure_basis.inv()
    records: list[dict[str, Any]] = []
    roots_m: list[sp.Matrix] = []
    roots_n: list[sp.Matrix] = []
    roots_fermat: list[sp.Matrix] = []
    closure_types: Counter[tuple[int, int]] = Counter()
    target_types: Counter[tuple[int, int]] = Counter()
    source_to_target_types: Counter[tuple[int, ...]] = Counter()

    for index, source in enumerate(active_records):
        root_l = sp.Matrix(source["primitive_integral_same_direction_coordinates"])
        orientation = int(sp.Rational(source["presentation_scale_from_primitive"]))
        if orientation not in (-1, 1):
            raise ArithmeticError("unexpected non-unit A16_6 presentation scale")
        raw_root_m = inverse_closure_basis * root_l
        if any(not value.is_Integer for value in raw_root_m):
            raise ArithmeticError(f"facet {index} is absent from the stable closure")
        closure_content = _content(raw_root_m)
        root_m = orientation * (raw_root_m / closure_content).applyfunc(int)
        div_m = _content(gram_m * root_m)
        root_n = (gram_m * root_m / div_m).applyfunc(int)
        root_fermat = transport * root_n

        square_l = int((root_l.T * gram_l * root_l)[0])
        div_l = _content(gram_l * root_l)
        square_m = int((root_m.T * gram_m * root_m)[0])
        square_n = int((root_n.T * gram_n * root_n)[0])
        div_n = _content(gram_n * root_n)
        reflection_m = _reflection(gram_m, root_m)
        reflection_n = _reflection(gram_n, root_n)
        reflection_fermat = _reflection(gram_fermat, root_fermat)
        factor_from_m = sp.Rational(2, div_m)
        factor_from_l = sp.Rational(2, div_m * closure_content)

        checks = {
            "source_square_archived": square_l == int(source["square"]),
            "source_pairing_column_and_orientation": (
                gram_l * (orientation * root_l)
                == sp.Matrix(source["fixed_pairing_column"])
            ),
            "primitive_in_M_N_and_Fermat": (
                _content(root_m)
                == _content(root_n)
                == _content(root_fermat)
                == 1
            ),
            "scaled_dual_square_law": square_n * div_m**2 == 2 * square_m,
            "exact_positive_halfspace_transport": (
                factor_from_m > 0
                and factor_from_l > 0
                and root_n.T * gram_n * gram_m
                == factor_from_m * root_m.T * gram_m
            ),
            "closure_reflection_integral": all(
                value.is_Integer for value in reflection_m
            ),
            "contragredient_reflection": reflection_n == reflection_m.T,
            "dual_and_Fermat_reflections_integral": (
                all(value.is_Integer for value in reflection_n)
                and all(value.is_Integer for value in reflection_fermat)
            ),
            "Fermat_reflection_intertwines": (
                reflection_fermat * transport == transport * reflection_n
            ),
            "Fermat_reflection_unimodular_involution_isometry": (
                reflection_fermat.det() == -1
                and reflection_fermat**2 == sp.eye(12)
                and reflection_fermat.T
                * gram_fermat
                * reflection_fermat
                == gram_fermat
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"facet {index} failed: "
                f"{[key for key, value in checks.items() if not value]}"
            )

        closure_types[(square_m, div_m)] += 1
        target_types[(square_n, div_n)] += 1
        source_to_target_types[
            (square_l, div_l, closure_content, square_m, div_m, square_n, div_n)
        ] += 1
        roots_m.append(root_m)
        roots_n.append(root_n)
        roots_fermat.append(root_fermat)
        records.append(
            {
                "facet_zero_based": index,
                "source_orbit_index_one_based": source["orbit_index"],
                "source_L_coordinates_unoriented": [int(value) for value in root_l],
                "source_presentation_orientation": orientation,
                "source_L_square": square_l,
                "source_L_divisibility": div_l,
                "coordinates_in_M_before_primitive_reduction_unoriented": [
                    int(value) for value in raw_root_m
                ],
                "M_coordinate_content": closure_content,
                "primitive_oriented_M_coordinates": [int(value) for value in root_m],
                "primitive_M_square": square_m,
                "primitive_M_divisibility": div_m,
                "primitive_oriented_N_coordinates": [int(value) for value in root_n],
                "primitive_oriented_Fermat_coordinates": [
                    int(value) for value in root_fermat
                ],
                "primitive_Fermat_square": square_n,
                "primitive_Fermat_divisibility": div_n,
                "positive_halfspace_factor_from_M": str(factor_from_m),
                "positive_halfspace_factor_from_original_L": str(factor_from_l),
                "checks": checks,
            }
        )

    pairwise_transport = all(
        (roots_n[left].T * gram_n * roots_n[right])[0]
        * _content(gram_m * roots_m[left])
        * _content(gram_m * roots_m[right])
        == 2 * (roots_m[left].T * gram_m * roots_m[right])[0]
        for left in range(29)
        for right in range(29)
    )
    roots_n_matrix = sp.Matrix.hstack(*roots_n).applyfunc(int)
    roots_fermat_matrix = sp.Matrix.hstack(*roots_fermat).applyfunc(int)
    hnf_n = hermite_normal_form(roots_n_matrix)
    hnf_fermat = hermite_normal_form(roots_fermat_matrix)
    carrier = sum(roots_n, sp.zeros(12, 1))
    carrier_fermat = transport * carrier
    carrier_pairings = [
        int((carrier.T * gram_n * root)[0]) for root in roots_n
    ]
    generation_checks = {
        "all_29_active_facets_transport": len(records) == 29,
        "all_archived_orientations_positive_unit": all(
            row["source_presentation_orientation"] == 1 for row in records
        ),
        "exact_pairwise_positive_Gram_rescaling": pairwise_transport,
        "primitive_N_roots_generate_full_lattice": hnf_n == sp.eye(12),
        "primitive_Fermat_roots_generate_full_fixed_divisor_lattice": (
            hnf_fermat == sp.eye(12)
        ),
        "sum_is_strict_timelike_interior": (
            int((carrier.T * gram_n * carrier)[0]) == 436
            and min(carrier_pairings) == 4
            and all(value > 0 for value in carrier_pairings)
            and carrier_fermat.T * gram_fermat * carrier_fermat
            == carrier.T * gram_n * carrier
        ),
        "closure_type_histogram": closure_types
        == Counter({(-2, 1): 17, (-1, 1): 7, (-4, 2): 5}),
        "target_type_histogram": target_types
        == Counter({(-4, 2): 17, (-2, 2): 7, (-2, 1): 5}),
    }
    if not all(generation_checks.values()):
        raise ArithmeticError(
            "generation/interior audit failed: "
            f"{[key for key, value in generation_checks.items() if not value]}"
        )

    vinberg_checks = {
        "literal_case_green": vinberg_case["all_checks_pass"] is True,
        "finite_volume_noncompact": (
            vinberg_case["finite_volume"] is True
            and vinberg_case["compact"] is False
            and vinberg_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
        ),
        "every_elliptic_edge_has_two_endpoints": (
            vinberg_case["bad_edge_count"] == 0
            and vinberg_case["elliptic_edge_count"] == 28798
        ),
        "all_parabolic_neighborhoods_good": (
            vinberg_case["bad_parabolic_neighborhood_count"] == 0
            and vinberg_case["connected_affine_component_count"] == 668
        ),
        "finite_and_ideal_vertices": (
            vinberg_case["finite_vertex_count"] == 3876
            and vinberg_case["ideal_vertex_count"] == 253
        ),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError(f"inherited Vinberg audit failed: {vinberg_checks}")

    checks = {
        "lattice_bridge": all(lattice_checks.values()),
        "all_facet_records": all(all(row["checks"].values()) for row in records),
        "full_generation_and_interior": all(generation_checks.values()),
        "literal_Vinberg_finite_volume": all(vinberg_checks.values()),
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "id": "A16-6-H11-29-FERMAT-CHAR3-SCALED-DUAL-CARRIER",
        "upstream": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[key],
            }
            for key, path in paths.items()
        },
        "theorem": (
            "The finite-volume noncompact A16_6 H11/29 chamber has an "
            "additional algebraic divisor-lattice carrier.  After its unique "
            "maximal reflection-stable odd index-eight closure M, the doubled "
            "dual M^vee(2) is explicitly the saturated C2-fixed divisor lattice "
            "U+E8(-2)+<-2>^2 on the supersingular Fermat quartic in "
            "characteristic three; all 29 positive halfspaces and reflections "
            "transport, and the primitive mirrors generate the whole fixed lattice."
        ),
        "stable_closure": {
            "operation": "unique maximal reflection-stable odd index-eight closure",
            "basis_columns_in_original_fixed_Q_space": _payload(closure_basis),
            "gram": _payload(gram_m),
            "determinant": int(gram_m.det()),
            "smith_invariants": _smith(gram_m),
            "parity": "odd",
            "primitive_mirror_types_square_divisibility": _histogram(closure_types),
        },
        "scaled_dual_bridge": {
            "operation": "N=M^vee(2), with Gram 2*Gram(M)^(-1)",
            "gram_N": _payload(gram_n),
            "N_determinant": int(gram_n.det()),
            "N_smith_invariants": _smith(gram_n),
            "N_parity": "even",
            "constructive_isometry_chain": {
                "SK_positive_gram": _payload(sk_gram),
                "dual_SK_gram": _payload(dual_sk_gram),
                "U_plus_A1x10_gram": _payload(a1_model),
                "U_plus_A1x10_basis_columns_in_N": _payload(a1_basis_in_n),
                "exceptional_U_split_basis": _payload(split),
                "exceptional_isotropic_vector_u": [int(value) for value in split[:, 0]],
                "isotropic_partner_v": [int(value) for value in split[:, 1]],
                "positive_complement_gram": _payload(positive_complement),
                "positive_complement_to_E8_2_A1_2_squared": _payload(
                    complement_isometry
                ),
                "standard_basis_columns_in_N": _payload(standard_basis_in_n),
            },
            "target": {
                "variety": "supersingular Fermat quartic in characteristic 3",
                "divisor_lattice": "NS(X)^C2=U+E8(-2)+<-2>^2",
                "archived_basis": (
                    "saturated fixed-divisor basis inside the published "
                    "rank-22 Neron-Severi basis"
                ),
                "existing_Fermat_section_wall_count": 206,
                "known_full_reflection_chamber_wall_count": 14,
                "archived_fixed_gram": _payload(gram_fermat),
            },
            "N_to_Fermat_fixed_isometry": _payload(transport),
            "isometry_determinant": int(transport.det()),
            "primitive_target_mirror_types_square_divisibility": _histogram(
                target_types
            ),
            "source_to_closure_to_target_type_histogram": _histogram(
                source_to_target_types
            ),
            "positive_halfspace_identity": (
                "For oriented primitive q in M, p=Gram(M)q/div_M(q) and "
                "Phi(x)=Gram(M)x satisfy (Phi(x),p)_N="
                "(2/div_M(q))*(x,q)_M."
            ),
            "checks": lattice_checks,
        },
        "facet_records": records,
        "full_lattice_generation": {
            "N_root_matrix_HNF": _payload(hnf_n),
            "Fermat_root_matrix_HNF": _payload(hnf_fermat),
            "generation_index": 1,
            "strict_interior_N_coordinates": [int(value) for value in carrier],
            "strict_interior_Fermat_coordinates": [
                int(value) for value in carrier_fermat
            ],
            "strict_interior_square": int((carrier.T * gram_n * carrier)[0]),
            "minimum_wall_pairing": min(carrier_pairings),
            "checks": generation_checks,
        },
        "inherited_literal_Vinberg_audit": {
            "criterion": vinberg_case["criterion"],
            "dimension": vinberg_case["dimension"],
            "facet_count": vinberg_case["facet_count"],
            "finite_vertex_count": vinberg_case["finite_vertex_count"],
            "ideal_vertex_count": vinberg_case["ideal_vertex_count"],
            "elliptic_edge_count": vinberg_case["elliptic_edge_count"],
            "connected_affine_component_count": vinberg_case[
                "connected_affine_component_count"
            ],
            "bad_edge_count": vinberg_case["bad_edge_count"],
            "bad_parabolic_neighborhood_count": vinberg_case[
                "bad_parabolic_neighborhood_count"
            ],
            "classification": vinberg_case["classification"],
            "checks": vinberg_checks,
        },
        "lineage_and_scope": {
            "prior_full_rank_complex_K3_carrier": (
                "the determinant-16 U+D8(-1)+A1(-1)^2 Brandhorst--Hofmann "
                "very-general Neron--Severi lattice"
            ),
            "new_positive_gain": (
                "an arithmetically different determinant-1024 carrier on an "
                "actual characteristic-three Fermat fixed divisor lattice, "
                "with coordinate-level scaled-dual wall transport"
            ),
            "target_already_has_other_chambers": (
                "the same fixed lattice already carries the Fermat 206-wall "
                "section and Alexeev's 14-wall full reflection chamber"
            ),
            "known_abstract_polytope": True,
            "not_a_new_abstract_polytope": True,
            "not_a_new_lattice_or_maximal_reflection_group": True,
            "not_the_complete_Fermat_nef_chamber": (
                "the imported arrangement has 29 mixed mirrors, including "
                "seventeen primitive square-minus-four divisor cuts; extra "
                "roots may refine it"
            ),
            "not_an_automorphism_fixed_chamber_claim": True,
            "minus_two_divisors_not_asserted_irreducible": True,
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
        raise ArithmeticError("A16 Fermat scaled-dual certificate failed")
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists():
            raise ArithmeticError(f"missing archive: {args.output}")
        if args.output.read_text(encoding="utf-8") != serialized:
            raise ArithmeticError("archived A16 Fermat scaled-dual payload changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()
