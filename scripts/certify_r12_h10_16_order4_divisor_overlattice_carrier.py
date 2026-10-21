#!/usr/bin/env python3
"""Certify an order-four divisor-overlattice carrier for R12 H10/16.

The archived C2b action on

    M = U + D8(-1) + A1(-1)^2

does not extend to a K3 automorphism with ``M`` primitive Picard.  It does,
however, have a unique nontrivial action-stable even isotropic discriminant
overgroup.  That index-two overgroup is

    Mbar = U + E8(-1) + A1(-1)^2,

with the order-two Picard action induced by the actual order-four
Brandhorst--Hofmann 0.4.4.1 / explicit Weierstrass family.  Consequently the
same surface contains an invariant finite-index divisor sublattice carrying
the archived complete 19-root arrangement and its H10/16 fixed section.

The 19-root chamber is a coarse divisor-root chamber (a union of full nef
chambers), because the saturated Picard lattice has additional roots.  Its
fixed section is nevertheless literally the frozen finite-volume Coxeter
polytope.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    induced_row_action,
    saturated_row_kernel,
)
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _q_histogram,
)
from scripts.certify_r12_h10_16_nikulin_obstruction import (
    _bits_vector,
    _combine_bits,
    _discriminant_action,
    _discriminant_numerator_mod_eight,
    _mod_two_kernel_basis,
    _rank_bits,
    _vector_bits,
)
from scripts.certify_roulleau_r14_c2d_special_k3_order4_lift import (
    smith_diagonal,
)
from scripts.certify_roulleau_ud8d4_bh_k3_lifts import (
    apply_discriminant_map,
    apply_reduced_fixed_map,
    discriminant_map,
    glue_pairs,
    verify_fixed_discriminant_isometry,
)
from scripts.certify_ue8a1x2_order4_weierstrass_family import (
    replay_brandhorst_hofmann_zip,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
FOLD = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
OLD_BOUNDARY = ROOT / "data/enumerations/r12_h10_16_nikulin_obstruction.json"
FAMILY = ROOT / "data/enumerations/ue8a1x2_order4_weierstrass_family.json"
OUTPUT = (
    ROOT
    / "data/enumerations/r12_h10_16_order4_divisor_overlattice_carrier.json"
)

EXPECTED_SOURCE_SHA256 = {
    "fold": "51751779daa723e24dcaf27ec5dd69463dc68e662da28751506672defe7cbe63",
    "old_boundary": "58b5e005bc7612e060728dcc75f2eed2f52cadf826ea767a345b2b54c0be38f9",
    "explicit_family": "c933994f12b272de2286513cf38d88ce2a5dec718367569d3babe2d0f4bf8431",
}

BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(value)
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def matrix_digest(matrix: sp.Matrix) -> str:
    raw = json.dumps(matrix_payload(matrix), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def payload_digest(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    raw = json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def negative_cartan(rank: int, edges: Iterable[tuple[int, int]]) -> sp.Matrix:
    gram = -2 * sp.eye(rank)
    for left, right in edges:
        gram[left, right] = gram[right, left] = 1
    return gram


def standard_row39_pair() -> tuple[sp.Matrix, sp.Matrix]:
    """Return U+E8(-1)+A1(-1)^2 and the A1 swap, in row convention."""

    hyperbolic_plane = sp.Matrix([[0, 1], [1, 0]])
    e8 = negative_cartan(
        8,
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)],
    )
    a1_pair = sp.diag(-2, -2)
    gram = sp.diag(hyperbolic_plane, e8, a1_pair)
    action = sp.diag(sp.eye(10), sp.Matrix([[0, 1], [1, 0]]))
    return gram, action


def parity_kernel_basis(functional: Sequence[int]) -> sp.Matrix:
    """Basis of {x in Z^n : functional dot x = 0 mod 2}."""

    bits = [int(value) & 1 for value in functional]
    if not any(bits):
        return sp.eye(len(bits))
    pivot = bits.index(1)
    rows: list[list[int]] = []
    for index in range(len(bits)):
        if index == pivot:
            continue
        row = [0] * len(bits)
        row[index] = 1
        if bits[index]:
            row[pivot] = 1
        rows.append(row)
    last = [0] * len(bits)
    last[pivot] = 2
    rows.append(last)
    return sp.Matrix(rows)


def span_bits(generators: Sequence[int]) -> frozenset[int]:
    return frozenset(
        _combine_bits(generators, coefficients)
        for coefficients in range(1 << len(generators))
    )


def all_binary_subspaces(dimension: int) -> set[frozenset[int]]:
    """Enumerate every subspace of F_2^dimension (67 when dimension=4)."""

    vectors = range(1 << dimension)
    return {
        span_bits(generators)
        for generator_count in range(dimension + 1)
        for generators in itertools.product(vectors, repeat=generator_count)
    }


def independent_basis(subspace: Iterable[int]) -> list[int]:
    basis: list[int] = []
    for vector in sorted(subspace):
        if _rank_bits([*basis, vector]) > len(basis):
            basis.append(vector)
    return basis


def overgroup_basis(
    rank: int, discriminant_basis_bits: Sequence[int], subgroup_basis: Sequence[int]
) -> sp.Matrix:
    """Rows in the original lattice coordinates for a 2-elementary overgroup."""

    glue_rows = []
    for coefficients in subgroup_basis:
        bits = _combine_bits(discriminant_basis_bits, coefficients)
        glue_rows.append([(bits >> index) & 1 for index in range(rank)])
    generators = 2 * sp.eye(rank)
    if glue_rows:
        generators = generators.col_join(sp.Matrix(glue_rows))
    return hermite_normal_form(generators.T).T / 2


def pair_data(gram: sp.Matrix, action: sp.Matrix) -> dict[str, Any]:
    fixed_basis = saturated_row_kernel(action - sp.eye(gram.rows))
    fixed_gram = sp.simplify(fixed_basis * gram * fixed_basis.T)
    coinvariant_basis = saturated_row_kernel(gram * fixed_basis.T)
    coinvariant_gram = sp.simplify(
        coinvariant_basis * gram * coinvariant_basis.T
    )
    coinvariant_action = induced_row_action(coinvariant_basis, action)
    return {
        "fixed_basis": fixed_basis,
        "fixed_gram": fixed_gram,
        "coinvariant_basis": coinvariant_basis,
        "coinvariant_gram": coinvariant_gram,
        "coinvariant_action": coinvariant_action,
        "glue": glue_pairs(fixed_basis, coinvariant_basis, gram),
    }


def exact_pair_conjugacy(
    source_gram: sp.Matrix,
    source_action: sp.Matrix,
    target_gram: sp.Matrix,
    target_action: sp.Matrix,
    reduced_fixed_map: sp.Matrix,
) -> dict[str, Any]:
    """Certify an integral action-pair isometry by fixed/coinvariant glue."""

    source = pair_data(source_gram, source_action)
    target = pair_data(target_gram, target_action)
    fixed_certificate = verify_fixed_discriminant_isometry(
        source["fixed_gram"], target["fixed_gram"], reduced_fixed_map
    )
    if not (
        source["coinvariant_gram"] == sp.Matrix([[-4]])
        and target["coinvariant_gram"] == sp.Matrix([[-4]])
        and source["coinvariant_action"] == sp.Matrix([[-1]])
        and target["coinvariant_action"] == sp.Matrix([[-1]])
    ):
        raise AssertionError("the rank-one rootless coinvariant changed")

    coinvariant_isometry = sp.Matrix([[1]])
    coinvariant_map = discriminant_map(
        coinvariant_isometry,
        source["coinvariant_gram"],
        target["coinvariant_gram"],
    )
    target_coin_moduli, _, _ = _discriminant_data(target["coinvariant_gram"])
    target_graph = {coin: fixed for fixed, coin in target["glue"]}
    all_glue_carried = True
    for source_fixed, source_coin in source["glue"]:
        target_fixed = apply_reduced_fixed_map(
            source_fixed,
            source["fixed_gram"],
            target["fixed_gram"],
            reduced_fixed_map,
        )
        target_coin = apply_discriminant_map(
            coinvariant_map, source_coin, target_coin_moduli
        )
        if target_graph.get(target_coin) != target_fixed:
            all_glue_carried = False
            break
    fixed_moduli, _, _ = _discriminant_data(source["fixed_gram"])
    fixed_length = sum(modulus > 1 for modulus in fixed_moduli)
    stable_lifting = source["fixed_gram"].rows >= fixed_length + 2
    if not (all_glue_carried and stable_lifting):
        raise AssertionError("the fixed/coinvariant action pairs do not glue")
    return {
        "source_fixed_rank_signature_determinant": [
            source["fixed_gram"].rows,
            list(exact_inertia(source["fixed_gram"])),
            int(source["fixed_gram"].det()),
        ],
        "target_fixed_rank_signature_determinant": [
            target["fixed_gram"].rows,
            list(exact_inertia(target["fixed_gram"])),
            int(target["fixed_gram"].det()),
        ],
        "source_fixed_smith": smith_diagonal(source["fixed_gram"]),
        "target_fixed_smith": smith_diagonal(target["fixed_gram"]),
        "source_fixed_q_histogram": _q_histogram(source["fixed_gram"]),
        "target_fixed_q_histogram": _q_histogram(target["fixed_gram"]),
        "fixed_discriminant_isometry": fixed_certificate,
        "fixed_rank_at_least_discriminant_length_plus_two": stable_lifting,
        "nikulin_stable_lifting_applies": True,
        "coinvariant_gram": matrix_payload(source["coinvariant_gram"]),
        "coinvariant_action": matrix_payload(source["coinvariant_action"]),
        "coinvariant_integral_isometry": [[1]],
        "coinvariant_discriminant_map": matrix_payload(coinvariant_map),
        "source_glue": [
            [list(fixed), list(coin)] for fixed, coin in source["glue"]
        ],
        "target_glue": [
            [list(fixed), list(coin)] for fixed, coin in target["glue"]
        ],
        "all_glue_elements_carried": True,
        "full_integral_lattice_action_pairs_conjugate": True,
    }


def build_certificate(external_zip: Path | None = None) -> dict[str, Any]:
    source_hashes = {
        "fold": sha256(FOLD),
        "old_boundary": sha256(OLD_BOUNDARY),
        "explicit_family": sha256(FAMILY),
    }
    if source_hashes != EXPECTED_SOURCE_SHA256:
        raise AssertionError("a pinned H10/16 or row39 source changed")
    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    old_boundary = json.loads(OLD_BOUNDARY.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    if not (
        fold["all_checks_pass"]
        and old_boundary["all_checks_pass"]
        and family["all_checks_pass"]
    ):
        raise AssertionError("an upstream certificate failed")

    stock = next(row for row in fold["stocks"] if row["stock_key"] == "R12_UD8A1x2")
    survivor = next(row for row in stock["survivors"] if row["class_index"] == 3)
    class_record = next(row for row in stock["classes"] if row["class_index"] == 3)
    geometry = survivor["geometry"]
    volume = geometry["vinberg_volume_criterion"]
    m_gram = sp.Matrix(stock["basis_gram"])
    wall_rows = sp.Matrix(stock["wall_coordinate_rows"])
    permutation = class_record["representative_generators_one_based"][0]
    m_action_column = sp.zeros(12)
    for source_index in range(12):
        m_action_column[:, source_index] = wall_rows[
            permutation[source_index] - 1, :
        ].T
    m_action = m_action_column.T
    archived_action = sp.Matrix(
        survivor["neron_severi_action"]["generator_action_matrices"][0]
    )
    wall_gram = wall_rows * m_gram * wall_rows.T
    if not (
        m_action_column == archived_action
        and m_action * m_gram * m_action.T == m_gram
        and m_action**2 == sp.eye(12)
        and stock["wall_count"] == 19
        and all(wall_gram[index, index] == -2 for index in range(19))
        and survivor["all_checks_pass"]
        and all(geometry["checks"].values())
        and geometry["facet_count"] == 16
        and volume["decision"] == "finite_volume_noncompact"
        and volume["finite_volume_bad_edges"] == []
        and volume["edge_extension_histogram"] == {"2": 426}
        and volume["parabolic_neighborhood_audit_count"] == 14
    ):
        raise AssertionError("the archived H10/16 chamber changed")

    # Complete inverse-overextension census in the four-bit discriminant form.
    discriminant_basis = _mod_two_kernel_basis(m_gram)
    finite_action = _discriminant_action(discriminant_basis, m_action_column)
    finite_action = finite_action.applyfunc(lambda value: int(value) & 1)

    def q_numerator(coefficients: int) -> int:
        return _discriminant_numerator_mod_eight(
            m_gram, _combine_bits(discriminant_basis, coefficients)
        )

    def image_coefficients(coefficients: int) -> int:
        image = finite_action * _bits_vector(coefficients, 4)
        return _vector_bits([int(value) for value in image])

    subspaces = all_binary_subspaces(4)
    isotropic_subspaces = [
        subspace
        for subspace in subspaces
        if all(q_numerator(element) == 0 for element in subspace)
    ]
    stable_subspaces = [
        subspace
        for subspace in subspaces
        if {image_coefficients(element) for element in subspace} == set(subspace)
    ]
    stable_isotropic_subspaces = [
        subspace
        for subspace in isotropic_subspaces
        if {image_coefficients(element) for element in subspace} == set(subspace)
    ]
    nontrivial = [subspace for subspace in stable_isotropic_subspaces if len(subspace) > 1]
    if not (
        discriminant_basis == [344, 600, 1112, 2058]
        and finite_action
        == sp.Matrix(
            [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
        )
        and len(subspaces) == 67
        and len(isotropic_subspaces) == 4
        and len(stable_subspaces) == 27
        and len(stable_isotropic_subspaces) == 2
        and nontrivial == [frozenset({0, 1})]
    ):
        raise AssertionError("the complete stable-isotropic census changed")

    subgroup_basis = independent_basis(nontrivial[0])
    over_basis = overgroup_basis(12, discriminant_basis, subgroup_basis)
    over_gram = sp.simplify(over_basis * m_gram * over_basis.T)
    over_action = sp.simplify(over_basis * m_action * over_basis.inv())
    standard_gram, standard_action = standard_row39_pair()
    over_to_row39 = exact_pair_conjugacy(
        over_gram,
        over_action,
        standard_gram,
        standard_action,
        sp.Matrix([[1]]),
    )
    if not (
        abs(over_basis.det()) == sp.Rational(1, 2)
        and all(value.q == 1 for value in over_gram)
        and all(over_gram[index, index] % 2 == 0 for index in range(12))
        and over_gram.det() == -4
        and exact_inertia(over_gram) == (1, 11, 0)
        and smith_diagonal(over_gram) == [*[1] * 10, 2, 2]
        and _q_histogram(over_gram) == {"0": 1, "1": 1, "3/2": 2}
        and all(value.q == 1 for value in over_action)
        and over_action * over_gram * over_action.T == over_gram
        and over_action**2 == sp.eye(12)
        and sp.factor(over_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 11 * (sp.Symbol("lambda") + 1)
        and over_to_row39["full_integral_lattice_action_pairs_conjugate"]
    ):
        raise AssertionError("the unique even overgroup is not row39")

    # A direct index-two divisor sublattice in the standard row39 marking.
    # Coordinates are U, the displayed E8 simple basis, and the two A1 roots;
    # requiring the seventh E8 coefficient (global coordinate 8) to be even
    # selects one particularly small witness.
    parity = [0] * 12
    parity[8] = 1
    divisor_basis = parity_kernel_basis(parity)
    divisor_gram = sp.simplify(divisor_basis * standard_gram * divisor_basis.T)
    divisor_action = induced_row_action(divisor_basis, standard_action)
    divisor_to_archive = exact_pair_conjugacy(
        divisor_gram,
        divisor_action,
        m_gram,
        m_action,
        sp.Matrix([[1, 1, 1], [0, 1, 1], [0, 0, 1]]),
    )
    if not (
        abs(divisor_basis.det()) == 2
        and divisor_gram.det() == -16
        and exact_inertia(divisor_gram) == (1, 11, 0)
        and smith_diagonal(divisor_gram) == [*[1] * 8, 2, 2, 2, 2]
        and _q_histogram(divisor_gram)
        == {"0": 4, "1": 4, "1/2": 2, "3/2": 6}
        and divisor_action * divisor_gram * divisor_action.T == divisor_gram
        and divisor_action**2 == sp.eye(12)
        and divisor_to_archive["full_integral_lattice_action_pairs_conjugate"]
    ):
        raise AssertionError("the direct row39 divisor sublattice failed")

    family_handoff = family["brandhorst_hofmann_handoff"]
    if not (
        family["family"]["global_order"] == 4
        and family_handoff["id"] == "0.4.4.1"
        and family_handoff["full_NS_action_pairs_integrally_conjugate"]
        and family["neron_severi"]["rank_signature_determinant"]
        == [12, [1, 11, 0], -4]
        and family["picard_action"]["coinvariant"]["rootless"]
        and family["picard_action"]["coinvariant"]["gram"] == [[-4]]
    ):
        raise AssertionError("the explicit order-four family handoff changed")
    external = (
        replay_brandhorst_hofmann_zip(external_zip) if external_zip else None
    )

    checks = {
        "archived_complete_19_root_chamber": stock["wall_count"] == 19,
        "archived_all_wall_norms_minus_two": all(
            wall_gram[index, index] == -2 for index in range(19)
        ),
        "archived_c2b_action_reconstructed": m_action_column == archived_action,
        "complete_four_bit_subspace_census": len(subspaces) == 67,
        "unique_nontrivial_stable_even_isotropic_overgroup": len(nontrivial) == 1,
        "overgroup_is_even_index_two": abs(over_basis.det()) == sp.Rational(1, 2),
        "overgroup_is_row39_action_pair": over_to_row39[
            "full_integral_lattice_action_pairs_conjugate"
        ],
        "direct_divisor_sublattice_is_archived_c2b_pair": divisor_to_archive[
            "full_integral_lattice_action_pairs_conjugate"
        ],
        "actual_order_four_weierstrass_carrier": family["family"][
            "global_order"
        ]
        == 4,
        "full_picard_coinvariant_rootless": family["picard_action"][
            "coinvariant"
        ]["rootless"],
        "strict_fixed_section_interior_exists": geometry["checks"][
            "strict_timelike_interior_exists"
        ],
        "vinberg_finite_volume_exact": volume["decision"]
        == "finite_volume_noncompact",
        "vinberg_no_bad_audits": not volume["finite_volume_bad_edges"],
    }
    if not all(checks.values()):
        raise AssertionError(f"a positive carrier check failed: {checks}")

    certificate: dict[str, Any] = {
        "schema": "r12_h10_16_order4_divisor_overlattice_carrier_v1",
        "all_checks_pass": True,
        "theorem": (
            "The very general explicit BH 0.4.4.1 K3 family with "
            "NS=U+E8(-1)+A1(-1)^2 and order-four automorphism sigma contains "
            "a sigma-stable index-two divisor sublattice isometric, as an "
            "action pair, to the archived C2b lattice "
            "U+D8(-1)+A1(-1)^2.  Its complete 19-root divisor chamber has "
            "the finite-volume noncompact H10/16 fixed Coxeter section."
        ),
        "source_sha256": source_hashes,
        "case_id": "R12-D8A1x2-C2b-H10-16",
        "old_primitive_boundary_retained": {
            "primitive_picard_lift_exists": False,
            "reason": old_boundary["k3_consequence"]["reason"],
            "new_escape": (
                "M is not primitive in NS(X): its unique action-stable even "
                "index-two saturation is the actual row39 Picard lattice"
            ),
            "contradiction": False,
        },
        "inverse_overextension_census": {
            "discriminant_group": "A_M=(Z/2)^4",
            "half_vector_basis_bits": discriminant_basis,
            "action_on_coefficient_columns": matrix_payload(finite_action),
            "quadratic_classes": [
                {
                    "coefficient_bits_as_integer": coefficients,
                    "coefficient_vector": [
                        (coefficients >> index) & 1 for index in range(4)
                    ],
                    "bGb_mod_8": q_numerator(coefficients),
                    "action_image_bits_as_integer": image_coefficients(coefficients),
                }
                for coefficients in range(16)
            ],
            "all_subspace_count": len(subspaces),
            "dimension_histogram": {
                str(dimension): sum(
                    len(subspace) == 1 << dimension for subspace in subspaces
                )
                for dimension in range(5)
            },
            "totally_q_isotropic_subspace_count": len(isotropic_subspaces),
            "action_stable_subspace_count": len(stable_subspaces),
            "action_stable_totally_q_isotropic_subspace_count": len(
                stable_isotropic_subspaces
            ),
            "nontrivial_action_stable_even_overgroup_count": len(nontrivial),
            "unique_nontrivial_subgroup_elements": sorted(nontrivial[0]),
            "unique_nontrivial_subgroup_basis_coefficients": subgroup_basis,
            "unique_glue_half_vector_bits": _combine_bits(
                discriminant_basis, subgroup_basis[0]
            ),
            "conclusion": (
                "there is exactly one nontrivial h-stable even isotropic "
                "overgroup; it is index two and no further even h-stable "
                "overextension exists"
            ),
        },
        "unique_saturation": {
            "basis_rows_in_archived_M_coordinates": matrix_payload(over_basis),
            "gram_sha256": matrix_digest(over_gram),
            "action_sha256": matrix_digest(over_action),
            "rank_signature_determinant": [12, [1, 11, 0], -4],
            "smith": smith_diagonal(over_gram),
            "q_histogram_mod_2Z": _q_histogram(over_gram),
            "isometry_class": "U+E8(-1)+A1(-1)^2",
            "index_over_M": 2,
            "induced_action_characteristic_polynomial": "(x-1)^11*(x+1)",
            "exact_action_pair_conjugacy_to_row39": over_to_row39,
        },
        "direct_row39_divisor_sublattice": {
            "standard_row39_basis_order": (
                "U, displayed E8 simple roots, two A1 roots"
            ),
            "parity_functional": parity,
            "description": (
                "the coefficient of the seventh displayed E8 simple root "
                "(global zero-based coordinate 8) is even"
            ),
            "basis_rows_in_standard_row39_coordinates": matrix_payload(
                divisor_basis
            ),
            "index_in_NS": 2,
            "gram_sha256": matrix_digest(divisor_gram),
            "action_sha256": matrix_digest(divisor_action),
            "rank_signature_determinant": [12, [1, 11, 0], -16],
            "smith": smith_diagonal(divisor_gram),
            "q_histogram_mod_2Z": _q_histogram(divisor_gram),
            "exact_action_pair_conjugacy_to_archived_C2b": divisor_to_archive,
        },
        "explicit_algebraic_carrier": {
            "family_certificate": str(FAMILY.relative_to(ROOT)),
            "equation": family["family"]["affine_equation"],
            "A": family["family"]["A"],
            "C": family["family"]["C"],
            "automorphism": family["family"]["automorphism"],
            "automorphism_global_order": 4,
            "automorphism_picard_order": 2,
            "very_general_NS": "U+E8(-1)+A1(-1)^2",
            "brandhorst_hofmann_id": family_handoff["id"],
            "brandhorst_hofmann_row_one_based": family_handoff["row_one_based"],
            "brandhorst_hofmann_row_sha256_with_newline": family_handoff[
                "row_sha256_with_newline"
            ],
            "brandhorst_hofmann_zip_sha256": BH_ZIP_SHA256,
            "external_database_replay_performed": external is not None,
            "external_database_replay": external,
        },
        "divisor_chamber_handoff": {
            "ambient_picard_is_saturation_of_M": True,
            "archived_M_root_count": stock["wall_count"],
            "all_archived_roots_remain_algebraic_minus_two_classes": True,
            "archived_complete_M_reflection_chamber_retained": True,
            "full_NS_nef_chamber_may_be_a_strict_refinement": True,
            "scope": (
                "the 19 walls define an effective divisor-root subsystem "
                "chamber, a union of full NS Weyl chambers; "
                "this is not claimed to be the full nef cone of the row39 surface"
            ),
            "fixed_kahler_selection": (
                "the full Picard coinvariant is <-4>, hence rootless.  No "
                "Picard root is orthogonal to the entire sigma-fixed positive "
                "space.  Choose a generic point in the archived strict fixed "
                "interior avoiding the locally finite extra-root hyperplanes; "
                "the NS Weyl chamber containing it is sigma-stable and can be "
                "marked nef"
            ),
            "effectivity": (
                "with that invariant ample class, Riemann--Roch chooses the "
                "archived positive sign of every one of the 19 algebraic "
                "minus-two roots, so all are effective divisor classes"
            ),
        },
        "hyperbolic_polytope": {
            "id": "R12-D8A1x2-C2b-H10-16",
            "dimension": geometry["hyperbolic_dimension"],
            "facets": geometry["facet_count"],
            "classification": volume["decision"],
            "finite_vertices": geometry["finite_vertex_count"],
            "cusps": geometry["ideal_vertex_count"],
            "hyperideal_vertices": geometry["hyperideal_vertex_count"],
            "elliptic_edges": volume["elliptic_edge_count"],
            "edge_extension_histogram": volume["edge_extension_histogram"],
            "endpoint_type_histogram": volume["endpoint_type_histogram"],
            "parabolic_neighborhood_audit_count": volume[
                "parabolic_neighborhood_audit_count"
            ],
            "bad_edge_count": len(volume["finite_volume_bad_edges"]),
            "canonical_full_audit_sha256": volume[
                "canonical_full_audit_sha256"
            ],
            "simple_roots_complete": volume["simple_roots_complete"],
        },
        "novelty_scope": (
            "This is a new project-level algebraic divisor-overlattice "
            "realization of the known finite-volume noncompact H10/16 chamber. "
            "It is not a new abstract Coxeter polytope, arithmetic class, full "
            "nef chamber, or worldwide-priority claim."
        ),
        "checks": checks,
    }
    certificate["payload_sha256"] = payload_digest(certificate)
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--brandhorst-hofmann-zip", type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(args.brandhorst_hofmann_zip)
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if certificate != archived:
            raise AssertionError("the archived positive carrier certificate changed")
        print(f"verified {args.output}")
        return
    args.output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")
    print(certificate["payload_sha256"])


if __name__ == "__main__":
    main()
