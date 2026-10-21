#!/usr/bin/env python3
"""Certify two genuine K3 divisor-sublattice carriers of H10/45.

The archived H10/45 action lives on

    M = U + D4(-1)^2 + A1(-1)^2.

It was previously not known to come from a K3 automorphism.  The key point is
that ``M`` need not be the whole Neron--Severi lattice.  In the canonical
Brandhorst--Hofmann Picard bases, the following stable finite-index
sublattices carry the exact H10/45 pair:

* index 2 in row 0.4.4.2, by making coordinate 10 even;
* index 4 in row 0.4.4.1, by making coordinates 10 and 12 even.

Both rows are genuine purely non-symplectic order-four K3 actions and both
have explicit Weierstrass families elsewhere in this repository.  This
script checks the sublattices, the fixed/coinvariant discriminant forms, the
primitive eigen-glue, the Nikulin lifting hypotheses for exact integral pair
conjugacy, the finite overextension census, and the frozen Vinberg audit.

The resulting 45-wall object is an effective divisor-subsystem chamber.  It
is not asserted to be the full nef chamber of the carrier surface: extra
roots of the finite-index ambient Neron--Severi lattice can refine it.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence
import zipfile

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _q_histogram,
)
from scripts.certify_r12_h10_45_e7_embedding_search_boundary import (
    discriminant_action,
    fixed_q2_histogram,
    gf2_nullspace,
)
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    saturated_row_kernel,
    smith_diagonal,
)
from scripts.certify_roulleau_rank12_bh_k3_lifts import (
    BH_ZIP_SHA256,
    MEMBERS,
    NS_GRAMS,
    PURE_ACTION,
    _derive_external_row,
    _parse_bh_row,
    sha256_bytes,
)
from scripts.certify_roulleau_ud8d4_bh_k3_lifts import (
    apply_discriminant_map,
    apply_reduced_fixed_map,
    discriminant_map,
    exact_order,
    glue_pairs,
    verify_fixed_discriminant_isometry,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data"
    / "enumerations"
    / "roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "r12_h10_45_bh_divisor_sublattice_carrier.json"
)
UE8_FAMILY = ROOT / "data" / "enumerations" / "ue8a1x2_order4_weierstrass_family.json"
UD8_FAMILY = ROOT / "data" / "enumerations" / "ud8a1x2_order4_weierstrass_family.json"

EXPECTED_SOURCE_SHA256 = (
    "448bcc47937a2763bbd4ea97466fbc48924dc862e88ba8719e5ab98b0cc4c916"
)
BH_MEMBER = "K3Groups/results/purely_ns/order04.txt"
BH_MEMBER_SHA256 = MEMBERS["pure04"]["sha256"]
BH_ROW_COUNT = MEMBERS["pure04"]["rows"]
BH_ROWS = {
    "0.4.4.1": {
        "one_based": 39,
        "sha256_with_newline": (
            "03ef01ed30b98e3d0fa419cb8275577c90594868fcd3400c48413555b1248452"
        ),
    },
    "0.4.4.2": {
        "one_based": 40,
        "sha256_with_newline": (
            "bbba67834e034344032bbf137a1f2f1499f06f3adac123ce2b580e0bf7614849"
        ),
    },
}

TARGET_CASE = "R12-D4x2A1x2-S6transposition-H10-45"
TARGET_SMITH = [*[1] * 6, *[2] * 6]
TARGET_Q = {"0": 16, "1": 16, "1/2": 12, "3/2": 20}
TARGET_FIXED_SMITH = [*[1] * 6, *[2] * 4, 4]
TARGET_FIXED_Q = {
    "0": 8,
    "1": 8,
    "1/2": 8,
    "1/4": 8,
    "3/2": 8,
    "3/4": 4,
    "5/4": 8,
    "7/4": 12,
}

# This is a map A_F(target) -> A_F(carrier) on the nontrivial Smith
# coordinates (2,2,2,2,4).  All 64 elements are checked below.
FIXED_DISCRIMINANT_MAP = sp.Matrix(
    [
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 1, 1],
        [0, 0, 1, 0, 0],
        [0, 0, 2, 0, 1],
    ]
)

CARRIERS: tuple[dict[str, Any], ...] = (
    {
        "id": "BH-0.4.4.2-index2-H10-45",
        "bh_id": "0.4.4.2",
        "ambient_lattice": "U+D8(-1)+A1(-1)^2",
        "ambient_determinant": -16,
        "index": 2,
        "even_coordinates_one_based": [10],
        "family_path": UD8_FAMILY,
        "family_lattice": "U+D8(-1)+A1(-1)^2",
    },
    {
        "id": "BH-0.4.4.1-index4-H10-45",
        "bh_id": "0.4.4.1",
        "ambient_lattice": "U+E8(-1)+A1(-1)^2",
        "ambient_determinant": -4,
        "index": 4,
        "even_coordinates_one_based": [10, 12],
        "family_path": UE8_FAMILY,
        "family_lattice": "U+E8(-1)+A1(-1)^2",
    },
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def matrix_digest(matrix: sp.Matrix) -> str:
    return hashlib.sha256(
        json.dumps(matrix_payload(matrix), separators=(",", ":")).encode()
    ).hexdigest()


def _target_pair(source: dict[str, Any]) -> dict[str, Any]:
    gram = sp.Matrix(source["coordinate_reconstruction"]["basis_gram"])
    action = sp.eye(12)
    action[1, 1] = action[2, 2] = 0
    action[1, 2] = action[2, 1] = 1
    fixed = saturated_row_kernel(action - sp.eye(12))
    coinvariant = saturated_row_kernel(action + sp.eye(12))
    fixed_gram = fixed * gram * fixed.T
    coinvariant_gram = coinvariant * gram * coinvariant.T
    half_basis = gf2_nullspace(gram)
    discr_action = discriminant_action(action, half_basis)
    if not (
        gram.det() == -64
        and exact_inertia(gram) == (1, 11, 0)
        and smith_diagonal(gram) == TARGET_SMITH
        and _q_histogram(gram) == TARGET_Q
        and action * gram * action.T == gram
        and exact_order(action, 2)
        and fixed_gram.det() == 64
        and exact_inertia(fixed_gram) == (1, 10, 0)
        and smith_diagonal(fixed_gram) == TARGET_FIXED_SMITH
        and _q_histogram(fixed_gram) == TARGET_FIXED_Q
        and coinvariant_gram == sp.Matrix([[-4]])
        and fixed_q2_histogram(discr_action, half_basis, gram)
        == {"0": 8, "1": 8, "2": 8, "3": 8}
    ):
        raise AssertionError("the archived H10/45 lattice-action pair changed")
    return {
        "gram": gram,
        "action": action,
        "fixed": fixed,
        "coinvariant": coinvariant,
        "fixed_gram": fixed_gram,
        "coinvariant_gram": coinvariant_gram,
        "glue": glue_pairs(fixed, coinvariant, gram),
    }


def _inclusion_basis(even_coordinates_one_based: Sequence[int]) -> sp.Matrix:
    basis = sp.eye(12)
    for coordinate in even_coordinates_one_based:
        basis[coordinate - 1, coordinate - 1] = 2
    return basis


def _family_record(specification: dict[str, Any]) -> dict[str, Any]:
    family = json.loads(specification["family_path"].read_text(encoding="utf-8"))
    exactness = family["neron_severi"]["generic_exactness"]
    if not (
        family["all_checks_pass"]
        and family["brandhorst_hofmann_handoff"]["id"]
        == specification["bh_id"]
        and family["neron_severi"]["lower_bound_lattice"]
        == specification["family_lattice"]
        and exactness["parameter_moduli_dimension"] == 4
        and "very general" in exactness["scope"]
        and family["family"]["global_order"] == 4
        and family["family"]["holomorphic_two_form"].startswith("dt")
    ):
        raise AssertionError("an explicit order-four Weierstrass carrier changed")
    return {
        "archive": str(specification["family_path"].relative_to(ROOT)),
        "payload_sha256": family["payload_sha256"],
        "equation": family["family"]["affine_equation"],
        "automorphism": family["family"]["automorphism"],
        "very_general_NS": specification["family_lattice"],
        "very_general_NS_exact": True,
        "global_automorphism_order": 4,
        "brandhorst_hofmann_id": specification["bh_id"],
        "formula_automorphism_role": (
            "supplies the full order-four K3 Hodge action; the simultaneous "
            "effectivity theorem may replace it by a Weyl-conjugate Torelli "
            "automorphism on the same surface"
        ),
    }


def _carrier_record(
    specification: dict[str, Any], target: dict[str, Any]
) -> dict[str, Any]:
    ambient_gram = sp.Matrix(NS_GRAMS[specification["bh_id"]])
    ambient_action = sp.Matrix(PURE_ACTION)
    inclusion = _inclusion_basis(specification["even_coordinates_one_based"])
    sub_gram = sp.simplify(inclusion * ambient_gram * inclusion.T)
    sub_action = sp.simplify(inclusion * ambient_action * inclusion.inv())
    ambient_coinvariant = saturated_row_kernel(ambient_action + sp.eye(12))
    ambient_coinvariant_gram = (
        ambient_coinvariant * ambient_gram * ambient_coinvariant.T
    )
    fixed = saturated_row_kernel(sub_action - sp.eye(12))
    coinvariant = saturated_row_kernel(sub_action + sp.eye(12))
    fixed_gram = fixed * sub_gram * fixed.T
    coinvariant_gram = coinvariant * sub_gram * coinvariant.T

    fixed_certificate = verify_fixed_discriminant_isometry(
        target["fixed_gram"], fixed_gram, FIXED_DISCRIMINANT_MAP
    )
    coin_map = discriminant_map(
        sp.eye(1), target["coinvariant_gram"], coinvariant_gram
    )
    target_coin_moduli, _, _ = _discriminant_data(coinvariant_gram)
    carrier_glue = glue_pairs(fixed, coinvariant, sub_gram)
    carrier_graph = {coin: fixed_class for fixed_class, coin in carrier_glue}
    glue_carried = True
    for source_fixed, source_coin in target["glue"]:
        target_fixed = apply_reduced_fixed_map(
            source_fixed,
            target["fixed_gram"],
            fixed_gram,
            FIXED_DISCRIMINANT_MAP,
        )
        target_coin = apply_discriminant_map(
            coin_map, source_coin, target_coin_moduli
        )
        glue_carried &= carrier_graph.get(target_coin) == target_fixed

    fixed_length = len([value for value in smith_diagonal(fixed_gram) if value > 1])
    checks = {
        "ambient_even_hyperbolic": (
            exact_inertia(ambient_gram) == (1, 11, 0)
            and all(ambient_gram[index, index] % 2 == 0 for index in range(12))
        ),
        "ambient_determinant": (
            ambient_gram.det() == specification["ambient_determinant"]
        ),
        "ambient_action_isometry": (
            ambient_action * ambient_gram * ambient_action.T == ambient_gram
            and exact_order(ambient_action, 2)
        ),
        "ambient_coinvariant_rootless": (
            ambient_coinvariant_gram == sp.Matrix([[-4]])
        ),
        "stable_sublattice_index": (
            abs(int(inclusion.det())) == specification["index"]
            and all(sp.Rational(value).q == 1 for value in sub_action)
            and inclusion * ambient_action == sub_action * inclusion
        ),
        "target_lattice_invariants": (
            sub_gram.det() == -64
            and exact_inertia(sub_gram) == (1, 11, 0)
            and smith_diagonal(sub_gram) == TARGET_SMITH
            and _q_histogram(sub_gram) == TARGET_Q
        ),
        "target_action_character": (
            exact_order(sub_action, 2)
            and sub_action * sub_gram * sub_action.T == sub_gram
            and sp.factor(sub_action.charpoly().as_expr())
            == (sp.Symbol("lambda") - 1) ** 11
            * (sp.Symbol("lambda") + 1)
        ),
        "target_fixed_lattice": (
            fixed_gram.det() == 64
            and exact_inertia(fixed_gram) == (1, 10, 0)
            and smith_diagonal(fixed_gram) == TARGET_FIXED_SMITH
            and _q_histogram(fixed_gram) == TARGET_FIXED_Q
        ),
        "target_coinvariant": coinvariant_gram == sp.Matrix([[-4]]),
        "fixed_discriminant_isometry": fixed_certificate["bijective"]
        and fixed_certificate["quadratic_form_preserved_mod_2Z"],
        "all_eigen_glue_carried": glue_carried,
        "nikulin_rank_gate": fixed_gram.rows >= fixed_length + 2,
    }
    if not all(checks.values()):
        raise AssertionError(
            f"the divisor-sublattice carrier {specification['bh_id']} failed: {checks}"
        )

    return {
        "carrier_id": specification["id"],
        "brandhorst_hofmann": {
            "id": specification["bh_id"],
            "effective_group": "C4",
            "Picard_image": "C2",
            "member": BH_MEMBER,
            "row_one_based": BH_ROWS[specification["bh_id"]]["one_based"],
            "row_sha256_with_newline": BH_ROWS[specification["bh_id"]][
                "sha256_with_newline"
            ],
        },
        "ambient_NS": {
            "lattice": specification["ambient_lattice"],
            "rank_signature_determinant": [
                12,
                [1, 11],
                int(ambient_gram.det()),
            ],
            "gram_sha256": matrix_digest(ambient_gram),
            "Picard_coinvariant_lattice": "<-4>",
            "Picard_coinvariant_rootless": True,
        },
        "divisor_sublattice": {
            "lattice": "U+D4(-1)^2+A1(-1)^2",
            "ambient_index": specification["index"],
            "ambient_quotient": (
                "Z/2" if specification["index"] == 2 else "(Z/2)^2"
            ),
            "canonical_BH_coordinates_required_even_one_based": specification[
                "even_coordinates_one_based"
            ],
            "basis_in_ambient_NS": matrix_payload(inclusion),
            "gram": matrix_payload(sub_gram),
            "gram_sha256": matrix_digest(sub_gram),
            "action": matrix_payload(sub_action),
            "action_sha256": matrix_digest(sub_action),
            "rank_signature_determinant": [12, [1, 11], -64],
            "smith": TARGET_SMITH,
            "discriminant_q_histogram": TARGET_Q,
            "nonprimitive_full_rank_divisor_sublattice": True,
        },
        "exact_H10_45_action_pair": {
            "fixed_rank_signature_determinant": [11, [1, 10], 64],
            "fixed_smith": TARGET_FIXED_SMITH,
            "fixed_discriminant_q_histogram": TARGET_FIXED_Q,
            "coinvariant_lattice": "<-4>",
            "eigenstack_index": len(carrier_glue),
            "fixed_discriminant_isometry": fixed_certificate,
            "coinvariant_discriminant_map": matrix_payload(coin_map),
            "all_eigen_glue_pairs_carried": glue_carried,
            "indefinite_fixed_rank_at_least_discriminant_length_plus_two": True,
            "nikulin_uniqueness_and_surjectivity_apply": True,
            "full_integral_lattice_action_pairs_conjugate": True,
        },
        "explicit_weierstrass_carrier": _family_record(specification),
        "checks": checks,
        "all_checks_pass": True,
    }


def _add(
    left: Sequence[int], right: Sequence[int], moduli: Sequence[int]
) -> tuple[int, ...]:
    return tuple((left[index] + right[index]) % moduli[index] for index in range(len(moduli)))


def _closure(
    generators: Iterable[Sequence[int]], moduli: Sequence[int]
) -> frozenset[tuple[int, ...]]:
    zero = tuple(0 for _ in moduli)
    result = {zero}
    for generator_sequence in generators:
        generator = tuple(generator_sequence)
        old = tuple(result)
        multiples = [zero]
        current = zero
        for _ in range(math.prod(moduli)):
            current = _add(current, generator, moduli)
            if current == zero:
                break
            multiples.append(current)
        result = {
            _add(left, right, moduli) for left in old for right in multiples
        }
    return frozenset(result)


def _binary_discriminant_isometry_exists(
    source_gram: sp.Matrix, target_gram: sp.Matrix
) -> bool:
    """Decide isometry of the small 2-elementary discriminant forms exactly."""

    source_moduli, _, source_q = _discriminant_data(source_gram)
    target_moduli, _, target_q = _discriminant_data(target_gram)
    source_indices = [i for i, modulus in enumerate(source_moduli) if modulus > 1]
    target_indices = [i for i, modulus in enumerate(target_moduli) if modulus > 1]
    if not (
        len(source_indices) == len(target_indices)
        and all(source_moduli[i] == 2 for i in source_indices)
        and all(target_moduli[i] == 2 for i in target_indices)
    ):
        return False
    dimension = len(source_indices)
    reduced_elements = list(itertools.product((0, 1), repeat=dimension))

    def q_table(
        discriminant_gram: sp.Matrix,
        indices: Sequence[int],
        rank: int,
    ) -> dict[tuple[int, ...], sp.Rational]:
        result = {}
        for reduced in reduced_elements:
            full = sp.zeros(rank, 1)
            for position, index in enumerate(indices):
                full[index] = reduced[position]
            result[reduced] = sp.Rational(
                (full.T * discriminant_gram * full)[0]
            ) % 2
        return result

    source_values = q_table(source_q, source_indices, source_gram.rows)
    target_values = q_table(target_q, target_indices, target_gram.rows)
    for columns in itertools.product(reduced_elements, repeat=dimension):
        images = []
        preserves_q = True
        for source in reduced_elements:
            image = tuple(
                sum(source[column] * columns[column][row] for column in range(dimension))
                % 2
                for row in range(dimension)
            )
            images.append(image)
            if source_values[source] != target_values[image]:
                preserves_q = False
                break
        if preserves_q and len(set(images)) == 2**dimension:
            return True
    return False


def _overextension_census(target: dict[str, Any]) -> dict[str, Any]:
    """Enumerate every small even fixed-lattice extension compatible with glue."""

    gram = target["fixed_gram"]
    moduli, left_inverse, discr_gram = _discriminant_data(gram)
    elements = list(itertools.product(*(range(modulus) for modulus in moduli)))
    gamma = next(
        fixed_class
        for fixed_class, coin_class in target["glue"]
        if any(coin_class)
    )

    def q(element: Sequence[int]) -> sp.Rational:
        vector = sp.Matrix(element)
        return sp.Rational((vector.T * discr_gram * vector)[0]) % 2

    def bilinear(left: Sequence[int], right: Sequence[int]) -> sp.Rational:
        return sp.Rational(
            (sp.Matrix(left).T * discr_gram * sp.Matrix(right))[0]
        ) % 1

    isotropic_perpendicular = [
        element
        for element in elements
        if q(element) == 0 and bilinear(element, gamma) == 0
    ]
    subgroups: set[frozenset[tuple[int, ...]]] = set()
    for element in isotropic_perpendicular:
        subgroup = _closure([element], moduli)
        if len(subgroup) in (2, 4) and all(q(value) == 0 for value in subgroup):
            subgroups.add(subgroup)
    for left, right in itertools.combinations(isotropic_perpendicular, 2):
        subgroup = _closure([left, right], moduli)
        if (
            len(subgroup) == 4
            and all(q(value) == 0 for value in subgroup)
            and all(
                bilinear(x, y) == 0 for x in subgroup for y in subgroup
            )
        ):
            subgroups.add(subgroup)

    def extended_ambient(
        subgroup: frozenset[tuple[int, ...]],
    ) -> sp.Matrix:
        rows = [sp.eye(12).row(index) for index in range(12)]
        for element in subgroup:
            if not any(element):
                continue
            fixed_coordinates = gram.inv() * left_inverse * sp.Matrix(element)
            rows.append(fixed_coordinates.T * target["fixed"])
        denominator = math.lcm(
            *(sp.Rational(value).q for row in rows for value in row)
        )
        integral_generators = sp.Matrix.vstack(
            *(denominator * row for row in rows)
        ).applyfunc(int)
        basis = hermite_normal_form(integral_generators.T).T / denominator
        return sp.simplify(basis * target["gram"] * basis.T)

    histogram: Counter[
        tuple[int, int, tuple[tuple[str, int], ...], bool]
    ] = Counter()
    serialized = []
    for subgroup in sorted(subgroups, key=lambda value: (len(value), sorted(value))):
        ambient = extended_ambient(subgroup)
        reference_id = "0.4.4.2" if len(subgroup) == 2 else "0.4.4.1"
        reference = sp.Matrix(NS_GRAMS[reference_id])
        same_discriminant_form = _binary_discriminant_isometry_exists(
            ambient, reference
        )
        key = (
            len(subgroup),
            abs(int(ambient.det())),
            tuple(sorted(_q_histogram(ambient).items())),
            same_discriminant_form,
        )
        histogram[key] += 1
        serialized.append(
            {
                "order": len(subgroup),
                "elements_in_fixed_smith_coordinates": [
                    list(element) for element in sorted(subgroup)
                ],
                "even_hyperbolic": (
                    exact_inertia(ambient) == (1, 11, 0)
                    and all(ambient[index, index] % 2 == 0 for index in range(12))
                ),
                "reference_discriminant_form_isometric": same_discriminant_form,
                "reference_lattice": (
                    "U+D8(-1)+A1(-1)^2"
                    if len(subgroup) == 2
                    else "U+E8(-1)+A1(-1)^2"
                ),
            }
        )

    counts = Counter(len(subgroup) for subgroup in subgroups)
    checks = {
        "fixed_discriminant_group_order": math.prod(moduli) == 64,
        "order_two_count": counts[2] == 7,
        "order_four_count": counts[4] == 3,
        "every_order_two_extension_has_row_0_4_4_2_form": all(
            determinant == 16
            and dict(q_values) == _q_histogram(sp.Matrix(NS_GRAMS["0.4.4.2"]))
            and same_discriminant_form
            for (order, determinant, q_values, same_discriminant_form), count in histogram.items()
            if order == 2 and count
        ),
        "every_order_four_extension_has_row_0_4_4_1_form": all(
            determinant == 4
            and dict(q_values) == _q_histogram(sp.Matrix(NS_GRAMS["0.4.4.1"]))
            and same_discriminant_form
            for (order, determinant, q_values, same_discriminant_form), count in histogram.items()
            if order == 4 and count
        ),
        "all_extensions_even_hyperbolic": all(
            row["even_hyperbolic"] for row in serialized
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"the fixed-overextension census changed: {checks}")
    subgroup_digest = hashlib.sha256(
        json.dumps(serialized, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "criterion": (
            "an even isotropic subgroup H of A_F, orthogonal to the fixed "
            "eigen-glue class gamma, enlarges F and hence the full action pair; "
            "the ambient divisor-sublattice index is |H|"
        ),
        "fixed_discriminant_moduli": moduli,
        "fixed_discriminant_group_order": 64,
        "nonzero_eigen_glue_fixed_class": list(gamma),
        "admissible_order_two_subgroups": counts[2],
        "admissible_order_four_subgroups": counts[4],
        "all_order_two_extensions_land_in": "U+D8(-1)+A1(-1)^2",
        "all_order_four_extensions_land_in": "U+E8(-1)+A1(-1)^2",
        "landing_identification": (
            "exact discriminant-form isometry plus the indefinite rank-at-least-"
            "length-plus-two uniqueness theorem"
        ),
        "subgroup_list_sha256": subgroup_digest,
        "checks": checks,
        "all_checks_pass": True,
    }


def replay_brandhorst_hofmann_zip(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    if sha256_bytes(payload) != BH_ZIP_SHA256:
        raise AssertionError("the pinned Brandhorst--Hofmann ZIP changed")
    rows = []
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        raw = archive.read(BH_MEMBER)
        if sha256_bytes(raw) != BH_MEMBER_SHA256:
            raise AssertionError("the pinned pure order-four member changed")
        lines = raw.decode("utf-8").splitlines()
        if len(lines) != BH_ROW_COUNT:
            raise AssertionError("the pure order-four row count changed")
        for identifier, specification in BH_ROWS.items():
            line = lines[specification["one_based"] - 1]
            if (
                sha256_bytes((line + "\n").encode())
                != specification["sha256_with_newline"]
            ):
                raise AssertionError(f"Brandhorst--Hofmann row {identifier} changed")
            parsed = _parse_bh_row(line)
            replay = _derive_external_row(line, 4)
            if not (
                parsed["id"] == identifier
                and replay["ns_gram"] == sp.Matrix(NS_GRAMS[identifier])
                and replay["ns_generator"] == sp.Matrix(PURE_ACTION)
            ):
                raise AssertionError(f"row {identifier} failed exact K3 replay")
            rows.append(
                {
                    "id": identifier,
                    "row_one_based": specification["one_based"],
                    "generic_NS_gram_sha256": matrix_digest(replay["ns_gram"]),
                    "generic_transcendental_rank": replay[
                        "transcendental_gram"
                    ].rows,
                    "full_K3_lattice_replayed": True,
                }
            )
    return {
        "zip_sha256": BH_ZIP_SHA256,
        "member": BH_MEMBER,
        "member_sha256": BH_MEMBER_SHA256,
        "row_count": BH_ROW_COUNT,
        "rows": rows,
    }


def build_certificate(bh_zip: Path | None = None) -> dict[str, Any]:
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise AssertionError("the pinned H10/45 source changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(row for row in source["cases"] if row["id"] == TARGET_CASE)
    vinberg = case["vinberg_finite_volume_audit"]
    target = _target_pair(source)

    wall_rows = sp.Matrix(source["coordinate_reconstruction"]["wall_coordinate_rows"])
    wall_gram = wall_rows * target["gram"] * wall_rows.T
    invariant_interior = sp.Matrix([[3, 3, 3, 0, 0, 2, 2, 6, 6, 3, 1, 7]])
    invariant_interior_pairings = list(
        invariant_interior * target["gram"] * wall_rows.T
    )
    wall_orbits = case["subgroup"]["wall_orbits_one_based"]
    permutation = list(range(59))
    for orbit in wall_orbits:
        if len(orbit) == 2:
            left, right = orbit[0] - 1, orbit[1] - 1
            permutation[left], permutation[right] = right, left
    if not (
        source["all_checks_pass"]
        and case["all_checks_pass"]
        and wall_rows.rows == 59
        and all(wall_gram[index, index] == -2 for index in range(59))
        and all(
            wall_rows.row(index) * target["action"]
            == wall_rows.row(permutation[index])
            for index in range(59)
        )
        and invariant_interior * target["action"] == invariant_interior
        and (invariant_interior * target["gram"] * invariant_interior.T)[0]
        == 40
        and min(invariant_interior_pairings) == 1
        and all(value > 0 for value in invariant_interior_pairings)
        and vinberg["all_checks_pass"]
        and vinberg["dimension"] == 10
        and vinberg["facet_count"] == 45
        and vinberg["elliptic_edge_count"] == 7332
        and vinberg["finite_vertex_count"] == 1048
        and vinberg["ideal_vertex_count"] == 93
        and vinberg["bad_edge_count"] == 0
        and vinberg["bad_parabolic_neighborhood_count"] == 0
    ):
        raise AssertionError("the frozen H10/45 chamber or Vinberg audit changed")

    carriers = [_carrier_record(specification, target) for specification in CARRIERS]
    census = _overextension_census(target)
    external = replay_brandhorst_hofmann_zip(bh_zip) if bh_zip else None
    certificate: dict[str, Any] = {
        "schema": "coxeter.r12_h10_45.bh_divisor_sublattice_carrier.v1",
        "artifact_id": "R12-H10-45-BH-DIVISOR-SUBLATTICE-CARRIER",
        "theorem": (
            "There are explicit complex projective K3 surfaces admitting purely "
            "non-symplectic automorphisms tau of order four whose divisor "
            "lattices contain a tau-stable index-two or index-four copy of "
            "U+D4(-1)^2+A1(-1)^2.  The restricted Picard involution is "
            "integrally conjugate to the archived H10/45 action.  Its 59-root "
            "divisor-subsystem chamber therefore has the certified finite-volume "
            "noncompact 45-facet fixed Coxeter section in H10.  The displayed "
            "Weierstrass automorphism sigma supplies the Hodge action; tau is "
            "its Weyl conjugate furnished on the same surface by strong Torelli."
        ),
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": EXPECTED_SOURCE_SHA256,
            "target_case": TARGET_CASE,
        },
        "archived_target_pair": {
            "lattice": "U+D4(-1)^2+A1(-1)^2",
            "rank_signature_determinant": [12, [1, 11], -64],
            "smith": TARGET_SMITH,
            "discriminant_q_histogram": TARGET_Q,
            "action_characteristic": "(lambda-1)^11*(lambda+1)",
            "fixed_rank_signature_determinant": [11, [1, 10], 64],
            "fixed_smith": TARGET_FIXED_SMITH,
            "fixed_discriminant_q_histogram": TARGET_FIXED_Q,
            "fixed_discriminant_q2_histogram": {
                "0": 8,
                "1": 8,
                "2": 8,
                "3": 8,
            },
            "coinvariant_lattice": "<-4>",
            "eigenstack_index": len(target["glue"]),
            "wall_root_count": 59,
            "all_wall_roots_have_square_minus_two": True,
            "all_walls_permuted": True,
            "invariant_integral_interior_vector": [
                int(value) for value in invariant_interior
            ],
            "invariant_interior_square": 40,
            "invariant_interior_wall_pairing_minimum": 1,
            "invariant_interior_wall_pairing_maximum": int(
                max(invariant_interior_pairings)
            ),
        },
        "carriers": carriers,
        "finite_index_carrier_criterion": census,
        "algebro_geometric_effectivity": {
            "actual_K3_automorphisms": True,
            "explicit_Weierstrass_families": True,
            "argument": [
                "Transport the archived invariant 59-root chamber through the exact integral action-pair conjugacy.",
                "Its archived integral interior vector is already Picard-invariant, has square 40, and pairs positively with all 59 roots; perturb it inside the fixed cone away from the locally finite extra-root hyperplanes.",
                "The resulting ambient Neron--Severi Weyl chamber C is invariant because the ambient Picard coinvariant is the rootless lattice <-4>.",
                "Choose w in the Neron--Severi Weyl group carrying C to the actual ample chamber.  Products of square-minus-two reflections act trivially on the discriminant group, so w extends as the identity on the transcendental lattice.",
                "The conjugated order-four Hodge isometry w*sigma^*w^{-1} preserves the ample chamber; strong Torelli gives an automorphism tau of the same explicit K3 surface.",
                "The 59 w-transported roots pair positively with an ample class, so Riemann--Roch makes all of them effective simultaneously.",
            ],
            "all_59_transported_roots_are_algebraic_divisor_classes": True,
            "all_59_effective_signs_can_be_chosen_simultaneously": True,
            "invariant_ample_class_exists": True,
            "archived_invariant_interior_vector_has_square_40": True,
            "archived_invariant_interior_pairs_positively_with_all_59_roots": True,
            "ambient_Picard_coinvariant_is_rootless_minus_four": True,
            "Weyl_conjugator_acts_trivially_on_NS_discriminant": True,
            "Weyl_conjugator_extends_as_identity_on_transcendental_lattice": True,
            "strong_Torelli_produces_effective_chamber_automorphism": True,
            "formula_sigma_need_not_equal_effective_tau": True,
            "scope": (
                "effective invariant divisor-subsystem chamber; extra roots of "
                "the ambient NS can refine it, so it is not claimed to be the "
                "full nef chamber of the carrier K3"
            ),
        },
        "vinberg_handoff": {
            "dimension": 10,
            "facets": 45,
            "classification": "finite_volume_noncompact",
            "finite_volume": True,
            "compact": False,
            "elliptic_edges": 7332,
            "finite_vertices": 1048,
            "cusps": 93,
            "connected_affine_components_audited": vinberg[
                "connected_affine_component_count"
            ],
            "bad_edge_endpoints": 0,
            "bad_affine_zero_neighborhoods": 0,
            "criterion": vinberg["criterion"],
        },
        "brandhorst_hofmann": {
            "article": "https://doi.org/10.1017/fms.2023.50",
            "database": "https://doi.org/10.5281/zenodo.10410278",
            "archive_version": "2.1",
            "zip_sha256": BH_ZIP_SHA256,
            "member": BH_MEMBER,
            "member_sha256": BH_MEMBER_SHA256,
            "external_full_replay_performed": external is not None,
            "external_replay": external,
        },
        "relation_to_previous_no_go": (
            "The primitive E7(-2)/Nikulin obstruction remains valid.  These "
            "carriers escape it because the H10/45 lattice is a nonprimitive "
            "finite-index divisor sublattice of the full Picard lattice and the "
            "geometric automorphism has order four."
        ),
        "novelty_scope": (
            "new project-level explicit algebraic realizations of the already "
            "certified H10/45 polytope; no worldwide-new abstract Coxeter "
            "polytope or arithmetic commensurability class is claimed"
        ),
        "all_checks_pass": all(row["all_checks_pass"] for row in carriers)
        and census["all_checks_pass"],
    }
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(payload).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--brandhorst-hofmann-zip", type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(args.brandhorst_hofmann_zip)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise AssertionError("the archived H10/45 divisor carrier changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "carriers": len(certificate["carriers"]),
                "all_checks_pass": certificate["all_checks_pass"],
                "payload_sha256": certificate["certificate_payload_sha256"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
