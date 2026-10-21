#!/usr/bin/env python3
"""Certify a named rational U+E8+D4 K3 with a cyclic order-12 action.

The surface is the minimal resolution of

    y^2 = x^3 + t^3 (t^2-1)(t^2-4).

Point counts at p=7 prove its geometric Picard lattice is exactly
U+E8(-1)+D4(-1).  The surface carries compatible automorphisms of orders
4, 6, and 12.  Their exact Brandhorst--Hofmann rows and the corresponding
full/triality Coxeter chambers are certified below.

The default replay uses frozen exact point counts and BH row matrices.
``--recount`` independently recomputes the four finite-field counts using
the optional python-flint dependency.  ``--brandhorst-hofmann-zip`` scans
the three complete pinned classification files and replays the target rows.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import itertools
import json
from pathlib import Path
from typing import Any
import zipfile

import sympy as sp

from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    saturated_row_kernel,
    smith_diagonal,
)
from scripts.certify_roulleau_ue8d4_generic_ns_k3_automorphisms import (
    BH_DATA_DOI,
    BH_DATA_VERSION,
    BH_README_MEMBER,
    BH_README_SHA256,
    BH_ZIP_MD5,
    BH_ZIP_SHA256,
    CASES,
    _bh_d4_decomposition,
    derive_bh_row,
    discriminant_audit,
    matrix_digest,
)
from scripts.certify_ue8d4_order6_weierstrass_family import (
    D4,
    D4_TRIALITY,
    E8,
    U,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/ue8d4_named_c12_member.json"
ORDER6_SOURCE = (
    ROOT / "data/enumerations/roulleau_ue8d4_generic_ns_k3_automorphisms.json"
)
STOCK_SOURCE = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"

P = 7
ROOTS_OVER_Q = (-2, -1, 1, 2)
ROOTS_MOD_7 = (5, 6, 1, 2)

NS_GRAM = sp.diag(U, E8, D4)
NS_IDENTITY = sp.eye(14)
NS_TRIALITY = sp.diag(sp.eye(10), D4_TRIALITY)

# Exact BH 0.4.3.1 generic NS Gram.  The action is the identity.  The basis
# below is already block-orthogonal as U, an even unimodular negative rank-8
# block, and a rank-4 block; ORDER4_D4_CHANGE turns the last block into the
# standard D4(-1) Gram.
ORDER4_NS_GRAM = sp.Matrix(
    [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, -2, -1, 1, 1, 1, -1, 1, 1, 0, 0, 0, 0, 0],
        [0, -1, -2, 1, 1, 1, -1, 1, 1, 0, 0, 0, 0, 0],
        [0, 1, 1, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, -2, -1, 1, -1, -1, 0, 0, 0, 0, 0],
        [0, 1, 1, -1, -1, -2, 1, -1, -1, 0, 0, 0, 0, 0],
        [0, -1, -1, 0, 1, 1, -2, 1, 1, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, -1, -1, 1, -2, -1, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, -1, -1, 1, -1, -2, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, -1, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -2, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, -2, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
)
ORDER4_D4_CHANGE = sp.Matrix(
    [[-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, -1], [0, 1, -1, 0]]
)

POINT_COUNTS = (
    {
        "n": 1,
        "q": 7,
        "raw_finite_weierstrass_sum": 56,
        "i0_star_split": False,
        "i0_star_resolution_correction": 7,
        "ii_star_fiber_points": 64,
        "resolved_surface_points": 127,
        "h2_trace": 77,
        "known_ns_normalized_trace": 11,
        "known_ns_trace": 77,
        "quotient_trace": 0,
    },
    {
        "n": 2,
        "q": 49,
        "raw_finite_weierstrass_sum": 2422,
        "i0_star_split": False,
        "i0_star_resolution_correction": 49,
        "ii_star_fiber_points": 442,
        "resolved_surface_points": 2913,
        "h2_trace": 511,
        "known_ns_normalized_trace": 11,
        "known_ns_trace": 539,
        "quotient_trace": -28,
    },
    {
        "n": 3,
        "q": 343,
        "raw_finite_weierstrass_sum": 117992,
        "i0_star_split": True,
        "i0_star_resolution_correction": 1372,
        "ii_star_fiber_points": 3088,
        "resolved_surface_points": 122452,
        "h2_trace": 4802,
        "known_ns_normalized_trace": 14,
        "known_ns_trace": 4802,
        "quotient_trace": 0,
    },
    {
        "n": 4,
        "q": 2401,
        "raw_finite_weierstrass_sum": 5779942,
        "i0_star_split": False,
        "i0_star_resolution_correction": 2401,
        "ii_star_fiber_points": 21610,
        "resolved_surface_points": 5803953,
        "h2_trace": 39151,
        "known_ns_normalized_trace": 11,
        "known_ns_trace": 26411,
        "quotient_trace": 12740,
    },
)

BH_TARGETS: dict[int, dict[str, Any]] = {
    4: {
        "id": "0.4.3.1",
        "member": "K3Groups/results/purely_ns/order04.txt",
        "member_sha256": "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8",
        "member_row_count": 79,
        "row_one_based": 25,
        "row_sha256_with_newline": "a4ad1c37a6864c2d6fc37d1139a99ba3a6f8a889b9a4a5c387ba2420dad00c4b",
        "ns_gram_digest": "f217e1cec1bb4cd615471c7d3f18287fd1b41c690c10207d8347e9c463f74320",
        "ns_action_digest": "e4185c4fdc9dc1c675059566ff40bcb9ab8d2b86c414709914eac5bd258ba48a",
        "transcendental_gram_digest": "d86d496503981b1e11a670c61bcf95ffb1153edc702cde97e5afdd307a27b46c",
        "transcendental_characteristic": "(lambda**2 + 1)**4",
        "full_characteristic": "(lambda - 1)**14*(lambda**2 + 1)**4",
        "unique_rank14_det_minus4_hit": [25, "0.4.3.1", "(lambda - 1)**14"],
    },
    6: {
        "id": "0.6.3.1",
        "member": "K3Groups/results/purely_ns/order06.txt",
        "member_sha256": "45074265c680c21efafe549b03e34cfde10467f10e2a6cd73ef3a9eb19401f63",
        "member_row_count": 150,
        "row_one_based": 57,
        "row_sha256_with_newline": "aacb4c5b6aa9e9f622921803a076c4b0baa0ff40a19d3f4cfd377325bd56b344",
        "ns_gram_digest": "1e97b94f50ac272a897006ffd419af6d2bb2d4b8909e8364041e410f3431782e",
        "ns_action_digest": "2adeb56426c8bddffbfccf5c22a3413d92c6f37fe229e9694818ee280f548cf9",
        "transcendental_gram_digest": "c4c0fdd31752d2907ed9a66b370f7b01f1fc19eced6caaa7d82675c32a657f3b",
        "transcendental_characteristic": "(lambda**2 - lambda + 1)**4",
        "full_characteristic": "(lambda - 1)**12*(lambda**2 - lambda + 1)**4*(lambda**2 + lambda + 1)",
        "unique_rank14_det_minus4_hit": [57, "0.6.3.1", "(lambda - 1)**12*(lambda**2 + lambda + 1)"],
    },
    12: {
        "id": "0.12.1.1",
        "member": "K3Groups/results/purely_ns/order12.txt",
        "member_sha256": "d37bc4bbe2ae99c635cedab379909c0436a720ed06e28cdaa0f8973bf5b48c0b",
        "member_row_count": 48,
        "row_one_based": 13,
        "row_sha256_with_newline": "8e3c52ca3b9809d1b7873adb119b4c34022f9fd9d61f3a43ab7b9e0057dc7bb3",
        "ns_gram_digest": "1e97b94f50ac272a897006ffd419af6d2bb2d4b8909e8364041e410f3431782e",
        "ns_action_digest": "2adeb56426c8bddffbfccf5c22a3413d92c6f37fe229e9694818ee280f548cf9",
        "transcendental_gram_digest": "acd6d3844019f1172db3433e1190b2c913762121aff756277d58d0412d76f5a0",
        "transcendental_characteristic": "(lambda**4 - lambda**2 + 1)**2",
        "full_characteristic": "(lambda - 1)**12*(lambda**2 + lambda + 1)*(lambda**4 - lambda**2 + 1)**2",
        "unique_rank14_det_minus4_hit": [13, "0.12.1.1", "(lambda - 1)**12*(lambda**2 + lambda + 1)"],
    },
}


def matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if sp.Rational(value).q == 1 else str(value) for value in matrix.row(i)]
        for i in range(matrix.rows)
    ]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload_digest(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    encoded = json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _finite_field_recount(n: int) -> dict[str, Any]:
    try:
        from flint import fq_default_ctx
    except ImportError as exc:  # pragma: no cover - dependency-specific
        raise RuntimeError(
            "--recount requires the optional dependency python-flint>=0.9"
        ) from exc

    ctx = fq_default_ctx(P, n, "u", fq_type="FQ_ZECH")
    field = [ctx(list(coefficients)) for coefficients in itertools.product(range(P), repeat=n)]
    q = P**n
    if len(set(field)) != q:
        raise AssertionError("finite-field enumeration is not bijective")
    zero = ctx.zero()
    squares = {value * value for value in field}
    cubes = Counter(value * value * value for value in field)

    def quadratic_character(value: Any) -> int:
        if value == zero:
            return 0
        return 1 if value in squares else -1

    def cubic_character_sum(b_value: Any) -> int:
        return sum(
            multiplicity * quadratic_character(cube + b_value)
            for cube, multiplicity in cubes.items()
        )

    b_histogram: Counter[Any] = Counter()
    one = ctx.one()
    four = ctx(4)
    for t_value in field:
        b_value = t_value**3 * (t_value**2 - one) * (t_value**2 - four)
        b_histogram[b_value] += 1
    cached_sums = {b: cubic_character_sum(b) for b in b_histogram}
    raw = sum(
        multiplicity * (q + 1 + cached_sums[b])
        for b, multiplicity in b_histogram.items()
    )

    # At t=0 the first blow-up is Y^2=t(X^3+4).  The three outer
    # nonidentity arms are indexed by X^3=-4=3.  They are a 3-cycle unless
    # this cubic splits.
    arm_root_count = cubes[ctx(3)]
    if arm_root_count not in (0, 3):
        raise AssertionError("unexpected I0* arm splitting type")
    split = arm_root_count == 3
    i0_correction = (4 if split else 1) * q
    ii_star_points = 9 * q + 1
    total = raw + i0_correction + ii_star_points
    h2_trace = total - 1 - q**2
    normalized_ns_trace = 14 if split else 11
    ns_trace = normalized_ns_trace * q
    return {
        "n": n,
        "q": q,
        "raw_finite_weierstrass_sum": raw,
        "i0_star_split": split,
        "i0_star_resolution_correction": i0_correction,
        "ii_star_fiber_points": ii_star_points,
        "resolved_surface_points": total,
        "h2_trace": h2_trace,
        "known_ns_normalized_trace": normalized_ns_trace,
        "known_ns_trace": ns_trace,
        "quotient_trace": h2_trace - ns_trace,
    }


def point_count_audit(recount: bool = False) -> dict[str, Any]:
    rows = [dict(row) for row in POINT_COUNTS]
    for row in rows:
        q = row["q"]
        split = row["i0_star_split"]
        expected_correction = (4 if split else 1) * q
        expected_ns_factor = 14 if split else 11
        if not (
            q == P ** row["n"]
            and split == (row["n"] % 3 == 0)
            and row["i0_star_resolution_correction"] == expected_correction
            and row["ii_star_fiber_points"] == 9 * q + 1
            and row["resolved_surface_points"]
            == row["raw_finite_weierstrass_sum"]
            + expected_correction
            + 9 * q
            + 1
            and row["h2_trace"] == row["resolved_surface_points"] - 1 - q**2
            and row["known_ns_normalized_trace"] == expected_ns_factor
            and row["known_ns_trace"] == expected_ns_factor * q
            and row["quotient_trace"] == row["h2_trace"] - row["known_ns_trace"]
        ):
            raise AssertionError("a frozen point-count identity failed")
    if recount:
        replayed = [_finite_field_recount(n) for n in range(1, 5)]
        if replayed != rows:
            raise AssertionError("the independent finite-field recount changed")

    traces = [row["quotient_trace"] for row in rows]
    elementary = [sp.Integer(1)]
    for degree in range(1, 5):
        numerator = sum(
            (-1) ** (power - 1)
            * elementary[degree - power]
            * traces[power - 1]
            for power in range(1, degree + 1)
        )
        if numerator % degree:
            raise AssertionError("Newton identity lost integrality")
        elementary.append(sp.Integer(numerator // degree))
    if elementary != [1, 0, 14, 0, -3087]:
        raise AssertionError("the first half of the Frobenius polynomial changed")

    # There is no sign choice here.  Over F_7 the order-six automorphism is
    # defined because mu_3 is contained in F_7.  The quotient is the sum of
    # two four-dimensional primitive-sixth-root eigenspaces.  Cup product
    # pairs them, so if A is Frobenius on one, the other block is
    # q^2*(A^-1)^t.  Their determinant product is q^8.
    coefficients = elementary + [
        P**2 * elementary[3],
        P**4 * elementary[2],
        P**6 * elementary[1],
        P**8,
    ]
    T, x = sp.symbols("T x")
    polynomial = sp.Poly(
        sum((-1) ** degree * coefficients[degree] * T ** (8 - degree) for degree in range(9)),
        T,
    )
    expected = sp.Poly(T**8 + 14 * T**6 - 3087 * T**4 + 33614 * T**2 + P**8, T)
    if polynomial != expected:
        raise AssertionError("the Frobenius polynomial changed")
    normalized_primitive = sp.Poly(sp.expand(polynomial.as_expr().subs(T, P * x) / P**7), x)
    expected_normalized = sp.Poly(7 * x**8 + 2 * x**6 - 9 * x**4 + 2 * x**2 + 7, x)
    if normalized_primitive != expected_normalized:
        raise AssertionError("the normalized Frobenius polynomial changed")

    possible_orders = [m for m in range(1, 31) if sp.totient(m) <= 8]
    if possible_orders != [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 24, 30]:
        raise AssertionError("the degree-at-most-eight cyclotomic list changed")
    cyclotomic_gcd_degrees = {
        str(order): sp.gcd(normalized_primitive, sp.Poly(sp.cyclotomic_poly(order, x), x)).degree()
        for order in possible_orders
    }
    if any(cyclotomic_gcd_degrees.values()):
        raise AssertionError("a normalized Frobenius eigenvalue became cyclotomic")

    u = sp.symbols("u")
    algebraic_integer_obstruction = sp.Poly(7 * u**2 + 2 * u - 23, u)
    if not (
        algebraic_integer_obstruction.is_irreducible
        and sp.discriminant(algebraic_integer_obstruction.as_expr(), u) == 648
    ):
        raise AssertionError("the algebraic-integer obstruction changed")

    return {
        "prime": P,
        "good_reduction": {
            "roots_mod_7": list(ROOTS_MOD_7),
            "pairwise_distinct_and_nonzero": len(set(ROOTS_MOD_7)) == 4 and 0 not in ROOTS_MOD_7,
            "fiber_types_preserved": ["I0*", "II", "II", "II", "II", "II*"],
        },
        "fiber_resolution_formula": {
            "raw_i0_star_cusp_points": "q+1",
            "nonsplit_i0_star_points": "2q+1; correction q",
            "split_i0_star_points": "5q+1; correction 4q",
            "ii_star_points": "9q+1",
            "arm_cubic": "X^3+4=0, or X^3=3 modulo 7",
            "frobenius_arm_action": "3-cycle unless 3 divides n; identity otherwise",
        },
        "counts": rows,
        "quotient_dimension": 8,
        "quotient_power_traces_n1_to_n4": traces,
        "determinant_sign_proof": (
            "Frobenius commutes with zeta over F7; W=V_lambda+V_lambda^-1 with "
            "both dimensions 4, and cup product makes the second block "
            "7^2*(A^-1)^t, so det(F|W)=7^8"
        ),
        "frobenius_polynomial": str(polynomial.as_expr()),
        "normalized_primitive_polynomial": str(normalized_primitive.as_expr()),
        "cyclotomic_orders_exhausted": possible_orders,
        "cyclotomic_gcd_degrees": cyclotomic_gcd_degrees,
        "algebraic_integer_exclusion": {
            "substitution": "u=x^2+x^(-2)",
            "equation": "7*u^2+2*u-23=0",
            "roots": ["(-1+9*sqrt(2))/7", "(-1-9*sqrt(2))/7"],
            "reason": (
                "a root of unity makes u an algebraic integer, but neither root "
                "lies in O_Q(sqrt(2))=Z[sqrt(2)]"
            ),
        },
        "independent_python_flint_recount_available": True,
    }


def standard_action_audit() -> dict[str, Any]:
    if not (
        NS_GRAM.det() == -4
        and exact_inertia(NS_GRAM) == (1, 13, 0)
        and smith_diagonal(NS_GRAM) == [*[1] * 12, 2, 2]
        and NS_TRIALITY * NS_GRAM * NS_TRIALITY.T == NS_GRAM
        and NS_TRIALITY**3 == sp.eye(14)
    ):
        raise AssertionError("the standard U+E8+D4 action changed")

    coordinates = sp.eye(14)
    fixed_d4 = sp.Matrix([[0, 1, 0, 0], [1, 0, 1, 1]])
    coinvariant_d4 = sp.Matrix([[1, 0, -1, 0], [0, 0, 1, -1]])
    fixed_basis = coordinates[list(range(10)), :].col_join(
        fixed_d4 * coordinates[list(range(10, 14)), :]
    )
    coinvariant_basis = coinvariant_d4 * coordinates[list(range(10, 14)), :]
    fixed_gram = fixed_basis * NS_GRAM * fixed_basis.T
    coinvariant_gram = coinvariant_basis * NS_GRAM * coinvariant_basis.T
    glued_basis = fixed_basis.col_join(coinvariant_basis)
    glue_smith = smith_diagonal(glued_basis)
    outer_root_coordinates = coordinates[10, :] * glued_basis.inv()
    if not (
        fixed_gram.det() == -3
        and exact_inertia(fixed_gram) == (1, 11, 0)
        and coinvariant_gram == sp.Matrix([[-4, 2], [2, -4]])
        and abs(int(glued_basis.det())) == 3
        and glue_smith == [*[1] * 13, 3]
        and list(outer_root_coordinates)
        == [*[0] * 10, 0, sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(1, 3)]
    ):
        raise AssertionError("the triality fixed/coinvariant glue changed")

    identity_discriminant = discriminant_audit(NS_GRAM, NS_IDENTITY)
    triality_discriminant = discriminant_audit(NS_GRAM, NS_TRIALITY)
    for audit in (identity_discriminant, triality_discriminant):
        audit.pop("_action")
    if not (
        identity_discriminant["q_histogram_mod_2Z"] == {"0": 1, "1": 3}
        and identity_discriminant["nonzero_action_cycle_lengths"] == [1, 1, 1]
        and triality_discriminant["q_histogram_mod_2Z"] == {"0": 1, "1": 3}
        and triality_discriminant["nonzero_action_cycle_lengths"] == [3]
    ):
        raise AssertionError("the discriminant actions changed")

    return {
        "neron_severi": {
            "class": "U+E8(-1)+D4(-1)",
            "rank_signature_determinant": [14, [1, 13, 0], -4],
            "smith": [*[1] * 12, 2, 2],
            "discriminant_q_histogram_mod_2Z": {"0": 1, "1": 3},
        },
        "identity_action": {
            "fixed_lattice": "all of NS",
            "coinvariant_lattice": "0",
            "discriminant": identity_discriminant,
        },
        "triality_action": {
            "picard_characteristic": "(lambda-1)^12*(lambda^2+lambda+1)",
            "fixed_lattice": {
                "class": "U+E8(-1)+A2(-1)",
                "rank_signature_determinant": [12, [1, 11, 0], -3],
                "gram": matrix_payload(fixed_gram),
            },
            "coinvariant_lattice": {
                "class": "A2(-2)",
                "gram": matrix_payload(coinvariant_gram),
                "rootless": True,
            },
            "fixed_plus_coinvariant_index": 3,
            "quotient_smith": glue_smith,
            "glue_generator": (
                "outer D4 root = (fixed outer-sum + 2*c1 + c2)/3 "
                "modulo the displayed fixed/coinvariant bases"
            ),
            "discriminant": triality_discriminant,
        },
    }


def bh_offline_audit() -> dict[str, Any]:
    order6_case = next(case for case in CASES if case["bh_id"] == "0.6.3.1")
    order6_gram = sp.Matrix(order6_case["ns_gram"])
    order6_action = sp.Matrix(order6_case["ns_action"])
    order6_block = _bh_d4_decomposition(order6_case)

    identity = sp.eye(14)
    coordinate_rows = sp.eye(14)
    order4_block_basis = coordinate_rows[[0, 13, *range(1, 9)], :].col_join(
        ORDER4_D4_CHANGE * coordinate_rows[[9, 10, 11, 12], :]
    )
    order4_block_gram = order4_block_basis * ORDER4_NS_GRAM * order4_block_basis.T
    order4_e8 = order4_block_gram[2:10, 2:10]
    if not (
        abs(int(order4_block_basis.det())) == 1
        and order4_block_gram[:2, :2] == U
        and order4_block_gram[:10, 10:] == sp.zeros(10, 4)
        and order4_block_gram[10:, :10] == sp.zeros(4, 10)
        and order4_block_gram[10:, 10:] == D4
        and exact_inertia(order4_e8) == (0, 8, 0)
        and order4_e8.det() == 1
        and all(order4_e8[index, index] % 2 == 0 for index in range(8))
        and matrix_digest(ORDER4_NS_GRAM) == BH_TARGETS[4]["ns_gram_digest"]
        and matrix_digest(identity) == BH_TARGETS[4]["ns_action_digest"]
        and matrix_digest(order6_gram) == BH_TARGETS[6]["ns_gram_digest"]
        and matrix_digest(order6_action) == BH_TARGETS[6]["ns_action_digest"]
        and BH_TARGETS[12]["ns_gram_digest"] == BH_TARGETS[6]["ns_gram_digest"]
        and BH_TARGETS[12]["ns_action_digest"] == BH_TARGETS[6]["ns_action_digest"]
        and order6_block["d4_action"] == D4_TRIALITY
        and order6_block["block_basis_determinant"] == 1
    ):
        raise AssertionError("a frozen BH integral pair decomposition changed")

    identity_disc = discriminant_audit(ORDER4_NS_GRAM, identity)
    triality_disc = discriminant_audit(order6_gram, order6_action)
    for audit in (identity_disc, triality_disc):
        audit.pop("_action")
    if not (
        identity_disc["nonzero_action_cycle_lengths"] == [1, 1, 1]
        and triality_disc["nonzero_action_cycle_lengths"] == [3]
        and identity_disc["q_histogram_mod_2Z"]
        == triality_disc["q_histogram_mod_2Z"]
        == {"0": 1, "1": 3}
    ):
        raise AssertionError("a BH discriminant action changed")

    rows = []
    for order in (4, 6, 12):
        target = dict(BH_TARGETS[order])
        target["automorphism_order"] = order
        target["generic_ns_rank_signature_determinant"] = [14, [1, 13, 0], -4]
        target["generic_ns_smith"] = [*[1] * 12, 2, 2]
        target["complete_file_uniqueness_criterion"] = (
            "unique row whose generic NS has rank 14 and determinant -4"
        )
        target["fixed_coinvariant_glue_pair"] = (
            "identity/full NS/trivial glue"
            if order == 4
            else "U+E8+A2 fixed, A2(-2) coinvariant, index-3 glue, literal D4 triality"
        )
        rows.append(target)
    return {
        "article": "https://doi.org/10.1017/fms.2023.50",
        "database": f"https://doi.org/{BH_DATA_DOI}",
        "database_version": BH_DATA_VERSION,
        "database_zip_sha256": BH_ZIP_SHA256,
        "rows": rows,
        "order4_exact_block_decomposition": {
            "block_basis_determinant": int(order4_block_basis.det()),
            "blocks": "U + E8(-1) + D4(-1)",
            "d4_change_of_basis": matrix_payload(ORDER4_D4_CHANGE),
            "e8_reason": "even negative-definite unimodular rank 8, hence E8(-1)",
            "action": "identity on every block",
        },
        "order6_order12_exact_pair_equality": {
            "same_ns_gram_entrywise": True,
            "same_ns_action_entrywise": True,
            "block_basis_determinant": order6_block["block_basis_determinant"],
            "blocks": "U+E8(-1) [identity] + D4(-1) [outer 3-cycle]",
            "d4_action": matrix_payload(order6_block["d4_action"]),
        },
        "discriminant_actions": {
            "order4": identity_disc,
            "orders6_and12": triality_disc,
        },
        "full_external_scan_available": True,
    }


def replay_bh_zip(path: Path) -> None:
    payload = path.read_bytes()
    if not (
        hashlib.md5(payload).hexdigest() == BH_ZIP_MD5  # noqa: S324; provenance only
        and hashlib.sha256(payload).hexdigest() == BH_ZIP_SHA256
    ):
        raise AssertionError("the pinned Brandhorst--Hofmann ZIP changed")
    order6_case = next(case for case in CASES if case["bh_id"] == "0.6.3.1")
    expected_triality_gram = sp.Matrix(order6_case["ns_gram"])
    expected_triality_action = sp.Matrix(order6_case["ns_action"])
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        readme = archive.read(BH_README_MEMBER)
        if not (
            hashlib.sha256(readme).hexdigest() == BH_README_SHA256
            and b"each line represents a moduli space parametrizing pairs (X,G)" in readme
        ):
            raise AssertionError("the BH database semantics changed")
        for order in (4, 6, 12):
            target = BH_TARGETS[order]
            member = archive.read(target["member"])
            if hashlib.sha256(member).hexdigest() != target["member_sha256"]:
                raise AssertionError("a complete BH order file changed")
            lines = member.decode().splitlines()
            if len(lines) != target["member_row_count"]:
                raise AssertionError("a complete BH row count changed")
            hits = []
            selected = None
            for row_number, line in enumerate(lines, 1):
                row = derive_bh_row(line, order)
                gram = row["ns_gram"]
                action = row["ns_action"]
                if gram.rows == 14 and gram.det() == -4:
                    hits.append(
                        [row_number, row["id"], str(sp.factor(action.charpoly().as_expr()))]
                    )
                if row["id"] == target["id"]:
                    selected = (row_number, line, row)
            if hits != [target["unique_rank14_det_minus4_hit"]] or selected is None:
                raise AssertionError("a BH uniqueness scan changed")
            row_number, line, row = selected
            expected_gram = ORDER4_NS_GRAM if order == 4 else expected_triality_gram
            expected_action = sp.eye(14) if order == 4 else expected_triality_action
            if not (
                row_number == target["row_one_based"]
                and hashlib.sha256((line + "\n").encode()).hexdigest()
                == target["row_sha256_with_newline"]
                and row["ns_gram"] == expected_gram
                and row["ns_action"] == expected_action
                and matrix_digest(row["transcendental_gram"])
                == target["transcendental_gram_digest"]
                and str(sp.factor(row["transcendental_action"].charpoly().as_expr()))
                == target["transcendental_characteristic"]
                and str(sp.factor(row["k3_action"].charpoly().as_expr()))
                == target["full_characteristic"]
            ):
                raise AssertionError("an exact BH target replay changed")


def build_certificate(recount: bool = False) -> dict[str, Any]:
    point_counts = point_count_audit(recount)
    lattice_actions = standard_action_audit()
    bh = bh_offline_audit()

    order6_source = json.loads(ORDER6_SOURCE.read_text(encoding="utf-8"))
    order6_case = next(
        case for case in order6_source["cases"] if case["brandhorst_hofmann_row"]["id"] == "0.6.3.1"
    )
    hyperbolic = order6_case["hyperbolic_polytope"]
    stock = json.loads(STOCK_SOURCE.read_text(encoding="utf-8"))
    stock_case = next(case for case in stock["stocks"] if case["key"] == "R14_UE8D4")
    if not (
        order6_source["all_checks_pass"]
        and hyperbolic["classification"] == "finite_volume_noncompact"
        and hyperbolic["dimension"] == 11
        and hyperbolic["facets"] == 13
        and hyperbolic["bad_edges"] == 0
        and stock["all_checks_pass"]
        and stock_case["hyperbolic_dimension"] == 13
        and stock_case["wall_count"] == 15
        and stock_case["source_completeness_method"]
        == "direct_curve_enumeration_and_completeness_test"
    ):
        raise AssertionError("an upstream Coxeter handoff changed")

    checks = {
        "global_k3_fiber_euler_sum": 6 + 4 * 2 + 10 == 24,
        "named_roots_distinct": len(set(ROOTS_OVER_Q)) == 4 and 0 not in ROOTS_OVER_Q,
        "point_counts_reconstruct_degree8_weil_polynomial": point_counts[
            "frobenius_polynomial"
        ]
        == "T**8 + 14*T**6 - 3087*T**4 + 33614*T**2 + 5764801",
        "no_cyclotomic_quotient_eigenvalue": not any(
            point_counts["cyclotomic_gcd_degrees"].values()
        ),
        "reduction_picard_rank_is_14_by_tate": True,
        "specialization_bounds_char0_picard_rank_by_14": True,
        "displayed_trivial_lattice_has_rank_14": NS_GRAM.rows == 14,
        "no_same_rank_even_overlattice": (1 - 13) % 8 != 0,
        "exact_char0_ns": lattice_actions["neron_severi"]["class"]
        == "U+E8(-1)+D4(-1)",
        "bh_three_unique_exact_matches": len(bh["rows"]) == 3,
        "vinberg_triality_finite_volume": hyperbolic["classification"]
        == "finite_volume_noncompact",
    }
    if not all(checks.values()):
        raise AssertionError(f"a named-member check failed: {checks}")

    payload: dict[str, Any] = {
        "schema": "ue8d4-named-c12-member-v1",
        "theorem": (
            "The minimal resolution over Qbar of y^2=x^3+t^3(t^2-1)(t^2-4) "
            "has NS=U+E8(-1)+D4(-1).  It carries a cyclic order-12 action whose "
            "order-4, order-6, and order-12 cyclic subactions realize BH rows "
            "0.4.3.1, 0.6.3.1, and 0.12.1.1 and the known full H13/15 and "
            "triality H11/13 Coxeter chambers."
        ),
        "surface": {
            "field": "Q",
            "affine_equation": "y^2 = x^3 + t^3*(t^2-1)*(t^2-4)",
            "family_parameters": {"a": -1, "b": 2, "c": -2},
            "fibers": [
                {"place": "0", "kodaira": "I0*", "root_lattice": "D4(-1)", "euler": 6},
                {"place": "+/-1,+/-2", "kodaira": "four II", "root_lattice": "0", "euler": 8},
                {"place": "infinity", "kodaira": "II*", "root_lattice": "E8(-1)", "euler": 10},
            ],
            "trivial_lattice": "U+E8(-1)+D4(-1)",
        },
        "symmetric_c12_subfamily": {
            "equation": "y^2=x^3+t^3*(t^2-r^2)*(t^2-s^2)",
            "open_condition": "r*s*(r^2-s^2) != 0",
            "moduli_dimension": 1,
            "nonisotrivial_parameter": "the ratio r/s, modulo a finite symmetry",
            "generic_exact_ns_argument": (
                "The displayed divisors give rho>=14.  Primitive order 12 makes "
                "rank(T) divisible by phi(12)=4, so a jump has rho>=18; its "
                "primitive eigenspace has projective period dimension 0, whereas "
                "the symmetric family is genuinely one-dimensional.  Thus very "
                "general rho=14, and the determinant-4 saturation gate gives exact NS."
            ),
        },
        "automorphisms": {
            "relations": {
                "zeta": "(x,y,t)->(omega*x,-y,t), order 6",
                "iota": "(x,y,t)->(-x,i*y,-t), order 4",
                "h=zeta*iota": "(x,y,t)->(-omega*x,-i*y,-t), order 12",
                "common_square_relation": "iota^2=zeta^3=(x,-y,t)",
                "commute": True,
                "cyclic_group": "<zeta,iota>=<h>=C12; zeta=h^10 and iota=h^3",
            },
            "two_form_multipliers": {
                "zeta": "-omega, primitive order 6",
                "iota": "-i, primitive order 4",
                "h": "i*omega, primitive order 12",
            },
            "resolved_fiber_action": {
                "i0_star_blowup_chart": "x=t*X, y=t*Y; singular arm points satisfy X^3+4=0",
                "zeta": "X->omega*X, a 3-cycle on the outer D4 arms",
                "iota": "X->X and Y->-i*Y, fixing every D4 arm",
                "h": "the same outer D4 3-cycle as zeta",
                "ii_star": "the E8 affine diagram is rigid, so every component class is fixed",
                "U": "fiber and zero-section classes are fixed",
                "four_II_fibers": "iota and h swap +1 with -1 and +2 with -2; these add no root classes",
            },
            "picard_actions": {
                "order4_iota": "identity on NS",
                "order6_zeta": "identity on U+E8(-1), outer triality on D4(-1)",
                "order12_h": "identity on U+E8(-1), outer triality on D4(-1)",
            },
        },
        "point_count_and_tate_certificate": point_counts,
        "picard_rank_conclusion": {
            "tate_source": "https://arxiv.org/abs/1301.6326",
            "reduction_geometric_picard_rank": 14,
            "specialization_injection": "NS(X_Qbar) injects into NS(X_F7bar)",
            "characteristic_zero_picard_rank": 14,
            "saturation": (
                "a proper index-2 even enlargement would be even unimodular of "
                "signature (1,13), impossible because 1-13 is not divisible by 8"
            ),
            "exact_ns": "U+E8(-1)+D4(-1)",
        },
        "fixed_coinvariant_glue": lattice_actions,
        "brandhorst_hofmann": bh,
        "coxeter_handoff": {
            "order4_full_chamber": {
                "dimension": 13,
                "facets": 15,
                "status": "known complete U+E8(-1)+D4(-1) K3 nef Coxeter chamber",
                "finite_volume_noncompact": True,
                "source_completeness": stock_case["source_completeness_method"],
            },
            "orders6_and12_triality_chamber": {
                "dimension": hyperbolic["dimension"],
                "facets": hyperbolic["facets"],
                "classification": hyperbolic["classification"],
                "finite_vertices": hyperbolic["finite_vertices"],
                "cusps": hyperbolic["cusps"],
                "elliptic_edges": hyperbolic["elliptic_edges"],
                "bad_edges": hyperbolic["bad_edges"],
                "parabolic_neighborhood_audits": hyperbolic["parabolic_neighborhood_audits"],
            },
            "novelty": (
                "The abstract chambers are known.  The positive result is a named "
                "rational K3, an arithmetic proof of its exact NS, and simultaneous "
                "equation-level C4/C6/C12 realizations."
            ),
        },
        "provenance": {
            "order6_source": str(ORDER6_SOURCE.relative_to(ROOT)),
            "order6_source_sha256": sha256_file(ORDER6_SOURCE),
            "stock_source": str(STOCK_SOURCE.relative_to(ROOT)),
            "stock_source_sha256": sha256_file(STOCK_SOURCE),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["payload_sha256"] = payload_digest(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recount", action="store_true")
    parser.add_argument("--brandhorst-hofmann-zip", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.brandhorst_hofmann_zip is not None:
        replay_bh_zip(args.brandhorst_hofmann_zip)
    payload = build_certificate(recount=args.recount)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise AssertionError(f"archive mismatch: {args.output}")
        print(f"archive OK: {args.output}")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()
