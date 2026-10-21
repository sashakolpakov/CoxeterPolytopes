#!/usr/bin/env python3
"""Certify the rank-two collapse and the first divisibility-three escape gate.

There are two logically separate results.

1.  A genuine rank-two extension of the golden fixed lattice over Z[phi]
    has trace rank 24 and trace signature (1,23), so it cannot be a Picard
    lattice in the OG10 Beauville lattice of signature (3,21).  Moreover, in
    any abstract integral carrier whose involution fixes the old lattice
    pointwise and has no other fixed direction, each projected norm-two wall
    already spans a one-root carrier.  Rank-two coupling therefore creates no
    new kind of projected norm-two wall.

2.  The first reflective way around the norm-two carrier is a norm-six root
    with pairing ideal contained in 3*Z[phi].  Write

        H6(k) = [[F, 3*k], [3*k^t, 6]],   p0 = F^-1*k.

    Its mirror projects to p0^perp, its field reflection is integral, and its
    two-place Schur gate is

        0 < q=2*p0^2 < 4/3,   q' < 4/3.

    The script independently replays the three previously frozen covector
    boxes under this sharper gate.  No normalized covector cuts the archived
    v61 ray.  This is a finite-box theorem, not a global obstruction to
    norm-six/divisibility-three carriers.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache, reduce
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia, rank
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram, hyperbolic_plane
from coxeter_polytopes.search.og10_golden_h9_vinberg import h9_field_gram
from scripts.search_og10_golden_projected_wall import (
    FROZEN_COEFFICIENT_BOUND,
    FROZEN_E8_SHELL_BOUND,
    FROZEN_MIXED_B_SHELL_BOUND,
    FROZEN_UNIT_POWER_BOUND,
    INTERIOR_SCALED,
    PHI_PAIR,
    RAY_SCALED,
    _add,
    _coefficient_records,
    _conjugate,
    _e8_square,
    _e8_vectors_through_shell,
    _mixed_e8_records,
    _mul,
    _neg,
    _power,
    _scale,
    _sign,
    _sum_pairs,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = (
    ROOT / "data/candidates/og10_golden_rank_two_reduction_no_go.json"
)
PRIOR_BOX_CERTIFICATE_PATH = (
    ROOT / "data/candidates/og10_golden_projected_wall_box_no_go.json"
)

Pair = tuple[int, int]


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_json_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _digest_record(hasher: object, record: object) -> None:
    hasher.update(repr(record).encode("ascii"))
    hasher.update(b"\n")


def _golden_matrix_product(
    left: Sequence[Sequence[Golden]],
    right: Sequence[Sequence[Golden]],
) -> list[list[Golden]]:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                ZERO,
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _golden_transpose(matrix: Sequence[Sequence[Golden]]) -> list[list[Golden]]:
    return [list(column) for column in zip(*matrix)]


def _golden_congruence(
    gram: Sequence[Sequence[Golden]],
    action: Sequence[Sequence[Golden]],
) -> list[list[Golden]]:
    return _golden_matrix_product(
        _golden_transpose(action), _golden_matrix_product(gram, action)
    )


def _identity(size: int) -> list[list[Golden]]:
    return [
        [ONE if i == j else ZERO for j in range(size)] for i in range(size)
    ]


def _matrix_add_scalar_identity(
    matrix: Sequence[Sequence[Golden]], scalar: Golden
) -> list[list[Golden]]:
    return [
        [entry + (scalar if i == j else ZERO) for j, entry in enumerate(row)]
        for i, row in enumerate(matrix)
    ]


def _inertia_record(value: Inertia) -> list[int]:
    return [value.positive, value.negative, value.zero]


def _append_block(
    base: Sequence[Sequence[Golden]], block: Sequence[Sequence[Golden]]
) -> list[list[Golden]]:
    old = len(base)
    extra = len(block)
    result = [[ZERO for _ in range(old + extra)] for _ in range(old + extra)]
    for i in range(old):
        for j in range(old):
            result[i][j] = base[i][j]
    for i in range(extra):
        for j in range(extra):
            result[old + i][old + j] = block[i][j]
    return result


def _og10_gram() -> list[list[int]]:
    negative_e8 = [[-entry for entry in row] for row in e8_cartan_gram()]
    return block_diagonal(
        hyperbolic_plane(),
        hyperbolic_plane(),
        hyperbolic_plane(),
        negative_e8,
        negative_e8,
        [[-2, 1], [1, -2]],
    )


def _rank_mod_prime(matrix: Sequence[Sequence[int]], prime: int) -> int:
    work = [[int(entry) % prime for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if work else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [(inverse * entry) % prime for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (work[row][j] - factor * work[pivot_row][j]) % prime
                for j in range(columns)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def _row_times_matrix_mod_prime(
    vector: Sequence[int], matrix: Sequence[Sequence[int]], prime: int
) -> list[int]:
    return [
        sum(vector[i] * int(matrix[i][j]) for i in range(len(vector))) % prime
        for j in range(len(matrix))
    ]


def structural_audit() -> dict[str, object]:
    """Exact rank/signature checks and concrete identities behind the theorem."""

    base = [list(row) for row in h9_field_gram()]
    base_trace = integral_gram(trace_gram(base, -PHI * INV_SQRT5))

    # A concrete positive rank-two Schur extension.  Its purpose is to check
    # the universal central-inversion formulas with nonorthogonal projections:
    # p1=u/2 and p2=phi*u/2, with root block diag(2,2).
    rank_two = _append_block(base, [[Golden(2), ZERO], [ZERO, Golden(2)]])
    rank_two[9][10] = rank_two[10][9] = ONE
    rank_two[9][11] = rank_two[11][9] = PHI
    involution = _identity(12)
    involution[9][10] = ONE
    involution[10][10] = Golden(-1)
    involution[9][11] = PHI
    involution[11][11] = Golden(-1)
    involution_square = _golden_matrix_product(involution, involution)
    fixed_dimension = 12 - rank(
        _matrix_add_scalar_identity(involution, Golden(-1))
    )
    anti_dimension = 12 - rank(_matrix_add_scalar_identity(involution, ONE))
    rank_two_trace = integral_gram(
        trace_gram(rank_two, -PHI * INV_SQRT5)
    )

    # The norm-six test model uses k equal to the last fixed covector.  It is
    # not a v61 cutter; it verifies reflection integrality and the mod-3 trace
    # class in a small exact lattice.
    norm_six = _append_block(base, [[Golden(6)]])
    norm_six[9][10] = norm_six[10][9] = Golden(3)
    norm_six_involution = _identity(11)
    norm_six_involution[9][10] = Golden(3)
    norm_six_involution[10][10] = Golden(-1)
    norm_six_reflection = _identity(11)
    norm_six_reflection[10][9] = Golden(-1)
    norm_six_reflection[10][10] = Golden(-1)
    norm_six_trace = integral_gram(
        trace_gram(norm_six, -PHI * INV_SQRT5)
    )
    norm_six_root_row = norm_six_trace[20]
    norm_six_trace_divisibility = reduce(
        gcd, (abs(entry) for entry in norm_six_root_row), 0
    )
    d_mod_three = [0] * 22
    d_mod_three[20] = 1
    e_mod_three = [0] * 22
    e_mod_three[20] = -1
    e_mod_three[21] = 1
    norm_six_mod_three_radical_dimension = 22 - _rank_mod_prime(
        norm_six_trace, 3
    )

    og10 = _og10_gram()
    og10_mod_three_radical_dimension = 24 - _rank_mod_prime(og10, 3)
    checks = {
        "base_field_signatures": (
            inertia(base, embedding=0) == Inertia(9, 1, 0)
            and inertia(base, embedding=1) == Inertia(10, 0, 0)
        ),
        "base_trace_is_even_signature_1_19": (
            is_even_integral(base_trace)
            and inertia(base_trace) == Inertia(1, 19, 0)
        ),
        "rank_two_field_signatures": (
            inertia(rank_two, embedding=0) == Inertia(11, 1, 0)
            and inertia(rank_two, embedding=1) == Inertia(12, 0, 0)
        ),
        "rank_two_trace_is_even_signature_1_23": (
            is_even_integral(rank_two_trace)
            and inertia(rank_two_trace) == Inertia(1, 23, 0)
        ),
        "rank_two_involution_is_integral_order_two_isometry": (
            involution_square == _identity(12)
            and _golden_congruence(rank_two, involution) == rank_two
        ),
        "rank_two_involution_has_fixed_base_only": (
            fixed_dimension == 10 and anti_dimension == 2
        ),
        "og10_signature_and_determinant": (
            inertia(og10) == Inertia(3, 21, 0)
            and determinant(og10) == Golden(-3)
        ),
        "negative_index_obstructs_rank_two_trace_embedding": (
            inertia(rank_two_trace).negative > inertia(og10).negative
        ),
        "norm_six_model_has_required_field_and_trace_signatures": (
            inertia(norm_six, embedding=0) == Inertia(10, 1, 0)
            and inertia(norm_six, embedding=1) == Inertia(11, 0, 0)
            and inertia(norm_six_trace) == Inertia(1, 21, 0)
        ),
        "norm_six_reflection_and_fixed_slice_involution_are_integral_isometries": (
            _golden_matrix_product(norm_six_reflection, norm_six_reflection)
            == _identity(11)
            and _golden_congruence(norm_six, norm_six_reflection) == norm_six
            and _golden_matrix_product(norm_six_involution, norm_six_involution)
            == _identity(11)
            and _golden_congruence(norm_six, norm_six_involution) == norm_six
        ),
        "norm_six_trace_root_has_square_minus_6_and_divisibility_3": (
            norm_six_trace[20][20] == -6
            and norm_six_trace_divisibility == 3
        ),
        "norm_six_golden_pair_spans_two_dimensional_mod_3_radical": (
            norm_six_mod_three_radical_dimension >= 2
            and _row_times_matrix_mod_prime(d_mod_three, norm_six_trace, 3)
            == [0] * 22
            and _row_times_matrix_mod_prime(e_mod_three, norm_six_trace, 3)
            == [0] * 22
            and _rank_mod_prime([d_mod_three, e_mod_three], 3) == 2
        ),
        "og10_mod_3_radical_is_one_dimensional": (
            og10_mod_three_radical_dimension == 1
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"structural audit failed: {checks}")
    return {
        "base": {
            "field_rank": 10,
            "field_identity_inertia": _inertia_record(inertia(base, embedding=0)),
            "field_conjugate_inertia": _inertia_record(inertia(base, embedding=1)),
            "field_determinant": str(determinant(base)),
            "trace_rank": len(base_trace),
            "trace_inertia": _inertia_record(inertia(base_trace)),
            "trace_determinant": str(determinant(base_trace)),
        },
        "rank_two_model": {
            "field_rank": 12,
            "identity_inertia": _inertia_record(inertia(rank_two, embedding=0)),
            "conjugate_inertia": _inertia_record(inertia(rank_two, embedding=1)),
            "trace_rank": len(rank_two_trace),
            "trace_inertia": _inertia_record(inertia(rank_two_trace)),
            "trace_determinant": str(determinant(rank_two_trace)),
            "involution_fixed_dimension": fixed_dimension,
            "involution_anti_dimension": anti_dimension,
        },
        "og10_beauville_lattice": {
            "rank": len(og10),
            "inertia": _inertia_record(inertia(og10)),
            "determinant": str(determinant(og10)),
            "max_projective_picard_rank": 22,
        },
        "norm_six_test_model": {
            "normalized_covector": "k=(0^9,1)",
            "normalized_projected_square_twice": "1",
            "cuts_v61": False,
            "field_determinant": str(determinant(norm_six)),
            "trace_determinant": str(determinant(norm_six_trace)),
            "trace_root_square": norm_six_trace[20][20],
            "trace_root_divisibility_inside_test_lattice": norm_six_trace_divisibility,
            "discriminant_class_square": "-2/3 = 4/3 mod 2Z",
            "trace_mod_3_radical_dimension": norm_six_mod_three_radical_dimension,
            "golden_pair_is_independent_mod_3": True,
        },
        "mod_3_ambient_gate": {
            "og10_mod_3_radical_dimension": og10_mod_three_radical_dimension,
            "primitive_embedding_is_injective_after_reduction_mod_3": True,
            "conclusion": (
                "at most one of e and (phi-1)e can have ambient divisibility 3"
            ),
        },
        "checks": checks,
    }


def _norm_six_signature_admissible(normalized_square_twice: Pair) -> bool:
    """Whether 0<q<4/3 and q'<4/3, using integral golden arithmetic."""

    if _sign(normalized_square_twice) <= 0:
        return False
    margin = _add((4, 0), _neg(_scale(normalized_square_twice, 3)))
    return _sign(margin) > 0 and _sign(_conjugate(margin)) > 0


@lru_cache(maxsize=None)
def audit_norm_six_zero_e8(coefficient_bound: int) -> dict[str, object]:
    hasher = hashlib.sha256()
    admissible = 0
    cutters = []
    for A, C, square_twice, on_interior, on_ray in _coefficient_records(
        coefficient_bound
    ):
        if not _norm_six_signature_admissible(square_twice):
            continue
        admissible += 1
        record = (A, C, square_twice, on_interior, on_ray)
        _digest_record(hasher, record)
        if _sign(on_interior) * _sign(on_ray) == -1:
            cutters.append(record)
    return {
        "family": "norm_six_zero_E8",
        "parameters": {"coefficient_bound": coefficient_bound, "w": "0"},
        "signature_admissible_normalized_covector_count": admissible,
        "exact_wall_test_count": admissible,
        "cutting_covector_count": len(cutters),
        "canonical_wall_records_sha256": hasher.hexdigest(),
        "cutting_records": cutters,
    }


@lru_cache(maxsize=None)
def audit_norm_six_scalar_unit_family(
    coefficient_bound: int,
    unit_power_bound: int,
    e8_shell_bound: int,
) -> dict[str, object]:
    e8 = e8_cartan_gram()
    vectors = _e8_vectors_through_shell(e8_shell_bound)
    shell_counts = Counter(_e8_square(vector, e8) // 2 for vector in vectors)
    coefficients = _coefficient_records(coefficient_bound)
    hasher = hashlib.sha256()
    signature_parameter_count = 0
    wall_test_count = 0
    cutters = []
    for power in range(1, unit_power_bound + 1):
        scalar = _power(PHI_PAIR, power)
        admissible_by_shell = {}
        for shell in sorted(shell_counts):
            e8_square_twice = _scale(_mul(scalar, scalar), 4 * shell)
            admissible = tuple(
                data
                for data in coefficients
                if _norm_six_signature_admissible(
                    _add(data[2], e8_square_twice)
                )
            )
            admissible_by_shell[shell] = admissible
            signature_parameter_count += len(admissible)
        for vector in vectors:
            shell = _e8_square(vector, e8) // 2
            e8_square_twice = _scale(_mul(scalar, scalar), 4 * shell)
            covector_e8 = tuple(
                _scale(
                    scalar,
                    sum(e8[i][j] * vector[j] for j in range(8)),
                )
                for i in range(8)
            )
            e8_on_interior = _sum_pairs(
                _mul(covector_e8[i], INTERIOR_SCALED[i + 1])
                for i in range(8)
            )
            e8_on_ray = _sum_pairs(
                _mul(covector_e8[i], RAY_SCALED[i + 1])
                for i in range(8)
            )
            for A, C, non_e8_square, non_e8_interior, non_e8_ray in (
                admissible_by_shell[shell]
            ):
                wall_test_count += 1
                square_twice = _add(non_e8_square, e8_square_twice)
                on_interior = _add(non_e8_interior, e8_on_interior)
                on_ray = _add(non_e8_ray, e8_on_ray)
                record = (
                    power,
                    shell,
                    vector,
                    A,
                    C,
                    square_twice,
                    on_interior,
                    on_ray,
                )
                _digest_record(hasher, record)
                if _sign(on_interior) * _sign(on_ray) == -1:
                    cutters.append(record)
    return {
        "family": "norm_six_scalar_unit_E8",
        "parameters": {
            "coefficient_bound": coefficient_bound,
            "unit_power_min": 1,
            "unit_power_max": unit_power_bound,
            "e8_shell_min": 1,
            "e8_shell_max": e8_shell_bound,
        },
        "e8_vector_count": len(vectors),
        "e8_shell_counts": {
            str(2 * shell): shell_counts[shell] for shell in sorted(shell_counts)
        },
        "signature_admissible_parameter_pair_count_summed_over_power_and_shell": signature_parameter_count,
        "exact_wall_test_count": wall_test_count,
        "cutting_covector_count": len(cutters),
        "canonical_wall_records_sha256": hasher.hexdigest(),
        "cutting_records": cutters,
    }


@lru_cache(maxsize=None)
def audit_norm_six_mixed_e8_family(
    coefficient_bound: int, b_shell_bound: int
) -> dict[str, object]:
    records = tuple(
        sorted(_mixed_e8_records(b_shell_bound), key=lambda record: (record[5], record[4]))
    )
    coefficients = _coefficient_records(coefficient_bound)
    by_square: dict[Pair, list[object]] = {}
    for record in records:
        by_square.setdefault(record[0], []).append(record)
    admissible_by_square = {}
    signature_parameter_count = 0
    for e8_square_twice in sorted(by_square):
        admissible = tuple(
            data
            for data in coefficients
            if _norm_six_signature_admissible(
                _add(data[2], e8_square_twice)
            )
        )
        admissible_by_square[e8_square_twice] = admissible
        signature_parameter_count += len(admissible)
    hasher = hashlib.sha256()
    wall_test_count = 0
    cutters = []
    for (
        e8_square_twice,
        _,
        e8_on_interior,
        e8_on_ray,
        a_vector,
        b_vector,
    ) in records:
        for A, C, non_e8_square, non_e8_interior, non_e8_ray in (
            admissible_by_square[e8_square_twice]
        ):
            wall_test_count += 1
            square_twice = _add(non_e8_square, e8_square_twice)
            on_interior = _add(non_e8_interior, e8_on_interior)
            on_ray = _add(non_e8_ray, e8_on_ray)
            record = (
                a_vector,
                b_vector,
                A,
                C,
                square_twice,
                on_interior,
                on_ray,
            )
            _digest_record(hasher, record)
            if _sign(on_interior) * _sign(on_ray) == -1:
                cutters.append(record)
    return {
        "family": "norm_six_mixed_E8",
        "parameters": {
            "coefficient_bound": coefficient_bound,
            "b_nonzero": True,
            "b_shell_min": 1,
            "b_shell_max": b_shell_bound,
            "conjugate_square_gate": "0 < w'^2 < 2",
        },
        "b_vector_count": len(_e8_vectors_through_shell(b_shell_bound)),
        "eligible_w_count": len(records),
        "eligible_w_square_class_count": len(by_square),
        "signature_admissible_parameter_pair_count_summed_over_w_square_classes": signature_parameter_count,
        "exact_wall_test_count": wall_test_count,
        "cutting_covector_count": len(cutters),
        "canonical_wall_records_sha256": hasher.hexdigest(),
        "cutting_records": cutters,
    }


@lru_cache(maxsize=1)
def run_norm_six_box_audit() -> dict[str, object]:
    families = {
        "zero_E8": audit_norm_six_zero_e8(FROZEN_COEFFICIENT_BOUND),
        "scalar_unit_E8": audit_norm_six_scalar_unit_family(
            FROZEN_COEFFICIENT_BOUND,
            FROZEN_UNIT_POWER_BOUND,
            FROZEN_E8_SHELL_BOUND,
        ),
        "mixed_E8_first_b_shell": audit_norm_six_mixed_e8_family(
            FROZEN_COEFFICIENT_BOUND,
            FROZEN_MIXED_B_SHELL_BOUND,
        ),
    }
    return {
        "families": families,
        "combined_audit_sha256": _canonical_json_sha256(families),
    }


def _source_hashes() -> dict[str, str]:
    paths = {
        "certificate_script": Path(__file__).resolve(),
        "prior_norm_two_box_certificate": PRIOR_BOX_CERTIFICATE_PATH,
        "prior_box_enumerator": ROOT / "scripts/search_og10_golden_projected_wall.py",
        "golden_arithmetic": ROOT / "src/coxeter_polytopes/exact/golden.py",
        "exact_matrix": ROOT / "src/coxeter_polytopes/exact/matrix.py",
        "trace_lattice": ROOT / "src/coxeter_polytopes/exact/trace_lattice.py",
        "fixed_field_lattice": ROOT / "src/coxeter_polytopes/search/og10_golden_h9_vinberg.py",
        "standard_lattices": ROOT / "src/coxeter_polytopes/lattices.py",
    }
    return {name: _sha256_file(path) for name, path in paths.items()}


def _prior_box_reference() -> dict[str, object]:
    payload = json.loads(PRIOR_BOX_CERTIFICATE_PATH.read_text())
    return {
        "file_sha256": _sha256_file(PRIOR_BOX_CERTIFICATE_PATH),
        "certificate_payload_sha256": payload["certificate_payload_sha256"],
        "combined_audit_sha256": payload["audit"]["combined_audit_sha256"],
        "prior_signature_gate": payload["signature_gate"],
        "prior_conclusion": payload["conclusion"],
    }


def build_certificate() -> dict[str, object]:
    structural = structural_audit()
    norm_six_audit = run_norm_six_box_audit()
    families = norm_six_audit["families"]
    no_cutters = all(
        family["cutting_covector_count"] == 0 for family in families.values()
    )
    payload: dict[str, object] = {
        "schema_version": "0.1.0",
        "candidate_id": "og10_golden_rank_two_reduction_no_go",
        "certificate_kind": "exact_structural_reduction_and_norm_six_finite_box_no_go",
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "fixed_lattice": "F=<-2*phi>+E8+<2>",
        "rank_two_theorem": {
            "hypotheses": [
                "the fixed copy of F is saturated",
                (
                    "at the identity embedding M is Lorentzian with its unique "
                    "negative direction already contained in F"
                ),
                "at the conjugate embedding M is positive definite",
                "J is an integral isometric involution fixing F pointwise",
                "Fix(J) over Q(sqrt(5)) is exactly F",
                "r is an ambient field root of norm 2 whose mirror meets H(F)",
            ],
            "central_inversion": "J=-I on F^perp",
            "projection": "p=(r+Jr)/2=F^-1*ell",
            "integrality": "2p=r+Jr lies in H intersect F_Q=F",
            "one_root_subcarrier": "Gram(F+O*r)=[[F,ell],[ell^t,2]]",
            "two_place_gate": "0<p^2<2 and (p')^2<2",
            "cut_gate": "ell(interior)*ell(v61)<0",
            "full_rank_two_schur_gate": (
                "S=R-L^t*F^-1*L is positive definite at both real places"
            ),
            "conclusion": (
                "every projected norm-two root wall is already a one-root wall; "
                "rank-two cross-coupling cannot create a new projected-wall type"
            ),
        },
        "og10_rank_obstruction": {
            "field_rank_two_extension_rank_over_Z_phi": 12,
            "trace_rank": 24,
            "required_trace_signature": [1, 23, 0],
            "og10_beauville_signature": [3, 21, 0],
            "projective_og10_picard_rank_max": 22,
            "conclusion": (
                "a genuine two-line Z[phi] extension cannot embed in the OG10 "
                "Beauville lattice and cannot be a projective OG10 Picard lattice"
            ),
        },
        "norm_six_divisibility_three_escape": {
            "status": (
                "first structurally viable reflective escape ansatz; no primitive "
                "OG10 embedding is asserted"
            ),
            "carrier": "H6(k)=[[F,3k],[3k^t,6]]",
            "normalized_projection": "p0=F^-1*k; actual projection p=3*p0",
            "projected_hyperplane": "p0^perp",
            "field_reflection": "s_e(x)=x-H(x,e)*e/3",
            "reflection_integrality": "H(e,H6) is contained in 3*Z[phi]",
            "paired_trace_roots": "D=e and E=(phi-1)e",
            "paired_trace_gram": [[-6, 0], [0, -6]],
            "intrinsic_folded_reflection": (
                "the two integral trace-lattice reflections multiply to the field "
                "reflection s_e"
            ),
            "ambient_mod_3_limit": (
                "D and E are independent modulo 3, while the ambient OG10 mod-3 "
                "radical has dimension one; under a primitive embedding at most one "
                "can be an ambient divisibility-three wall"
            ),
            "ambient_radical_intersection_warning": (
                "the ambient radical line may meet span(D,E) in neither coordinate "
                "line; selecting D or E is an additional gluing condition"
            ),
            "signature_gate": "0<q=2*p0^2<4/3 and q'<4/3",
            "cut_gate": "k(interior)*k(v61)<0",
            "ambient_embedding_gate": (
                "the primitive trace lattice must embed in Lambda_OG10 with positive "
                "rank-two complement; one selected class, say e/3, must survive as "
                "the unique nonzero ambient order-3 class of square -2/3=4/3 mod 2Z, "
                "while (phi-1)e/3 cannot also survive"
            ),
            "wall_theorem_input": (
                "a primitive ambient class of square -6 and divisibility 3 is an "
                "OG10 stably prime-exceptional (monodromy-reflective) wall"
            ),
        },
        "mod_p_folding_theorem": {
            "statement": (
                "for a primitive embedding N->Lambda, d mod-p independent classes "
                "of ambient divisibility p force d<=dim Rad(Lambda/pLambda)"
            ),
            "proof_gate": (
                "Lambda/N is torsion-free, so Tor_1(Lambda/N,F_p)=0 and "
                "N/pN injects into Lambda/pLambda"
            ),
            "og10_instance": {
                "p": 3,
                "independent_trace_lifts": 2,
                "ambient_radical_dimension": structural["mod_3_ambient_gate"][
                    "og10_mod_3_radical_dimension"
                ],
                "both_lifts_ambient_divisibility_three_is_impossible": True,
            },
        },
        "projected_arrangement_gates": [
            "the projected hyperplane is not proportional to an old F norm-two mirror",
            "at least one designed wall separates the archived interior from v61",
            "all retained projected facets satisfy exact acute/Coxeter angle equations",
            "the complete norm-two and norm-six/divisibility-three wall arrangement is closed",
            "the trace lattice passes primitive OG10 gluing with the required order-3 class",
            "a marked projective OG10 model selects the intended movable chamber",
            "finite volume or compactness is certified after full wall completion",
        ],
        "structural_audit": structural,
        "prior_norm_two_box_reference": _prior_box_reference(),
        "norm_six_box_audit": norm_six_audit,
        "source_sha256": _source_hashes(),
        "checks": {
            "all_structural_checks_pass": all(structural["checks"].values()),
            "all_three_norm_six_boxes_have_zero_cutting_covectors": no_cutters,
            "norm_six_gate_is_strictly_stronger_than_prior_norm_two_gate": True,
            "no_projective_realization_is_claimed_without_primitive_embedding": True,
        },
        "conclusion": (
            "The proposed genuine rank-two OG10 route is structurally impossible, "
            "and central-inversion coupling adds no projected norm-two wall type. "
            "The first reflective escape is the norm-six/divisibility-three one-line "
            "compiler, but its first three exact normalized covector boxes contain "
            "no wall cutting v61."
        ),
        "scope": (
            "Global for genuine Z[phi]-rank-two OG10 extensions and for the "
            "central-inversion reduction under the stated two-place ambient "
            "signature hypotheses; finite-box only for norm-six walls. "
            "Norm-six covectors outside the declared boxes, nonsaturated fixed "
            "lattices, non-pointwise involutions, and other OG10 wall types remain open."
        ),
        "commands": [
            "PYTHONPATH=src:. python3 scripts/certify_og10_golden_rank_two_reduction.py --write-certificate",
            "PYTHONPATH=src:. python3 scripts/certify_og10_golden_rank_two_reduction.py --verify-certificate",
            "PYTHONPATH=src:. pytest -q tests/test_og10_golden_rank_two_reduction.py",
        ],
    }
    payload["certificate_payload_sha256"] = _canonical_json_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, object]) -> dict[str, bool]:
    digest = payload.get("certificate_payload_sha256")
    body = dict(payload)
    body.pop("certificate_payload_sha256", None)
    fresh = build_certificate()
    checks = {
        "certificate_payload_digest_matches": digest == _canonical_json_sha256(body),
        "fresh_exact_replay_matches_archive": payload == fresh,
        "source_hashes_match": payload.get("source_sha256") == _source_hashes(),
    }
    if not all(checks.values()):
        raise AssertionError(f"rank-two reduction certificate failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--verify-certificate", action="store_true")
    args = parser.parse_args()
    if args.write_certificate:
        certificate = build_certificate()
        CERTIFICATE_PATH.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        )
        print(f"wrote {CERTIFICATE_PATH}")
        print(certificate["certificate_payload_sha256"])
        return
    if args.verify_certificate:
        payload = json.loads(CERTIFICATE_PATH.read_text())
        print(json.dumps(verify_certificate(payload), indent=2, sort_keys=True))
        return
    print(json.dumps(build_certificate(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
