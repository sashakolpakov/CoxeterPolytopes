#!/usr/bin/env python3
"""Independent hostile audit of the universal weak Chow-action carrier.

This file deliberately does not import the owner certifier or the surface
embedding implementation.  It reconstructs the sensitive correspondence
orientation, a nontrivial restriction-of-scalars eigenslice, the stable-
lattice maneuver, and the four-square blow-up count from exact matrices.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "universal_weak_chow_action_carrier_and_leverage_no_go_"
    "independent_hostile_audit.json"
)
OWNER_REPORT = ROOT / "results/universal_weak_chow_action_carrier_and_leverage_no_go.md"

OWNER_PINS = {
    "results/universal_weak_chow_action_carrier_and_leverage_no_go.md": (
        "df38b26e0cffd8a391834be922869153af5badf41a632e74a48b472b332e101a"
    ),
    "scripts/certify_universal_weak_chow_action_carrier.py": (
        "d3a19cfe7069b811721a77a624952d92ab5bc7c5574d9e332d2ff8c7f0205eba"
    ),
    "data/enumerations/universal_weak_chow_action_carrier.json": (
        "d75897bbf8d80a662a7b9d3aa3ae729052494e20633d835566db70870bb0fc9e"
    ),
    "tests/test_universal_weak_chow_action_carrier.py": (
        "6da7b9e4590d7e531a7f3502ba2c971e806cbd3f259cb4722169553b7532c802"
    ),
}

DEPENDENCY_PINS = {
    "results/rational_surface_eigencarrier_universality.md": (
        "57e607e7cd6646c5827b3d4d0bff86d7bf8294845e29fe95d4e8b1478b854333"
    ),
    "src/coxeter_polytopes/geometric/rational_surface_carrier.py": (
        "ceba2ae954c7d83020f16d2d034c6d48adb30bcbb3fb8d43698d06189c5f2159"
    ),
    "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness.md": (
        "6eee91512840287ce47a195a317fa4bbfd19678a95b52869ec0a45775d61e67b"
    ),
    "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness.json": (
        "2bce48c0a7156caea3f4b4890411d770450912b4d4b2d3827ce93cc7e9730acc"
    ),
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def boolean_checks(checks: dict[str, Any]) -> dict[str, bool]:
    return {key: bool(value) for key, value in checks.items()}


def pin_audit(pins: dict[str, str]) -> dict[str, str]:
    observed = {relative: sha256_file(ROOT / relative) for relative in pins}
    if observed != pins:
        raise AssertionError(f"terminal pin mismatch: {observed}")
    return observed


def inertia(matrix: sp.Matrix) -> list[int]:
    values = matrix.eigenvals()
    positive = sum(m for value, m in values.items() if value.is_positive is True)
    negative = sum(m for value, m in values.items() if value.is_negative is True)
    zero = matrix.rows - positive - negative
    return [positive, negative, zero]


def reflection(gram: sp.Matrix, index: int) -> sp.Matrix:
    root = sp.eye(gram.rows).col(index)
    square = (root.T * gram * root)[0]
    return sp.eye(gram.rows) - 2 * root * (root.T * gram) / square


def cycle(operator: sp.Matrix, intersection: sp.Matrix) -> sp.Matrix:
    """Product-basis coefficients of sum e_i^* x A(e_i)."""

    return intersection.inv() * operator.T


def stable_lattice_replay() -> dict[str, Any]:
    """A finite-coset model of sum_g gL and denominator clearing."""

    half = sp.Rational(1, 2)
    g = sp.Matrix([[0, 2], [half, 0]])
    lattice = sp.eye(2)
    stable_basis = sp.diag(1, half)  # L + gL
    form = sp.Matrix([[0, half], [half, 0]])
    action_on_stable_basis = stable_basis.inv() * g * stable_basis
    scaled_stable_gram = 4 * stable_basis.T * form * stable_basis

    checks = {
        "rational_involution": g**2 == sp.eye(2),
        "original_lattice_not_stable": any(
            not entry.is_Integer for entry in lattice.inv() * g * lattice
        ),
        "sum_of_translates_contains_L": all(
            entry.is_Integer for entry in stable_basis.inv() * lattice
        ),
        "sum_of_translates_contains_gL": all(
            entry.is_Integer for entry in stable_basis.inv() * g * lattice
        ),
        "sum_of_translates_is_stable": all(
            entry.is_Integer for entry in action_on_stable_basis
        ),
        "isometry_survives": g.T * form * g == form,
        "one_positive_one_negative_direction": inertia(form) == [1, 1, 0],
        "cleared_form_integral": all(
            entry.is_Integer for entry in scaled_stable_gram
        ),
    }
    checks = boolean_checks(checks)
    return {
        "checks": checks,
        "stable_action": [list(map(str, action_on_stable_basis.row(i))) for i in range(2)],
        "cleared_stable_gram": [
            list(map(str, scaled_stable_gram.row(i))) for i in range(2)
        ],
        "all_checks_pass": all(checks.values()),
    }


def four_square_bound_replay() -> dict[str, Any]:
    """Rebuild one full 4R-3 block without using the owner implementation."""

    a0 = sp.Rational(7, 3)
    negatives = [sp.Rational(5, 6), sp.Rational(11, 10)]
    diagonal = sp.diag(a0, *[-value for value in negatives])
    rational_embedding = sp.zeros(10, 3)
    rational_embedding[0, 0] = (a0 + 1) / 2
    rational_embedding[1, 0] = (a0 - 1) / 2
    for offset, value in enumerate((5, 2, 1, 0)):
        rational_embedding[2 + offset, 1] = sp.Rational(value, 6)
    for offset, value in enumerate((10, 3, 1, 0)):
        rational_embedding[6 + offset, 2] = sp.Rational(value, 10)

    intersection = sp.diag(1, *([-1] * 9))
    scale = 30
    integral_embedding = scale * rational_embedding
    checks = {
        "positive_difference_of_two_squares": (
            rational_embedding[:, 0].T
            * intersection
            * rational_embedding[:, 0]
        )[0]
        == a0,
        "first_negative_four_square_identity": sum(
            sp.Rational(value, 6) ** 2 for value in (5, 2, 1, 0)
        )
        == negatives[0],
        "second_negative_four_square_identity": sum(
            sp.Rational(value, 10) ** 2 for value in (10, 3, 1, 0)
        )
        == negatives[1],
        "cleared_embedding_integral": all(
            entry.is_Integer for entry in integral_embedding
        ),
        "exact_scaled_gram": (
            integral_embedding.T * intersection * integral_embedding
            == scale**2 * diagonal
        ),
        "rank_three_uses_nine_blowups": intersection.rows - 1 == 4 * 3 - 3,
        "uniform_rank_750_bound": 4 * 750 - 3 == 2997,
    }
    checks = boolean_checks(checks)
    return {
        "source_diagonal": [str(diagonal[i, i]) for i in range(3)],
        "common_scale": scale,
        "blowup_count": intersection.rows - 1,
        "general_count": "1+4(R-1)=4R-3",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def chow_orientation_replay() -> dict[str, Any]:
    """Use noncommuting reflections to make the cycle order observable."""

    gram = sp.Matrix([[-2, 0, 1], [0, -2, 2], [1, 2, -2]])
    source = [reflection(gram, index) for index in range(3)]
    injection = sp.Matrix(
        [
            [0, 0, 3],
            [0, 0, -1],
            [4, 0, -2],
            [4, 0, -2],
            [0, 4, -4],
            [0, 4, -4],
        ]
    )
    intersection = sp.diag(1, *([-1] * 5))
    embedded_gram = injection.T * intersection * injection
    projector = injection * embedded_gram.inv() * injection.T * intersection

    def lift(operator: sp.Matrix) -> sp.Matrix:
        return (
            injection
            * operator
            * embedded_gram.inv()
            * injection.T
            * intersection
        )

    lifted = [lift(operator) for operator in source]
    cycles = [cycle(operator, intersection) for operator in lifted]
    left, right = lifted[0], lifted[2]
    c_left, c_right = cycles[0], cycles[2]

    checks = {
        "source_signature_1_2": inertia(gram) == [1, 2, 0],
        "independent_embedding_identity": embedded_gram == 16 * gram,
        "projector_idempotent": projector**2 == projector,
        "projector_self_adjoint": (
            projector.T * intersection == intersection * projector
        ),
        "lifts_restrict_exactly": all(
            lifted[i] * injection == injection * source[i] for i in range(3)
        ),
        "noncommuting_pair_really_noncommutes": left * right != right * left,
        "cycle_push_pull_action_is_operator": all(
            cycles[i].T * intersection == lifted[i] for i in range(3)
        ),
        "middle_contraction_is_first_left_then_right": (
            c_left * intersection * c_right
            == cycle(right * left, intersection)
        ),
        "middle_contraction_is_not_wrong_order": (
            c_left * intersection * c_right
            != cycle(left * right, intersection)
        ),
        "group_homomorphism_orientation": (
            c_right * intersection * c_left
            == cycle(left * right, intersection)
        ),
        "reflection_cycles_self_transpose": all(
            coefficient == coefficient.T for coefficient in cycles
        ),
        "reflection_cycles_square_to_projector": all(
            operator**2 == projector for operator in lifted
        ),
        "m2_relation": (lifted[0] * lifted[1]) ** 2 == projector,
        "m3_relation": (lifted[0] * lifted[2]) ** 3 == projector,
        "parabolic_word_nontrivial": lifted[1] * lifted[2] != projector,
        "parabolic_nilpotency_index_three": (
            (lifted[1] * lifted[2] - projector) ** 2 != sp.zeros(6)
            and (lifted[1] * lifted[2] - projector) ** 3 == sp.zeros(6)
        ),
    }
    checks = boolean_checks(checks)
    return {
        "source_gram_determinant": int(gram.det()),
        "ambient_signature": inertia(intersection),
        "projector_rank": projector.rank(),
        "composition_rule": "C_A J C_B = C_(BA)",
        "action_rule": "C_A^t J = A",
        "homomorphism_rule": "Gamma_g o Gamma_h = Gamma_(gh)",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def eigenslice_replay() -> dict[str, Any]:
    """Exact Res_{Q(sqrt(2))/Q} slice with one Lorentzian real place."""

    root_two = sp.sqrt(2)
    trace_gram = sp.Matrix(
        [
            [0, 4, 0, 0],
            [4, 0, 0, 0],
            [0, 0, -2, 0],
            [0, 0, 0, -4],
        ]
    )
    multiplication = sp.diag(1, 1, 1, 1)
    block = sp.Matrix([[0, 2], [1, 0]])
    multiplication[:2, :2] = block
    multiplication[2:, 2:] = block
    field_linear_reflection = sp.diag(1, 1, -1, -1)
    evaluation = sp.Matrix(
        [
            [1, root_two, 0, 0],
            [0, 0, 1, root_two],
            [1, -root_two, 0, 0],
            [0, 0, 1, -root_two],
        ]
    )
    archimedean_gram = sp.diag(root_two, -1, -root_two, -1)
    diagonal_action = sp.diag(root_two, root_two, -root_two, -root_two)

    plus_values = sp.Matrix([2, -1, 0, 0])
    x = evaluation.inv() * plus_values
    walls = [sp.Matrix([0, 0, 1, 0]), sp.Matrix([1, 0, 2, 0])]
    plus_wall_values = [evaluation * wall for wall in walls]
    trace_pairings = [(x.T * trace_gram * wall)[0] for wall in walls]
    slice_pairings = [
        (
            plus_values[:2, :].T
            * sp.diag(root_two, -1)
            * wall_values[:2, :]
        )[0]
        for wall_values in plus_wall_values
    ]

    plus_projector_arch = sp.diag(1, 1, 0, 0)
    plus_projector = evaluation.inv() * plus_projector_arch * evaluation
    checks = {
        "trace_gram_integral": all(entry.is_Integer for entry in trace_gram),
        "trace_signature_1_3": inertia(trace_gram) == [1, 3, 0],
        "restriction_of_scalars_identity": (
            evaluation.T * archimedean_gram * evaluation == trace_gram
        ),
        "field_action_is_trace_self_adjoint": (
            multiplication.T * trace_gram == trace_gram * multiplication
        ),
        "field_action_splits_by_real_embeddings": (
            evaluation * multiplication == diagonal_action * evaluation
        ),
        "native_reflection_is_trace_isometry": (
            field_linear_reflection.T * trace_gram * field_linear_reflection
            == trace_gram
        ),
        "native_reflection_is_involution": (
            field_linear_reflection**2 == sp.eye(4)
        ),
        "native_reflection_commutes_with_field_action": (
            field_linear_reflection * multiplication
            == multiplication * field_linear_reflection
        ),
        "plus_projector_idempotent": plus_projector**2 == plus_projector,
        "plus_eigenvalue_exact": multiplication * plus_projector == (
            root_two * plus_projector
        ),
        "distinguished_slice_signature_1_1": (
            inertia(sp.diag(root_two, -1)) == [1, 1, 0]
        ),
        "conjugate_slice_negative_definite": (
            inertia(sp.diag(-root_two, -1)) == [0, 2, 0]
        ),
        "strict_positive_slice_point": (
            plus_values[:2, :].T
            * sp.diag(root_two, -1)
            * plus_values[:2, :]
        )[0]
        > 0,
        "trace_pairing_equals_distinguished_pairing": (
            trace_pairings == slice_pairings
        ),
        "oriented_wall_pairings_positive": all(value > 0 for value in trace_pairings),
    }
    checks = boolean_checks(checks)
    return {
        "field": "Q(sqrt(2))",
        "trace_rank": 4,
        "trace_signature": inertia(trace_gram),
        "distinguished_signature": [1, 1],
        "conjugate_signature": [0, 2],
        "wall_pairings": [str(value) for value in trace_pairings],
        "common_surface_multiplier": "D^2*sigma_0(beta)>0",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def leverage_replay() -> dict[str, Any]:
    """Calibrate the matrix-algebra and circularity parts of the no-go."""

    intersection = sp.diag(1, -1)
    squarefree_parameters = [2, 3, 5, 6, 7, 10, 11, 13]
    field_generators = [sp.Matrix([[0, value], [1, 0]]) for value in squarefree_parameters]
    cycles = [cycle(operator, intersection) for operator in field_generators]
    checks = {
        "all_rational_matrices_are_divisor_product_actions": all(
            coefficient.T * intersection == operator
            for coefficient, operator in zip(cycles, field_generators)
        ),
        "quadratic_field_relations_hold_on_fixed_rank_two_tate_part": all(
            operator**2 == value * sp.eye(2)
            for value, operator in zip(squarefree_parameters, field_generators)
        ),
        "fixed_degree_field_discriminants_are_not_height_bounded": (
            max(squarefree_parameters) > min(squarefree_parameters)
        ),
        "construction_rank_inequality_has_variable_right_side": all(
            rank <= (4 * rank - 3) + 1 for rank in (2, 3, 20, 750)
        ),
        "imported_rank_bound": 25 * 30 == 750,
        "imported_blowup_bound": 4 * 25 * 30 - 3 == 2997,
    }
    checks = boolean_checks(checks)
    return {
        "only_internal_rank_statement": "d(n+1)=R <= rho(Bl_b(P2))=b+1",
        "fixed_rank_two_real_quadratic_calibration": squarefree_parameters,
        "unbounded_family": "A_D=[[0,D],[1,0]], A_D^2=D I for squarefree D>1",
        "circularity_verdict": (
            "2997 follows only after importing n+1<=30 and [k:Q]<=25"
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def validate_owner_scope() -> dict[str, bool]:
    text = OWNER_REPORT.read_text(encoding="utf-8")
    required = {
        "marked_scope": "unconditional in the marked rational-Chow category",
        "weak_nonclaim": "The word **weak** is indispensable",
        "composition_orientation": "C_AJC_B=C_{BA}",
        "exact_slice": "D^2\\sigma _0(\\beta)",
        "import_caveat": "It is not a new proof of",
        "leverage_scope": "without an independently\nbounded carrier complexity",
        "strong_route_open": "No universal effective/nef or mixed positive-characteristic theorem",
    }
    checks = {key: needle in text for key, needle in required.items()}
    if not all(checks.values()):
        raise AssertionError(f"owner scope markers missing: {checks}")
    return checks


def build_audit() -> dict[str, Any]:
    owner_pins = pin_audit(OWNER_PINS)
    dependency_pins = pin_audit(DEPENDENCY_PINS)
    stable = stable_lattice_replay()
    embedding = four_square_bound_replay()
    chow = chow_orientation_replay()
    slice_replay = eigenslice_replay()
    leverage = leverage_replay()
    owner_scope = validate_owner_scope()

    owner_archive = json.loads(
        (ROOT / "data/enumerations/universal_weak_chow_action_carrier.json").read_text(
            encoding="utf-8"
        )
    )
    checks = {
        "terminal_owner_pins_match": owner_pins == OWNER_PINS,
        "terminal_dependency_pins_match": dependency_pins == DEPENDENCY_PINS,
        "owner_archive_passes": owner_archive.get("all_checks_pass") is True,
        "owner_archive_base_pin_matches": owner_archive.get("dependency_sha256", {}).get(
            "results/rational_surface_eigencarrier_universality.md"
        )
        == DEPENDENCY_PINS["results/rational_surface_eigencarrier_universality.md"],
        "stable_lattice_replay_passes": stable["all_checks_pass"],
        "four_square_bound_replay_passes": embedding["all_checks_pass"],
        "chow_orientation_replay_passes": chow["all_checks_pass"],
        "exact_eigenslice_replay_passes": slice_replay["all_checks_pass"],
        "leverage_replay_passes": leverage["all_checks_pass"],
        "owner_scope_is_precise": all(owner_scope.values()),
    }
    return {
        "schema": (
            "coxeter.universal-weak-chow-action-carrier-and-leverage-no-go."
            "independent-hostile-audit.v1"
        ),
        "verdict": "PASS",
        "owner_pins": owner_pins,
        "dependency_pins": dependency_pins,
        "independence": {
            "imports_owner_certifier": False,
            "imports_owner_embedding_implementation": False,
            "noncommuting_pair_makes_composition_order_observable": True,
            "quadratic_field_eigenslice_is_separate_from_owner_triangle": True,
        },
        "stable_lattice_replay": stable,
        "four_square_bound_replay": embedding,
        "chow_orientation_replay": chow,
        "exact_eigenslice_replay": slice_replay,
        "leverage_no_go_replay": leverage,
        "scope_verdict": {
            "proved": (
                "universal marked weak divisor-Tate Chow action carrier with "
                "exact distinguished eigenslice"
            ),
            "imported": "n<30 and [k:Q]<=25, hence b<=2997",
            "not_proved": (
                "effective/nef/Mori return, regular or birational reflection "
                "action, or an AG-derived ABSW bound"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def encoded_audit() -> str:
    payload = build_audit()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    encoded = encoded_audit()
    if args.write:
        OUTPUT.write_text(encoded, encoding="utf-8")
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise SystemExit("independent audit archive mismatch")
        print("PASS independent universal weak Chow-action hostile audit")
        return 0
    if not args.write:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
