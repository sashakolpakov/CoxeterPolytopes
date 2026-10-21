"""Certify the universal weak rational-Chow reflection-action addendum.

The general construction is proved in the accompanying report.  This replay
uses the arithmetic (2,3,infinity) triangle.  It verifies the Lorentzian
embedding, orthogonal Chow projector, all three lifted side reflections,
literal divisor-product cycle composition, the finite Coxeter relations, and
the nontrivial parabolic product.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

from coxeter_polytopes.geometric.rational_surface_carrier import (
    embed_lorentzian_lattice,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/universal_weak_chow_action_carrier.json"
BASE_REPORT = ROOT / "results/rational_surface_eigencarrier_universality.md"
SCHEMA = "coxeter-polytopes.universal-weak-chow-action-carrier.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix_rows(matrix: sp.MatrixBase) -> list[list[str]]:
    return [[str(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _reflection(gram: sp.Matrix, index: int) -> sp.Matrix:
    root = sp.eye(gram.rows).col(index)
    square = (root.T * gram * root)[0]
    return sp.eye(gram.rows) - 2 * root * (root.T * gram) / square


def _cycle_coefficients(operator: sp.Matrix, intersection: sp.Matrix) -> sp.Matrix:
    """Coefficients of sum e_i^* x A(e_i)."""

    return intersection.inv() * operator.T


def build_certificate() -> dict[str, object]:
    # Root-basis Gram for the finite-area (2,3,infinity) triangle.  All roots
    # have square -2; pairings 0, 1, and 2 give labels 2, 3, and infinity.
    gram = sp.Matrix([[-2, 0, 1], [0, -2, 2], [1, 2, -2]])
    reflections = [_reflection(gram, i) for i in range(3)]

    eigenvalues = gram.eigenvals()
    positive_count = sum(
        multiplicity
        for value, multiplicity in eigenvalues.items()
        if value.is_positive is True
    )
    negative_count = sum(
        multiplicity
        for value, multiplicity in eigenvalues.items()
        if value.is_negative is True
    )
    source_checks = {
        "gram_symmetric": gram == gram.T,
        "gram_determinant_is_2": gram.det() == 2,
        "gram_signature_is_1_2": positive_count == 1 and negative_count == 2,
        "reflections_integral": all(
            all(entry.is_Integer for entry in reflection)
            for reflection in reflections
        ),
        "reflections_are_isometries": all(
            reflection.T * gram * reflection == gram
            for reflection in reflections
        ),
        "reflections_are_involutions": all(
            reflection**2 == sp.eye(3) for reflection in reflections
        ),
        "m2_relation": (reflections[0] * reflections[1]) ** 2 == sp.eye(3),
        "m3_relation": (reflections[0] * reflections[2]) ** 3 == sp.eye(3),
    }

    parabolic = reflections[1] * reflections[2]
    parabolic_nilpotent = parabolic - sp.eye(3)
    source_checks.update(
        {
            "parallel_product_is_nontrivial": parabolic != sp.eye(3),
            "parallel_product_is_unipotent": (
                sp.expand(parabolic.charpoly().as_expr() - (sp.Symbol("lambda") - 1) ** 3)
                == 0
            ),
            "parallel_product_nilpotency_index_3": (
                parabolic_nilpotent**2 != sp.zeros(3)
                and parabolic_nilpotent**3 == sp.zeros(3)
            ),
        }
    )

    embedded = embed_lorentzian_lattice(gram)
    injection = sp.Matrix(embedded.embedding)
    intersection = sp.Matrix(embedded.ambient_gram)
    embedded_gram = injection.T * intersection * injection
    projector = (
        injection
        * embedded_gram.inv()
        * injection.T
        * intersection
    )

    def lift(operator: sp.Matrix) -> sp.Matrix:
        return (
            injection
            * operator
            * embedded_gram.inv()
            * injection.T
            * intersection
        )

    lifted = [lift(reflection) for reflection in reflections]
    projector_cycle = _cycle_coefficients(projector, intersection)
    lifted_cycles = [_cycle_coefficients(operator, intersection) for operator in lifted]

    all_cycle_compositions = True
    for left_index, left in enumerate(lifted):
        for right_index, right in enumerate(lifted):
            left_cycle = lifted_cycles[left_index]
            right_cycle = lifted_cycles[right_index]
            composed_cycle = _cycle_coefficients(right * left, intersection)
            if left_cycle * intersection * right_cycle != composed_cycle:
                all_cycle_compositions = False

    complement_basis = sp.Matrix.hstack(*(injection.T * intersection).nullspace())
    complement_gram = complement_basis.T * intersection * complement_basis
    complement_negative = complement_gram.is_negative_definite is True

    carrier_checks = {
        "embedding_rank_is_3": injection.rank() == 3,
        "embedding_identity": embedded_gram == embedded.vector_scale**2 * gram,
        "blowup_count_is_9": embedded.blowup_count == 9,
        "projector_rank_is_3": projector.rank() == 3,
        "projector_idempotent": projector**2 == projector,
        "projector_self_adjoint": projector.T * intersection == intersection * projector,
        "projector_cycle_idempotent": (
            projector_cycle * intersection * projector_cycle == projector_cycle
        ),
        "orthogonal_complement_negative_definite": complement_negative,
        "lift_restricts_to_native_action": all(
            lifted[i] * injection == injection * reflections[i] for i in range(3)
        ),
        "lifted_reflections_self_adjoint": all(
            operator.T * intersection == intersection * operator for operator in lifted
        ),
        "lifted_reflections_square_to_projector": all(
            operator**2 == projector for operator in lifted
        ),
        "lifted_cycles_self_transpose": all(
            cycle == cycle.T for cycle in lifted_cycles
        ),
        "all_pair_cycle_composition_identities": all_cycle_compositions,
        "lifted_m2_relation": (lifted[0] * lifted[1]) ** 2 == projector,
        "lifted_m3_relation": (lifted[0] * lifted[2]) ** 3 == projector,
    }

    lifted_parabolic = lifted[1] * lifted[2]
    carrier_checks.update(
        {
            "lifted_parallel_product_nontrivial_on_summand": (
                lifted_parabolic != projector
            ),
            "lifted_parallel_product_cubes_nilpotently_about_projector": (
                (lifted_parabolic - projector) ** 2 != sp.zeros(intersection.rows)
                and (lifted_parabolic - projector) ** 3
                == sp.zeros(intersection.rows)
            ),
        }
    )

    native_interior = gram.inv() * sp.ones(3, 1)
    ambient_interior = injection * native_interior
    wall_divisors = [injection.col(i) for i in range(3)]
    wall_pairings = [
        (ambient_interior.T * intersection * divisor)[0]
        for divisor in wall_divisors
    ]
    wall_checks = {
        "transported_wall_gram_exact": (
            sp.Matrix.hstack(*wall_divisors).T
            * intersection
            * sp.Matrix.hstack(*wall_divisors)
            == embedded.vector_scale**2 * gram
        ),
        "strict_interior_positive_square": (
            (ambient_interior.T * intersection * ambient_interior)[0] > 0
        ),
        "strict_interior_all_wall_pairings_positive": all(
            pairing > 0 for pairing in wall_pairings
        ),
        "common_wall_pairing_scale": all(
            pairing == embedded.vector_scale**2 for pairing in wall_pairings
        ),
    }

    published_degree_bound = 25
    published_dimension_plus_one_bound = 30
    trace_rank_bound = published_degree_bound * published_dimension_plus_one_bound
    blowup_bound = 4 * trace_rank_bound - 3
    numerical_checks = {
        "published_degree_bound": published_degree_bound,
        "published_n_plus_one_bound": published_dimension_plus_one_bound,
        "trace_rank_bound": trace_rank_bound,
        "uniform_blowup_bound": blowup_bound,
        "uniform_blowup_bound_is_2997": blowup_bound == 2997,
        "uniform_picard_rank_bound": blowup_bound + 1,
        "construction_only_inequality": "d(n+1)=R <= rho(X)",
        "classical_bounds_imported_not_derived": True,
        "weak_carrier_alone_bounds_dimension_or_field": False,
    }

    checks = {
        key: bool(value)
        for key, value in {**source_checks, **carrier_checks, **wall_checks}.items()
    }
    return {
        "schema": SCHEMA,
        "dependency_sha256": {
            str(BASE_REPORT.relative_to(ROOT)): _sha256(BASE_REPORT),
        },
        "scope": {
            "input": "marked admissible lattice-form arithmetic Coxeter datum",
            "finite_volume_needed_for_carrier": False,
            "reflection_maximality_needed_for_carrier": False,
            "carrier": "divisor Tate summand on a rational blowup with rational Chow action",
            "effectivity_or_nef_return_claimed": False,
            "regular_or_birational_surface_action_claimed": False,
        },
        "triangle_calibration": {
            "name": "arithmetic (2,3,infinity) triangle",
            "root_gram": _matrix_rows(gram),
            "source_reflections": [_matrix_rows(matrix) for matrix in reflections],
            "vector_scale": embedded.vector_scale,
            "surface": f"Bl_{embedded.blowup_count}(P^2)",
            "ambient_rank": intersection.rows,
            "projector_rank": projector.rank(),
            "checks": checks,
            "all_checks_pass": all(checks.values()),
        },
        "unconditional_published_bound_corollary": numerical_checks,
        "all_checks_pass": all(checks.values())
        and numerical_checks["uniform_blowup_bound_is_2997"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_certificate()
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"

    if args.write:
        ARCHIVE.write_text(serialized)
    if args.check:
        if not ARCHIVE.exists():
            raise SystemExit(f"missing archive: {ARCHIVE}")
        if ARCHIVE.read_text() != serialized:
            raise SystemExit("archive mismatch")
        if not result["all_checks_pass"]:
            raise SystemExit("certificate checks failed")
        print("PASS universal weak Chow-action carrier certificate")
        return
    if not args.write:
        print(serialized, end="")


if __name__ == "__main__":
    main()
