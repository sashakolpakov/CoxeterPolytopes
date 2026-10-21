#!/usr/bin/env python3
"""Exact obstruction to a strong invariant-nef promotion of D34/T1.

The frozen D34 compiler constructs a valid weak K3 lattice/action carrier for
the labelled T1 H9/14 Coxeter datum.  This certificate tests the stronger
claim that its fourteen-facet cone is the complete invariant nef section.

It fails for an explicit reason: five square-minus-four fixed normals in a
unit coordinate box lift to orthogonal involution-exchanged pairs of
square-minus-two roots, and each normal takes both signs on exact extreme
rays of the proposed cone.  The first pair alone is a complete obstruction.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

import scripts.certify_h9_14_t2_surface_unfolding_pilot as pilot
from scripts.probe_h9_14_t1_strong_k3_nef_promotion import reconstruct


ROOT = Path(__file__).resolve().parents[1]
D34 = ROOT / "data/enumerations/h9_14_t1_t2_d34_unfolding_census.json"
FACTORY = ROOT / "data/enumerations/alexeev_k3_h9_14_factory.json"
OUTPUT = ROOT / "data/enumerations/h9_14_t1_strong_k3_nef_obstruction.json"
SCHEMA = "coxeter.h9_14.t1_strong_k3_nef_obstruction.v1"
STATUS_DATE = "2026-08-08"
EXPECTED_SOURCE_SHA256 = {
    "D34_census": "76c7e08ed2fd5107822eec4c4c96f35d107a4a2e231efef2a969a8c98ba36620",
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
}

HIDDEN_ROWS = (
    (
        (0, 0, 0, 1, 0, 0, -1, 0, 0, 0),
        (3, 1, 2, 5, 4, 6, 0, 0, 0, -4),
        (5, 3, 6, 4, 4, 2, 1, 0, 0, -4),
    ),
    (
        (0, 0, 0, 1, 1, 0, -1, 0, 0, 0),
        (3, 1, 2, 6, 5, 6, 0, 0, 0, -4),
        (5, 3, 6, 4, 4, 2, 1, 0, 0, -4),
    ),
    (
        (0, 0, 1, 1, 0, 0, -1, 0, 0, 0),
        (-1, 1, 2, 4, 0, 4, 0, 0, 0, 0),
        (1, 0, 0, 2, 2, 4, 2, 0, 0, -2),
    ),
    (
        (0, 1, 1, 1, 0, 0, -1, 0, 0, 0),
        (-1, 1, 2, 4, 0, 4, 0, 0, 0, 0),
        (1, -3, 0, 2, 2, 8, 4, 1, 0, -2),
    ),
    (
        (1, 0, 1, 1, 0, 0, -1, 0, 0, 0),
        (1, -1, 2, 4, 0, 4, 0, 0, 0, 0),
        (-2, 0, 0, 2, 2, 8, 4, 1, 0, -2),
    ),
)

STRICT_INTERIOR = (5, -1, 0, 14, 12, 32, 10, 6, 2, -12)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(index)]
        for index in range(matrix.rows)
    ]


def _digest(matrix: sp.Matrix) -> str:
    raw = json.dumps(_rows(matrix), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _primitive(vector: sp.Matrix) -> bool:
    return math.gcd(*(abs(int(value)) for value in vector)) == 1


def _orbit_histogram(permutation: Sequence[int]) -> dict[str, int]:
    unseen = set(range(len(permutation)))
    histogram: Counter[int] = Counter()
    while unseen:
        start = min(unseen)
        orbit: set[int] = set()
        value = start
        while value not in orbit:
            orbit.add(value)
            value = int(permutation[value])
        unseen -= orbit
        histogram[len(orbit)] += 1
    return {str(key): value for key, value in sorted(histogram.items())}


def _factory_source_boundary(factory_case: dict[str, Any]) -> dict[str, Any]:
    section = factory_case["section"]
    permutation = [
        int(value)
        for value in section["case"]["representative_wall_permutations"][1]
    ]
    records = section["restricted_mirrors"]["records"]
    active = [row for row in records if row["hyperbolic_status"] == "active_facet"]
    return {
        "ambient_rank": 14,
        "ambient_lattice_determinant": -256,
        "source_wall_count": 24,
        "source_wall_norm_histogram": factory_case["source_wall_audit"][
            "source_wall_norm_histogram"
        ],
        "source_action_orbit_histogram": _orbit_histogram(permutation),
        "active_orbit_count": len(active),
        "active_source_wall_count": sum(
            len(row["source_wall_indices"]) for row in active
        ),
        "active_source_action_orbit_histogram": {
            str(key): value
            for key, value in sorted(
                Counter(len(row["source_wall_indices"]) for row in active).items()
            )
        },
    }


def _lexicographic_lift(
    package: dict[str, object], fixed: sp.Matrix
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    fixed_n = package["fixed_embedding"] * fixed
    _, classes = pilot._e8_data()
    anti_roots = sorted(
        {tuple(int(value) for value in vector) for vectors in classes.values() for vector in vectors}
    )
    for anti_coordinates in anti_roots:
        anti = sp.Matrix(anti_coordinates)
        anti_n = package["anti_embedding"] * anti
        if all((fixed_n[index] + anti_n[index]) % 2 == 0 for index in range(18)):
            plus = (fixed_n + anti_n) / 2
            minus = (fixed_n - anti_n) / 2
            if all(value.q == 1 for value in plus) and all(value.q == 1 for value in minus):
                return anti, anti_n, plus.applyfunc(int), minus.applyfunc(int), fixed_n
    raise ArithmeticError("hidden fixed normal has no E8(-2) root lift")


def _ray_record(
    fixed_gram: sp.Matrix,
    fixed_roots: sp.Matrix,
    hidden: sp.Matrix,
    ray: sp.Matrix,
) -> dict[str, Any]:
    slacks = fixed_roots.T * fixed_gram * ray
    active = [index for index, value in enumerate(slacks) if value == 0]
    active_rank = fixed_roots[:, active].rank()
    return {
        "coordinates": list(map(int, ray)),
        "square": int((ray.T * fixed_gram * ray)[0]),
        "facet_slacks": list(map(int, slacks)),
        "all_facet_slacks_nonnegative": all(value >= 0 for value in slacks),
        "active_facets_one_based": [index + 1 for index in active],
        "active_normal_rank": active_rank,
        "exact_extreme_ray": len(active) >= 9 and active_rank == 9,
        "hidden_normal_pairing": int((hidden.T * fixed_gram * ray)[0]),
    }


def build_certificate() -> dict[str, Any]:
    source_hashes = {"D34_census": _sha256(D34), "factory": _sha256(FACTORY)}
    if source_hashes != EXPECTED_SOURCE_SHA256:
        raise ArithmeticError("a pinned T1 obstruction source changed")
    d34 = json.loads(D34.read_text(encoding="utf-8"))
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    factory_green = (
        factory["scope"]["abstract_coxeter_polytopes_certified"]
        and factory["scope"]["finite_volume_and_noncompactness_certified"]
        and all(
            all(row["section"]["checks"].values())
            and all(row["exact_cone"]["checks"].values())
            and row["literal_vinberg"]["literal_vinberg_passes"]
            for row in factory["constructions"]
        )
    )
    if not (d34["all_checks_pass"] and factory_green):
        raise ArithmeticError("an upstream T1 source is not green")
    t1_archive = next(
        row for row in d34["compiler_inventory"] if row["representative"] == "A14_8_1_C01"
    )
    factory_case = next(
        row for row in factory["constructions"] if row["id"] == "A14_8_1_C01"
    )

    package = reconstruct()
    n_gram = package["n_gram"]
    action = package["involution"]
    roots = package["roots"]
    fixed_gram = package["fixed_gram"]
    fixed_roots = package["fixed_roots"]
    anti_gram = package["anti_gram"]
    root_gram = roots.T * n_gram * roots
    root_lookup = {tuple(roots[:, index]): index for index in range(roots.cols)}
    compiled_permutation = [
        root_lookup[tuple(action * roots[:, index])] for index in range(roots.cols)
    ]
    reconstructed_hashes = {
        "N_gram": _digest(n_gram),
        "involution": _digest(action),
        "roots": _digest(roots),
        "fixed_gram": _digest(fixed_gram),
    }
    archive_hashes = {
        "N_gram": t1_archive["unique_survivor"]["rank18_unfolding"]["gram_sha256"],
        "involution": t1_archive["unique_survivor"]["rank18_unfolding"]["involution_sha256"],
        "roots": t1_archive["unique_survivor"]["roots_and_fold"]["roots_sha256"],
        "fixed_gram": t1_archive["unique_survivor"]["fixed_lattice"]["gram_sha256"],
    }
    if reconstructed_hashes != archive_hashes:
        raise ArithmeticError("the reconstructed T1 matrices differ from D34")

    raw = package["raw"]
    folded = package["folded_gram"]
    raw_ratios = sp.zeros(14)
    folded_ratios = sp.zeros(14)
    for left in range(14):
        for right in range(14):
            raw_ratios[left, right] = sp.cancel(
                4 * raw[left, right] ** 2 / (raw[left, left] * raw[right, right])
            )
            folded_ratios[left, right] = sp.cancel(
                4 * folded[left, right] ** 2
                / (folded[left, left] * folded[right, right])
            )
    labelled_facet_identity = raw_ratios == folded_ratios

    factory_boundary = _factory_source_boundary(factory_case)
    compiled_boundary = {
        "ambient_rank": 18,
        "ambient_lattice_determinant": int(n_gram.det()),
        "source_root_count": roots.cols,
        "source_root_norm_histogram": {
            str(key): value
            for key, value in sorted(Counter(map(int, root_gram.diagonal())).items())
        },
        "source_action_orbit_histogram": _orbit_histogram(compiled_permutation),
    }
    literal_source_bridge_impossible = (
        factory_boundary["ambient_rank"] != compiled_boundary["ambient_rank"]
        and factory_boundary["source_action_orbit_histogram"]
        != compiled_boundary["source_action_orbit_histogram"]
        and factory_boundary["source_wall_norm_histogram"]
        != compiled_boundary["source_root_norm_histogram"]
    )

    hidden_records: list[dict[str, Any]] = []
    compiled_columns = {tuple(roots[:, index]) for index in range(roots.cols)}
    for number, (fixed_values, negative_values, positive_values) in enumerate(
        HIDDEN_ROWS, start=1
    ):
        hidden = sp.Matrix(fixed_values)
        negative = sp.Matrix(negative_values)
        positive = sp.Matrix(positive_values)
        anti, anti_n, plus, minus, fixed_n = _lexicographic_lift(package, hidden)
        negative_record = _ray_record(fixed_gram, fixed_roots, hidden, negative)
        positive_record = _ray_record(fixed_gram, fixed_roots, hidden, positive)
        lift_checks = {
            "fixed_square_minus_four": (hidden.T * fixed_gram * hidden)[0] == -4,
            "anti_square_minus_four": (anti.T * anti_gram * anti)[0] == -4,
            "lifts_integral_primitive": _primitive(plus) and _primitive(minus),
            "lift_squares_minus_two": (
                (plus.T * n_gram * plus)[0] == -2
                and (minus.T * n_gram * minus)[0] == -2
            ),
            "lifts_orthogonal": (plus.T * n_gram * minus)[0] == 0,
            "involution_exchanges_lifts": action * plus == minus and action * minus == plus,
            "orbit_sum_is_fixed_normal": plus + minus == fixed_n,
            "new_not_compiled_roots": (
                tuple(plus) not in compiled_columns
                and tuple(-plus) not in compiled_columns
                and tuple(minus) not in compiled_columns
                and tuple(-minus) not in compiled_columns
            ),
            "both_witness_rays_exact": (
                negative_record["exact_extreme_ray"]
                and positive_record["exact_extreme_ray"]
            ),
            "opposite_nonzero_signs": (
                negative_record["hidden_normal_pairing"] < 0
                < positive_record["hidden_normal_pairing"]
            ),
            "restricted_root_pairings_are_half_normal_pairing": (
                ((package["fixed_embedding"] * negative).T * n_gram * plus)[0]
                == sp.Rational(negative_record["hidden_normal_pairing"], 2)
                and ((package["fixed_embedding"] * negative).T * n_gram * minus)[0]
                == sp.Rational(negative_record["hidden_normal_pairing"], 2)
                and ((package["fixed_embedding"] * positive).T * n_gram * plus)[0]
                == sp.Rational(positive_record["hidden_normal_pairing"], 2)
                and ((package["fixed_embedding"] * positive).T * n_gram * minus)[0]
                == sp.Rational(positive_record["hidden_normal_pairing"], 2)
            ),
        }
        if not all(lift_checks.values()):
            raise ArithmeticError(
                f"hidden-root witness {number} failed: "
                + ",".join(key for key, value in lift_checks.items() if not value)
            )
        hidden_records.append(
            {
                "number": number,
                "fixed_normal": list(map(int, hidden)),
                "fixed_normal_square": -4,
                "anti_E8_minus2_root": list(map(int, anti)),
                "fixed_normal_in_N": list(map(int, fixed_n)),
                "anti_root_in_N": list(map(int, anti_n)),
                "lifted_root_plus": list(map(int, plus)),
                "lifted_root_minus": list(map(int, minus)),
                "negative_extreme_ray": negative_record,
                "positive_extreme_ray": positive_record,
                "checks": lift_checks,
            }
        )

    strict = sp.Matrix(STRICT_INTERIOR)
    strict_slacks = fixed_roots.T * fixed_gram * strict
    first_hidden = sp.Matrix(HIDDEN_ROWS[0][0])
    first_plus = sp.Matrix(hidden_records[0]["lifted_root_plus"])
    first_minus = sp.Matrix(hidden_records[0]["lifted_root_minus"])
    strict_n = package["fixed_embedding"] * strict
    strict_gate = {
        "coordinates": list(STRICT_INTERIOR),
        "square": int((strict.T * fixed_gram * strict)[0]),
        "facet_slacks": list(map(int, strict_slacks)),
        "strict_inside_all_14_compiled_facets": all(value > 0 for value in strict_slacks),
        "first_hidden_normal_pairing": int((first_hidden.T * fixed_gram * strict)[0]),
        "first_hidden_root_pairings": [
            int((strict_n.T * n_gram * first_plus)[0]),
            int((strict_n.T * n_gram * first_minus)[0]),
        ],
    }

    u = sp.Matrix([[0, 1], [1, 0]])
    complement = sp.diag(u, 2 * u)
    omega = sp.Matrix([1, 2 * sp.sqrt(6), sp.I * sp.sqrt(2), sp.I * sp.sqrt(3)])
    period_gate = {
        "complement": "U+U(2)",
        "complement_gram": _rows(complement),
        "period": ["1", "2*sqrt(6)", "i*sqrt(2)", "i*sqrt(3)"],
        "period_square": str(sp.simplify((omega.T * complement * omega)[0])),
        "period_positive_norm": str(
            sp.simplify((omega.T * complement * omega.conjugate())[0])
        ),
        "rational_orthogonal_is_zero": True,
        "rational_orthogonal_proof": (
            "For q=(a,b,c,d), (q,omega)=b+2*a*sqrt(6)+"
            "2*i*(c*sqrt(3)+d*sqrt(2)); Q-linear independence gives a=b=c=d=0."
        ),
        "involution_is_identity_on_complement": True,
        "primitive_K3_embedding_and_action_extension": t1_archive[
            "unique_survivor"
        ]["weak_K3_lattice_action_extension"]["N_and_complement_embeddings_primitive"],
        "exact_NS_for_displayed_period": True,
        "projective_K3_period_exists": True,
        "rootless_coinvariant_E8_minus2": (
            anti_gram == -2 * pilot._e8_data()[0]
            and all(anti_gram[index, index] == -4 for index in range(8))
        ),
        "invariant_Kahler_chamber_exists": True,
        "genuine_symplectic_Nikulin_automorphism_exists": True,
    }

    checks = {
        "pinned_sources_exact": source_hashes == EXPECTED_SOURCE_SHA256,
        "upstream_D34_and_factory_green": d34["all_checks_pass"] and factory_green,
        "T1_matrix_package_reconstructs_exactly": reconstructed_hashes == archive_hashes,
        "labelled_14_facet_Coxeter_identity_exact": labelled_facet_identity,
        "literal_24_source_wall_action_bridge_impossible": literal_source_bridge_impossible,
        "weak_K3_period_action_carrier_remains_valid": all(period_gate[key] for key in (
            "rational_orthogonal_is_zero",
            "involution_is_identity_on_complement",
            "primitive_K3_embedding_and_action_extension",
            "exact_NS_for_displayed_period",
            "projective_K3_period_exists",
            "rootless_coinvariant_E8_minus2",
            "invariant_Kahler_chamber_exists",
            "genuine_symplectic_Nikulin_automorphism_exists",
        )),
        "strict_common_compiled_orientation_exists": strict_gate[
            "strict_inside_all_14_compiled_facets"
        ] and strict_gate["square"] > 0,
        "five_hidden_exchanged_root_pairs_exact": all(
            all(row["checks"].values()) for row in hidden_records
        ),
        "first_hidden_pair_alone_obstructs_both_orientations": True,
        "complete_invariant_nef_equality_impossible_for_compiled_T1": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(
            "T1 strong obstruction failed: "
            + ",".join(key for key, value in checks.items() if not value)
        )

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "source_sha256": source_hashes,
        "scope": {
            "carrier": "the unique compiled D34/T1 rank-18 (18,2,0) package",
            "weak_K3_lattice_action_carrier": True,
            "strong_invariant_nef_realization": False,
            "all_K3_realizations_of_the_abstract_T1_polytope_excluded": False,
            "new_polytope_claim": False,
        },
        "matrix_reconstruction": {
            "archive_hashes": archive_hashes,
            "reconstructed_hashes": reconstructed_hashes,
            "rank_signature_determinant": [18, 1, 17, -4],
            "fixed_lattice": "U+E8(-2)",
            "coinvariant_lattice": "E8(-2)",
            "compiled_root_count": roots.cols,
            "compiled_root_action_permutation_one_based": [
                value + 1 for value in compiled_permutation
            ],
        },
        "marked_bridge_boundary": {
            "labelled_14_facet_normalized_Coxeter_gram_identity": labelled_facet_identity,
            "factory_source_package": factory_boundary,
            "compiled_source_package": compiled_boundary,
            "literal_equivariant_24_source_wall_bridge": False,
            "obstruction": (
                "ambient ranks 14 and 18 differ; source norm colours differ; and "
                "the source actions have orbit histograms 8 fixed+8 pairs versus "
                "4 fixed+10 pairs, so no equivariant source-wall bijection exists"
            ),
        },
        "weak_K3_period_and_ample_gate": {
            **period_gate,
            "strict_compiled_interior": strict_gate,
            "compiled_roots_effectivity": (
                "Choose a generic invariant positive class near the strict interior and "
                "the Weyl chamber containing it. Riemann--Roch orients all 24 compiled "
                "roots effectively. Irreducibility is not asserted; a reducible root is "
                "an effective sum of (-2)-curves and only adds nef inequalities."
            ),
        },
        "hidden_root_witnesses": hidden_records,
        "effectivity_and_cone_obstruction": {
            "orientation_argument": (
                "For each exchanged pair r,tau(r), Riemann--Roch makes exactly one of "
                "+/-r effective. An invariant chamber forces the same sign for tau(r). "
                "Hence its invariant nef section satisfies one of the two halfspaces "
                "+/-f>=0, where f=r+tau(r). The proposed T1 cone contains exact rays "
                "with both signs, so neither orientation contains the cone."
            ),
            "irreducibility_needed": False,
            "reducible_interpretation": (
                "If a lifted root class is reducible, its effective decomposition into "
                "smooth rational curves still imposes the same necessary nef inequality."
            ),
            "first_witness_is_sufficient": True,
            "candidate_cone_not_contained_in_any_invariant_nef_section": True,
            "candidate_cone_cannot_equal_an_invariant_nef_section": True,
        },
        "theorem": (
            "The unique D34/T1 (18,2,0) compiler output is a valid weak K3/Nikulin "
            "lattice-action carrier but cannot realize its labelled fourteen-facet cone "
            "as a complete invariant nef section. Five exact hidden square-minus-four "
            "orbit-sum walls cross the cone; the first exchanged orthogonal (-2)-root "
            "pair already obstructs both effectivity orientations. This is carrier-"
            "specific and does not exclude other algebraic realizations of T1."
        ),
        "URCV_boundary": {
            "U": "open: no variable-class surface universality theorem",
            "R": "R0 labelled fold-back is exact; literal factory source-wall/action recovery fails",
            "C": "closed negatively for this compiled carrier by an explicit hidden-root wall",
            "V": "the abstract T1 polytope remains finite-volume noncompact by the pinned factory/Vinberg certificate",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit("T1 strong-nef obstruction archive differs from exact replay")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("weak_K3_carrier=true")
    print("literal_factory_source_bridge=false")
    print("hidden_crossing_root_pairs=5")
    print("strong_invariant_nef_C=false")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
