#!/usr/bin/env python3
"""Compose the variable D34 denominator with the strong T2/T1 geometry.

This is a theorem synthesizer, not another search.  It verifies that the
surface-free S16 short4/long10 survivor is exactly the marked T2 package with
the existing strong K3 realization, and that the analogous S17 package is
exactly the carrier killed by the exchanged-root obstruction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import scripts.certify_h9_14_d34_cut_denominator as denominator


ROOT = Path(__file__).resolve().parents[1]
VARIABLE = ROOT / "data/enumerations/h9_14_d34_cut_denominator.json"
T2_PILOT = ROOT / "data/enumerations/h9_14_t2_surface_unfolding_pilot.json"
T2_BRIDGE = ROOT / "data/enumerations/a18_c01_t2_marked_bridge_hostile_audit.json"
PROMOTION = ROOT / "data/enumerations/a18_c01_c07_k3_nef_promotions.json"
T1_OBSTRUCTION = ROOT / "data/enumerations/h9_14_t1_strong_k3_nef_obstruction.json"
OUTPUT = ROOT / "data/enumerations/h9_14_s16_strong_transfer_corollary.json"
SCHEMA = "coxeter.h9_14_s16_strong_transfer_corollary.v1"

SOURCE_HASHES = {
    "variable_denominator": "fa89247c46f7063cba4b26904b1ff9c749647bf3f001354eb107237391cad1f9",
    "T2_surface_unfolding": "e72a7b5690c16c8c405ddf8679cef4bafc18cf93f59f38f5e635acc9131d60c2",
    "T2_marked_bridge": "92396280c6b276979a005f9732408944c97fc71c08998d1f127833c7bacd6595",
    "strong_K3_promotion": "2ab7acc379e48ac55944e92b8d614c0992d7c482af1300c670f2bf8e899ba1b1",
    "T1_exchanged_root_obstruction": "4eb338784ec452cb08fb4a75b6f4be4c43155d1e6fce2a2259bbe39f410037c1",
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _jhash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _load() -> tuple[dict[str, Any], ...]:
    paths = (VARIABLE, T2_PILOT, T2_BRIDGE, PROMOTION, T1_OBSTRUCTION)
    for (name, expected), path in zip(SOURCE_HASHES.items(), paths, strict=True):
        if _sha(path) != expected:
            raise ArithmeticError(f"pinned synthesis input changed: {name}")
    return tuple(json.loads(path.read_text()) for path in paths)


def _support(payload: dict[str, Any], support_id: str) -> dict[str, Any]:
    return next(
        row for row in payload["support_censuses"] if row["support_id"] == support_id
    )["passing_diagram"]


def _generated_gram(record: dict[str, Any]) -> list[list[int]]:
    edges = denominator._zero_based(denominator.SUPPORTS[record["support_id"]])
    mask = int(record["m4_cut_mask_hex"], 16)
    side = sum(
        1 << (vertex - 1)
        for vertex in record["canonical_cut_side_vertices_one_based"]
    )
    gram, _squares = denominator._gram_lane(edges, mask, side, 0)
    return denominator._rows(gram)


def build_certificate() -> dict[str, Any]:
    variable, pilot, bridge, promotion, obstruction = _load()
    if not all(
        payload["all_checks_pass"]
        for payload in (variable, pilot, bridge, promotion, obstruction)
    ):
        raise ArithmeticError("an upstream theorem archive is not green")
    s16 = _support(variable, "S16")
    s17 = _support(variable, "S17")
    s16_lanes = s16["scaled_dual_lanes"]
    s17_lanes = s17["scaled_dual_lanes"]
    s16_package = s16["unique_minimal_involutive_weak_K3_package"]
    s17_package = s17["unique_minimal_involutive_weak_K3_package"]

    generated_t2 = _generated_gram(s16)
    pilot_t2 = [
        [int(value) for value in row]
        for row in pilot["surface_free_logical_input"]["facet_gram"]
    ]
    pilot_survivor = pilot["unique_survivor"]
    exact_package_hashes = {
        "fixed_gram": _jhash(pilot_survivor["fixed_lattice"]["gram"]),
        "rank18_gram": _jhash(pilot_survivor["rank18_unfolding"]["gram"]),
        "rank18_involution": _jhash(
            pilot_survivor["rank18_unfolding"]["involution"]
        ),
        "ordered_24_roots": _jhash(
            pilot_survivor["source_roots_and_fold"]["source_root_coordinates"]
        ),
    }
    variable_package_hashes = {
        "fixed_gram": s16_package["fixed_lattice"]["gram_sha256"],
        "rank18_gram": s16_package["rank18_unfolding"]["gram_sha256"],
        "rank18_involution": s16_package["rank18_unfolding"]["involution_sha256"],
        "ordered_24_roots": s16_package["roots_and_fold"]["roots_sha256"],
    }
    c01 = next(case for case in promotion["cases"] if case["id"] == "A18_4_0_C01")

    obstruction_hashes = obstruction["matrix_reconstruction"]["archive_hashes"]
    variable_t1_hashes = {
        "N_gram": s17_package["rank18_unfolding"]["gram_sha256"],
        "fixed_gram": s17_package["fixed_lattice"]["gram_sha256"],
        "involution": s17_package["rank18_unfolding"]["involution_sha256"],
        "roots": s17_package["roots_and_fold"]["roots_sha256"],
    }

    s16_checks = {
        "variable_S16_has_one_Lorentzian_diagram_orbit": (
            next(
                row for row in variable["support_censuses"] if row["support_id"] == "S16"
            )["rank10_lorentzian_corank4_diagram_orbit_count"] == 1
        ),
        "surface_free_short4_long10_lane_selected": (
            s16_lanes[0]["normal_square_histogram"] == {"2": 4, "4": 10}
            and s16_lanes[0]["stable_even_overlattice_candidate_count"] == 16
            and s16_lanes[0]["weak_K3_survivor_count"] == 1
        ),
        "scaled_dual_lane_not_silently_promoted": (
            s16_lanes[1]["normal_square_histogram"] == {"2": 10, "4": 4}
            and s16_lanes[1]["weak_K3_survivor_count"] == 0
        ),
        "generated_Gram_equals_T2_pilot_entry_for_entry": generated_t2 == pilot_t2,
        "variable_package_equals_T2_pilot_in_exact_coordinates": (
            variable_package_hashes == exact_package_hashes
        ),
        "marked_bridge_is_unimodular_and_exact": (
            bridge["marked_isometry"]["all_three_equalities_exact"]
            and bridge["marked_isometry"]["determinant"] == 1
            and bridge["marked_isometry"]["inverse_integral"]
            and bridge["checks"]["fourteen_folded_sums_and_gram_exact"]
            and bridge["checks"]["full_48_wall_chamber_action_exact"]
        ),
        "intrinsic_strong_K3_moduli_stratum_nonempty": (
            c01["period_and_torelli"]["picard_lattice_is_exactly_N"]
            and c01["period_and_torelli"]["invariant_kahler_chamber_exists"]
            and "there exists a projective complex K3 surface"
            in c01["period_and_torelli"]["global_torelli_conclusion"]
        ),
        "complete_invariant_nef_section_C": (
            "Nef(X)" in c01["nef_and_quotient_cone"]["invariant_nef_section_equality"]
            and c01["nef_and_quotient_cone"]["all_restricted_long_walls_redundant"]
        ),
        "fresh_literal_V": (
            s16["fresh_cone_and_volume_certificate"]["all_14_facets_irredundant"]
            and s16["fresh_cone_and_volume_certificate"]["literal_vinberg"]["decision"]
            == "finite_volume_noncompact"
            and s16["fresh_cone_and_volume_certificate"]["literal_vinberg"]["bad_edge_count"]
            == 0
            and s16["fresh_cone_and_volume_certificate"]["literal_vinberg"]["bad_parabolic_neighborhood_count"]
            == 0
        ),
    }
    s17_checks = {
        "variable_S17_has_one_Lorentzian_diagram_orbit": (
            next(
                row for row in variable["support_censuses"] if row["support_id"] == "S17"
            )["rank10_lorentzian_corank4_diagram_orbit_count"] == 1
        ),
        "surface_free_short4_long10_lane_has_five_candidates_one_weak": (
            s17_lanes[0]["normal_square_histogram"] == {"2": 4, "4": 10}
            and s17_lanes[0]["stable_even_overlattice_candidate_count"] == 5
            and s17_lanes[0]["weak_K3_survivor_count"] == 1
        ),
        "variable_T1_package_is_exact_obstruction_input": (
            variable_t1_hashes == obstruction_hashes
        ),
        "weak_K3_surface_and_genuine_action_exist": (
            obstruction["weak_K3_period_and_ample_gate"]
            ["exact_NS_for_displayed_period"]
            and obstruction["weak_K3_period_and_ample_gate"]
            ["genuine_symplectic_Nikulin_automorphism_exists"]
        ),
        "compiled_strong_moduli_strata_empty": (
            obstruction["effectivity_and_cone_obstruction"]
            ["candidate_cone_cannot_equal_an_invariant_nef_section"]
            and obstruction["checks"]
            ["complete_invariant_nef_equality_impossible_for_compiled_T1"]
        ),
        "fresh_polytope_V_nonempty": (
            s17["fresh_cone_and_volume_certificate"]["literal_vinberg"]["decision"]
            == "finite_volume_noncompact"
            and s17["fresh_cone_and_volume_certificate"]["literal_vinberg"]["bad_edge_count"]
            == 0
            and s17["fresh_cone_and_volume_certificate"]["literal_vinberg"]["bad_parabolic_neighborhood_count"]
            == 0
        ),
    }
    if not all(s16_checks.values()) or not all(s17_checks.values()):
        raise ArithmeticError("the S16 transfer or S17 countercalibration is not exact")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "source_sha256": SOURCE_HASHES,
        "strong_marked_S16_class": {
            "name": "D16_plus",
            "surface_free_definition": {
                "support": "S16 with the fixed fourteen facet labels",
                "relations": "exactly eight m4 support edges, all other support edges m3, all nonedges m2",
                "primitive_norm_square_alphabet": [2, 4],
                "orientation_selector": "exactly ten norm-square-4 facets and four norm-square-2 facets",
                "Gram_rank_inertia": [10, 1, 9, 4],
                "positive_sheet": "the component containing the generated strict cdd interior",
                "does_not_name_a_Gram_surface_or_polytope": True,
            },
            "finite_A_D": {
                "compiler_lift_alphabet": "minimal fixed-root/orthogonal-exchanged-pair lift with coinvariant E8(-2)",
                "diagram_orbits": 1,
                "marked_integral_Gram_lanes_after_orientation_selector": 1,
                "stable_even_fixed_overlattice_packages": 16,
                "weak_K3_packages": 1,
                "strong_nonempty_packages": 1,
            },
            "intrinsic_moduli_stratum": (
                "marked projective K3 surfaces with exact NS=N, the specified "
                "symplectic Nikulin action, the marked effective root orbits, "
                "and equality of the fourteen-wall section with Nef(X) intersect Fix(g); "
                "finite hyperbolic volume is not part of the definition"
            ),
            "checks": s16_checks,
            "obligations": {
                "U": (
                    "closed for D16_plus: the variable census forces the one marked "
                    "T2 Gram and package, and the exact bridge transports every "
                    "realization of that marked Gram, modulo labelled orientation-"
                    "preserving integral isometry, to the nonempty K3 stratum"
                ),
                "R": (
                    "closed: the unimodular bridge transports the exact form, action, "
                    "ordered 24 roots, fourteen orbit sums, and full 48-wall action"
                ),
                "C": "closed by the independent complete invariant-nef-section theorem",
                "V": (
                    "closed afresh by generated strict interior, fourteen irredundancy "
                    "witnesses, 683 elliptic edges, 90 finite and 20 ideal endpoints, "
                    "and zero bad full parabolic neighborhoods"
                ),
            },
            "marked_equivalence": (
                "Pol_fv^mark(D16_plus) is nonempty iff M_D16_plus^alg is nonempty"
            ),
            "marked_equivalence_proved": True,
            "unmarked_discrete_Gram_image_count": 1,
            "unmarked_polytope_type": "the known finite-volume noncompact T2 H9/14 polytope",
        },
        "bare_S16_corollary": {
            "class": "the same S16 datum allowing both scaled-dual norm orientations",
            "statement": (
                "the unmarked polytope existence problem has one congruence type and "
                "is nonempty iff the D16_plus strong K3 stratum is nonempty"
            ),
            "marked_warning": (
                "formal marked U is not asserted for the bare class: the complementary "
                "short10/long4 integral lane has no minimal weak K3 package"
            ),
        },
        "S17_countercalibration": {
            "name": "D17_plus",
            "surface_free_definition": (
                "replace S16 by S17 in D16_plus, retaining the short4/long10 selector"
            ),
            "finite_A_D": {
                "diagram_orbits": 1,
                "stable_even_fixed_overlattice_packages": 5,
                "weak_K3_packages": 1,
                "strong_nonempty_packages": 0,
            },
            "checks": s17_checks,
            "weak_obligations": {
                "U_weak": "closed for the unique compiled package",
                "R0": "exact labelled fourteen-facet fold-back",
                "V": "fresh finite-volume noncompact certificate with 112 finite and 15 ideal endpoints",
            },
            "strong_failure": (
                "four arithmetic failures in D17_plus and the exchanged-root "
                "obstruction on its sole weak survivor make every strong minimal-"
                "E8(-2) stratum empty; the complementary lane has two further "
                "arithmetic failures, although the T1 polytope exists"
            ),
            "minimal_carrier_equivalence_false": True,
            "scope": (
                "this disproves formal U only for the stated minimal K3 carrier alphabet; "
                "other K3 lift alphabets and other algebraic carrier types remain open"
            ),
        },
        "milestone_assessment": {
            "M3b": (
                "closed for the sharply scoped oriented D16_plus class by explicit "
                "separate U/R/C/V proofs"
            ),
            "M3c": (
                "the existence and discrete-Gram-image consequence is closed (one T2 "
                "image), but a full marking/monodromy-orbit classification of the K3 "
                "moduli stratum is not supplied and the broad milestone should remain open"
            ),
            "broad_surface_softening": "open beyond the D16_plus support/orientation denominator",
        },
        "theorem": (
            "The variable surface-free oriented S16 denominator satisfies the full "
            "marked U/R/C/V transfer theorem and has one unmarked Gram image, T2. "
            "The parallel S17 denominator is an exact countercalibration for the "
            "minimal E8(-2) K3 alphabet: its T1 polytope exists and passes V, but all "
            "strong carrier strata are empty."
        ),
        "checks": {
            "all_sources_green": True,
            "S16_full_marked_transfer": all(s16_checks.values()),
            "S17_exact_minimal_carrier_countercalibration": all(s17_checks.values()),
            "scaled_dual_and_broad_scope_warnings_recorded": True,
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("S16 transfer corollary archive differs from exact replay")
        print("archive_check=PASS")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"wrote={args.output}")
    print("S16_marked_URCV=PASS; unmarked_Gram_images=1")
    print("S17_minimal_strong_equivalence=FALSE; polytope=NONEMPTY; strong_moduli=EMPTY")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
