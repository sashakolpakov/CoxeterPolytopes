import json
from functools import lru_cache
from pathlib import Path

from scripts.certify_raw_adjunction_log_stage0 import (
    ARCHIVE,
    FIXTURES,
    INPUT_SCHEMA,
    build_calibration,
    build_fixtures_from_sources,
    solve_raw_input,
)


@lru_cache(maxsize=1)
def calibration():
    return build_calibration()


def cases_by_id():
    return {case["case_id"]: case for case in calibration()["cases"]}


def test_fixture_projection_is_raw_and_reproducible():
    frozen = json.loads(FIXTURES.read_text(encoding="utf-8"))
    assert frozen == build_fixtures_from_sources()
    for case in frozen["cases"]:
        assert set(case["raw_input"]) == {"schema", "picard_gram", "wall_vectors"}
        assert case["raw_input"]["schema"] == INPUT_SCHEMA
        forbidden = {
            "canonical_class",
            "boundary_indices_zero_based",
            "boundary_weights",
            "surface",
            "characteristic",
            "finite_volume",
            "carrier",
        }
        assert forbidden.isdisjoint(case["raw_input"])


def test_all_six_calibrations_and_mutations_pass():
    certificate = calibration()
    assert certificate["status"] == "PASS"
    assert certificate["all_checks_pass"]
    assert len(certificate["cases"]) == 6
    assert all(case["all_checks_pass"] for case in certificate["cases"])
    assert all(
        all(case["mutation_checks"].values()) for case in certificate["cases"]
    )


def test_toric_h6_recovers_one_weak_branch_and_rejects_false_doubling():
    cert = cases_by_id()["toric_h6_9"]["certificate"]
    assert len(cert["canonical_candidates"]) == 2
    assert cert["accepted_candidate_count"] == 1
    accepted = cert["accepted_candidates"][0]
    assert accepted["branch"] == "weak_minus_one_minus_two"
    assert accepted["canonical_class"] == [-3, 1, 1, 1, 1, 1, 1]
    assert accepted["physical_square_histogram"] == {"-2": 6, "-1": 3}
    rejected = next(candidate for candidate in cert["canonical_candidates"] if not candidate["accepted"])
    assert rejected["physical_square_histogram"] == {"-4": 3, "-2": 6}
    assert not rejected["log_equation"]


def test_four_coble_rows_infer_boundary_and_K_without_markings():
    expected = {
        "coble_type_ii_h10_30": ({"-4": 1, "-2": 12, "-1": 17}, 1),
        "coble_type_i_h11_26": ({"-4": 2, "-2": 10, "-1": 14}, 2),
        "coble_type_vi_char3_h14_60": ({"-4": 5, "-2": 10, "-1": 45}, 5),
        "coble_type_mii_char3_h17_48": ({"-4": 8, "-2": 24, "-1": 16}, 8),
    }
    rows = cases_by_id()
    for case_id, (histogram, boundary_count) in expected.items():
        cert = rows[case_id]["certificate"]
        assert cert["accepted_candidate_count"] == 1
        candidate = cert["accepted_candidates"][0]
        assert candidate["branch"] == "coble_index_two_log_boundary"
        assert candidate["physical_square_histogram"] == histogram
        assert len(candidate["boundary_indices_zero_based_inferred"]) == boundary_count
        assert candidate["ramification_indices_inferred"] == [2] * boundary_count
        boundary_scales = [
            candidate["physical_scales"][index]
            for index in candidate["boundary_indices_zero_based_inferred"]
        ]
        if case_id == "coble_type_ii_h10_30":
            assert boundary_scales == [2]  # C=2r with primitive r^2=-1
            assert candidate["all_reflections_integral_on_NS"]
        else:
            assert boundary_scales == [1] * boundary_count  # primitive B^2=-4
            assert not candidate["all_reflections_integral_on_NS"]
        assert candidate["boundary_pairwise_disjoint"]
        assert candidate["log_residual_K_plus_Delta"] == [0] * cert["validation"]["rank"]
        assert candidate["reflection_module"][
            "all_wall_reflections_preserve_inferred_adapted_module"
        ]


def test_h10_alternative_characteristic_K_is_not_misclassified():
    cert = cases_by_id()["coble_type_ii_h10_30"]["certificate"]
    assert len(cert["canonical_candidates"]) == 2
    rejected = next(candidate for candidate in cert["canonical_candidates"] if not candidate["accepted"])
    assert rejected["physical_square_histogram"] == {"-4": 17, "-2": 12, "-1": 1}
    assert rejected["standard_boundary_shapes"]
    assert not rejected["log_equation"]


def test_maximal_F4_triangle_survives_arithmetic_extension_but_fails_log_equation():
    cert = cases_by_id()["blowup_F4_maximal_triangle_negative_control"]["certificate"]
    assert cert["accepted_candidate_count"] == 0
    assert len(cert["canonical_candidates"]) == 1
    candidate = cert["canonical_candidates"][0]
    assert candidate["canonical_class"] == [-2, -6, 1]
    assert candidate["physical_square_histogram"] == {"-4": 1, "-1": 2}
    assert candidate["standard_boundary_shapes"]
    assert candidate["boundary_pairwise_disjoint"]
    assert candidate["reflection_module"][
        "all_wall_reflections_preserve_inferred_adapted_module"
    ]
    assert candidate["reflection_module"]["ell_scaled_form_is_integral"]
    assert not candidate["all_reflections_integral_on_NS"]
    assert not candidate["log_equation"]
    assert candidate["log_residual_K_plus_Delta"] == ["-3/2", -6, 1]


def test_finite_exhaustion_bound_and_global_characteristic_gate():
    for case in calibration()["cases"]:
        cert = case["certificate"]
        finite = cert["finite_enumeration"]
        assert finite["basis_pairing_assignments_checked"] == finite["enumeration_upper_bound"]
        assert finite["enumeration_upper_bound"] <= 2 ** cert["validation"]["rank"]
        for candidate in cert["canonical_candidates"]:
            gram = cert["input"]["picard_gram"]
            canonical = candidate["canonical_class"]
            covector = [
                sum(gram[i][j] * canonical[j] for j in range(len(canonical)))
                for i in range(len(canonical))
            ]
            assert all((covector[i] - gram[i][i]) % 2 == 0 for i in range(len(canonical)))
            assert candidate["canonical_square"] == 10 - cert["validation"]["rank"]


def test_nonunimodular_and_nonfullspan_inputs_are_rejected():
    nonunimodular = {
        "schema": INPUT_SCHEMA,
        "picard_gram": [[2, 0], [0, -2]],
        "wall_vectors": [[0, 1], [1, 1]],
    }
    assert solve_raw_input(nonunimodular)["status"] == "INVALID_INPUT"

    nonfullspan = {
        "schema": INPUT_SCHEMA,
        "picard_gram": [[1, 0], [0, -1]],
        "wall_vectors": [[0, 1]],
    }
    invalid = solve_raw_input(nonfullspan)
    assert invalid["status"] == "INVALID_INPUT"
    assert not invalid["validation"]["full_span"]


def test_schema_forbids_supplied_canonical_or_boundary_data():
    frozen = json.loads(FIXTURES.read_text(encoding="utf-8"))["cases"][0]["raw_input"]
    with_K = dict(frozen, canonical_class=[0] * len(frozen["picard_gram"]))
    invalid = solve_raw_input(with_K)
    assert invalid["status"] == "INVALID_INPUT"
    assert invalid["validation"]["extra_keys"] == ["canonical_class"]


def test_archive_is_exact_rebuild():
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == calibration()


def test_input_schema_file_has_closed_key_set():
    schema_path = Path("data/fixtures/raw_adjunction_log_stage0_input.schema.json")
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    assert schema["additionalProperties"] is False
    assert set(schema["required"]) == {"schema", "picard_gram", "wall_vectors"}
    assert set(schema["properties"]) == {"schema", "picard_gram", "wall_vectors"}
