from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path

import pytest

from scripts.certify_h9_14_record30_a2_schedule_no_go import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_EXCEPTIONAL_SUMMARY_SHA256,
    EXPECTED_INDEPENDENT_INCIDENCE_SHA256,
    EXPECTED_LABEL_SUMMARY_SHA256,
    EXPECTED_PRODUCTION_SHA256,
    build,
    verify,
)
from scripts.independent_h9_14_record30_a2_exceptional_whole_nerves import (
    _exact_sign,
    _field_norm,
    _residue,
)


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ARCHIVE_FILE_SHA256 = (
    "06759a184708d9c2b8a4ab15fc772fc2a4b7e66ac275d2ff537bebf92cb9ebe6"
)
EXPECTED_TARGET_COUNTS = (
    1,
    1,
    3,
    0,
    0,
    1,
    11,
    2,
    0,
    0,
    1,
    4,
    1,
    0,
    37,
    3,
    3,
    0,
    1,
    5,
    3,
)
EXPECTED_EXCEPTIONAL_TRANSCRIPTS = {
    "orbit7_target1": (
        "fc3207d7b856db2cbb2828e6827eb57e558136cd0bc268796cff6fb69c4c1825"
    ),
    "orbit10_target0": (
        "5f4bd85088843d7cfb97432046a26268bfaa9296165f2ac202ba31a668afabd5"
    ),
}


@pytest.fixture(scope="session")
def certificate() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_frozen_schedule_certificate_and_scope(certificate) -> None:
    assert sha256(ARCHIVE.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_FILE_SHA256
    assert verify(certificate) == {
        "record_index": 30,
        "excluded_compact_schedule": 2,
        "positive_A4_orbit_count": 21,
        "pure_target_orbit_count": 77,
        "sparse_label_no_go_target_count": 74,
        "exceptional_target_count": 3,
        "a2_compact_schedule_excluded": True,
        "whole_record30_shadow_excluded": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    scope = certificate["scope"]
    assert scope["compact_only"] is True
    assert scope["whole_record30_shadow_excluded"] is False
    assert scope["finite_volume_noncompact_chamber_excluded"] is False
    assert scope["new_hyperbolic_coxeter_polytope"] is False


def test_production_and_explicit_incidence_replays_cover_same_77_orbits(
    certificate,
) -> None:
    assert certificate["production_nauty_transcript_sha256"] == (
        EXPECTED_PRODUCTION_SHA256
    )
    assert certificate["independent_explicit_incidence_transcript_sha256"] == (
        EXPECTED_INDEPENDENT_INCIDENCE_SHA256
    )
    comparison = certificate["incidence_cross_audit"]
    assert comparison["production_and_independent_common_counters_identical"]
    assert comparison["production_and_independent_target_orbit_sets_identical"]
    assert comparison["positive_A4_orbit_count"] == 21
    assert comparison["target_orbit_count"] == 77
    rows = certificate["production_nauty_incidence_replay"]["orbit_audits"]
    assert tuple(row["pure_128_facet_target_orbit_count"] for row in rows) == (
        EXPECTED_TARGET_COUNTS
    )
    independent = certificate["independent_explicit_incidence_replay"]
    assert independent["graph_automorphism_count"] == 1_935_360
    assert independent["positive_A4_orbit_count"] == 21
    assert independent["pure_128_facet_target_orbit_count"] == 77
    orbit14 = independent["orbit_audits"][14]
    assert orbit14["stabilizer_order"] == 2_880
    assert orbit14["search_node_count"] == 369_607
    assert orbit14["literal_target_count_with_fixed_A4_pair"] == 4_290
    assert orbit14["target_orbit_size_histogram"] == {
        "60": 17,
        "90": 1,
        "120": 14,
        "180": 1,
        "240": 1,
        "360": 3,
    }


def test_complete_sparse_label_replay_has_exactly_three_survivors(
    certificate,
) -> None:
    replay = certificate["independent_all_target_sparse_label_replay"]
    assert replay["summary_sha256"] == EXPECTED_LABEL_SUMMARY_SHA256
    assert replay["imports_coxeter_polytopes"] is False
    assert replay["complete_order_four_lanner_placement_count"] == 108
    assert replay["complete_order_five_lanner_placement_count"] == 420
    assert replay["target_orbit_count"] == 77
    assert replay["sparse_label_no_go_target_count"] == 74
    assert replay["terminal_sparse_state_count_histogram"] == {
        "0": 74,
        "648": 2,
        "11664": 1,
    }
    assert replay["sparse_label_survivors"] == [
        {
            "positive_A4_orbit_index": 1,
            "orbit_target_index": 0,
            "terminal_sparse_state_count": 11_664,
        },
        {
            "positive_A4_orbit_index": 7,
            "orbit_target_index": 1,
            "terminal_sparse_state_count": 648,
        },
        {
            "positive_A4_orbit_index": 10,
            "orbit_target_index": 0,
            "terminal_sparse_state_count": 648,
        },
    ]


def test_exceptional_whole_nerve_blocks_have_exact_nonzero_determinants(
    certificate,
) -> None:
    replay = certificate["independent_exceptional_exact_whole_nerve_replay"]
    assert replay["summary_sha256"] == EXPECTED_EXCEPTIONAL_SUMMARY_SHA256
    assert replay["imports_coxeter_polytopes"] is False
    assert replay["second_modulus_check"]["proof_role"] == (
        "redundant replay check only"
    )
    assert replay["exact_coefficient_basis"] == [
        "1",
        "sqrt(2)",
        "sqrt(5)",
        "sqrt(10)",
    ]
    for case in replay["cases"]:
        assert case["actual_rank_nine_vertex_count"] == 128
        assert case["terminal_sparse_state_count"] == 648
        assert case["nerve_automorphism_order"] == 1_152
        assert case["sparse_state_orbit_count"] == 6
        assert case["sparse_state_orbit_sizes"] == [72, 144, 144, 72, 144, 72]
        assert case["connected_shell_counts_per_sparse_state"] == {
            "2": 56,
            "3": 128,
            "4": 38,
        }
        assert case["first_empty_connected_shell_size"] == 5
        assert case["exact_quotient_principal_block_count"] == 1_332
        assert case["determinant_sign_histogram"] == {"1": 1_332}
        assert case["zero_exact_field_norm_count"] == 0
        assert case["zero_second_prime_residue_count"] == 0
        assert case["exact_determinant_transcript_sha256"] == (
            EXPECTED_EXCEPTIONAL_TRANSCRIPTS[case["case"]]
        )
        rows = case["exact_determinant_rows"]
        assert len(rows) == 1_332
        for row in rows:
            coefficients = tuple(
                row[
                    "determinant_coefficients_1_sqrt2_sqrt5_sqrt10"
                ]
            )
            assert _field_norm(coefficients) == row["determinant_field_norm"]
            assert row["determinant_field_norm"] != 0
            assert _exact_sign(coefficients) == row["determinant_sign"] == 1
            assert _residue(coefficients) == row["second_prime_residue"] != 0
        assert case[
            "all_literal_assignments_have_rank_at_least_eleven_by_symmetry_transport"
        ]
        assert case["whole_nerve_compact_coxeter_impossible"]


def test_independent_sources_do_not_import_project_label_code() -> None:
    paths = (
        ROOT / "scripts/independent_h9_14_record30_a2_orbit14_labels.py",
        ROOT / "scripts/independent_h9_14_record30_a2_all_target_labels.py",
        ROOT / "scripts/independent_h9_14_record30_a2_exceptional_whole_nerves.py",
    )
    for path in paths:
        source = path.read_text(encoding="utf-8")
        assert "from coxeter_polytopes" not in source
        assert "import coxeter_polytopes" not in source


@pytest.mark.skipif(
    os.environ.get("COXETER_RUN_SLOW_EXACT_REPLAYS") != "1",
    reason="set COXETER_RUN_SLOW_EXACT_REPLAYS=1 for the approximately 8-minute replay",
)
def test_fresh_full_schedule_replay(certificate) -> None:
    assert json.loads(json.dumps(build(), sort_keys=True)) == json.loads(
        json.dumps(certificate, sort_keys=True)
    )
