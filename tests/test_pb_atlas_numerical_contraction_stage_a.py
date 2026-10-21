from __future__ import annotations

import copy
import json
from pathlib import Path

import jsonschema
import pytest

from scripts.certify_pb_atlas_numerical_stage_a import build_certificate
from scripts.freeze_pb_atlas_stage_a_fixtures import build_fixtures
from scripts.pb_atlas_contraction_core import (
    _permutation_lift,
    enumerate_frame_orbits,
    parse_physical_branch,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_a_inputs.json"
ARCHIVE = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
INPUT_SCHEMA = ROOT / "schemas/pb-atlas-physical-branch.schema.json"
OUTPUT_SCHEMA = ROOT / "schemas/pb-atlas-numerical-stage-a.schema.json"


def _fixtures() -> dict:
    return json.loads(FIXTURES.read_text(encoding="utf-8"))


def _archive() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_frozen_fixture_extractor_is_exact() -> None:
    assert build_fixtures() == _fixtures()


def test_every_arrays_only_input_validates_against_schema() -> None:
    schema = json.loads(INPUT_SCHEMA.read_text(encoding="utf-8"))
    for case in _fixtures()["cases"]:
        jsonschema.validate(case["input"], schema)
        branch = parse_physical_branch(case["input"])
        assert branch.wall_sublattice_index == 1
        assert branch.wall_count == len(case["input"]["wall_vectors"])


def test_parser_rejects_labels_but_accepts_nonsaturated_full_span() -> None:
    payload = copy.deepcopy(_fixtures()["cases"][0]["input"])
    payload["known_surface_name"] = "forbidden routing datum"
    with pytest.raises(ValueError, match="accepts arrays only"):
        parse_physical_branch(payload)

    payload = copy.deepcopy(_fixtures()["cases"][0]["input"])
    # Full-rank, primitive, acute, index-two subset found without using an
    # endpoint witness.  The marked half-boundary remains at position zero.
    subset = (0, 4, 5, 6, 10, 18, 19, 22, 26, 27, 28)
    payload["wall_vectors"] = [payload["wall_vectors"][index] for index in subset]
    payload["boundary_indices_zero_based"] = [0]
    branch = parse_physical_branch(payload)
    assert branch.wall_sublattice_index == 2


def test_frozen_archive_and_output_packets_validate() -> None:
    archive = _archive()
    schema = json.loads(OUTPUT_SCHEMA.read_text(encoding="utf-8"))
    assert archive["all_checks_pass"] is True
    assert archive["status"] == "PASS_PB_ATLAS_NUMERICAL_STAGE_A"
    assert len(archive["cases"]) == 4
    for case in archive["cases"]:
        assert case["all_checks_pass"] is True
        assert all(case["expected_checks"].values())
        assert all(case["mutation_checks"].values())
        certificate = case["certificate"]
        jsonschema.validate(certificate, schema)
        for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits"):
            for orbit in certificate[family]:
                program = orbit["program"]
                assert orbit["orbit_size"] == len(orbit["member_masks"])
                assert program["stage_B_payload_complete"] is True
                assert len(program["all_wall_coordinates_in_endpoint_blowup_basis"]) == certificate[
                    "input"
                ]["wall_count"]
                assert set(program["contracted_wall_indices_zero_based"]) == set(
                    program["sequential_contraction_order_zero_based"]
                )


def test_exact_calibration_counts_and_new_high_rank_endpoints() -> None:
    cases = _archive()["cases"]
    assert cases[0]["certificate"]["P2_labelled_count"] == 573
    assert cases[1]["certificate"]["P2_labelled_count"] == 52
    assert cases[2]["certificate"]["hirzebruch_labelled_histogram"] == {
        "F_0": 105,
        "F_1": 310,
        "F_2": 120,
        "F_3": 20,
    }
    assert cases[2]["certificate"]["P2_labelled_count"] == 5
    assert cases[3]["certificate"]["hirzebruch_labelled_histogram"] == {"F_0": 1}
    assert cases[3]["certificate"]["P2_labelled_count"] == 0
    assert _archive()["aggregate_terminal_histogram"] == {
        "F_0": 1029,
        "F_1": 2586,
        "F_2": 564,
        "F_3": 44,
        "F_4": 12,
        "P2": 630,
    }
    hostile = _archive()["nonsaturated_symmetry_hostile_fixture"]
    assert hostile["all_checks_pass"] is True
    assert all(hostile["checks"].values())
    assert hostile["integrally_liftable_group_order"] < hostile[
        "raw_marked_Gram_group_order"
    ]


def test_h17_orbit_quotient_replays_every_depth_without_surface_labels() -> None:
    fixture = _fixtures()["cases"][3]
    frozen = _archive()["cases"][3]["certificate"]
    branch = parse_physical_branch(fixture["input"])
    automorphisms = frozen["symmetry_quotient"][
        "automorphism_permutations_zero_based"
    ]
    enumeration = enumerate_frame_orbits(branch, automorphisms=automorphisms)
    assert [row["labelled_frame_count"] for row in enumeration["level_census"]] == fixture[
        "expected"
    ]["labelled_frame_counts_by_depth"]
    assert [row["orbit_count"] for row in enumeration["level_census"]] == [
        1,
        1,
        3,
        7,
        25,
        62,
        144,
        240,
        329,
        317,
        243,
        124,
        52,
        12,
        4,
        1,
        1,
        0,
    ]


def test_nontrivial_wall_gram_symmetries_have_integral_picard_lifts() -> None:
    fixtures = _fixtures()
    archive = _archive()
    for fixture, frozen in zip(fixtures["cases"], archive["cases"]):
        branch = parse_physical_branch(fixture["input"])
        permutations = frozen["certificate"]["symmetry_quotient"][
            "automorphism_permutations_zero_based"
        ]
        nontrivial = next(
            permutation
            for permutation in permutations
            if permutation != list(range(branch.wall_count))
        )
        lift = _permutation_lift(branch, nontrivial)
        assert lift is not None
        assert abs(int(lift.det())) == 1
        assert lift.T * branch.gram * lift == branch.gram
        assert lift * branch.canonical == branch.canonical


@pytest.mark.slow
def test_full_frozen_certificate_replay() -> None:
    assert build_certificate() == _archive()
