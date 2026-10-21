from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys

import sympy as sp

from scripts.certify_coble_finite_geometry_raw_recognizer import (
    build_manifest_certificate,
    recognize_wall_data,
)
from scripts.generate_coble_finite_geometry_raw_fixtures import build_fixtures


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "data/fixtures/coble_finite_geometry_raw"
MANIFEST = FIXTURES / "manifest.json"
ARCHIVE = ROOT / "data/enumerations/coble_finite_geometry_raw_recognizer.json"
SCRIPT = ROOT / "scripts/certify_coble_finite_geometry_raw_recognizer.py"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _positive_payloads() -> list[dict]:
    return [_load(FIXTURES / f"positive_{index:02d}.json") for index in range(1, 6)]


def test_frozen_fixtures_are_mechanical_normalizations_of_the_five_rows() -> None:
    files, manifest = build_fixtures()
    assert _load(MANIFEST) == manifest
    for name, payload in files.items():
        assert _load(FIXTURES / name) == payload
        assert set(payload) == {
            "schema",
            "lattice_gram",
            "primitive_wall_vectors",
            "orientation",
        }
    assert len([case for case in manifest["cases"] if case["expected_accepted"]]) == 5


def test_all_five_structural_cores_are_accepted_from_raw_arrays() -> None:
    results = [recognize_wall_data(payload) for payload in _positive_payloads()]
    assert all(result["accepted"] for result in results)
    assert {result["finite_geometry"]["core"] for result in results} == {
        "polar_chord",
        "grid_projectivity",
        "mobius_S_3_4_10",
        "desargues_triangle_cut",
        "tetrahedral_4_plus_6",
    }
    assert all(result["independence"]["archive_provenance_used"] is False for result in results)
    assert all(result["independence"]["stored_output_Gram_comparison"] is False for result in results)


def test_nontrivial_unimodular_basis_change_and_wall_permutation_are_invariant() -> None:
    for payload in _positive_payloads():
        baseline = recognize_wall_data(payload)
        gram = sp.Matrix(payload["lattice_gram"])
        rank = gram.rows
        change = sp.eye(rank)
        change[0, 1] = 1
        change[2, 0] = -1
        assert int(change.det()) == 1
        inverse = change.inv()
        transformed = copy.deepcopy(payload)
        new_gram = change.T * gram * change
        transformed["lattice_gram"] = [
            [int(value) for value in new_gram.row(index)] for index in range(rank)
        ]
        transformed["primitive_wall_vectors"] = [
            [int(value) for value in inverse * sp.Matrix(wall)]
            for wall in reversed(payload["primitive_wall_vectors"])
        ]
        replay = recognize_wall_data(transformed)
        assert replay["accepted"]
        assert replay["finite_geometry"]["core"] == baseline["finite_geometry"]["core"]
        assert replay["finite_geometry"] == baseline["finite_geometry"]
        assert replay["input_summary"] == baseline["input_summary"]


def test_tetrahedral_core_reconstructs_all_five_intrinsic_frames() -> None:
    result = recognize_wall_data(_positive_payloads()[4])
    assert result["accepted"]
    design = result["finite_geometry"]
    assert design["core"] == "tetrahedral_4_plus_6"
    assert design["plane_frame_count"] == 5
    assert design["frame_center_frequency_histogram"] == {"0": 35, "2": 10, "5": 6}
    assert design["pairwise_frame_intersection_histogram"] == {"7": 10}
    assert design["root_kind_histogram_per_frame"] == {
        "conic": 4,
        "line": 12,
        "quartic": 4,
    }
    assert design["colored_design_isomorphism_count_per_frame"] == 24
    assert all(
        frame["all_full_wall_coordinates_integral"]
        and frame["exceptional_basis_determinant_absolute"] == 1
        for frame in design["frame_records"]
    )


def test_decoy_archive_type_hash_and_stored_gram_metadata_are_ignored() -> None:
    for payload in _positive_payloads():
        payload = copy.deepcopy(payload)
        baseline = recognize_wall_data(payload)
        payload.update(
            {
                "archive_name": "deliberately-wrong-row.json",
                "known_type_label": "deliberately wrong",
                "source_hash": "not-a-digest",
                "stored_output_gram": [[999]],
            }
        )
        assert recognize_wall_data(payload) == baseline


def test_acute_orientation_mutations_preserve_crude_histograms_but_are_rejected() -> None:
    certificate = build_manifest_certificate(MANIFEST)
    mutations = [
        row for row in certificate["case_results"] if row.get("preserves_crude_signature")
    ]
    assert len(mutations) == 5
    assert all(row["crude_signature"] == row["baseline_crude_signature"] for row in mutations)
    assert all(row["recognition"]["accepted"] is False for row in mutations)
    assert all(row["recognition"]["stage"] == "orientation" for row in mutations)
    scope = certificate["hostile_mutation_scope"]
    assert scope["orientation_preserving_deeper_incidence_mutation_frozen"] is False


def test_schema_primitivity_and_orientation_are_real_gates() -> None:
    nonprimitive = _load(FIXTURES / "negative_nonprimitive.json")
    assert recognize_wall_data(nonprimitive)["stage"] == "primitivity"
    missing_orientation = copy.deepcopy(_positive_payloads()[0])
    del missing_orientation["orientation"]
    assert recognize_wall_data(missing_orientation)["stage"] == "orientation"


def test_genuine_single_input_cli() -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(FIXTURES / "positive_05.json")],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    payload = json.loads(completed.stdout)
    assert payload["accepted"]
    assert payload["finite_geometry"]["core"] == "tetrahedral_4_plus_6"


def test_frozen_certificate_matches_live_manifest_replay() -> None:
    assert build_manifest_certificate(MANIFEST) == _load(ARCHIVE)


def test_five_gate_rating_is_scoped_thin_and_not_global() -> None:
    payload = _load(ARCHIVE)
    assert payload["recognized_core_class"]["name"] == "C_FG"
    gates = payload["five_gate_rating"]
    assert gates["C_independent_P_only_class"]["verdict"] == "PASS_SCOPED"
    assert gates["Car_nonempty_for_every_recognized_core"]["verdict"] == "PASS_SCOPED_ROOT_SEED"
    assert gates["Return_exact_same_P"]["verdict"] == "CONDITIONAL_PASS"
    assert gates["AG_independent_classification"]["verdict"] == "CONDITIONAL_PASS"
    assert gates["Pushback_new_theorem_on_every_member"]["verdict"] == "FAIL_NONTHIN_NOVELTY"
    assert gates["overall"].startswith("L3_PLUS_CLOSED_GENERIC_RECOGNITION_GATE")
    layers = payload["four_layer_theorem"]
    assert set(layers) == {
        "i_recognized_grammar",
        "ii_uniform_explicit_realization",
        "iii_exact_same_P_return",
        "iv_Kondo_pushback",
    }
