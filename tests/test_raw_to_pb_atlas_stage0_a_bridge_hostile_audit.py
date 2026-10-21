"""Hostile provenance and algebra checks for the raw-to-PB bridge."""

import copy
import json

import pytest
import sympy as sp

import scripts.certify_raw_adjunction_log_stage0 as raw_module
import scripts.certify_raw_to_pb_atlas_stage0_a_bridge as bridge_module
from scripts.certify_raw_adjunction_log_stage0 import solve_raw_input
from scripts.certify_raw_to_pb_atlas_stage0_a_bridge import (
    CASE_TO_ALIAS,
    PB_FIXTURES,
    RAW_FIXTURES,
    candidate_to_pb_input,
)
from scripts.pb_atlas_contraction_core import (
    INPUT_SCHEMA as PB_INPUT_SCHEMA,
    parse_physical_branch,
)


def _raw_cases():
    return {
        case["case_id"]: case
        for case in json.loads(RAW_FIXTURES.read_text(encoding="utf-8"))["cases"]
    }


def _pb_cases():
    return {
        case["calibration_alias"]: case["input"]
        for case in json.loads(PB_FIXTURES.read_text(encoding="utf-8"))["cases"]
    }


def _unique_candidate(raw_input):
    certificate = solve_raw_input(raw_input)
    assert certificate["accepted_candidate_count"] == 1
    return certificate["accepted_candidates"][0]


def test_boundary_weight_is_exactly_the_physical_scale_and_gives_minus_two_K():
    for case_id in CASE_TO_ALIAS:
        raw_input = _raw_cases()[case_id]["raw_input"]
        candidate = _unique_candidate(raw_input)
        gram = sp.Matrix(raw_input["picard_gram"])
        walls = [sp.Matrix(row) for row in raw_input["wall_vectors"]]
        canonical = sp.Matrix(candidate["canonical_class"])
        boundary = candidate["boundary_indices_zero_based_inferred"]
        scales = candidate["physical_scales"]
        weights = [scales[index] for index in boundary]

        physical_boundary = [scales[index] * walls[index] for index in boundary]
        assert sum(physical_boundary, sp.zeros(gram.rows, 1)) == -2 * canonical
        assert all(int((curve.T * gram * curve)[0]) == -4 for curve in physical_boundary)

        payload = candidate_to_pb_input(raw_input, candidate)
        assert payload["boundary_indices_zero_based"] == boundary
        assert payload["boundary_weights"] == weights
        assert sum(
            (
                weight * walls[index]
                for index, weight in zip(
                    payload["boundary_indices_zero_based"], payload["boundary_weights"]
                )
            ),
            sp.zeros(gram.rows, 1),
        ) == -2 * canonical


def test_four_payloads_are_reassembled_without_the_bridge_converter():
    """Rebuild each payload fieldwise, independently of candidate_to_pb_input."""

    raw_cases = _raw_cases()
    native = _pb_cases()
    for case_id, alias in CASE_TO_ALIAS.items():
        raw_input = raw_cases[case_id]["raw_input"]
        candidate = _unique_candidate(raw_input)
        boundary = candidate["boundary_indices_zero_based_inferred"]
        independently_assembled = {
            "schema": PB_INPUT_SCHEMA,
            "picard_gram": copy.deepcopy(raw_input["picard_gram"]),
            "wall_vectors": copy.deepcopy(raw_input["wall_vectors"]),
            "boundary_indices_zero_based": list(boundary),
            "boundary_weights": [
                candidate["physical_scales"][index] for index in boundary
            ],
            "canonical_class": list(candidate["canonical_class"]),
        }
        independently_assembled = parse_physical_branch(
            independently_assembled
        ).input_payload
        assert independently_assembled == native[alias]
        assert independently_assembled["canonical_class"] == candidate["canonical_class"]
        assert independently_assembled["boundary_indices_zero_based"] == boundary
        assert independently_assembled["boundary_weights"] == [
            candidate["physical_scales"][index] for index in boundary
        ]


def test_converter_has_no_runtime_dataflow_from_native_or_expected_fixtures(monkeypatch):
    class Forbidden:
        def __fspath__(self):  # pragma: no cover - only failure path
            raise AssertionError("forbidden fixture path access")

    raw_cases = _raw_cases()
    precomputed = []
    for case_id in CASE_TO_ALIAS:
        raw_input = raw_cases[case_id]["raw_input"]
        precomputed.append((raw_input, _unique_candidate(raw_input)))

    # Neither conversion nor Raw Stage 0 is allowed to consult the frozen PB
    # payload or the source file used to project the calibration fixtures.
    monkeypatch.setattr(bridge_module, "PB_FIXTURES", Forbidden())
    monkeypatch.setattr(bridge_module, "_load", lambda path: (_ for _ in ()).throw(
        AssertionError("converter attempted fixture load")
    ))
    monkeypatch.setattr(raw_module, "PB_FIXTURES", Forbidden())
    for raw_input, candidate in precomputed:
        replayed = solve_raw_input(raw_input)
        assert replayed["accepted_candidates"] == [candidate]
        assert candidate_to_pb_input(raw_input, candidate)["canonical_class"] == candidate[
            "canonical_class"
        ]


def test_poisoning_native_derived_fields_does_not_change_construction_but_breaks_equality():
    raw_cases = _raw_cases()
    native = _pb_cases()
    for case_id, alias in CASE_TO_ALIAS.items():
        raw_input = raw_cases[case_id]["raw_input"]
        derived = candidate_to_pb_input(raw_input, _unique_candidate(raw_input))
        poisoned = copy.deepcopy(native[alias])
        poisoned["canonical_class"][0] += 2
        poisoned["boundary_indices_zero_based"] = list(
            reversed(poisoned["boundary_indices_zero_based"])
        )
        poisoned["boundary_weights"][0] += 1
        assert candidate_to_pb_input(raw_input, _unique_candidate(raw_input)) == derived
        assert poisoned != derived


def test_gram_and_walls_have_shared_fixture_provenance_not_independent_evidence():
    raw_fixtures = json.loads(RAW_FIXTURES.read_text(encoding="utf-8"))
    provenance_key = "data/fixtures/pb_atlas_numerical_stage_a_inputs.json"
    assert provenance_key in raw_fixtures["provenance"]

    native = _pb_cases()
    for case_id, alias in CASE_TO_ALIAS.items():
        raw_input = _raw_cases()[case_id]["raw_input"]
        assert raw_input["picard_gram"] == native[alias]["picard_gram"]
        assert raw_input["wall_vectors"] == native[alias]["wall_vectors"]


def test_each_derived_weight_is_semantically_enforced_by_the_stage_a_equation():
    for case_id in CASE_TO_ALIAS:
        raw_input = _raw_cases()[case_id]["raw_input"]
        payload = candidate_to_pb_input(raw_input, _unique_candidate(raw_input))
        for position in range(len(payload["boundary_weights"])):
            mutated = copy.deepcopy(payload)
            mutated["boundary_weights"][position] += 1
            with pytest.raises(ValueError, match="weighted boundary sum is not -2K"):
                parse_physical_branch(mutated)
