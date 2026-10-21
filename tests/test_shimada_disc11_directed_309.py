import json
from pathlib import Path

import pytest

from coxeter_polytopes.exact.golden import Golden
from scripts.probe_shimada_disc11_hyperideal_escape import (
    CHECKPOINT_SCHEMA_VERSION,
    FROZEN_CHECKPOINT_PAYLOAD_SHA256,
    FROZEN_DIRECTED_OUTPUT_HASHES,
    FROZEN_TERMINAL_ACTIVE_FACETS,
    FROZEN_TERMINAL_NORM,
    FROZEN_TERMINAL_PAIRINGS_SHA256,
    FROZEN_TERMINAL_RAY,
    FROZEN_TERMINAL_ROOTS_SHA256,
    FROZEN_TERMINAL_WALL_COUNT,
    SECOND_EXPANDED_ROOTS_SHA256,
    _canonical_sha,
    _decoded_vector,
    _extract_checkpoint,
)


ARTIFACT = Path("data/candidates/shimada_disc11_directed_309.json")


def _artifact() -> dict:
    return json.loads(ARTIFACT.read_text())


def test_frozen_directed_checkpoint_hash_chain() -> None:
    artifact = _artifact()
    payload, digest = _extract_checkpoint(artifact)
    assert artifact["status"] == (
        "verified_exact_facet_lower_bound_with_surviving_hyperideal_ray"
    )
    relation_counts = artifact["theorem"]["coxeter_relation_counts"]
    assert relation_counts == {
        "m2": 2680,
        "m3": 776,
        "m5": 460,
        "ultraparallel": 43670,
    }
    assert sum(relation_counts.values()) == 309 * 308 // 2
    assert digest == FROZEN_CHECKPOINT_PAYLOAD_SHA256
    assert payload["schema_version"] == CHECKPOINT_SCHEMA_VERSION
    assert payload["base"] == {
        "wall_count": 189,
        "roots_sha256": SECOND_EXPANDED_ROOTS_SHA256,
    }

    steps = payload["steps"]
    assert len(steps) == 30
    assert [step["iteration"] for step in steps] == list(range(3, 33))
    assert tuple(step["output"]["roots_sha256"] for step in steps) == (
        FROZEN_DIRECTED_OUTPUT_HASHES
    )
    previous_hash = SECOND_EXPANDED_ROOTS_SHA256
    for offset, step in enumerate(steps):
        assert step["input_wall_count"] == 189 + 4 * offset
        assert step["input_roots_sha256"] == previous_hash
        assert step["output"]["wall_count"] == 193 + 4 * offset
        assert step["output"]["new_distinct_wall_count"] == 4
        assert step["output"]["pairwise_acute"] is True
        assert step["output"]["target_ray_orbit_surviving_after_wall"] == 0
        assert step["new_wall"]["V4_orbit_size"] == 4
        assert step["new_wall"]["physical_NS_square"] == "-2"
        assert step["new_wall"]["field_root_norm"] == {"a": "2", "b": "0"}
        assert Golden.from_json(
            step["new_wall"]["target_ray_root_pairing"]
        ).sign() > 0
        face = step["new_wall"]["strict_timelike_face_witness"]
        assert face["new_root_pairing_zero"] is True
        assert face["strict_on_all_prior_walls"] is True
        assert Golden.from_json(face["square"]).sign() < 0
        location = step["shimada_location"]
        assert location["kind"] == "root_boundary"
        assert len(location["auxiliary_gallery"]) == location[
            "auxiliary_gallery_length"
        ]
        assert _canonical_sha(location["auxiliary_gallery"]) == location[
            "auxiliary_gallery_sha256"
        ]
        previous_hash = step["output"]["roots_sha256"]


def test_frozen_terminal_ray_is_exactly_hyperideal() -> None:
    payload, _ = _extract_checkpoint(_artifact())
    terminal = payload["terminal"]
    assert terminal["wall_count"] == FROZEN_TERMINAL_WALL_COUNT
    assert terminal["roots_sha256"] == FROZEN_TERMINAL_ROOTS_SHA256
    assert terminal["escape_witness_count"] == 1
    witness = terminal["escape_search"]["exact_noncompact_witnesses"][0]
    assert witness["kind"] == "hyperideal"
    assert _decoded_vector(witness["ray"]) == FROZEN_TERMINAL_RAY
    assert Golden.from_json(witness["norm"]) == FROZEN_TERMINAL_NORM
    assert FROZEN_TERMINAL_NORM.sign() > 0
    assert tuple(witness["active_facets"]) == FROZEN_TERMINAL_ACTIVE_FACETS
    assert witness["all_wall_pairings_nonpositive"] is True
    assert witness["wall_pairings_sha256"] == FROZEN_TERMINAL_PAIRINGS_SHA256


def test_checkpoint_seal_rejects_tampering() -> None:
    artifact = _artifact()
    artifact["checkpoint"]["payload"]["terminal"]["wall_count"] = 310
    with pytest.raises(ValueError, match="payload hash mismatch"):
        _extract_checkpoint(artifact)
