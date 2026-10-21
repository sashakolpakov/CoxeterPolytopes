import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from scripts.probe_shimada_disc11_hyperideal_escape import (
    FROZEN_DIRECTED_OUTPUT_HASHES,
    FROZEN_TERMINAL_ROOTS_SHA256,
    _extract_checkpoint,
)


ARTIFACT = Path("data/candidates/shimada_disc11_directed_309_plus5.json")
REPLAY_AUDIT = Path("data/candidates/shimada_disc11_directed_329_replay_audit.json")
PAYLOAD_SHA256 = "c527c7eed01f2defbd3f46adbc5c7df8384e1616e9b0bbdb05500f6adec336e4"
TERMINAL_ROOTS_SHA256 = "eb11488e15a822878b2d1f9b41a0039fdc289b6d887d8fb4c829206a3d6b6f16"
TERMINAL_PAIRINGS_SHA256 = "8c031d3c5f1fedd21e50f4fb080a11cbfb1bf8a57a79340ca6646590d0074578"
NEW_OUTPUT_HASHES = (
    "0255b9b41468b3cc8df241b40dc6645e504656e74165e7412508c2c0228dd751",
    "d5ee045260b0034e3370561b3cacd5cd74ec418da12ad0bfa09fb408fbd265f5",
    "464bfb55546dbbc1f357ab1681b75889364aab38ba3bb36dfa31d4c3fdaf3f94",
    "340d0e8e3f99e388592b3c57bf2725e0244fdf83c533a426dccac0a6b6ce1060",
    TERMINAL_ROOTS_SHA256,
)


def _artifact() -> dict:
    return json.loads(ARTIFACT.read_text())


def test_329_checkpoint_extends_the_frozen_309_chain() -> None:
    payload, digest = _extract_checkpoint(_artifact())
    assert digest == PAYLOAD_SHA256
    steps = payload["steps"]
    assert len(steps) == 35
    assert [step["iteration"] for step in steps] == list(range(3, 38))
    assert tuple(step["output"]["roots_sha256"] for step in steps[:30]) == (
        FROZEN_DIRECTED_OUTPUT_HASHES
    )
    assert steps[30]["input_roots_sha256"] == FROZEN_TERMINAL_ROOTS_SHA256
    assert tuple(step["output"]["roots_sha256"] for step in steps[30:]) == (
        NEW_OUTPUT_HASHES
    )
    assert [step["output"]["wall_count"] for step in steps[30:]] == [
        313,
        317,
        321,
        325,
        329,
    ]
    assert [step["shimada_location"]["auxiliary_gallery_length"] for step in steps[30:]] == [
        138,
        105,
        64,
        108,
        55,
    ]
    for step in steps[30:]:
        assert step["output"]["new_distinct_wall_count"] == 4
        assert step["output"]["pairwise_acute"] is True
        assert step["output"]["target_ray_orbit_surviving_after_wall"] == 0
        assert step["new_wall"]["V4_orbit_size"] == 4
        assert step["new_wall"]["physical_NS_square"] == "-2"
        assert step["new_wall"]["field_root_norm"] == {"a": "2", "b": "0"}
        assert step["new_wall"]["strict_timelike_face_witness"][
            "strict_on_all_prior_walls"
        ] is True


def test_329_terminal_outer_cone_still_has_an_exact_hyperideal_ray() -> None:
    payload, _ = _extract_checkpoint(_artifact())
    terminal = payload["terminal"]
    assert terminal["wall_count"] == 329
    assert terminal["roots_sha256"] == TERMINAL_ROOTS_SHA256
    assert terminal["escape_witness_count"] == 1
    witness = terminal["escape_search"]["exact_noncompact_witnesses"][0]
    assert witness["kind"] == "hyperideal"
    assert Golden.from_json(witness["norm"]).sign() > 0
    assert witness["all_wall_pairings_nonpositive"] is True
    assert witness["wall_pairings_sha256"] == TERMINAL_PAIRINGS_SHA256
    assert witness["active_facets"] == [4, 9, 13, 16, 20, 38, 58, 62, 193]


def test_full_external_atlas_replay_covers_all_35_galleries() -> None:
    audit = json.loads(REPLAY_AUDIT.read_text())
    assert audit["verified"] is True
    assert audit["payload_sha256"] == PAYLOAD_SHA256
    assert audit["replayed_orbits"] == 35
    assert audit["replayed_galleries"] == 35
    assert audit["replayed_gallery_steps"] == 2687
    assert audit["frozen_309_prefix_match"] is True
    assert audit["terminal_equals_frozen_309"] is False
    assert audit["terminal_wall_count"] == 329
    assert audit["terminal_roots_sha256"] == TERMINAL_ROOTS_SHA256
    assert audit["terminal_escape_kind"] == "hyperideal"
    relations = audit["terminal_coxeter_relation_counts"]
    assert relations == {
        "m2": 2816,
        "m3": 816,
        "m5": 472,
        "ultraparallel": 49852,
    }
    assert sum(relations.values()) == 329 * 328 // 2
