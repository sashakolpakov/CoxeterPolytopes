import json
from collections import Counter
from pathlib import Path

from scripts.certify_borcherds_char2_all_cyclic_elements import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "borcherds_char2_all_cyclic_element_completion.json"
)


def _payload() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_full_element_class_denominator() -> None:
    payload = _payload()
    reconstruction = payload["reconstruction"]
    assert reconstruction["full_order"] == 241920
    assert reconstruction["geometric_subgroup_order"] == 120960
    assert reconstruction["class_count"] == 28
    classes = payload["classes"]
    assert len(classes) == 28
    assert sum(record["class_size"] for record in classes) == 241920
    assert [record["stable_index"] for record in classes] == list(range(1, 29))
    assert len(
        {record["invariant_class_fingerprint_sha256"] for record in classes}
    ) == 28
    assert Counter(record["order"] for record in classes) == {
        1: 1,
        2: 4,
        3: 3,
        4: 3,
        5: 1,
        6: 6,
        7: 1,
        8: 4,
        10: 1,
        12: 1,
        14: 1,
        15: 1,
        21: 1,
    }
    assert Counter(record["geometric"] for record in classes) == {
        True: 16,
        False: 12,
    }


def test_rank_gate_leaves_exactly_four_nonidentity_casts() -> None:
    payload = _payload()
    gate = payload["rank_gate"]
    assert gate["nonidentity_high_class_count"] == 4
    assert gate["nonidentity_high_class_names"] == ["2a", "2c", "2d", "3b"]
    classes = {record["class_name"]: record for record in payload["classes"]}
    assert {
        name: (classes[name]["fixed_rank"], classes[name]["hyperbolic_dimension"])
        for name in ("2a", "2c", "2d", "3b")
    } == {
        "2a": (14, 13),
        "2c": (15, 14),
        "2d": (11, 10),
        "3b": (12, 11),
    }
    assert max(
        record["fixed_rank"]
        for record in payload["classes"]
        if record["order"] > 1
        and record["class_name"] not in {"2a", "2c", "2d", "3b"}
    ) <= 10


def test_all_four_high_classes_replay_current_casts_exactly() -> None:
    payload = _payload()
    classes = {record["class_name"]: record for record in payload["classes"]}
    expected = {
        "2a": (50, "SS2-2A-H13-50"),
        "2c": (98, "SS2-2C-H14-98"),
        "2d": (47, "SS2-2D-H10-47"),
        "3b": (22, "SS2-3B-H11-22"),
    }
    for name, (facets, identifier) in expected.items():
        screen = classes[name]["high_rank_screen"]
        assert screen["negative_facet_count"] == facets
        assert screen["all_negative_restrictions_are_genuine_facets"] is True
        assert screen["coxeter"]["all_pairs_are_exact_coxeter_relations"] is True
        assert screen["all_facet_reflections_integral"] is True
        assert screen["coxeter_reflective_candidate"] is True
        assert screen["deduplication"] == {
            "current_cast_id": identifier,
            "current_cast_facet_gram_sha256": screen["facet_gram_sha256"],
            "exact_facet_gram_match": True,
            "conclusion": "same pinned chamber, not an additional cyclic cast",
        }
    assert classes["2a"]["geometric"] is True
    assert classes["2c"]["geometric"] is False
    assert classes["2d"]["geometric"] is True
    assert classes["3b"]["geometric"] is True


def test_full_cyclic_certificate_rebuilds_exactly() -> None:
    assert build_certificate() == _payload()
