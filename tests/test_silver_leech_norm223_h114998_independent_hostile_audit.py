import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "scripts/independent_hostile_audit_silver_leech_norm223_h114998.py"
)
ARTIFACT = (
    ROOT
    / "data/enumerations/"
    "silver_leech_norm223_h114998_independent_hostile_audit.json"
)


def load_audit():
    specification = importlib.util.spec_from_file_location(
        "silver_leech_norm223_h114998_hostile", SCRIPT
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


AUDIT = load_audit()


def frozen() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_independent_full_binary_replay_and_limiting_locus() -> None:
    payload = AUDIT.encode(AUDIT.build_payload())
    census = payload["independent_binary_replay"]
    assert census["raw_type_count_in_complete_prefix"] == 2875
    assert census["short_combination_exclusion_count_in_complete_prefix"] == 20
    assert census["surviving_type_count"] == 2855
    assert census["limiting_surviving_type_count"] == 216
    assert census["independently_tightened_component_count"] == 2
    assert [row["triple"] for row in census["maximizers"]] == [
        [115270, 57364, 81316],
        [115266, 57366, 81316],
    ]
    assert [row["shifted_difference_square"] for row in census["maximizers"]] == [
        2,
        0,
    ]
    cross = payload["production_cross_check"]
    assert cross["sha256"] == (
        "6250394c4c24b129bb3fff6839ba5a020708985b5c96756c903f71a158285c4f"
    )
    assert cross["consulted_only_after_independent_replay"] is True
    assert cross["verdict_and_key_exact_fields_agree"] is True
    assert payload == frozen()


def test_explicit_root_and_complete_chronological_acceptance() -> None:
    payload = frozen()
    root = payload["limiting_root"]
    assert root["d_e_gram"] == [[4, -1], [-1, 6]]
    assert root["u_v_gram"] == [115266, 57366, 81316]
    assert root["direction_pairing"] == {"a": "679", "b": "479"}
    assert root["root_square"] == {"a": "2", "b": "0"}
    assert root["primitive"] is True
    assert root["reflection_integral"] is True

    chronology = payload["earlier_chronology"]
    assert chronology["component_count"] == 34
    assert chronology["binary_type_count"] == 2639
    assert chronology["all_strictly_negative"] is True
    assert chronology["tight_control_ordinal"] == 565
    assert chronology["tight_threshold_gap"] == {"a": "63/4", "b": "45/4"}
    assert chronology["tight_rejection_discriminant"] == {"a": "57", "b": "41"}

    endpoint = payload["endpoint_chronology"]
    assert endpoint["height_zero_empty"] is True
    assert endpoint["same_height_q1_empty"] is True
    assert endpoint["same_height_q2_complete_trace_lemma"][
        "positive_distinct_pairing_impossible"
    ] is True
    assert endpoint["same_height_acute"] is True
    assert payload["verdict"]["chronologically_vinberg_accepted"] is True
    assert payload["verdict"]["new_compact_example"] is False
