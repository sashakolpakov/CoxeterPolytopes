from __future__ import annotations

from collections import Counter
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FAMILY_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_a0_independent_audit.json"
)
DOWNSTREAM_ARCHIVE = (
    ROOT
    / "data/enumerations/h9_14_record18_a0_downstream_independent_audit.json"
)
DOWNSTREAM_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_a0_downstream_audit.py"
)


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def test_exact_literal_and_normalized_family_equality() -> None:
    record = load(FAMILY_ARCHIVE)
    independent = record["independent_survivors"]
    production = record["production_survivors"]
    independent_literal = {
        tuple(row["missing_five_masks"]) for row in independent
    }
    production_literal = {
        tuple(row["missing_five_masks"]) for row in production
    }
    independent_normalized = Counter(
        tuple(row["canonical_family_key"]) for row in independent
    )
    production_normalized = Counter(
        tuple(row["canonical_family_key"]) for row in production
    )
    assert len(independent_literal) == len(production_literal) == 196
    assert independent_literal == production_literal
    assert independent_normalized == production_normalized
    assert sorted(independent_normalized.values()) == [16, 36, 144]
    assert record["literal_family_sets_equal"] is True
    assert record["canonical_normalized_multisets_equal"] is True


def test_frozen_downstream_terminals() -> None:
    record = load(DOWNSTREAM_ARCHIVE)
    assert record["canonical_class_count"] == 3
    assert record["class_multiplicities"] == [16, 36, 144]
    assert record["all_classes_pass_connected_mod2_homology_audit"] is True
    assert record["all_classes_have_unrestricted_compact_coxeter_no_go"] is True
    expected = {
        36: [420, 96, 6, 0],
        144: [420, 96, 6, 0],
        16: [420, 96, 6, 4, 0],
    }
    for audit in record["classes"]:
        topology = audit["topology"]
        labels = audit["coxeter_labels"]
        assert topology["face_vector"] == [
            14,
            86,
            310,
            735,
            1204,
            1372,
            1052,
            495,
            110,
        ]
        assert topology["ridge_degree_histogram"] == {"2": 495}
        assert topology["passes_connected_homology_audit"] is True
        assert labels["terminal_state_count"] == 0
        assert [
            stage["passing_state_count"] for stage in labels["stages"]
        ] == expected[audit["raw_family_multiplicity"]]


def test_no_project_import_downstream_replay() -> None:
    specification = importlib.util.spec_from_file_location(
        "record18_a0_downstream_independent", DOWNSTREAM_SOURCE
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    replay = module.build_record(FAMILY_ARCHIVE)
    frozen = load(DOWNSTREAM_ARCHIVE)
    replay["source"] = frozen["source"]
    assert json.loads(json.dumps(replay, sort_keys=True)) == frozen
