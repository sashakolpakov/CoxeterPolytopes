from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.json"
)
HIGH_A = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json"
)
REPORT = ROOT / "results/gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.md"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("partition33_cc_tower", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_archive_covers_every_canonical_black_record() -> None:
    module = load_module()
    data = archived(ARCHIVE)
    expected = [tuple(value) for value in module.canonical_black_records(0)]
    assert len(expected) == len(data["rows"]) == 1292
    assert data["prime_histogram"] == {"101": 1277, "109": 15}
    actual = []
    for row in data["rows"]:
        record = (row["a_lambda"], row["b"], row["c"], row["d"], row["e"])
        actual.append(record)
        assert sum(record) == 29
        assert record[1] <= record[2] <= record[3] <= record[4]
        assert row["base_dimension"] == 0
        assert row["base_quotient_dimension"] == 24
        assert row["successive_endpoint_standard_basis"] == [1]
        assert row["singular_stdout"] == "BASE_DIM\n0\n24\nUNIT"
        assert row["structural_homogenization_degrees"] == [4, 3, 2, 1, 29]
        assert row["uniform_partition_infinity_lemma_applied"] is True
    assert actual == expected


def test_predetermined_tranche_and_high_a_overlap_are_exact() -> None:
    module = load_module()
    data = archived(ARCHIVE)
    high_a = archived(HIGH_A)
    tranche = [
        value
        for value in module.canonical_black_records(0)
        if 17 <= value[0] <= 20
    ]
    assert len(tranche) == data["predetermined_tranche"]["records"] == 41
    assert data["predetermined_tranche"]["passed_before_full_CC_run"] is True
    assert high_a["black_degree_records_eliminated"] == 12
    assert data["previously_eliminated_CC_overlap_from_high_a_theorem"] == 12
    assert data["CC_marked_passports_certified"] - 12 == 1280


def test_fast_representative_towers_recompute() -> None:
    module = load_module()
    for weights in (
        (1, 1, 1, 1, 25),
        (17, 3, 3, 3, 3),
        (25, 1, 1, 1, 1),
    ):
        row = module.certify_one(101, weights)
        assert row is not None
        assert row["base_quotient_dimension"] == 24
        assert row["successive_endpoint_standard_basis"] == [1]


def test_projective_gate_uses_structural_degree_and_collision() -> None:
    data = archived(ARCHIVE)
    certificate = data["projective_infinity_certificate"]
    assert certificate["structural_degrees"] == [4, 3, 2, 1, 29]
    assert certificate["valid_for_all_selected_primes"] is True
    assert "t=0,s=1" in certificate["lemma"]
    assert "all five black labels" in certificate["lemma"]


def test_complete_CC_and_raw_residual_scope() -> None:
    data = archived(ARCHIVE)
    assert data["CC_black_degree_records_certified"] == 1292
    assert data["new_marked_passports_eliminated"] == 1280
    assert data["partition_black_degree_systems_remaining"] == 1280
    assert data["partition_marked_passports_remaining"] == 3840
    assert data["five_black_status_free_necessary_systems_remaining"] == 9742
    assert data["five_black_ww_marked_passports_remaining"] == 28388
    assert data["five_black_marked_passports_remaining"] == 129788
    assert data["global_marked_passports_eliminated"] == 108602
    assert data["global_marked_passports_remaining"] == 42502076
    assert data["complete_CC_no_go_claimed"] is True
    assert not data["complete_partition_no_go_claimed"]
    assert not data["complete_five_black_ww_no_go_claimed"]
    assert not data["complete_degree_29_no_go_claimed"]


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "Partition-`33` `CC` theorem" in text
    assert "1,292" in text and "1,280" in text
    assert "3,840" in text and "42,502,076" in text
    assert "global_degree_29_status = OPEN" in text
