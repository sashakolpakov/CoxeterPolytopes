from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_high_a.py"
)
QUOTIENT = ROOT / "data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json"
ARCHIVE = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json"
)
REPORT = ROOT / "results/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.md"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("partition33_high_a", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_partition_33_endpoint_multiplier_is_four() -> None:
    quotient = archived(QUOTIENT)
    row = next(
        value
        for value in quotient["white_partitions"]
        if value["partition"] == "(3^2,1^23)"
    )
    assert row["two_ordered_white_endpoint_choices"] == {
        "SS": 1,
        "SC": 1,
        "CS": 1,
        "CC": 1,
    }


def test_archive_covers_exactly_all_a_at_least_21_records() -> None:
    module = load_module()
    data = archived(ARCHIVE)
    expected = [
        tuple(value)
        for value in module.canonical_black_records(0)
        if value[0] >= 21
    ]
    assert len(expected) == len(data["rows"]) == 12
    assert data["prime_histogram"] == {"101": 11, "131": 1}
    actual = []
    for row in data["rows"]:
        record = (row["a_lambda"], row["b"], row["c"], row["d"], row["e"])
        actual.append(record)
        assert record[0] >= 21 and sum(record) == 29
        assert record[1] <= record[2] <= record[3] <= record[4]
        assert row["homogenization_degrees"] == [4, 3, 2, 1, 29, 29, 29]
        assert row["affine_standard_basis"] == [1]
        assert row["singular_stdout"] == "UNIT"
        assert row["uniform_infinity_lemma_applied"] is True
        assert len(row["equation_sha256"]) == 64
    assert actual == expected


def test_fast_representative_certificates_recompute() -> None:
    module = load_module()
    for weights in ((21, 1, 1, 1, 5), (21, 2, 2, 2, 2), (25, 1, 1, 1, 1)):
        row = module.certify_one(101, weights)
        assert row is not None
        assert row["affine_standard_basis"] == [1]


def test_uniform_infinity_metadata_has_both_collision_branches() -> None:
    data = archived(ARCHIVE)
    certificate = data["projective_infinity_certificate"]
    assert certificate["equation_degrees"] == [4, 3, 2, 1, 29, 29, 29]
    assert certificate["valid_for_all_selected_primes"] is True
    assert "For t!=s" in certificate["lemma"]
    assert "For t=s" in certificate["lemma"]


def test_exact_residual_and_open_scope() -> None:
    data = archived(ARCHIVE)
    assert data["black_degree_records_eliminated"] == 12
    assert data["marked_passports_eliminated"] == 48
    assert data["partition_black_degree_systems_remaining"] == 1280
    assert data["partition_marked_passports_remaining"] == 5120
    assert data["five_black_status_free_necessary_systems_remaining"] == 9742
    assert data["five_black_ww_marked_passports_remaining"] == 29668
    assert data["five_black_marked_passports_remaining"] == 131068
    assert data["global_marked_passports_remaining"] == 42503356
    assert data["global_degree_29_status"] == "OPEN"
    assert not data["complete_partition_no_go_claimed"]
    assert not data["complete_five_black_ww_no_go_claimed"]
    assert not data["complete_degree_29_no_go_claimed"]
    assert not data["degree_29_witness_found"]


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "a>=21" in text
    assert "12" in text and "48" in text
    assert "1,280" in text and "5,120" in text
    assert "global_degree_29_status = OPEN" in text
