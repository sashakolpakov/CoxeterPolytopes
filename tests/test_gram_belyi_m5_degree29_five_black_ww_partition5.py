from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition5.py"
)
QUOTIENT = (
    ROOT / "data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition5.json"
)
REPORT = ROOT / "results/gram_belyi_m5_degree29_five_black_ww_partition5.md"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_m5_degree29_five_black_ww_partition5", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_partition_5_has_exactly_three_ordered_endpoint_markings() -> None:
    quotient = archived(QUOTIENT)
    row = next(
        value
        for value in quotient["white_partitions"]
        if value["partition"] == "(5,1^24)"
    )
    assert row["two_ordered_white_endpoint_choices"] == {
        "SS": 1,
        "SC": 1,
        "CS": 1,
        "CC": 0,
    }


def test_every_black_record_has_an_exact_affine_certificate() -> None:
    module = load_module()
    data = archived(ARCHIVE)
    expected = [tuple(value) for value in module.canonical_black_records(0)]
    assert data["status"] == "EXACT_FIVE_BLACK_WW_PARTITION_5_NO_GO"
    assert len(expected) == len(data["rows"]) == 1292
    assert sum(data["prime_histogram"].values()) == 1292
    actual = []
    for row in data["rows"]:
        record = (row["a_lambda"], row["b"], row["c"], row["d"], row["e"])
        actual.append(record)
        assert sum(record) == 29
        assert row["b"] <= row["c"] <= row["d"] <= row["e"]
        assert row["derivative_homogenization_degrees"] == [4, 3, 2, 1]
        assert row["endpoint_homogenization_degrees"] == [29, 29]
        assert row["affine_standard_basis"] == [1]
        assert row["singular_stdout"] == "UNIT"
        assert row["uniform_infinity_lemma_applied"] is True
        assert len(row["equation_sha256"]) == 64
    assert actual == expected


def test_representative_affine_certificates_recompute() -> None:
    module = load_module()
    for weights in (
        (1, 1, 1, 1, 25),
        (5, 6, 6, 6, 6),
        (1, 4, 7, 8, 9),
        (25, 1, 1, 1, 1),
    ):
        row = module.certify_one(101, weights)
        assert row is not None
        assert row["affine_standard_basis"] == [1]
        assert row["uniform_infinity_lemma_applied"] is True


def test_bad_reduction_is_inconclusive_and_next_prime_certifies() -> None:
    module = load_module()
    weights = (1, 2, 6, 6, 14)
    assert module.certify_one(101, weights) is None
    row = module.certify_one(109, weights)
    assert row is not None and row["prime"] == 109


def test_uniform_projective_infinity_metadata_is_scoped() -> None:
    data = archived(ARCHIVE)
    certificate = data["projective_infinity_certificate"]
    assert certificate["homogenized_variables"] == ["u", "v", "w", "z", "t", "s"]
    assert certificate["derivative_equation_degrees"] == [4, 3, 2, 1]
    assert certificate["endpoint_equation_degrees"] == [29, 29]
    assert certificate["valid_for_all_selected_primes"] is True
    assert "ord_t(D_infinity)=4 forces n_t=5" in certificate["lemma"]


def test_scoped_count_boundary_remains_open() -> None:
    data = archived(ARCHIVE)
    assert data["exact_partition_systems_eliminated"] == 1292
    assert data["eliminated_five_black_marked_passports"] == 3876
    assert data["eliminated_five_black_binary_status_master_systems"] == 0
    assert data["remaining_five_black_binary_status_master_systems"] == 22068
    assert data["ww_binary_status_master_systems_before_collapse"] == 5168
    assert data["ww_status_free_master_systems_after_collapse"] == 1292
    assert data["status_free_master_systems_after_collapse_by_endpoint_block"] == {
        "BW": 4225,
        "WB": 4225,
        "WW": 1292,
        "total": 9742,
    }
    assert data["remaining_five_black_status_free_necessary_systems"] == 9742
    assert data["remaining_five_black_ww_marked_passports"] == 29716
    assert data["remaining_five_black_marked_passports"] == 131116
    assert data["global_marked_passports_eliminated"] == 107274
    assert data["global_marked_passports_remaining"] == 42503404
    assert data["global_degree_29_status"] == "OPEN"
    assert not data["degree_29_witness_found"]
    assert not data["complete_five_black_ww_no_go_claimed"]
    assert not data["complete_five_black_no_go_claimed"]
    assert not data["complete_degree_29_no_go_claimed"]
    assert not data["broad_plane_tree_search_used"]


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "partition `(5,1^24)` theorem" in text
    assert "1,292" in text
    assert "3,876" in text
    assert "29,716" in text
    assert "global degree-`29` problems" in text
    assert "remain open" in text
