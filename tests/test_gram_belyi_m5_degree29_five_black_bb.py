from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_gram_belyi_m5_degree29_five_black_bb.py"
QUOTIENT = (
    ROOT / "data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json"
)
ARCHIVE = ROOT / "data/enumerations/gram_belyi_m5_degree29_five_black_bb.json"
REPORT = ROOT / "results/gram_belyi_m5_degree29_five_black_bb.md"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_m5_degree29_five_black_bb", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_five_black_quotient_is_exactly_the_claimed_denominator() -> None:
    data = archived(QUOTIENT)
    assert data["status"] == "EXACT_FIVE_BLACK_PASSPORT_MASTER_QUOTIENT"
    assert data["raw_marked_passports_by_endpoint_colors"] == {
        "BB": 55_055,
        "BW": 50_700,
        "WB": 50_700,
        "WW": 33_592,
        "total": 190_047,
    }
    assert data["binary_status_master_systems_by_endpoint_colors"] == {
        "BB": 11_011,
        "BW": 8_450,
        "WB": 8_450,
        "WW": 5_168,
        "total": 33_079,
    }
    assert data["black_count_crosscheck"]["agree"] is True
    assert not data["five_black_no_go_claimed"]


def test_every_BB_record_has_affine_and_infinity_certificates() -> None:
    module = load_module()
    data = archived(ARCHIVE)
    expected = [tuple(value) for value in module.canonical_black_records(2)]
    assert data["status"] == "EXACT_FIVE_BLACK_BB_NO_GO"
    assert len(expected) == len(data["rows"]) == 11_011
    assert data["prime_histogram"] == {"101": 10_589, "109": 401, "131": 21}
    actual = []
    for row in data["rows"]:
        record = (
            row["a_lambda"],
            row["b_x0"],
            row["c_x1"],
            row["d"],
            row["e"],
        )
        actual.append(record)
        assert sum(record) == 29 and row["d"] <= row["e"]
        assert all(
            degree <= bound
            for degree, bound in zip(
                row["critical_equation_total_degrees"], [27, 27, 26]
            )
        )
        assert row["homogenization_degrees"] == [27, 27, 26]
        assert row["affine_standard_basis"] == [1]
        assert row["singular_stdout"] == "UNIT"
        assert row["infinity_chart_gcd"] == [1]
        assert any(row["infinity_point_u0_v1_values"])
        assert len(row["critical_equation_sha256"]) == 64
    assert actual == expected


def test_representative_projective_certificates_recompute() -> None:
    module = load_module()
    for weights in ((1, 1, 1, 1, 25), (5, 6, 7, 5, 6), (25, 1, 1, 1, 1)):
        row = module.certify_one(101, weights)
        assert row is not None
        assert row["affine_standard_basis"] == [1]
        assert row["infinity_chart_gcd"] == [1]
        assert any(row["infinity_point_u0_v1_values"])


def test_bad_reduction_is_inconclusive_and_next_prime_certifies() -> None:
    module = load_module()
    weights = (7, 16, 2, 1, 3)
    assert module.certify_one(101, weights) is None
    row = module.certify_one(109, weights)
    assert row is not None and row["prime"] == 109


def test_scoped_count_boundary_remains_open() -> None:
    data = archived(ARCHIVE)
    assert data["eliminated_five_black_master_systems"] == 11_011
    assert data["eliminated_five_black_marked_passports"] == 55_055
    assert data["remaining_five_black_master_systems"] == 22_068
    assert data["remaining_five_black_marked_passports"] == 134_992
    assert data["global_marked_passports_eliminated"] == 103_398
    assert data["global_marked_passports_remaining"] == 42_507_280
    assert data["global_degree_29_status"] == "OPEN"
    assert not data["degree_29_witness_found"]
    assert not data["complete_five_black_no_go_claimed"]
    assert not data["complete_degree_29_no_go_claimed"]
    assert not data["broad_plane_tree_search_used"]


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "Five-black `BB` theorem" in text
    assert "11,011" in text
    assert "55,055" in text
    assert "22,068" in text
    assert "global degree-`29` problems remain open" in text
