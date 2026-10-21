from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_gram_belyi_m5_degree29_reduction.py"
ARCHIVE = ROOT / "data/enumerations/gram_belyi_m5_degree29_reduction.json"
REPORT = ROOT / "results/gram_belyi_m5_degree29_reduction.md"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_m5_degree29_reduction", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_archive_is_exact_replay() -> None:
    module = load_module()
    data = archived()
    assert data == module.build_certificate()
    assert data["all_checks_pass"] is True


def test_passport_denominator_and_low_vertex_partition() -> None:
    data = archived()
    census = data["passport_census"]
    assert census["ordinary_passport_pair_count"] == 235_465
    assert census["marked_ramification_passport_count"] == 42_610_678
    low = census["marked_by_lambda_color_vertex_count_and_endpoint_colors"]
    assert low["1"] == {"BB": 0, "BW": 0, "WB": 0, "WW": 1, "total": 1}
    assert low["2"] == {
        "BB": 0,
        "BW": 56,
        "WB": 56,
        "WW": 84,
        "total": 196,
    }
    assert low["3"] == {
        "BB": 756,
        "BW": 1512,
        "WB": 1512,
        "WW": 1372,
        "total": 5152,
    }


def test_every_bounded_trial_has_an_exact_modular_obstruction() -> None:
    data = archived()
    reductions = data["exact_reductions"]
    expected = {
        "two_black_one_endpoint_black": 56,
        "two_black_both_endpoints_white": 28,
        "three_black_both_endpoints_black": 378,
        "three_black_exactly_one_endpoint_black": 756,
        "Galois_symmetric_q_power_ansatz": 210,
    }
    for name, count in expected.items():
        block = reductions[name]
        assert block["trial_count"] == count
        assert len(block["rows"]) == count
        assert sum(block["witness_prime_histogram"].values()) == count
        assert len(block["canonical_rows_sha256"]) == 64
        for row in block["rows"]:
            assert row["prime"] in (101, 109, 131, 139, 149, 151)
            if "gcd" in row:
                assert row["gcd"] == [1]


def test_scope_is_a_reduction_not_a_global_no_go() -> None:
    data = archived()
    theorem = data["theorem"]
    claims = data["claim_boundary"]
    assert theorem["degree_29_global_status"] == "OPEN"
    assert theorem["marked_ramification_passports_eliminated"] == 3_977
    assert theorem["marked_ramification_passports_remaining"] == 42_606_701
    assert theorem["three_black_remaining_marked_passport_count"] == 1_372
    assert theorem["four_or_more_black_remaining_marked_passport_count"] == 42_605_329
    assert not claims["degree_29_witness_found"]
    assert not claims["complete_degree_29_no_go_claimed"]
    assert not claims["three_black_WW_coefficient_systems_solved"]
    assert not claims["four_or_more_black_strata_solved"]
    assert claims["ordinary_passports_are_not_dessins"]
    assert not claims["modular_timeout_used_as_evidence"]


def test_report_scope_and_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "sharpness remains open" in text
    assert "points `x=0,1` are white" in text
    assert "not `42,610,678` dessins" in text
    assert "complete_degree_29_no_go_claimed = false" in text
