from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_gram_belyi_m5_degree29_four_black_reduction.py"
ARCHIVE = ROOT / "data/enumerations/gram_belyi_m5_degree29_four_black_reduction.json"
REPORT = ROOT / "results/gram_belyi_m5_degree29_four_black_reduction.md"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_m5_degree29_four_black_reduction", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_exact_combinatorial_quotient() -> None:
    module = load_module()
    assert len(module.four_black_degree_records()) == 640
    assert len(module.one_black_endpoint_records()) == 1_729
    assert len(module.both_black_endpoint_records()) == 3_276

    data = archived()
    quotient = data["passport_system_quotient"]
    assert quotient["BB"] == {"marked_passports": 9_828, "coefficient_systems": 3_276}
    assert quotient["BW"] == {"marked_passports": 12_103, "coefficient_systems": 3_458}
    assert quotient["WB"] == {"marked_passports": 12_103, "coefficient_systems": 3_458}
    assert quotient["WW"] == {"marked_passports": 8_960, "coefficient_systems": 2_560}
    assert sum(row["marked_passports"] for row in (quotient["BB"], quotient["BW"], quotient["WB"], quotient["WW"])) == 42_994
    assert quotient["at_least_one_white_endpoint_critical"] == {
        "coefficient_systems": 5_378,
        "marked_passports": 20_872,
    }


def test_every_frozen_resultant_has_full_degree_and_gcd_one() -> None:
    data = archived()
    bb = data["eliminated_blocks"]["both_endpoints_black"]
    critical = data["eliminated_blocks"]["one_black_other_white_critical"]
    assert bb["coefficient_system_count"] == 3_276
    assert bb["marked_passport_count"] == 9_828
    assert bb["prime_histogram"] == {"101": 3_156, "109": 106, "131": 14}
    assert critical["coefficient_system_count"] == 3_458
    assert critical["marked_passport_count"] == 13_832
    assert critical["prime_histogram"] == {"101": 3_285, "109": 164, "131": 9}
    for row in bb["rows"]:
        assert row["remainder_degrees_R1_R2"] == [27, 27]
        assert row["gcd"] == [1]
    for row in critical["rows"]:
        assert row["cleared_remainder_degrees_R1_R2"] == [55, 54]
        assert row["exceptional_branch_resultant_mod_prime"] != 0
        assert row["gcd"] == [1]


def test_parameter_remainder_matches_direct_specialization() -> None:
    module = load_module()
    prime = 101
    _sqrt5, lam = module.lambda_mod(prime)
    ring = module.FpPoly(prime)
    for weights in ((1, 1, 1, 26), (7, 8, 9, 5), (26, 1, 1, 1)):
        parameter = module.parameter_remainder(
            prime, ([lam], [0], [1], [0, 1]), weights
        )
        for t_value in (2, 17, 73):
            direct = module.scalar_remainder(
                prime, (lam, 0, 1, t_value), weights
            )
            assert [ring.evaluate(coefficient, t_value) for coefficient in parameter] == direct


def test_interpolated_critical_equations_are_exact_on_all_residues() -> None:
    module = load_module()
    prime = 101
    _sqrt5, lam = module.lambda_mod(prime)
    ring = module.FpPoly(prime)
    for alpha, weights in (
        (lam, (1, 1, 1, 26)),
        ((1 - lam) % prime, (7, 8, 9, 5)),
    ):
        first, second, exceptional = (
            module.critical_endpoint_equations_by_interpolation(
                prime, alpha, weights
            )
        )
        assert [len(first) - 1, len(second) - 1] == [55, 54]
        assert exceptional != 0
        assert ring.gcd(first, second) == [1]


def test_exact_remainder_and_global_claim_boundary() -> None:
    data = archived()
    assert data["status"] == "EXACT_QUANTIFIED_FOUR_BLACK_REDUCTION"
    assert data["global_degree_29_status"] == "OPEN"
    assert data["eliminated_four_black_coefficient_systems"] == 6_734
    assert data["eliminated_four_black_marked_passports"] == 23_660
    assert data["remaining_four_black_coefficient_systems"] == 6_018
    assert data["remaining_four_black_marked_passports"] == 19_334
    assert data["global_marked_passports_eliminated_through_four_black"] == 29_009
    assert data["global_marked_passports_remaining"] == 42_581_669
    remaining = data["remaining_blocks"]
    assert remaining == {
        "one_black_other_white_simple": {
            "coefficient_systems": 3_458,
            "marked_passports": 10_374,
        },
        "both_white_at_least_one_critical": {
            "coefficient_systems": 1_920,
            "marked_passports": 7_040,
        },
        "both_white_simple": {
            "coefficient_systems": 640,
            "marked_passports": 1_920,
        },
    }
    assert not data["degree_29_witness_found"]
    assert not data["complete_four_black_no_go_claimed"]
    assert not data["complete_degree_29_no_go_claimed"]
    assert not data["broad_plane_tree_search_used"]


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "Four-black partial theorem" in text
    assert "23,660" in text
    assert "6,018" in text
    assert "sharpness remains open" in text
    assert "not a claim that every master system constructs" in text
