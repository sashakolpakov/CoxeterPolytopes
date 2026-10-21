from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/certify_gram_belyi_m5_degree29_four_black_one_black_simple.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_four_black_one_black_simple.json"
)
REPORT = (
    ROOT
    / "results/gram_belyi_m5_degree29_four_black_one_black_simple.md"
)


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_m5_degree29_four_black_one_black_simple", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_final_four_black_block_is_completely_represented() -> None:
    data = archived()
    assert data["status"] == "EXACT_FOUR_BLACK_COMPLETE_NO_GO"
    assert data["black_endpoint_degree_record_count"] == 1_729
    assert data["endpoint_orientations"] == ["BW", "WB"]
    assert data["broad_orientation_system_count"] == 3_458
    assert data["remaining_simple_endpoint_status_system_count_before"] == 3_458
    assert data["remaining_simple_marked_passport_count_before"] == 10_374
    assert data["eliminated_simple_endpoint_status_system_count"] == 3_458
    assert data["eliminated_simple_marked_passport_count"] == 10_374


def test_every_resultant_component_and_endpoint_orbit_is_certified() -> None:
    data = archived()
    assert len(data["rows"]) == 1_729
    assert data["prime_histogram"] == {
        "101": 1_529,
        "109": 174,
        "131": 23,
        "139": 2,
        "149": 1,
    }
    orbit_count = 0
    for row in data["rows"]:
        prime = row["prime"]
        lam = row["lambda_mod_prime"]
        endpoints = row["endpoint_coordinates_mod_prime"]
        assert (29 * lam * lam - 14 * lam + 1) % prime == 0
        assert endpoints["BW"] * lam % prime == 1
        assert endpoints["WB"] * (1 - lam) % prime == 1
        assert row["critical_equation_total_degrees"] == [28, 27]
        assert row["critical_equation_bidegrees"] == [[27, 27], [27, 27]]
        assert row["critical_resultant_degree"] == 756
        resultant = row["critical_resultant_monic_coefficients"]
        assert len(resultant) == 757 and resultant[-1] == 1
        assert row["u_factor_count"] == len(row["components"])
        assert sum(
            component["u_factor_degree"] * component["u_factor_exponent"]
            for component in row["components"]
        ) == 756
        for component in row["components"]:
            assert max(component["specialized_R1_R2_degrees"]) == 27
            assert sum(
                degree * exponent
                for degree, exponent in component["v_factor_degrees"]
            ) == component["common_v_degree"]
            checks = component["endpoint_nonzero_BW_WB"]
            assert len(checks) == len(component["v_factor_degrees"])
            assert all(value == [True, True] for value in checks)
            orbit_count += len(checks)
        assert row["all_BW_WB_endpoint_differences_nonzero"] is True
    assert orbit_count == data["total_critical_orbits_checked"] == 28_159


def test_canonical_degree_records_and_orientation_count_recompute() -> None:
    module = load_module()
    records = list(module.one_black_endpoint_records())
    assert len(records) == 1_729
    assert len(set(records)) == 1_729
    assert all(sum(record) == 29 and record[2] <= record[3] for record in records)
    assert 2 * len(records) == 3_458


def test_representative_tower_certificate_recomputes_both_orientations() -> None:
    module = load_module()
    for degrees in ((1, 1, 1, 26), (1, 2, 8, 18), (26, 1, 1, 1)):
        row = module.certify_one(101, degrees)
        assert row is not None
        assert row["critical_resultant_degree"] == 756
        assert row["all_BW_WB_endpoint_differences_nonzero"] is True
        assert all(
            check == [True, True]
            for component in row["components"]
            for check in component["endpoint_nonzero_BW_WB"]
        )


def test_four_black_completion_and_global_claim_boundary() -> None:
    data = archived()
    assert data["four_black_coefficient_systems_remaining"] == 0
    assert data["four_black_marked_passports_remaining"] == 0
    assert data["global_marked_passports_eliminated"] == 48_343
    assert data["global_marked_passports_remaining"] == 42_562_335
    assert data["global_degree_29_status"] == "OPEN"
    assert not data["degree_29_witness_found"]
    assert data["complete_four_black_no_go_claimed"]
    assert not data["complete_degree_29_no_go_claimed"]
    assert not data["broad_plane_tree_search_used"]


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "Four-black completion theorem" in text
    assert "28,159" in text
    assert "at least five distinct black" in text
    assert "42,562,335" in text
    assert "global degree-`29` sharpness remains open" in text
