from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_gram_belyi_m5_degree29_four_black_ww_simple.py"
ARCHIVE = ROOT / "data/enumerations/gram_belyi_m5_degree29_four_black_ww_simple.json"
REPORT = ROOT / "results/gram_belyi_m5_degree29_four_black_ww_simple.md"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_m5_degree29_four_black_ww_simple", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archived() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_all_WW_endpoint_statuses_are_represented() -> None:
    data = archived()
    assert data["status"] == "EXACT_FOUR_BLACK_WW_NO_GO"
    assert data["black_degree_master_system_count"] == 640
    assert data["represented_endpoint_statuses"] == ["SS", "SC", "CS", "CC"]
    assert data["represented_endpoint_status_system_count"] == 2_560
    assert data["represented_marked_passport_count"] == 8_960
    assert data["eliminated_black_degree_master_system_count"] == 640
    assert data["eliminated_endpoint_status_system_count"] == 2_560
    assert data["eliminated_marked_passport_count"] == 8_960
    assert data["remaining_black_degree_master_system_count"] == 0
    assert data["remaining_endpoint_status_system_count"] == 0
    assert data["remaining_marked_passport_count"] == 0


def test_every_resultant_component_and_endpoint_orbit_is_certified() -> None:
    data = archived()
    assert len(data["rows"]) == 640
    assert data["prime_histogram"] == {"101": 572, "109": 58, "131": 10}
    orbit_count = 0
    for row in data["rows"]:
        prime = row["prime"]
        q = row["endpoint_ratio_q_mod_prime"]
        assert (q * q + 12 * q + 16) % prime == 0
        assert row["critical_equation_total_degrees"] == [28, 27]
        assert row["critical_equation_bidegrees"] == [[27, 27], [27, 27]]
        assert row["critical_resultant_degree"] == 756
        resultant = row["critical_resultant_monic_coefficients"]
        assert len(resultant) == 757 and resultant[-1] == 1
        assert row["u_factor_count"] == len(row["components"])
        for component in row["components"]:
            assert max(component["specialized_R1_R2_degrees"]) == 27
            assert sum(
                degree * exponent
                for degree, exponent in component["v_factor_degrees"]
            ) == component["common_v_degree"]
            assert all(degree == 0 for degree in component["endpoint_gcd_degrees"])
            assert len(component["v_factor_degrees"]) == len(
                component["endpoint_gcd_degrees"]
            )
            orbit_count += len(component["endpoint_gcd_degrees"])
        assert row["all_endpoint_gcds_are_one"] is True
    assert orbit_count == data["total_endpoint_orbits_checked"] == 10_696


def test_canonical_remainder_matches_direct_specialization() -> None:
    module = load_module()
    prime = 101
    _sqrt5, lam = module.lambda_mod(prime)
    ring = module.FpPoly(prime)

    for weights in ((1, 1, 1, 26), (7, 8, 9, 5), (26, 1, 1, 1)):
        _context, remainder = module.canonical_remainder(prime, lam, weights)
        for u_value, v_value in ((2, 3), (17, 73)):
            evaluated = []
            for coefficient in remainder:
                value = 0
                for (u_power, v_power), scalar in coefficient.terms():
                    value += (
                        int(scalar)
                        * pow(u_value, u_power, prime)
                        * pow(v_value, v_power, prime)
                    )
                evaluated.append(value % prime)

            roots = (0, 1, u_value, v_value)
            factors = [[-root % prime, 1] for root in roots]
            product_value = [1]
            for factor, exponent in zip(factors, weights):
                product_value = ring.mul(
                    product_value, ring.power(factor, exponent)
                )
            derivative_residual = [0]
            for omitted, weight in enumerate(weights):
                term = [1]
                for index, factor in enumerate(factors):
                    if index != omitted:
                        term = ring.mul(term, factor)
                derivative_residual = ring.add(
                    derivative_residual, ring.scale(term, weight)
                )
            _quotient, direct = ring.divmod(product_value, derivative_residual)
            direct += [0] * (3 - len(direct))
            assert evaluated == direct


def test_representative_tower_certificate_recomputes() -> None:
    module = load_module()
    row = module.certify_one(101, (1, 1, 1, 26))
    assert row is not None
    assert row["critical_resultant_degree"] == 756
    assert row["endpoint_orbit_count"] == 37
    assert all(
        degree == 0
        for component in row["components"]
        for degree in component["endpoint_gcd_degrees"]
    )


def test_updated_claim_boundary_remains_global_open() -> None:
    data = archived()
    assert data["four_black_coefficient_systems_remaining"] == 3_458
    assert data["four_black_marked_passports_remaining"] == 10_374
    assert data["global_marked_passports_eliminated"] == 37_969
    assert data["global_marked_passports_remaining"] == 42_572_709
    assert data["global_degree_29_status"] == "OPEN"
    assert not data["degree_29_witness_found"]
    assert not data["complete_four_black_no_go_claimed"]
    assert not data["complete_degree_29_no_go_claimed"]
    assert not data["broad_plane_tree_search_used"]


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "Four-black `WW` theorem" in text
    assert "10,696" in text
    assert "8,960" in text
    assert "global degree-`29` sharpness remains open" in text
    assert "only remaining" in text
