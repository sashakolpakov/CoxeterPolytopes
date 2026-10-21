from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_gram_belyi_m5_degree29_three_black_ww.py"
ARCHIVE = ROOT / "data/enumerations/gram_belyi_m5_degree29_three_black_ww.json"
REPORT = ROOT / "results/gram_belyi_m5_degree29_three_black_ww.md"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_m5_degree29_three_black_ww", SCRIPT
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
    assert data == module.build_archive()
    module.verify_archive(data)


def test_census_collapses_1372_records_to_784_systems() -> None:
    data = archived()
    assert data["marked_black_degree_records"] == 196
    assert data["marked_three_black_WW_records_before"] == 1_372
    assert data["distinct_three_black_WW_coefficient_systems"] == 784
    assert data["distinct_endpoint_critical_systems"] == 588
    assert data["endpoint_critical_marked_records"] == 980
    assert data["distinct_simple_simple_systems"] == 196
    assert data["simple_simple_marked_records"] == 392


def test_all_endpoint_critical_resultants_are_separated() -> None:
    data = archived()
    both = data["both_endpoints_critical"]
    one = data["exactly_one_endpoint_required_critical"]
    assert len(both) == 196
    assert len(one) == 392
    assert data["endpoint_critical_prime_histogram"] == {"101": 525, "109": 63}
    for row in both:
        assert row["equation_degrees"] == [2, 29 - row["a"]]
        assert row["gcd"] == [1]
    for row in one:
        assert row["equation_degrees"] == [29 - row["a"], 58]
        assert row["exceptional_branch_resultant_mod_prime"] != 0
        assert row["gcd"] == [1]


def test_all_simple_simple_extension_gcds_are_one() -> None:
    module = load_module()
    data = archived()
    rows = data["simple_simple"]
    assert len(rows) == 196
    assert data["simple_simple_prime_histogram"] == {
        "101": 173,
        "109": 22,
        "131": 1,
    }
    for row in rows:
        prime = row["prime"]
        q = row["endpoint_ratio_q_mod_prime"]
        assert (q * q + 12 * q + 16) % prime == 0
        assert row["shabat_equation_degree"] == 28
        assert len(row["shabat_equation_coefficients_mod_prime"]) == 29
        factors = row["irreducible_factor_checks"]
        assert sum(
            factor["factor_degree"] * factor["factor_exponent"]
            for factor in factors
        ) == 28
        assert all(factor["H_degrees"] == [29, 29] for factor in factors)
        assert all(factor["gcd_degree"] == 0 for factor in factors)
        assert row["all_extension_gcds_are_one"] is True

        ring = module.FpPoly(prime)
        reconstructed = [1]
        for factor in factors:
            reconstructed = ring.mul(
                reconstructed,
                ring.power(
                    factor["factor_coefficients"], factor["factor_exponent"]
                ),
            )
        shabat = row["shabat_equation_coefficients_mod_prime"]
        shabat_monic = ring.scale(shabat, pow(shabat[-1], -1, prime))
        assert reconstructed == shabat_monic


def test_elimination_formulas_against_direct_specialization() -> None:
    """Check the cleared formulas independently at ordinary residue points."""

    module = load_module()
    prime = 101
    _sqrt5, lam = module.lambda_mod(prime)
    ring = module.FpPoly(prime)
    inverse29 = pow(29, -1, prime)

    def f_scalar(x_value: int, r_value: int, s_value: int, a: int, b: int, c: int) -> int:
        return (
            pow((x_value - lam) % prime, a, prime)
            * pow((x_value - r_value) % prime, b, prime)
            * pow((x_value - s_value) % prime, c, prime)
        ) % prime

    for a, b, c in ((1, 1, 27), (7, 8, 14), (19, 5, 5)):
        # Both endpoints critical: the first returned equation is
        # (a+b)B and the second is (a+b)^c(F(1)-F(0)).
        first, second, _bounds = module.both_endpoints_critical_equations(
            prime, lam, a, b, c
        )
        denominator = (a + b) % prime
        for r_value in (2, 17):
            s_num = (29 - (b + c) * lam - (a + c) * r_value) % prime
            s_value = s_num * pow(denominator, -1, prime) % prime
            A = (
                (b + c) * lam + (a + c) * r_value + (a + b) * s_value
            ) % prime
            B = (
                a * r_value * s_value
                + b * lam * s_value
                + c * lam * r_value
            ) % prime
            assert A == 29
            assert ring.evaluate(first, r_value) == denominator * B % prime
            assert ring.evaluate(second, r_value) == (
                pow(denominator, c, prime)
                * (
                    f_scalar(1, r_value, s_value, a, b, c)
                    - f_scalar(0, r_value, s_value, a, b, c)
                )
            ) % prime

        # One endpoint critical: check the eliminated root, the other
        # critical point, and both cleared value equations directly.
        for endpoint in (0, 1):
            first, second, c_poly, _bounds, _exceptional = (
                module.endpoint_critical_equations(
                    prime, lam, a, b, c, endpoint
                )
            )
            for r_value in (2, 17):
                C = ring.evaluate(c_poly, r_value)
                if C == 0:
                    continue
                if endpoint == 0:
                    K = c * lam * r_value % prime
                else:
                    K = (
                        29
                        - (b + c) * lam
                        - (a + c) * r_value
                        + c * lam * r_value
                    ) % prime
                s_value = -K * pow(C, -1, prime) % prime
                A = (
                    (b + c) * lam
                    + (a + c) * r_value
                    + (a + b) * s_value
                ) % prime
                B = (
                    a * r_value * s_value
                    + b * lam * s_value
                    + c * lam * r_value
                ) % prime
                assert (29 * endpoint * endpoint - A * endpoint + B) % prime == 0
                z = (A * inverse29 - endpoint) % prime
                f_endpoint = f_scalar(endpoint, r_value, s_value, a, b, c)
                assert ring.evaluate(first, r_value) == (
                    pow(C, c, prime)
                    * (
                        f_scalar(1 - endpoint, r_value, s_value, a, b, c)
                        - f_endpoint
                    )
                ) % prime
                assert ring.evaluate(second, r_value) == (
                    pow(29 * C, 29, prime)
                    * (f_scalar(z, r_value, s_value, a, b, c) - f_endpoint)
                ) % prime

        # Canonical R_0,R_1 agree with direct division after specializing t.
        r0, r1 = module.canonical_critical_remainder(prime, a, b, c)
        for t_value in (2, 17):
            polynomial = ring.mul(
                ring.mul(
                    ring.power([0, 1], a), ring.power([-1, 1], b)
                ),
                ring.power([-t_value, 1], c),
            )
            A = (a + c + (a + b) * t_value) % prime
            D = [a * t_value, -A, 29]
            _quotient, remainder = ring.divmod(polynomial, D)
            expected = ring.trim(
                [ring.evaluate(r0, t_value), ring.evaluate(r1, t_value)]
            )
            assert remainder == expected


def test_claim_boundary_remains_global_open() -> None:
    data = archived()
    assert data["status"] == "EXACT_THREE_BLACK_WW_NO_GO"
    assert data["global_degree_29_status"] == "OPEN"
    assert data["eliminated_three_black_WW_marked_records"] == 1_372
    assert data["remaining_three_black_WW_marked_records"] == 0
    assert (
        data["remaining_degree_29_marked_records_with_at_least_four_black_vertices"]
        == 42_605_329
    )
    assert not data["degree_29_witness_found"]
    assert not data["complete_degree_29_no_go_claimed"]
    assert data["unresolved_three_black_WW_systems"] == []


def test_report_scope_and_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "Three-black theorem" in text
    assert "all `196` equations pass" in text
    assert "global degree-`29` sharpness remains open" in text
    assert "not a complete" in text
