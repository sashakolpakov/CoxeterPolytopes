from pathlib import Path
import subprocess
import sys

from scripts.certify_rational_single_chain_first_mixed_ray_bound import (
    build_certificate,
    certify_f2_minus_four_low_a,
    certify_plane_minus_four_low_R,
    f2_canonical,
    f2_first_mixed_ray,
    f2_pair,
    f2_standard_walls,
    plane_canonical,
    plane_first_mixed_ray,
    plane_pair,
    plane_standard_walls,
)


ROOT = Path(__file__).parents[1]
ARCHIVE = (
    ROOT / "data/enumerations/rational_single_chain_first_mixed_ray_bound.json"
)
REPORT = ROOT / "results/rational_single_chain_first_mixed_ray_bound.md"
CERTIFIER = (
    ROOT / "scripts/certify_rational_single_chain_first_mixed_ray_bound.py"
)


def test_plane_first_mixed_ray_is_spacelike_and_checks_every_standard_wall() -> None:
    for m in range(12, 33):
        ray = plane_first_mixed_ray(m)
        pairings = [plane_pair(ray, wall) for wall in plane_standard_walls(m)]
        expected = [0] * (m + 1)
        expected[1] = 1
        expected[12] = 2
        assert pairings == expected
        assert plane_pair(ray, ray) == -4
        assert plane_pair(plane_canonical(m), ray) == 4


def test_f2_first_mixed_ray_is_spacelike_and_checks_every_standard_wall() -> None:
    for m in range(11, 33):
        ray = f2_first_mixed_ray(m)
        pairings = [f2_pair(ray, wall) for wall in f2_standard_walls(m)]
        expected = [0] * (m + 2)
        expected[0] = 1
        expected[12] = 2
        assert pairings == expected
        assert f2_pair(ray, ray) == -4
        assert f2_pair(f2_canonical(m), ray) == 4


def test_plane_pairing_identity_after_adjunction_substitution() -> None:
    samples = (
        (8, 1, (4, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 0)),
        (9, 2, (5, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1)),
        (10, 4, (6, 4, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1)),
    )
    for d, q, multiplicities in samples:
        # Adjust the final tail entry so only adjunction, not geometry, is tested.
        values = list(multiplicities)
        values[-1] += 3 * d + q - 2 - sum(values)
        assert values[-1] >= 0
        m = len(values)
        curve = (d, *(-value for value in values))
        direct = plane_pair(curve, plane_first_mixed_ray(m))
        tail = sum(values[12:])
        closed = d - values[0] - 2 * q + 4 + 2 * tail
        assert direct == closed


def test_f2_pairing_identity_after_adjunction_substitution() -> None:
    samples = (
        (2, 5, 1, (2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        (3, 7, 2, (2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)),
        (4, 9, 4, (2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0)),
    )
    for a, b, q, multiplicities in samples:
        values = list(multiplicities)
        values[-1] += 2 * b + q - 2 - sum(values)
        assert values[-1] >= 0
        m = len(values)
        curve = (a, b, *(-value for value in values))
        direct = f2_pair(curve, f2_first_mixed_ray(m))
        tail = sum(values[11:])
        closed = a - 2 * q + 4 + 2 * tail
        assert direct == closed


def test_parity_free_low_pairing_case_splits_have_no_integer_solution() -> None:
    plane = certify_plane_minus_four_low_R()
    ruled = certify_f2_minus_four_low_a()
    assert plane["tested_R_values"] == [0, 1, 2, 3]
    assert ruled["tested_a_values"] == [0, 1, 2, 3]
    assert plane["solutions"] == 0
    assert ruled["solutions"] == 0
    assert all(row["square_solutions"] == 0 for row in plane["rows"])
    assert all(row["square_solutions"] == 0 for row in ruled["rows"])


def test_scope_is_exact_and_covers_adapted_hosts_without_parity() -> None:
    payload = build_certificate(32)
    assert payload["status"] == "PASS_EXACT_CONDITIONAL_DEPTH_BOUND"
    scope = payload["scope"]
    assert scope["adapted_half_boundary_host_covered_on_the_exact_frames"] is True
    assert scope["arbitrary_single_chain_covered"] is False
    assert scope["reflection_classification_or_census_used"] is False
    conclusion = payload["finite_volume_conclusion"]
    assert conclusion["plane_chain_length_at_most"] == 11
    assert conclusion["F2_chain_length_at_most"] == 10


def test_report_has_no_stale_doubling_or_satellite_claim() -> None:
    report = REPORT.read_text(encoding="utf-8")
    assert "C=2r" not in report
    assert "full-NS integrality gives" not in report
    assert "second centre is satellite" not in report
    assert "fibre-following free infinitely-near point" in report
    assert "never divides" in report
    assert "without a divisibility assumption" in report


def test_archive_replays_exactly() -> None:
    subprocess.run(
        [
            sys.executable,
            str(CERTIFIER),
            "--check",
            str(ARCHIVE),
        ],
        cwd=ROOT,
        check=True,
    )


def test_certifier_is_independent_of_root_enumerators_and_classifications() -> None:
    source = CERTIFIER.read_text(encoding="utf-8")
    assert "certify_rational_i1n_inverse_design" not in source
    assert "enumerate_simple_roots" not in source
    assert "nikulin" not in source.lower()
    assert "vinberg" not in source.lower()
