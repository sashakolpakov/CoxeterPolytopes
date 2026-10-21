"""Independent hostile checks for the first-mixed-ray depth theorem.

The arithmetic below does not import the producer certifier.  In particular,
the primitive square-minus-four case is checked without a divisibility or
full-Neron--Severi hypothesis.
"""

from __future__ import annotations

import hashlib
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "rational_single_chain_first_mixed_ray_bound_hostile_audit.md"
)
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "rational_single_chain_first_mixed_ray_bound.json"
)
CERTIFIER = (
    ROOT
    / "scripts"
    / "certify_rational_single_chain_first_mixed_ray_bound.py"
)
SOURCE_HASHES = {
    "results/rational_single_chain_first_mixed_ray_bound.md":
        "57bef687ef27540c92fe277e722200f6906391a2584dbd62c8a680e1367bdbd1",
    "scripts/certify_rational_single_chain_first_mixed_ray_bound.py":
        "f9c4c4bcaa9dfb0d758daff31e6ae0ad6399614d11b047cd8ca571c834564abb",
    "data/enumerations/rational_single_chain_first_mixed_ray_bound.json":
        "25e0ba8392719b68e4e7e0d328bd6795419c48de620d7e84e796749e2acee078",
    "tests/test_rational_single_chain_first_mixed_ray_bound.py":
        "db8e9edcb8c8af5edb44c0ed11f98ec37d93c86ddf58e0f22191d6d5559c435c",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def p2_pair(left, right):
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def p2_mixed_ray(m: int):
    return (7, -3, *((-2,) * 11), *((0,) * (m - 12)))


def p2_standard_walls(m: int):
    walls = [(1, -1, -1, -1, *((0,) * (m - 3)))]
    for index in range(1, m):
        wall = [0] * (m + 1)
        wall[index], wall[index + 1] = 1, -1
        walls.append(tuple(wall))
    endpoint = [0] * (m + 1)
    endpoint[m] = 1
    walls.append(tuple(endpoint))
    return tuple(walls)


def f2_pair(left, right):
    return (
        -2 * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def f2_mixed_ray(m: int):
    return (4, 9, *((-2,) * 11), *((0,) * (m - 11)))


def f2_standard_walls(m: int):
    section = (1, 0, *((0,) * m))
    fibre = (0, 1, -1, -1, *((0,) * (m - 2)))
    walls = [section, fibre]
    for index in range(1, m):
        wall = [0] * (m + 2)
        wall[index + 1], wall[index + 2] = 1, -1
        walls.append(tuple(wall))
    endpoint = [0] * (m + 2)
    endpoint[-1] = 1
    walls.append(tuple(endpoint))
    return tuple(walls)


def monotone_partitions(length, cap, total, square_total):
    """Yield nonincreasing nonnegative tuples with prescribed first moments."""
    if length == 0:
        if total == 0 and square_total == 0:
            yield ()
        return
    if total < 0 or square_total < 0:
        return
    upper = min(cap, total, isqrt(square_total))
    for value in range(upper, -1, -1):
        for tail in monotone_partitions(
            length - 1,
            value,
            total - value,
            square_total - value * value,
        ):
            yield (value, *tail)


def test_both_mixed_rays_are_spacelike_and_satisfy_every_standard_wall():
    for m in range(12, 41):
        ray = p2_mixed_ray(m)
        pairings = [p2_pair(ray, wall) for wall in p2_standard_walls(m)]
        assert p2_pair(ray, ray) == -4
        assert p2_pair((-3, *((1,) * m)), ray) == 4
        assert pairings[1] == 1
        assert pairings[12] == 2
        assert all(value >= 0 for value in pairings)
        assert sum(value != 0 for value in pairings) == 2

    for m in range(11, 41):
        ray = f2_mixed_ray(m)
        pairings = [f2_pair(ray, wall) for wall in f2_standard_walls(m)]
        assert f2_pair(ray, ray) == -4
        assert f2_pair((-2, -4, *((1,) * m)), ray) == 4
        assert pairings[0] == 1
        assert pairings[12] == 2
        assert all(value >= 0 for value in pairings)
        assert sum(value != 0 for value in pairings) == 2


def test_plane_primitive_minus_four_T_zero_bad_R_cases_are_exhausted():
    # For T=sum_{i>12} mu_i=0 and R=d-mu_1, adjunction and C^2=-4 give
    #   sum_{2..12} mu_i   = 2d+R+2,
    #   sum_{2..12} mu_i^2 = 2dR-R^2+4.
    # A negative pairing requires R<=3.  The line wall additionally gives
    # mu_2+mu_3<=R.  The finite search below is exhaustive because the
    # latter inequality bounds every remaining multiplicity by R and hence
    # 2d+R+2<=11R, so d<=14 for R<=3.
    survivors = []
    for R in range(4):
        for degree in range(1, 15):
            mu_1 = degree - R
            if mu_1 < 0:
                continue
            for tail in monotone_partitions(
                11,
                mu_1,
                2 * degree + R + 2,
                2 * degree * R - R * R + 4,
            ):
                if tail[0] + tail[1] <= R:
                    survivors.append((R, degree, mu_1, tail))
    assert survivors == []


def test_F2_primitive_minus_four_T_zero_bad_a_cases_are_exhausted():
    # With T=sum_{i>11} mu_i=0, adjunction and C^2=-4 give the two targets
    # below.  A negative pairing requires a<=3.  The fibre wall bounds every
    # multiplicity by a, so the finite b range is exhaustive.
    survivors = []
    for a in range(4):
        for b in range(2 * a, 17):
            for multiplicities in monotone_partitions(
                11,
                a,
                2 * b + 2,
                -2 * a * a + 2 * a * b + 4,
            ):
                if multiplicities[0] + multiplicities[1] <= a:
                    survivors.append((a, b, multiplicities))
    assert survivors == []


def test_tail_sum_one_forces_exact_nonnegative_boundary_values():
    # A nonincreasing integral sequence with positive tail sum one has a 1
    # in the first tail position and hence 1 in every earlier position.
    plane = (1,) * 13 + (0, 0)
    assert sum(plane[12:]) == 1
    assert plane[1] + plane[2] == 2
    assert (2 - 4 + 2 * sum(plane[12:])) == 0

    ruled = (1,) * 12 + (0, 0)
    assert sum(ruled[11:]) == 1
    assert ruled[0] + ruled[1] == 2
    assert (2 - 4 + 2 * sum(ruled[11:])) == 0


def test_scope_hypotheses_have_exact_formal_countermodels_when_removed():
    m = 12
    w_plane = p2_mixed_ray(m)
    k_plane = (-3, *((1,) * m))
    alpha_0 = p2_standard_walls(m)[0]

    # A primitive conic with eight successive simple multiplicities has the
    # correct rational (-4) numerics and negative mixed pairing, but meets
    # the marked line negatively.  This is numerical scope evidence, not an
    # assertion that the two curves coexist geometrically.
    no_line = (2, *((-1,) * 8), *((0,) * 4))
    assert p2_pair(no_line, no_line) == -4
    assert p2_pair(k_plane, no_line) == 2
    assert p2_pair(w_plane, no_line) == -3
    assert p2_pair(alpha_0, no_line) == -1

    # Moving one simple multiplicity past a zero repairs line
    # nonintersection but violates the exceptional-chain proximity order.
    mu = [1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]
    no_proximity = (2, *[-value for value in mu])
    assert p2_pair(no_proximity, no_proximity) == -4
    assert p2_pair(k_plane, no_proximity) == 2
    assert p2_pair(w_plane, no_proximity) == -3
    assert p2_pair(alpha_0, no_proximity) == 0
    assert mu[1] < mu[2]

    m = 11
    w_ruled = f2_mixed_ray(m)
    k_ruled = (-2, -4, *((1,) * m))
    section, fibre_wall = f2_standard_walls(m)[:2]

    # This class respects the section and proximity inequalities but not
    # the marked fibre wall.
    no_fibre = (1, 2, *((-1,) * 6), *((0,) * 5))
    assert f2_pair(no_fibre, no_fibre) == -4
    assert f2_pair(k_ruled, no_fibre) == 2
    assert f2_pair(w_ruled, no_fibre) == -3
    assert f2_pair(section, no_fibre) == 0
    assert f2_pair(fibre_wall, no_fibre) == -1

    # This one respects both global walls but violates monotonicity at the
    # second exceptional place.
    mu = [1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    no_proximity = (1, 2, *[-value for value in mu])
    assert f2_pair(no_proximity, no_proximity) == -4
    assert f2_pair(k_ruled, no_proximity) == 2
    assert f2_pair(w_ruled, no_proximity) == -3
    assert f2_pair(section, no_proximity) == 0
    assert f2_pair(fibre_wall, no_proximity) == 0
    assert mu[1] < mu[2]


def test_frozen_producer_bytes_archive_and_scope_replay_exactly():
    assert {
        relative: sha256(ROOT / relative)
        for relative in SOURCE_HASHES
    } == SOURCE_HASHES
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert payload["status"] == "PASS_EXACT_CONDITIONAL_DEPTH_BOUND"
    scope = payload["scope"]
    assert scope["adapted_half_boundary_host_covered_on_the_exact_frames"] is True
    assert scope["arbitrary_single_chain_covered"] is False
    assert scope["reflection_classification_or_census_used"] is False
    conclusion = payload["finite_volume_conclusion"]
    assert conclusion["plane_chain_length_at_most"] == 11
    assert conclusion["F2_chain_length_at_most"] == 10


def test_hostile_report_records_the_proof_and_keeps_the_scope_exact():
    report = AUDIT.read_text(encoding="utf-8")
    required = [
        "**PASS.**",
        r"\(m\leq11\)",
        r"\(m\leq10\)",
        r"\boxed{\dim P\leq11}",
        "parity-free",
        "adapted",
        "interior null direction",
        "open ideal-boundary patch",
        "fibre-following free infinitely-near",
        "not satellite",
        "ten ones",
        "formal countermodels",
        "does not force either marked frame",
        "not an existence or sharpness theorem",
    ]
    for phrase in required:
        assert phrase in report

    forbidden_controls = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden_controls for byte in AUDIT.read_bytes())

    producer = (
        ROOT / "results/rational_single_chain_first_mixed_ray_bound.md"
    ).read_text(encoding="utf-8")
    assert "fibre-following free infinitely-near point" in producer
    assert "so the point is not satellite" in producer
    assert "never divides" in producer
    assert "without a divisibility assumption" in producer
    assert "C=2r" not in producer
    assert "full-NS integrality gives" not in producer
