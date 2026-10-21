"""Independent hostile checks for the marked plane single-chain gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "rational_single_chain_standard_simplex_gate_hostile_audit.md"
)
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "rational_single_chain_standard_simplex_gate.json"
)
SOURCE_HASHES = {
    "results/rational_single_chain_standard_simplex_gate.md":
        "ccf0eac1ed58d7c92ae76043780cf6492e1ef261e7785470d55c330a45e4a5f2",
    "scripts/certify_rational_single_chain_standard_simplex_gate.py":
        "e40d400819c5d3f1d2499c325b3ee34cc7e3c7b870691c9d1cf798fdb44fa430",
    "data/enumerations/rational_single_chain_standard_simplex_gate.json":
        "bce15a2ad61ee32591aa912321c12b17bdf3e60f4a277ac53999dbbd9e780f38",
    "tests/test_rational_single_chain_standard_simplex_gate.py":
        "93a843b8b39d31ba26e6d82ed598eb1cd12b8c2781be0dcab0ad931d7cc2b832",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pairing(left, right):
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def standard_roots(m: int):
    alpha_0 = (1, -1, -1, -1, *((0,) * (m - 3)))
    roots = [alpha_0]
    for index in range(1, m):
        root = [0] * (m + 1)
        root[index], root[index + 1] = 1, -1
        roots.append(tuple(root))
    last = [0] * (m + 1)
    last[m] = 1
    roots.append(tuple(last))
    return tuple(roots)


def closed_form_duals(m: int):
    h = (1, *((0,) * m))
    rays = [h]
    for index in range(1, m + 1):
        ray = [min(index, 3)] + [0] * m
        ray[1:index + 1] = [-1] * index
        rays.append(tuple(ray))
    return tuple(rays)


def canonical(m: int):
    return (-3, *((1,) * m))


def f2_pairing(left, right):
    return (
        -2 * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def f2_roots(m: int):
    section = (1, 0, *((0,) * m))
    fibre_root = (0, 1, -1, -1, *((0,) * (m - 2)))
    roots = [section, fibre_root]
    for index in range(1, m):
        root = [0] * (m + 2)
        root[index + 1], root[index + 2] = 1, -1
        roots.append(tuple(root))
    endpoint = [0] * (m + 2)
    endpoint[-1] = 1
    roots.append(tuple(endpoint))
    return tuple(roots)


def f2_closed_form_duals(m: int):
    fibre = (0, 1, *((0,) * m))
    section_dual = (1, 2, *((0,) * m))
    rays = [fibre, section_dual]
    for index in range(1, m + 1):
        coefficient = min(index, 2)
        ray = [coefficient, 2 * coefficient, *([0] * m)]
        ray[2:index + 2] = [-1] * index
        rays.append(tuple(ray))
    return tuple(rays)


def f2_canonical(m: int):
    return (-2, -4, *((1,) * m))


def test_producer_pins_and_archive_scope_are_stable():
    assert {
        relative: sha256(ROOT / relative)
        for relative in SOURCE_HASHES
    } == SOURCE_HASHES
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert payload["status"] == "PASS_EXACT_SCOPED_NO_GO"
    scope = payload["scope"]
    assert scope["alpha_0_not_universal"] is True
    assert scope["geometric_effectivity_of_doubled_numerical_wall_claimed"] is False
    assert scope["uniform_chain_length_bound_claimed"] is False
    assert scope["classification_or_census_used"] is False


def test_dual_basis_inverse_diagonal_and_visible_threshold_independently():
    finite = []
    for m in range(3, 41):
        roots = standard_roots(m)
        rays = closed_form_duals(m)
        assert [
            [pairing(ray, root) for root in roots]
            for ray in rays
        ] == [
            [int(row == column) for column in range(m + 1)]
            for row in range(m + 1)
        ]
        squares = [pairing(ray, ray) for ray in rays]
        assert squares == [
            1,
            *[min(index, 3) ** 2 - index for index in range(1, m + 1)],
        ]
        if all(square >= 0 for square in squares):
            finite.append(m)
    assert finite == list(range(3, 10))


def test_tail_arithmetic_genus_and_physical_adjunction_identity():
    for m in range(10, 33):
        rays = closed_form_duals(m)
        k_class = canonical(m)
        for j in range(10, m + 1):
            ray = rays[j]
            assert pairing(ray, ray) == 9 - j
            assert pairing(k_class, ray) == j - 9
            assert 1 + (
                pairing(ray, ray) + pairing(k_class, ray)
            ) // 2 == 1

            # The identity before imposing smooth-rational adjunction is
            # C.v_j=(-K).C+sum_{i>j} mu_i.
            multiplicities = [((2 * i + j) % 4) for i in range(1, m + 1)]
            curve = (7, *[-value for value in multiplicities])
            tail_sum = sum(multiplicities[j:])
            assert pairing(curve, ray) == (
                -pairing(k_class, curve) + tail_sum
            )

    # Three numerical smooth-rational adjunction classes reproduce the sign
    # table for the complete (-1),(-2),(-4) alphabet at j=10.
    m = 10
    v_10 = closed_form_duals(m)[10]
    examples = (
        (-1, (1, -1, -1, *((0,) * 8)), 1),
        (-2, (1, -1, -1, -1, *((0,) * 7)), 0),
        (-4, (6, *((-2,) * 10)), -2),
    )
    for expected_square, curve, expected_pairing in examples:
        assert pairing(curve, curve) == expected_square
        assert pairing(curve, curve) + pairing(canonical(m), curve) == -2
        assert pairing(curve, v_10) == expected_pairing


def reflection_is_integral(vector) -> bool:
    rank = len(vector)
    for index in range(rank):
        basis = [0] * rank
        basis[index] = 1
        scalar = Fraction(pairing(basis, vector), 2)
        image = [
            Fraction(value) + scalar * root_value
            for value, root_value in zip(basis, vector)
        ]
        if any(coordinate.denominator != 1 for coordinate in image):
            return False
    return True


def test_unimodular_norm_minus_four_integrality_is_exactly_doubling():
    checked = 0
    for vector in itertools.product(range(-2, 3), repeat=5):
        if vector == (0,) * 5 or pairing(vector, vector) != -4:
            continue
        checked += 1
        integral = reflection_is_integral(vector)
        even = all(coordinate % 2 == 0 for coordinate in vector)
        assert integral is even
        if even:
            half = tuple(coordinate // 2 for coordinate in vector)
            assert pairing(half, half) == -1
    assert checked > 100

    # Even nonnegative tail multiplicities with total at most one vanish.
    for multiplicities in itertools.product(range(0, 5, 2), repeat=4):
        if sum(multiplicities) <= 1:
            assert multiplicities == (0, 0, 0, 0)


def independent_augmented_rays(m: int):
    rays = closed_form_duals(m)
    v_10 = rays[10]
    coefficients = [pairing(ray, v_10) for ray in rays]
    assert coefficients == [
        3, 2, 4, 6, 5, 4, 3, 2, 1, 0, *((-1,) * (m - 9))
    ]
    surviving = [
        rays[index]
        for index, coefficient in enumerate(coefficients)
        if coefficient >= 0
    ]
    mixed = []
    for left, positive in enumerate(coefficients):
        if positive <= 0:
            continue
        for right, negative in enumerate(coefficients):
            if negative >= 0:
                continue
            # Here negative=-1.  These are the unique two-coordinate
            # solutions of positive*t_left+negative*t_right=0.
            mixed.append(tuple(
                rays[left][coordinate] + positive * rays[right][coordinate]
                for coordinate in range(m + 1)
            ))
    return coefficients, tuple((*surviving, *mixed))


def test_one_wall_extreme_rays_and_first_mixed_escape_independently():
    finite = []
    for m in range(10, 33):
        coefficients, rays = independent_augmented_rays(m)
        assert len(rays) == 10 + 9 * (m - 9)
        assert len(set(rays)) == len(rays)
        squares = [pairing(ray, ray) for ray in rays]
        if all(square >= 0 for square in squares):
            finite.append(m)

        # Independently verify the closed square formula for every mixed ray.
        duals = closed_form_duals(m)
        offset = 10
        position = offset
        for i, positive in enumerate(coefficients):
            if positive <= 0:
                continue
            for j in range(10, m + 1):
                assert rays[position] == tuple(
                    a + positive * b for a, b in zip(duals[i], duals[j])
                )
                assert pairing(rays[position], rays[position]) == (
                    pairing(duals[i], duals[i])
                    + positive ** 2 * (11 - j)
                )
                position += 1
        assert position == len(rays)

    assert finite == [10, 11]
    v_1 = closed_form_duals(12)[1]
    v_12 = closed_form_duals(12)[12]
    first_mixed = tuple(a + 2 * b for a, b in zip(v_1, v_12))
    assert pairing(first_mixed, first_mixed) == -4


def test_F2_duals_inverse_diagonal_and_visible_threshold_independently():
    finite = []
    for m in range(2, 41):
        roots = f2_roots(m)
        rays = f2_closed_form_duals(m)
        assert [
            [f2_pairing(ray, root) for root in roots]
            for ray in rays
        ] == [
            [int(row == column) for column in range(m + 2)]
            for row in range(m + 2)
        ]
        squares = [f2_pairing(ray, ray) for ray in rays]
        assert squares == [
            0,
            2,
            *[
                2 * min(index, 2) ** 2 - index
                for index in range(1, m + 1)
            ],
        ]
        if all(square >= 0 for square in squares):
            finite.append(m)
    assert finite == list(range(2, 9))


def test_F2_escape_adjunction_and_doubled_stage_nine_wall():
    for m in range(9, 33):
        rays = f2_closed_form_duals(m)
        k_class = f2_canonical(m)
        for j in range(9, m + 1):
            ray = rays[j + 1]
            assert f2_pairing(ray, ray) == 8 - j
            assert f2_pairing(k_class, ray) == j - 8
            assert 1 + (
                f2_pairing(ray, ray) + f2_pairing(k_class, ray)
            ) // 2 == 1

        u_9 = rays[10]
        if m == 9:
            assert u_9 == tuple(-value for value in k_class)
        wall = tuple(2 * value for value in u_9)
        assert f2_pairing(wall, wall) == -4
        assert f2_pairing(k_class, wall) == 2
        assert 1 + (
            f2_pairing(wall, wall) + f2_pairing(k_class, wall)
        ) // 2 == 0
        assert [
            f2_pairing(wall, rays[j + 1])
            for j in range(9, m + 1)
        ] == [-2] * (m - 8)


def independent_F2_augmented_rays(m: int):
    rays = f2_closed_form_duals(m)
    u_9 = rays[10]
    coefficients = [f2_pairing(ray, u_9) for ray in rays]
    assert coefficients == [
        2, 4, 3, 6, 5, 4, 3, 2, 1, 0,
        *((-1,) * (m - 8)),
    ]
    surviving = [
        ray
        for ray, coefficient in zip(rays, coefficients)
        if coefficient >= 0
    ]
    mixed = []
    for left, positive in enumerate(coefficients):
        if positive <= 0:
            continue
        for right, negative in enumerate(coefficients):
            if negative >= 0:
                continue
            mixed.append(tuple(
                a + positive * b
                for a, b in zip(rays[left], rays[right])
            ))
    return coefficients, tuple((*surviving, *mixed))


def test_F2_one_wall_extremes_and_first_mixed_escape_independently():
    finite = []
    for m in range(9, 33):
        coefficients, rays = independent_F2_augmented_rays(m)
        assert len(rays) == 10 + 9 * (m - 8)
        assert len(set(rays)) == len(rays)
        squares = [f2_pairing(ray, ray) for ray in rays]
        if all(square >= 0 for square in squares):
            finite.append(m)

        duals = f2_closed_form_duals(m)
        position = 10
        for i, positive in enumerate(coefficients):
            if positive <= 0:
                continue
            for j in range(9, m + 1):
                assert rays[position] == tuple(
                    a + positive * b
                    for a, b in zip(duals[i], duals[j + 1])
                )
                assert f2_pairing(rays[position], rays[position]) == (
                    f2_pairing(duals[i], duals[i])
                    + positive ** 2 * (10 - j)
                )
                position += 1
        assert position == len(rays)

    assert finite == [9, 10]
    fibre = f2_closed_form_duals(11)[0]
    u_11 = f2_closed_form_duals(11)[12]
    first_mixed = tuple(a + 2 * b for a, b in zip(fibre, u_11))
    assert f2_pairing(first_mixed, first_mixed) == -4


def test_archive_replays_and_hostile_report_keeps_the_scope_narrow():
    subprocess.run(
        [
            sys.executable,
            str(
                ROOT
                / "scripts"
                / "certify_rational_single_chain_standard_simplex_gate.py"
            ),
            "--max-m",
            "32",
            "--check",
            str(ARCHIVE),
        ],
        cwd=ROOT,
        check=True,
    )
    text = AUDIT.read_text(encoding="utf-8")
    required = [
        "**PASS.**",
        r"\(m\leq9\)",
        "arithmetic-genus",
        r"does not assert that \(v_j\) is effective",
        "forces the tail sum to vanish",
        "direct full-NS lane",
        "support must contain exactly one",
        r"\(v_1+2v_{12}\)",
        "does not prove a chain-depth bound",
        "not an effectivity theorem",
        "fibre-following",
        r"visible threshold is \(m=8\)",
        r"\(f+2u_{11}\)",
        "no missing higher-support extreme type",
        "complete-nef/maximal counterfamily",
    ]
    for phrase in required:
        assert phrase in text

    forbidden_controls = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden_controls for byte in AUDIT.read_bytes())

    producer = (
        ROOT / "results/rational_single_chain_standard_simplex_gate.md"
    ).read_text(encoding="utf-8")
    assert "The second centre\nis the fibre-following free" in producer
    assert "so it is not satellite" in producer
    assert "a free chain on a fibre" not in producer
    f2 = json.loads(ARCHIVE.read_text(encoding="utf-8"))[
        "F2_minimal_base_companion"
    ]
    assert "fibre-following free point" in f2["marking"]
    assert f2["visible_simplex"]["finite_volume_exactly_m_at_most"] == 8
    assert f2["one_wall_augmented_cone"]["finite_volume_exactly_m"] == [9, 10]
    assert f2["geometric_effectivity_of_doubled_numerical_wall_claimed"] is False
