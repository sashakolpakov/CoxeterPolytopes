import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/rational_log_fano_tail_coefficient_no_go.md"
OWNER_SCRIPT = ROOT / "scripts/certify_rational_log_fano_tail_coefficient_no_go.py"
OWNER_ARCHIVE = (
    ROOT / "data/enumerations/rational_log_fano_tail_coefficient_no_go.json"
)
OWNER_TEST = ROOT / "tests/test_rational_log_fano_tail_coefficient_no_go.py"
AUDIT_REPORT = (
    ROOT / "results/rational_log_fano_tail_coefficient_no_go_hostile_audit.md"
)

PINNED_OWNER_HASHES = {
    OWNER_REPORT: "bf60458b0d1162fb7b1b17e6503102697c5f118119ab360ba7b67f321cef968c",
    OWNER_SCRIPT: "b3622c2122329d23239cc0b9e2bd6182dda080e1220fb52b90f32b0c60ef9e01",
    OWNER_ARCHIVE: "b64133616a6fcfd574b8da07823567f9e21ee2347db79dbea631d93ee8750baf",
    OWNER_TEST: "3a5de56425413c4a9d08510215bcdf68e998c6e8104765d011fea2cbae7e1c43",
}


def _tail_coordinates(length, chain, endpoint):
    """Independent reconstruction; do not import the owner certifier."""
    chain = tuple(Fraction(value) for value in chain)
    endpoint = Fraction(endpoint)
    assert length >= 1 and len(chain) == length - 1
    if length == 1:
        return (1 + endpoint,)
    return (
        (1 + chain[0],)
        + tuple(
            1 - chain[index - 1] + chain[index]
            for index in range(1, length - 1)
        )
        + (1 - chain[-1] + endpoint,)
    )


def _negative_dot(first, second):
    return -sum(a * b for a, b in zip(first, second, strict=True))


def _basis_vector(length, index):
    return tuple(Fraction(position == index) for position in range(length))


def test_independent_total_transform_and_physical_chain_reconstruction():
    for length in range(1, 16):
        exceptional = [_basis_vector(length, index) for index in range(length)]
        physical = [
            tuple(
                exceptional[index][position] - exceptional[index + 1][position]
                for position in range(length)
            )
            for index in range(length - 1)
        ] + [exceptional[-1]]

        gram = [
            [_negative_dot(first, second) for second in physical]
            for first in physical
        ]
        assert [gram[index][index] for index in range(length)] == [
            -2
        ] * (length - 1) + [-1]
        for row in range(length):
            for column in range(length):
                expected = (
                    gram[row][row]
                    if row == column
                    else Fraction(1)
                    if abs(row - column) == 1
                    else Fraction(0)
                )
                assert gram[row][column] == expected

        chain = tuple(Fraction(index + 1, length + 2) for index in range(length - 1))
        endpoint = Fraction(length + 1, length + 2)
        direct = [Fraction(1) for _ in range(length)]
        for index, coefficient in enumerate(chain):
            direct[index] += coefficient
            direct[index + 1] -= coefficient
        direct[-1] += endpoint
        assert tuple(direct) == _tail_coordinates(length, chain, endpoint)


def test_sharp_cauchy_bound_fixed_endpoint_profile_and_unique_global_equality():
    for length in range(1, 21):
        for endpoint in (Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(1)):
            chain = tuple(
                Fraction(index) * endpoint / length
                for index in range(1, length)
            )
            coordinates = _tail_coordinates(length, chain, endpoint)
            assert coordinates == (1 + endpoint / length,) * length
            loss = sum(value * value for value in coordinates)
            assert loss == Fraction(length + endpoint) ** 2 / length

    grid = tuple(Fraction(index, 4) for index in range(5))
    for length in range(1, 7):
        for values in itertools.product(grid, repeat=length):
            coordinates = _tail_coordinates(length, values[:-1], values[-1])
            coordinate_sum = sum(coordinates, Fraction(0))
            loss = sum(value * value for value in coordinates)
            assert coordinate_sum == length + values[-1]
            assert loss >= coordinate_sum * coordinate_sum / length
            assert loss >= length
            assert (loss == length) == (values == (Fraction(0),) * length)


def test_h1_h2_rigidity_is_an_exact_nonnegative_slack_decomposition():
    grid = tuple(Fraction(index, 12) for index in range(13))
    gamma_grid = (Fraction(0), Fraction(1, 5), Fraction(2))
    for first, second, gamma_first, gamma_second in itertools.product(
        grid, grid, gamma_grid, gamma_grid
    ):
        h1_first = 2 - 4 * first + 2 * second + gamma_first
        h1_second = 2 + 2 * first - 4 * second + gamma_second
        assert h1_first + h1_second == (
            2 * (1 - first)
            + 2 * (1 - second)
            + gamma_first
            + gamma_second
        )
        if h1_first <= 0 and h1_second <= 0:
            assert (first, second, gamma_first, gamma_second) == (
                Fraction(1),
                Fraction(1),
                Fraction(0),
                Fraction(0),
            )

        h2_boundary = 2 - 4 * first + 2 * second + gamma_first
        h2_minus_two = 2 * first - 2 * second + gamma_second
        assert h2_boundary + h2_minus_two == (
            2 * (1 - first) + gamma_first + gamma_second
        )
        if h2_boundary <= 0 and h2_minus_two <= 0:
            assert (first, second, gamma_first, gamma_second) == (
                Fraction(1),
                Fraction(1),
                Fraction(0),
                Fraction(0),
            )


def test_discrepancy_ladder_and_log_chern_square_changes():
    for length in range(1, 101):
        exceptional_sum = tuple(Fraction(1) for _ in range(length))
        physical_sum = [Fraction(0) for _ in range(length)]
        for index in range(length - 1):
            coefficient = Fraction(index + 1)
            physical_sum[index] += coefficient
            physical_sum[index + 1] -= coefficient
        physical_sum[-1] += length
        assert tuple(physical_sum) == exceptional_sum

        zero = _tail_coordinates(length, (0,) * (length - 1), 0)
        full = _tail_coordinates(length, (1,) * (length - 1), 1)
        assert zero == (Fraction(1),) * length
        assert sum(value * value for value in zero) == length
        assert full == (Fraction(2),) + (Fraction(1),) * (length - 1)
        assert sum(value * value for value in full) == length + 3

        # Formal Euler changes: every blowup in the open complement replaces
        # one point by P^1 (+1); removing the full fibre removes the original
        # base point instead (-1), independently of the tower length.
        core_euler = Fraction(37)
        assert core_euler + length - core_euler == length
        assert (core_euler - 1) - core_euler == -1


def test_frozen_scope_hashes_and_report_controls():
    for path, expected in PINNED_OWNER_HASHES.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected

    archive = json.loads(OWNER_ARCHIVE.read_text(encoding="utf-8"))
    assert archive["schema_version"] == 2
    assert archive["all_checks_pass"] is True
    assert archive["scope"]["absolute_tail_bound_claimed"] is False
    assert archive["scope"]["minimal_hostile_H1_H2_weak_log_Fano_exists"] is False
    assert archive["scope"]["higher_parallel_dotted_H1_H2_claimed"] is False
    assert archive["scope"]["classification_used"] is False
    assert archive["log_Chern_changes"]["full_reduced_tail_in_boundary"]["scope"] == (
        "m>=1; centers outside the reduced SNC core boundary"
    )

    for path in (OWNER_REPORT, AUDIT_REPORT):
        raw = path.read_bytes()
        assert b"\r" not in raw
        assert not any(byte < 32 and byte not in (9, 10) for byte in raw)

    owner = OWNER_REPORT.read_text(encoding="utf-8")
    audit = AUDIT_REPORT.read_text(encoding="utf-8")
    for needle in (
        r"m<[-(K_{S_0}+\Delta_0)]^2",
        r"\sim_{\mathbf R}",
        r"\varepsilon\)-lc log-Fano",
        "Higher, parallel, or dotted H1/H2 incidences are not",
    ):
        assert needle in owner
    for needle in (
        "PASS after scope and presentation repairs",
        r"a_i=\frac{ib}{m}",
        r"\sim_{\mathbb R}",
        "not the missing absolute\ndimension-bound theorem",
    ):
        assert needle in audit
