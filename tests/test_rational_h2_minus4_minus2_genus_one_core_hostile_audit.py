"""Independent hostile checks for the minimal H2 adjoint-core theorem."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_h2_minus4_minus2_genus_one_core.md"
AUDIT = ROOT / "results/rational_h2_minus4_minus2_genus_one_core_hostile_audit.md"
SCRIPT = ROOT / "scripts/certify_rational_h2_minus4_minus2_genus_one_core.py"
ARCHIVE = ROOT / "data/enumerations/rational_h2_minus4_minus2_genus_one_core.json"
OWNER_TEST = ROOT / "tests/test_rational_h2_minus4_minus2_genus_one_core.py"

OWNER_HASHES = {
    REPORT: "1464b327d588f9a5e8c5b2a6809667194d27c09a69aad5805250df4053c9634b",
    SCRIPT: "f75dc268fbaeea9efd23f39ebca8c7a29ac9903a0f2d20422c6399cad2f09adf",
    ARCHIVE: "d240a14f49df3672ce604aae0f7ddb67cf3700eca966b96b3a83e59917a5b7e5",
    OWNER_TEST: "f0153a604090ac2a65e2179183a075ca6a4093ad68b2f5b370a02ac0621fb80d",
}


def _det2(matrix: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def _solve2(
    matrix: tuple[tuple[int, int], tuple[int, int]],
    rhs: tuple[int, int],
) -> tuple[Fraction, Fraction]:
    determinant = _det2(matrix)
    return (
        Fraction(rhs[0] * matrix[1][1] - matrix[0][1] * rhs[1], determinant),
        Fraction(matrix[0][0] * rhs[1] - rhs[0] * matrix[1][0], determinant),
    )


def _pair(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return -2 * sum(a * b for a, b in zip(left, right, strict=True))


def _reflect(vector: tuple[int, ...], root: tuple[int, ...]) -> tuple[Fraction, ...]:
    coefficient = Fraction(2 * _pair(vector, root), _pair(root, root))
    return tuple(Fraction(x) - coefficient * r for x, r in zip(vector, root, strict=True))


def _tail_roots(length: int) -> tuple[tuple[int, ...], ...]:
    roots = []
    for index in range(length - 1):
        root = [0] * length
        root[index], root[index + 1] = 1, -1
        roots.append(tuple(root))
    terminal = [0] * length
    terminal[-1] = 1
    roots.append(tuple(terminal))
    return tuple(roots)


def test_pair_discrepancy_and_core_square_are_recomputed_independently():
    pair = ((-4, 2), (2, -2))
    assert _det2(pair) == 4
    assert _solve2(pair, (2, 0)) == (Fraction(-1), Fraction(-1))

    d_square = -4 - 2 + 2 * 2
    k_dot_d = (2 * 0 - 2 - (-4)) + (2 * 0 - 2 - (-2))
    assert (d_square, k_dot_d) == (-2, 2)
    assert 1 + Fraction(d_square + k_dot_d, 2) == 1
    assert 10 + 2 * k_dot_d + d_square == 12


def test_sharp_equality_blowup_arithmetic_is_exact():
    base_rank = 2
    base_squares = (2, 2)
    private_blowups = (6, 4)
    assert tuple(a - b for a, b in zip(base_squares, private_blowups, strict=True)) == (-4, -2)
    assert base_rank + sum(private_blowups) == 12


def test_free_tail_gram_and_reflections_without_owner_helpers():
    for length in (1, 2, 3, 7, 19):
        roots = _tail_roots(length)
        gram = tuple(tuple(_pair(a, b) for b in roots) for a in roots)
        expected = [[0] * length for _ in range(length)]
        for index in range(length):
            expected[index][index] = -4 if index + 1 < length else -2
            if index + 1 < length:
                expected[index][index + 1] = expected[index + 1][index] = 2
        assert gram == tuple(tuple(row) for row in expected)

        basis = tuple(
            tuple(int(i == j) for i in range(length)) for j in range(length)
        )
        for index, root in enumerate(roots[:-1]):
            images = tuple(_reflect(vector, root) for vector in basis)
            assert images[index] == tuple(Fraction(x) for x in basis[index + 1])
            assert images[index + 1] == tuple(Fraction(x) for x in basis[index])
        terminal_images = tuple(_reflect(vector, roots[-1]) for vector in basis)
        assert terminal_images[-1] == tuple(-Fraction(x) for x in basis[-1])


def test_archive_scope_and_conditional_bounds_are_frozen():
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert payload["pair_arithmetic"]["adjoint_square_formula"] == "12-rho"
    assert payload["scope"] == {
        "core_only": True,
        "original_rank_bound_claimed": False,
        "pair": "B^2=-4, R^2=-2, B.R=2 on a smooth rational surface",
        "reflection_maximal_counterfamily_claimed": False,
    }
    assert payload["conditional_total_marked_degree_bound"] == {
        "premise": (
            "every contracted tail component is a physical facet and "
            "T is the total integral ample-marking degree of all physical facets"
        ),
        "tail_length_bound": "m<=T",
        "picard_rank_bound": "rho(S)<=12+T",
        "hyperbolic_dimension_bound": "dim<=11+T",
        "per_facet_degree_bound_alone_is_sufficient": False,
    }
    assert payload["conditional_trace_bound"] == {
        "rank": "d(n+1)",
        "bound_if_the_trace_block_descends": "d(n+1)<=12",
        "rational_field_bound": "n<=11",
        "descent_is_automatic": False,
    }


def test_owner_hashes_and_report_transport_hygiene():
    for path, expected in OWNER_HASHES.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected

    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = raw.decode("utf-8")
    for lost_delimiter_signature in (
        "Let (S) be",
        "exceptional (B_m) tail",
        "on (mathbf P^1)",
        "orbifold-(c_2)",
        "by (operatorname {NS}(S))",
    ):
        assert lost_delimiter_signature not in text
    for required in (
        r"\omega_D\simeq\mathcal O_D",
        r"\boxed{\rho(S_0)\leq12}",
        r"m\leq T",
        "bound on each individual facet degree",
        r"\boxed{d(n+1)\leq12}",
        "Descent is not automatic",
    ):
        assert required in text

    audit = AUDIT.read_text(encoding="utf-8")
    assert "hostile verdict is **PASS**" in audit
    assert "not a PASS for an absolute original-rank" in audit
