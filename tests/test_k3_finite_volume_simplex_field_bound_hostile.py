"""Independent hostile replay of the finite-volume K3 simplex theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_finite_volume_simplex_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_finite_volume_simplex_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_finite_volume_simplex_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_finite_volume_simplex_field_bound.py"
AUDIT = ROOT / "results/k3_finite_volume_simplex_field_bound_hostile_audit.md"
PINS = {
    REPORT: "80e14c2f63a7ce211bc25697ccd8fd694d498f9647b20979405cba1b66716051",
    CERTIFIER: "665b679b7bd5bedbc76cc2fdb2afadea9b0fc506967f9b44e8dd77dfeeb26009",
    ARCHIVE: "ce4f721c1720f0e9cb3d11c327551048418a8708bf51095e82f63d1716d0e5d0",
    OWNER_TESTS: "501f16d39b709ec394752c80fb02b4de1eb72a28f998eb06de5065b39f515dd8",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multiply(left, right):
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def trace(matrix):
    return sum(matrix[i][i] for i in range(len(matrix)))


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_projective_module_cycle_trace_formula_replays_in_a_free_frame():
    cartan = (
        (Fraction(2), Fraction(-1, 2), Fraction(-2, 3), Fraction(-3, 4)),
        (Fraction(-4, 5), Fraction(2), Fraction(-5, 6), Fraction(-6, 7)),
        (Fraction(-7, 8), Fraction(-8, 9), Fraction(2), Fraction(-9, 10)),
        (Fraction(-10, 11), Fraction(-11, 12), Fraction(-12, 13), Fraction(2)),
    )
    projections = tuple(
        tuple(
            tuple(cartan[i][column] if row == i else Fraction(0) for column in range(4))
            for row in range(4)
        )
        for i in range(4)
    )
    for cycle in ((0, 1), (0, 1, 2), (0, 1, 2, 3)):
        product = projections[cycle[0]]
        for index in cycle[1:]:
            product = multiply(product, projections[index])
        expected = Fraction(1)
        for position, i in enumerate(cycle):
            expected *= cartan[i][cycle[(position + 1) % len(cycle)]]
        assert trace(product) == expected


def test_symmetric_cycle_square_bound_for_elliptic_and_parabolic_pairs():
    for pair_values in (
        (Fraction(0), Fraction(1), Fraction(3)),
        (Fraction(4), Fraction(4), Fraction(4)),
        (Fraction(1, 4), Fraction(15, 4), Fraction(2)),
    ):
        assert all(0 <= value <= 4 for value in pair_values)
        cycle_square = pair_values[0] * pair_values[1] * pair_values[2]
        assert cycle_square <= 4**3
    # An ultraparallel pair would have p>4 and is precisely the excluded case.
    assert Fraction(17, 4) > 4


def test_closed_walk_factorization_needs_only_simple_cycles():
    weights = {
        (0, 1): Fraction(2),
        (1, 2): Fraction(3),
        (2, 1): Fraction(5),
        (1, 3): Fraction(7),
        (3, 0): Fraction(11),
    }
    walk = (0, 1, 2, 1, 3, 0)
    full = Fraction(1)
    for left, right in zip(walk, walk[1:]):
        full *= weights[(left, right)]
    first_cycle = weights[(1, 2)] * weights[(2, 1)]
    second_cycle = weights[(0, 1)] * weights[(1, 3)] * weights[(3, 0)]
    assert full == first_cycle * second_cycle


def test_primitive_grid_and_uniform_exponent_table_are_independent():
    assert len(tuple(combinations(range(6), 2))) == 15
    rows = []
    for module_rank in range(3, 21):
        degree = 20 // module_rank
        exponent = 0 if degree < 2 else (module_rank + 5) * degree * (degree - 1)
        rows.append((module_rank, degree, exponent))
    assert rows[:8] == [
        (3, 6, 240),
        (4, 5, 180),
        (5, 4, 120),
        (6, 3, 66),
        (7, 2, 24),
        (8, 2, 26),
        (9, 2, 28),
        (10, 2, 30),
    ]
    assert max(exponent for _, _, exponent in rows) == 240
    assert 2**240 == 256**30


def test_archive_replays_and_owner_scope_is_unchanged():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["uniform_binary_discriminant_exponent"] == 240
    assert payload["uniform_field_discriminant_bound_expression"] == "256^30"
    conclusions = payload["conclusions"]
    assert conclusions["finite_field_set"] is True
    assert conclusions["ideal_finite_volume_simplices_covered"] is True
    assert conclusions["imposed_mixed_height_needed"] is False
    assert conclusions["angle_classification_used"] is False
    assert conclusions["reflection_maximality_used_for_bound"] is False
    assert conclusions["carrier_existence_proved"] is False
    assert conclusions["general_higher_facet_polytope_closed"] is False


def test_report_and_audit_pin_finite_volume_and_hyperideal_boundary():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "Every two facets of an ordinary finite-volume simplex meet at a finite or ideal ridge",
        "Every proper rank-two subdiagram of a finite-volume Coxeter simplex is elliptic or parabolic",
        r"&=a_{i_1i_2}a_{i_2i_3}\cdots a_{i_\ell i_1}",
        "Every closed directed walk decomposes into simple directed cycles",
        r"2^{(m+5)d(d-1)}",
        "finite-volume simplices with ideal faces",
        "truncated hyperideal polytopes",
    ):
        assert phrase in report
    for phrase in (
        "PASS for ordinary finite-volume simplices",
        "normalized two-normal Gram block is positive semidefinite",
        "canonical trace of an endomorphism of a finite projective module",
        r"maximum is \(240\)",
        "gratuitous correspondence-field extension is not covered",
    ):
        assert phrase in audit


def test_tex_transport_and_no_classification_import_are_clean():
    combined = REPORT.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    normalized = " ".join(combined.split())
    assert "No Coxeter-simplex classification, angle list" in normalized
    assert "No list of simplex diagrams or dihedral angles is being used" in normalized
