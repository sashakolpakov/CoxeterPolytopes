"""Independent hostile checks for exceptional-forest comparability."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "rational_exceptional_forest_comparability_hostile_audit.md"
)
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "rational_exceptional_forest_comparability.json"
)
SOURCE_HASHES = {
    "results/rational_exceptional_forest_comparability_theorem.md":
        "3b3760c2636827b7f0f8cf6db506c46fc3da87bd7724f25c789dd26554541450",
    "scripts/certify_rational_exceptional_forest_comparability.py":
        "161198bdbcb82bf63390da8c28f789c5eac5dba70ecc305a5005812822803ba0",
    "data/enumerations/rational_exceptional_forest_comparability.json":
        "10c89c66639e112ce57b2ab36f4f33feeb413a13fc21f04d475e813d5fa4d4eb",
    "tests/test_rational_exceptional_forest_comparability.py":
        "c3d3b565b6cac782e4eba01902cfaedcfad75779d78aee7e53097f15462f99f6",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_equal_weights(proximity, left, right):
    children = [[] for _ in proximity]
    for child, predecessors in enumerate(proximity):
        for predecessor in predecessors:
            children[predecessor].append(child)

    def evaluate(slacks):
        weights = [0] * len(proximity)
        for vertex in range(len(proximity) - 1, -1, -1):
            weights[vertex] = slacks[vertex] + sum(
                weights[child] for child in children[vertex]
            )
        return weights

    slacks = [1] * len(proximity)
    slacks[left] = slacks[right] = 0
    base = evaluate(slacks)
    target = max(base[left], base[right]) + 7
    slacks[left] = target - base[left]
    slacks[right] = target - base[right]
    weights = evaluate(slacks)
    return weights, [
        weights[vertex] - sum(weights[child] for child in children[vertex])
        for vertex in range(len(proximity))
    ]


def test_source_pins_and_archive_scope():
    assert {
        relative: sha256(ROOT / relative)
        for relative in SOURCE_HASHES
    } == SOURCE_HASHES
    archive = json.loads(ARCHIVE.read_text())
    assert archive["payload"]["verdict"] == "PROVED"
    scope = archive["payload"]["scope"]
    assert scope["odd_parity_proves_maximality"] is False
    assert scope["depth_bound_obtained"] is False


def test_independent_satellite_weight_construction():
    proximity = (
        (),
        (0,),
        (0,),
        (1, 0),
        (2, 0),
        (3, 1),
        (4, 2),
        (5, 3),
        (6, 4),
    )
    for left, right in ((1, 2), (3, 4), (5, 6), (7, 8)):
        weights, margins = independent_equal_weights(
            proximity, left, right
        )
        assert weights[left] == weights[right]
        assert min(weights) > 0
        assert min(margins) > 0

    # Comparable ancestors cannot be equal under strict proximity margins.
    weights, _ = independent_equal_weights(proximity, 1, 2)
    assert weights[0] > weights[1] > weights[3] > weights[5] > weights[7]


def test_difference_reflection_and_ample_mirror_are_exact():
    rank = 8
    left, right = 2, 6
    root = [0] * rank
    root[left], root[right] = 1, -1

    def pair(first, second):
        return -sum(a * b for a, b in zip(first, second))

    assert pair(root, root) == -2

    def reflect(vector):
        return [
            value + pair(vector, root) * root_value
            for value, root_value in zip(vector, root)
        ]

    basis = [
        [int(row == column) for row in range(rank)]
        for column in range(rank)
    ]
    assert reflect(basis[left]) == basis[right]
    assert reflect(basis[right]) == basis[left]
    weights = [19, 11, 7, 5, 3, 2, 7, 1]
    ample_exceptional_part = [-value for value in weights]
    assert pair(ample_exceptional_part, root) == 0


def test_even_parity_formula_and_odd_integrality_failure():
    # Exceptional coefficients suffice because the core is root-orthogonal.
    root = (1, -1, 0, 0)

    def surface_pair(first, second):
        return -sum(a * b for a, b in zip(first, second))

    even_boundary = (3, 1, 0, 0)
    odd_boundary = (2, 1, 0, 0)
    assert surface_pair(even_boundary, root) == -2
    assert surface_pair(odd_boundary, root) == -1

    even_shift = [
        Fraction(surface_pair(even_boundary, root), 2) * value
        for value in root
    ]
    odd_shift = [
        Fraction(surface_pair(odd_boundary, root), 2) * value
        for value in root
    ]
    assert all(value.denominator == 1 for value in even_shift)
    assert any(value.denominator == 2 for value in odd_shift)
    assert Fraction(surface_pair(odd_boundary, root), 2).denominator == 2


def test_hostile_audit_states_sharp_gain_and_limit():
    text = AUDIT.read_text()
    required = [
        "**PASS.**",
        "A satellite point lies",
        "forces strict decrease",
        r"\(\pi\)-ample",
        "parity is an if-and-only-if test",
        "Odd parity does not prove maximality",
        r"at most \(2^b\) vertices",
        "arbitrary depth, as the affine",
    ]
    for phrase in required:
        assert phrase in text
