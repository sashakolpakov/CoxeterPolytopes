"""Standalone hostile scope replay for the rank-four K3 field no-go."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rank4_k3_correspondence_infinite_field_no_go.md"
CERTIFIER = (
    ROOT
    / "scripts/certify_rank4_k3_correspondence_infinite_field_no_go.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rank4_k3_correspondence_infinite_field_no_go.json"
)
OWNER_TESTS = (
    ROOT / "tests/test_rank4_k3_correspondence_infinite_field_no_go.py"
)
AUDIT = (
    ROOT
    / "results/rank4_k3_correspondence_infinite_field_no_go_hostile_scope_audit.md"
)
PINS = {
    REPORT: "3406e7d4f43ac1c96c13bd93dd8f6eae7e90564efcf06667262de4c155dd26c8",
    CERTIFIER: "cb70ef2f39da8309cac648a87fd0731187812208a955dd3ba587a1a8c83e82f9",
    ARCHIVE: "d9792dd05caa49335b593ae621296f1c48721328e02b56a6823d715a04c5cd88",
    OWNER_TESTS: "3d67e694fd4a8192b9da7c90b2499e115141d8b596f0ce30a143d1278f8268ff",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def transpose(matrix):
    return tuple(tuple(row[column] for row in matrix) for column in range(len(matrix[0])))


def multiply(left, right):
    columns = transpose(right)
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def block_diagonal(*blocks):
    size = sum(len(block) for block in blocks)
    result = [[0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        for row in range(len(block)):
            for column in range(len(block)):
                result[offset + row][offset + column] = block[row][column]
        offset += len(block)
    return tuple(tuple(row) for row in result)


def determinant(matrix):
    work = [[Fraction(entry) for entry in row] for row in matrix]
    value = Fraction(1)
    for column in range(len(work)):
        pivot = next(row for row in range(column, len(work)) if work[row][column])
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            value = -value
        diagonal = work[column][column]
        value *= diagonal
        for row in range(column + 1, len(work)):
            ratio = work[row][column] / diagonal
            for index in range(column, len(work)):
                work[row][index] -= ratio * work[column][index]
    assert value.denominator == 1
    return int(value)


def owner_matrices(a: int):
    c = ((0, -1), (1, 2 * a))
    b0 = ((-1, 0), (0, 1))
    ba = ((1 - 2 * a, -1), (-1, -1))
    g = block_diagonal(
        tuple(tuple(2 * entry for entry in row) for row in b0),
        tuple(tuple(2 * entry for entry in row) for row in ba),
    )
    operator = block_diagonal(c, c)
    return c, b0, ba, g, operator


def test_owner_bytes_are_frozen():
    assert {path: sha256(path) for path in PINS} == PINS


def test_self_adjointness_determinant_and_polynomial_rebuild_independently():
    identity = tuple(
        tuple(1 if row == column else 0 for column in range(4))
        for row in range(4)
    )
    for a in range(2, 513):
        c, b0, ba, g, operator = owner_matrices(a)
        assert multiply(transpose(c), b0) == multiply(b0, c)
        assert multiply(transpose(c), ba) == multiply(ba, c)
        assert multiply(transpose(operator), g) == multiply(g, operator)
        square = multiply(operator, operator)
        polynomial = tuple(
            tuple(
                square[row][column]
                - 2 * a * operator[row][column]
                + identity[row][column]
                for column in range(4)
            )
            for row in range(4)
        )
        assert polynomial == ((0, 0, 0, 0),) * 4
        assert ba[0][0] < 0
        assert determinant(ba) == 2 * (a - 1) > 0
        assert determinant(g) == -32 * (a - 1)


def test_K3_embedding_gram_and_primitivity_rebuild_independently():
    u = ((0, 1), (1, 0))
    target = block_diagonal(u, u, u, ((-2,),))
    for a in range(2, 513):
        vectors = (
            (1, -1, 0, 0, 0, 0, 0),
            (0, 0, 1, 1, 0, 0, 0),
            (0, 0, 0, 0, 1, 2 - 2 * a, 1),
            (0, 0, 0, 0, 0, 0, 1),
        )
        columns = transpose(vectors)
        pullback = multiply(multiply(transpose(columns), target), columns)
        assert pullback == owner_matrices(a)[3]
        pivot = tuple(
            tuple(columns[row][column] for column in range(4))
            for row in (0, 2, 4, 6)
        )
        assert determinant(pivot) == 1


def test_archive_scope_forbids_maximal_same_P_promotion():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = envelope["payload"]
    scope = payload["scope"]
    assert payload["family"]["NS_rank"] == 4
    assert payload["family"]["primitive_full_trace_block"] is True
    assert scope["exact_complete_restricted_root_component_returned"] is True
    false_keys = (
        "prescribed_input_chamber_returned",
        "restricted_component_Coxeter_claimed",
        "restricted_component_finite_sided_claimed",
        "restricted_component_finite_volume_claimed",
        "reflection_saturation_claimed",
        "reflection_maximality_claimed",
        "exact_maximal_same_P_infinite_field_family_claimed",
    )
    assert all(scope[key] is False for key in false_keys)
    assert scope["strong_exact_maximal_field_set_gate"] == "open"


def test_Pell_samples_do_not_conflate_parameters_with_fields():
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    rows = payload["field_samples"]
    for row in rows:
        assert (
            row["a"] * row["a"]
            - row["squarefree_field_parameter"] * row["Pell_y"] ** 2
            == 1
        )
    mapping = {row["a"]: row["squarefree_field_parameter"] for row in rows}
    assert mapping[2] == mapping[7] == 3
    assert payload["distinct_sample_field_count"] > 50
    assert (
        payload["field_infinitude_proof"]["bounded_degree_alone_implies_finite_fields"]
        is False
    )


def test_report_and_scope_audit_state_the_nonclaim_literally():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "does **not** prove that the restricted eigen-nef component is finite-volume",
        "reflection-maximal same-\\(P\\) eigen-nef return could still impose",
        "strongest exact-maximal field-set gate remains open",
    ):
        assert phrase in report
    for phrase in (
        "PASS at the stated weak boundary",
        "not a second provenance claim",
        "The following stronger statement is **not certified**",
        "prescribed same-P eigen-nef packages occur for infinitely many fields",
    ):
        assert phrase in audit


def test_certifier_archive_replay_passes_as_a_separate_process():
    completed = subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "PASS archive replay" in completed.stdout
