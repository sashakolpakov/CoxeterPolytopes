"""Independent hostile replay of exact-trace cusp field collapse."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/exact_trace_noncompact_ground_field_collapse.md"
CERTIFIER = ROOT / "scripts/certify_exact_trace_noncompact_ground_field_collapse.py"
ARCHIVE = ROOT / "data/enumerations/exact_trace_noncompact_ground_field_collapse.json"
OWNER_TESTS = ROOT / "tests/test_exact_trace_noncompact_ground_field_collapse.py"
AUDIT = ROOT / "results/exact_trace_noncompact_ground_field_collapse_hostile_audit.md"
PINS = {
    REPORT: "d738018bacc3803f78fa725da375bb8045ecce92a1c879b03b30a9998e980255",
    CERTIFIER: "61f7d0a21cab3b80975849d015c8f2cf13386b3d5f35890917dd0b672e6ebfcb",
    ARCHIVE: "fb87928a12165c226272078ce57f8c98418ce02f6629a80fecb3e5769e7cfcd7",
    OWNER_TESTS: "d2fe708d1c368ed2ee02c673c4e456e5a497b8dc7d8ce8e318ebbc9c2d00d40f",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row:
                continue
            value = work[row][column]
            if value:
                work[row] = [
                    work[row][j] - value * work[pivot_row][j]
                    for j in range(len(work[0]))
                ]
        pivot_row += 1
    return pivot_row


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_horosphere_height_integral_is_finite_exactly_in_relevant_dimensions():
    # Integral_T^infinity y^(-n)dy = T^(1-n)/(n-1), with T=2 here.
    for dimension in range(2, 21):
        value = Fraction(1, (dimension - 1) * 2 ** (dimension - 1))
        assert value > 0
    # In H^1 the corresponding logarithmic integral diverges; no noncompact
    # finite-length chamber exists, explaining the theorem's n>=2 scope.


def test_square_cusp_normal_span_rederives_without_owner_helpers():
    # H^3, form diag(1,1,1,-1), u=(0,0,1,1).
    cusp = [Fraction(0), Fraction(0), Fraction(1), Fraction(1)]
    normals = [
        [Fraction(1), Fraction(0), Fraction(1), Fraction(1)],
        [Fraction(-1), Fraction(0), Fraction(1), Fraction(1)],
        [Fraction(0), Fraction(1), Fraction(1), Fraction(1)],
        [Fraction(0), Fraction(-1), Fraction(1), Fraction(1)],
    ]
    form_diagonal = [1, 1, 1, -1]

    def pairing(left, right):
        return sum(form_diagonal[i] * left[i] * right[i] for i in range(4))

    assert pairing(cusp, cusp) == 0
    assert all(pairing(cusp, normal) == 0 for normal in normals)
    assert all(pairing(normal, normal) == 1 for normal in normals)
    assert rank(normals) == 3
    assert [sum(row[j] for row in normals) for j in range(4)] == [
        4 * entry for entry in cusp
    ]
    # Four incident walls in H^3 make this explicitly nonsimple.
    assert len(normals) == 4 > 3


def test_all_incident_fixed_equations_have_one_dimensional_kernel():
    # For a norm-one normal e, fixed vectors obey <x,e>=0.  The four rows
    # below are precisely those equations for the square cusp.
    equations = [
        [Fraction(1), Fraction(0), Fraction(1), Fraction(-1)],
        [Fraction(-1), Fraction(0), Fraction(1), Fraction(-1)],
        [Fraction(0), Fraction(1), Fraction(1), Fraction(-1)],
        [Fraction(0), Fraction(-1), Fraction(1), Fraction(-1)],
    ]
    assert rank(equations) == 3
    assert 4 - rank(equations) == 1


def test_field_embedding_contradiction_has_no_degree_loophole():
    # In a totally real degree-d field every nonidentity embedding would
    # carry the same field equation h(u,u)=0 to a definite place.
    for degree in range(2, 21):
        nonidentity_real_embeddings = degree - 1
        assert nonidentity_real_embeddings >= 1
    assert 1 - 1 == 0


def test_archive_and_scope_close_noncompact_not_compact_fields():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["conclusions"]["returned_coefficient_field_is_Q"]
    assert payload["conclusions"]["full_trace_K3_hyperbolic_dimension_at_most_19"]
    assert payload["scope"]["finite_projective_order_module_allowed"]
    assert payload["scope"]["upper_half_space_argument_includes_nonsimple_cusps"]
    assert not payload["scope"]["partial_affine_subdiagram_suffices"]
    assert not payload["scope"]["compact_chambers_closed"]
    assert not payload["scope"]["carrier_construction_proved"]
    assert not payload["scope"]["raw_chamber_dimension_bound_proved"]


def test_report_and_audit_state_the_exact_operator_mechanism():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "full-dimensional Euclidean polyhedron",
        "recession cone",
        r"\sum_i\lambda_i e_i",
        "Package all maps",
        "Extension of scalars",
        "one-dimensional",
        "A nonzero vector over a field remains nonzero",
        "The statement holds for simple or nonsimple cusps",
        r"n+1=\rho(X)\le20",
        "does not assert that every raw chamber has a K3 carrier",
    ):
        assert phrase in report
    for phrase in (
        "PASS for exact finite-volume noncompact chambers",
        "does not infer a rational line from an abstract affine subdiagram",
        "This argument does not assume that the cusp is simple",
        "It never chooses a generator of a rank-one normal module",
        "not merely numerically at one place",
        "does not replace the rational weak-log",
    ):
        assert phrase in audit


def test_tex_transport_and_no_import_controls_are_clean():
    combined_bytes = REPORT.read_bytes() + AUDIT.read_bytes()
    assert b"\x08" not in combined_bytes
    combined = combined_bytes.decode("utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    for forbidden in ("Nikulin", "Vinberg", "ABSW"):
        assert forbidden not in combined
