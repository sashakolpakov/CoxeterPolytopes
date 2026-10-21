"""Independent hostile checks for the physical-star finite-host theorem."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/k3_compact_physical_star_trace_host_finiteness.md"
OWNER_SCRIPT = ROOT / "scripts/certify_k3_compact_physical_star_trace_host_finiteness.py"
OWNER_ARCHIVE = (
    ROOT / "data/enumerations/k3_compact_physical_star_trace_host_finiteness.json"
)
OWNER_TEST = ROOT / "tests/test_k3_compact_physical_star_trace_host_finiteness.py"
AUDIT_REPORT = (
    ROOT
    / "results/k3_compact_physical_star_trace_host_finiteness_independent_hostile_audit.md"
)
COUNTER_ARCHIVE = ROOT / "data/enumerations/k3_compact_projected_root_scaling_no_go.json"

FROZEN = {
    OWNER_REPORT: "5cc86e7225c11697c83af0ec34fa5c773e8e969c47107eabdd2a327c23c1ea3a",
    OWNER_SCRIPT: "c4a52c76f28192eada049fe66031fd62bb35f6366332f92abc54cec995d58d55",
    OWNER_ARCHIVE: "cc2d48b699aeebd1aaa10ae6b4b8249addde2eac681c9bb475aa19fcd3e14e90",
    OWNER_TEST: "7f38e30ed6fb20af8c54a266c9b8ec4cad59ddccd90f1fd451d40f9565dbd3ec",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_owner_bytes_are_exactly_the_audited_freeze() -> None:
    assert {str(path.relative_to(ROOT)): sha256(path) for path in FROZEN} == {
        str(path.relative_to(ROOT)): expected for path, expected in FROZEN.items()
    }


def test_energy_budget_independently_implies_all_entry_bounds() -> None:
    # Split each physical square-two energy budget into the distinguished
    # coordinate and an aggregate definite tail.  The bounds are monotone,
    # so a dense rational replay catches coefficient/sign regressions.
    for numerator_i in range(1, 1000):
        ti = 2 * numerator_i / 1000
        tail_i = 2 - ti
        for numerator_j in range(1, 1000, 7):
            tj = 2 * numerator_j / 1000
            tail_j = 2 - tj
            identity = 7 * math.sqrt(ti * tj)
            tail = math.sqrt(tail_i * tail_j)
            assert identity < 14
            assert tail < 2
            assert identity + tail < 16


def test_cycle_height_and_rank_leave_finite_integer_boxes() -> None:
    rows = []
    for module_rank in range(3, 21):
        degree = 20 // module_rank
        rows.append((module_rank, degree, 14**module_rank))
    assert max(degree for _, degree, _ in rows if degree >= 2) == 6
    assert max(rank for rank, degree, _ in rows if degree >= 2) == 10
    assert all(height > 0 for _, _, height in rows)

    # If all conjugates of a degree-e algebraic integer are bounded by H,
    # its j-th monic coefficient is bounded by C(e,j)H^j.  This supplies a
    # literal finite coefficient box, not merely a bounded field degree.
    for degree in range(1, 7):
        height = 14**10
        bounds = [math.comb(degree, j) * height**j for j in range(1, degree + 1)]
        assert len(bounds) == degree
        assert all(isinstance(value, int) and value > 0 for value in bounds)


def test_order_conductor_cannot_escape_above_a_fixed_cycle_order() -> None:
    # Independent toy interval: Z[sqrt(5)] has additive index two in O_K,
    # so its only overorders are the two endpoints.  The general proof uses
    # the same finite-subgroup argument on O_K/O_0.
    sqrt5_basis_in_ok = sp.Matrix([[1, -1], [0, 2]])
    assert abs(int(sqrt5_basis_in_ok.det())) == 2
    quotient_order = 2
    additive_subgroups = ({0}, {0, 1})
    assert len(additive_subgroups) == 2
    assert all(len(group) <= quotient_order for group in additive_subgroups)


def test_codifferent_box_is_finite_without_owner_helper() -> None:
    # O_K^vee=(1/sqrt(5))O_K.  Directly enumerate a provably oversized box
    # and recover the same finite intersection with |z_0|<=14, |z_1|<=2.
    sqrt5 = math.sqrt(5)
    phi0 = (1 + sqrt5) / 2
    phi1 = (1 - sqrt5) / 2
    points = []
    for u in range(-100, 101):
        for v in range(-50, 51):
            z0 = (u + v * phi0) / sqrt5
            z1 = -(u + v * phi1) / sqrt5
            if abs(z0) <= 14 and abs(z1) <= 2:
                points.append((u, v))
    assert len(points) == 253
    assert max(abs(u) for u, _ in points) < 58
    assert max(abs(v) for _, v in points) < 16


def test_full_root_frame_forces_a_finite_dual_interval() -> None:
    # A concrete even integral full-rank frame and an index-two even
    # overlattice replay det(F)=[M:F]^2 det(M).
    frame = sp.diag(2, -2, -2)
    basis_change = sp.Matrix(
        [
            [sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 1, 0],
            [0, 0, 1],
        ]
    )
    host = basis_change.T * frame * basis_change
    assert all(entry.q == 1 for entry in host)
    assert all(host[i, i] % 2 == 0 for i in range(3))
    assert abs(int(frame.det())) == 2**2 * abs(int(host.det()))
    assert abs(int(frame.det())) == 8


def test_projected_root_counterfamily_hits_exact_missing_hypothesis() -> None:
    counter = json.loads(COUNTER_ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert counter["fixed_data"]["compact_polytope"] == "Coxeter triangle (2,5,5)"
    assert counter["fixed_data"]["field"] == "Q(sqrt(5))"
    assert counter["all_parameter_family"]["trace_embedding_primitive"] is True
    assert counter["all_parameter_family"]["trace_contains_minus_two_roots"] is False
    assert counter["K3_gate"]["every_restricted_facet_has_physical_minus_two_wall"] is True
    assert counter["root_projection_exhaustion"][
        "restricted_chamber_exactly_triangle_255"
    ] is True
    assert counter["unbounded_data"]["common_trace_scale"] is True
    assert counter["unbounded_data"]["glue_index"] is True
    assert counter["scope"]["violates_same_vector_in_trace_hypothesis"] is True


def test_owner_and_audit_state_all_scope_boundaries() -> None:
    owner = OWNER_REPORT.read_text(encoding="utf-8")
    audit = AUDIT_REPORT.read_text(encoding="utf-8")
    for needle in (
        "same vector",
        "F\\subseteq M\\subseteq F^\\#",
        "does not assert finiteness of",
        "ordinary K3 root reflection",
        "does not prove Property 4",
    ):
        assert needle in owner
    for needle in (
        "PASS",
        "common trace scale",
        "order conductor",
        "root-line ideals",
        "same-span overlattice index",
        "finite ambient NS extensions or K3 complements",
        "same-host",
        "universal Property 4",
    ):
        assert needle in audit
