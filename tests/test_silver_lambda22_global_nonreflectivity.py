import json
import subprocess
import sys
from pathlib import Path

from coxeter_polytopes.classification.silver_lambda22_global_nonreflective import (
    BASEPOINT,
    FRONTIER_HEIGHT_BOUND,
    ORBIT_ONE,
    SEPARATING_CONTROL_BOUND,
    SEPARATING_HEIGHT_BOUND,
    TIME_ZERO_NORM_CASES,
    infinite_order_audit,
    inverse_isometry,
    isometry,
    isometry_audit,
    no_time_zero_mirror_audit,
    orbit_recurrence,
    separator_implies_positive_frontier_pairing,
)
from coxeter_polytopes.exact.silver import Silver


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_global_nonreflectivity.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_g_is_an_integral_isometry_with_integral_inverse() -> None:
    assert isometry_audit()
    assert isometry(BASEPOINT) == ORBIT_ONE
    assert inverse_isometry(ORBIT_ONE) == BASEPOINT


def test_g_has_expanding_plane_eigenvalue_and_infinite_order() -> None:
    assert infinite_order_audit()
    assert orbit_recurrence(0) == (Silver(1), Silver(0))
    assert orbit_recurrence(1) == (Silver(5, 4), Silver(2, 1))
    assert orbit_recurrence(2) == (Silver(113, 80), Silver(36, 26))
    assert orbit_recurrence(3) == (Silver(2405, 1700), Silver(774, 547))


def test_separator_height_is_strictly_inside_complete_frontier() -> None:
    assert SEPARATING_CONTROL_BOUND == Silver(8, 6)
    assert SEPARATING_CONTROL_BOUND**2 == Silver(136, 96)
    assert SEPARATING_HEIGHT_BOUND == Silver(56, 40)
    assert SEPARATING_HEIGHT_BOUND < FRONTIER_HEIGHT_BOUND == Silver(150)
    assert separator_implies_positive_frontier_pairing()


def test_complete_norm_alphabet_excludes_time_zero_mirrors() -> None:
    assert len(TIME_ZERO_NORM_CASES) == 8
    assert no_time_zero_mirror_audit()


def test_frozen_payload_states_the_global_conclusion() -> None:
    data = payload()
    assert data["scope"].startswith("Global theorem")
    assert data["frontier_dependency"][
        "base_H14_and_post_H14_complete_through_150"
    ] is True
    assert data["chamber_argument"] == {
        "v0_and_gv0_same_regular_reflection_chamber": True,
        "g_permutes_complete_primitive_reflective_root_set": True,
        "g_normalizes_full_reflection_group": True,
        "g_stabilizes_chamber": True,
        "ambient_integral_isometry_group_discrete": True,
        "infinite_cyclic_chamber_stabilizer": True,
        "finite_volume_chamber_has_finite_stabilizer": True,
    }
    assert data["conclusion"] == {
        "reflection_chamber_has_infinite_volume": True,
        "reflection_subgroup_has_infinite_index": True,
        "host_is_reflective": False,
        "global_nonreflectivity_theorem": True,
    }


def test_frozen_payload_verifies() -> None:
    process = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_silver_lambda22_global_nonreflectivity.py"),
            "--verify",
        ],
        cwd=ROOT,
        env={"PYTHONPATH": str(ROOT / "src")},
        text=True,
        capture_output=True,
        check=False,
    )
    assert process.returncode == 0, process.stdout + process.stderr
    assert process.stdout.strip() == "PASS"
