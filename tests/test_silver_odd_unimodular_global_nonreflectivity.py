import json
import subprocess
import sys
from pathlib import Path

from coxeter_polytopes.classification.silver_odd_unimodular_global_nonreflective import (
    CATALOGUE_HOSTS,
    SEPARATOR_HEIGHT_BOUND,
    boost_audit,
    plane_isometry_audit,
)
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_s26 import minimal_direction_audit_through


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_odd_unimodular_global_nonreflectivity.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_exact_boost_and_separator() -> None:
    assert boost_audit()
    assert plane_isometry_audit()
    assert SEPARATOR_HEIGHT_BOUND == Silver(543, 384)


def test_universal_frontier_is_strict_at_exact_cutoff() -> None:
    audit = minimal_direction_audit_through(SEPARATOR_HEIGHT_BOUND)
    assert (
        audit.control_count,
        audit.distinct_height_count,
        audit.nonempty_control_count,
        audit.gated_binary_type_count,
        audit.nonempty_height_count,
    ) == (139, 108, 32, 297, 26)
    assert audit.strict_hyperideal
    assert audit.worst_row.witness_margin == Silver(-51, 37) / 55


def test_frozen_scope_is_universal_but_not_field_wide() -> None:
    data = payload()
    assert data["hypotheses"] == {
        "spatial_lattice": "positive-definite integral unimodular",
        "minimum": 3,
        "chosen_direction_square": 3,
        "rank_restriction": None,
        "shell_or_orbit_assumption": False,
    }
    assert data["conclusion"][
        "universal_minimum_three_unimodular_family_nonreflective"
    ] is True
    assert data["conclusion"]["field_wide_nonexistence"] is False
    assert data["conclusion"]["named_catalogue_instances_nonreflective"] == [
        name for name, _ in CATALOGUE_HOSTS
    ]
    corollary = data["published_classification_corollary"]
    assert corollary["rank_counts"] == {
        "23": 1,
        "24": 1,
        "25": 0,
        "26": 1,
        "27": 3,
        "28": 38,
        "29": 10092,
    }
    assert corollary["total"] == 10136
    assert "not all arithmetic hosts" in corollary["scope"]
    local_proof = data["complete_reflective_arithmetic"]["proof"]
    assert "(y,L)" in local_proof and "pL" in local_proof
    assert "S26" not in local_proof


def test_frozen_chamber_to_infinite_volume_bridge_is_explicit() -> None:
    chamber = payload()["chamber_logic"]
    assert chamber["reflection_subgroup_normal"] is True
    assert chamber["reflection_subgroup_chamber_stabilizer_trivial"] is True
    assert chamber["boost_powers_inject_into_Gamma_mod_W"] is True
    assert chamber["conjugate_form_positive_definite"] is True
    assert chamber["ambient_arithmetic_quotient_cocompact"] is True
    assert chamber["reflection_chamber_infinite_volume"] is True
    assert chamber["reflection_subgroup_infinite_index"] is True


def test_frozen_payload_verifies() -> None:
    process = subprocess.run(
        [
            sys.executable,
            str(
                ROOT
                / "scripts/audit_silver_odd_unimodular_global_nonreflectivity.py"
            ),
            "--verify",
        ],
        cwd=ROOT,
        env={"PYTHONPATH": str(ROOT / "src")},
        text=True,
        capture_output=True,
        check=False,
    )
    assert process.returncode == 0, process.stdout + process.stderr
    assert process.stdout.strip() == "PASS_GLOBAL_ODD_UNIMODULAR_NONREFLECTIVITY"
