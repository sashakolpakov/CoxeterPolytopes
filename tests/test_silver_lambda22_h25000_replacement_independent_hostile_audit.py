"""Frozen checks for the clean-room silver Lambda22 H25000 replay."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/silver_lambda22_h25000_replacement_independent_hostile_audit.json"
)
SOURCE = (
    ROOT
    / "scripts/independent_hostile_audit_silver_lambda22_h25000_replacement.py"
)


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_clean_room_frontier_agrees_exactly() -> None:
    frozen = payload()
    frontier = frozen["frontier"]
    assert frozen["verdict"] == "PASS"
    assert frontier["control_count"] == 301
    assert frontier["retained_type_count"] == 5_575
    assert frontier["locally_admissible_extension_count_with_multiplicity"] == 2_663_705
    assert frontier["nonnegative_type_count"] == 2
    assert frontier["retained_type_stream_sha256"] == (
        "21240a46b23221f1c167ffd5c6267c71eb328a57b0a067fa70e488b3955f2e5f"
    )
    assert frozen["source"]["all_frontier_counts_and_digest_agree"] is True


def test_clean_room_control_box_and_threat_reductions() -> None:
    frozen = payload()
    proof = frozen["cutoff_semantics"]["coefficient_box_proof"]
    assert proof == {
        "identity_embedding_absolute_bound": "<1118",
        "conjugate_embedding_absolute_bound": "<3",
        "derived_rational_coefficient_bound": "<=560",
        "derived_radical_coefficient_bound": "<=396",
        "enumerated_symmetric_coefficient_box": 561,
    }
    first, second = frozen["critical_reductions"]
    assert (first["forced_relation"], first["short_vector_norm"]) == ("v=40d", 4)
    assert (first["gate_sublattice"], first["gate_sublattice_minimum"]) == (
        "S2",
        12,
    )
    assert (second["forced_relation"], second["short_vector_norm"]) == ("v=60d", 6)
    assert (second["gate_sublattice"], second["gate_sublattice_minimum"]) == (
        "S3",
        24,
    )
    assert all(row["empty"] for row in (first, second))


def test_auditor_imports_no_project_search_code_and_keeps_bounded_scope() -> None:
    source = SOURCE.read_text()
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    scope = payload()["logical_scope"]
    assert scope["all_controls_reconstructed"] is True
    assert scope["all_integral_PSD_extensions_enumerated"] is True
    assert scope["bounded_prefix_only"] is True
    assert scope["vinberg_termination"] is False
    assert scope["global_nonreflectivity"] is False
