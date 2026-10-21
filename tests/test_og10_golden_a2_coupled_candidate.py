import json
from pathlib import Path

from scripts.certify_og10_golden_a2_coupled_candidate import (
    certify,
    projected_wall_y_values,
)


def test_a2_coupled_og10_exact_carrier_and_volume_screen() -> None:
    result = certify()
    archived = json.loads(
        Path(
            "data/candidates/og10_golden_a2_coupled_candidate.json"
        ).read_text()
    )

    assert result == archived
    assert result["verified"] is True
    assert all(result["checks"].values())

    field = result["field_lattice"]
    assert field["determinant"] == "-2 - 6*phi"
    assert field["determinant_norm"] == -20
    assert field["identity_signature"] == [10, 1]
    assert field["conjugate_signature"] == [11, 0]

    trace = result["trace_lattice"]
    assert trace["signature"] == [1, 21]
    assert trace["determinant"] == -20
    assert trace["smith_nonunit_invariants"] == [2, 10]

    embedding = result["primitive_og10_embedding"]
    assert embedding["status"] == (
        "proved_by_explicit_gluing_and_integral_isometry"
    )
    assert embedding["positive_complement_gram"] == [[2, 0], [0, 30]]
    assert embedding["n_embedding_smith_diagonal"] == [1] * 6
    assert embedding["k_embedding_smith_diagonal"] == [1] * 2
    assert embedding["full_ambient_type"] == "3U + 2E8(-1) + A2(-1)"

    local = result["local_fold_diagnostic"]
    assert local["root_pairing"] == "-1*phi"
    assert local["root_pair_coxeter_label"] == 5
    assert local["fixed_projection_norm"] == "1 - 1/2*phi"
    assert local["folded_fixed_norm_two_root"][:3] == ["-1", "0", "1*phi"]

    projected = result["projected_wall_theorem"]
    assert projected["possible_middle_coordinates"] == [
        "-1",
        "-1 + 1*phi",
        "0",
        "1 - 1*phi",
        "1",
    ]
    assert len(projected_wall_y_values()) == 5
    assert "if_and_only_if" in projected["arrangement_conclusion"]

    diagnostics = result["fixed_symmetry_lift_diagnostics"]
    assert diagnostics["bugaenko_reversal_extension_status"] == (
        "both_extensions_are_Z[phi]_integral_isometric_involutions"
    )
    assert diagnostics["pell_fold_extension_status"] == (
        "both_field_isometric_involutions_are_nonintegral"
    )
    assert [len(entries) for entries in diagnostics["pell_fold_nonintegral_entries"]] == [
        2,
        2,
    ]

    decision = result["norm_two_volume_decision"]
    assert decision["ambient_hyperbolic_dimension"] == 10
    assert decision["ambient_norm_two_chamber_volume"] == "infinite"
    assert "mixed_Kahler_wall_shell_audit" in decision["scope"]

