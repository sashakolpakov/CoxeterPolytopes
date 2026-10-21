import json
from pathlib import Path

from scripts.certify_roulleau_ue8d4_generic_ns_k3_automorphisms import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/roulleau_ue8d4_generic_ns_k3_automorphisms.json"
)


def test_archive_rebuilds_exactly_offline() -> None:
    assert build_certificate() == json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_two_genuine_generic_ns_realizations() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert payload["summary"] == {
        "compact": 0,
        "dimensions": [11, 12],
        "finite_volume_noncompact": 2,
        "genuine_k3_automorphism_realizations": 2,
        "new_abstract_polytope_claims": 0,
    }
    assert [
        (
            row["brandhorst_hofmann_row"]["id"],
            row["brandhorst_hofmann_row"]["automorphism_order"],
            row["exact_picard_action_match"]["picard_action_order"],
            row["exact_picard_action_match"]["archived_type_id"],
        )
        for row in payload["cases"]
    ] == [
        ("0.6.3.1", 6, 3, "R14_UE8D4_C3_S3"),
        ("0.16.0.1", 16, 2, "R14_UE8D4_C2"),
    ]
    assert all(
        row["brandhorst_hofmann_row"]["generic_picard_lattice"][
            "isometry_class"
        ]
        == "U+E8(-1)+D4(-1)"
        and row["brandhorst_hofmann_row"]["generic_picard_lattice"][
            "rank_signature_determinant"
        ]
        == [14, [1, 13, 0], -4]
        for row in payload["cases"]
    )


def test_exact_action_and_kahler_gates() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert payload["chamber_transport_criterion"][
        "diagram_action_on_discriminant_is_faithful"
    ]
    for row in payload["cases"]:
        match = row["exact_picard_action_match"]
        discr = match["discriminant_form_and_action"]
        assert match["equivariant_orthogonal_decomposition"][
            "actions_conjugate_integrally"
        ]
        assert discr["same_q_action_class"]
        assert discr["archived"]["q_histogram_mod_2Z"] == {"0": 1, "1": 3}
        assert discr["brandhorst_hofmann"]["q_histogram_mod_2Z"] == {
            "0": 1,
            "1": 3,
        }
        assert match["coinvariant_lattice"]["rootless"]
        assert row["kahler_and_torelli"][
            "database_row_is_effective_genuine_automorphism_moduli"
        ]
        assert row["kahler_and_torelli"]["generic_picard_is_exactly_rank_14"]


def test_vinberg_finite_volume_handoff() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert [
        (
            row["hyperbolic_polytope"]["dimension"],
            row["hyperbolic_polytope"]["facets"],
            row["hyperbolic_polytope"]["finite_vertices"],
            row["hyperbolic_polytope"]["cusps"],
            row["hyperbolic_polytope"]["elliptic_edges"],
        )
        for row in payload["cases"]
    ] == [(11, 13, 27, 1, 162), (12, 14, 35, 2, 234)]
    assert all(
        row["hyperbolic_polytope"]["classification"]
        == "finite_volume_noncompact"
        and row["hyperbolic_polytope"]["bad_edges"] == 0
        and row["hyperbolic_polytope"]["full_reflection_chamber_index"] == 1
        for row in payload["cases"]
    )
