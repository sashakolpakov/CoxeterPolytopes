import json

from scripts.certify_h9_14_record35_coxeter_labels import OUTPUT, verify


EXPECTED_SHA256 = (
    "fcb5b0374c1ee68b32193491dcb1d905a7ac39ef61e7c56a7bbfe6bf6b56d23c"
)


def test_frozen_record35_coxeter_label_certificate_replays_exactly() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "record_index": 35,
        "joined_common_core_states": 4_620,
        "rank_five_survivors": 96,
        "compact_coxeter_realization_exists": False,
        "certificate_sha256": EXPECTED_SHA256,
    }


def test_certificate_states_the_unrestricted_scope_and_actual_vertices() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert "unrestricted" in manifest["theorem"]["proof_scope"]
    gate = manifest["census"]["three_support_common_core_cluster"][
        "reusable_common_four_core_gate"
    ]
    assert gate["unrestricted_external_label_no_go"] is True
    assert gate["natural_S3_x_S4_symmetry_preserves_restrictions"] is True
    assert len(gate["actual_vertex_restrictions"]) == 29
    assert [
        (
            stage["order"],
            stage["restriction_count"],
            stage["cumulative_passing_state_count"],
        )
        for stage in gate["cumulative_restriction_stages"]
    ] == [(4, 13, 4_620), (5, 12, 96), (6, 4, 0)]
    assert "supplied actual vertex list" in gate["reuse_scope"]
