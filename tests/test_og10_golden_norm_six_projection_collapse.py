import json

from coxeter_polytopes.lattices import e8_cartan_gram
from scripts.independent_og10_golden_norm_six_projection_collapse import (
    CERTIFICATE_PATH,
    E8_GRAM,
    build_certificate,
    q_rectangle,
    verify_certificate,
)


def test_norm_six_projection_collapse_replays() -> None:
    archived = json.loads(CERTIFICATE_PATH.read_text())
    assert archived == build_certificate()
    assert all(verify_certificate(archived).values())
    assert q_rectangle() == [(1, 0)]


def test_parity_subset_and_v61_signs_are_exact() -> None:
    payload = build_certificate()
    audit = payload["frozen_71_root_prefix_audit"]
    assert audit["parity_subset_indices"] == [8, 9, 61, 70]
    assert audit["parity_subset_indices_through_61"] == [8, 9]
    assert all(row["norm"] == ["2", "0"] for row in audit["rows"])
    assert all(row["on_v61_sign"] == -1 for row in audit["rows"])


def test_literal_e8_matrix_matches_authoritative_project_matrix() -> None:
    assert all(len(row) == 8 for row in E8_GRAM)
    assert E8_GRAM == tuple(tuple(row) for row in e8_cartan_gram())
