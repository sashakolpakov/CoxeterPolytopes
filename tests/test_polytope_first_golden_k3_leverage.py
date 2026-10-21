import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from scripts.certify_polytope_first_golden_k3_leverage import (
    EXPECTED_TRACE_SMITH,
    OUTPUT,
    SEEDS,
    build,
    field_gram,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]


def test_eight_seed_grams_are_exact_compact_lanner_data() -> None:
    payload = build()
    rows = payload["lanner_seeds"]
    assert len(rows) == 8
    assert [row["name"] for row in rows] == [seed[0] for seed in SEEDS]
    assert sorted({row["dimension"] for row in rows}) == [2, 3, 4]
    for row in rows:
        dimension = row["dimension"]
        assert row["identity_inertia"] == [dimension, 1, 0]
        assert row["conjugate_inertia"] == [dimension + 1, 0, 0]
        assert row["vertex_link_inertias"] == [
            [dimension, 0, 0]
        ] * (dimension + 1)
        assert Golden.from_json(row["field_determinant"]).sign(0) < 0
        assert Golden.from_json(row["conjugate_field_determinant"]).sign(0) > 0


def test_trace_smith_stable_range_and_periods() -> None:
    rows = build()["lanner_seeds"]
    for row in rows:
        assert tuple(row["trace_smith"]) == EXPECTED_TRACE_SMITH[row["name"]]
        assert row["trace_rank"] == 2 * (row["dimension"] + 1)
        assert row["k3_complement_rank"] == 20 - 2 * row["dimension"]
        assert row["k3_complement_rank"] >= row["trace_rank"] + 2
        assert row["k3_period_dimension"] == 18 - 2 * row["dimension"]
        assert row["checks"]["trace_even_integral"] is True
        assert row["checks"]["trace_signature"] is True
        assert row["checks"]["nikulin_stable_range"] is True


def test_paired_root_arithmetic_is_a1_plus_a1() -> None:
    paired = build()["paired_root_arithmetic"]
    assert paired["root_trace_square"] == -2
    assert paired["unit_partner_trace_square"] == -2
    assert paired["trace_pairing"] == 0
    assert paired["unit_partner"] == "phi-1"


def test_gram_constructor_keeps_omitted_edges_orthogonal() -> None:
    gram = field_gram(3, ((0, 1, 5), (1, 2, 3), (2, 3, 5)))
    assert gram[0][2] == 0
    assert gram[0][3] == 0
    assert gram[1][3] == 0
    assert gram[0][1] == -Golden(0, 1)
    assert gram[1][2] == -1


def test_high_dimensional_scope_is_not_overclaimed() -> None:
    payload = build()
    boundary = payload["high_dimensional_boundary"]
    assert boundary["bugaenko_H8"] == {
        "trace_rank": 18,
        "discriminant_length": 2,
        "stable_bound_right_side": 2,
        "period_dimension": 2,
        "complete_chamber": "compact_11_facets",
    }
    assert boundary["silver_H9"]["complete_chamber_volume"] == "infinite"
    assert boundary["full_trace_quadratic_k3_ceiling"] == 9
    assert "not arithmetic Coxeter polytopes themselves" in payload[
        "theorem_scope"
    ]["nonexistence_scope"]
    embedding_scope = payload["theorem_scope"]["embedding_orbit_scope"]
    assert "complement isometry classes" in embedding_scope
    assert "anti-isometry double cosets" in embedding_scope
    wall_scope = payload["theorem_scope"]["trace_even_wall_scope"]
    assert "five-shell list" in wall_scope
    assert "restricted to field-even packages" in wall_scope
    compactness_scope = payload["theorem_scope"]["compactness_bound_scope"]
    assert "sigma(t)" in compactness_scope
    assert "not uniformly 2" in compactness_scope


def test_frozen_certificate_replay() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(archived)
    assert summary["seed_count"] == 8
    assert summary["dimensions"] == [2, 3, 4]
    assert summary["full_trace_quadratic_k3_ceiling"] == 9
