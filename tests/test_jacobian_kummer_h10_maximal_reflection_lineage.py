from __future__ import annotations

import json
from pathlib import Path

from coxeter_polytopes.search.jk_h10_maximal_vinberg import (
    divisibility,
    pair,
    run_search,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data" / "enumerations" / "jacobian_kummer_h10_maximal_reflection_lineage.json"


def test_exact_vinberg_stopping_and_root_types() -> None:
    search = run_search(3, attempt_stopping=True)
    assert len(search.roots) == 13
    assert search.stopping.fired
    assert search.stopping.finite_vertex_count == 29
    assert search.stopping.cusp_count == 3
    assert search.stopping.hyperideal_ray_count == 0
    assert search.stopping.opposite_sheet_ray_count == 0
    assert sorted((pair(root, root), divisibility(root)) for root in search.roots) == [
        (-4, 2),
    ] * 11 + [(-4, 4), (-2, 1)]


def test_archived_exact_corpus_identification() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"]
    assert payload["vinberg_maximal_integral_reflection_chamber"]["facet_count"] == 13
    matches = payload["coxiter_snapshot"]["exact_matches"]
    assert [match["path"] for match in matches] == [
        "graphs/group units qf/10-2,1,1,1,1,1,1,1,1,1,1.coxiter"
    ]
    assert payload["lineage_decision"]["induced_84_wall_group_is_a_finite_index_reflection_subgroup"]


def test_reflection_doubling_archive_has_no_84_wall_h10_record() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    archive = payload["coxiter_snapshot"]["reflection_doubling_archive"]
    assert archive["coxiter_file_count"] == 2_959
    assert archive["h10_facet_count_histogram"] == {
        "12": 3,
        "13": 4,
        "14": 3,
        "15": 1,
        "16": 3,
        "19": 1,
        "20": 3,
        "28": 2,
        "44": 2,
    }
    assert archive["h10_84_facet_records"] == []
