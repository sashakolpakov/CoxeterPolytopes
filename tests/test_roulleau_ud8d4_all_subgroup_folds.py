import json
from pathlib import Path

from scripts.certify_roulleau_ud8d4_all_subgroup_folds import (
    _ambient_reconstruction,
    _source_gram,
)
from scripts.certify_roulleau_ud8d4_fold_lineage import build_certificate
from scripts.search_k3_high_rank_involution_folds import (
    _group_checks,
    _weighted_graph_automorphisms,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = json.loads(
    (ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json").read_text()
)
LINEAGE_PATH = ROOT / "data/enumerations/roulleau_ud8d4_fold_lineage.json"
LINEAGE = json.loads(LINEAGE_PATH.read_text())


def test_source_graph_and_integral_ns_reconstruction() -> None:
    gram = _source_gram()
    basis_gram, wall_rows = _ambient_reconstruction(gram)
    group, _ = _weighted_graph_automorphisms(gram)
    assert gram.rank() == 14
    assert exact_inertia(gram) == (1, 13, 6)
    assert basis_gram.det() == -16
    assert wall_rows * basis_gram * wall_rows.T == gram
    assert len(group) == 240
    assert all(_group_checks(group).values())


def test_complete_subgroup_and_five_type_denominator() -> None:
    assert CONSTRUCTION["schema"] == "coxeter.roulleau_ud8d4_all_subgroup_folds.v1"
    assert CONSTRUCTION["all_checks_pass"] is True
    assert CONSTRUCTION["weighted_diagram_group"]["structure"] == "C2 x S5"
    assert CONSTRUCTION["weighted_diagram_group"]["id_group"] == [240, 189]
    assert CONSTRUCTION["subgroup_census"]["class_count_including_identity_and_full_group"] == 57
    assert CONSTRUCTION["subgroup_census"]["proper_fixed_rank_at_least_11_count"] == 13
    assert [
        (record["hyperbolic_dimension"], record["facet_count"])
        for record in CONSTRUCTION["abstract_diagram_types"]
    ] == [(12, 18), (11, 16), (11, 16), (10, 14), (10, 14)]
    assert all(
        record["geometry"]["finite_volume_certified"]
        and record["geometry"]["compact"] is False
        and record["geometry"]["checks"]["every_orbit_wall_is_irredundant"]
        and record["geometry"]["checks"]["all_full_parabolic_neighborhoods_pass"]
        for record in CONSTRUCTION["abstract_diagram_types"]
    )


def test_exact_corpus_lineage_is_one_known_and_four_bounded_nonmatches() -> None:
    match_counts = [
        record["lineage"]["exact_CoxIter_match_count"]
        for record in CONSTRUCTION["abstract_diagram_types"]
    ]
    assert match_counts == [0, 0, 0, 0, 1]
    known = CONSTRUCTION["abstract_diagram_types"][4]["lineage"]["coxiter"]
    assert [
        row["path"] for row in known["ordinary_exact_matches"]
    ] == ["graphs/group units qf/10-3,1,1,1,1,1,1,1,1,1,1.coxiter"]
    assert known["doubling_exact_matches"] == []


def test_stable_overlattice_and_index_lineage_replays() -> None:
    assert build_certificate() == LINEAGE
    assert LINEAGE["all_checks_pass"] is True
    assert LINEAGE["summary"]["indices"] == {
        "R14_UD8D4_TYPE_1": 8060,
        "R14_UD8D4_TYPE_2": None,
        "R14_UD8D4_TYPE_3": None,
        "R14_UD8D4_TYPE_4": 992,
        "R14_UD8D4_TYPE_5": 1,
    }
    assert [
        record["maximal_overlattice_count"] for record in LINEAGE["types"]
    ] == [3, 1, 3, 1, 1]

