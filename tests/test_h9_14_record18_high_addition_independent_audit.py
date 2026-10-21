from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HIGH_ARCHIVE = (
    ROOT
    / "data/enumerations/h9_14_record18_high_addition_independent_audit.json"
)
A1_DOWNSTREAM_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_a1_downstream_independent_audit.json"
)
A1_DOWNSTREAM_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_a1_downstream_audit.py"
)
A2_DOWNSTREAM_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_a2_downstream_independent_audit.json"
)
A2_DOWNSTREAM_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_a2_downstream_audit.py"
)


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def test_high_addition_solver_independent_census() -> None:
    record = load(HIGH_ARCHIVE)
    assert "128-bit" in record["backend"]
    cases = {
        tuple(case["representative_A4_masks"]): case
        for case in record["cases"]
    }
    assert set(cases) == {
        (53, 197),
        (53, 54),
        (53, 86),
        (53, 198),
        (23, 101),
        (53, 58),
        (53, 90),
        (53, 202),
        (23, 57),
        (23, 105),
        (15, 53),
        (23, 27),
        (23, 43),
        (23, 45),
        (23,),
        (15,),
        (53,),
    }
    for key in ((53, 54), (23, 27), (23,), (15,)):
        case = cases[key]
        assert case["scheduled_ridge_terminal_count"] == 0
        assert case["literal_pure_target_count"] == 0
        assert case["compact_incidence_unsat"] is True
    target = cases[(53,)]
    assert target["scheduled_ridge_terminal_count"] == 8
    assert target["literal_pure_target_count"] == 8
    assert target["canonical_target_count"] == 8
    assert sorted(target["canonical_target_orbit_sizes"]) == [
        20,
        40,
        60,
        120,
        180,
        180,
        240,
        240,
    ]
    assert target["literal_target_orbit_weight_sum"] == 1080
    assert cases[(53, 198)]["canonical_target_orbit_sizes"] == [6]
    assert cases[(53, 202)]["canonical_target_orbit_sizes"] == [1, 1]
    for key, case in cases.items():
        if key in {(53,), (53, 198), (53, 202)}:
            continue
        assert case["canonical_target_count"] == 0


def test_a1_downstream_topology_and_label_terminals() -> None:
    record = load(A1_DOWNSTREAM_ARCHIVE)
    assert record["raw_target_weight"] == 1080
    assert record["canonical_class_count"] == 8
    assert record["class_orbit_sizes"] == [20, 40, 60, 120, 180, 180, 240, 240]
    assert record["all_classes_pass_connected_mod2_homology_audit"] is True
    assert record["all_classes_have_unrestricted_compact_coxeter_no_go"] is True
    observed = sorted(
        (
            audit["literal_orbit_size"],
            tuple(
                stage["passing_state_count"]
                for stage in audit["coxeter_labels"]["stages"]
            ),
        )
        for audit in record["classes"]
    )
    assert observed == sorted(
        (
            (180, (108, 84, 8, 0)),
            (180, (108, 84, 24, 16, 72, 0)),
            (240, (108, 84, 24, 16, 0)),
            (240, (108, 84, 8, 0)),
            (40, (108, 0)),
            (120, (108, 84, 8, 0)),
            (60, (108, 84, 8, 0)),
            (20, (108, 0)),
        )
    )
    for audit in record["classes"]:
        topology = audit["topology"]
        assert topology["face_vector"] == [
            14,
            86,
            310,
            734,
            1198,
            1358,
            1036,
            486,
            108,
        ]
        assert topology["ridge_degree_histogram"] == {"2": 486}
        assert topology["passes_connected_homology_audit"] is True
        assert audit["coxeter_labels"]["terminal_state_count"] == 0


def test_a1_downstream_replay_matches_frozen_archive() -> None:
    specification = importlib.util.spec_from_file_location(
        "record18_a1_downstream_independent", A1_DOWNSTREAM_SOURCE
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    replay = module.build_record(HIGH_ARCHIVE)
    frozen = load(A1_DOWNSTREAM_ARCHIVE)
    replay["source"] = frozen["source"]
    assert json.loads(json.dumps(replay, sort_keys=True)) == frozen


def test_a2_downstream_topology_and_label_terminals() -> None:
    record = load(A2_DOWNSTREAM_ARCHIVE)
    assert record["a2_representative_count"] == 14
    assert record["a2_add4_target_raw_weight"] == 8
    assert record["a2_add4_target_subgroup_class_count"] == 3
    assert record["a2_add4_literal_family_count_audited"] == 8
    assert record["a2_compact_branch_excluded"] is True
    assert record["all_target_classes_pass_connected_mod2_homology_audit"] is True
    assert (
        record["all_target_classes_have_unrestricted_compact_coxeter_no_go"]
        is True
    )
    observed = sorted(
        (
            tuple(audit["A4_masks"]),
            tuple(
                stage["passing_state_count"]
                for stage in audit["coxeter_labels"]["stages"]
            ),
        )
        for audit in record["classes"]
    )
    assert observed == sorted(
        [((53, 198), (108, 84, 4, 28, 196, 0))] * 6
        + [((53, 202), (108, 84, 8, 16, 0))] * 2
    )
    for audit in record["classes"]:
        topology = audit["topology"]
        assert topology["face_vector"] == [
            14,
            86,
            310,
            733,
            1192,
            1344,
            1020,
            477,
            106,
        ]
        assert topology["ridge_degree_histogram"] == {"2": 477}
        assert topology["passes_connected_homology_audit"] is True
        assert audit["coxeter_labels"]["terminal_state_count"] == 0


def test_a2_downstream_replay_matches_frozen_archive() -> None:
    specification = importlib.util.spec_from_file_location(
        "record18_a2_downstream_independent", A2_DOWNSTREAM_SOURCE
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    replay = module.build_record(HIGH_ARCHIVE)
    frozen = load(A2_DOWNSTREAM_ARCHIVE)
    replay["source"] = frozen["source"]
    replay["raw_family_source"] = frozen["raw_family_source"]
    assert json.loads(json.dumps(replay, sort_keys=True)) == frozen
