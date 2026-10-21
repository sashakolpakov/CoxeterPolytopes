from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    lanner_placement_codes,
)
from scripts.explore_h9_14_record0_a1_dual_facets import build_dual_model
from scripts.explore_h9_14_record0_a1_joint_label_cegar import (
    _branch_stabilizer_actions,
    _digest,
    _move,
    _move_code_on_support,
    a4_label_orbits,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import vertices


EXPECTED_SOURCE1_A4_PARTITION_DIGEST = (
    "df9a46605bdf4899a5994ef4e8ce9a8002f9c1cc0d832341c547eefba433f938"
)


def test_source1_all_108_A4_placements_have_one_orbit_representative() -> None:
    orbits = a4_label_orbits(1)
    assert len(orbits) == 19
    assert sorted(map(len, orbits)) == [3, 3, *([6] * 17)]
    flattened = [code for orbit in orbits for code in orbit]
    assert len(flattened) == 108
    assert len(set(flattened)) == 108
    assert set(flattened) == set(lanner_placement_codes(4))
    assert _digest(orbits) == EXPECTED_SOURCE1_A4_PARTITION_DIGEST
    assert all(
        sum(code in orbit for orbit in orbits) == 1
        for code in lanner_placement_codes(4)
    )


def test_exact_code_stabilizer_actions_preserve_every_branch_model_ground() -> None:
    data = build_dual_model(1, symmetry_breaking=False)
    four = tuple(data["problem"]["fours"])[0]
    support = vertices(four)
    facets = tuple(data["facets"])
    facet_set = set(facets)
    for orbit in a4_label_orbits(1):
        representative = orbit[0]
        actions = _branch_stabilizer_actions(four, representative)
        # Seven adjacent transpositions generate the outside S8.  The two
        # size-three placement orbits have one additional internal involution.
        assert len(actions) == (8 if len(orbit) == 3 else 7)
        for action in actions:
            assert _move(four, action) == four
            assert (
                _move_code_on_support(representative, support, action)
                == representative
            )
            assert {_move(facet, action) for facet in facets} == facet_set


def test_frozen_branch0_checkpoint_has_branch_not_source_scope() -> None:
    path = (
        "data/enumerations/"
        "h9_14_record0_a1_source1_branch00_joint_checkpoint.json"
    )
    document = json.loads(open(path, encoding="utf-8").read())
    assert document["a4_placement_partition"][
        "coverage_is_disjoint_and_complete"
    ] is True
    assert document["a4_placement_partition"]["literal_placement_count"] == 108
    assert document["a4_placement_partition"]["orbit_count"] == 19
    assert document["a4_placement_partition"][
        "partition_digest"
    ] == EXPECTED_SOURCE1_A4_PARTITION_DIGEST
    assert document["a4_branch"]["orbit_index"] == 0
    assert document["a4_branch"]["placement_orbit_size"] == 6
    assert document["seeded_abstract_nonfinite_pattern_count"] == 23
    assert document["seeded_global_monotone_face_cut_count"] == 1769801
    assert document["terminal_solver_status"] == "INFEASIBLE"
    assert document["a4_branch_compact_coxeter_no_go"] is True
    assert document["fixed_source_compact_coxeter_no_go"] is False


def test_every_frozen_source1_branch_checkpoint_is_consistent() -> None:
    paths = sorted(
        Path("data/enumerations").glob(
            "h9_14_record0_a1_source1_branch??_joint_checkpoint.json"
        )
    )
    assert paths
    seen: set[int] = set()
    for path in paths:
        document = json.loads(path.read_text(encoding="utf-8"))
        branch = document["a4_branch"]["orbit_index"]
        assert branch not in seen
        seen.add(branch)
        assert path.name == (
            "h9_14_record0_a1_source1_"
            f"branch{branch:02d}_joint_checkpoint.json"
        )
        assert document["source_orbit_index"] == 1
        assert document["a4_placement_partition"][
            "coverage_is_disjoint_and_complete"
        ] is True
        assert document["a4_placement_partition"][
            "partition_digest"
        ] == EXPECTED_SOURCE1_A4_PARTITION_DIGEST
        assert document["seeded_abstract_nonfinite_pattern_count"] == 23
        assert document["seeded_global_monotone_face_cut_count"] == 1769801
        assert document["installed_abstract_nonfinite_pattern_count"] == 23
        assert document["installed_global_monotone_face_cut_count"] == 1769801
        assert document["completed_cut_round_count"] == 0
        assert document["terminal_solver_status"] == "INFEASIBLE"
        assert document["stopped_reason"] == (
            "terminal_joint_incidence_label_infeasible"
        )
        assert document["compact_label_survivor"] is None
        assert document["a4_branch_compact_coxeter_no_go"] is True
        assert document["fixed_source_compact_coxeter_no_go"] is False
