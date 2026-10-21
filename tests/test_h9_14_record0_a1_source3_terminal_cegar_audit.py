from __future__ import annotations

import ast
from collections import Counter
from functools import lru_cache
import json

from scripts.audit_h9_14_record0_a1_source3_terminal_cegar import (
    CHECKPOINT_SHA256,
    EXPECTED_CORE_COUNT,
    EXPECTED_ORBIT_IMAGE_CLAUSE_COUNT,
    OUTPUT,
    ROOT,
    SUPPORT_SHA256,
    build_certificate,
)


EXPECTED_MODEL_SHA256 = (
    "7485e42496c5f9f8f03cedde468a9f7180b98c9ee92ab6f86a23e68a0b735ec3"
)


@lru_cache(maxsize=1)
def _structural_replay() -> dict:
    # This still rebuilds all 41 label proofs, every source-group orbit, all
    # 427,140 cuts, and the exact terminal model.  Only the costly second
    # CP-SAT infeasibility solve is omitted from the routine regression gate.
    return build_certificate(solve_terminal=False)


@lru_cache(maxsize=1)
def _archived() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def _without_core_times(payload: dict) -> list[dict]:
    rows = []
    for source in payload["local_label_logic"]["core_summaries"]:
        row = dict(source)
        row.pop("orbit_rebuild_wall_time_seconds", None)
        rows.append(row)
    return rows


def test_all_41_local_proofs_and_full_source_orbits_replay() -> None:
    rebuilt = _structural_replay()
    archived = _archived()
    assert rebuilt["pinned_inputs"]["terminal_checkpoint"]["sha256"] == CHECKPOINT_SHA256
    assert rebuilt["pinned_inputs"]["nerve_support_schedule"]["sha256"] == SUPPORT_SHA256
    assert rebuilt["base_model"] == archived["base_model"]
    assert _without_core_times(rebuilt) == _without_core_times(archived)
    rows = rebuilt["local_label_logic"]["core_summaries"]
    assert len(rows) == EXPECTED_CORE_COUNT
    assert all(row["terminal_sparse_label_state_count"] == 0 for row in rows)
    assert all(row["all_checks_pass"] for row in rows)
    assert all(all(row["local_checks"].values()) for row in rows)
    assert all(all(row["orbit_checks"].values()) for row in rows)
    assert sum(row["full_source_core_orbit_size"] for row in rows) == (
        EXPECTED_ORBIT_IMAGE_CLAUSE_COUNT
    )


def test_exact_terminal_model_rebuild_matches_frozen_proto() -> None:
    rebuilt = _structural_replay()
    archived = _archived()
    assert rebuilt["terminal_model"] == archived["terminal_model"]
    terminal = rebuilt["terminal_model"]
    assert terminal["terminal_model_proto_sha256"] == EXPECTED_MODEL_SHA256
    assert terminal["installed_clause_count"] == EXPECTED_ORBIT_IMAGE_CLAUSE_COUNT
    assert terminal["terminal_model_variable_count"] == 10_602
    assert terminal["terminal_model_constraint_count"] == 458_793
    assert terminal["terminal_model_constraint_type_histogram"] == {
        "bool_or": 437_124,
        "linear": 20_122,
        "table": 1_547,
    }
    assert terminal["cut_arity_histogram"] == {
        "5": 5_040,
        "6": 56_784,
        "7": 261_492,
        "8": 55_440,
        "9": 48_384,
    }


def test_terminal_infeasibility_and_scope_are_frozen_exactly() -> None:
    archived = _archived()
    theorem = archived["theorem"]
    solver = archived["terminal_solver"]
    assert archived["all_checks_pass"]
    assert all(archived["checks"].values())
    assert solver["status"] == "INFEASIBLE" and solver["infeasible"]
    assert solver["num_search_workers"] == 1
    assert solver["time_limit_seconds"] is None
    assert solver["solver_wall_time_seconds"] > 500
    assert theorem["record_index"] == 0
    assert theorem["schedule_a"] == 1
    assert theorem["source_orbit_index"] == 3
    assert theorem["representative_A4"] == [2, 3, 4, 13]
    assert theorem["compact_only"]
    assert theorem["finite_volume_noncompact_checked"] is False
    assert theorem["not_global_H9_14"]
    assert archived["decision"]["fixed_source_compact_coxeter_no_go"]
    assert "source orbit 3 only" in archived["decision"]["scope"]


def test_source_is_no_trust_replay_without_exploratory_imports() -> None:
    script = ROOT / "scripts/audit_h9_14_record0_a1_source3_terminal_cegar.py"
    tree = ast.parse(script.read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    assert not any(name.startswith("scripts.explore_h9_14") for name in imports)
    assert Counter(
        row["full_source_core_orbit_size"]
        for row in _archived()["local_label_logic"]["core_summaries"]
    ) == Counter(
        {
            84: 1,
            336: 2,
            1008: 2,
            2016: 16,
            2520: 2,
            4032: 1,
            5040: 2,
            10080: 5,
            20160: 8,
            80640: 2,
        }
    )

