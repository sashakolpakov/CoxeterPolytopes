import json

from coxeter_polytopes.combinatorics.h9_14_record6 import (
    EXPECTED_NUMERICAL_BRANCHES,
    audit_ten_branches,
    derive_numerical_branches,
    verify_target_graph,
)
from scripts.certify_h9_14_record6 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)
from scripts.independent_h9_14_record6_audit import (
    EXPECTED_SUMMARY_SHA256,
    audit as hostile_audit,
)


def test_record6_exact_ten_branch_exclusion() -> None:
    assert verify_target_graph() is True
    assert derive_numerical_branches() == EXPECTED_NUMERICAL_BRANCHES
    branches = audit_ten_branches()
    assert len(branches) == 10
    assert all(branch.excluded for branch in branches)
    by_name = {branch.name: branch for branch in branches}
    assert (
        by_name["EEE_z0_a0_0"].link_compatible_candidates,
        by_name["EEE_z0_a0_0"].independent_triple_compatible_candidates,
        by_name["EEE_z0_a0_0"].expected_a0_compatible_candidates,
        by_name["EEE_z0_a0_0"].isolate_four_histogram,
    ) == (54, 36, 0, ((1, 36),))
    assert by_name["FHH_z1_a0_0"].missing_triple_histogram == ((1, 1440),)
    assert by_name["HHH_z3_a0_0"].missing_triple_histogram == ((1, 780),)
    assert by_name["EHH_z0_a0_1"].missing_triple_histogram == ((1, 60),)


def test_frozen_record6_incidence_no_go_certificate() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(manifest) == {
        "fixed_base_link_type_count": 42,
        "numerical_branch_count": 10,
        "excluded_branch_count": 10,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert manifest["theorem_status"] == "certified_incidence_impossible"
    assert all(row["excluded"] for row in manifest["branch_audits"])


def test_record6_independent_hostile_replay() -> None:
    summary = hostile_audit()
    assert summary["imports_production_module"] is False
    assert summary["production_counts_match"] is True
    assert summary["all_ten_branches_excluded"] is True
    assert summary["summary_sha256"] == EXPECTED_SUMMARY_SHA256
