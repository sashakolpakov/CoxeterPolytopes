from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_p17_p19_three_branch_property4_escape import ARCHIVE, build


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/p17_p19_three_branch_property4_escape.md"


def frozen() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text())


def test_live_certificate_matches_archive() -> None:
    assert build() == frozen()


def test_all_three_current_transfer_branches_are_empty() -> None:
    payload = frozen()
    assert payload["status"] == "PASS_CURRENT_THREE_BRANCH_COMPILER_EMPTY_ON_P17_P18_P19"
    assert [row["n"] for row in payload["rows"]] == [17, 18, 19]
    for row in payload["rows"]:
        assert row["direct_rational_full_nef"] == "EMPTY_BY_MARKING_INDEPENDENT_MMP_ABEL"
        assert row["weak_log_rational"]["status"] == "EMPTY_AT_ADJUNCTION_LINEAR_SYSTEM"
        assert row["literal_full_rank_K3"]["status"] == "EMPTY_AT_PHYSICAL_ROOT_NORMALIZATION"


def test_weak_log_boundary_adjunction_has_no_candidate() -> None:
    expected_markings = {17: 8, 18: 24, 19: 120}
    expected_colourings = {17: 4, 18: 32768, 19: 32}
    expected_first_failure = {17: {"19": 4}, 18: {"20": 32768}, 19: {"21": 32}}
    for row in frozen()["rows"]:
        assert row["canonical_marking_count"] == expected_markings[row["n"]]
        weak = row["weak_log_rational"]
        assert weak["colorings_tested"] == expected_colourings[row["n"]]
        assert weak["consistent_colorings"] == []
        assert len(weak["basis_facet_labels"]) == row["rho"]
        assert weak["max_initial_equations_satisfied"] == row["rho"]
        assert weak["first_failed_facet_histogram"] == expected_first_failure[row["n"]]


def test_mixed_norm_classes_poison_literal_k3_physical_roots() -> None:
    for row in frozen()["rows"]:
        k3 = row["literal_full_rank_K3"]
        assert k3["odd_host"] is True
        assert k3["facet_norm_square_classes"] == [1, 2]
        assert k3["rational_common_similarity_obstruction"] == "2 is not a square in Q"
        assert row["norm_one_facet_labels"]
        assert row["norm_two_facet_labels"]


def test_scope_does_not_promote_branch_emptiness_to_polytope_no_go() -> None:
    scope = frozen()["scope"]
    assert scope["proved"].startswith("emptiness of the three currently declared")
    assert "nonexistence of P17, P18, or P19 as arithmetic Coxeter chambers" in scope["not_claimed"]
    assert "nonexistence of every possible carrier on a nonrational surface or higher-dimensional variety" in scope[
        "not_claimed"
    ]
    text = REPORT.read_text()
    assert "not a no-surface theorem" in text
    assert "another carrier species" in text
    assert "uncontrolled" in text


def test_report_states_the_exact_property4_consequence() -> None:
    text = REPORT.read_text()
    assert "finite-volume complete chamber in the pinned host" in text
    assert "nonemptiness of the current three carrier branches" in text
    assert "same-module/same-line physical-root contract" in text
    assert "positive-characteristic obstruction" in text
    assert "not by itself a strict" in text
