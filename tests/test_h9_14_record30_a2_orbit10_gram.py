from __future__ import annotations

from fractions import Fraction
import json

import sympy as sp

from scripts.certify_h9_14_record30_a2_orbit10_gram import (
    ARCHIVE,
    _finite_labels,
    _forest_matching_expansion,
    _twice_gram,
    build,
    verify,
)


EXPECTED_SHA256 = (
    "a3b5996c71cc616fb18dd99f6c54befab01196d9b7dfb23b964ceb2d7f7eef17"
)


def test_record30_a2_orbit10_exact_gram_obstruction_replays() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(archived) == {
        "rank_nine_vertices": 128,
        "principal_order": list(range(11)),
        "principal_determinant": "7*(3 - sqrt(5))",
        "full_gram_rank_lower_bound": 11,
        "H9_gram_rank_upper_bound": 10,
        "this_exact_labeling_excluded": True,
        "certificate_sha256": EXPECTED_SHA256,
    }
    assert build() == archived


def test_decisive_minor_has_two_independent_exact_derivations() -> None:
    principal = _twice_gram(tuple(range(11)), _finite_labels())
    assert sp.simplify(principal.det() - 7 * (3 - sp.sqrt(5))) == 0
    matching_determinant, matching_sums = _forest_matching_expansion()
    assert matching_determinant == (Fraction(21), Fraction(-7))
    assert matching_sums == (
        (Fraction(1), Fraction(0)),
        (Fraction(10), Fraction(1)),
        (Fraction(35), Fraction(7)),
        (Fraction(99, 2), Fraction(29, 2)),
        (Fraction(51, 2), Fraction(19, 2)),
        (Fraction(7, 2), Fraction(3, 2)),
    )
    assert 7 * (3 - sp.sqrt(5)) > 0


def test_scope_does_not_exclude_the_nerve_or_schedule() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived["local_replay"] == {
        "all_order_four_and_five_minimal_nonfaces_lanner": True,
        "all_rank_nine_vertex_diagrams_finite": True,
        "bad_higher_lanner_supports": [],
        "bad_rank_nine_vertices": [],
        "warning": (
            "Passing these local tests is necessary but not sufficient; "
            "the global Gram-rank test below is decisive."
        ),
    }
    scope = archived["scope"]
    assert scope["single_sparse_terminal_state"] is True
    assert scope["all_values_of_the_three_dotted_entries"] is True
    assert scope["other_sparse_terminal_states_excluded"] is False
    assert scope["other_cross_edge_placements_excluded"] is False
    assert scope["whole_orbit10_nerve_excluded"] is False
    assert scope["whole_record30_a2_schedule_excluded"] is False
    assert scope["finite_volume_noncompact_statement"] is False
