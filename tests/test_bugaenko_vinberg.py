import copy
import json
from pathlib import Path

import pytest

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.known_examples.bugaenko import (
    bugaenko_h8_certificate,
    bugaenko_simple_roots,
)
from coxeter_polytopes.search.bugaenko_vinberg import (
    admissible_control_coordinates,
    build_bugaenko_vinberg_transcript,
    e8_roots_from_simple_reflections,
    run_bugaenko_vinberg_search,
)
from coxeter_polytopes.verify.vinberg import (
    VinbergTranscriptError,
    verify_bugaenko_vinberg_completion,
    verify_bugaenko_vinberg_transcript,
)


TRANSCRIPT = (
    Path(__file__).parents[1]
    / "data"
    / "known_examples"
    / "bugaenko_vinberg_transcript.json"
)


def _archived_transcript() -> dict:
    return json.loads(TRANSCRIPT.read_text())


def test_height_zero_e8_enumeration_is_exhaustive() -> None:
    roots = e8_roots_from_simple_reflections()
    assert len(roots) == 240
    assert len(set(roots)) == 240
    assert all(
        all(coefficient >= 0 for coefficient in root)
        or all(coefficient <= 0 for coefficient in root)
        for root in roots
    )


def test_all_control_coordinates_through_stopping_height() -> None:
    assert admissible_control_coordinates(7) == (
        Golden(1),
        PHI,
        1 + PHI,
        2 * PHI,
        1 + 2 * PHI,
        2 + 2 * PHI,
        1 + 3 * PHI,
        2 + 3 * PHI,
    )


def test_exact_vinberg_search_reproduces_eleven_walls_without_root_input() -> None:
    search = run_bugaenko_vinberg_search()
    assert [len(level.roots) for level in search.levels] == [1, 1, 0, 0, 0, 0, 0, 1]
    assert len(search.roots) == 11
    assert set(search.roots) == {tuple(root) for root in bugaenko_simple_roots()}
    assert [level.height for level in search.levels] == sorted(
        level.height for level in search.levels
    )


def test_archived_transcript_is_a_clean_reproducible_run() -> None:
    archived = _archived_transcript()
    assert archived == build_bugaenko_vinberg_transcript()
    search = verify_bugaenko_vinberg_transcript(archived)
    assert len(search.roots) == 11
    assert sum(sum(level.visited_by_depth) for level in search.levels) == 21_547


def test_transcript_and_independent_geometry_prove_vinberg_completion() -> None:
    verification = verify_bugaenko_vinberg_completion(
        _archived_transcript(), bugaenko_h8_certificate()
    )
    assert len(verification.search.roots) == 11
    assert verification.chamber.facet_count == 11
    assert len(verification.chamber.finite_vertices) == 41
    assert verification.chamber.cusps == ()
    assert verification.diagram.decision == "compact"
    assert verification.diagram.elliptic_edge_count == 164
    assert verification.diagram.finite_vertex_count == 41


def test_transcript_verifier_rejects_a_falsified_empty_height() -> None:
    falsified = copy.deepcopy(_archived_transcript())
    falsified["levels"][2]["solution_count"] = 1
    with pytest.raises(VinbergTranscriptError, match="solution_count"):
        verify_bugaenko_vinberg_transcript(falsified)
