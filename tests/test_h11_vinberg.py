import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.search.h11_vinberg import (
    D,
    admissible_control_coordinates,
    build_h11_vinberg_transcript,
    h11_basis_gram,
    run_h11_vinberg_extension,
    run_h11_vinberg_search,
    standard_e8_roots_from_simple_reflections,
)
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from coxeter_polytopes.verify.vinberg import verify_h11_vinberg_transcript


TRANSCRIPT = (
    Path(__file__).parents[1] / "data" / "candidates" / "h11_vinberg_prefix.json"
)
BOUND20_TRANSCRIPT = (
    Path(__file__).parents[1]
    / "data"
    / "candidates"
    / "h11_vinberg_prefix_bound20.json"
)
BOUND5_HYPERIDEAL_WITNESS = (0, 1, 2, 4, 5, 6, 8, 10, 12)
BOUND20_HYPERIDEAL_WITNESS = (0, 3, 4, 5, 6, 7, 8, 11, 85)


def _decode_vector(value):
    return tuple(Golden.from_json(entry) for entry in value)


def test_control_scan_is_exhaustive_through_bound_five() -> None:
    assert admissible_control_coordinates(5) == (
        PHI - 1,
        Golden(1),
        PHI,
        1 + PHI,
        2 * PHI,
        2 + PHI,
        1 + 2 * PHI,
    )


def test_height_zero_system_and_bounded_wall_counts() -> None:
    assert len(standard_e8_roots_from_simple_reflections()) == 240
    search = run_h11_vinberg_search(5)
    assert search.cone_root_count == 242
    assert len(search.cone_roots) == 9
    assert [len(level.roots) for level in search.levels] == [1, 1, 1, 0, 1, 0, 3]
    assert len(search.roots) == 16
    assert [level.height for level in search.levels] == sorted(
        level.height for level in search.levels
    )
    assert search.levels[0].height == D * (PHI - 1) ** 2 / 4


def test_archived_bounded_transcript_is_reproducible() -> None:
    archived = json.loads(TRANSCRIPT.read_text())
    assert archived == build_h11_vinberg_transcript(5)
    assert len(verify_h11_vinberg_transcript(archived).roots) == 16


def test_exact_search_can_resume_from_a_certified_prefix() -> None:
    prefix = run_h11_vinberg_search(5)
    full = run_h11_vinberg_search(6)
    extension = run_h11_vinberg_extension(
        prefix.roots[len(prefix.cone_roots) :],
        start_after=prefix.candidate_coordinates[-1],
        identity_bound=6,
    )
    assert extension.candidate_coordinates == full.candidate_coordinates[
        len(prefix.candidate_coordinates) :
    ]
    assert extension.levels == full.levels[len(prefix.levels) :]
    assert extension.new_roots == full.roots[len(prefix.roots) :]


def test_exact_hyperideal_extreme_ray_blocks_vinberg_stopping() -> None:
    search = run_h11_vinberg_search(5)
    witness = verify_hyperideal_extreme_ray(
        search.roots,
        h11_basis_gram(),
        BOUND5_HYPERIDEAL_WITNESS,
        dimension=9,
    )
    assert witness.norm == Golden(561, 0) / 8 - Golden(693, 0) * PHI / 16
    assert witness.norm.sign(0) > 0
    assert all(pairing.sign(0) <= 0 for pairing in witness.wall_pairings)


def test_archived_bound_twenty_prefix_has_an_exact_hyperideal_obstruction() -> None:
    archived = json.loads(BOUND20_TRANSCRIPT.read_text())
    roots = tuple(
        _decode_vector(root)
        for root in (
            *archived["height_zero"]["simple_roots"],
            *(root for level in archived["levels"] for root in level["roots"]),
        )
    )
    assert archived["control_coordinate_scan"]["identity_bound"] == 20
    assert [level["solution_count"] for level in archived["levels"]] == [
        1, 1, 1, 0, 1, 0, 3, 5, 1, 5, 0, 1, 5, 2, 1, 1, 9, 2, 4, 8,
        0, 0, 9, 0, 3, 11, 0, 6, 6,
    ]
    assert sum(sum(level["visited_by_depth"]) for level in archived["levels"]) == 65_621_346
    assert archived["output"]["wall_count"] == len(roots) == 95
    gram = h11_basis_gram()
    assert all(_pair(root, gram, root) == 2 for root in roots)
    assert all(
        _pair(left, gram, right).sign(0) <= 0
        for i, left in enumerate(roots)
        for right in roots[i + 1 :]
    )
    witness = verify_hyperideal_extreme_ray(
        roots, gram, BOUND20_HYPERIDEAL_WITNESS, dimension=9
    )
    assert witness.norm == Golden(81) / 28 - Golden(25) * PHI / 14
    assert all(pairing.sign(0) <= 0 for pairing in witness.wall_pairings)
