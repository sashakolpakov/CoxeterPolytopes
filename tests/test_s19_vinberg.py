from coxeter_polytopes.search.s19_vinberg import (
    dominant_a2_vectors,
    dominant_a2_divisible_three_vectors,
    dominant_e8_vectors,
    enumerate_level_candidates,
    height_zero_roots,
    height_zero_simple_roots,
    pair,
    root_divisibility,
    build_s19_vinberg_transcript,
    run_s19_vinberg_search,
)


def test_height_zero_system_is_exact_a1_2e8_g2() -> None:
    all_roots = height_zero_roots()
    simple = height_zero_simple_roots()
    assert len(all_roots) == 494
    assert len(simple) == 19
    assert sum(pair(root, root) == -2 for root in all_roots) == 488
    assert sum(pair(root, root) == -6 for root in all_roots) == 6
    assert sum(pair(root, root) == -2 for root in simple) == 18
    assert sum(pair(root, root) == -6 for root in simple) == 1
    assert {(-pair(root, root), root_divisibility(root)) for root in all_roots} == {
        (2, 1),
        (6, 3),
    }


def test_dominant_ade_enumerators_reconstruct_low_shells() -> None:
    e8 = dominant_e8_vectors(6)
    a2 = dominant_a2_vectors(6)
    assert [(vector.square, vector.weights) for vector in e8] == [
        (0, (0, 0, 0, 0, 0, 0, 0, 0)),
        (2, (0, 0, 0, 0, 0, 0, 1, 0)),
        (4, (1, 0, 0, 0, 0, 0, 0, 0)),
        (6, (0, 0, 0, 0, 0, 1, 0, 0)),
    ]
    assert [(vector.square, vector.weights) for vector in a2] == [
        (0, (0, 0)),
        (2, (1, 1)),
        (6, (0, 3)),
    ]
    assert [(vector.square, vector.weights) for vector in dominant_a2_divisible_three_vectors(24)] == [
        (0, (0, 0)),
        (6, (0, 3)),
        (18, (3, 3)),
        (24, (0, 6)),
    ]


def test_first_mixed_vinberg_levels_are_exhaustive_and_acute() -> None:
    search = run_s19_vinberg_search(3)
    assert [(level.control, level.norm_abs) for level in search.levels] == [
        (1, 2),
        (3, 6),
        (2, 2),
        (3, 2),
    ]
    assert search.levels[0].candidate_count == 3
    assert len(search.levels[0].roots) == 3
    assert enumerate_level_candidates(1, 6) == ()
    assert len(search.roots) >= 22
    assert all(
        pair(left, right) >= 0
        for i, left in enumerate(search.roots)
        for right in search.roots[i + 1 :]
    )


def test_prefix_through_control_six_reaches_known_norm_two_wall_count() -> None:
    search = run_s19_vinberg_search(6)
    assert [
        (level.control, level.norm_abs, level.candidate_count, len(level.roots))
        for level in search.levels
    ] == [
        (1, 2, 3, 3),
        (3, 6, 1, 0),
        (2, 2, 8, 0),
        (3, 2, 11, 0),
        (6, 6, 4, 0),
        (4, 2, 37, 1),
        (5, 2, 56, 0),
        (9, 6, 7, 0),
        (6, 2, 148, 2),
    ]
    assert len(search.roots) == 25
    assert sum(pair(root, root) == -2 for root in search.roots) == 24
    assert sum(pair(root, root) == -6 for root in search.roots) == 1

    # The very next norm-six level contributes two further walls.  Hence the
    # height-18 prefix is genuinely incomplete and Vinberg stopping cannot
    # have fired there.
    next_long = tuple(
        candidate.root
        for candidate in enumerate_level_candidates(12, 6)
        if all(pair(candidate.root, wall) >= 0 for wall in search.roots)
    )
    assert len(next_long) == 2

    transcript = build_s19_vinberg_transcript(6)
    assert transcript["output"]["wall_count_by_abs_norm"] == {"2": 24, "6": 1}
    assert transcript["stopping"]["attempted"] is False
    assert transcript["stopping"]["fired"] is False
