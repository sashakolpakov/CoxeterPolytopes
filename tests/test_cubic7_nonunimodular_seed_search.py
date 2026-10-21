from collections import Counter
from pathlib import Path
from fractions import Fraction

from coxeter_polytopes.exact.cubic7 import THETA7, Cubic7
from coxeter_polytopes.search.cubic7_nonunimodular import (
    DELTA7,
    L2_SEED,
    L7_SEED,
    P13_R10_SEED,
    P13_R7_SEED,
    P13_R8_SEED,
    P13_SEEDS,
    enumerate_cone_rays_incremental,
    enumerate_height_zero_mirrors,
    enumerate_nonunimodular_root_pool,
    height_zero_is_b3_times_a1,
    nonunimodular_diagonal_pair,
    run_nonunimodular_vinberg_search,
)
from coxeter_polytopes.verify.cubic7_alvin import parse_alvin_roots_text


_ALVIN_FIXTURES = Path(__file__).with_name("fixtures") / "cubic7_alvin"


def test_local_norm_manifests_are_totally_positive_and_squarefree():
    assert DELTA7 == 2 - THETA7
    assert DELTA7.norm() == 7
    assert all(DELTA7.sign(embedding) > 0 for embedding in range(3))
    assert L2_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert L7_SEED.canonical_mirror_norms == (
        Cubic7(1),
        Cubic7(2),
        DELTA7,
        2 * DELTA7,
    )


def test_height_zero_mirror_exhaustion_is_b3_times_a1_for_both_seeds():
    for seed in (L2_SEED, L7_SEED):
        roots = enumerate_height_zero_mirrors(seed)
        assert len(roots) == 20
        assert height_zero_is_b3_times_a1(seed)
        assert all(root[0] == 0 for root in roots)

    assert Counter(
        nonunimodular_diagonal_pair(L2_SEED, root, root)
        for root in enumerate_height_zero_mirrors(L2_SEED)
    ) == Counter({Cubic7(2): 14, Cubic7(1): 6})
    assert Counter(
        nonunimodular_diagonal_pair(L7_SEED, root, root)
        for root in enumerate_height_zero_mirrors(L7_SEED)
    ) == Counter({Cubic7(2): 12, Cubic7(1): 6, DELTA7: 2})


def test_rank_five_seed_prefixes_are_exhaustive_through_height_twenty():
    expected = {
        "L2": {
            "controls": ((Cubic7(1), 1), (Cubic7(2), 7)),
            "pool": Counter({Cubic7(2): 10, Cubic7(1): 1}),
            "accepted": Counter({Cubic7(2): 4}),
            "walls": 8,
        },
        "L7": {
            "controls": (
                (Cubic7(1), 1),
                (Cubic7(2), 7),
                (DELTA7, 8),
                (2 * DELTA7, 20),
            ),
            "pool": Counter(
                {Cubic7(2): 11, 2 * DELTA7: 2, Cubic7(1): 1, DELTA7: 1}
            ),
            "accepted": Counter({Cubic7(2): 3, 2 * DELTA7: 1}),
            "walls": 8,
        },
    }

    for seed in (L2_SEED, L7_SEED):
        search = run_nonunimodular_vinberg_search(seed, 20)
        datum = expected[seed.seed_id]
        assert search.pool.control_counts_by_norm == datum["controls"]
        assert Counter(root.norm for root in search.pool.roots) == datum["pool"]
        assert Counter(root.norm for root in search.accepted_roots) == datum["accepted"]
        assert len(search.walls) == datum["walls"]
        assert all(
            nonunimodular_diagonal_pair(seed, left, right).sign(0) <= 0
            for i, left in enumerate(search.walls)
            for right in search.walls[i + 1 :]
        )
        for root in search.pool.roots:
            assert nonunimodular_diagonal_pair(seed, root.vector, root.vector) == root.norm
            assert (Cubic7(20) - root.height).sign(0) >= 0
            assert root.control.sign(0) > 0


def test_l7_spatial_r1_wall_is_enumerated_and_selected_exactly():
    """Guard the nonintegral R1 branch in AlVin's recursive coordinate box.

    AlVin sorts the diagonal so that ``2-theta`` is its first spatial
    coefficient and prints this wall as
    ``([0,1,-1], [0,-1,-1], 0, 0, 0)``.  In this backend the exceptional
    coefficient is kept last, so the equivalent vector has ``1+theta`` in
    its last coordinate.
    """

    target = (
        Cubic7(-3, 1, 2),  # [0,1,-1] = -3 + theta + 2 theta^2
        Cubic7(0),
        Cubic7(0),
        Cubic7(0),
        Cubic7(1, 1),  # [0,-1,-1] = 1 + theta
    )
    expected_norm = 2 * DELTA7
    expected_height = Cubic7(Fraction(-1, 2), 1, Fraction(1, 2))

    assert nonunimodular_diagonal_pair(L7_SEED, target, target) == expected_norm
    pool = enumerate_nonunimodular_root_pool(L7_SEED, 2)
    root = next(candidate for candidate in pool.roots if candidate.vector == target)
    assert root.norm == expected_norm
    assert root.height == expected_height

    search = run_nonunimodular_vinberg_search(L7_SEED, 2)
    assert target in tuple(candidate.vector for candidate in search.accepted_roots)


def test_incremental_exact_ray_census_matches_archived_brute_counts():
    expected = {"L2": (9, 0, 6), "L7": (9, 0, 7)}
    for seed in (L2_SEED, L7_SEED):
        census = enumerate_cone_rays_incremental(
            run_nonunimodular_vinberg_search(seed, 20)
        )
        assert (
            census.finite_count,
            census.cusp_count,
            census.hyperideal_count,
        ) == expected[seed.seed_id]


def test_same_height_vinberg_packets_are_selected_simultaneously():
    # Both seeds have a genuine multiplicity-two accepted height by H=100.
    # This guards against an arbitrary secondary vector order choosing one
    # representative and suppressing another at the same Vinberg height.
    for seed in (L2_SEED, L7_SEED):
        search = run_nonunimodular_vinberg_search(seed, 100)
        multiplicities = Counter(root.height for root in search.accepted_roots)
        assert max(multiplicities.values()) >= 2
        for height, multiplicity in multiplicities.items():
            if multiplicity < 2:
                continue
            packet = tuple(root for root in search.accepted_roots if root.height == height)
            assert all(
                nonunimodular_diagonal_pair(seed, left.vector, right.vector).sign(0)
                <= 0
                for index, left in enumerate(packet)
                for right in packet[index + 1 :]
            )


def test_v3_alvin_inert2_ten_wall_diagnostic_matches_exact_l2_prefix():
    """Freeze the bounded live diagnostic without claiming termination.

    AlVin was deliberately run with ``-maxv 10`` on the spatial-R1-fixed
    image.  The exact backend exhausts through height 100, well beyond these
    ten walls.  Walls 6 and 7 have the same exact Vinberg height, so their
    harmless secondary-order swap is compared as an unordered packet.
    """

    transcript = (_ALVIN_FIXTURES / "inert2_v3_maxv10.alvin.txt").read_text()
    assert "Algorithm ended" not in transcript
    assert "The algorithm did not terminate" in transcript
    assert "Possible values for (e,e): 1, 2, 4" in transcript

    live = parse_alvin_roots_text(transcript)
    exact = run_nonunimodular_vinberg_search(L2_SEED, 100).walls[:10]
    assert len(live) == len(exact) == 10
    assert live[:5] == exact[:5]
    assert set(live[5:7]) == set(exact[5:7])
    assert live[7:] == exact[7:]

    live_norms = tuple(nonunimodular_diagonal_pair(L2_SEED, root, root) for root in live)
    assert all(norm in L2_SEED.canonical_mirror_norms for norm in live_norms)
    assert Cubic7(4) not in live_norms
    assert L2_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))


def test_v3_alvin_p13_r7_ten_wall_diagnostic_matches_exact_prefix():
    """Freeze AlVin's bounded p13-r7 prefix, including its height tie.

    The tenth printed root is the third member of an exact three-wall
    same-height packet.  ``-maxv 10`` truncates inside that packet, so this
    fixture checks enumeration agreement only and makes no termination or
    finite-volume claim.
    """

    transcript = (_ALVIN_FIXTURES / "p13_r7_v3_maxv10.alvin.txt").read_text()
    assert "Algorithm ended" not in transcript
    assert "The algorithm did not terminate" in transcript
    assert "Unknown weight" not in transcript

    live = parse_alvin_roots_text(transcript)
    search = run_nonunimodular_vinberg_search(P13_R7_SEED, 300)
    exact = search.walls
    assert len(live) == 10
    assert live[:9] == exact[:9]
    assert live[9] == exact[11]

    metadata = {root.vector: root for root in search.accepted_roots}
    packet = exact[9:12]
    assert len(packet) == 3
    assert len({metadata[root].height for root in packet}) == 1
    assert all(
        nonunimodular_diagonal_pair(P13_R7_SEED, left, right).sign(0) <= 0
        for index, left in enumerate(packet)
        for right in packet[index + 1 :]
    )
    assert live[9] in packet


def test_p13_seed_norm_classes_and_height_zero_systems_are_exact():
    expected_coefficients = (
        Cubic7(0, -1, 2),
        Cubic7(1, -1, 1),
        Cubic7(1, 0, 1),
    )
    assert P13_SEEDS == (P13_R7_SEED, P13_R10_SEED, P13_R8_SEED)
    for seed, coefficient in zip(P13_SEEDS, expected_coefficients, strict=True):
        assert seed.last_coefficient == coefficient
        assert coefficient.norm() == 13
        assert all(coefficient.sign(embedding) > 0 for embedding in range(3))
        assert seed.canonical_mirror_norms == (
            Cubic7(1),
            Cubic7(2),
            coefficient,
            2 * coefficient,
        )
        roots = enumerate_height_zero_mirrors(seed)
        assert len(roots) == 20
        assert height_zero_is_b3_times_a1(seed)
        assert Counter(
            nonunimodular_diagonal_pair(seed, root, root) for root in roots
        ) == Counter({Cubic7(2): 12, Cubic7(1): 6, coefficient: 2})


def test_p13_prefixes_are_exhaustive_and_open_through_height_twenty():
    expected = {
        "P13-r7": (12, 3, 7, (7, 0, 5)),
        "P13-r10": (10, 3, 7, (7, 0, 5)),
        "P13-r8": (11, 3, 7, (7, 0, 5)),
    }
    for seed in P13_SEEDS:
        search = run_nonunimodular_vinberg_search(seed, 20)
        pool_count, accepted_count, wall_count, ray_counts = expected[seed.seed_id]
        assert len(search.pool.roots) == pool_count
        assert len(search.accepted_roots) == accepted_count
        assert len(search.walls) == wall_count
        assert Counter(root.norm for root in search.accepted_roots) == Counter(
            {Cubic7(2): 3}
        )
        census = enumerate_cone_rays_incremental(search)
        assert (
            census.finite_count,
            census.cusp_count,
            census.hyperideal_count,
        ) == ray_counts
        assert all(
            nonunimodular_diagonal_pair(seed, left, right).sign(0) <= 0
            for index, left in enumerate(search.walls)
            for right in search.walls[index + 1 :]
        )
