from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_seed_manifest import (
    alvin_rc7_to_cubic7,
    cubic7_to_alvin_rc7,
    initial_cubic7_diagonal_seeds,
)


def test_cubic7_alvin_basis_conversion_is_exact():
    values = (
        Cubic7(1),
        Cubic7(0, 1),
        Cubic7(0, 0, 1),
        Cubic7(3, -4, 5),
    )
    assert cubic7_to_alvin_rc7(Cubic7(1)) == (-1, -1, -1)
    assert cubic7_to_alvin_rc7(Cubic7(0, 1)) == (1, 0, 0)
    for value in values:
        assert alvin_rc7_to_cubic7(cubic7_to_alvin_rc7(value)) == value


def test_initial_cubic7_diagonal_seed_box_is_exact_and_nonunimodular():
    seeds = initial_cubic7_diagonal_seeds()
    assert len(seeds) == 13
    assert len({seed.seed_id for seed in seeds}) == len(seeds)
    assert [seed.radical_norm for seed in seeds] == [
        1,
        7,
        8,
        13,
        13,
        13,
        56,
        91,
        91,
        91,
        104,
        104,
        104,
    ]
    assert all(seed.field_norm == seed.radical_norm for seed in seeds)
    assert all(
        seed.coefficient.sign(embedding) > 0
        for seed in seeds
        for embedding in range(3)
    )
    assert all(seed.field_norm != 1 for seed in seeds if not seed.calibration)

    # The three norm-13 elements generate the three different primes above
    # 13: theta has residues 7, 10 and 8, respectively.
    split = seeds[3:6]
    residues = []
    for seed in split:
        value = seed.coefficient
        roots = tuple(
            root
            for root in (7, 8, 10)
            if (
                int(value.a)
                + int(value.b) * root
                + int(value.c) * root * root
            )
            % 13
            == 0
        )
        residues.append(roots)
    assert residues == [(7,), (10,), (8,)]


def test_manifest_emits_rank_five_and_six_alvin_forms():
    seeds = initial_cubic7_diagonal_seeds()
    assert seeds[0].alvin_qf(5) == "[-1,0,0],1,1,1,[-1,-1,-1]"
    assert seeds[1].alvin_qf(5) == "[-1,0,0],1,1,1,[-3,-2,-2]"
    assert seeds[2].alvin_qf(6) == "[-1,0,0],1,1,1,1,[-2,-2,-2]"
    for seed in seeds:
        assert len(seed.alvin_qf(5).split(",")) >= 5
        assert "-k rc7" in seed.alvin_command(5)
