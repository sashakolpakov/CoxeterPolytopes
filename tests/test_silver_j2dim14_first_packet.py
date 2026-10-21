from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_j2dim14 import (
    FIVE_P7_SEED,
    FIVE_SQRT2_SEED,
    FIVE_TWO_SQRT2_SEED,
    FIVE_UNIT_SEED,
    J2DIM14_NORM_TEN_SHELL_COUNT,
    P7_SEED,
    SQRT2_SEED,
    TWO_SQRT2_SEED,
    UNIT_SEED,
    enumerate_controls,
    gated_gram_triples,
    height_zero_gated_gram_triples,
    j2dim14_catalogue_audit,
    norm_classes,
    unit_first_packet_audit,
)


def test_catalogue_and_height_zero_local_gates():
    assert j2dim14_catalogue_audit()
    seeds = (
        UNIT_SEED, SQRT2_SEED, TWO_SQRT2_SEED, P7_SEED,
        FIVE_UNIT_SEED, FIVE_SQRT2_SEED, FIVE_TWO_SQRT2_SEED, FIVE_P7_SEED,
    )
    assert [len(norm_classes(seed)) for seed in seeds] == [4, 4, 4, 8, 4, 4, 4, 8]
    assert all(
        not height_zero_gated_gram_triples(norm_class)
        for seed in seeds
        for norm_class in norm_classes(seed)
    )


def test_unit_first_packet_and_hyperideal_witness():
    audit = unit_first_packet_audit()
    assert audit.control_count == 12
    assert audit.distinct_height_count == 8
    assert audit.earlier_control_count == 10
    assert audit.wall_count == J2DIM14_NORM_TEN_SHELL_COUNT == 374_976
    assert audit.dual_square.numerator == 8 and audit.dual_square.denominator == 5
    assert audit.maximum_shell_pairing == 2
    assert audit.witness_lorentz_norm == Silver(13, -5) / 5
    assert audit.wall_margin == -1 + 2 * SQRT2
    assert audit.strict_hyperideal


def test_other_first_necessary_heights():
    expected = {
        SQRT2_SEED: 96 + 68 * SQRT2,
        TWO_SQRT2_SEED: 96 + 68 * SQRT2,
        P7_SEED: 124 + 88 * SQRT2,
    }
    for seed, height in expected.items():
        live = [
            row for row in enumerate_controls(seed, Silver(300))
            if gated_gram_triples(row)
        ]
        assert live and live[0].height == height
