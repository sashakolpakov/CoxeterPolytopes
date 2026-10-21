from fractions import Fraction

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, Silver
from coxeter_polytopes.search.niemeier_host_triage import (
    cartan_gram,
    highest_root_marks,
)
from coxeter_polytopes.search.silver_d12_squared import (
    D12_COMPONENT,
    SQRT2_D12_SQUARED_SEED,
    TWO_SQRT2_D12_SQUARED_SEED,
    UNIT_D12_SQUARED_SEED,
    d12_squared_basis_audit,
    d12_squared_basis_glue_audit,
    d12_squared_glue_membership,
    first_packet_geometry,
    height_zero_seed_dominant_roots,
    run_d12_squared_seed_prefix,
    strict_hyperideal_witness,
)


def _residue_weight(index: int):
    return tuple(SILVER_ONE if position == index else SILVER_ZERO for position in range(12))


def test_vector_spinor_graph_glue_and_unimodular_basis_are_exact():
    assert d12_squared_basis_glue_audit()
    assert d12_squared_basis_audit()
    assert highest_root_marks(D12_COMPONENT) == (1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1)
    assert cartan_gram(D12_COMPONENT)[9][10] == -1
    assert cartan_gram(D12_COMPONENT)[9][11] == -1
    vector = _residue_weight(0)
    spinor = _residue_weight(10)
    cospinor = _residue_weight(11)
    assert d12_squared_glue_membership(vector, spinor)
    assert d12_squared_glue_membership(spinor, vector)
    assert d12_squared_glue_membership(cospinor, cospinor)
    assert not d12_squared_glue_membership(vector, vector)


def test_all_three_height_zero_censuses_are_exactly_d12_squared():
    for seed in (
        UNIT_D12_SQUARED_SEED,
        SQRT2_D12_SQUARED_SEED,
        TWO_SQRT2_D12_SQUARED_SEED,
    ):
        roots = height_zero_seed_dominant_roots(seed)
        assert len(roots) == 2
        assert all(root.norm == 2 and root.control == 0 for root in roots)
        assert all(not root.glue_nontrivial for root in roots)


def test_unit_seed_three_nonempty_packets_and_deep_witness():
    prefix = run_d12_squared_seed_prefix(UNIT_D12_SQUARED_SEED, 30)
    nonempty = [packet for packet in prefix.packets if packet.roots]
    assert [len(packet.roots) for packet in nonempty] == [2, 8, 89]
    assert [packet.rejected_count for packet in nonempty] == [0, 0, 4]
    assert [sum(root.glue_nontrivial for root in packet.roots) for packet in nonempty] == [0, 4, 33]
    witness = strict_hyperideal_witness(UNIT_D12_SQUARED_SEED)
    assert witness.weight == Fraction(8, 5)
    assert witness.norm == Silver(Fraction(39, 25), -1)
    assert witness.all_positive_walls_strict and witness.hyperideal
    assert len(witness.positive_wall_pairings) == 99


def test_sqrt2_seed_two_nonempty_packets_and_strict_witness():
    prefix = run_d12_squared_seed_prefix(SQRT2_D12_SQUARED_SEED, 40)
    nonempty = [packet for packet in prefix.packets if packet.roots]
    assert [len(packet.roots) for packet in nonempty] == [2, 17]
    assert [sum(root.glue_nontrivial for root in packet.roots) for packet in nonempty] == [0, 5]
    witness = strict_hyperideal_witness(SQRT2_D12_SQUARED_SEED)
    assert (witness.dynkin_index, witness.weight) == (2, Fraction(7, 10))
    assert witness.norm == Silver(Fraction(147, 100), -1)
    assert max(witness.positive_wall_pairings) == Silver(Fraction(-3, 5), Fraction(2, 5))
    assert witness.all_positive_walls_strict and witness.hyperideal


def test_two_sqrt2_seed_is_wall_rich_but_still_hyperideal():
    prefix = run_d12_squared_seed_prefix(TWO_SQRT2_D12_SQUARED_SEED, 27)
    nonempty = [packet for packet in prefix.packets if packet.roots]
    assert [len(packet.roots) for packet in nonempty] == [4, 8, 15]
    assert [packet.rejected_count for packet in nonempty] == [0, 0, 2]
    assert [sum(root.glue_nontrivial for root in packet.roots) for packet in nonempty] == [0, 2, 5]
    witness = strict_hyperideal_witness(TWO_SQRT2_D12_SQUARED_SEED)
    assert (witness.dynkin_index, witness.weight) == (4, Fraction(7, 9))
    assert witness.norm == Silver(Fraction(245, 81), -2)
    assert max(witness.positive_wall_pairings) == Silver(Fraction(-1, 9), Fraction(-1, 9))
    assert witness.all_positive_walls_strict and witness.hyperideal


def test_unit_affine_product_geometry_has_hyperideal_rays():
    geometry = first_packet_geometry()
    assert (
        geometry.extreme_ray_count,
        geometry.finite_count,
        geometry.cusp_count,
        geometry.hyperideal_count,
    ) == (169, 3, 0, 166)
