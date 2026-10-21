from fractions import Fraction

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, Silver
from coxeter_polytopes.search.silver_a12_squared import (
    a12_residue,
    a12_squared_basis_audit,
    a12_squared_glue_membership,
    first_packet_geometry,
    height_zero_dominant_roots,
    run_prefix,
    strict_hyperideal_witness,
)


def _fundamental(index):
    return tuple(SILVER_ONE if i == index else SILVER_ZERO for i in range(12))


def test_prime13_graph_and_basis():
    assert a12_squared_basis_audit()
    assert a12_residue(_fundamental(4)) == (5, 0)
    assert a12_residue(_fundamental(0)) == (1, 0)
    assert a12_squared_glue_membership(_fundamental(4), _fundamental(0))
    assert not a12_squared_glue_membership(_fundamental(0), _fundamental(0))


def test_height_zero_and_affine_geometry():
    roots = height_zero_dominant_roots()
    assert len(roots) == 2 and all(root.norm == 2 for root in roots)
    geometry = first_packet_geometry()
    assert (geometry.extreme_ray_count, geometry.finite_count, geometry.cusp_count, geometry.hyperideal_count) == (169,1,0,168)


def test_deep_prefix_and_strict_witness():
    prefix = run_prefix(150,None)
    nonempty = [packet for packet in prefix.packets if packet.roots]
    assert [len(packet.roots) for packet in nonempty] == [2,4,20,30,164,28,32]
    assert [packet.rejected_count for packet in nonempty] == [0,0,4,10,16,20,12]
    assert len(prefix.accepted_roots) == 280
    assert all(not root.glue_nontrivial for root in prefix.accepted_roots)
    witness = strict_hyperideal_witness()
    assert witness.weight == Fraction(5,3)
    assert witness.norm == Silver(Fraction(61,39),-1)
    assert max(witness.pairings) == Silver(0,Fraction(-1,3))
    assert witness.all_positive_walls_strict and witness.hyperideal
