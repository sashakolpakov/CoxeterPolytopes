from fractions import Fraction
from coxeter_polytopes.exact.silver import SILVER_ONE,SILVER_ZERO,Silver
from coxeter_polytopes.search.silver_a11_d7_e6 import *

def fundamental(rank,index):return tuple(SILVER_ONE if i==index else SILVER_ZERO for i in range(rank))
def test_cyclic_glue_and_basis():
    assert basis_audit()
    assert glue_membership(fundamental(11,0),fundamental(7,5),fundamental(6,0))
    assert a11_residue(fundamental(11,0))==(1,0)
    assert d7_residue(fundamental(7,5))==(1,0)
    assert e6_residue(fundamental(6,0))==(1,0)
def test_height_zero_geometry_and_packets():
    assert len(height_zero())==3 and all(r.norm==2 for r in height_zero())
    g=first_geometry();assert (g.rays,g.finite,g.cusps,g.hyperideal)==(672,3,0,669)
    p=run_prefix(15);nonempty=[x for x in p.packets if x.roots]
    assert [len(x.roots) for x in nonempty]==[3,35]
    assert sum(r.glue_nontrivial for r in nonempty[1].roots)==22
    assert sum(r.three_factor_glue for r in nonempty[1].roots)==12
def test_strict_witness():
    w=strict_witness();assert w.weight==Fraction(5,3)
    assert w.norm==Silver(Fraction(167,108),-1)
    assert max(w.pairings)==Silver(Fraction(-1,3),Fraction(-1,3))
    assert w.strict and w.hyperideal and len(w.pairings)==38
