from fractions import Fraction
from coxeter_polytopes.exact.silver import SILVER_ONE,SILVER_ZERO,Silver
from coxeter_polytopes.search.silver_d10_e7_squared import *
def f(rank,i):return tuple(SILVER_ONE if j==i else SILVER_ZERO for j in range(rank))
def test_binary_glue_basis():
 assert basis_audit()
 assert glue_membership(f(10,8),f(7,5),(SILVER_ZERO,)*7)
 assert glue_membership(f(10,9),(SILVER_ZERO,)*7,f(7,5))
def test_packets_geometry():
 assert len(height_zero())==3 and all(r.norm==2 for r in height_zero())
 g=first_geometry();assert(g.rays,g.finite,g.cusps,g.hyperideal)==(704,4,0,700)
 p=run_prefix(15);non=[x for x in p.packets if x.roots];assert[len(x.roots) for x in non]==[3,16];assert sum(r.glue_nontrivial for r in non[1].roots)==7;assert sum(r.three_factor_glue for r in non[1].roots)==3
def test_witness():
 w=strict_witness();assert w.weight==Fraction(5,3);assert w.norm==Silver(Fraction(16,9),-1);assert max(w.pairings)==Silver(Fraction(-1,3),Fraction(-1,3));assert w.strict and w.hyperideal
