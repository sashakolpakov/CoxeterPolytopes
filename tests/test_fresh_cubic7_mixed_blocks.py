"""Regression gates for the changed A3/dihedral host."""
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import fresh_cubic7_mixed_blocks_independent_verify_20261013 as audit


def test_extra_a3_norm_four_reflection_is_primitive():
    f=audit.f
    r=(f.ZERO,*(f.ZERO for _ in range(6)),f.ONE,f.ZERO,f.neg(f.ONE))
    assert audit.norm(r)==f.FOUR
    assert all(f.integral(f.scale(x,f.F(1,2))) for x in audit.cov(r))
    assert not all(f.integral(f.scale(x,f.F(1,2))) for x in r)


def test_arithmetic_acute_two_pi_over_seven_is_rejected():
    f=audit.f
    a=(f.ZERO,f.ONE,f.ZERO,*(f.ZERO for _ in range(7)))
    y=(f.ONE,f.ZERO)
    found=False
    for _ in range(14):
        b=(f.ZERO,*y,*(f.ZERO for _ in range(7)))
        p=f.dot(a,audit.cov(b))
        if p==f.neg(f.THETA):
            assert audit.norm(a)==audit.norm(b)==f.TWO
            assert f.sign(p)<0 and f.mul(p,p)==f.mul(f.THETA,f.THETA)
            assert not audit.compatible(a,(b,))
            found=True
        y=(f.sub(f.mul(f.C,y[0]),y[1]),y[0])
    assert found


def test_norm_eight_gate_and_unique_integral_overclosure():
    assert audit.modular_norm8()['admissible_residue_count']==8
    assert audit.odd_overclosure()['integral_time_residue_classes_mod2']==[[0,0,0]]
