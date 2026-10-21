"""Missing-label and dependent-radical regressions for the cap graph audit."""
from fractions import Fraction as F
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import (
    Z,O,P,add,scale,bsign,angle_polynomials,finite_label)


def test_beyond_golden_angles_and_wrong_conjugate():
    polys=angle_polynomials()
    assert finite_label((scale(O,2),O,Z,Z),scale(O,2),Z,polys)==8
    assert finite_label((scale(O,2),scale(O,-1),Z,Z),scale(O,2),Z,polys) is None
    d=add(scale(O,9),scale(P,-3))
    assert finite_label((add(scale(O,2),scale(P,F(1,2))),scale(O,F(1,2)),Z,Z),d,Z,polys)==15


def test_dependent_radicals_and_zero_extension():
    assert bsign((Z,O,scale(O,-1),Z),scale(O,2),scale(O,2))==0
    assert bsign((O,scale(O,-999),Z,Z),Z,Z)==1
