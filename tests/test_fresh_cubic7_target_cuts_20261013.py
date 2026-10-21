"""Independent enumeration comparisons and a forbidden acute-angle regression."""
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from fresh_cubic7_dihedral_h9_20261013 import alphabet,K,T,O,Z,C,BETA,Host,mv,pair
from fresh_cubic7_target_alphabet_20261013 import thin_alphabet
from fresh_cubic7_asymmetric_cuts_20261013 import compatible

@pytest.mark.parametrize('caps',[(K(31),O,O),(K(100)/T,-O/T,-O/T),
                                 (2*BETA+200,2*BETA,2*BETA),(K(4),K(4),K(4))])
def test_thin_enumeration_matches_complete_box(caps):
    assert thin_alphabet(caps)==alphabet(caps)

def test_integral_acute_mirror_pair_can_have_forbidden_angle():
    host=Host(2*T,4,True);g=host.gram
    old=host.initial()[0]
    bad=(Z,1+T,C,*([Z]*7))
    assert pair(bad,g,bad)==2
    assert pair(bad,g,old)==-T
    assert not compatible(bad,K(2),(mv(g,old),),(K(2),))
    good=(Z,C,O,*([Z]*7))
    assert pair(good,g,good)==2
    assert compatible(good,K(2),(mv(g,old),),(K(2),))

def test_projection_obstruction_certificate():
    from fresh_cubic7_block_symmetry_obstruction_20261013 import certificate,verify
    d=certificate();verify(d)
    d['positive_coefficients'][0]['a']='1'
    with pytest.raises(AssertionError):verify(d)
