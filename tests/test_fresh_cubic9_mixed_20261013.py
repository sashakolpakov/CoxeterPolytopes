"""Focused gates for the changed conductor-nine mixed-block host."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import fresh_cubic9_mixed_search_20261013 as m

def test_primitive_norm_four_is_not_discarded():
    roots,_=m.candidates(m.T*m.T-2,m.K(4))
    assert len(roots)==1
    r=roots[0]
    assert m.pair(r,m.G,r)==4
    assert all((v/2).is_integral for v in m.mv(m.G,r))
    assert not all((v/2).is_integral for v in r)

def test_full_zero_census_includes_all_safe_classes():
    counts=[len(m.candidates(m.Z,q)[0]) for q in m.NORMS]
    assert counts==[4,0,0,0,2,0,0,0]
    assert len(m.initial())==9

def test_wrong_cubic_conjugate_angle_is_rejected():
    wrong=4+m.T-m.T*m.T  # 4cos^2(2pi/9), not pi/m for integral m.
    assert wrong.sign()>0 and (wrong-4).sign()<0
    assert wrong not in m.ALLOWED
    assert 2+m.T in m.ALLOWED  # pi/18.
