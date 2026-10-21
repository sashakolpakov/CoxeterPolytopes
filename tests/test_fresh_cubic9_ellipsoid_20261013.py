"""Compare the trace ellipsoid against an independent complete box census."""
from pathlib import Path
import sys
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import fresh_cubic9_dihedral_independent_audit_20261013 as f
from fresh_cubic9_ellipsoid_audit_20261013 import ellipsoid,box_alphabet

@pytest.mark.parametrize('cap',[
    f.ONE,f.scale(f.ONE,4),f.mul(f.THETA,f.THETA),f.BETA,
    f.add(f.TWO,f.THETA),f.mul(f.add(f.ONE,f.THETA),f.add(f.ONE,f.THETA)),
])
def test_exact_ellipsoid_matches_unpruned_box(cap):
    assert set(ellipsoid(cap))==set(box_alphabet((cap,)*3))
