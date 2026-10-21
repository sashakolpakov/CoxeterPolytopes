"""Field arithmetic, complete small alphabets, and geometric certificate regressions."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from coxeter_polytopes.exact.cubic9 import Cubic9 as K,THETA9 as T,cubic9_determinant,cubic9_inertia
from fresh_cubic9_alphabet_20261013 import alphabet
import fresh_cubic9_dihedral_independent_audit_20261013 as f
from fresh_cubic9_enumeration_audit_20261013 import alphabet as independent_alphabet

def raw(x):return x.a,x.b,x.c

def test_field_against_generic_fraction_reduction():
    xs=[K(*v) for v in product(range(-1,2),repeat=3)]
    for x in xs:
        assert raw(x.conjugate(1))==raw(K(x.a-2*x.b+2*x.c,-x.b-x.c,x.b))
        assert x.conjugate(3)==x and x.trace()==3*x.a+6*x.c
        assert x.norm()==f.field_norm(raw(x))
        if x:assert raw(x.inverse())==f.inverse(raw(x))
        for e in range(3):assert x.sign(e)==f.sign(raw(x),e)
        for y in xs:assert raw(x*y)==f.mul(raw(x),raw(y))

@pytest.mark.parametrize('caps',[(K(1),)*3,(K(4),)*3,(K(30)/T,-K(1)/T,-K(1)/T)])
def test_thin_alphabet_matches_independent_unpruned_box(caps):
    assert {raw(x) for x in alphabet(caps)}==set(independent_alphabet(tuple(raw(x) for x in caps)))

def test_signature_and_known_compact_triangle():
    g=((-T,K(0),K(0)),(K(0),K(2),-T),(K(0),-T,K(2)))
    assert cubic9_determinant(g)==-T*(4-T*T)
    assert cubic9_inertia(g).negative==1
    assert all(cubic9_inertia(g,embedding=e).positive==3 for e in (1,2))
    path=Path(__file__).resolve().parents[1]/'results/fresh_cubic9_dihedral_triangle_h30_20261013.json'
    d=json.loads(path.read_text());dominant,_=f.binary_census();out=f.verify(d,dominant)
    assert out['verdict']=='EXACT_COMPACT_KNOWN_2_4_18_TRIANGLE'
    d['geometry']['positive_affine_normal_dependence'][0]['a']='999'
    with pytest.raises(AssertionError):f.verify(d,dominant)
