"""Adversarial arithmetic and archived-witness checks for the new H9 scout."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import itertools

ROOT=Path(__file__).resolve().parents[1]
def load(name):
    spec=spec_from_file_location(name,ROOT/'scripts'/f'{name}.py')
    module=module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module
S=load('search_quadratic_diagonal_h9_2026_10_13')
V=load('verify_quadratic_diagonal_h9_2026_10_13')

def test_field_arithmetic_against_existing_independent_classes():
    from coxeter_polytopes.exact.golden import Golden
    from coxeter_polytopes.exact.silver import Silver
    for d,T in [(2,Silver),(5,Golden)]:
        K=S.Field(d)
        for x in itertools.product(range(-6,7),repeat=2):
            v=T(*x)
            for p in (0,1):assert K.sign(x,p)==v.sign(p)
            y=(x[1]-1,x[0]+2)
            expected=v*T(*y)
            assert K.mul(x,y)==(expected.a,expected.b)

def test_low_dimensional_search_count_against_cartesian_enumeration():
    for d,alpha in [(2,(1,1)),(3,(0,1)),(5,(0,1))]:
        K=S.Field(d);n=2;H=3;count=0
        for q in (1,2):
            # Separate rectangular brute force is larger than these controls.
            controls=[]
            for x in itertools.product(range(-12,13),repeat=2):
                ax2=K.mul(alpha,K.sq(x))
                if K.sign(x)>0 and K.sign(S.sub((q*H,0),ax2))>=0 and K.sign(S.add((q,0),ax2),1)>=0:controls.append(x)
            assert set(controls)==set(K.controls(alpha,q,H))
            for x in controls:
                target=S.add((q,0),K.mul(alpha,K.sq(x)))
                alphabet=K.spatial(target)
                for y,z in itertools.product(alphabet,repeat=2):
                    if K.sign(S.sub(y,z))>=0 and S.add(K.sq(y),K.sq(z))==target:count+=1
        assert S.search(d,alpha,n,H)['dominant_root_count']==count

def test_all_frozen_prefixes_have_exact_spacelike_witnesses():
    assert V.verify_all()==132
