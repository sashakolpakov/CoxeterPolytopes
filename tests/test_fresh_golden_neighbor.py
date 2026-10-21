from fractions import Fraction
import json
from pathlib import Path
import pytest
from coxeter_polytopes.exact.golden import Golden as G,PHI
from scripts.fresh_golden_neighbor_verify import verify

ROOT=Path(__file__).resolve().parents[1]

def form(matrix,left,right):
    return sum((left[i]*matrix[i][j]*right[j] for i in range(2) for j in range(2)),G(0))

def test_unique_glue_diagonalizes_to_unimodular_binary_block():
    gram=((-2*PHI,G(0)),(G(0),G(2)))
    w=(-PHI/2,G(Fraction(-1,2)))
    f=((PHI-1)/2,PHI/2)
    assert form(gram,w,w)==-PHI
    assert form(gram,f,f)==1
    assert form(gram,w,f)==0
    assert w[0]*f[1]-w[1]*f[0]==G(Fraction(-1,2))
    # Reduction of the integral-overlattice condition modulo the inert2.
    values=[G(a,b) for a in range(2) for b in range(2)]
    isotropic=[(a,b) for a in values for b in values if a or b
        if all(v.denominator==1 and int(v)%2==0 for v in ((b*b-PHI*a*a).a,(b*b-PHI*a*a).b))]
    assert len(isotropic)==3 # one projective F4 line
    assert all(a for a,b in isotropic)

def test_second_binary_variant_is_only_an_integral_basis_change():
    gram=((-PHI,G(0)),(G(0),G(2)))
    w=(2*PHI,PHI*PHI); f=(PHI,PHI)
    assert form(gram,w,w)==-2*PHI
    assert form(gram,f,f)==1 and form(gram,w,f)==0
    assert w[0]*f[1]-w[1]*f[0]==1

def test_generic_generator_calibrates_against_prior_engine():
    a=json.loads((ROOT/'results/fresh_golden_neighbor_calibration.json').read_text())
    b=json.loads((ROOT/'results/golden_h9_odd_a1_full_norms_h30_20261013.json').read_text())
    assert a['completed'] and a['wall_count']==17
    assert a['all_walls']==b['all_walls']

@pytest.mark.parametrize('name',['sqrt5_h30','sqrt5_h100','3phi_h30','1p4phi_h30','m1p3phi_h30'])
def test_standalone_escape_replays(name):
    verify(ROOT/f'results/fresh_golden_neighbor_{name}.json')
