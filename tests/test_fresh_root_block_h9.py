"""Structural and independent-certificate checks for the new E7+A2 lane."""
import json
from pathlib import Path
import pytest
from coxeter_polytopes.exact.golden import Golden as G, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.search.bugaenko_vinberg import _pair, _fraction_inverse
from scripts.fresh_root_block_h9_search import (
    host_gram, norm_candidates, spatial_gram, zero_roots,
)
from scripts.fresh_root_block_h9_verify import verify

ROOT=Path(__file__).resolve().parents[1]

def test_root_block_and_discriminant_exponent():
    spatial=spatial_gram()
    assert determinant(tuple(tuple(G(x) for x in row) for row in spatial))==6
    assert all((6*x).denominator==1 for row in _fraction_inverse(spatial) for x in row)
    for alpha in (PHI,2*PHI):
        gram=host_gram(alpha)
        assert determinant(gram)==-6*alpha
        assert (inertia(gram,embedding=0).positive,inertia(gram,embedding=0).negative)==(9,1)
        assert inertia(gram,embedding=1).positive==10
    assert norm_candidates(PHI)==tuple(map(G,(1,2,3,4,6,12)))
    assert norm_candidates(2*PHI)==tuple(map(G,(2,4,6,12)))

def test_full_zero_system_includes_the_a2_diagram_reflection():
    gram=host_gram(PHI); roots=zero_roots()
    assert len(roots)==9
    assert [_pair(r,gram,r) for r in roots]==[G(2)]*8+[G(6)]
    assert _pair(roots[-2],gram,roots[-1])==-3
    assert all(_pair(r,gram,s).sign()<=0 for i,r in enumerate(roots) for s in roots[:i])
    omitted_a2_simple=tuple(G(-int(j==9)) for j in range(10))
    assert _pair(omitted_a2_simple,gram,roots[-1])==3
    for r in roots:
        q=_pair(r,gram,r)
        for i in range(10):
            e=tuple(G(int(j==i)) for j in range(10))
            assert (2*_pair(r,gram,e)/q).is_integral

def test_spatial_norm_two_root_closure_has_e7_and_a2_components():
    gram=spatial_gram()
    roots={tuple(int(i==j) for j in range(9)) for i in range(9)}
    pending=list(roots)
    while pending:
        r=pending.pop()
        for i in range(9):
            pairing=sum(gram[i][j]*r[j] for j in range(9))
            s=tuple(int(r[j]-pairing*int(j==i)) for j in range(9))
            if s not in roots:roots.add(s);pending.append(s)
    assert len(roots)==132
    assert sum(all(v==0 for v in r[7:]) for r in roots)==126
    assert sum(all(v==0 for v in r[:7]) for r in roots)==6

@pytest.mark.parametrize('name,walls',[('phi',29),('2phi',30)])
def test_standalone_escape_replay(name,walls):
    path=ROOT/f'results/fresh_root_block_h9_{name}_h30.json'
    data=json.loads(path.read_text())
    assert data['completed'] and data['wall_count']==walls
    verify(path)

def test_standalone_replay_rejects_forged_escape(tmp_path):
    data=json.loads((ROOT/'results/fresh_root_block_h9_phi_h30.json').read_text())
    data['hyperideal_witness']['norm']={'a':'1','b':'0'}
    bad=tmp_path/'bad.json';bad.write_text(json.dumps(data))
    with pytest.raises(AssertionError):verify(bad)
