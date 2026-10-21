"""Golden H3/A6 enumeration and hostile certificate regressions."""
from copy import deepcopy
from fractions import Fraction
import pytest
from scripts.fresh_h3_blocks_h9_a6_search import G, PHI, vectors
from scripts.fresh_h3_blocks_h9_independent_verify import (
    H3,A6,O,P,scale,add,brute_block,finite_roots,read,matrix,escape,
)

def independent(x):return (Fraction(x.a),Fraction(x.b))

@pytest.mark.parametrize('which,B', [(0,H3),(1,A6)])
@pytest.mark.parametrize('q', [1,2,4,7,14,28])
def test_schur_enumerator_matches_independent_zero_height_box(which,B,q):
    found,_=vectors(which,G(q),G(q))
    converted={(independent(e),tuple(map(independent,r))) for e,r in found}
    assert converted==set(brute_block(B,scale(O,q),scale(O,q)))

@pytest.mark.parametrize('which,B', [(0,H3),(1,A6)])
def test_positive_height_block_matches_independent_cartesian_box(which,B):
    target=G(2)+2*PHI
    found,_=vectors(which,target,G(1))
    converted={(independent(e),tuple(map(independent,r))) for e,r in found}
    assert converted==set(brute_block(B,add(scale(O,2),scale(P,2)),O))

def test_full_finite_root_orbits_have_expected_counts():
    assert finite_roots(H3,30)==30
    assert finite_roots(A6,42)==42

@pytest.mark.parametrize('name', ['phi','2phi'])
def test_forged_spacelike_escape_norm_is_rejected(name):
    payload=read('results/fresh_h3_blocks_h9_a6_'+name+'_h30_20261013.json')
    witness=deepcopy(payload['hyperideal_witness'])
    witness['norm']['a']='999'
    with pytest.raises(AssertionError):
        escape(witness,matrix(payload['all_walls']),matrix(payload['gram']))
