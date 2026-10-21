"""Bound completeness and algebraic obstruction regressions."""
from copy import deepcopy
import pytest
from coxeter_polytopes.exact.golden import Golden as G
from scripts.fresh_rank3_h9_search import alphabet
from scripts.fresh_rank3_h9_verify import verify_obstruction
from scripts.fresh_h4_blocks_independent_verify import Z,O,scale,matrix
from scripts.fresh_multibasepoint_h9_verify import read

@pytest.mark.parametrize('bound,cap',[(10,G(1)/14),(20,G(2)/14),(8,G(2))])
def test_ternary_alphabet_matches_independent_box(bound,cap):
    expected={G(a,b) for a in range(-30,31) for b in range(-30,31)
              if G(a,b).sign()>=0 and (bound-G(a,b)).sign()>=0 and (cap-G(a,b)**2).sign(1)>=0}
    assert set(alphabet(bound,cap))==expected

def obstruction_data():
    cones=read('results/fresh_rank3_h9_cones_20261013.json');old=read(cones['base_source'])
    z=(Z,*(scale(O,k) for k in (7,14,20,16,12,8,4,10)),Z)
    return cones,matrix(old['all_walls']),matrix(old['gram']),matrix(old['standard_block_swap_involution']),z

def test_fixed_slice_obstruction():assert verify_obstruction(*obstruction_data())

def test_forged_obstruction_is_rejected():
    data=list(obstruction_data());data[0]=deepcopy(data[0])
    data[0]['slice_obstruction']['identities'][0]['nonnegative_coefficients'][0]['a']='-7'
    with pytest.raises(AssertionError):verify_obstruction(*data)
