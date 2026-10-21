"""Independent arithmetic and rejection checks for the silver F4 audit."""
from copy import deepcopy
import pytest
from coxeter_polytopes.exact.silver import Silver as S
from coxeter_polytopes.search.silver_path_host import bounded_silver_integers
from scripts.fresh_silver_f4_h9_20261013 import block_vectors
from scripts.fresh_silver_f4_independent_verify import (
    O,P,ETA,add,sub,scale,parse,matrix,read,bounded,angle_ok,zero_blocks,check_escape,verify_case,
)

@pytest.mark.parametrize('caps',[(S(1),S(2)),(S(2,1),S(4,-2)),(S(7,-2),S(2,1))])
def test_source_embedding_box_against_independent_bound(caps):
    expected=set(bounded(tuple((x.a,x.b) for x in caps)))
    assert {(x.a,x.b) for x in bounded_silver_integers(caps)}==expected

@pytest.mark.parametrize('q',[S(1),S(2,1),S(2),S(4,2),S(4)])
def test_zero_block_recursion_against_direct_product(q):
    actual,_=block_vectors(q,q)
    assert {((e.a,e.b),tuple((x.a,x.b) for x in r)) for e,r in actual}==set(zero_blocks((q.a,q.b)))

def test_non_coxeter_three_pi_over_eight_rejected():
    assert angle_ok(ETA)
    assert not angle_ok(sub(scale(O,2),P))

def test_forged_escape_norm_rejected():
    d=read('results/fresh_silver_f4_h9_unit_h30_20261013.json')
    w=deepcopy(d['geometry']['hyperideal_witness']);w['norm']['a']='110'
    with pytest.raises(AssertionError):check_escape(w,matrix(d['all_walls']),matrix(d['gram']))

def test_compact_calibration_has_all_five_finite_vertices():
    assert verify_case('h4_calibration')['finite_vertices']==5

def test_old_stabilizer_reflection_is_already_wall46():
    from scripts.fresh_silver_f4_stabilizer_check import verify
    d=verify()
    assert d['wall_index']==46 and d['old_fixed_witness_on_mirror']
