"""Independent replay and adversarial checks for compatible H9 cut unions."""
import copy
from pathlib import Path
import pytest
from scripts.fresh_h4_blocks_independent_verify import Z,O,P,H4,sub,neg,mul,div,scale,pair,matrix,parse
from scripts.fresh_multibasepoint_h9_verify import verify as verify_merge,valid_coxeter_value,read
from scripts.fresh_multibasepoint_h9_targeted_verify import verify as verify_targeted,plane

ROOT=Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('name,walls',[('2phi',364),('sqrt5',85)])
def test_complete_merged_certificate_and_all_angles(name,walls):
    result=verify_merge(ROOT/f'results/fresh_multibasepoint_h9_{name}_20261013.json')
    assert result['walls']==walls

def test_all_targeted_cut_orbits_and_final_escape():
    result=verify_targeted(ROOT/'results/fresh_multibasepoint_h9_2phi_targeted_rounds_20261013.json')
    assert result['walls']==372 and result['cut_orbits']==4 and result['new_walls']==8
    assert result['new_finite_angle_counts']=={2:62,4:2}

def test_acute_integral_h4_pair_at_2pi_over5_is_rejected():
    # These are actual norm-two integral H4 roots. Acuteness alone is insufficient.
    r=(O,Z,Z,Z);s=(neg(P),neg(P),Z,Z)
    q=pair(r,H4,r);other=pair(s,H4,s);product=pair(r,H4,s)
    assert q==other==scale(O,2) and product==sub(O,P)
    ratio=div(scale(mul(product,product),4),mul(q,other))
    assert ratio==sub(scale(O,2),P) and not valid_coxeter_value(ratio)

def test_rejects_forged_targeted_cut_norm():
    data=read('results/fresh_multibasepoint_h9_2phi_escape_plane_20261013.json')
    base=read(data['source']);data=copy.deepcopy(data);data['cuts'][0]['norm']={'a':'2','b':'0'}
    with pytest.raises(AssertionError):
        plane(data,matrix(base['all_walls']),matrix(base['gram']),base['hyperideal_witness'],tuple(map(parse,base['strict_interior'])))

def test_rejects_forged_local_obstruction():
    data=read('results/fresh_multibasepoint_h9_sqrt5_escape_plane_20261013.json')
    base=read(data['source']);data=copy.deepcopy(data)
    data['norm_searches'][0]['local_obstructions'][0]['required_square_residue']=0
    with pytest.raises(AssertionError):
        plane(data,matrix(base['all_walls']),matrix(base['gram']),base['hyperideal_witness'],tuple(map(parse,base['strict_interior'])))
