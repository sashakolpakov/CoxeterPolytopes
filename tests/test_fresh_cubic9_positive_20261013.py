"""Regression gates for the changed coefficient and targeted A4 completion."""
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import fresh_cubic9_positive_20261013 as p

def test_norm_pool_changes_with_time_parity():
    odd=p.Host(p.T,2)
    assert p.O in odd.norms and p.K(4) not in odd.norms
    from fresh_cubic9_mixed_search_20261013 import NORMS,candidates
    assert p.K(4) in NORMS
    assert candidates(p.T*p.T-2,p.K(4))[0]  # Previously found primitive norm4.

def test_two_weight_slice_is_saturated_not_only_five_weight_rays():
    weights=(p.O,p.K(2),p.Z,p.Z)
    y=p.mv(p.AINV,weights)
    assert y==(p.K(2),p.K(3),p.K(2),p.O)
    assert all(t.is_integral for t in y)
    assert not (weights[0]/5).is_integral

def test_new_wall_cuts_the_original_end_weight_escape():
    base=json.loads((ROOT/'results/fresh_cubic9_positive_20261013_h9_theta_h30.json').read_text())
    cuts=json.loads((ROOT/'results/fresh_cubic9_positive_20261013_weights.json').read_text())
    h=p.Host(p.T,2);v=tuple(map(p.K.from_json,base['geometry']['hyperideal_witness']['vector']))
    walls=[tuple(map(p.K.from_json,row))for row in base['all_walls']]
    cov=[p.mv(h.G,r)for r in walls];norms=[p.dot(r,c)for r,c in zip(walls,cov)]
    roots=[tuple(map(p.K.from_json,row['vector']))for row in cuts['selected_new_roots']]
    assert any(p.pair(v,h.G,r).sign()>0 for r in roots)
    assert all(p.compatible(r,p.pair(r,h.G,r),cov,norms) for r in roots)
