#!/usr/bin/env python3
"""Large-component-only savings cannot lower the current rank65 finite maximum."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_rank65_small_component_obstruction_20261017'
BUDGET='fresh_a2_d11_full_coupon_forest_budget_20261017'
RANK17='fresh_finite_d17_long_pair_support_20261017'


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in(BUDGET,RANK17)}
    for s,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    row=docs[BUDGET]['certificate']['dimension_maxima'][65]
    assert row['doubled_cost']==974 and row['finite_curvature_coefficient']=='480'
    assert docs[BUDGET]['doubled_representative_cross_formula']==(
        'a*(2*d10+3*tD11to16+2*a17)+2*a1*(d11+d12)+2*a4*d13+4*a2*d11')
    d17=next(r for r in docs[RANK17]['records']if r['kind']=='D')
    assert len(d17['retained_pairs'])==d17['new_cost']==127
    d16=16*15//2;assert d16==120
    cost=127+3*d16;assert cost==487
    assert 17+3*16==65
    # No current cross coupon joins two components of ranks16 or17.
    assert all(min(a,b)>5 for a in(16,17)for b in(16,17))
    a1_four=4*d16+4*F(3,2);assert a1_four==486
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in(BUDGET,RANK17)for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_RANK65_LARGE_COMPONENT_ONLY_REDUCTION_OBSTRUCTION',
        rank=65,current_exact_maximum='487',curvature='480',
        unchanged_attainer=['D17','D16','D16','D16'],attainer_cost='487',
        second_forest=['A1','D16','D16','D16','D16'],second_cost='486',
        optimistic_zero_cost_rank_threshold=18,optimistic_maximum='487',
        necessary_deficits=dict(D17_plus_three_D16='7',A1_plus_four_D16='6'),
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        proof='Setting every payment involving a whole component of rank>=18 to zero cannot increase the old exact maximum487. The unchanged valid typed forest D17+3D16 still attains487, so the optimistic maximum remains487. This is an abstract finite-budget obstruction, not an actual hyperbolic vertex existence claim.',
        scope='The present rule cannot reach finite allowance480 through reductions confined to components of rank>=18. A local-cost solution must lower D16 or D17 same-component cost; changing small-component coupons alone cannot change the first witness. Actual-vertex restrictions or justified global compensation are alternative routes. No nonexistence theorem or new ambient bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['unchanged_attainer'],out['optimistic_maximum'],flush=True)
