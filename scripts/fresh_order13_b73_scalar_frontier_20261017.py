#!/usr/bin/env python3
"""Complete ambient65 scalar domain after the proved order13 B refinement."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order13_b73_scalar_frontier_20261017'
BASE='fresh_ideal5_weighted_dimension65_20261017'
LOCAL='fresh_order13_affine_b_fork_projection_transfer_20261017'
FINITE='fresh_a2_d11_full_coupon_forest_budget_20261017'


def profiles(left,minimum=1,prefix=(),low=False,high=False):
    if not left:yield prefix;return
    for r in range(minimum,left+1):
        if(r<10 and low)or(r>=18 and high):continue
        yield from profiles(left-r,r,prefix+(r,),low or r<10,high or r>=18)


def cost(r):
    return F(r,2)if r<=10 else {11:F(16,3),12:F(74,13),13:F(6),14:F(32,5),15:F(107,16),16:F(120,17)}[r]if r<=16 else F(163,23)


def scalar(rs):
    h=[1]
    for r in rs:
        new=[0]*(len(h)+r)
        for j,v in enumerate(h):
            for k in range(r+1):new[j+k]+=v
        h=new
    d=sum(rs);V=prod(r+1 for r in rs);G=sum(v*min(j,d-j)for j,v in enumerate(h))
    assert sum(h)==V and h==h[::-1]
    margin=F(5*G-d*V,2*V)-sum(map(cost,rs),F())
    return dict(ranks=list(rs),link_vertices=V,G=G,h_polynomial=h,margin=str(margin))


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in(BASE,LOCAL,FINITE)}
    for s,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[LOCAL]['new_order13_costs']==dict(A=65,B=73,C=73,D=74)
    assert docs[LOCAL]['finite_weights_unchanged'] and not docs[LOCAL]['added_finite_coupons']
    rows=[scalar(rs)for rs in profiles(64)];assert len(rows)==1650
    oldnegative=[r for r in rows if F(r['margin'])-F(r['ranks'].count(12),13)<0]
    negative=[r for r in rows if F(r['margin'])<0]
    minimum=min(rows,key=lambda r:F(r['margin']))
    assert len(oldnegative)==16 and len(negative)==7
    assert(minimum['margin'],minimum['ranks'])==('-35835/454597',[10,10,12,16,16])
    scenarios=[]
    for drop15 in(0,1,2):
        drops={14:drop15,16:1}
        values=[dict(ranks=r['ranks'],margin=str(F(r['margin'])+sum(F(r['ranks'].count(k)*v,k+1)for k,v in drops.items())))for r in rows]
        best=min(values,key=lambda r:F(r['margin']))
        scenarios.append(dict(additional_total_factor_cost_drops={str(k+1):v for k,v in drops.items()},
            negative_profiles=sum(F(r['margin'])<0 for r in values),minimum=best,rows=values))
    assert scenarios[1]['negative_profiles']==1 and scenarios[1]['minimum']==dict(ranks=[10,10,10,10,10,14],margin='-23/14641')
    assert scenarios[2]['negative_profiles']==0 and scenarios[2]['minimum']==dict(ranks=[9,10,10,10,10,15],margin='61/1936')
    max65=docs[FINITE]['certificate']['dimension_maxima'][65]
    assert max65['doubled_cost']==974 and max65['finite_slack']=='-7'
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in(BASE,LOCAL,FINITE)for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths.append(ROOT/'results/fresh_distance14_unrestricted_dimension156_20261017.md')
    return dict(status='PASS_ORDER13_B73_COMPLETE_AMBIENT65_SCALAR_DIAGNOSTIC',
        ambient_dimension=65,complete_profiles=1650,rows=rows,old_negative_profiles=16,
        current_negative_profiles=7,current_negative_rows=negative,minimum=minimum,
        conditional_scenarios=scenarios,finite_maximum='487',finite_curvature='480',finite_slack='-7',
        concentration='Arbitrary repetitions10..17; at most one1..9 and one>=18. Only linear low-cost and constant high-cost moves are used; no duplicate13 move.',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        scope='Proved order13 change reduces16 negatives to7. Additional order15/17 drops are conditional geometric targets, not proved here. Finite deficit7 remains. No new ambient theorem or example.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['complete_profiles'],out['current_negative_profiles'],out['minimum']['margin'],flush=True)
