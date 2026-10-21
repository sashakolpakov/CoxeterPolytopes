#!/usr/bin/env python3
"""Saved-row exact frontier after the proved order16 D105 transfer.

No convolution, forest optimization, or new local geometry is performed.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order16_d105_ambient61_frontier_20261018'
FRONT='fresh_ambient61_payment_frontier_20261018'
LOCAL='fresh_order16_d_distance9_half_transfer_20261018'
FAMILIES='fresh_sparse_coupon_order11_13_dimension62_20261018'
STEMS=(FRONT,LOCAL,FAMILIES)


def domain():
    out=[]
    def rest(left,minimum,prefix):
        if not left:out.append(prefix);return
        for r in range(minimum,left+1):rest(left-r,r,prefix+[r])
    rest(60,9,[])
    for r in range(1,9):rest(60-r,9,[r])
    return sorted(out)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[LOCAL]['new_order16_maximum']=='105'
    assert docs[LOCAL]['new_family_totals']==dict(A='80',B='105',C='103',D='105')
    assert docs[LOCAL]['new_finite_coupons']==[] and docs[LOCAL]['all_finite_weights_unchanged']
    family_rows=[r for r in docs[FAMILIES]['new_factor_interfaces'] if r['order']==11]
    families={r['kind']:r['total_cost'] for r in family_rows}
    assert len(family_rows)==4 and families==dict(A='99/2',B='54',C='54',D='53')
    old=docs[FRONT]['ideal'];expected=domain()
    assert len(expected)==2365 and sorted(r['ranks'] for r in old['rows'])==expected
    rows=[]
    for r in old['rows']:
        saving=F(r['ranks'].count(15),16)
        margin=F(r['normalized_margin'])+saving
        assert F(r['absolute_margin'])==r['link_vertices']*F(r['normalized_margin'])
        rows.append(dict(ranks=r['ranks'],link_vertices=r['link_vertices'],
            old_normalized_margin=r['normalized_margin'],proved_saving=str(saving),
            normalized_margin=str(margin),absolute_margin=str(r['link_vertices']*margin)))
    negative=[r for r in rows if F(r['normalized_margin'])<0]
    zero=[r for r in rows if F(r['normalized_margin'])==0]
    assert len(negative)==4 and not zero
    threshold=max(-F(r['normalized_margin'])*11/r['ranks'].count(10) for r in negative)
    assert threshold==F(1159,11616)
    # One rank10 A/D factor already improves the envelope by at least1/11.
    assert all(r['ranks'].count(10)>0 and F(r['normalized_margin'])+F(1,11)>0 for r in negative)
    conditional=[F(r['normalized_margin'])+F(r['ranks'].count(10),44) for r in rows]
    assert min(conditional)==F(3,250) and all(v>0 for v in conditional)
    assert min(v*r['link_vertices'] for v,r in zip(conditional,rows))==64
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_PROVED_ORDER16_D105_AMBIENT61_FRONTIER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        ambient_dimension=61,link_dimension=60,complete_profiles=2365,rows=rows,
        negative_count=4,zero_count=0,negative_rows=negative,
        normalized_minimum=str(min(F(r['normalized_margin']) for r in rows)),
        absolute_minimum=str(min(F(r['absolute_margin']) for r in rows)),
        current_order11_family_costs=families,
        rank10_nonmaximal_family_clears_each_negative=True,
        strict_order11_total_saving_threshold=str(threshold),
        conditional_order11_saving='1/4',conditional_normalized_minimum='3/250',
        conditional_absolute_minimum='64',conditional_geometry_proved=False,
        no_new_finite_cost=True,new_dimension_bound=None,
        scope='Bounded exact update of every frozen low9-only terminal row after the proved D105 transfer. The independent sorted-part domain is reconstructed, but old convolution, curvature, concentration and type bounds remain inherited. Four negative rows remain; only all-B/C rank10 factors can attain these envelope obstructions. Their physical feasibility is neither asserted nor excluded. No new order11 theorem or finite budget follows.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['negative_count'],out['normalized_minimum'],flush=True)
