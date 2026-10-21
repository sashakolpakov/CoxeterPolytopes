#!/usr/bin/env python3
"""Compose the order14 B/C deletions and close every ambient62 cusp margin.

The complete preceding scalar rows are inherited, not reconvolved. This
checks their exact cost update and independent complete profile domain.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order14_ambient62_ideal_margin_closure_20261018'
B='fresh_order14_b_terminal_residual_transfer_20261018'
C='fresh_order14_c_inner_endpoint_projection_transfer_20261018'
PRE='fresh_ambient62_three_transfer_frontier_20261018'
CURRENT='fresh_sparse_coupon_order11_13_dimension62_20261018'
STEMS=(B,C,PRE,CURRENT)


def high_parts(total,minimum=9):
    if not total:yield ();return
    for r in range(minimum,total+1):
        for tail in high_parts(total-r,r):yield (r,)+tail


def summary(rows):
    n=min(F(r['normalized_margin']) for r in rows)
    a=min(F(r['absolute_margin']) for r in rows)
    return dict(normalized_minimum=str(n),
        normalized_minimizers=[r['ranks'] for r in rows if F(r['normalized_margin'])==n],
        absolute_minimum=str(a),
        absolute_minimizers=[r['ranks'] for r in rows if F(r['absolute_margin'])==a],
        negative_count=sum(F(r['normalized_margin'])<0 for r in rows),
        zero_count=sum(F(r['normalized_margin'])==0 for r in rows))


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    b,c,old=docs[B],docs[C],docs[PRE]
    assert docs[CURRENT]['dimension_bound']==old['inherited_completed_ceiling']==62
    assert b['old_cost']==c['old_cost']==84 and b['new_cost']==82 and c['new_cost']==83
    assert b['deleted_pairs']==[[9,12],[9,13]] and c['deleted_pair']==[1,12]
    for d in (b,c):
        assert d['finite_weights_unchanged'] and not d['added_finite_coupons']
        assert d['triangle_minimum']==2 and d['quadrilateral_minimum']==1
        retained=set(map(tuple,d['retained_unit_pairs']))
        assert set(map(tuple,d['short8_pairs_retained']))<=retained
        assert set(map(tuple,d['basic_pairs_retained']))<=retained
        assert set(map(tuple,d['triangle_pairs_retained']))<=retained
    families={'A':70,'B':82,'C':83,'D':82};assert max(families.values())==83
    assert b['new_order14_family_costs']==dict(A=70,B=82,C=84,D=82)
    assert c['new_order14_family_costs']==dict(A='70',B='84',C='83',D='82')
    profiles=set(high_parts(61))
    for low in range(1,9):profiles.update((low,)+p for p in high_parts(61-low))
    previous=old['ideal']['current_rows']
    assert len(profiles)==len(previous)==old['ideal']['complete_profiles']==2641
    assert {tuple(r['ranks']) for r in previous}==profiles
    assert len({tuple(r['ranks']) for r in previous})==len(previous)
    prior_summary=summary(previous);saved_summary=old['ideal']['current_summary']
    assert prior_summary['negative_count']==saved_summary['negative']==3
    assert prior_summary['zero_count']==saved_summary['zero']==0
    assert prior_summary['normalized_minimum']==saved_summary['minimum']['normalized_margin']
    assert prior_summary['absolute_minimum']==saved_summary['minimum_absolute']['absolute_margin']
    assert sorted((r for r in previous if F(r['normalized_margin'])<0),key=lambda r:r['ranks'])==sorted(saved_summary['negative_rows'],key=lambda r:r['ranks'])
    rows=[]
    for r in previous:
        ranks=r['ranks'];V=prod(a+1 for a in ranks);assert V==r['link_vertices']
        count=ranks.count(13);saving=F(count,14)
        margin=F(r['normalized_margin'])+saving
        cost=F(r['normalized_cost'])-saving
        assert F(r['absolute_margin'])==V*F(r['normalized_margin'])
        rows.append(dict(ranks=ranks,link_vertices=V,rank13_factors=count,
            normalized_cost=str(cost),normalized_margin=str(margin),absolute_margin=str(V*margin),
            normalized_cost_saving=str(saving)))
    result=summary(rows)
    assert result==dict(normalized_minimum='2741/371712',
        normalized_minimizers=[[10,10,11,15,15]],absolute_minimum='73',
        absolute_minimizers=[[61]],negative_count=0,zero_count=0)
    assert old['finite']['exact_maximum']=='441' and old['finite']['curvature']=='434'
    assert old['finite']['finite_slack']=='-7' and old['strict_global_constant']=='4805'
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
                           for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_AMBIENT62_ORDER14_COMPLETE_IDEAL_MARGIN_CLOSURE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        ambient_dimension=62,inherited_unrestricted_ceiling=62,new_unrestricted_dimension_bound=None,
        order14_family_costs=families,old_order14_maximum='84',new_order14_maximum='83',
        old_rank13_envelope='6',new_rank13_envelope='83/14',
        finite_weights_unchanged=True,new_finite_coupons=[],
        triangle_minimum=2,quadrilateral_minimum=1,
        complete_profiles=2641,
        concentration_domain='Sum61; at most one part1..8; arbitrary repeated parts>=9; no high-part concentration.',
        previous_summary=summary(previous),current_summary=result,rows=rows,
        all_actual_cusp_absolute_surplus='73',
        finite_maximum='441',finite_curvature='434',finite_slack='-7',
        finite_actual_attaining_forest=old['finite']['actual_attaining_forest'],
        strict_global_constant='4805',necessary_vertex_inequality='7*v_f >=4805+73*v_i',
        all_ideal_class_dimension_bound=61,
        actual_corner_priority='Inherited C11/B12/D13/cycle11/BCD11 priorities first; current short8 route remains independent. In the remaining selected branch every changed B14/C14 occurrence is excluded by its complete local projection or residual proof; all other selected and nonselected payments remain. D17/B16 long-orbit proofs have the same preserved short-path outputs.',
        scope='Completed simultaneous local B/C14 rule and exact2641-row cost update. Curvature rows, all-rank low9 concentration, factor fibers, finite441 and global constant are inherited completed inputs. No Bellman or historical geometric census replay; finite deficit7 prevents a new unrestricted bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_profiles'],out['current_summary'],flush=True)
