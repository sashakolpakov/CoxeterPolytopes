#!/usr/bin/env python3
"""Complete rank62 finite budget after proved D17=120, no other deletion.

Reuses the frozen exact certificate engine. Generation builds the changed
tables once; --verify checks saved Bellman upper/attainment certificates.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_ambient62_payment_frontier_20261018 as B

ROOT=B.ROOT
TAG='fresh_finite_d17_cost120_rank62_budget_20261018'
LOCAL='fresh_finite_d17_five_unit_reduction_20261018'
THREE='fresh_ambient62_three_transfer_frontier_20261018'
STEMS=(LOCAL,THREE,B.TAG,B.FINITE,
       'fresh_order13_d_center_pair_zero_transfer_20261018',
       'fresh_order17_e6_unit_coupon_zero_transfer_20261018')
CHANGES=dict(D17='5',D13_zero_coupons=True)


def actual_base(t):
    return B.original_base(t)-(10 if t==('D',17) else 0)


def cross(a,b,upper):
    value=(B.upper_cross if upper else B.actual_cross)(a,b)
    for small,large in ((a,b),(b,a)):
        if large==('D',10):
            value+=2 if small==('A',2) else 1 if small==('A',3) else 0
    if not upper and {a,b}=={('E',6),('D',11)}:
        value+=4
    return value


def cost(names,upper=False):
    ts=[(s[0],int(s[1:])) for s in names]
    return sum(actual_base(t) for t in ts)+sum(cross(a,b,upper)
                for i,a in enumerate(ts) for b in ts[i+1:])


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for p,h in doc['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    local=docs[LOCAL]
    assert local['old_cost']=='125' and local['new_cost']=='120' and local['total_saving']=='5'
    assert local['new_finite_coupons']==[] and local['all_other_pair_weights_unchanged']
    assert docs[THREE]['finite']['exact_maximum']=='441'
    assert docs[THREE]['finite']['inherited_finite_scenario_index']==1
    assert docs[STEMS[-1]]['finite_coupon']['new_total']=='4'
    scenario=B.finite_scenario(docs[B.FINITE],CHANGES,
                               None if saved is None else saved['certificate'])
    assert scenario['conditional_changes']==CHANGES
    # The reused engine's actual-E6 value is the old2. Its upper value is
    # already4. Recount the true current cost, adding the missing doubled4.
    comparisons=[]
    for row in scenario['alternatives']:
        count=Counter(row['forest'])
        actual=row['actual_doubled']+4*count['E6']*count['D11']
        assert actual==cost(row['forest'])<=row['doubled_upper']==cost(row['forest'],True)
        comparisons.append(dict(counts=row['counts'],forest=row['forest'],
                                actual_doubled=actual,upper_doubled=row['doubled_upper']))
    types=([('A',r) for r in range(1,63)]+[('B',r) for r in range(2,63)]+
           [('D',r) for r in range(4,63)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    for a in types:
        assert actual_base(B.representative(a))>=actual_base(a)
        for b in types:
            assert cross(B.representative(a),B.representative(b),True)>=cross(a,b,False)
    maximum=max(r['upper_doubled'] for r in comparisons)
    assert str(F(maximum,2))==scenario['exact_maximum']
    attainers=sorted({tuple(r['forest']) for r in comparisons if r['actual_doubled']==r['upper_doubled']==maximum})
    assert attainers
    violations=sorted((r for r in comparisons if r['actual_doubled']>868),
                      key=lambda r:(-r['actual_doubled'],r['forest']))
    # Expand only the actual A16/D16 type ties at an attaining forest;
    # each expanded forest is literally recounted, not assumed equivalent.
    expanded=set()
    for names in attainers:
        for replace in range(names.count('D16')+1):
            rest=list(names)
            for _ in range(replace):rest.remove('D16');rest.append('A16')
            if cost(rest)==maximum:expanded.add(tuple(sorted(rest)))
    paths=[Path(__file__),Path(B.__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
             for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE_D17_COST120_COMPLETE_RANK62_BUDGET',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        ambient_dimension=62,finite_maximum=str(F(maximum,2)),finite_curvature='434',
        finite_slack=str(F(434)-F(maximum,2)),actual_attaining_forests=attainers,
        expanded_A16_D16_actual_attainers=sorted(expanded),
        proved_component_change=dict(D17_old='125',D17_new='120',new_coupons=[]),
        retained_component_costs=dict(A15='105',D15='105',A16='114',B16='104',D16='114',A17='114'),
        E6_D11_actual_total='4',D13_zero_increments=dict(A2_D10='1',A3middle_D10='1/2'),
        simultaneous_type_pairs=len(types)**2,certificate=scenario,
        current_actual_comparisons=comparisons,actual_representative_forests_above434=violations,
        complete_large_vectors=len(comparisons),interaction_keys=len(scenario['tables']),
        no_other_conditional_component_drop=True,new_dimension_bound=None,
        scope='Only proved D17 reduction5 is applied to the complete current three-transfer finite rule. Current E6/D11 total4 and D13 zero increments are included. Exact upper and attaining certificates cover every rank62 forest; no unproved rank15/A16/D16 reductions, geometric census, global model optimization or ambient theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    saved=json.loads(path.read_text()) if args.verify else None
    out=json.loads(json.dumps(run(saved)))
    if args.verify:assert out==saved
    else:
        assert not path.exists(),'Existing output; use --verify.'
        path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'maximum',out['finite_maximum'],'slack',out['finite_slack'],
          'attainers',out['actual_attaining_forests'],flush=True)
