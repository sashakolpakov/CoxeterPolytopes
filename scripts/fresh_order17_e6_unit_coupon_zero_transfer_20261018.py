#!/usr/bin/env python3
"""Delete four affine D16 halves using whole E6/D11 unit coupons.

Exact saved-interface composition and pointwise finite-budget absorption.
No geometric census or Bellman optimization is repeated here.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order17_e6_unit_coupon_zero_transfer_20261018'
OLD='fresh_order17_e6_half_coupon_transfer_20261018'
GEO='fresh_order17_affine_d_three_quarter_transfer_20261017'
CURRENT='fresh_sparse_coupon_order11_13_dimension62_20261018'
BASE='fresh_order11_order17_dimension63_20261018'
BUDGET='fresh_sparse_coupon_rank63_forest_budget_20261018'
DZERO='fresh_order13_d_center_pair_zero_transfer_20261018'
BZERO='fresh_order16_affine_b_singleton_projection_transfer_20261018'
STEMS=(OLD,GEO,CURRENT,BASE,BUDGET,DZERO,BZERO,
       'fresh_order11_c_endpoint_zero_transfer_20261018',
       'fresh_order12_b_center_terminal_transfer_20261018',
       'fresh_order13_d_center_pair_half_transfer_20261018',
       'fresh_cycle11_half_pair_transfer_20261018',
       'fresh_order11_bcd_half_pair_transfer_20261018',
       'fresh_depth11_finite8_zero_transfer_20261018',
       'fresh_deep_finite8_cycle_interface_20261017',
       'fresh_finite7_ideal5_triangle_payments_20261017')
CHANGE={(3,15),(3,16),(9,13),(9,14)}


def distances(adj,start):
    d={start:0};q=deque([start])
    while q:
        i=q.popleft()
        for j in adj[i]:
            if j not in d:d[j]=d[i]+1;q.append(j)
    return d


def role_check(gram,role):
    es,ds=role['small_nodes'],role['large_nodes']
    assert role['small_type']=='E6' and role['large_type']=='D11'
    assert len(es)==6 and len(ds)==11 and role['whole_components']
    assert not set(es)&set(ds) and all(gram[i][j]==0 for i in es for j in ds)
    output={}
    for name,ids,mark in [('E6',es,role['small_mark']),('D11',ds,role['large_mark'])]:
        adj={i:[j for j in ids if i!=j and gram[i][j]] for i in ids}
        assert all(gram[i][j]==-1 for i in ids for j in adj[i])
        forks=[i for i in ids if len(adj[i])==3];assert len(forks)==1
        fork=forks[0];dist=distances(adj,fork);assert set(dist)==set(ids)
        leaves=[i for i in ids if len(adj[i])==1]
        lengths=sorted(dist[i] for i in leaves)
        assert lengths==([1,2,2] if name=='E6' else [1,1,8])
        eligible=sorted(i for i in leaves if dist[i]==(2 if name=='E6' else 1))
        assert len(eligible)==2 and mark in eligible
        output[name]=dict(nodes=ids,fork=fork,leaf_arm_lengths=lengths,
                          eligible_marks=eligible,actual_mark=mark)
    return output


def absorption(budget):
    allowance=budget['representative_allowance']
    assert allowance['types']==['E6','D11'] and allowance['actual_total']=='2'
    assert allowance['upper_total']=='4'
    rows=[];changed=0
    for r in budget['type_pair_checks']:
        pair={tuple(t) for t in r['types']}
        extra=4 if pair=={('E',6),('D',11)} else 0
        new=r['actual_doubled']+extra
        assert new<=r['upper_doubled']
        if extra:
            changed+=1;assert new==r['upper_doubled']==8
            rows.append(dict(types=r['types'],old_actual_doubled=r['actual_doubled'],
                             new_actual_doubled=new,unchanged_upper_doubled=r['upper_doubled']))
    assert len(budget['type_pair_checks'])==36864 and changed==2
    checks=[]
    for r in budget['alternatives']:
        c=Counter(r['formal_upper_attainer']);extra=4*c['E6']*c['D11']
        new=r['actual_doubled_at_that_forest']+extra
        assert new<=r['doubled_upper']
        checks.append(dict(counts=r['counts'],new_actual_doubled=new,
                           unchanged_upper_doubled=r['doubled_upper']))
    assert len(checks)==1329 and max(r['unchanged_upper_doubled'] for r in checks)==898
    witness=budget['actual_attaining_forest']
    assert witness==['D17','D17','D29'] and budget['actual_attainer_cost']=='449'
    assert 'E6' not in witness and 'D11' not in witness
    return dict(ordered_type_pairs_checked=36864,changed_ordered_type_pairs=rows,
                saved_formal_attainers_checked=1329,formal_attainer_comparisons=checks,
                unchanged_representative_E6_D11_allowance='4',
                exact_rank63_maximum_for_this_upgrade_only='449',attaining_forest=witness,
                Bellman_recomputed=False,
                excludes_separate_D13_zero_coupon_increments=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    old,geo=docs[OLD],docs[GEO]
    assert docs[CURRENT]['dimension_bound']==62
    assert old['new_cost']=='118' and old['new_ideal_weight']=='1/2'
    assert old['inherited_complete_patterns']==140 and old['inherited_affected_maps']==96
    assert old['inherited_proper_exclusions']==64 and set(map(tuple,old['changed_pairs']))==CHANGE
    assert geo['outcomes']==dict(excluded_projection_above_norm=32,
        excluded_multiple_core_contact_roots=16,excluded_forced_lorentz_index=16,
        sealed_E6_D11_quarter_coupon_maps=32)
    weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values(),F())==118
    for p in CHANGE:assert weights[p]==F(1,2);weights[p]=F()
    assert sum(weights.values(),F())==116
    target=old['target'];adj={i:[j for j,m in enumerate(row) if m] for i,row in enumerate(target)}
    assert all(distances(adj,i)[j]==10 for i,j in CHANGE)
    for mp in geo['automorphisms']:assert {tuple(sorted(mp[i] for i in p)) for p in CHANGE}==CHANGE
    assert all(weights[tuple(p)]==1 for p in geo['basic_pairs_retained']+geo['triangle_pairs_retained'])
    short8=[(i,j) for i in range(17) for j in range(i+1,17) if distances(adj,i)[j]<=8]
    assert all(weights[p]==1 for p in short8)
    rows=[]
    for case in geo['cases']:
        if case['index'] not in geo['whole_affine_contexts']:continue
        proof=case['proof'];g=[[F(x) for x in r] for r in proof['core_gram']]
        z=list(map(F,proof['radical']));assert len(g)==19 and proof['signature']==[17,1,1]
        assert all(sum(g[i][j]*z[j] for j in range(19))==0 for i in range(19))
        aff=proof['affine_nodes'];ker=list(map(F,proof['positive_affine_kernel']))
        assert len(aff)==17 and all(x>0 for x in ker)
        assert all(sum(g[i][j]*ker[k] for k,j in enumerate(aff))==0 for i in aff)
        assert proof['allowed_contacts']==[3,13] and z[3]*z[13]<0 and proof['both_contacts_shortcut']==2
        new=[];role=None;original=None;coupon_corner=None
        for c in proof['corners']:
            oi,of=F(c['used_ideal_payment']),F(c['used_new_finite_coupon'])
            if oi:assert oi==F(3,4) and not of;original=c['marks']
            if of:
                assert of==F(1,4) and len(c['new_coupon_roles'])==1
                role=role_check(g,c['new_coupon_roles'][0]);coupon_corner=c['marks']
            assert all(not x['used_in_this_proof'] for x in c.get('unused_same_component_payments',[]))
            new.append(dict(marks=c['marks'],ideal_payment='0',finite_coupon=str(int(bool(of))),
                            total=str(int(bool(of)))))
        assert original is not None and coupon_corner is not None and original!=coupon_corner
        assert role is not None and len(proof['maps'])==8 and sum(F(c['total']) for c in new)==1
        saved=next(r for r in old['changed_cases'] if r['case']==case['index'])
        assert role==saved['coupon_role']
        rows.append(dict(case=case['index'],maps=8,coupon_role=role,original_marks=original,
                         distinct_unit_marks=coupon_corner,corners=new,total='1'))
    assert [r['case'] for r in rows]==[1981,2374,2549,2551]
    assert docs[DZERO]['new_affine_total']=='73' and docs[BZERO]['new_order16_maximum']=='106'
    short_interfaces=docs[BASE]['retained_small_tree_short8_interfaces']
    assert len(short_interfaces)==18
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
                           for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER17_D116_E6_D11_UNIT_COUPON_ZERO_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        inherited_complete_patterns=140,inherited_affected_maps=96,inherited_proper_exclusions=64,
        changed_cases=rows,target=target,changed_pairs=sorted(CHANGE),old_ideal_weight='1/2',new_ideal_weight='0',
        old_cost='118',new_cost='116',weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(weights.items())],
        retained_short8_pairs=short8,inherited_small_tree_short8_interfaces=short_interfaces,
        finite_coupon=dict(small_type='whole E6',small_mark='endpoint of either length2 arm',
            large_type='whole D11',large_mark='short fork leaf',old_weight='1/2',new_weight='1',
            marked_pairs_per_component_pair=4,old_total='2',new_total='4',added_total='2'),
        new_order17_costs=dict(A=85,B=111,C=108,D=116),new_order17_envelope='116/17',
        finite_budget_absorption=absorption(docs[BUDGET]),uses_finite_same_component_units=False,
        triangle_minimum=2,quadrilateral_minimum=1,
        actual_corner_priority=['C11 endpoint zero','B12 center-terminal zero','D13 center half or upgraded zero',
            'cycle11 distance5 half','coupled BCD11 halves','changed order17 D orbit','no changed order17 corner'],
        compatible_with_B16_singleton_deletions=True,compatible_with_D13_center_zero=True,
        new_dimension_bound=None,
        scope='Complete local zero-transfer composition: inherited arbitrary-support proper exclusions, four rechecked radical/kernel and whole finite role interfaces, distinct full finite units, independent current short8 priority. Pointwise E6/D11 budget absorption is proved without rerunning census or Bellman tables; separate D13 zero increments and global frontier remain separate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'ideal',out['new_cost'],'finite increment',out['finite_coupon']['added_total'],flush=True)
