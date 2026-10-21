#!/usr/bin/env python3
"""Retain the complete face proof after lowering one unused broad role.

The universal proof is in the report and its bound provenance inventory.
This checker verifies intrinsic multiplicities, frozen output interfaces
and simultaneous actual/representative costs. It runs no forest DP.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_d12_half_face_composition_20261018 as B
import fresh_d12_type_split_ambient58_forest_frontier_20261018 as C
import fresh_small_d15_fork_coupon_provenance_20261018 as I

ROOT=B.ROOT
TAG='fresh_small_d15_center_three_eighth_retention_20261018'
PINS={
 f'scripts/{I.TAG}.py':'a83e813335848fe314109140298f545d969592faa3285aa25be618eccb840602',
 f'results/{I.TAG}.json':'3e535d1a66fa3547fcc5b82d6e991d0778b994188147067b43dcf681b094318a',
 f'results/{I.TAG}.md':'750390cd694464a2820154444cef4c3e96865915dd53408c9d82b73dc0327e9d',
 f'scripts/{B.TAG}.py':'0f66fca4eb6084d013ff46a24b491d10dd272332969705c86de2f874a9819155',
 f'results/{B.TAG}.json':'a2ab5675edc3808e8ebee40770c7cb68acb8df7b9a7832391c2eb9b01e71016b',
 f'results/{B.TAG}.md':'9f697fed9259db6a96ca94376ea69691794f9376995bdedbeb968a7851ba0a47',
}


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def debit(a,b):
    return sum(4*s[1]for s,t in((a,b),(b,a))if 3<=s[1]<=5 and t==('D',15))
def actual_base(t):return B.actual_base(t)
def upper_base(t):return C.upper_base(t)
def actual_cross(a,b):return B.actual_cross(a,b)-debit(a,b)
def upper_cross(a,b):return C.upper_cross(a,b)-debit(a,b)
def context(forest,upper=False):return C.context(forest,upper)


def run():
    for path,h in PINS.items():assert digest(path)==h,path
    inv=json.loads((ROOT/f'results/{I.TAG}.json').read_text())
    prior=json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    for d in(inv,prior):
        for path,h in d['sources'].items():assert digest(path)==h,path
    assert inv['status']=='PASS_BOUNDED_SMALL_D15_FORK_COUPON_PROVENANCE'
    assert inv['baseline']==B.TAG and prior['ideal_weights_changed'] is False
    old=inv['literal_inventory'];later=inv['later_finite_inventory']
    assert(old['records'],old['corners'],later['corner_count'])==(216,832,564)
    assert old['small_rank3_to5_D15_corners']==later['D15_cross_corners']==0
    assert len(old['D15_cross_corners'])==10
    assert all(r['small_type']=='A1'for r in old['D15_cross_corners'])
    assert old['role_histogram']=={'fork_center':2,'short_fork_leaf':4,'long_arm':4}
    assert len(later['packets'])==22
    # A simple D15 tree has one intrinsic center and a two-node leaf orbit.
    edges=[(i,i+1)for i in range(12)]+[(12,13),(12,14)]
    adj={i:[]for i in range(15)}
    for a,b in edges:adj[a].append(b);adj[b].append(a)
    centers=[i for i in adj if len(adj[i])==3]
    assert centers==[12] and [i for i in adj[12]if len(adj[i])==1]==[13,14]
    rules=[]
    for rank in (3,4,5):
        old_weights=[[i,j,'1/2']for i in range(rank)for j in(12,13,14)]
        new_weights=[[i,j,'3/8'if j==12 else'1/2']for i in range(rank)for j in(12,13,14)]
        before=sum(F(w)for _,_,w in old_weights);after=sum(F(w)for _,_,w in new_weights)
        assert before==F(3*rank,2) and after==F(11*rank,8) and before-after==F(rank,8)
        assert sorted((i,27-j if j in(13,14)else j,w)for i,j,w in new_weights)==sorted(map(tuple,new_weights))
        rules.append(dict(small_rank=rank,old_weights=old_weights,new_weights=new_weights,
            old_component_pair_cost=before,new_component_pair_cost=after,saving=before-after))
    changed=[]
    for a in C.TYPES:
        assert actual_base(a)==C.actual_base(a)<=upper_base(C.representative(a))
        for b in C.TYPES:
            assert B.actual_cross(a,b)==C.actual_cross(a,b)
            x=actual_cross(a,b);u=upper_cross(C.representative(a),C.representative(b))
            assert 0<=x==actual_cross(b,a)<=B.actual_cross(a,b) and x<=u
            assert u<=C.upper_cross(C.representative(a),C.representative(b))
            if debit(a,b):changed.append(dict(a=a,b=b,old=F(B.actual_cross(a,b),32),new=F(x,32),saving=F(debit(a,b),32)))
    assert len(C.TYPES)==177 and len(changed)==22
    assert actual_cross(('A',1),('D',15))==50
    for small in C.TYPES:
        if small[1]==2:assert actual_cross(small,('D',15))==96
        if 3<=small[1]<=5:
            assert actual_cross(small,('D',15))==44*small[1]
            assert actual_cross(small,('A',15))==16*small[1]
    paths=set(PINS)|set(inv['sources'])|set(prior['sources'])|{f'scripts/{TAG}.py'}
    paths.update(f'{folder}/{C.TAG}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return I.enc(dict(status='PASS_SMALL_D15_CENTER_THREE_EIGHTH_FACE_RETENTION',
        sources={p:digest(p)for p in sorted(paths)},baseline=B.TAG,provenance=I.TAG,scale=32,
        changed_rule=dict(small='whole finite connected component of rank3,4,5; any mark',
            large='whole finite D15',role='unique fork center',canonical_mark=12,old='1/2',new='3/8'),
        component_pair_rules=rules,changed_ordered_type_pairs=changed,simultaneous_type_comparisons=len(C.TYPES)**2,
        retained_D15_roles=dict(A1_total='25/16',rank2_total='3',short_leaves='1/2',A1_node11='1/16'),
        literal_interface=dict(records=216,corners=832,changed_coupon_debits=0,all_ten_D15_small_sides='A1'),
        later_finite_interface=dict(packets=22,corners=564,changed_coupon_debits=0),
        output_routes=inv['universal_output_inventory'],
        same_component_weights_changed=False,ideal_weights_changed=False,contextual_rules_changed=False,
        new_coupon_cost='0',triangle_minimum='2',quadrilateral_minimum='1',
        residual_item_delta='-4*r*c15 at ranks3,4; rank5 fixed contribution -20*c5*c15, in units1/32',
        forest_DP_performed=False,new_forest_maximum_claimed=False,completed_unrestricted_bound=58,new_dimension_bound=None,
        scope='Narrow whole-face proof retention by explicit universal and replacement-output inventory. Saved-corner checks supplement inherited whole-component proofs. No new geometry census, forest optimization or dimension bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['simultaneous_type_comparisons'],flush=True)
