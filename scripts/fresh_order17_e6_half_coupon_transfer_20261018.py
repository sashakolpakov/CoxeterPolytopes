#!/usr/bin/env python3
"""Refine four ideal D16 weights to halves using the same sealed finite coupon.

No geometric census is rerun: frozen complete support and sealing are inputs.
This script independently checks literal roles, weight interfaces and totals.
"""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order17_e6_half_coupon_transfer_20261018'
OLD='fresh_order17_affine_d_three_quarter_transfer_20261017'
CURRENT='fresh_deep_finite8_dimension64_20261018'
CHANGE={(3,15),(3,16),(9,13),(9,14)}

def distances(adj,start):
    d={start:0};q=deque([start])
    while q:
        i=q.popleft()
        for j in adj[i]:
            if j not in d:d[j]=d[i]+1;q.append(j)
    return d

def role_check(gram,role):
    es,ds=role['small_nodes'],role['large_nodes'];em,dm=role['small_mark'],role['large_mark']
    assert role['small_type']=='E6'and role['large_type']=='D11'and len(es)==6 and len(ds)==11
    assert role['whole_components']and not set(es)&set(ds)
    assert all(gram[i][j]==0 for i in es for j in ds)
    output={}
    for name,ids,mark in [('E6',es,em),('D11',ds,dm)]:
        adj={i:[j for j in ids if i!=j and gram[i][j]]for i in ids}
        assert all(gram[i][j]==-1 for i in ids for j in adj[i])
        forks=[i for i in ids if len(adj[i])==3];assert len(forks)==1
        fork=forks[0];dist=distances(adj,fork);assert set(dist)==set(ids)
        leaves=[i for i in ids if len(adj[i])==1]
        lengths=sorted(dist[i]for i in leaves);assert lengths==([1,2,2]if name=='E6'else[1,1,8])
        eligible=sorted(i for i in leaves if dist[i]==(2 if name=='E6'else 1))
        assert len(eligible)==2 and mark in eligible
        output[name]=dict(nodes=ids,fork=fork,leaf_arm_lengths=lengths,eligible_marks=eligible,actual_mark=mark)
    return output

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in(OLD,CURRENT)}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    old=docs[OLD];assert old['new_cost']=='119'and old['new_changed_weight']=='3/4'
    assert set(map(tuple,old['changed_pairs']))==CHANGE
    assert old['complete_eligible_patterns']==140 and old['affected_maps']==96
    expected=Counter(excluded_projection_above_norm=32,excluded_multiple_core_contact_roots=16,
        excluded_forced_lorentz_index=16,sealed_E6_D11_quarter_coupon_maps=32)
    assert old['outcomes']==expected and old['whole_affine_contexts']==[1981,2374,2549,2551]
    weights={tuple(r['pair']):F(r['weight'])for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values(),F())==119
    for pair in CHANGE:assert weights[pair]==F(3,4);weights[pair]=F(1,2)
    target=old['target'];adj={i:[j for j,m in enumerate(row)if m]for i,row in enumerate(target)}
    for mp in old['automorphisms']:assert {tuple(sorted(mp[i]for i in p))for p in CHANGE}==CHANGE
    assert all(distances(adj,i)[j]==10 for i,j in CHANGE)
    assert all(weights[tuple(p)]==1 for p in old['basic_pairs_retained']+old['triangle_pairs_retained'])
    assert all(weights[i,j]==1 for i in range(17)for j in range(i+1,17)if distances(adj,i)[j]<=8)
    rows=[];proper=0;seen=[]
    for case in old['cases']:
        proof=case['proof']
        if case['index']not in old['whole_affine_contexts']:
            assert all(r['outcome'].startswith('excluded_')for r in proof['maps']);proper+=len(proof['maps']);continue
        seen.append(case['index']);g=[[F(x)for x in r]for r in proof['core_gram']]
        assert len(g)==19 and proof['signature']==[17,1,1]and proof['both_contacts_shortcut']==2
        z=list(map(F,proof['radical']));assert all(sum(g[i][j]*z[j]for j in range(19))==0 for i in range(19))
        aff=proof['affine_nodes'];ker=list(map(F,proof['positive_affine_kernel']))
        assert len(aff)==17 and all(x>0 for x in ker)and all(sum(g[i][j]*ker[k]for k,j in enumerate(aff))==0 for i in aff)
        assert proof['allowed_contacts']==[3,13]and z[3]*z[13]<0
        new=[];old_total=F();new_total=F();role=None
        for c in proof['corners']:
            oi,of=F(c['used_ideal_payment']),F(c['used_new_finite_coupon']);old_total+=oi+of
            if oi:
                assert oi==F(3,4)and of==0 and c['kind']=='affineD16_order17';ni,nf=F(1,2),F()
            elif of:
                assert of==F(1,4)and len(c['new_coupon_roles'])==1
                role=role_check(g,c['new_coupon_roles'][0]);ni,nf=F(),F(1,2)
            else:ni=nf=F()
            assert all(not x['used_in_this_proof']for x in c.get('unused_same_component_payments',[]))
            new_total+=ni+nf;new.append(dict(marks=c['marks'],ideal_weight=str(ni),finite_coupon=str(nf),total=str(ni+nf)))
        assert old_total==new_total==1 and role is not None and len(proof['maps'])==8
        rows.append(dict(case=case['index'],changed_maps=8,coupon_role=role,corners=new,total='1'))
    assert proper==64 and seen==old['whole_affine_contexts']and sum(weights.values(),F())==118
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in(OLD,CURRENT)for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER17_D118_E6_D11_HALF_COUPON',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        inherited_complete_patterns=140,inherited_affected_maps=96,inherited_proper_exclusions=64,changed_cases=rows,
        target=target,changed_pairs=sorted(CHANGE),old_ideal_weight='3/4',new_ideal_weight='1/2',old_cost='119',new_cost='118',
        weighted_pairs=[dict(pair=p,weight=str(w))for p,w in sorted(weights.items())],
        finite_coupon=dict(small_type='whole E6',small_mark='endpoint of either length2 arm',large_type='whole D11',large_mark='short fork leaf',
            old_weight='1/4',new_weight='1/2',marked_pairs_per_component_pair=4,old_total='1',new_total='2',added_total='1'),
        new_order17_costs=dict(A=85,B=111,C=108,D=118),new_order17_envelope='118/17',
        uses_finite_same_component_units=False,other_ideal_weights_unchanged=True,triangle_minimum=2,quadrilateral_minimum=1,
        compatible_with_deep11_zero=True,compatible_with_finite_D16_D17_exchange=True,
        finite_budget_proved_here=False,new_dimension_bound=None,
        scope='Exact weight refinement of the frozen complete order17 transfer. Literal19-root radical/kernel and coupon roles rechecked; whole-component sealing and complete140-pattern/96-map coverage inherited. New finite obligation is1 additional unit per wholeE6-D11 pair. No new geometric census or global dimension theorem.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'ideal',out['new_cost'],'finite increment',out['finite_coupon']['added_total'],flush=True)
