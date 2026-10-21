#!/usr/bin/env python3
"""Remove an unused broad D16 coupon allowance and rebalance D114.

All historical complete geometry remains bound input. Reconstruct the216
current literal corners and four seals, certify the changed finite problem,
and update every row of the complete current ideal domain exactly.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_d15_balanced_ambient60_frontier_20261018 as P
import fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018 as M

ROOT=P.ROOT;B=P.B;Z=P.Z.Z;L=P.L;T=L.T;Q=L.Q;A=L.A;D=L.D;N=60
TAG='fresh_small_d16_coupon_and_ideal17_rebalance_20261018'
PROVENANCE=(
 'fresh_broad_finite_half_coupons_20261017',
 'fresh_finite8_a17_d11_coupon_budget_20261017',
 'fresh_ideal7_finite_cross_payments_20261017',
 'fresh_order16_finite_cross_payments_20261017',
 'fresh_finite8_ideal5_combined_payments_20261017',
 'fresh_finite8_cycle_replacement_20261017',
 'fresh_finite8_ideal5_local_replacement_20261017',
 'fresh_cycle14_ideal5_transfer_20261017',
 'fresh_affine_d13_long_pair_transfer_20261017',
 'fresh_affine_d14_long_pair_transfer_20261017',
 'fresh_finite_d17_d5_d12_unit_transfer_20261018',
 'fresh_cycle11_half_pair_transfer_20261018',
 'fresh_order11_bcd_half_pair_transfer_20261018',
 'fresh_order11_bc_zero_transfer_20261018',
 'fresh_order12_b_center_terminal_transfer_20261018',
 'fresh_order13_d_center_pair_half_transfer_20261018',
 'fresh_order13_d_center_pair_zero_transfer_20261018',
 'fresh_large_finite_d_leaf_half_transfer_20261017',
 'fresh_finite_rank15_existing_half_transfer_20261018',
 'fresh_finite_a15_endpoint_coupon_zero_20261018',
 'fresh_finite_b15_d15_existing_payment_transfer_20261018',
 'fresh_finite_d17_a4_d13_half_transfer_20261018',
 'fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018')
STEMS=tuple(dict.fromkeys((P.TAG,M.TAG,*P.STEMS,*M.STEMS,*PROVENANCE)))
FINAL_MD={P.TAG:'738c5ff6f3d2e3c28fb70c83179ac97e86d512059d60d4a71ed1b1cd6f38bb48'}


def sparse_debit(a,b):
    return sum(12*s[1] for s,t in ((a,b),(b,a)) if 2<=s[1]<=5 and t==('D',16))


def actual_base(t):return P.actual_base(t)
def upper_base(t):return P.upper_base(t)
def actual_cross(a,b):
    return P.actual_cross(a,b)-sparse_debit(a,b)-(2 if {a,b}=={('D',5),('D',12)} else 0)
def upper_cross(a,b):
    return P.upper_cross(a,b)-sparse_debit(a,b)-(2 if {a,b}=={('A',5),('D',12)} else 0)
def item(r,key):
    return P.item(r,key[:7])-(12*r*key[7] if 2<=r<=5 else 0)-(2*(key[1]-key[3]) if r==5 else 0)


def value(forest,upper=False):
    base,cross=(upper_base,upper_cross) if upper else (actual_base,actual_cross)
    return sum(base(t) for t in forest)+sum(cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])


def table(key):
    weights={r:item(r,key) for r in Z.ITEMS};v=[0]+[-10**9]*N
    for r,w in weights.items():
        for s in range(r,N+1):v[s]=max(v[s],v[s-r]+w)
    last=[0]+[next(r for r,w in weights.items() if r<=s and v[s]==v[s-r]+w) for s in range(1,N+1)]
    return dict(key=list(key),values=v,last=last)


def check_table(t,key):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==61 and v[0]==last[0]==0
    weights={r:item(r,key) for r in Z.ITEMS};count=0
    for s in range(1,61):
        for r,w in weights.items():
            if r<=s:assert v[s]>=v[s-r]+w;count+=1
        r=last[s];assert r in weights and r<=s and v[s]==v[s-r]+weights[r]
    return count


def finite(saved=None):
    rep=Z.B.B.representative;types=[t for t in Z.B.B.types() if t[1]<=60]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)
        for b in types:
            assert 0<=actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)
    assert actual_cross(('D',5),('D',12))==upper_cross(('A',5),('D',12))==90
    for r in range(2,6):assert actual_cross(('A',r),('D',16))==upper_cross(('A',r),('D',16))==0
    assert actual_cross(('A',1),('D',16))==12
    dom=[(cs,used,(*key,cs[7])) for cs,used,key in Z.large_domain()];assert len(dom)==1044
    keys=sorted({key for cs,used,key in dom})
    tabs={tuple(t['key']):t for t in saved['tables']} if saved else {k:table(k) for k in keys}
    assert len(keys)==1042 and set(tabs)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tabs[k],k) for k in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];items=0
    for cs,used,key in dom:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in fixed)==used and key[7]==cs[7]
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.ITEMS,residual):
            assert item(r,key)==upper_base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=N-used;tab=tabs[key];v=sum(upper_base(t) for t in fixed)+tab['values'][rem];forest=fixed[:]
        while rem:
            r=tab['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N and value(forest,True)==v and value(forest)<=v
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
                         octupled_upper=v,octupled_actual=value(forest)))
    maximum=max(r['octupled_upper'] for r in rows)
    witness=[('D',17),('D',17),('D',26)];assert maximum==3260==value(witness)
    first=[('D',16)]*3+[('A',5)]*2+[('A',2)]
    second=[('A',1)]*4+[('D',5)]*4+[('D',12)]*3
    assert P.actual_value(first)==3258 and value(first)==2826
    assert P.actual_value(second)==3248 and value(second)==3224
    assert sum(r['octupled_upper']>3240 for r in rows)==6
    return dict(exact_maximum='815/2',curvature='405',finite_deficit='5/2',complete_fixed_vectors=len(dom),
        interaction_keys=len(keys),simultaneous_type_pairs=len(types)**2,residual_item_checks=items,
        Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=[tabs[k] for k in keys],alternatives=rows,above_curvature=sum(r['octupled_upper']>3240 for r in rows),
        actual_attainer=['D17','D17','D26'],D5_D12_actual_cost='45/4',
        small_D16_obstruction=dict(forest=['D16']*3+['A5']*2+['A2'],old='1629/4',new='1413/4'),
        D5_D12_obstruction=dict(forest=['A1']*4+['D5']*4+['D12']*3,old='406',new='403'),
        finite_cost_unit='1/8',new_residual_tables=True,replay_optimizes=False,
        actual_rank15_costs=dict(A15='101',B15='97',D15='201/2'),rank15_representative_upper='101',
        formal_forest_attainment_not_geometric_realization=True)


def seals(docs):
    old={r['index']:r for r in docs[M.O.TAG]['cases']};over=docs[M.O.OVER];proofs=[]
    weights={tuple(r['pair']):F(r['weight']) for r in docs[B.TAG]['local_proofs']['ideal_weighted_pairs']}
    assert sum(weights.values())==114
    for pair in M.CHANGE:assert weights[pair]==0;weights[pair]=F(1,8)
    assert sum(weights.values())==F(229,2)
    for frozen in docs[M.TAG]['selected_proofs']:
        i=frozen['index'];c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        proof=M.inspect(c,core,old[i]);assert proof=={k:v for k,v in frozen.items() if k!='maps'}
        corners=[]
        for corner in proof['corners']:
            ip=F(1,8) if corner['kind']=='affineD16_order17' else F()
            fp=F(7,8) if F(corner['finite_coupon']) else F()
            corners.append(dict(marks=corner['marks'],ideal_payment=ip,finite_coupon=fp))
        assert sum(c['ideal_payment']+c['finite_coupon'] for c in corners)==1
        proofs.append(dict(index=i,old_seal=proof,corners=corners,total='1'))
    assert [p['index'] for p in proofs]==M.EXPECTED
    assert docs[M.TAG]['complete_eligible_patterns']==140 and docs[M.TAG]['affected_maps']==32
    assert sum(len(r['embeddings']) for r in docs[M.TAG]['reconstructed_support'])==488
    role=docs[B.H.TAG]['disjoint_coupon_roles']
    assert role['this_packet']==[[0,9]] and role['preceding_D114']==[[0,10],[0,11]]
    assert role['finite_D17']==[[3,9],[3,10],[3,11],[4,9],[4,10],[4,11]]
    return dict(proofs=proofs,complete_inherited_selected_patterns=140,complete_inherited_maps=488,affected_maps=32,
        ideal_changed_pairs=sorted(M.CHANGE),ideal_new_weight='1/8',ideal_new_cost='229/2',
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        finite_coupon=dict(small='whole_D5',small_role='unique_long_arm_endpoint',small_marks=[0],
            large='whole_D12',large_role='short_fork_leaf',large_marks=[10,11],old_weight='1',new_weight='7/8',
            component_pair_delta='-1/4'),disjoint_coupon_roles=role,
        historical_map_search_repeated=False,triangle_minimum=2,quadrilateral_minimum=1)


def literal(docs):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in B.C.R.B.C.L.STEMS}
    prev={(r['source'],B.C.R.B.C.L.N.normalize(r['case'])):r for r in docs[P.TAG]['combined_literal_interface']['rows']}
    target=docs[M.TAG]['target'];rows=[];changed=[];d16=[];credits=debits=sparse=0
    for stem,index,g,white,cs in T.literal_inputs(inputs):
        before=prev[stem,B.C.R.B.C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));ip=cp=sp=F();images=[];roles=[]
            sig=Q.inertia(A.sub(g,ids))
            if sig==[len(ids),0,0]:
                blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
                if aa is not bb:
                    for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                        if la['type']=='D16':
                            d16.append(dict(source=stem,case=index,marks=marks,small=sm,large=la))
                            assert sm['type']=='A1'
                            if 2<=len(sm['nodes'])<=5 and la['roles'].get(str(lm)) in ('fork_center','short_fork_leaf'):sp=-F(1,2)
                        if sm['type']=='D5' and la['type']=='D12':
                            si=T.mark_pairs(g,sm['nodes'],[m],5);li=T.mark_pairs(g,la['nodes'],[lm],12)
                            if si==[(0,)] and set(li)=={(10,),(11,)}:
                                cp=-F(1,8);roles.append(dict(small_images=si,large_images=li))
            elif sig[1]==0:
                shared=next((p for p in A.components(g,ids) if set(marks)<=set(p)),None)
                if shared and len(shared)==17:
                    graph=T.T.labels(g,shared)
                    if graph is not None:
                        images=sorted({tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in M.P.embeddings(graph,target)})
                        if set(images)&M.CHANGE:assert set(images)<=M.CHANGE;ip=F(1,8)
            v=F(before['corners'][j]['retained_payment_lower_bound'])+ip+cp+sp;assert v>=0
            credits+=ip>0;debits+=cp<0;sparse+=sp<0
            corners.append(dict(marks=marks,canonical_images=images,ideal_increment=ip,coupon_delta=cp,
                removed_D16_coupon_delta=sp,coupon_roles=roles,retained_payment_lower_bound=v))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['ideal_increment'] or c['coupon_delta'] or c['removed_D16_coupon_delta'] for c in corners):changed.append(row)
    assert len(rows)==216 and len(d16)==10 and sparse==0
    assert credits==debits==4
    return dict(configurations=216,ideal_credits=credits,coupon_debits=debits,removed_D16_coupon_occurrences=sparse,
        separate_D16_corners=d16,rows=rows,changed=changed,all_totals_at_least_one=True)


def ideal(previous):
    domain=sorted((low,)+p if low else p for low in range(9) for p in Z.high_parts(59-low))
    assert len(domain)==2113 and domain==sorted(tuple(r['ranks']) for r in previous['rows'])
    rows=[];limits=[]
    for old in previous['rows']:
        row=dict(old);count=row['ranks'].count(16);delta=F(count,34)
        row['normalized_cost']=str(F(row['normalized_cost'])+delta)
        row['normalized_margin']=str(F(row['normalized_margin'])-delta)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['long_short_rebalance_normalized_cost_increase']=str(delta);rows.append(row)
        if count:limits.append((F(old['normalized_margin'])*17/count,old['ranks']))
    summary=Z.summarize(rows)
    assert summary['negative_count']==summary['zero_count']==0
    assert summary['normalized_minimum']=='31/6000' and summary['absolute_minimum']=='54'
    assert min(limits)==(F(897,1360),[9,9,9,16,16])
    return dict(complete_profiles=2113,rows=rows,summary=summary,old_order17_cost='114',new_order17_cost='229/2',
        strict_order17_total_increase_threshold='897/1360',threshold_profile=[9,9,9,16,16],
        low9_concentration_only=True,positive_terminal_absolute_surplus_extends_to_all_actual_cusps=True,
        domain_independently_enumerated=True,convolutions_inherited_without_replay=True)


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,(s,path)
    for s,h in FINAL_MD.items():assert sha256((ROOT/f'results/{s}.md').read_bytes()).hexdigest()==h,s
    assert docs[P.TAG]['finite']['exact_maximum']=='815/2'
    assert docs[P.TAG]['ideal']['order17_maximum']=='114'
    paths={Path(__file__),Path(P.__file__),Path(M.__file__),Path(M.O.__file__),Path(T.__file__),Path(D.__file__),Path(Z.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return P.L.enc(dict(status='PASS_SMALL_D16_COUPON_IDEAL17_REBALANCE',sources=sources,
        completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        deleted_coupon=dict(small_whole_ranks=[2,3,4,5],large='whole_D16',large_roles=['fork_center','short_fork_leaf'],
            old_weight='1/2',new_weight='0',A1_D16_unchanged=True,all_other_cross_component_weights_unchanged_except_D5_D12=True),
        local_proofs=seals(docs),combined_literal_interface=literal(docs),
        finite=finite(saved['finite'] if saved else None),ideal=ideal(docs[P.TAG]['ideal']),
        triangle_minimum=2,quadrilateral_minimum=1,new_unrestricted_dimension_bound=None,
        necessary_ambient60_inequality=dict(finite_vertex_coefficient='5/2',constant='4500',ideal_vertex_coefficient='54'),
        all_ideal_dimension_bound=59,scope='Sparse removal of an unused broad D16 coupon enlargement with complete proof-interface inventory in the report, all216 current literal corners, four inherited radical seals, complete changed finite certificates and all2113 ideal row updates. Formal maximum407.5 remains above405; no59 bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    if not args.verify:assert not path.exists(),'Frozen certificate exists; use --verify.'
    out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'above',out['finite']['above_curvature'],
          'keys',out['finite']['interaction_keys'],'ideal',out['ideal']['summary'],flush=True)
