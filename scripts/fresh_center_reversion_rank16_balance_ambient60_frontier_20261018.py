#!/usr/bin/env python3
"""Restore ideal17 center halves, then rebalance sealed A/D16--17 corners.

Two exact complete forest certificates and one reconstructed ideal domain.
The historical geometry is inherited; eight existing seals are recomputed.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_order12_order16_depth8_ambient60_frontier_20261018 as C
import fresh_order17_d5_d12_center_quarter_transfer_20261018 as H
import fresh_a16_a17_full_payment_exchange_20261018 as A
import fresh_d16_d17_eighth_balance_20261018 as D

ROOT=C.ROOT;Z=C.Z;T=A.T;Q=A.Q;P=A.P;DIM=60
TAG='fresh_center_reversion_rank16_balance_ambient60_frontier_20261018'
STEMS=tuple(dict.fromkeys((C.TAG,H.TAG,A.TAG,D.TAG,H.M.TAG,*C.STEMS,*H.STEMS,*A.STEMS,*D.STEMS)))
FINAL_MD={C.TAG:'9ba5373780bcbf2cc907b99e12447144108011ab7806b59c037262ec292a9e3a',
 H.TAG:'ec94fd0c6f5ae5ba87d289bc6621b4361b308b5c1d887d8de129a7c8b194a5bb',
 A.TAG:'ced5478a43bd7245e29feb3006aa2cb36a3f684ff84c932435168b202f2669de',
 D.TAG:'3786fc72577499c9b1e8d8aea6d5d10f515906928d4de2b27f4147fab2b19e2d'}
SHIFT={('A',16):2,('A',17):-2,('D',16):4,('D',17):-4}
PAIRS={'A16':A.ZERO,'A17':A.RESTORE,'D16':D.UP,'D17':D.DOWN}
DELTA={'A16':F(1,8),'A17':F(-1,8),'D16':F(1,4),'D17':F(-1,4)}
NEW_WEIGHT={'A16':F(1,8),'A17':F(7,8),'D16':F(7,8),'D17':F(1,8)}


def base(t,stage):return C.base(t)+(SHIFT.get(t,0) if stage else 0)
def actual_cross(a,b):return C.actual_cross(a,b)-(2 if {a,b}=={('D',5),('D',12)} else 0)
def upper_cross(a,b):return C.upper_cross(a,b)-(2 if {a,b}=={('A',5),('D',12)} else 0)
def item(r,key,stage):return C.item(r,key)-(2*(key[1]-key[3]) if r==5 else 0)-(4 if stage and r==17 else 0)


def table(key,stage):
    values=[0]+[-10**9]*DIM
    for r in Z.Z.ITEMS:
        w=item(r,key,stage)
        for n in range(r,DIM+1):values[n]=max(values[n],values[n-r]+w)
    last=[0]+[next(r for r in Z.Z.ITEMS if r<=n and values[n]==values[n-r]+item(r,key,stage)) for n in range(1,DIM+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(t,key,stage):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==61 and v[0]==last[0]==0
    count=0
    for n in range(1,61):
        for r in Z.Z.ITEMS:
            if r<=n:assert v[n]>=v[n-r]+item(r,key,stage);count+=1
        r=last[n];assert r in Z.Z.ITEMS and r<=n and v[n]==v[n-r]+item(r,key,stage)
    return count


def actual_value(forest,stage):
    return sum(base(t,stage) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])


def finite(stage,saved=None):
    rep=Z.Z.B.B.representative;types=[t for t in Z.Z.B.B.types() if t[1]<=DIM]
    for a in types:
        assert base(rep(a),stage)>=base(a,stage)
        for b in types:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)
    assert actual_cross(('D',5),('D',12))==upper_cross(('A',5),('D',12))==92
    domain=Z.Z.large_domain();keys=sorted({key for cs,used,key in domain})
    assert len(domain)==1044 and len(keys)==792
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {key:table(key,stage) for key in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tables[key],key,stage) for key in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in Z.Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];items=0
    for cs,used,key in domain:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in fixed)==used and key[1]-key[3]==cs[3]
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.Z.ITEMS,residual):
            assert item(r,key,stage)==base(t,stage)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=DIM-used;t=tables[key];value=sum(base(x,stage) for x in fixed)+t['values'][rem]
        forest=fixed[:]
        while rem:
            r=t['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==DIM
        upper=sum(base(t,stage) for t in forest)+sum(upper_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert upper==value and actual_value(forest,stage)<=upper
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            octupled_upper=value,octupled_actual=actual_value(forest,stage)))
    maximum=max(r['octupled_upper'] for r in rows);assert maximum==(3264 if stage else 3268)
    witness=[('D',15)]*4 if stage else [('D',17),('D',17),('D',26)]
    assert sum(r for k,r in witness)==60 and actual_value(witness,stage)==maximum
    former=[('A',1)]*4+[('D',5)]*4+[('D',12)]*3;assert actual_value(former,stage)==3248
    return dict(stage='center_reversion_and_balance' if stage else 'center_reversion_only',
        exact_maximum=str(F(maximum,8)),curvature='405',finite_deficit=str(F(maximum,8)-405),
        complete_fixed_vectors=1044,interaction_keys=792,simultaneous_type_pairs=len(types)**2,
        residual_item_checks=items,Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=[tables[k] for k in keys],alternatives=rows,above_curvature=sum(r['octupled_upper']>3240 for r in rows),
        actual_attainer=[k+str(r) for k,r in witness],actual_attaining_cost=str(F(maximum,8)),
        former409_forest_cost='406',actual_D5_D12_cost='23/2',representative_A5_D12_allowance='23/2',
        costs={k+str(r):str(F(base((k,r),stage),8)) for k,r in SHIFT},
        finite_cost_unit='1/8',new_tables_generated=True,replay_optimizes=False,
        formal_forest_attainment_not_geometric_realization=True)


def local_proofs(docs):
    ideal={tuple(r['pair']):F(r['weight']) for r in docs[H.TAG]['weighted_pairs']}
    for p in H.CHANGE:assert ideal[p]==F(1,4);ideal[p]=F(1,2)
    assert sum(ideal.values())==114
    assert ideal=={tuple(r['pair']):F(r['weight']) for r in docs[H.M.TAG]['weighted_pairs']}
    role=docs[H.TAG]['disjoint_coupon_roles'];assert role['this_packet']==[[0,9]]
    assert role['preceding_D114']==[[0,10],[0,11]]
    assert role['finite_D17']==[[3,9],[3,10],[3,11],[4,9],[4,10],[4,11]]
    over=docs[H.O.OVER];old={r['index']:r for r in docs[H.O.TAG]['cases']};centers=[]
    for saved in docs[H.TAG]['selected_proofs']:
        i=saved['index'];c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        proof=H.inspect(c,core,old[i]);assert proof=={k:v for k,v in saved.items() if k!='maps'}
        corners=[]
        for r in proof['corners']:
            ip=F(1,2) if F(r['ideal_payment']) else F();fp=F(1,2) if F(r['finite_coupon']) else F()
            corners.append(dict(marks=r['marks'],ideal_payment=ip,finite_coupon=fp))
        assert sum(r['ideal_payment']+r['finite_coupon'] for r in corners)==1
        centers.append(dict(index=i,old_seal=proof,corners=corners,total='1'))
    assert [r['index'] for r in centers]==[1314,2529]
    weights={}
    for typ,stem in (('A16',A.TAG),('A17',A.TAG),('D16',D.TAG),('D17',D.TAG)):
        w={tuple(r['pair']):F(r['weight']) for r in docs[stem][typ+'_weighted_pairs']}
        assert sum(w.values())==F(C.base((typ[0],int(typ[1:]))),8)
        for p in PAIRS[typ]:assert w[p]+DELTA[typ]==NEW_WEIGHT[typ];w[p]=NEW_WEIGHT[typ]
        assert sum(w.values())==F(base((typ[0],int(typ[1:])),1),8)
        weights[typ]=w
    assert weights['A16'][1,13]==weights['A16'][2,14]==F(1,2) and weights['A17'][1,15]==1
    assert weights['D16'][1,11]==0 and weights['D16'][1,14]==weights['D16'][1,15]==F(1,2)
    # Complete saved A16 native support, including empty cases.
    sup=A.V.support('A',docs[A.T.TAG]);target=docs[A.A16]['target'];native=docs[A.A.OVER]
    assert A.V.validate_maps(sup,target,native)==238 and sup['complete_eligible_patterns']==394
    hits=[r['index'] for r in sup['rows'] if any(tuple(m['marked_images']) in A.ZERO for m in r['maps'])]
    assert hits==[1657,2067]
    ap=[]
    for old in docs[A.TAG]['selected_proofs']:
        i=old['index'];c=native['isolation']['cases'][i];core=native['isolation']['cores'][c['core_key']]
        proof=A.V.core_proof(c,core,'A',A.ZERO);assert A.V.enc(proof)==old['old_seal']
        corners=[]
        for r in proof['corners']:
            typ='A16' if r['reason']=='retained_selected_half' else 'A17' if r['reason']=='restored_whole_A17_half' else None
            value=NEW_WEIGHT[typ] if typ else F()
            if typ:assert set(r['canonical_images'])<=PAIRS[typ]
            else:assert r['payment']==0
            corners.append(dict(marks=r['marks'],type=typ,canonical_images=r['canonical_images'],payment=value))
        assert sum(r['payment'] for r in corners)==1
        ap.append(dict(index=i,old_seal=proof,corners=corners,total='1'))
    # The reduced D17 orbit has no original selected map; its only later
    # dependency is the four whole D16/D17 seals reconstructed below.
    target=docs[D.SUPPORT]['target'];sup=docs[D.SUPPORT]['complete_selected_support'];maps=0;hits=[]
    assert sup['complete_eligible_patterns']==140
    for row in sup['rows']:
        c=native['isolation']['cases'][row['index']];core=native['isolation']['cores'][c['core_key']]
        g=A.A.matrix(core['gram']);ns=c['marked_component'];graph=T.T.labels(g,ns);seen=set();pairs=set()
        for saved in row['maps']:
            mp=saved['images'];assert len(mp)==len(ns)==len(set(mp)) and tuple(mp) not in seen;seen.add(tuple(mp));maps+=1
            assert all(0<=x<17 for x in mp)
            assert all(graph[i][j]==target[mp[i]][mp[j]] for i in range(len(mp)) for j in range(len(mp)))
            pair=tuple(sorted(mp[ns.index(m)] for m in c['preferred_corner']));assert list(pair)==saved['marked_images'];pairs.add(pair)
            if pair in D.DOWN:hits.append(row['index'])
        assert sorted(pairs)==list(map(tuple,row['marked_images']))
    assert maps==168 and not hits
    dp=[]
    for old in docs[D.TAG]['reweighted_proofs']:
        i=old['index'];c=native['isolation']['cases'][i];core=native['isolation']['cores'][c['core_key']]
        proof=D.J.seal(c,core);assert D.J.enc(proof)==old['old_seal'];corners=[]
        for r in proof['corners']:
            typ='D16' if r['reason']=='selected_D16_half' else 'D17' if r['reason']=='restored_whole_D17_half' else None
            value=NEW_WEIGHT[typ] if typ else F()
            if typ:assert set(r['canonical_images'])==PAIRS[typ]
            else:assert r['payment']==0
            corners.append(dict(marks=r['marks'],type=typ,canonical_images=r['canonical_images'],payment=value))
        assert sum(r['payment'] for r in corners)==1
        dp.append(dict(index=i,old_seal=proof,corners=corners,total='1'))
    assert [r['index'] for r in dp]==[1644,1664,1931,2324]
    return dict(center_seals=centers,A_exchange_seals=ap,D_exchange_seals=dp,
        complete_A16_support=dict(patterns=394,maps=238,affected_contexts=[1657,2067]),
        complete_D17_support=dict(patterns=140,maps=168,reduced_pair_occurrences=0),
        order17_family_costs=dict(A='85',B='111',C='108',D='114'),ideal_pair_weight='1/2',
        ideal_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(ideal.items())],
        finite_weighted_pairs={k:[dict(pair=p,weight=v) for p,v in sorted(w.items())] for k,w in weights.items()},
        reduced_coupon=dict(small='whole_D5',small_role='unique_long_arm_endpoint',large='whole_D12',
            large_role='fork_center',old_weight='3/4',new_weight='1/2',eligible_marked_pairs=1,component_pair_change='-1/4'),
        disjoint_coupon_roles=role,triangle_minimum=2,quadrilateral_minimum=1,
        new_cross_component_coupons=[])


def literal(docs,stage):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in C.R.B.C.L.STEMS}
    previous={(r['source'],C.R.B.C.L.N.normalize(r['case'])):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    target=docs[H.TAG]['target'];rows=[];changed=[];counts=dict(ideal_credits=0,coupon_debits=0,finite_credits=0,finite_debits=0)
    for stem,index,g,white,cs in T.literal_inputs(inputs):
        before=previous[stem,C.R.B.C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));ip=F();cp=F();fp=F();images=[];typ=None;roles=[]
            sig=Q.inertia(A.A.sub(g,ids))
            if sig==[len(ids),0,0]:
                bs=A.D.finite_blocks(g,ids);aa,bb=[next(b for b in bs if m in b['nodes']) for m in marks]
                if aa is bb:
                    typ=aa['type']
                    if stage and typ in PAIRS:
                        images=T.mark_pairs(g,aa['nodes'],marks,len(aa['nodes']),typ[0])
                        if set(images)&PAIRS[typ]:assert set(images)<=PAIRS[typ];fp=DELTA[typ]
                else:
                    for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                        if sm['type']=='D5' and la['type']=='D12':
                            si=T.mark_pairs(g,sm['nodes'],[m],5);li=T.mark_pairs(g,la['nodes'],[lm],12)
                            if si==[(0,)] and li==[(9,)]:cp=F(-1,4);roles.append(dict(small_images=si,large_images=li))
            elif sig[1]==0:
                shared=next((p for p in A.A.components(g,ids) if set(marks)<=set(p)),None)
                if shared and len(shared)==17:
                    graph=T.T.labels(g,shared)
                    if graph is not None:
                        images=sorted({tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in P.embeddings(graph,target)})
                        if set(images)&H.CHANGE:assert set(images)<=H.CHANGE;ip=F(1,4)
            value=F(before['corners'][j]['retained_payment_lower_bound'])+ip+cp+fp;assert value>=0
            counts['ideal_credits']+=ip>0;counts['coupon_debits']+=cp<0;counts['finite_credits']+=fp>0;counts['finite_debits']+=fp<0
            corners.append(dict(marks=marks,type=typ,canonical_images=images,ideal_increment=ip,
                coupon_delta=cp,coupon_roles=roles,finite_delta=fp,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['ideal_increment'] or c['coupon_delta'] or c['finite_delta'] for c in corners):changed.append(row)
    assert len(rows)==216 and counts==dict(ideal_credits=2,coupon_debits=2,finite_credits=4 if stage else 0,finite_debits=0)
    return dict(configurations=216,counts=counts,changed=changed,rows=rows,all_totals_at_least_one=True)


def ideal(docs):
    previous=C.ideal(docs);assert previous==docs[C.TAG]['ideal'];rows=[]
    for old in previous['rows']:
        row=dict(old);penalty=F(row['ranks'].count(16),34)
        row['normalized_cost']=str(F(row['normalized_cost'])+penalty)
        row['normalized_margin']=str(F(row['normalized_margin'])-penalty)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['center_reversion_normalized_cost_increase']=str(penalty);rows.append(row)
    summary=Z.Z.summarize(rows)
    assert len(rows)==2113 and summary['negative_count']==summary['zero_count']==0
    assert summary['normalized_minimum']=='31/6000' and summary['absolute_minimum']=='54'
    return dict(complete_profiles=2113,rows=rows,summary=summary,low9_concentration_only=True,
        positive_terminal_absolute_surplus_extends_to_all_actual_cusps=True,order17_maximum='114')


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    for stem,h in FINAL_MD.items():assert sha256((ROOT/f'results/{stem}.md').read_bytes()).hexdigest()==h,stem
    assert docs[C.TAG]['finite']['exact_maximum']=='409' and docs[C.TAG]['completed_unrestricted_dimension_bound']==60
    paths={Path(__file__),Path(C.__file__),Path(H.__file__),Path(A.__file__),Path(D.__file__),
           Path(A.V.__file__),Path(D.J.__file__),Path(H.O.__file__),Path(T.__file__),Path(Z.Z.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    proofs=local_proofs(docs);ide=ideal(docs);stages=[]
    for stage in (0,1):
        stages.append(dict(finite=finite(stage,saved['stages'][stage]['finite'] if saved else None),literal_interface=literal(docs,stage)))
    return C.I16.enc(dict(status='PASS_CENTER_REVERSION_RANK16_BALANCE_AMBIENT60_FRONTIER',sources=sources,
        completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        local_proofs=proofs,ideal=ide,stages=stages,final_finite_maximum='408',final_finite_curvature='405',
        final_finite_costs=dict(A16='441/4',A17='467/4',D16='443/4',D17='467/4'),
        strict_global_constant='4500',constant_counted_once=True,
        necessary_ambient60_inequality=dict(finite_vertex_coefficient='3',constant='4500',ideal_vertex_coefficient='54'),
        all_ideal_dimension_bound=59,new_unrestricted_dimension_bound=None,new_finite_volume_example=False,
        scope='Compatible reversal to the proved ideal D114 center halves and two exact finite corner balances. Eight old whole-core seals reconstructed; complete saved native interfaces and216 current literals checked. Both1044/792 finite counts newly certified;2113 ideal convolutions reconstructed. Formal final408 remains above405; unrestricted60 unchanged.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if a.verify else None
    out=run(saved)
    if a.verify:assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',[s['finite']['exact_maximum'] for s in out['stages']],
          'ideal',out['ideal']['summary'],flush=True)
