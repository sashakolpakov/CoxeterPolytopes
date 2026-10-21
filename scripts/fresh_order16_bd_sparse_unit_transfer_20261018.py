#!/usr/bin/env python3
"""Order16 B102/D103 using forced terminals and a sparse existing allowance.

Complete native target domains, arbitrary-support completion proofs, literal
interfaces and pointwise finite coupon domination. No forest optimization.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_order16_affine_b_singleton_projection_transfer_20261018 as B
import fresh_order16_d_distance9_half_transfer_20261018 as H
import fresh_balanced_rank16_order14_ambient60_frontier_20261018 as C
import fresh_order13_bc_sparse_terminal_transfer_20261018 as BC

ROOT=B.ROOT;X=B.X;D=H.D;P=D.P;Q=D.Q;L=D.L;T=C.B.T;A=C.B.A
TAG='fresh_order16_bd_sparse_unit_transfer_20261018'
CHANGES={'B':{(0,10),(10,14),(10,15)},'D':{(2,14),(2,15),(9,12),(9,13)}}
B_EXPECT=[1019,1031,1036,1050,1266,1285,1614,1624,1647,1667,1906,1934,2299,2327]
D_EXPECT=[369,373,504,509,515,662,766,989,997,1031,1037,1051,1164,1177,1537,1614,
          1648,1668,1796,1830,1907,1912,1935,1947,1979,2035,2044,2189,2223,2300,
          2305,2328,2340,2372,2536,2540]
SEALED={1979,2372,2536,2540}
STEMS=(C.TAG,B.TAG,H.TAG,D.TAG,D.SMALL,D.OVER,X.TAG,BC.TAG,
       'fresh_order13_d_center_leaf_half_transfer_20261018',
       'fresh_order15_affine_b_residual_terminal_transfer_20261017',
       'fresh_distance9_supporting_extensions_review_20261017',
       'fresh_order11_bc_zero_transfer_20261018','fresh_depth10_finite7_zero_transfer_20261018',
       'fresh_finite7_ideal5_triangle_payments_20261017')


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def b_proof(c,core,maps,target):
    g,J,du,dv,ends=X.constraints(c,core);assert c['method']=='lorentz_residual'
    rows=[];base=c['marked_component']
    for mp in maps:
        image={t:base[i] for i,t in enumerate(mp)};entries=[]
        for new in sorted(set(range(16))-set(image)):
            for old,i in sorted(image.items()):
                label=target[new][old]
                if not label:continue
                assert label in (3,4) and i in c['possible_contacts']
                q2=F(1 if label==3 else 2)
                entries.append(dict(new_target=new,old_target=old,core_node=i,label=label,
                    physical_contact_square=q2,projection_lower_bound=q2*J[i][i]))
        bad=next((r for r in entries if r['projection_lower_bound']>2),None)
        if bad:
            assert bad['label']==4 and bad['projection_lower_bound'] in (F(4),F(8,3))
            rows.append(dict(map=mp,entries=entries,outcome='projection_above_norm_two',witness=bad));continue
        assert len(entries)==1
        e=entries[0];i=e['core_node']
        assert (e['new_target'],e['old_target'],e['label'])==(11,10,3)
        assert J[i][i]==2 and i in c['preferred_corner']
        assert 11 not in image and 12 not in image and target[11][12]==3
        same=next(end for end in c['end_sets'] if i in end)
        opposite=next(end for end in c['end_sets'] if i not in end)
        assert len(opposite)==1 and J[i][opposite[0]]==-2
        assert all(J[i][j]>0 for j in same)
        rows.append(dict(map=mp,entries=entries,outcome='saturated_entry_has_no_next_white',
            forced_old_contact=i,entry_target_node=11,next_target_node=12,next_label=3,
            residual_square=0,same_end=same,opposite_end=opposite,
            same_end_products={str(j):J[i][j] for j in same},opposite_product=-2,
            next_node_cannot_be_old_core=True))
    return dict(index=c['index'],gram=core['gram'],inverse=core['inverse'],
        complete_end_proof=ends,rows=rows)


def coupon(g,ids,marks):
    blocks=D.finite_blocks(g,ids);selected=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if selected[0] is selected[1]:return []
    out=[]
    for si,li in ((0,1),(1,0)):
        sm,la=selected[si],selected[li]
        if (sm['type'],la['type'])!=('D5','D11'):continue
        if la['roles'].get(str(marks[li]))!='short_fork_leaf':continue
        adj={i:[j for j in sm['nodes'] if i!=j and g[i][j]] for i in sm['nodes']}
        center=next(i for i in sm['nodes'] if len(adj[i])==3)
        leaves=[i for i in sm['nodes'] if len(adj[i])==1]
        long=[i for i in leaves if i not in adj[center]]
        assert len(long)==1
        if marks[si] not in long:continue
        assert all(not g[i][j] for i in sm['nodes'] for j in la['nodes'])
        out.append(dict(small='whole_D5',small_nodes=sm['nodes'],small_mark=marks[si],
            small_role='unique_long_arm_endpoint',large='whole_D11',large_nodes=la['nodes'],
            large_mark=marks[li],large_role='short_fork_leaf',old_weight='1/2',new_weight='1'))
    assert len(out)<=1
    return out


def d_proof(c,core,maps,target,old):
    g,con=D.core_constraints(c,core);rows=[]
    for mp in maps:
        r,completed=D.complete(c,core,g,con,target,mp)
        inherited=next(q for q in old['completions'] if any(z['map']==mp for z in q['embeddings']))
        # The predecessor merged identical literal Grams under target
        # automorphisms. Compare every field except the target ordering.
        ignored={'embeddings','target_to_literal','affine_nodes'}
        assert {k:v for k,v in enc(r).items() if k not in ignored}=={k:v for k,v in inherited.items() if k not in ignored}
        if 'affine_nodes' in r:assert sorted(r['affine_nodes'])==sorted(inherited['affine_nodes'])
        if r['outcome']=='sealed_completion':
            assert c['index'] in SEALED and not r['new_roots']
            assert r['signature']==[17,1,0] and len(completed)==18
            assert con['candidate_nodes']==[]
            corners=[]
            for cr in r['corners']:
                marks=cr['marks']
                cs=coupon(completed,sorted(set(r['white'])|set(marks)),marks) if cr['kind']=='finite' else []
                corners.append(dict(marks=marks,coupon=cs,payment=int(bool(cs))))
            assert sum(cr['payment'] for cr in corners)==1
            assert next(cr['marks'] for cr in corners if cr['payment'])!=c['preferred_corner']
            r['new_whole_corner_payments']=corners;r['new_total']='1'
            r['finite_same_component_units_used']=False
        else:assert c['index'] not in SEALED
        r['map']=mp;rows.append(r)
    return dict(index=c['index'],gram=core['gram'],inverse=core['inverse'],constraints=con,rows=rows)


def literal(current,targets):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in C.C.L.STEMS}
    before={(r['source'],C.C.L.N.normalize(r['case'])):r for r in current['combined_literal_interface']['rows']}
    rows=[];changed=[]
    for stem,index,g,white,cs in T.literal_inputs(inputs):
        old=before[stem,C.C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();inc=F();images={};coupons=[]
            sig=Q.inertia(A.sub(g,ids))
            if sig==[len(ids),0,0]:
                coupons=coupon(g,ids,marks);inc=F(len(coupons),2)
            elif sig[1]==0:
                shared=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if shared and len(shared)==16:
                    labels=T.T.labels(g,shared)
                    if labels is not None:
                        for kind,target in targets.items():
                            pairs=sorted({tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
                            images[kind]=pairs
                            if set(pairs)&CHANGES[kind]:
                                assert set(pairs)<=CHANGES[kind];debit=F(1) if kind=='B' else F(1,2)
            value=F(old['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert value>=0
            corners.append(dict(marks=marks,ideal_images=images,new_ideal_debit=debit,
                new_coupon_increment=inc,coupons=coupons,retained_payment_lower_bound=value))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(r['new_ideal_debit'] or r['new_coupon_increment'] for r in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changed,rows=rows)


def finite_comparison():
    rep=C.Z.Z.B.B.representative;types=[t for t in C.Z.Z.B.B.types() if t[1]<=60]
    comparisons=0
    for a in types:
        for b in types:
            increment=8 if {a,b}=={('D',5),('D',11)} else 0
            value=C.actual_cross(a,b)+increment
            assert value<=C.upper_cross(rep(a),rep(b));comparisons+=1
    assert C.actual_cross(('D',5),('D',11))==60
    assert C.upper_cross(('A',5),('D',11))==76
    return dict(complete_typed_comparisons=comparisons,cost_unit='1/8',
        old_actual_D5_D11='15/2',new_actual_D5_D11='17/2',unchanged_representative_A5_D11='19/2',
        residual_item_weights_unchanged=True,other_new_A4_B9_increment_not_counted_here=True,
        no_dynamic_program_or_new_forest_certificate=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,doc in docs.items():
        assert doc['status'].startswith('PASS'),s
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[C.TAG]['completed_unrestricted_dimension_bound']==60
    assert docs[BC.TAG]['simultaneous_order13_maximum']=='145/2'
    assert not docs[BC.TAG]['literal_interface']['changed']
    small=docs[D.SMALL];over=docs[D.OVER]['isolation'];records=[];targets={}
    patterns=[p for p in small['patterns'] if p['order_bound']>=16 and len(p['nodes'])<=16]
    assert len(patterns)==394
    oldd={r['index']:r for r in docs[D.TAG]['cases']}
    for kind,previous in (('B',docs[B.TAG]),('D',docs[H.TAG])):
        rec=next(r for r in small['records'] if (r['kind'],r['order'])==(kind,16));target=rec['target'];targets[kind]=target
        weights={tuple(r['pair']):F(r['weight']) for r in previous['weighted_pairs']} if kind=='D' else {p:F(int(p in set(map(tuple,previous['retained_unit_pairs'])))) for p in combinations(range(16),2)}
        assert sum(weights.values())==105
        for p in CHANGES[kind]:assert weights[p]==(1 if kind=='B' else F(1,2));weights[p]=F()
        assert sum(weights.values())==(102 if kind=='B' else 103)
        _,dist,features=P.distance_data(target);autos=P.embeddings(target,target)
        assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES[kind]}==CHANGES[kind] for mp in autos)
        assert all(dist[i][j]>=10 for i,j in CHANGES[kind])
        basic=L.S.basic_pairs(target);assert not basic&CHANGES[kind]
        assert all(weights[p]==1 for p in weights if dist[p[0]][p[1]]<=8)
        triangle={p for p in weights if dist[p[0]][p[1]]<=5 or (dist[p[0]][p[1]]<=8 and any(dist[i][f]<=5 for i in p for f in features))}
        assert all(weights[p]==1 for p in triangle)
        saved={r['case']:r for r in rec['case_embeddings']};support=[];proofs=[];counts=Counter()
        for p in patterns:
            maps=P.embeddings(p['labels'],target);images=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
            if maps:assert maps==saved[p['index']]['embeddings'] and images==list(map(tuple,saved[p['index']]['marked_images']))
            else:assert p['index'] not in saved
            support.append(dict(index=p['index'],map_count=len(maps),marked_images=images))
            chosen=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in CHANGES[kind]]
            if not chosen:continue
            c=over['cases'][p['index']];core=over['cores'][c['core_key']]
            proof=b_proof(c,core,chosen,target) if kind=='B' else d_proof(c,core,chosen,target,oldd[c['index']])
            proofs.append(proof);counts.update(r['outcome'] for r in proof['rows'])
        assert [p['index'] for p in proofs]==(B_EXPECT if kind=='B' else D_EXPECT)
        assert {r['index'] for r in support if r['map_count']}==set(saved)
        assert counts==(dict(projection_above_norm_two=14,saturated_entry_has_no_next_white=12) if kind=='B' else dict(excluded_projection_above_norm=64,excluded_multiple_core_contact_roots=56,excluded_lorentz_index=48,sealed_completion=32))
        records.append(dict(kind=kind,target=target,changed_pairs=sorted(CHANGES[kind]),old_cost='105',new_cost=str(sum(weights.values())),
            weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
            complete_patterns=394,successful_patterns=len(saved),complete_maps=sum(r['map_count'] for r in support),
            support=support,cases=proofs,affected_maps=sum(counts.values()),outcomes=dict(counts),retained_basic_pairs=sorted(basic),
            retained_triangle_pairs=sorted(triangle),automorphisms=autos))
    routing=Counter()
    for c in over['cases']:
        if c['method']=='finite_leaf9_or_ideal_feature5':routing['protected']+=1
        elif c['method']=='nine_context_adaptive_finite_payment':routing['adaptive']+=1
        elif c['component_order_bound']<16:routing['whole_bound_below16']+=1
        elif len(c['marked_component'])>16:routing['displayed_above16']+=1
        else:routing['eligible']+=1
    assert routing==dict(protected=74,adaptive=9,whole_bound_below16=2045,displayed_above16=92,eligible=394)
    paths={Path(__file__),Path(B.__file__),Path(H.__file__),Path(D.__file__),Path(X.__file__),Path(C.__file__),Path(P.__file__),Path(Q.__file__),Path(L.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_ORDER16_B102_D103_SPARSE_UNIT_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},records=records,
        complete_high_overlap_routing=dict(routing),new_order16_family_costs=dict(A=80,B=102,C=103,D=103),
        new_order16_maximum='103',new_normalized_rank15_envelope='103/16',
        upgraded_coupon=dict(small='whole_D5',small_role='unique_long_arm_endpoint',small_marks=[0],
            large='whole_D11',large_role='short_fork_leaf',large_marks=[9,10],old_weight='1/2',new_weight='1',
            eligible_marked_pairs=2,component_pair_cost_increase='1',old_component_pair_total='15/2',new_component_pair_total='17/2'),
        literal_interface=literal(docs[C.TAG],targets),finite_pointwise_comparison=finite_comparison(),
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='Complete394 native target searches per family;26 B forced-entry/saturation exclusions and200 D completions, including32 whole-factor maps with distinct D5-long-end/D11-shortleaf units. All current literal interfaces and finite pointwise allowance checked; no new forest, global scalar or ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_order16_family_costs'],'literal changes',len(out['literal_interface']['changed']),flush=True)
