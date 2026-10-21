#!/usr/bin/env python3
"""B12/order13 cost145/2 and C12/order13 cost72, with a sparse B9 coupon.

Complete necessary native maps, all74 unbounded protected bases and216
current literal interfaces. No finite budget or ambient bound is claimed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_order13_affine_b_singleton_transfer_20261017 as M
import fresh_b14_b15_singleton_completion_transfer_20261017 as X
import fresh_d17_fractional_ambient60_frontier_20261018 as C

ROOT=M.ROOT;P=M.P;Q=M.Q;T=C.L.T;A=C.L.A
TAG='fresh_order13_bc_sparse_terminal_transfer_20261018'
OLD_B='fresh_order13_affine_b_fork_projection_transfer_20261017'
OLD_C='fresh_ideal5_weighted_dimension71_20261017'
D13='fresh_order13_d_center_leaf_half_transfer_20261018'
STEMS=(C.TAG,M.DIAG,M.OVER,M.TAG,OLD_B,OLD_C,
       D13,
       'fresh_order14_c_distance10_terminal_projection_20261018',
       'fresh_order11_bc_zero_transfer_20261018','fresh_depth10_finite7_zero_transfer_20261018',
       'fresh_finite8_ideal5_local_replacement_20261017','fresh_finite8_cycle_replacement_20261017')
PAIRS={'B':(1,10),'C':(1,11)}
EXPECTED={'B':[157,174,952,2474],'C':[625,1556,1642]}


def proper(c,core,maps,target,kind):
    g=A.matrix(core['gram']);n=len(g);gauge=list(map(F,core['gauge']));J=Q.inverse(g)
    assert all(g[i][i]==2*gauge[i] for i in range(n)) and Q.inertia(g)==[n-1,1,0]
    assert M.encode(J)==core['inverse']
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]];black=u+v
    white=sorted(set(range(n))-set(black));assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,white+[b]))==[len(white)+1,0,0] for b in black)
    actual=Q.inspect_corner(g,gauge,black,*c['preferred_corner'])
    for k,val in actual.items():assert c[k]==val,(c['index'],k)
    ends=c['end_sets'];allowed=c['possible_contacts'];du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v)
    assert set(ends[0]).isdisjoint(ends[1]) and set(ends[0]+ends[1])==set(allowed)
    assert all(J[i][j]>0 for end in ends for i in end for j in end)
    assert all(min(du[i]+dv[j],dv[i]+du[j])+2<=8 for i in ends[0] for j in ends[1])
    base=c['marked_component'];other=sorted(set(white+c['preferred_corner'])-set(base))
    assert all(g[i][j]==0 for i in base for j in other)
    rows=[]
    for mp in maps:
        image={t:base[i] for i,t in enumerate(mp)};entries=[]
        for new in sorted(set(range(13))-set(image)):
            for old,node in image.items():
                bond=target[new][old]
                if not bond:continue
                assert bond in (3,4) and node in allowed
                physical=F(1 if bond==3 else 2);lower=physical*gauge[node]*J[node][node]
                entries.append(dict(new_target_node=new,old_target_node=old,core_node=node,native_label=bond,
                    physical_native_contact_square=physical,old_gauge_squared=gauge[node],
                    inverse_diagonal=J[node][node],projection_lower_bound=lower))
        assert entries and all(e['projection_lower_bound']>2 for e in entries)
        if kind=='B':assert len(entries)==1 and entries[0]['native_label']==3 and entries[0]['projection_lower_bound']==14
        else:assert len(entries)==2 and all(e['native_label']==4 for e in entries)
        rows.append(dict(map=mp,displayed_target_nodes=sorted(image),forced_entries=entries))
    return dict(index=c['index'],kind=kind,method=c['method'],core_key=c['core_key'],gram=core['gram'],gauge_squared=gauge,
        inverse=core['inverse'],allowed_contacts=allowed,end_sets=ends,other_old_corner_nodes=other,maps=rows)


def coupon(g,ids,marks):
    gauge=[g[i][i]/2 for i in range(len(g))];bs=M.blocks(g,ids,gauge)
    a,b=[next(b for b in bs if m in b['nodes']) for m in marks];out=[]
    if a is b:return out
    for small,large,m,lm in ((a,b,*marks),(b,a,marks[1],marks[0])):
        if small['type']=='A4' and large['type']=='B9' and small['roles'].get(str(m))=='endpoint' and large['roles'].get(str(lm))=='doubled_terminal':
            out.append(dict(small='whole_A4',small_nodes=small['nodes'],small_mark=m,small_role='endpoint',
                large='whole_B9',large_nodes=large['nodes'],large_mark=lm,large_role='doubled_terminal',payment='1/2'))
    assert len(out)<=1
    return out


def seal(c,core):
    p=M.inspect_case(c,core);g=A.matrix(core['gram']);corners=[]
    assert c['index'] in (174,2474)
    for old in p['corners']:
        marks=old['marks'];ids=sorted(set(p['white'])|set(marks))
        if old['kind']=='affineB12_order13':
            assert marks==c['preferred_corner'];pay=F(1,2);cs=[]
        else:cs=coupon(g,ids,marks);pay=F(len(cs),2)
        corners.append(dict(marks=marks,kind=old['kind'],coupon=cs,payment=pay))
    assert sum(r['payment'] for r in corners)==1 and sum(bool(r['coupon']) for r in corners)==1
    return dict(index=c['index'],inherited_full_core_seal=p,corners=corners,total='1')


def literal(current,targets):
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in C.L.STEMS}
    previous={(r['source'],C.L.N.normalize(r['case'])):r for r in current['combined_literal_interface']['rows']}
    rows=[];changed=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        old=previous[stem,C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();inc=F();images={};cp=[]
            sig=Q.inertia(A.sub(g,ids))
            if sig==[len(ids),0,0]:cp=coupon(g,ids,marks);inc=F(len(cp),2)
            elif sig[1]==0:
                shared=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if shared and len(shared)==13:
                    labels=T.T.labels(g,shared)
                    if labels is not None:
                        for kind,target in targets.items():
                            pairs=sorted({tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
                            images[kind]=pairs
                            if PAIRS[kind] in pairs:
                                assert pairs==[PAIRS[kind]];debit=F(1,2) if kind=='B' else F(1)
            value=F(old['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert value>=0
            corners.append(dict(marks=marks,ideal_images=images,new_debit=debit,new_coupon_increment=inc,coupons=cp,retained_payment_lower_bound=value))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(r['new_debit'] or r['new_coupon_increment'] for r in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changed,rows=rows)


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def run():
    data={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,d in data.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    diag=data[M.DIAG];over=data[M.OVER]['isolation'];assert len(diag['patterns'])==diag['eligible_patterns']==801
    records=[];targets={}
    for kind in ('B','C'):
        r=next(r for r in diag['records'] if r['kind']==kind);target=r['target'];targets[kind]=target;pair=PAIRS[kind]
        before=set(map(tuple,data[OLD_B]['retained_unit_pairs'])) if kind=='B' else set(map(tuple,next(r for r in data[OLD_C]['order13_interfaces']['records'] if r['kind']=='C')['retained_pairs']))
        assert len(before)==73 and pair in before
        weights={p:F(int(p in before)) for p in combinations(range(13),2)};weights[pair]=F(1,2) if kind=='B' else F()
        _,dist,features=P.distance_data(target);assert dist[pair[0]][pair[1]]==(9 if kind=='B' else 10)
        assert all(tuple(sorted(m[i] for i in pair))==pair for m in P.embeddings(target,target))
        assert all(weights[p]==1 for p in weights if dist[p[0]][p[1]]<=8)
        assert pair not in set(map(tuple,r['immediate_pairs']))
        triangle={p for p in weights if dist[p[0]][p[1]]<=5 or (dist[p[0]][p[1]]<=8 and any(dist[i][f]<=5 for i in p for f in features))}
        assert all(weights[p]==1 for p in triangle)
        saved={r['case']:r for r in r['case_embeddings']};support=[];proofs=[];changed_maps=0
        for p in diag['patterns']:
            maps=P.embeddings(p['labels'],target);pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
            if maps:assert maps==saved[p['index']]['embeddings'] and pairs==list(map(tuple,saved[p['index']]['marked_images']))
            else:assert p['index'] not in saved
            support.append(dict(index=p['index'],map_count=len(maps),marked_images=pairs))
            chosen=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks']))==pair]
            if not chosen:continue
            c=over['cases'][p['index']];core=over['cores'][c['core_key']];changed_maps+=len(chosen)
            result=seal(c,core) if c['method']=='radical' else proper(c,core,chosen,target,kind)
            proofs.append(dict(index=c['index'],maps=chosen,proof=result))
        assert [r['index'] for r in proofs]==EXPECTED[kind]
        assert {r['index'] for r in support if r['map_count']}==set(saved)
        assert changed_maps==(8 if kind=='B' else 6)
        protected=[]
        for c in over['cases']:
            if c['method']!='finite_leaf9_or_ideal_feature5':continue
            core=over['cores'][c['core_key']];g=A.matrix(core['gram']);base=c['marked_component'];labels=T.T.labels(g,base)
            maps=P.embeddings(labels,target) if labels is not None and len(base)<=13 else []
            images=sorted({tuple(sorted(mp[base.index(i)] for i in c['preferred_corner'])) for mp in maps})
            assert pair not in images
            protected.append(dict(index=c['index'],base_nodes=base,native_labels=labels,all_embeddings=maps,marked_images=images,
                excluded_by_order=len(base)>13,excluded_by_native_labels=labels is None))
        assert len(protected)==74
        records.append(dict(kind=kind,target=target,changed_pair=pair,old_cost='73',new_cost=str(sum(weights.values())),
            weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],complete_native_patterns=801,
            successful_patterns=len(saved),complete_native_maps=sum(r['map_count'] for r in support),affected_maps=changed_maps,
            support=support,cases=proofs,complete_protected_bases=protected,immediate_changed_pairs=[],triangle_pairs=sorted(triangle)))
    routing=Counter()
    for c in over['cases']:
        if c['method']=='finite_leaf9_or_ideal_feature5':routing['protected']+=1
        elif c['method']=='nine_context_adaptive_finite_payment':routing['adaptive']+=1
        elif c['component_order_bound']<=12:routing['whole_bound_at_most12']+=1
        elif len(c['marked_component'])>13:routing['displayed_order_above13']+=1
        else:routing['eligible']+=1
    assert sum(routing.values())==2614 and routing['eligible']==801
    assert C.actual_cross(('A',4),('B',9))==0
    assert data[D13]['new_cost']=='71' and not data[D13]['added_finite_coupons']
    assert next(r for r in data[OLD_C]['order13_interfaces']['records'] if r['kind']=='A')['cost']==65
    paths={Path(__file__),Path(M.__file__),Path(X.__file__),Path(C.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_ORDER13_B_SPARSE_HALF_C_TERMINAL_ZERO_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},records=records,
        new_costs=dict(B='145/2',C='72'),simultaneous_order13_family_costs=dict(A='65',B='145/2',C='72',D='71'),
        simultaneous_order13_maximum='145/2',complete_routing=dict(routing),literal_interface=literal(data[C.TAG],targets),
        new_coupon=dict(small='whole_A4',small_role='either_endpoint',large='whole_B9',large_role='unique_doubled_terminal',
            weight='1/2',eligible_marked_pairs=2,old_actual_component_pair_cost='0',new_actual_component_pair_cost='1',increment='1'),
        triangle_minimum=2,quadrilateral_minimum=1,finite_budget_proved=False,new_dimension_bound=None,
        scope='Complete801 native target searches per family, all74 protected native bases, physical gauge-aware first entries, two full radical seals and216 current literal interfaces. Historical all-label domain and universal short-path proofs are inherited. One actual A4/B9 coupon costs1; no forest or ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],'literal changes',len(out['literal_interface']['changed']),out['complete_routing'],flush=True)
