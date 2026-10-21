#!/usr/bin/env python3
"""Finite rank16/17 refinements with one precisely typed D14 coupon.

Local geometry and complete affected interfaces; no finite-forest maximum.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d17_long_pair_support_20261017 as T
import fresh_d16_partial_weight_transfer_20261017 as D

ROOT=T.ROOT;A=T.A;P=T.P;Q=D.Q
TAG='fresh_finite_d16_d17_long_pair_half_transfer_20261017'
NEW15='fresh_order15_affine_d_half_coupon_transfer_20261017'
NEW17='fresh_order17_affine_d_three_quarter_transfer_20261017'
BASE='fresh_ambient65_ideal_margin_closure_20261017'
FINITE='fresh_e6_d11_quarter_coupon_forest_budget_20261017'
EXTRA=(T.NEW_CYCLE,T.NEW_D,T.NEW_B,T.NEW_C,NEW15,NEW17,
    'fresh_order12_bc_unit_coupon_transfer_20261017',
    'fresh_order13_affine_b_fork_projection_transfer_20261017',
    'fresh_order14_affine_b_fork_projection_transfer_20261017',
    'fresh_order15_affine_b_residual_terminal_transfer_20261017',
    'fresh_order15_affine_c_terminal_projection_transfer_20261017',
    'fresh_large_finite_d_leaf_half_transfer_20261017')
STEMS=tuple(dict.fromkeys((T.TAG,A.TAG,*A.STEMS,*EXTRA,BASE,FINITE)))
EXPECTED={16:[1146,1161,1174,1185,1268,1909,2302,2491],
          17:[1017,1034,1048,1059,1264,1904,2297,2488]}


def target(n,kind='D'):
    edges=[(i,i+1,3) for i in range(n-1)]
    if kind=='B':edges[-1]=(n-2,n-1,4)
    if kind=='D':edges=[(i,i+1,3) for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)]
    return P.labelled(n,edges)


def changed(n):
    return {(n-15,j) for j in range(n-4,n)}


def mark_pairs(g,ids,marks,n,kind='D'):
    return sorted({tuple(m['marked_images']) for m in T.mapped(g,ids,marks,target(n,kind))})


def selected(over,n,kind='D'):
    rows=[];maps_count=0;pair_union=set();affected=[]
    for c in over['isolation']['cases']:
        if (c.get('component_order_bound') or 0)<n or len(c['marked_component'])>n:continue
        core=over['isolation']['cores'][c['core_key']];g=A.matrix(core['gram'])
        maps=T.mapped(g,c['marked_component'],c['preferred_corner'],target(n,kind))
        pairs=sorted({tuple(m['marked_images']) for m in maps})
        row=dict(index=c['index'],core_key=c['core_key'],method=c['method'],
                 component_order_bound=c['component_order_bound'],
                 base_nodes=c['marked_component'],marks=c['preferred_corner'],maps=maps,marked_images=pairs)
        rows.append(row);maps_count+=len(maps);pair_union.update(pairs)
        if kind=='D' and set(pairs)&changed(n):
            assert set(pairs)<=changed(n)
            affected.append(row)
    assert len(rows)==(394 if n==16 else 140)
    if kind=='D':assert [r['index'] for r in affected]==EXPECTED[n]
    return dict(rank=n,kind=kind,complete_eligible_patterns=len(rows),
                complete_maps=maps_count,marked_pair_union=sorted(pair_union),rows=rows,affected=affected)


def sealed(c,core,n):
    g=A.matrix(core['gram']);N=len(g)
    assert N==19 and core['gauge']==[1]*N and c['method']=='radical'
    assert c['core_key']=='(25, 25, ((9, 9),))' and Q.inertia(g)==[17,1,1]
    z=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*z[j] for j in range(N))==0 for i in range(N))
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    black=set(u+v);white=sorted(set(range(N))-black)
    assert len(black)==4 and Q.inertia(A.sub(g,white))==[15,0,0]
    for b in black:assert Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0]
    actual=Q.inspect_corner(g,[1]*N,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    assert len(c['marked_component'])==n
    assert set(mark_pairs(g,c['marked_component'],c['preferred_corner'],n))<=changed(n)
    assert Q.inertia(A.sub(g,c['marked_component']))==[n,0,0]
    allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);shortcuts=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2
            assert length<=8
            shortcuts.append(dict(nodes=[i,j],length=length))
    assert shortcuts
    # A radical of an actual indefinite Lorentz span represents the zero
    # vector. Thus any new root's contact vector annihilates z. All allowed
    # coefficients are nonzero, and mixed signs give the retained shortcut.
    corners=[];coupon_total=F();selected_total=F()
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        bs=D.finite_blocks(g,ids)
        if marks==c['preferred_corner']:
            selected_total+=F(1,2)
            corners.append(dict(marks=marks,components=bs,used_selected_half='1/2',used_coupon='0'))
            continue
        selected_blocks=[next(b for b in bs if m in b['nodes']) for m in marks]
        coupon=[]
        for small,large,sm,lm in ((selected_blocks[0],selected_blocks[1],marks[0],marks[1]),
                                  (selected_blocks[1],selected_blocks[0],marks[1],marks[0])):
            if small['type'] not in ('A1','A2') or large['type']!='D14':continue
            assert small['type']==('A2' if n==16 else 'A1')
            images=mark_pairs(g,large['nodes'],[lm],14)
            assert images and all(p[0] in (10,11,12,13) for p in images)
            new=all(p==(10,) for p in images)
            assert new or large['roles'][str(lm)] in ('fork_center','short_fork_leaf')
            coupon.append(dict(small_type=small['type'],small_nodes=small['nodes'],small_mark=sm,
                large_type='D14',large_nodes=large['nodes'],large_mark=lm,
                large_canonical_mark_images=images,new_coupon=new,weight='1/2',whole_components=True))
        assert len(coupon)<=1
        value=F(1,2) if coupon else F();coupon_total+=value
        corners.append(dict(marks=marks,components=bs,used_selected_half='0',
            used_coupon=str(value),coupon=coupon,other_same_component_units_used=False))
    assert selected_total==coupon_total==F(1,2)
    return dict(index=c['index'],core_key=c['core_key'],core_gram=core['gram'],
        radical=list(map(str,z)),white=white,allowed_contacts=allowed,
        all_mixed_sign_shortcuts=shortcuts,whole_core_sealed_unless_short8=True,
        corners=corners,selected_half='1/2',finite_coupon_half='1/2',total='1')


def literal_inputs(data):
    """Every full literal interface used by the prior rank17 proof, plus D16."""
    g17=A.matrix(data[A.OVER]['isolation']['cores']['(25, 25, ((9, 9),))']['gram'])
    g16=A.matrix(data[A.OLD16]['exceptional_D15']['gram'])
    for stem,g,cases in ((A.OLD16,g16,data[A.OLD16]['exceptional_D15']['cases']),
                         (A.OLD17,g17,data[A.OLD17]['order17']['replacement_cases'])):
        for r in cases:
            white=sorted(set(range(len(g)))-set(r['black']))
            yield stem,r['case'],g,white,[(c['marks'],int(c['finite'] and F(c['base'])==1),F(c['weight'])) for c in r['corners']]
    for stem,key,unitkey in ((A.D13,'extreme_equality','unit'),(A.D14,'equality','finite_small17_unit')):
        for i,r in enumerate(data[stem][key]['rows']):
            g=A.d13_gram(r) if stem==A.D13 else A.matrix(r['gram'])
            yield stem,i,g,r['white'],[(c['marks'],c[unitkey],F(c[unitkey])+F(c['half_coupon'])) for c in r['corners']]
    for r in data[T.PARTIAL]['cases']:
        yield T.PARTIAL,r['index'],A.matrix(r['gram']),r['white'],[
            (c['marks'],c['finite_unit'],F(c['finite_unit'])+F(c['finite_half'])+F(c['ideal_weight'])) for c in r['corners']]
    for r in data[D.TAG]['cases']:
        for j,p in enumerate(r['completions']):
            if 'corners' not in p:continue
            assert len(p['gram'])==18 and len(p['white'])==14
            assert all(len(b['nodes'])<=16 for c in p['corners'] for b in c.get('components',[]))
            yield D.TAG,(r['index'],j),A.matrix(p['gram']),p['white'],[
                (c['marks'],c['finite_unit'],F(c['finite_unit'])+F(c['finite_half'])+F(c['ideal_weight'])) for c in p['corners']]
    for j,r in enumerate(data[T.NEW_D]['two_A1_A3_spectators']['rows']):
        assert all(len(b['nodes'])<=16 for c in r['corners'] for b in c['components'])
        yield T.NEW_D,j,A.matrix(r['gram']),r['white'],[
            (c['marks'],c['unit'],F(c['unit'])+F(c['old_half'])) for c in r['corners']]


def literal_audit(data,ab_removed):
    old=data[T.TAG]['simultaneous_literal_totals'];oldkey={(r['source'],r['case']):r for r in old}
    rows=[];ab_hits=[]
    for stem,index,g,white,cs in literal_inputs(data):
        corners=[]
        for j,(marks,unit,weight) in enumerate(cs):
            ids=sorted(set(white)|set(marks));parts=A.components(g,ids)
            shared=next((p for p in parts if set(marks)<=set(p)),None)
            # First account for the completed simultaneous A/B/D17 deletion.
            key=(stem,tuple(index) if isinstance(index,list) else index)
            previous=F(oldkey[key]['corners'][j]['retained_weight']) if stem not in (D.TAG,T.NEW_D) else weight
            debit=F();images=[];ab_images=[]
            if unit and shared and len(shared) in (16,17):
                n=len(shared);images=mark_pairs(g,shared,marks,n)
                if set(images)&changed(n):
                    assert set(images)<=changed(n);debit=F(1,2)
                if n==16:
                    for kind in ('A','B'):
                        ps=mark_pairs(g,shared,marks,n,kind)
                        if ps:
                            ab_hits.append(dict(source=stem,case=index,marks=marks,kind=kind,images=ps))
                            ab_images.append(dict(kind=kind,images=ps))
                            if set(ps)&ab_removed[kind]:
                                assert set(ps)<=ab_removed[kind] and not debit
                                debit=F(1)
            value=previous-debit;assert value>=0
            corners.append(dict(marks=marks,old_current_weight=str(previous),finite_unit=unit,
                whole_D_images=images,whole_AB16_images=ab_images,debit=str(debit),retained_weight=str(value)))
        total=sum(F(c['retained_weight']) for c in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,retained_total=str(total)))
    assert len(rows)==216
    affected=[r for r in rows if any(F(c['debit']) for c in r['corners'])]
    assert len(affected)>=25
    return dict(complete_literal_configurations=216,affected_configurations=len(affected),rows=rows,
                affected_source_counts=dict(Counter(r['source'] for r in affected)),
                minimum_retained_total=str(min(F(r['retained_total']) for r in rows)),A16_B16_unit_images=ab_hits)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,doc in data.items():
        assert doc['status'].startswith('PASS'),s
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    over=data[A.OVER];records=[]
    for n in (16,17):
        support=selected(over,n);graph=target(n);autos=P.embeddings(graph,graph)
        adj,dist,_=P.distance_data(graph);change=changed(n)
        assert all(dist[a][b]>=11 for a,b in change)
        basic={p for p in combinations(range(n),2) if dist[p[0]][p[1]]<=8 or
               (dist[p[0]][p[1]]==9 and any(len(adj[i])==1 for i in p))}
        assert not(change&basic)
        assert all({tuple(sorted(m[i] for i in p)) for p in change}==change for m in autos)
        if n==16:paid=set(combinations(range(n),2))
        else:paid=set(map(tuple,next(r for r in data[T.TAG]['records'] if r['kind']=='D')['retained_pairs']))
        assert change<=paid and len(paid)==(120 if n==16 else 127)
        proofs=[sealed(over['isolation']['cases'][i],over['isolation']['cores'][over['isolation']['cases'][i]['core_key']],n)
                for i in EXPECTED[n]]
        assert sum(len(r['maps']) for r in support['affected'])==16
        records.append(dict(rank=n,target=graph,automorphisms=autos,changed_pairs=sorted(change),
            old_cost=len(paid),new_cost=str(F(len(paid))-2),basic_pairs_retained=sorted(basic),
            weighted_pairs=[dict(pair=p,weight='1/2' if p in change else '1') for p in sorted(paid)],
            selected_support=support,proofs=proofs))
    # Both adaptive bound16 corners still have marked distance9 in their
    # induced tree, so none of the present distance>=11 changes affects them.
    adaptive=[]
    for r in data[A.ADAPT]['cases']:
        p=r['payment']
        if 'whole_component_bound' in p:
            g=A.matrix(r['gram']);marks=p['corner'];z=next(z for z in r['all_corners'] if z['preferred_corner']==marks)
            ns=z['marked_component'];maps=T.mapped(g,ns,marks,target(16))
            assert all(tuple(m['marked_images']) not in changed(16) for m in maps)
            assert Q.graph_dist(A.sub(g,ns),[ns.index(marks[0])])[ns.index(marks[1])]==9
            adaptive.append(dict(index=r['index'],marks=marks,base_nodes=ns,maps=maps,marked_distance=9))
    assert len(adaptive)==2
    ab_records=[];ab_removed={}
    for kind in ('A','B'):
        support=selected(over,16,kind);graph=target(16,kind);adj,dist,_=P.distance_data(graph)
        pairs=set(combinations(range(16),2))
        basic={p for p in pairs if dist[p[0]][p[1]]<=8 or
               (dist[p[0]][p[1]]==9 and any(len(adj[i])==1 for i in p))}
        remove=pairs-basic-set(map(tuple,support['marked_pair_union']));autos=P.embeddings(graph,graph)
        assert len(remove)==(4 if kind=='A' else 16)
        assert all({tuple(sorted(m[i] for i in p)) for p in remove}==remove for m in autos)
        assert all(dist[a][b]>=10 for a,b in remove)
        ab_removed[kind]=remove
        ab_records.append(dict(kind=kind,rank=16,target=graph,automorphisms=autos,
            removed_pairs=sorted(remove),retained_pairs=sorted(pairs-remove),basic_pairs_retained=sorted(basic),
            old_cost=120,new_cost=120-len(remove),selected_support=support))
    literals=literal_audit(data,ab_removed)
    # Later D11 equality outputs are explicitly D13/D15/A16, never D16/17.
    typed=[]
    for stem in (T.NEW_D,'fresh_affine_d12_leaf_pair_transfer_20261017',T.NEW_B,'fresh_order12_bc_unit_coupon_transfer_20261017'):
        for path,c in A.all_dicts(data[stem]):
            if not('components' in c and 'marks' in c and c.get('finite_unit',c.get('unit',0))):continue
            for b in c['components']:
                if set(c['marks'])<=set(b['nodes']):
                    assert b['type'] not in ('D16','D17')
                    typed.append(dict(source=stem,path=path,type=b['type']))
    for stem in (NEW15,NEW17):
        assert data[stem]['used_finite_same_component_units'] is False
    paths=[Path(__file__),Path(T.__file__),Path(A.__file__),Path(D.__file__),Path(P.__file__),Path(Q.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE_D16_D17_LONG_PAIR_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        records=records,A16_B16_records=ab_records,complete_high_overlap_contexts=2614,
        literal_interfaces=literals,adaptive_bound16_interfaces=adaptive,later_typed_units=typed,
        added_coupon=dict(weight='1/2',small_components=['whole A1','whole A2'],small_mark='any node',
            large_component='whole D14',large_mark='the unique long-arm neighbor of its fork center',
            canonical_D14_node=10,added_A1_D14_cost='1/2',added_A2_D14_cost='1'),
        finite_maximum_proved_here=False,new_dimension_bound=None,
        new_finite_component_costs=dict(A16=116,B16=104,D16=118,D17=125),
        scope='Complete selected finite A16/B16/D16/D17 maps,16 sealed actual radical cores and32 affected D maps,216 simultaneous literal configurations, adaptive and later interfaces. Whole-component geometric reduction and universal short-path branches inherited. New D14 coupon requires a separate exact finite budget.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],[(r['rank'],r['new_cost']) for r in out['records']],flush=True)
