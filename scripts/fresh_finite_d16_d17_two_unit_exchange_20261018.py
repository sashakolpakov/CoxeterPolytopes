#!/usr/bin/env python3
"""Delete two finite D16 units and restore two finite D17 units.

Complete selected maps, radical sealing and simultaneous literal interfaces.
The accompanying written proof treats arbitrary additional supporting roots.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d16_d17_long_pair_half_transfer_20261017 as T

ROOT=T.ROOT; A=T.A; P=T.P; Q=T.Q; D=T.D
TAG='fresh_finite_d16_d17_two_unit_exchange_20261018'
CURRENT='fresh_deep_finite8_dimension64_20261018'
ZERO='fresh_depth11_finite8_zero_transfer_20261018'
DELETE={(0,11),(0,12)}
RESTORE={(1,12),(1,13)}
EXPECTED=[1681,1685,2500,2532]
STEMS=tuple(dict.fromkeys((T.TAG,*T.STEMS,CURRENT,ZERO)))


def normalize_case(x):
    return tuple(x) if isinstance(x,list) else x


def radical_exchange(c,core):
    g=A.matrix(core['gram']); N=len(g)
    assert N==19 and core['gauge']==[1]*N and c['method']=='radical'
    assert c['core_key']=='(25, 25, ((9, 9),))'
    assert Q.inertia(g)==[17,1,1]
    z=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*z[j] for j in range(N))==0 for i in range(N))
    u=c['left'][1]; v=[c['right_copy_map'][i] for i in c['right'][1]]
    black=set(u+v); white=sorted(set(range(N))-black)
    assert len(black)==4 and Q.inertia(A.sub(g,white))==[15,0,0]
    for b in black:
        assert Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0]
    actual=Q.inspect_corner(g,[1]*N,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    assert len(c['marked_component'])==c['component_order_bound']==16
    selected_images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],16)
    assert selected_images and set(selected_images)<=DELETE
    allowed=c['possible_contacts']; assert all(z[i] for i in allowed)
    du=Q.graph_dist(g,u); dv=Q.graph_dist(g,v); shortcuts=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2
            assert length<=8
            shortcuts.append(dict(nodes=[i,j],length=length))
    assert shortcuts
    corners=[]; alternate_total=F()
    for old in c['all_corners']:
        marks=old['preferred_corner']; ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids)
        shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
        images=[]; value=F()
        if shared and shared['type']=='D17':
            images=T.mark_pairs(g,shared['nodes'],marks,17)
            if set(images)&RESTORE:
                assert set(images)<=RESTORE and marks!=c['preferred_corner']
                value=F(1)
        alternate_total+=value
        corners.append(dict(marks=marks,components=blocks,whole_D17_images=images,
                            restored_unit_used=str(value)))
    assert alternate_total==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],
                radical=list(map(str,z)),white=white,selected_images=selected_images,
                allowed_contacts=allowed,mixed_sign_shortcuts=shortcuts,
                whole_core_sealed_unless_short8=True,corners=corners,
                deleted_selected_payment='0',alternate_restored_payment='1',total='1')


def literal_audit(data):
    old=data[T.TAG]['literal_interfaces']['rows']
    lookup={(r['source'],normalize_case(r['case'])):r for r in old}
    assert len(lookup)==216
    rows=[]
    for stem,index,g,white,cs in T.literal_inputs(data):
        baseline=lookup[(stem,normalize_case(index))]
        assert len(baseline['corners'])==len(cs)
        corners=[]
        for j,(marks,unit,_) in enumerate(cs):
            previous=F(baseline['corners'][j]['retained_weight'])
            assert marks==baseline['corners'][j]['marks']
            ids=sorted(set(white)|set(marks)); parts=A.components(g,ids)
            shared=next((part for part in parts if set(marks)<=set(part)),None)
            debit=credit=F(); images=[]
            if unit and shared and len(shared) in (16,17):
                n=len(shared); images=T.mark_pairs(g,shared,marks,n)
                if n==16 and set(images)&DELETE:
                    assert set(images)<=DELETE and previous>=1
                    debit=F(1)
                if n==17 and set(images)&RESTORE:
                    assert set(images)<=RESTORE and previous==0
                    assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
                    credit=F(1)
            retained=previous-debit+credit; assert retained>=0
            corners.append(dict(marks=marks,finite_unit=unit,whole_D_images=images,
                                old_current_weight=str(previous),debit=str(debit),
                                credit=str(credit),retained_weight=str(retained)))
        total=sum(F(c['retained_weight']) for c in corners); assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,retained_total=str(total)))
    assert len(rows)==216
    affected=[r for r in rows if any(F(c['debit']) or F(c['credit']) for c in r['corners'])]
    assert [(r['source'],normalize_case(r['case'])) for r in affected]==[
        (D.TAG,(i,0)) for i in (1975,1979,2505,2536)]
    assert all(r['retained_total']=='1' for r in affected)
    return dict(complete_configurations=216,affected_configurations=4,rows=rows,
                minimum_retained_total=str(min(F(r['retained_total']) for r in rows)))


def run():
    data={stem:json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem,doc in data.items():
        assert doc['status'].startswith('PASS'),stem
        for path,digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,(stem,path)
    over=data[A.OVER]
    support=T.selected(over,16)
    affected=[r for r in support['rows'] if set(map(tuple,r['marked_images']))&DELETE]
    assert [r['index'] for r in affected]==EXPECTED
    assert support['complete_eligible_patterns']==394 and support['complete_maps']==578
    assert all(set(map(tuple,r['marked_images']))<=DELETE for r in affected)
    proofs=[radical_exchange(over['isolation']['cases'][i],
              over['isolation']['cores'][over['isolation']['cases'][i]['core_key']]) for i in EXPECTED]
    rules=[]
    for n,change in ((16,DELETE),(17,RESTORE)):
        graph=T.target(n); adj,dist,_=P.distance_data(graph); autos=P.embeddings(graph,graph)
        assert all({tuple(sorted(m[i] for i in pair)) for pair in change}==change for m in autos)
        assert all(dist[a][b] in (11,12) for a,b in change)
        old=next(r for r in data[T.TAG]['records'] if r['rank']==n)
        weights={pair:F() for pair in combinations(range(n),2)}
        weights.update({tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']})
        assert sum(weights.values())==F(old['new_cost'])==(118 if n==16 else 125)
        assert all(weights[pair]==(1 if n==16 else 0) for pair in change)
        for pair in change: weights[pair]=F(0 if n==16 else 1)
        assert sum(weights.values())==(116 if n==16 else 127)
        basic={pair for pair in weights if dist[pair[0]][pair[1]]<=8 or
               (dist[pair[0]][pair[1]]==9 and any(len(adj[i])==1 for i in pair))}
        assert not(basic&change) and all(weights[pair]==1 for pair in basic)
        rules.append(dict(rank=n,changed_pairs=sorted(change),automorphisms=autos,
                          old_cost=old['new_cost'],new_cost=str(sum(weights.values())),
                          weighted_pairs=[dict(pair=pair,weight=str(w)) for pair,w in sorted(weights.items())],
                          basic_short8_and_leaf9_pairs=sorted(basic)))
    literals=literal_audit(data)
    adaptive=data[T.TAG]['adaptive_bound16_interfaces']
    assert len(adaptive)==2 and all(r['marked_distance']==9 for r in adaptive)
    assert all(tuple(m['marked_images']) not in DELETE for r in adaptive for m in r['maps'])
    assert all(r['type'] not in ('D16','D17') for r in data[T.TAG]['later_typed_units'])
    for stem in (T.NEW15,T.NEW17):
        assert data[stem]['used_finite_same_component_units'] is False
    paths=[Path(__file__),Path(T.__file__),Path(A.__file__),Path(P.__file__),Path(Q.__file__),Path(D.__file__)]
    paths += [ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE_D16_D17_TWO_UNIT_EXCHANGE',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                rules=rules,complete_high_overlap_contexts=2614,selected_D16_support=support,
                affected_selected_contexts=EXPECTED,radical_exchange_proofs=proofs,
                literal_interfaces=literals,adaptive_bound16_interfaces=adaptive,
                new_finite_component_costs=dict(A16=116,B16=104,D16=116,D17=127),
                all_coupons_unchanged=True,all_ideal_weights_unchanged=True,
                depth11_zero_compatible=True,finite_maximum_proved_here=False,
                new_dimension_bound=None,
                scope='Complete selected D16 maps and four all-label radical transfers to restored whole D17 units; all216 simultaneous literal interfaces and retained universal interfaces. Arbitrary supporting roots are handled by the written radical argument and completed short-path theorems. No ambient bound is claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify',action='store_true'); args=parser.parse_args()
    out=json.loads(json.dumps(run())); path=ROOT/f'results/{TAG}.json'
    if args.verify: assert out==json.loads(path.read_text())
    else: path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_finite_component_costs'],flush=True)
