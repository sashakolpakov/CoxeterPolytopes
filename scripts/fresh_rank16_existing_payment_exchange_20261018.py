#!/usr/bin/env python3
"""Bounded exact rank16 support diagnostic and existing-payment exchange.

Saved complete maps are consumed, not re-enumerated. No forest optimization.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d17_five_unit_reduction_20261018 as V
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

S=V.S;T=S.T;A=S.A;Q=S.Q;D=S.D;P=S.P;ROOT=S.ROOT
TAG='fresh_rank16_existing_payment_exchange_20261018'
AHALF={(0,11),(4,15)};DHALF={(2,13),(2,14),(2,15)};ARESTORE={(1,12),(4,15)}
APROOF=[1657,2067]
DPROPER=[1018,1035,1049,1265,1905,2298]
DRADICAL=[1305,1964,2357,2409,2421,2429]
STEMS=tuple(dict.fromkeys((*V.STEMS,V.TAG,T.TAG,T.T.TAG,V.COUPON,
    'fresh_order14_b_terminal_residual_transfer_20261018',
    'fresh_order14_c_inner_endpoint_projection_transfer_20261018')))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x


def support(kind,half):
    r=next(x for x in half['A16_B16_records']if x['kind']=='A')if kind=='A'else next(x for x in half['records']if x['rank']==16)
    return r['selected_support']


def validate_maps(sup,target,over):
    assert sup['complete_eligible_patterns']==394
    maps=0
    for r in sup['rows']:
        c=over['isolation']['cases'][r['index']];core=over['isolation']['cores'][c['core_key']];g=A.matrix(core['gram'])
        base=c['marked_component'];assert base==r['base_nodes'] and c['preferred_corner']==r['marks']
        graph=T.T.labels(g,base);assert graph is not None or not r['maps']
        seen=set();pairs=set()
        for mp in r['maps']:
            m=mp['images'];assert len(m)==len(base)and len(set(m))==len(m)and all(0<=i<16 for i in m)
            assert tuple(m)not in seen;seen.add(tuple(m))
            assert all(graph[i][j]==target[m[i]][m[j]]for i in range(len(base))for j in range(len(base)))
            pair=tuple(sorted(m[base.index(i)]for i in r['marks']));assert list(pair)==mp['marked_images'];pairs.add(pair);maps+=1
        assert sorted(pairs)==list(map(tuple,r['marked_images']))
    assert maps==sup['complete_maps']
    return maps


def core_proof(c,core,kind,changes):
    g=A.matrix(core['gram']);n=len(g);assert core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i]for i in c['right'][1]];black=set(u+v);white=sorted(set(range(n))-black)
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0]for b in black)
    selected=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],16,kind)
    assert selected and set(selected)<=changes
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    if c['method']=='radical':
        assert n==19 and Q.inertia(g)==[17,1,1]
        z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts'];assert all(z[i]for i in allowed)
        du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8;short.append(dict(nodes=[i,j],length=length))
        assert short
        sealing=dict(method='radical',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=short,
                     conclusion='Every extra white is orthogonal to the whole core unless a full-star short8 path exists.')
    else:
        assert c['index']in DPROPER and kind=='D' and len(c['marked_component'])==16
        gg,J,du,dv,end=X.constraints(c,core);assert gg==g and Q.inertia(g)==[17,1,0]
        corner=sorted(set(white)|set(c['preferred_corner']));assert corner==c['marked_component']
        remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
        assert remaining==c['omitted_black'] and len(remaining)==2
        assert sorted(J[i][i]for i in remaining)==[4,8]
        assert all(sum(i in ee for i in remaining)==1 for ee in c['end_sets'])
        sealing=dict(method='whole_selected_D16_and_Lorentz_projection',whole_selected_nodes=c['marked_component'],
            end_projection=end,remaining_contact_nodes=remaining,inverse_diagonals=[[i,J[i][i]]for i in remaining],
            conclusion='Whole D16 excludes every selected-corner contact. A remaining single-end contact projects above norm2; both ends give short8. Thus every extra white is core-orthogonal in the no-short8 branch.')
    corners=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks));assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
        payment=F();reason=None;images=[];roles=[]
        if marks==c['preferred_corner']:
            payment=F(1,2);reason='retained_selected_half';images=selected
        elif kind=='A'and aa is bb and aa['type']=='A17':
            images=T.mark_pairs(g,aa['nodes'],marks,17,'A');assert images and set(images)<=ARESTORE
            payment=F(1,2);reason='restored_whole_A17_half'
        elif kind=='D'and aa is not bb:
            for sm,la,m,nm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if la['type']!='D13' or len(sm['nodes'])>5:continue
                im=T.mark_pairs(g,la['nodes'],[nm],13);assert im and set(im)<={(10,),(11,),(12,)}
                assert not payment;payment=F(1,2);reason='existing_whole_small_D13_fork_half'
                roles=[S.marked_role(g,sm,m),S.marked_role(g,la,nm)]
        corners.append(dict(marks=marks,components=blocks,payment=payment,reason=reason,canonical_images=images,roles=roles))
    assert sum(r['payment']for r in corners)==1 and sum(bool(r['payment'])for r in corners)==2
    return dict(index=c['index'],kind=kind,core_key=c['core_key'],gram=core['gram'],selected_images=selected,
        white=white,edge_positive_definite=True,sealing=sealing,corners=corners,total='1')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    half=docs[T.TAG];sp=docs[S.TAG];over=docs[A.OVER];d17=docs[V.TAG]
    assert d17['new_cost']=='120'
    rules={r['kind']:r for r in sp['rules']if r['rank']==16}
    weights={k:{tuple(x['pair']):F(x['weight'])for x in r['weighted_pairs']}for k,r in rules.items()}
    oldA17=next(r for r in docs[T.T.TAG]['records']if r['kind']=='A')
    a17={p:F(int(list(p)in oldA17['retained_pairs']))for p in combinations(range(17),2)}
    assert all(a17[p]==0 for p in ARESTORE)and sum(a17.values())==114
    assert {tuple(sorted(16-i for i in p))for p in ARESTORE}==ARESTORE
    assert all(weights['A'][p]==1 for p in AHALF)and all(weights['D'][p]==1 for p in DHALF)
    sups={k:support(k,half)for k in ('A','D')};proofs=[];inventories=[]
    for kind in ('A','D'):
        target=rules[kind]['target'];sup=sups[kind];count=validate_maps(sup,target,over);assert count=={'A':238,'D':578}[kind]
        adj,dist,_=P.distance_data(target);autos=rules[kind]['automorphisms'];changes=AHALF if kind=='A'else DHALF
        assert all({tuple(sorted(m[i]for i in p))for p in changes}==changes for m in autos)
        assert all(dist[i][j]>=10 for i,j in changes)
        rows=[r for r in sup['rows']if set(map(tuple,r['marked_images']))&changes]
        assert sorted(r['index']for r in rows)==(APROOF if kind=='A'else sorted(DPROPER+DRADICAL))
        for r in rows:
            c=over['isolation']['cases'][r['index']];proofs.append(core_proof(c,over['isolation']['cores'][c['core_key']],kind,changes))
        seen=set();candidate=[]
        for p,w in sorted(weights[kind].items()):
            if p in seen or not w or dist[p[0]][p[1]]<10:continue
            orbit={tuple(sorted(m[i]for i in p))for m in autos};seen|=orbit
            hits=[r for r in sup['rows']if set(map(tuple,r['marked_images']))&orbit]
            candidate.append(dict(orbit=sorted(orbit),current_weight=w,half_debit_component_saving=F(len(orbit),2),
                selected_contexts=[dict(index=r['index'],method=r['method'],component_size=len(r['base_nodes']))for r in hits]))
        inventories.append(dict(kind=kind,patterns=394,maps=count,remaining_long_orbits=candidate))
    # Simultaneous debit from the latest complete literal record, including D17=120.
    previous={(r['source'],V.normalize(r['case'])):r for r in d17['literal_interfaces']['rows']};assert len(previous)==216
    literal=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        old=previous[stem,V.normalize(index)];corners=[]
        for j,(marks,*_)in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));delta=F();images=[];typ=None
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                shared=next((b for b in D.finite_blocks(g,ids)if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type']in ('A16','D16','A17'):
                    typ=shared['type'];images=T.mark_pairs(g,shared['nodes'],marks,int(typ[1:]),typ[0])
                    changed=AHALF if typ=='A16'else DHALF if typ=='D16'else ARESTORE
                    if set(images)&changed:
                        assert set(images)<=changed;delta=F(1,2)if typ=='A17'else F(-1,2);hits[typ]+=1
            val=F(old['corners'][j]['retained_payment'])+delta;assert val>=0
            corners.append(dict(marks=marks,type=typ,images=images,previous=old['corners'][j]['retained_payment'],delta=delta,retained=val))
        total=sum(r['retained']for r in corners);assert total>=1
        literal.append(dict(source=stem,case=index,corners=corners,total=total))
    assert len(literal)==216 and not hits
    paths=[Path(__file__),Path(V.__file__),Path(X.__file__),Path(S.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(P.__file__),Path(D.__file__)]
    paths +=[ROOT/f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return enc(dict(status='PASS_RANK16_EXISTING_PAYMENT_EXCHANGE_SUPPORT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        candidate_inventory=inventories,A16_half_pairs=sorted(AHALF),D16_half_pairs=sorted(DHALF),A17_restored_half_pairs=sorted(ARESTORE),
        old_costs=dict(A16='114',D16='114',A17='114'),new_costs=dict(A16='113',D16='225/2',A17='115'),
        selected_proofs=proofs,selected_contexts=14,whole_Lorentz_projection_cases=6,radical_cases=8,
        literal_configurations=216,literal_changed_corners=dict(hits),literal_rows=literal,new_cross_component_coupons=[],
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Consumes saved complete394-pattern/238-map A16 and394-pattern/578-map D16 domains, validates every saved map, and checks14 exact local exchange/sealing cases and all216 current literal interfaces. No original overlap/embedding census or finite forest optimization is rerun. The A16 saving restores two existing A17 long-pair halves; it is not a pure finite decrease. Written universal/current-priority composition supplies the full local face theorem; no ambient bound or new coupon is claimed.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:
        assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],'contexts',out['selected_contexts'],'literal',out['literal_configurations'],flush=True)
