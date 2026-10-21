#!/usr/bin/env python3
"""Whole finite A/B/D15 cost104 using only existing distinct halves."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_rank15_long_pair_support_scout_20261018 as R
import fresh_b14_b15_singleton_completion_transfer_20261017 as X
import fresh_rank16_existing_payment_exchange_20261018 as E

ROOT=R.ROOT;N=R.N;T=R.T;A=R.A;P=R.P;Q=R.Q;D=N.D;S=N.S
TAG='fresh_finite_rank15_existing_half_transfer_20261018'
CHANGES={'A':{(0,13),(1,14)},'B':{(0,14)},'D':{(0,13),(0,14)}}
EXPECTED={'A':[1156,1616],'B':[],'D':[1645,1665,1932,2325]}
STEMS=tuple(dict.fromkeys((*R.STEMS,R.TAG,E.TAG)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def validate_support(rec,over):
    assert rec['complete_eligible_patterns']==536
    target=rec['target'];count=0;hits=[]
    for row in rec['selected_patterns']:
        c=over['isolation']['cases'][row['index']]
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        base=c['marked_component'];graph=T.T.labels(g,base);seen=set();pairs=set()
        for m in row['maps']:
            mp=m['images'];assert len(mp)==len(base)==len(set(mp))
            assert all(0<=i<15 for i in mp) and tuple(mp) not in seen;seen.add(tuple(mp))
            assert all(graph[i][j]==target[mp[i]][mp[j]] for i in range(len(base)) for j in range(len(base)))
            pair=tuple(sorted(mp[base.index(i)] for i in c['preferred_corner']))
            assert list(pair)==m['marked_images'];pairs.add(pair);count+=1
        assert sorted(pairs)==list(map(tuple,row['marked_images']))
        if pairs&CHANGES[rec['kind']]:
            assert pairs<=CHANGES[rec['kind']];hits.append(row['index'])
    assert count==dict(A=482,B=240,D=870)[rec['kind']]
    assert hits==EXPECTED[rec['kind']]
    return dict(kind=rec['kind'],patterns=536,maps=count,affected_cases=hits)


def core_proof(c,core,kind,d16weights):
    g=A.matrix(core['gram']);n=len(g)
    assert core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    white=sorted(set(range(n))-set(u+v))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in u+v)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    assert len(c['marked_component'])==15
    images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],15,kind)
    assert images and set(images)<=CHANGES[kind]
    if kind=='A':
        assert c['method']=='radical' and n==19 and Q.inertia(g)==[17,1,1]
        z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
        du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
                short.append(dict(nodes=[i,j],length=length))
        assert short
        sealing=dict(method='zero_radical',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=short)
    else:
        gg,J,du,dv,end=X.constraints(c,core);assert gg==g and n==18
        remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
        assert sorted(J[i][i] for i in remaining)==[4,4,14]
        assert all(J[i][i]>2 for i in remaining)
        sealing=dict(method='whole_D15_and_positive_Lorentz_end_blocks',
            whole_selected_nodes=c['marked_component'],end_projection=end,
            remaining_contact_nodes=remaining,inverse_diagonals=[[i,J[i][i]] for i in remaining])
    corners=[];alternatives=0
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        bs=D.finite_blocks(g,ids);aa,bb=[next(b for b in bs if m in b['nodes']) for m in marks]
        payment=F();reason=None;roles=[];pairs=[]
        if marks==c['preferred_corner']:
            payment=F(1,2);reason='selected_rank15_half';pairs=images
        elif kind=='A' and aa is not bb:
            for sm,la,m,nm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if sm['type']=='A2' and la['type']=='A13':
                    small=T.mark_pairs(g,sm['nodes'],[m],2,'A')
                    large=T.mark_pairs(g,la['nodes'],[nm],13,'A')
                    assert small==[(0,),(1,)] and large==[(0,),(12,)]
                    payment=F(1,2);reason='existing_whole_A2_A13_endpoint_half'
                    roles=[S.marked_role(g,sm,m),S.marked_role(g,la,nm)];alternatives+=1
        elif kind=='D' and aa is bb and aa['type']=='D16':
            pairs=T.mark_pairs(g,aa['nodes'],marks,16)
            assert pairs==[(1,14),(1,15)]
            assert all(d16weights[p]==F(1,2) and p not in E.DHALF for p in pairs)
            payment=F(1,2);reason='retained_whole_D16_half';alternatives+=1
        corners.append(dict(marks=marks,components=bs,payment=payment,reason=reason,
                            canonical_images=pairs,roles=roles))
    assert alternatives==1 and sum(r['payment'] for r in corners)==1
    return dict(index=c['index'],kind=kind,core_key=c['core_key'],gram=core['gram'],
                white=white,selected_images=images,sealing=sealing,corners=corners,total='1')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    scout=docs[R.TAG];over=docs[A.OVER]
    assert scout['new_weight_rule_proved'] is False
    assert docs[E.TAG]['new_costs']==dict(A16='113',D16='225/2',A17='115')
    assert docs[E.TAG]['literal_changed_corners']=={}
    d16=next(r for r in docs[S.TAG]['rules'] if r['kind']=='D' and r['rank']==16)
    d16weights={tuple(r['pair']):F(r['weight']) for r in d16['weighted_pairs']}
    assert all(d16weights[p]==F(1,2) for p in ((1,14),(1,15)))
    coupon=next(c for c in docs[S.TAG]['coupons'] if c['name']=='A2_A13end')
    assert coupon['weight']=='1/2' and coupon['small_nodes']==[0,1] and coupon['large_nodes']==[0,12]
    rules=[];supports=[];proofs=[]
    for rec in scout['records']:
        kind=rec['kind'];target=rec['target'];adj,dist,_=P.distance_data(target)
        changes=CHANGES[kind];new=F() if kind=='B' else F(1,2)
        assert all({tuple(sorted(mp[i] for i in p)) for p in changes}==changes for mp in rec['automorphisms'])
        assert all(dist[i][j]>=10 for i,j in changes)
        weights={p:(new if p in changes else F(1)) for p in combinations(range(15),2)}
        assert sum(weights.values())==104
        assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
        supports.append(validate_support(rec,over))
        for i in EXPECTED[kind]:
            c=over['isolation']['cases'][i];proofs.append(core_proof(c,over['isolation']['cores'][c['core_key']],kind,d16weights))
        rules.append(dict(kind=kind,rank=15,target=target,automorphisms=rec['automorphisms'],
            changed_pairs=sorted(changes),new_weight=new,old_cost='105',new_cost='104',
            weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
            retained_basic_pairs=rec['retained_basic_pairs']))
    # Reconstruct all current literal corners, debit all rank15 changes at once,
    # and also apply the independently frozen rank16/A17 exchange.
    old={(r['source'],N.normalize(r['case'])):r for r in docs[N.TAG]['literal_interfaces']['rows']}
    literal=[];hits=Counter();exchange_hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=old[stem,N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks))
            debit=F();delta=F();typ=None;pairs=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                shared=next((b for b in D.finite_blocks(g,ids) if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type'] in ('A15','B15','D15','A16','D16','A17'):
                    typ=shared['type'];pairs=T.mark_pairs(g,shared['nodes'],marks,int(typ[1:]),typ[0])
                    if typ.endswith('15') and set(pairs)&CHANGES[typ[0]]:
                        assert set(pairs)<=CHANGES[typ[0]] and prev['corners'][j]['retained_payment']=='1'
                        debit=F(1) if typ=='B15' else F(1,2);hits[typ]+=1
                    elif typ in ('A16','D16','A17'):
                        change=E.AHALF if typ=='A16' else E.DHALF if typ=='D16' else E.ARESTORE
                        if set(pairs)&change:
                            assert set(pairs)<=change;delta=F(1,2) if typ=='A17' else F(-1,2);exchange_hits[typ]+=1
            value=F(prev['corners'][j]['retained_payment'])-debit+delta;assert value>=0
            corners.append(dict(marks=marks,type=typ,images=pairs,previous=prev['corners'][j]['retained_payment'],
                rank15_debit=debit,rank16_exchange_delta=delta,retained_payment=value))
        total=sum(r['retained_payment'] for r in corners);assert total>=1
        literal.append(dict(source=stem,case=index,corners=corners,total=total))
    assert len(literal)==216 and hits=={'D15':16} and not exchange_hits
    # Complete no-short8 routing: only the536 rank15-capable selected patterns
    # can use a long whole-rank15 pair; the83 protected/adaptive cases retain
    # their marked distance9/basic interfaces.
    cases=over['isolation']['cases'];assert len(cases)==2614
    low_bound=[c for c in cases if 0<(c.get('component_order_bound') or 0)<=14]
    protected=[c for c in cases if c.get('component_order_bound') is None]
    high_bound=[c for c in cases if (c.get('component_order_bound') or 0)>=15]
    too_large=[c for c in high_bound if len(c['marked_component'])>15]
    eligible=[c for c in cases if (c.get('component_order_bound') or 0)>=15 and len(c['marked_component'])<=15]
    assert [len(x) for x in (low_bound,protected,high_bound,too_large,eligible)]==[1691,83,840,304,536]
    assert Counter(c['method'] for c in protected)=={'finite_leaf9_or_ideal_feature5':74,'nine_context_adaptive_finite_payment':9}
    assert len(low_bound)+len(protected)+len(too_large)+len(eligible)==2614
    routing=dict(total=2614,bound_at_most14=1691,protected_or_adaptive=83,
                 bound_at_least15=840,displayed_larger_than15=304,eligible_rank15=536)
    paths=[Path(__file__),Path(R.__file__),Path(N.__file__),Path(X.__file__),Path(E.__file__),
           Path(S.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(P.__file__),Path(Q.__file__),Path(D.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return enc(dict(status='PASS_FINITE_RANK15_EXISTING_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        rules=rules,new_costs=dict(A15='104',B15='104',D15='104'),complete_selected_support=supports,
        selected_proofs=proofs,selected_contexts=6,radical_contexts=2,Lorentz_whole_component_contexts=4,
        complete_high_overlap_routing=routing,
        literal_configurations=216,literal_changed_corners=dict(hits),literal_rows=literal,
        simultaneous_rank16_exchange=docs[E.TAG]['new_costs'],rank16_literal_changes=dict(exchange_hits),
        preserved_D16_alternative_pairs=[[1,14],[1,15]],preserved_D16_alternative_weight='1/2',
        reused_A2_A13_coupon=coupon,new_cross_component_coupons=[],
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Complete frozen536-pattern native maps for each finite rank15 type, six arbitrary-support whole-core transfers, and all216 simultaneous current literal interfaces including the frozen rank16/A17 exchange. All nonselected/universal/current ideal routes are addressed in the written proof. No new coupon, finite optimization or dimension theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],'contexts',out['selected_contexts'],'literal',out['literal_configurations'],flush=True)
