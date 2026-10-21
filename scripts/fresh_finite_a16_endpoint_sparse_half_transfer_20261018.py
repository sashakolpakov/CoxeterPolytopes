#!/usr/bin/env python3
"""Finite A16 cost111 with two explicit sparse whole-component coupons."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d17_a4_d13_half_transfer_20261018 as E

R=E.R;V=R.E;T=E.T;A=E.A;Q=E.Q;D=E.D;P=E.P;S=E.S;ROOT=E.ROOT
TAG='fresh_finite_a16_endpoint_sparse_half_transfer_20261018'
DEPTH='fresh_depth8_finite8_zero_transfer_20261018'
HALF={(0,13),(2,15),(1,13),(2,14)}
CASES=[1013,1028,1138,1611]
STEMS=tuple(dict.fromkeys((*E.STEMS,E.TAG,DEPTH)))
COUPONS=[dict(name='A1_A13end',small='whole A1',small_nodes=[0],large='whole A13',large_nodes=[0,12],weight='1/2',component_pair_cost='1'),
         dict(name='A2_A14third',small='whole A2',small_nodes=[0,1],large='whole A14',large_nodes=[2,11],weight='1/2',component_pair_cost='2')]


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def seal(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert n==19 and core['gauge']==[1]*n and c['method']=='radical'
    assert Q.inertia(g)==[17,1,1]
    z=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    white=sorted(set(range(n))-set(u+v))
    assert len(white)==15 and Q.inertia(A.sub(g,white))==[15,0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0] for b in u+v)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    assert len(c['marked_component'])==16
    images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],16,'A')
    assert images and set(images)<=HALF
    allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
            short.append(dict(nodes=[i,j],length=length))
    assert short
    wanted='A1_A13end' if c['index'] in (1028,1611) else 'A2_A14third'
    spec=next(x for x in COUPONS if x['name']==wanted)
    corners=[];coupons=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[17,0,0]
        blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
        value=F();reason=None;roles=[]
        if marks==c['preferred_corner']:
            value=F(1,2);reason='selected_A16_half'
        elif aa is not bb:
            for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if 'whole '+sm['type']==spec['small'] and 'whole '+la['type']==spec['large']:
                    si=T.mark_pairs(g,sm['nodes'],[m],len(sm['nodes']),'A')
                    li=T.mark_pairs(g,la['nodes'],[lm],len(la['nodes']),'A')
                    assert si==[(i,) for i in spec['small_nodes']] and li==[(i,) for i in spec['large_nodes']]
                    assert not value
                    value=F(1,2);reason=wanted
                    roles=[S.marked_role(g,sm,m),S.marked_role(g,la,lm)]
                    coupons.append(dict(name=wanted,small_images=si,large_images=li,weight='1/2'))
        corners.append(dict(marks=marks,components=blocks,payment=value,reason=reason,roles=roles))
    assert len(coupons)==1 and sum(r['payment'] for r in corners)==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],radical=z,
        white=white,selected_images=images,allowed_contacts=allowed,mixed_sign_shortcuts=short,
        corners=corners,coupons=coupons,total='1',whole_core_sealed_unless_short8=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[E.TAG]['new_cost']=='118'
    assert docs[DEPTH]['feature_depth']==8 and docs[DEPTH]['quadrilateral_minimum']==1
    assert docs[V.TAG]['new_costs']==dict(A16='113',D16='225/2',A17='115')
    rec=next(r for r in docs[S.TAG]['rules'] if r['kind']=='A' and r['rank']==16)
    old={tuple(r['pair']):F(r['weight']) for r in rec['weighted_pairs']}
    assert sum(old.values())==114 and all(old[p]==1 for p in V.AHALF)
    old.update({p:F(1,2) for p in V.AHALF});assert sum(old.values())==113
    assert all(old[p]==1 for p in HALF)
    new={p:(F(1,2) if p in HALF else w) for p,w in old.items()}
    assert sum(new.values())==111
    target=rec['target'];adj,dist,_=P.distance_data(target)
    assert all(dist[i][j]>=10 for i,j in HALF)
    assert all(new[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    assert all({tuple(sorted(mp[i] for i in p)) for p in HALF}==HALF for mp in rec['automorphisms'])
    sup=V.support('A',docs[T.TAG]);over=docs[A.OVER]
    assert V.validate_maps(sup,target,over)==238
    affected=[r for r in sup['rows'] if set(map(tuple,r['marked_images']))&HALF]
    assert [r['index'] for r in affected]==CASES
    assert all(set(map(tuple,r['marked_images']))<=HALF and len(r['base_nodes'])==16 for r in affected)
    proofs=[seal(over['isolation']['cases'][i],over['isolation']['cores'][over['isolation']['cases'][i]['core_key']]) for i in CASES]
    for c in COUPONS:
        assert F(c['component_pair_cost'])==len(c['small_nodes'])*len(c['large_nodes'])*F(c['weight'])
    oldcoupons=docs[S.TAG]['coupons']
    assert next(c for c in oldcoupons if c['name']=='A2_A14second')['large_nodes']==[1,12]
    assert set(COUPONS[1]['large_nodes']).isdisjoint({1,12})
    # The older broad A13 role is its unique middle, disjoint from endpoints.
    assert set(COUPONS[0]['large_nodes']).isdisjoint({6})
    previous={(r['source'],E.N.normalize(r['case'])):r for r in docs[E.TAG]['literal_rows']}
    assert len(previous)==216
    literal=[];hits=0
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=previous[stem,E.N.normalize(index)];rows=[]
        assert len(cs)==len(prev['corners'])
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();images=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                b=next((b for b in D.finite_blocks(g,ids) if set(marks)<=set(b['nodes'])),None)
                if b and b['type']=='A16':
                    images=T.mark_pairs(g,b['nodes'],marks,16,'A')
                    changes={old[p]-new[p] for p in images};assert len(changes)==1;debit=changes.pop()
            value=F(prev['corners'][j]['retained_payment'])-debit;assert value>=0;hits+=bool(debit)
            rows.append(dict(marks=marks,whole_A16_images=images,debit=debit,retained_payment=value))
        total=sum(r['retained_payment'] for r in rows);assert total>=1
        literal.append(dict(source=stem,case=index,corners=rows,total=total))
    assert len(literal)==216 and hits==0
    for p in docs[R.TAG]['selected_proofs']+docs[V.TAG]['selected_proofs']+docs[E.TAG]['selected_proofs']:
        for c in p['corners']:
            if not F(c['payment']):continue
            for b in c['components']:
                if b['type']=='A16' and set(c['marks'])<=set(b['nodes']):
                    images=T.mark_pairs(A.matrix(p['gram']),b['nodes'],c['marks'],16,'A')
                    assert not set(images)&HALF
    paths=[Path(__file__),Path(E.__file__),Path(R.__file__),Path(V.__file__),Path(S.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)]
    paths +=[ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return enc(dict(status='PASS_FINITE_A16_ENDPOINT_SPARSE_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        old_cost='113',new_cost='111',target=target,automorphisms=rec['automorphisms'],half_pairs=sorted(HALF),
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(new.items())],retained_basic_pairs=rec['retained_basic_pairs'],
        complete_selected_support=sup,selected_contexts=4,affected_maps=sum(len(r['maps']) for r in affected),
        selected_proofs=proofs,new_cross_component_coupons=COUPONS,literal_rows=literal,literal_configurations=216,
        literal_debited_corners=hits,compatible_D17_cost='118',compatible_feature_depth=8,
        other_current_rank16_costs=dict(D16='225/2',A17='115'),all_other_finite_and_ideal_weights_unchanged=True,
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Complete saved394-pattern/238-map A16 domain, four full radical seals with exact whole-component roles, and all216 latest literal interfaces. Two new typed half coupons cost1 and2 per eligible component pair. Their budget must be included separately; no ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'A16',result['new_cost'],'contexts',result['selected_contexts'],'literal',result['literal_debited_corners'],flush=True)
