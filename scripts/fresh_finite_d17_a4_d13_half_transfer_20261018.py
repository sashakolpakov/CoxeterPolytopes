#!/usr/bin/env python3
"""Finite D17 cost120 ->118 using existing whole A4/D13 half coupons.

Reconstruct the complete selected maps and simultaneous current literal
debits. Arbitrary-support and whole-face coverage are proved in the report.
No finite-forest optimization or new ambient bound.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d17_five_unit_reduction_20261018 as N
import fresh_finite_rank15_existing_half_transfer_20261018 as R

ROOT=N.ROOT; T=N.T; A=N.A; Q=N.Q; D=N.D; S=N.S; P=N.P
TAG='fresh_finite_d17_a4_d13_half_transfer_20261018'
CURRENT='fresh_rank15_rank16_dimension61_20261018'
NEIGHBOR='fresh_finite8_a4_d13_neighbor_coupon_budget_20261017'
HALF={(3,13),(3,14),(3,15),(3,16)}
CASES=[1293,1944,2019,2032,2041,2049,2337,2508]
ROLES={1293:'fork_center',1944:'short_fork_leaf',2019:'fork_center',
       2032:'short_fork_leaf',2041:'short_fork_leaf',2049:'long_arm_neighbor',
       2337:'short_fork_leaf',2508:'long_arm_neighbor'}
STEMS=tuple(dict.fromkeys((*R.STEMS,R.TAG,N.TAG,CURRENT,NEIGHBOR)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def key(row):return row['source'],N.normalize(row['case'])


def seal(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert c['method']=='radical' and c['core_key']=='(25, 25, ((9, 9),))'
    assert n==19 and core['gauge']==[1]*n and Q.inertia(g)==[17,1,1]
    z=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    black=set(u+v);white=sorted(set(range(n))-black)
    assert len(white)==15 and Q.inertia(A.sub(g,white))==[15,0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0] for b in black)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[k]==c[k]
    assert len(c['marked_component'])==17
    images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],17)
    assert images and set(images)<=HALF
    allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2
            assert length<=8
            short.append(dict(nodes=[i,j],length=length))
    assert short
    corners=[];coupons=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[17,0,0]
        blocks=D.finite_blocks(g,ids)
        aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
        selected=marks==c['preferred_corner'];matches=[]
        if not selected and aa is not bb:
            for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if sm['type']!='A4' or la['type']!='D13':continue
                si=T.mark_pairs(g,sm['nodes'],[m],4,'A')
                li=T.mark_pairs(g,la['nodes'],[lm],13)
                assert si==[(0,),(3,)]
                role={((9,),):'long_arm_neighbor',((10,),):'fork_center',
                      ((11,),(12,)):'short_fork_leaf'}[tuple(li)]
                assert role==ROLES[c['index']]
                coupon=dict(small_type='whole A4',small_nodes=sm['nodes'],small_mark=m,
                    small_canonical_images=si,large_type='whole D13',large_nodes=la['nodes'],
                    large_mark=lm,large_canonical_images=li,large_role=role,
                    source=NEIGHBOR if role=='long_arm_neighbor' else N.COUPON,weight='1/2')
                matches.append(coupon);coupons.append(coupon)
        assert len(matches)<=1
        value=F(1,2) if selected or matches else F()
        corners.append(dict(marks=marks,components=blocks,selected=selected,
                            existing_coupons=matches,payment=value))
    assert len(coupons)==1 and sum(r['payment'] for r in corners)==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],radical=z,
        white=white,selected_images=images,allowed_contacts=allowed,mixed_sign_shortcuts=short,
        corners=corners,total='1',whole_components_sealed_unless_short8=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[CURRENT]['dimension_bound']==61
    assert docs[R.TAG]['new_costs']==dict(A15='104',B15='104',D15='104')
    assert docs[R.TAG]['simultaneous_rank16_exchange']==dict(A16='113',D16='225/2',A17='115')
    assert docs[NEIGHBOR]['component_count']['half_weight']=='1/2'
    assert docs[NEIGHBOR]['component_count']['cost_per_A4_D13_component_pair']=='1'
    assert 'endpoint of whole A4' in docs[NEIGHBOR]['additional_rule']
    assert 'non-leaf neighbor' in docs[NEIGHBOR]['additional_rule']
    assert [r['rank'] for r in docs[N.COUPON]['small_components']]==[1,2,3,4,5]
    assert next(r for r in docs[N.COUPON]['old_large_components'] if r['type']=='D13')['eligible_marks']==3
    old=docs[N.TAG];weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==120
    assert all(weights[p]==1 for p in HALF)
    new={p:(F(1,2) if p in HALF else w) for p,w in weights.items()}
    assert sum(new.values())==118
    target=T.target(17);adj,dist,_=P.distance_data(target);autos=P.embeddings(target,target)
    assert target==old['target'] and all({tuple(sorted(mp[i] for i in p)) for p in HALF}==HALF for mp in autos)
    assert all(dist[i][j]>=10 for i,j in HALF)
    assert all(new[tuple(p)]==1 for p in old['retained_basic_pairs'])
    over=docs[A.OVER];support=T.selected(over,17)
    assert enc(support)==old['complete_selected_support']
    assert support['complete_eligible_patterns']==140 and support['complete_maps']==168
    hits=[r for r in support['rows'] if set(map(tuple,r['marked_images']))&HALF]
    assert [r['index'] for r in hits]==CASES
    assert all(set(map(tuple,r['marked_images']))<=HALF and len(r['base_nodes'])==17 for r in hits)
    proofs=[seal(over['isolation']['cases'][i],over['isolation']['cores'][over['isolation']['cases'][i]['core_key']]) for i in CASES]
    latest={key(r):r for r in docs[R.TAG]['literal_rows']}
    assert len(latest)==216
    literal=[];affected=[];roles=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=latest[stem,N.normalize(index)];assert len(cs)==len(prev['corners'])
        corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks']
            ids=sorted(set(white)|set(marks));im=[];debit=F()
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                shared=next((b for b in D.finite_blocks(g,ids) if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type']=='D17':
                    im=T.mark_pairs(g,shared['nodes'],marks,17)
                    delta={weights[p]-new[p] for p in im};assert len(delta)==1
                    debit=delta.pop()
            value=F(prev['corners'][j]['retained_payment'])-debit
            assert value>=0
            corners.append(dict(marks=marks,whole_D17_images=im,old_payment=prev['corners'][j]['retained_payment'],
                debit=debit,retained_payment=value))
            if debit:roles[str(im)]+=1
        total=sum(r['retained_payment'] for r in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total=total);literal.append(row)
        if any(r['debit'] for r in corners):affected.append(row)
    assert len(literal)==216 and len(affected)==24
    assert Counter(r['source'] for r in affected)=={'fresh_affine_d13_long_pair_transfer_20261017':24}
    # The two intervening finite transfers use no D17 same-component role.
    for p in docs[R.TAG]['selected_proofs']+docs[R.E.TAG]['selected_proofs']:
        for c in p['corners']:
            if F(c['payment']):
                assert not any(b['type']=='D17' and set(c['marks'])<=set(b['nodes']) for b in c['components'])
    paths=[Path(__file__),Path(N.__file__),Path(R.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(P.__file__),Path(D.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return enc(dict(status='PASS_FINITE_D17_A4_D13_EXISTING_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        target=target,automorphisms=autos,half_pairs=sorted(HALF),old_cost='120',new_cost='118',total_saving='2',
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(new.items())],
        retained_basic_pairs=old['retained_basic_pairs'],complete_selected_support=support,
        selected_contexts=len(proofs),affected_maps=sum(len(r['maps']) for r in hits),selected_proofs=proofs,
        selected_coupon_roles=dict(Counter(ROLES.values())),literal_configurations=216,
        literal_changed_configurations=len(affected),literal_changed_roles=dict(roles),literal_rows=literal,
        literal_affected=affected,new_cross_component_coupons=[],all_other_pair_weights_unchanged=True,
        current_rank15_costs=docs[R.TAG]['new_costs'],current_rank16_exchange=docs[R.E.TAG]['new_costs'],
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Complete140 selected patterns/168 maps, eight arbitrary-support sealed radical contexts and all216 current literal configurations after the simultaneous rank15/rank16 exchange. Reuses existing wholeA4/D13 neighbor, fork-center and short-leaf halves. No new coupon, finite-forest optimization or ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'cost',out['new_cost'],'contexts',out['selected_contexts'],
          'maps',out['affected_maps'],'literal',out['literal_changed_configurations'],flush=True)
