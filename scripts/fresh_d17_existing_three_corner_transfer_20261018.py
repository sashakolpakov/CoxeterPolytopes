#!/usr/bin/env python3
"""Two whole-D17 reductions using three already assigned corner payments.

No new coupon or ideal change. Complete selected support, four radical
seals, and all current literal interfaces are checked exactly.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as X
import fresh_finite_d17_five_unit_reduction_20261018 as N
import fresh_center_reversion_rank16_balance_ambient60_frontier_20261018 as B
import fresh_finite_a16_endpoint_sparse_half_transfer_20261018 as E

ROOT=X.ROOT;T=X.T;A=X.A;Q=X.Q;D=X.D;P=N.P
TAG='fresh_d17_existing_three_corner_transfer_20261018'
CHANGES={(2,12):F(5,8),(2,14):F()}
DEBITS={(2,12):F(3,8),(2,14):F(1,2)}
CASES=[1017,1062,1264,2525]
STEMS=tuple(dict.fromkeys((*X.STEMS,X.TAG,N.TAG,B.TAG,E.TAG)))
enc=X.R.enc


def seal(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert n==19 and core['gauge']==[1]*n and c['method']=='radical'
    assert c['core_key']=='(25, 25, ((9, 9),))' and Q.inertia(g)==[17,1,1]
    z=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][j] for j in c['right'][1]]
    white=sorted(set(range(n))-set(u+v));assert len(white)==15
    assert Q.inertia(A.sub(g,white))==[15,0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0] for b in u+v)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[k]==c[k]
    images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],17,'D')
    assert len(c['marked_component'])==17 and len(images)==1 and images[0] in CHANGES
    allowed=c['possible_contacts'];assert all(z[j] for j in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2
            assert length<=8;short.append(dict(nodes=[i,j],length=length))
    assert short
    corners=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[17,0,0]
        blocks=D.finite_blocks(g,ids)
        aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
        payment=F();reason=None;roles=[]
        if marks==c['preferred_corner']:
            payment=CHANGES[images[0]];reason='changed_selected_D17'
        elif aa is not bb:
            for small,large,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if small['type']=='A1' and large['type']=='D5':
                    im=T.mark_pairs(g,large['nodes'],[lm],5,'D')
                    if im==[(0,)]:
                        payment=F(3,8);reason='existing_D5_long_A1';roles=im
                if images[0]==(2,14) and small['type']=='A1' and large['type']=='D14':
                    im=T.mark_pairs(g,large['nodes'],[lm],14,'D')
                    if im==[(11,)]:
                        payment=F(1,2);reason='existing_A1_D14_center';roles=im
                if images[0]==(2,14) and small['type']=='A2' and large['type']=='A14':
                    im=T.mark_pairs(g,large['nodes'],[lm],14,'A')
                    if im==[(2,),(11,)]:
                        payment=F(1,2);reason='existing_A2_A14_third';roles=im
        corners.append(dict(marks=marks,components=blocks,payment=payment,reason=reason,canonical_roles=roles))
    assert sum(r['payment'] for r in corners)==1
    assert sum(bool(r['payment']) for r in corners)==2
    assert sorted(r['payment'] for r in corners)==([F(),F(),F(3,8),F(5,8)] if images[0]==(2,12) else [F(),F(),F(1,2),F(1,2)])
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],radical=z,
        white=white,selected_images=images,allowed_contacts=allowed,mixed_sign_shortcuts=short,
        corners=corners,total='1',whole_core_sealed_unless_independent_short8=True)


def literal(docs):
    before={X.R.key(r['source'],r['case']):r for r in docs[X.TAG]['combined_literal_interface']['rows']}
    rows=[];hits=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=before[X.R.key(stem,index)];corners=[]
        assert len(cs)==len(prev['corners'])
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks));images=[];debit=F()
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                b=next((b for b in D.finite_blocks(g,ids) if b['type']=='D17' and set(marks)<=set(b['nodes'])),None)
                if b:
                    images=T.mark_pairs(g,b['nodes'],marks,17,'D')
                    values={DEBITS.get(p,F()) for p in images};assert len(values)==1;debit=values.pop()
            payment=F(prev['corners'][j]['retained_payment_lower_bound'])-debit;assert payment>=0
            corners.append(dict(marks=marks,D17_images=images,debit=debit,retained_payment_lower_bound=payment))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(r['debit'] for r in corners):hits.append(row)
    assert len(rows)==len(before)==216 and not hits
    return dict(configurations=216,changed=hits,rows=rows)


def actual_base(t):return X.actual_base(t)-(14 if t==('D',17) else 0)
def upper_base(t):return X.upper_base(t)-(14 if t==('D',17) else 0)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    weights={tuple(r['pair']):F(r['weight']) for r in docs[B.TAG]['local_proofs']['finite_weighted_pairs']['D17']}
    assert sum(weights.values())==F(X.actual_base(('D',17)),16)==F(467,4)
    for p,w in CHANGES.items():assert weights[p]-w==DEBITS[p];weights[p]=w
    assert sum(weights.values())==F(actual_base(('D',17)),16)==F(927,8)
    assert docs[X.TAG]['finite']['component_costs']['A1_D5']=='3/8'
    assert docs[E.TAG]['new_cross_component_coupons'][1]['weight']=='1/2'
    target=T.target(17);adj,dist,_=P.distance_data(target);autos=P.embeddings(target,target)
    assert all(dist[a][b]>=10 for a,b in CHANGES)
    assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES}==set(CHANGES) for mp in autos)
    over=docs[A.OVER];sup=T.selected(over,17)
    assert enc(sup)==docs[N.TAG]['complete_selected_support']
    assert sup['complete_eligible_patterns']==140 and sup['complete_maps']==168
    hits=[r for r in sup['rows'] if set(map(tuple,r['marked_images']))&CHANGES.keys()]
    assert [r['index'] for r in hits]==CASES
    assert all(len(r['base_nodes'])==17 and set(map(tuple,r['marked_images']))<=CHANGES.keys() for r in hits)
    proofs=[seal(over['isolation']['cases'][i],over['isolation']['cores'][over['isolation']['cases'][i]['core_key']]) for i in CASES]
    # Current positive same-component later balances consume the disjoint
    # D17 orbit at first node1 and short leaves15/16, not these node2 roles.
    assert not set(CHANGES)&B.PAIRS['D17']
    types=[t for t in X.P.Z.B.B.types() if t[1]<=59];rep=X.P.Z.B.B.representative
    assert all(upper_base(rep(t))>=actual_base(t)>=0 for t in types)
    witness=[('D',16),('D',16),('D',27)]
    attained=sum(actual_base(t) for t in witness)+sum(X.actual_cross(a,b) for i,a in enumerate(witness) for b in witness[i+1:])+X.R.context_cost(witness,False)
    assert sum(r for k,r in witness)==59 and F(attained,16)==F(docs[X.TAG]['finite']['exact_maximum'])==F(795,2)
    paths={Path(__file__),Path(X.__file__),Path(N.__file__),Path(B.__file__),Path(E.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_D17_EXISTING_THREE_CORNER_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=target,automorphisms=autos,changed_pairs=[dict(pair=p,new_weight=w,debit=DEBITS[p]) for p,w in CHANGES.items()],
        old_cost='467/4',new_cost='927/8',component_saving='7/8',
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        complete_selected_support=sup,selected_proofs=proofs,affected_maps=sum(len(r['maps']) for r in hits),
        combined_literal_interface=literal(docs),all_cross_component_costs_unchanged=True,ideal_weights_unchanged=True,
        triangle_minimum=2,quadrilateral_minimum=1,
        finite=dict(exact_maximum='795/2',curvature='783/2',actual_attainer=['D16','D16','D27'],
            proof='The frozen complete upper remains valid under pointwise decreases; the displayed actual formal forest is unchanged and attains it.',
            new_Bellman_run=False,representative_base_checks=len(types),remaining_alternatives_recounted=False),
        new_dimension_bound=None,new_finite_volume_example=False,
        scope='Fresh complete140-pattern/168-map finiteD17 selected reconstruction, four arbitrary-support radical seals, and all216 current literal interfaces. Existing whole-component coupons only; complete earlier geometry and current independent short/ideal routes are inherited. No new finite optimization or ideal replay.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_cost'],'maps',out['affected_maps'],'literal',out['combined_literal_interface']['configurations'],flush=True)
