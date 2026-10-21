#!/usr/bin/env python3
"""Bounded A17 alternate-corner and accounting pilot; no new theorem."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as X
import fresh_center_reversion_rank16_balance_ambient60_frontier_20261018 as B

ROOT=X.ROOT;T=X.T;A=X.A;D=X.D;Q=X.Q;P=X.R.P
TAG='fresh_a17_alternate_corner_pilot_20261018'
CHANGE={(2,12):F(),(4,14):F(),(2,13):F(1,2),(3,14):F(1,2),(2,14):F()}
STEMS=tuple(dict.fromkeys((*X.STEMS,X.TAG,B.TAG)))
enc=X.R.enc


def additions(g,blocks,marks):
    aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks];out=[]
    if aa is bb:return out
    for sm,la,m,n in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if sm['type']=='D5' and la['type']=='A12':
            si=T.mark_pairs(g,sm['nodes'],[m],5,'D');li=T.mark_pairs(g,la['nodes'],[n],12,'A')
            if si==[(3,),(4,)] and li==[(2,),(9,)]:out.append(dict(kind='D5_short_A12_third',weight=F(1),small_images=si,large_images=li))
        if sm['type']=='A1' and la['type']=='A14':
            li=T.mark_pairs(g,la['nodes'],[n],14,'A')
            if li in ([(2,),(11,)],[(3,),(10,)]):out.append(dict(kind='A1_A14_third_or_fourth',weight=F(1,2),large_images=li))
    assert len(out)<=1;return out


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    weights={tuple(r['pair']):F(r['weight']) for r in docs[B.TAG]['local_proofs']['finite_weighted_pairs']['A17']}
    assert sum(weights.values())==F(X.actual_base(('A',17)),16)==F(467,4)
    assert all(weights[p]==1 for p in CHANGE)
    assert not set(CHANGE)&{(0,13),(3,16),(1,12),(4,15),(1,15)}
    target=T.target(17,'A');_,dist,_=T.P.distance_data(target);autos=T.P.embeddings(target,target)
    assert all(dist[a][b]>=10 for a,b in CHANGE)
    assert all({tuple(sorted(m[i] for i in p)) for p in CHANGE}==set(CHANGE) for m in autos)
    over=docs[A.OVER];sup=T.selected(over,17,'A')
    hits=[r for r in sup['rows'] if set(map(tuple,r['marked_images']))&CHANGE.keys()]
    assert [r['index'] for r in hits]==[1009,1040,1044,2011,2057]
    proofs=[]
    for row in hits:
        c=over['isolation']['cases'][row['index']];core=over['isolation']['cores'][c['core_key']];g=A.matrix(core['gram']);n=len(g)
        assert n==19 and c['method']=='radical' and Q.inertia(g)==[17,1,1]
        z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
        u=c['left'][1];v=[c['right_copy_map'][j] for j in c['right'][1]];white=sorted(set(range(n))-set(u+v))
        assert Q.inertia(A.sub(g,white))==[15,0,0]
        assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0] for b in u+v)
        actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
        for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
        assert len(c['marked_component'])==17
        allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
        du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8;short.append(dict(nodes=[i,j],length=length))
        assert short;corners=[]
        for old in c['all_corners']:
            marks=old['preferred_corner'];ids=sorted(set(white)|set(marks));assert Q.inertia(A.sub(g,ids))==[17,0,0]
            bs=D.finite_blocks(g,ids);cp=additions(g,bs,marks);value=sum((r['weight'] for r in cp),F());images=[]
            if marks==c['preferred_corner']:
                images=T.mark_pairs(g,c['marked_component'],marks,17,'A');values={CHANGE[p] for p in images};assert len(values)==1;value+=values.pop()
            corners.append(dict(marks=marks,components=bs,new_coupons=cp,selected_images=images,hypothetical_payment=value))
        assert sum(r['hypothetical_payment'] for r in corners)==1
        proofs.append(dict(index=c['index'],gram=core['gram'],radical=z,white=white,allowed_contacts=allowed,
            mixed_sign_shortcuts=short,corners=corners,hypothetical_total='1'))
    previous={X.R.key(r['source'],r['case']):r for r in docs[X.TAG]['combined_literal_interface']['rows']}
    literal=[];changed=[];bad=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=previous[X.R.key(stem,index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();cp=[];images=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                bs=D.finite_blocks(g,ids);cp=additions(g,bs,marks)
                b=next((b for b in bs if b['type']=='A17' and set(marks)<=set(b['nodes'])),None)
                if b:
                    images=T.mark_pairs(g,b['nodes'],marks,17,'A');values={weights[p]-CHANGE.get(p,weights[p]) for p in images};assert len(values)==1;debit=values.pop()
            inc=sum((r['weight'] for r in cp),F());value=F(prev['corners'][j]['retained_payment_lower_bound'])-debit+inc
            corners.append(dict(marks=marks,A17_images=images,debit=debit,new_coupons=cp,increment=inc,hypothetical_payment=value))
        total=sum(r['hypothetical_payment'] for r in corners);r=dict(source=stem,case=index,corners=corners,hypothetical_total=total);literal.append(r)
        if any(c['debit'] or c['increment'] for c in corners):changed.append(r)
        if total<1 or any(c['hypothetical_payment']<0 for c in corners):bad.append(r)
    assert len(literal)==216
    rep=X.P.Z.B.B.representative;types=[t for t in X.P.Z.B.B.types() if t[1]<=59];fail=[]
    for a in types:
        for b in types:
            pair={a,b};inc=64 if pair=={('D',5),('A',12)} else 32 if pair=={('A',1),('A',14)} else 0
            actual=X.actual_cross(a,b)+inc;upper=X.upper_cross(rep(a),rep(b))
            if actual>upper:fail.append(dict(pair=[a,b],actual=actual,upper=upper))
    assert X.actual_cross(('A',1),('A',14))==0
    saving=sum(weights[p]-w for p,w in CHANGE.items());assert saving==4
    paths={Path(__file__),Path(X.__file__),Path(B.__file__),Path(T.__file__),Path(A.__file__),Path(D.__file__),Path(Q.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_A17_ALTERNATE_CORNER_PILOT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        proposed_pairs=[dict(pair=p,weight=w) for p,w in CHANGE.items()],hypothetical_component_saving=saving,
        old_A17_cost='467/4',hypothetical_A17_cost='451/4',selected_support=sup,selected_core_diagnostics=proofs,
        literal_records=literal,literal_changed=changed,literal_failures=bad,
        simultaneous_type_pairs=len(types)**2,domination_failures=fail,
        proposed_actual_costs=dict(D5_A12='27/4',A1_A14='2'),
        new_local_payment_theorem=False,finite_budget_performed=False,new_dimension_bound=None,
        scope='Complete selectedA17 map reconstruction and five exact radical/core-role diagnostics,216 current literal trials and full simultaneous cross-type cost check. A complete later-output/universal composition audit and independent mathematical review remain required. No ambient theorem or geometric realization.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'saving',out['hypothetical_component_saving'],'literal_failures',len(out['literal_failures']),
        'type_failures',len(out['domination_failures']),flush=True)
