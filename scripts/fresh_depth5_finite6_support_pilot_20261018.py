#!/usr/bin/env python3
"""Diagnostic for finite distance6 with feature depth5, not a theorem.

Complete retained native profiles and minimal-core screen; conditional
positive-block patches and unbounded Schur tests. Arbitrary-support patch
geometry and current-payment retention are explicitly unproved here.
"""
import argparse,json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import sympy as S
import fresh_depth8_finite7_profile_transfer_20261018 as P
import fresh_depth6_finite6_adaptive_pilot_20261018 as Z
import fresh_a13_endpoint_a1_context_retention_20261018 as CURRENT

ROOT=P.ROOT;TAG='fresh_depth5_finite6_support_pilot_20261018'
N,A,B=17,5,11
enc=Z.enc

def profiles():
    white=[i for i in range(N)if i not in(A,B)]
    allowed=[i for i in white if not(i>B and i-B<=4 or A<i<B and B-i<=4)]
    assert allowed==[0,1,2,3,4,6,16]
    rows=[];kept=[]
    for digits in product((0,1,2),repeat=len(allowed)):
        if sum(bool(q)for q in digits)>2:continue
        x={i:S.sqrt(q)for i,q in zip(allowed,digits)if q}
        edge=P.projection(x,white);adj=P.projection(x,[i for i in range(N)if i!=A])
        ok=P.field(2-edge)['sign']>0 and P.field(2-adj)['sign']>=0
        rows.append(dict(contacts={i:str(v)for i,v in x.items()},edge=P.field(edge),adjacent=P.field(adj),accepted=ok))
        if ok:kept.append(x)
    assert len(rows)==99 and len(kept)==17
    return kept,dict(candidate_count=99,accepted_count=17,allowed=allowed,candidates=rows,
        profiles=[{i:str(v)for i,v in x.items()}for x in kept],
        full_path_projections=[P.field(P.projection(x,list(range(N))))for x in kept],
        minimum_white_inverse_diagonal='5/6',three_contact_lower='5/2',large_label_strict_lower='25/12')

def minimal(cat):
    rows=[];survivors=[]
    for r in cat['diagrams']:
        n=r['order'];adj={i:{}for i in range(n)}
        for i,j,m in r['edges']:adj[i][j]=adj[j][i]=m
        for a,c in product(r['finite_deletions'],repeat=2):
            if a==c:continue
            dist={a:0};todo=deque([a])
            while todo:
                i=todo.popleft()
                for j in adj[i]:
                    if j!=c and j not in dist:dist[j]=dist[i]+1;todo.append(j)
            contacts=[i for i in dist if c in adj[i]]
            features=[i for i in dist if sum(j!=c for j in adj[i])>=3 or any(m!=3 and j!=c for j,m in adj[i].items())]
            fd=min((dist[i]for i in features),default=99);distance=min(dist[i]+1 for i in contacts)
            row=dict(diagram=r['index'],order=n,a=a,c=c,distance=distance,feature_depth=fd,contacts=contacts)
            rows.append(row)
            if fd>=5 and distance>=5:survivors.append(row)
    assert [(r['diagram'],r['a'],r['c'])for r in survivors]==[(25,7,1),(25,7,2),(25,8,0),(25,8,1),(25,8,2),(25,8,3)]
    return dict(ordered_pairs=rows,survivors=survivors,universal_distance_reduction_proved=False)

def positive_block(left,right):
    rows=[list(map(S.Rational,r))for r in P.path(N).tolist()];contacts=[{},{}]
    def attach(old):
        k=len(rows)
        for i,row in enumerate(rows):row.append(S.Rational(-int(i==old)))
        rows.append([S.Rational(-int(i==old))for i in range(k)]+[S.Rational(2)])
        return k
    def isolated():
        k=len(rows)
        for row in rows:row.append(S.Rational(0))
        rows.append([S.Rational(0)]*k+[S.Rational(2)]);return k
    for side,name in enumerate((left,right)):
        end=0 if side==0 else N-1
        if name.startswith('extend'):
            k=attach(end)
            if name=='extend2':attach(k)
            elif name=='extend1_isolated1':contacts[side][isolated()]=S.Rational(1)
            else:assert name=='extend1'
        elif name=='isolated1':contacts[side][isolated()]=S.Rational(1)
        elif name=='fork2':
            contacts[side][attach(end)]=S.Rational(1);attach(end)
        else:assert name=='none'
    g=S.Matrix(rows)
    determinants=[g[:i,:i].det()for i in range(1,len(rows)+1)]
    return g,contacts,determinants

def patches(profile,low):
    if not low:return('none',)
    if profile.get(1)==1:return('extend1','isolated1')
    if profile.get(0)==1:return('extend2','extend1_isolated1')
    if not profile.get(0,0)and not profile.get(1,0):return('fork2',)
    return()

def schur(ps):
    t,u,T,U=S.symbols('t u T U',nonnegative=True)
    rows=[];bad=[];strict=[];short=[];excluded=[];blocks={}
    for lp,rp in product(('none','extend1','isolated1','extend2','extend1_isolated1','fork2'),repeat=2):
        g,cs,piv=positive_block(lp,rp)
        blocks[lp,rp]=(g,cs,piv,g.inv()if all(z>0 for z in piv)else None)
    for ci,c in enumerate(ps):
        for di,rr in enumerate(ps):
            d={N-1-i:q for i,q in rr.items()}
            paths=[(i,j)for i in c for j in d if abs(i-j)+2<=5 and not any(min(i,j)<k<max(i,j)for k in(A,B))]
            if paths:short.append(dict(c=ci,d=di,witness=paths[0]));continue
            for tc,ud in product((0,1),repeat=2):
                clow=P.field(P.projection({**c,A:tc},list(range(N)))-2)['sign']<=0
                dlow=P.field(P.projection({**rr,A:ud},list(range(N)))-2)['sign']<=0
                if clow:assert tc==0
                if dlow:assert ud==0
                for lp,rp in product(patches(c,clow),patches(rr,dlow)):
                    g,extra,piv,J=blocks[lp,rp];key=dict(c=ci,d=di,t=tc,u=ud,c_patch=lp,d_patch=rp)
                    if J is None:
                        excluded.append(dict(**key,reason='selected_finite_block_not_positive',leading_determinants=list(map(str,piv))));continue
                    xc={**c,A:t,**extra[0]};yd={**d,B:u,**extra[1]}
                    x=S.Matrix([xc.get(i,0)for i in range(g.rows)]);y=S.Matrix([yd.get(i,0)for i in range(g.rows)])
                    edge_ids=[i for i in range(g.rows)if i not in(A,B)]
                    tests=[]
                    for vector,omit in((x,A),(y,B)):
                        edge=S.expand((vector.extract(edge_ids,[0]).T*g.extract(edge_ids,edge_ids).inv()*vector.extract(edge_ids,[0]))[0])
                        ids=[i for i in range(g.rows)if i!=omit]
                        adjacent=S.expand((vector.extract(ids,[0]).T*g.extract(ids,ids).inv()*vector.extract(ids,[0]))[0])
                        tests.append((P.field(2-edge),P.field(2-adjacent)))
                    if any(a['sign']<=0 or b['sign']<0 for a,b in tests):
                        excluded.append(dict(**key,reason='enlarged_edge_or_adjacent_corner_not_allowed',margins=tests));continue
                    q=S.expand((x.T*J*x)[0]);r=S.expand((y.T*J*y)[0]);s=S.expand((x.T*J*y)[0])
                    q0=P.field(q.subs(t,tc)-2);r0=P.field(r.subs(u,ud)-2)
                    if q0['sign']<=0 or r0['sign']<=0:
                        bad.append(dict(**key,reason='patch_projection_still_low',q_minus2=q0,r_minus2=r0));continue
                    poly=S.Poly(S.expand(((q-2)*(r-2)-s*s).subs({t:1+T if tc else 0,u:1+U if ud else 0})),T,U)
                    terms=[dict(powers=list(m),value=P.field(v))for m,v in poly.terms()]
                    constant=P.field(poly.eval({T:0,U:0}))
                    row=dict(**key,q_minus2=q0,r_minus2=r0,coefficients=terms,constant=constant)
                    rows.append(row)
                    if constant['sign']>0 and all(z['value']['sign']>=0 for z in terms):strict.append(row)
                    else:bad.append(row)
    return dict(conditional_branches=rows,branch_count=len(rows),strict_branches=len(strict),
        remaining_branches=bad,short_profile_pairs=short,conditionally_excluded=excluded,
        block_models=[dict(left=k[0],right=k[1],gram=[[str(x)for x in r]for r in z[0].tolist()],
            black_patch_contacts=[{i:str(v)for i,v in c.items()}for c in z[1]],
            positive_leading_determinants=list(map(str,z[2])),positive=z[3]is not None)for k,z in blocks.items()],
        arbitrary_support_patch_geometry_proved=False,current_payment_retention_proved=False)

def run():
    stems=(P.TAG,Z.CAT,CURRENT.TAG);docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for name,h in d['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(s,name)
    ps,pr=profiles();print('PROFILES',len(ps),flush=True)
    small=minimal(docs[Z.CAT]);print('MINIMAL_CORE_SURVIVORS',len(small['survivors']),flush=True)
    sc=schur(ps)
    paths={f'scripts/{TAG}.py',f'scripts/{Z.TAG}.py'}|{f'{folder}/{s}.{ext}'for s in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_DEPTH5_FINITE6_CONDITIONAL_SUPPORT_DIAGNOSTIC',
        sources={name:sha256((ROOT/name).read_bytes()).hexdigest()for name in sorted(paths)},
        retained_path=dict(n=N,a=A,b=B),profiles=pr,minimal_core_screen=small,schur=sc,
        conjectural_new_pair_counts=dict(A='max(0,n-16)',B='max(0,n-17)',D='max(0,n-18)'),
        universal_local_theorem_proved=False,new_dimension_bound=None,
        scope='Complete retained native profiles and minimal-core screen; conditional enlarged-block Schur diagnostic only. Patch completeness/arbitrary-contact geometry and current payments remain unproved.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'branches',out['schur']['branch_count'],'strict',out['schur']['strict_branches'],
        'remaining',len(out['schur']['remaining_branches']),flush=True)
