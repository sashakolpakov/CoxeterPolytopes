#!/usr/bin/env python3
"""Exact depth6/7 finite7 profile diagnostic, not yet a universal theorem.

Central-contact reduction and non-strict surviving branches are open inputs.
This finite screen records every candidate; no dimension claim follows.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import sympy as S
import fresh_depth8_finite7_profile_transfer_20261018 as P

ROOT=P.ROOT
TAG='fresh_depth6_depth7_finite7_profile_pilot_20261018'
CURRENT='fresh_d15_balanced_ambient60_frontier_20261018'
SQ=S.sqrt(2)

def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x

def profiles(h):
    n=2*h+8;a=h;b=h+7;white=[i for i in range(n)if i not in(a,b)]
    allowed=[]
    for i in white:
        if i>b and i-b<=5 or a<i<b and b-i<=5:continue
        if i<a and P.field(2-P.inv(h,i,i))['sign']<=0:continue
        if i>b and P.field(2-P.inv(h,i-b-1,i-b-1))['sign']<=0:continue
        allowed.append(i)
    corner=[i for i in range(n)if i!=a]
    rows=[];kept=[]
    for digits in product((0,1,2),repeat=len(allowed)):
        if sum(bool(x)for x in digits)>2:continue
        x={i:S.sqrt(q)for i,q in zip(allowed,digits)if q}
        edge=P.projection(x,white);adj=P.projection(x,corner)
        ok=P.field(2-edge)['sign']>0 and P.field(2-adj)['sign']>=0
        rows.append(dict(contacts={str(i):str(q)for i,q in x.items()},
            edge_projection=P.field(edge),corner_projection=P.field(adj),accepted=ok))
        if ok:kept.append(x)
    assert len(rows)==1+2*len(allowed)+2*len(allowed)*(len(allowed)-1)
    return kept,dict(depth=h,n=n,a=a,b=b,allowed=allowed,candidates=rows,
        candidate_count=len(rows),accepted_profiles=[{str(i):str(q)for i,q in x.items()}for x in kept])

def branches(h,ps):
    n=2*h+8;a=h;b=h+7
    t,u,T,U=S.symbols('t u T U',nonnegative=True)
    J=S.Matrix(n,n,lambda i,j:P.inv(n,i,j))
    rows=[];short=[];outer=[];problem=[];positive=0
    for ci,c in enumerate(ps):
        for di,rr in enumerate(ps):
            d={n-1-i:q for i,q in rr.items()}
            shortcut=[(i,j)for i in c for j in d if abs(i-j)+2<=6
                      and not any(min(i,j)<k<max(i,j)for k in(a,b))]
            if shortcut:
                short.append(dict(c=ci,d=di,witness=shortcut[0]));continue
            x=S.Matrix([c.get(i,0)+(t if i==a else 0)for i in range(n)])
            y=S.Matrix([d.get(i,0)+(u if i==b else 0)for i in range(n)])
            q=S.expand((x.T*J*x)[0]);r=S.expand((y.T*J*y)[0]);s=S.expand((x.T*J*y)[0])
            det=S.expand((q-2)*(r-2)-s*s)
            for tc,ud in product((0,1),repeat=2):
                q0=q.subs(t,tc);r0=r.subs(u,ud)
                if P.field(q0-2)['sign']<=0 or P.field(r0-2)['sign']<=0:
                    outer.append(dict(c=ci,d=di,t=tc,u=ud,q=P.field(q0),r=P.field(r0)));continue
                poly=S.Poly(S.expand(det.subs({t:1+T if tc else 0,u:1+U if ud else 0})),T,U)
                coeff=[dict(powers=list(m),value=P.field(v))for m,v in poly.terms()]
                constant=P.field(poly.eval({T:0,U:0}))
                row=dict(c=ci,d=di,t=tc,u=ud,q=P.field(q0),r=P.field(r0),coefficients=coeff,constant=constant)
                rows.append(row)
                if constant['sign']>0 and all(z['value']['sign']>=0for z in coeff):positive+=1
                else:problem.append(row)
    return dict(strict_branches=positive,remaining_branches=problem,all_branches=rows,
        short_profile_pairs=short,unproved_outer_or_zero_projection_branches=outer)

def central_extension(h):
    n=2*h+8;a=h;b=h+7;c=n;d=n+1
    g=S.diag(P.path(n),S.Matrix([[2,0],[0,2]]))
    g[a+1,c]=g[c,a+1]=-1;g[b-1,d]=g[d,b-1]=-1
    j=g.inv();left=list(range(a+1))+[c];right=list(range(b+5,n))
    tail=j.extract(right,right)
    schur=j.extract(left,left)-j.extract(left,right)*tail.inv()*j.extract(right,left)
    pivots=[tail[:i,:i].det()for i in range(1,len(right)+1)]
    assert all(P.field(x)['sign']>0for x in pivots)
    assert all(P.field(x)['sign']>0for x in schur)
    assert P.field(schur[-1,-1]-2)['sign']>0
    q,r,s=(P.inv(n,i,k)for i,k in((a+1,a+1),(b-1,b-1),(a+1,b-1)))
    two=S.Matrix([[2-q,-s],[-s,2-r]])
    assert two.det()<0
    return dict(gram=[[str(x)for x in row]for row in g.tolist()],
        full_inertia=[n+1,1,0],two_root_schur=[[str(x)for x in row]for row in two.tolist()],
        no_short6_contact_domain_if_extra_white_touches_c=dict(left=left,right=right,
            excluded_middle_cycle_positions=list(range(a+1,b))),
        eliminated_tail_pivots=list(map(str,pivots)),remaining_schur=[[str(x)for x in row]for row in schur.tolist()],
        forced_unit_c_projection_lower=str(schur[-1,-1]),
        scope='Conditional linear-algebra seal of mark c after excluding middle contacts by positive-edge acyclicity and contacts near adjacent b,d by short6. Reflection gives d. This diagnostic does not establish all other branches or the central-contact lemma.')

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in(P.TAG,CURRENT)}
    for s,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[CURRENT]['finite']['exact_maximum']=='815/2'
    rows=[]
    for h in(6,7):
        ps,dom=profiles(h);br=branches(h,ps);ext=central_extension(h)
        rows.append(dict(depth=h,domain=dom,schur=br,central_extension=ext))
        print('DEPTH_PROFILE',h,'profiles',len(ps),'strict',br['strict_branches'],
              'remaining',len(br['remaining_branches']),'outer',len(br['unproved_outer_or_zero_projection_branches']),
              'endpoint_lower',ext['forced_unit_c_projection_lower'],flush=True)
    paths={Path(__file__),Path(P.__file__),Path(P.D.B.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}'for s in(P.TAG,CURRENT)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_DEPTH6_DEPTH7_PROFILE_AND_CONDITIONAL_ENDPOINT_PILOT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        rows=rows,new_universal_payment_rule=False,new_forest_count=False,new_dimension_bound=None,
        scope='Complete native retained-white profile screens and exact polynomial branch diagnostics for two proposed depths; conditional central-endpoint projection seals. No universal central-contact reduction, all surviving-branch proof, new cost rule or dimension theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
