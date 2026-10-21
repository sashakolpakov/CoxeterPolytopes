#!/usr/bin/env python3
"""Complete retained-profile diagnostic for deep finite distance six.

This pilot does not prove a new payment rule. It records every Schur
branch, including nonpositive branches and the required smaller-core
and current-payment interfaces that still need a universal proof.
"""
import argparse,json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as S
import fresh_depth8_finite7_profile_transfer_20261018 as P
import fresh_depth6_finite7_adaptive_transfer_20261018 as OLD

ROOT=P.ROOT;TAG='fresh_depth6_finite6_adaptive_pilot_20261018'
CURRENT='fresh_order15_d16_d17_ambient59_frontier_20261018'
CAT=OLD.CAT
N,A,B=19,6,12
SQ=S.sqrt(2)
enc=OLD.enc


def profiles():
    white=[i for i in range(N) if i not in(A,B)]
    allowed=[i for i in white if not(i>B and i-B<=4 or A<i<B and B-i<=4)]
    rows=[];kept=[]
    for digits in product((0,1,2),repeat=len(allowed)):
        if sum(bool(x) for x in digits)>2:continue
        x={i:S.sqrt(q) for i,q in zip(allowed,digits) if q}
        edge=P.projection(x,white);adj=P.projection(x,[i for i in range(N) if i!=A])
        ok=P.field(2-edge)['sign']>0 and P.field(2-adj)['sign']>=0
        rows.append(dict(contacts={i:str(v) for i,v in x.items()},edge=P.field(edge),adjacent=P.field(adj),accepted=ok))
        if ok:kept.append(x)
    assert len(rows)==1+2*len(allowed)+2*len(allowed)*(len(allowed)-1)
    low=[]
    for i,x in enumerate(kept):
        q=P.projection(x,list(range(N)))
        if P.field(q-2)['sign']<=0:low.append(dict(profile=i,contacts={j:str(v) for j,v in x.items()},projection=P.field(q)))
    assert min(P.inv(len(block),i,i) for block in (list(range(6)),list(range(5))) for i in range(len(block)))==S.Rational(5,6)
    return kept,dict(candidate_count=len(rows),accepted_count=len(kept),allowed=allowed,candidates=rows,
        profiles=[{i:str(v) for i,v in x.items()} for x in kept],low_projection_profiles=low,
        smallest_white_inverse_diagonal='5/6',three_contacts_projection_lower='5/2',label_at_least5_projection_strict_lower='25/12')


def small_cores(cat):
    rows=[];survivors=[]
    for r in cat['diagrams']:
        n=r['order'];g={i:{} for i in range(n)}
        for i,j,m in r['edges']:g[i][j]=g[j][i]=m
        for a,c in product(r['finite_deletions'],repeat=2):
            if a==c:continue
            dd={a:0};todo=deque([a])
            while todo:
                i=todo.popleft()
                for j in g[i]:
                    if j!=c and j not in dd:dd[j]=dd[i]+1;todo.append(j)
            contacts=[i for i in dd if c in g[i]];assert contacts
            features=[i for i in dd if sum(j!=c for j in g[i])>=3 or any(m!=3 and j!=c for j,m in g[i].items())]
            fd=min((dd[i] for i in features),default=99);dist=min(dd[i]+1 for i in contacts)
            row=dict(diagram=r['index'],order=n,a=a,c=c,component_distances=dd,
                component_contact_nodes=contacts,distance=dist,branch_or_double_depth=fd)
            rows.append(row)
            if fd>=6 and dist>=6:survivors.append(row)
    assert [(r['diagram'],r['a'],r['c']) for r in survivors]==[(25,8,1),(25,8,2)]
    return dict(ordered_pairs=rows,survivors=survivors,
        scope='Exact catalogue screen under depth>=6 and distance>=6. Derivation of this distance condition and both adaptive patch interfaces is a separate proof obligation.')


def schur(ps):
    t,u,T,U=S.symbols('t u T U',nonnegative=True)
    rows=[];bad=[];excluded=[];short=[];positive=[]
    def inner(x,y,n,shift):
        return S.expand(sum(v*w*P.inv(n,i+shift,j+shift) for i,v in x.items() for j,w in y.items()))
    for ci,c in enumerate(ps):
        for di,rr in enumerate(ps):
            d={N-1-i:q for i,q in rr.items()}
            paths=[(i,j) for i in c for j in d if abs(i-j)+2<=5 and not any(min(i,j)<k<max(i,j) for k in(A,B))]
            if paths:short.append(dict(c=ci,d=di,witness=paths[0]));continue
            x={**c,A:t};y={**d,B:u};oq=inner(x,x,N,0);orr=inner(y,y,N,0)
            for tc,ud in product((0,1),repeat=2):
                clow=P.field(oq.subs(t,tc)-2)['sign']<=0;dlow=P.field(orr.subs(u,ud)-2)['sign']<=0
                if clow and (tc or c!={1:1}) or dlow and (ud or d!={N-2:1}):
                    excluded.append(dict(c=ci,d=di,t=tc,u=ud,c_low=clow,d_low=dlow));continue
                for cp,dp in product(('extend','isolated') if clow else ('none',),('extend','isolated') if dlow else ('none',)):
                    shift=int(cp=='extend');n=N+shift+int(dp=='extend')
                    q=inner(x,x,n,shift)+S.Rational(int(cp=='isolated'),2)
                    r=inner(y,y,n,shift)+S.Rational(int(dp=='isolated'),2)
                    ss=inner(x,y,n,shift);q0=P.field(q.subs(t,tc)-2);r0=P.field(r.subs(u,ud)-2)
                    assert q0['sign']>0 and r0['sign']>0
                    poly=S.Poly(S.expand(((q-2)*(r-2)-ss**2).subs({t:1+T if tc else 0,u:1+U if ud else 0})),T,U)
                    coeff=[dict(powers=list(m),value=P.field(v)) for m,v in poly.terms()]
                    const=P.field(poly.eval({T:0,U:0}))
                    row=dict(c=ci,d=di,t=tc,u=ud,c_patch=cp,d_patch=dp,
                        q_minus2_lower=q0,r_minus2_lower=r0,coefficients=coeff,constant=const)
                    rows.append(row)
                    if const['sign']>0 and all(z['value']['sign']>=0 for z in coeff):positive.append(row)
                    else:bad.append(row)
    return dict(branch_count=len(rows),strict_branches=len(positive),remaining_branches=bad,
        branches=rows,conditional_small_core_exclusions=excluded,short_profile_pairs=short,
        scope='Unbounded opposite products0 or1+nonnegative parameter. Adaptive exclusions remain conditional until their geometry is proved.')


def central_extension_diagnostic():
    c,d=N,N+1;g=S.diag(P.path(N),S.diag(2,2));g[A+1,c]=g[c,A+1]=g[B-1,d]=g[d,B-1]=-1
    inv=g.inv();left=list(range(A+1))+[c];right=list(range(B+4,N));tail=inv.extract(right,right)
    k=inv.extract(left,left)-inv.extract(left,right)*tail.inv()*inv.extract(right,left)
    assert all(tail[:i,:i].det()>0 for i in range(1,len(right)+1)) and all(x>0 for x in k)
    assert k[-1,-1]==S.Rational(13,18)<2
    # A finite opposite corner with both marks deep would have outward
    # simple arms. Their hypothetical exact principals have two forks.
    arms=[]
    for left_arm,right_arm in product(range(1,8),repeat=2):
        # central path has five whites; two old arms have seven roots.
        h=P.path(5);h[0,0]-=S.Rational(7,8)+S.Rational(left_arm,left_arm+1)
        h[-1,-1]-=S.Rational(7,8)+S.Rational(right_arm,right_arm+1)
        inner=h.extract([1,2,3],[1,2,3]);ends=[0,4]
        z=h.extract(ends,ends)-h.extract(ends,[1,2,3])*inner.inv()*h.extract([1,2,3],ends)
        sig='two_negative' if z[0,0]<0 and z.det()>0 else 'not_excluded'
        arms.append(dict(c_arm_roots=left_arm,d_arm_roots=right_arm,
            schur=[[str(x) for x in row] for row in z.tolist()],classification=sig))
    return dict(gram=[[str(x) for x in row] for row in g.tolist()],
        free_contact_projection_lower=str(k[-1,-1]),single_endpoint_seal_fails=True,
        two_outward_arms=arms,
        scope='Hypothetical induced arm diagrams only. Existence/disjointness/contact exclusions under arbitrary extra whites and retained opposite-corner payment remain unproved.')


def run():
    stems=(P.TAG,OLD.TAG,CAT,CURRENT);docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for s,d in docs.items():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    ps,pr=profiles();sc=schur(ps)
    paths={Path(__file__),Path(P.__file__),Path(OLD.__file__)}
    paths.update(ROOT/f'{f}/{s}.{e}' for s in stems for f,e in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_DEPTH6_FINITE6_ADAPTIVE_DIAGNOSTIC',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        retained_path=dict(n=N,a=A,b=B),profiles=pr,minimal_core_screen=small_cores(docs[CAT]),schur=sc,
        central_extension=central_extension_diagnostic(),universal_payment_rule_proved=False,
        finite_forest_count_performed=False,new_dimension_bound=None,
        scope='Complete finite retained profile and conditional adaptive Schur diagnostic; outstanding universal and current-payment steps are not asserted.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run()
    p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'profiles',out['profiles']['accepted_count'],'branches',out['schur']['branch_count'],
        'remaining',len(out['schur']['remaining_branches']),flush=True)
    for row in out['schur']['remaining_branches']:print(row,flush=True)
