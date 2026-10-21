#!/usr/bin/env python3
"""Depth-six finite-seven transfer: complete profiles and adaptive small cores.

Published minimal-diagram classification and current short-six/triangle
payment theorems remain explicit inherited inputs. No forest optimization.
"""
import argparse,json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as S
import fresh_depth8_finite7_profile_transfer_20261018 as P

ROOT=P.ROOT
TAG='fresh_depth6_finite7_adaptive_transfer_20261018'
CAT='fresh_high_rank_quasilanner_census_20261017'
CURRENT='fresh_small_d16_coupon_and_ideal17_rebalance_20261018'
N,A,B=20,6,13
SQ=S.sqrt(2)

def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x

def profiles():
    white=[i for i in range(N) if i not in(A,B)]
    allowed=[i for i in white if not(i>B and i-B<=5 or A<i<B and B-i<=5)]
    assert allowed==[0,1,2,3,4,5,7,19]
    rows=[];kept=[]
    for digits in product((0,1,2),repeat=len(allowed)):
        if sum(bool(x) for x in digits)>2:continue
        x={i:S.sqrt(q) for i,q in zip(allowed,digits) if q}
        edge=P.projection(x,white);adj=P.projection(x,[i for i in range(N) if i!=A])
        ok=P.field(2-edge)['sign']>0 and P.field(2-adj)['sign']>=0
        rows.append(dict(contacts={i:str(v) for i,v in x.items()},edge=P.field(edge),adjacent=P.field(adj),accepted=ok))
        if ok:kept.append(x)
    assert len(rows)==129 and len(kept)==18
    assert kept==[{}, {19:1},{19:SQ},{7:1},{7:1,19:1},{7:SQ},
                 {5:1},{5:1,19:1},{5:1,7:1},{5:SQ},{4:1},{3:1},{2:1},
                 {1:1},{0:1},{0:1,19:1},{0:1,7:1},{0:SQ}]
    assert P.inv(6,0,0)==S.Rational(6,7)
    return kept,dict(candidates=rows,candidate_count=129,accepted_count=18,
        profiles=[{i:str(v) for i,v in x.items()} for x in kept],allowed=allowed,
        three_contacts_projection_lower='18/7',label_at_least5_projection_strict_lower='15/7')

def small_cores(cat):
    rows=[];survivors=[]
    for r in cat['diagrams']:
        n=r['order'];g={i:{} for i in range(n)}
        for i,j,m in r['edges']:g[i][j]=g[j][i]=m
        assert set(r['finite_deletions'])=={i for i,z in enumerate(r['deletions']) if all(t['kind']=='finite' for t in z)}
        for a,c in product(r['finite_deletions'],repeat=2):
            if a==c:continue
            dd={a:0};todo=deque([a])
            while todo:
                i=todo.popleft()
                for j in g[i]:
                    if j!=c and j not in dd:dd[j]=dd[i]+1;todo.append(j)
            contacts=[i for i in dd if c in g[i]]
            assert contacts
            dist=min(dd[i]+1 for i in contacts)
            features=[i for i in dd if sum(j!=c for j in g[i])>=3 or any(m!=3 and j!=c for j,m in g[i].items())]
            fd=min((dd[i] for i in features),default=99)
            row=dict(diagram=r['index'],order=n,a=a,c=c,component_distances=dd,
                component_contact_nodes=contacts,distance=dist,branch_or_double_depth=fd)
            rows.append(row)
            if fd>=6 and dist>=6:survivors.append(row)
    assert len(cat['diagrams'])==26
    assert [(r['diagram'],r['a'],r['c']) for r in survivors]==[(25,8,1),(25,8,2)]
    core=cat['diagrams'][25]
    assert core['order']==10 and core['index']==25
    kg=S.Matrix(10,10,lambda i,j:2 if i==j else 0)
    for i,j,m in core['edges']:kg[i,j]=kg[j,i]=-1 if m==3 else -SQ
    patches=[]
    for which,cm,mp in (
        ('extend',1,{0:1,1:20,2:0,3:22,4:2,5:3,6:4,7:5,8:6,9:7}),
        ('isolated',2,{0:1,1:0,2:20,3:22,4:2,5:3,6:4,7:5,8:6,9:7})):
        g=S.diag(P.path(20),S.diag(2,2,2));g[1,20]=g[20,1]=-1
        j=0 if which=='extend' else 20;g[j,22]=g[22,j]=-1
        assert g.extract([mp[i] for i in range(10)],[mp[i] for i in range(10)])==kg
        patches.append(dict(kind=which,c_node=cm,map_to_retained_path=mp,
            extra_white=22,extra_white_contact=j,contact_magnitude='1',old_path_contact_of_c=1))
    # Strict W+c excludes an isolated patch joining the opposite old endpoint.
    assert P.inv(6,1,1)==S.Rational(10,7)
    assert 1/(2-P.inv(6,5,5))==S.Rational(7,8)
    assert S.Rational(10,7)+S.Rational(7,8)>2
    return dict(all_ordered_finite_deletion_pairs=rows,survivors=survivors,patches=patches,
        catalogue_diagrams=26,lower_order_exclusion='A connected diagram on at most five nodes has distance at most four.',
        isolated_patch_forbidden_far_endpoint_projection='129/56',
        mutual_patch_orthogonality='A contact gives an adjacent c-to-d white path of length at most five unless both are selected-component extensions; then it gives a cycle in the selected finite component.')

def schur(ps):
    t,u,T,U=S.symbols('t u T U',nonnegative=True)
    rows=[];special=[];excluded=[];short=[];positive=[]
    def inner(x,y,n,shift):
        return S.expand(sum(v*w*P.inv(n,i+shift,j+shift) for i,v in x.items() for j,w in y.items()))
    for ci,c in enumerate(ps):
        for di,rr in enumerate(ps):
            d={19-i:q for i,q in rr.items()}
            paths=[(i,j) for i in c for j in d if abs(i-j)+2<=6 and not any(min(i,j)<k<max(i,j) for k in(A,B))]
            if paths:short.append(dict(c=ci,d=di,witness=paths[0]));continue
            x={**c,A:t};y={**d,B:u}
            oq=inner(x,x,20,0);orr=inner(y,y,20,0)
            for tc,ud in product((0,1),repeat=2):
                clow=P.field(oq.subs(t,tc)-2)['sign']<=0
                dlow=P.field(orr.subs(u,ud)-2)['sign']<=0
                if clow and (tc or c!={1:1}) or dlow and (ud or d!={18:1}):
                    excluded.append(dict(c=ci,d=di,t=tc,u=ud,c_low=clow,d_low=dlow));continue
                for cp,dp in product(('extend','isolated') if clow else ('none',),('extend','isolated') if dlow else ('none',)):
                    shift=int(cp=='extend');n=20+shift+int(dp=='extend')
                    q=inner(x,x,n,shift)+S.Rational(int(cp=='isolated'),2)
                    r=inner(y,y,n,shift)+S.Rational(int(dp=='isolated'),2)
                    ss=inner(x,y,n,shift)
                    q0=P.field(q.subs(t,tc)-2);r0=P.field(r.subs(u,ud)-2)
                    assert q0['sign']>0 and r0['sign']>0
                    poly=S.Poly(S.expand(((q-2)*(r-2)-ss**2).subs({t:1+T if tc else 0,u:1+U if ud else 0})),T,U)
                    coeff=[dict(powers=list(m),value=P.field(v)) for m,v in poly.terms()]
                    const=poly.eval({T:0,U:0})
                    row=dict(c=ci,d=di,t=tc,u=ud,c_patch=cp,d_patch=dp,
                        q_minus2_lower=q0,r_minus2_lower=r0,coefficients=coeff,constant=P.field(const))
                    rows.append(row)
                    if P.field(const)['sign']>0:
                        assert all(z['value']['sign']>=0 for z in coeff);positive.append(const)
                    else:
                        assert (ci,di,tc,ud,cp,dp)==(3,3,0,0,'none','none') and poly.as_expr()==-S.Rational(4,7)
                        special.append(row)
    assert len(rows)==766 and len(positive)==765 and len(special)==1 and len(excluded)==261
    least=positive[0]
    for v in positive[1:]:
        if P.field(v-least)['sign']<0:least=v
    return dict(branch_count=len(rows),positive_branches=len(positive),minimum_positive_constant=P.field(least),
        branches=rows,central_exception=special,small_core_excluded=excluded,short_profile_pairs=short,
        arbitrary_opposite_products='Each product magnitude is zero or 1 plus an arbitrary nonnegative parameter.')

def endpoint_seal():
    c,d=20,21
    g=S.diag(P.path(20),S.diag(2,2));g[7,c]=g[c,7]=g[12,d]=g[d,12]=-1
    j=g.inv();left=list(range(7))+[c];right=[18,19]
    tail=j.extract(right,right)
    k=j.extract(left,left)-j.extract(left,right)*tail.inv()*j.extract(right,left)
    pivots=[tail[:i,:i].det() for i in range(1,3)]
    assert all(x>0 for x in pivots) and all(x>0 for x in k) and k[-1,-1]==S.Rational(5,2)
    q=P.inv(20,7,7);s=P.inv(20,7,12)
    block=S.Matrix([[2-q,-s],[-s,2-q]])
    assert block.det()==-S.Rational(4,7)
    return dict(gram=[[str(x) for x in row] for row in g.tolist()],inertia=[21,1,0],
        two_root_schur=[[str(x) for x in row] for row in block.tolist()],
        allowed_contact_domain_for_extra_white_touching_c=dict(left=left,right=right),
        eliminated_tail_positive_leading_determinants=list(map(str,pivots)),
        remaining_schur=[[str(x) for x in row] for row in k.tolist()],
        norm2_projection_lower='5/2',both_marked_degrees=1,
        retained_whole_alternatives=['finite_A8','finite_D9','affine_D9_order10','affine_E8_order9'])

def counts():
    rows=[]
    for kind in 'ABD':
        for n in range(1 if kind=='A' else 2 if kind=='B' else 4,81):
            edges=[(i,i+1,3) for i in range(n-1)] if kind!='D' else [(i,i+1,3) for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)]
            if kind=='B':edges[-1]=(n-2,n-1,4)
            adj={i:[] for i in range(n)}
            for i,j,m in edges:adj[i].append(j);adj[j].append(i)
            features={i for i in adj if len(adj[i])!=2}|{i for a,b,m in edges if m!=3 for i in(a,b)}
            dist=[]
            for i in adj:
                dd={i:0};todo=deque([i])
                while todo:
                    v=todo.popleft()
                    for z in adj[v]:
                        if z not in dd:dd[z]=dd[v]+1;todo.append(z)
                dist.append(dd)
            depth=[min(dist[i][f] for f in features) for i in adj]
            pairs=[(i,j) for i,j in combinations(range(n),2) if dist[i][j]==7 and min(depth[i],depth[j])>=6]
            assert len(pairs)==max(0,n-{'A':19,'B':20,'D':21}[kind])
            rows.append(dict(kind=kind,rank=n,pairs=pairs))
    return dict(formulas=dict(A='max(0,n-19)',B='max(0,n-20)',D='max(0,n-21)'),rows=rows)

def run():
    stems=(P.TAG,CAT,CURRENT,P.D.TAG,P.D.B.TAG)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[CURRENT]['finite']['above_curvature']==6
    paths={Path(__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in stems for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    ps,pr=profiles()
    return enc(dict(status='PASS_DEPTH6_FINITE7_ADAPTIVE_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        profiles=pr,minimal_core_reduction=small_cores(docs[CAT]),schur=schur(ps),endpoint_seal=endpoint_seal(),counts=counts(),
        triangle_minimum=2,quadrilateral_minimum=1,new_coupons=False,new_forest_count=False,new_dimension_bound=None,
        scope='Universal local depth-six transfer with two classified adaptive supporting-white patches, an exhaustive unbounded Schur certificate and arbitrary-contact endpoint seals. Current short-six/triangle composition and minimal-diagram classification are explicit proof inputs. No finite-forest or global dimension computation.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'branches',out['schur']['branch_count'],'minimum',out['schur']['minimum_positive_constant'],flush=True)
