#!/usr/bin/env python3
"""Independent SymPy reconstruction of the 69 supporting-extension contexts.

The literal context list is a supplied finite domain, not an independent
exhaustion of all high-core overlaps. No producer arithmetic is imported.
"""
import argparse
from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance9_supporting_extensions_review_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
SEEDS=[[[18, [1, 3]], [25, [1, 2]], [[6, 9], [7, 8]], [1, 9]], [[18, [1, 3]], [25, [1, 3]], [[6, 9], [7, 8]], [1, 9]], [[18, [1, 3]], [25, [2, 3]], [[6, 9], [7, 8]], [1, 10]], [[22, [1, 2]], [25, [1, 2]], [[6, 9], [7, 8], [8, 7]], [1, 10]], [[22, [1, 2]], [25, [1, 3]], [[6, 9], [7, 8], [8, 7]], [1, 10]], [[22, [1, 2]], [25, [2, 3]], [[6, 9], [7, 8], [8, 7]], [1, 11]], [[22, [1, 3]], [25, [1, 2]], [[6, 9], [7, 8], [8, 7]], [1, 10]], [[22, [1, 3]], [25, [1, 3]], [[6, 9], [7, 8], [8, 7]], [1, 10]], [[22, [1, 3]], [25, [2, 3]], [[6, 9], [7, 8], [8, 7]], [1, 11]], [[22, [1, 5]], [25, [1, 2]], [[3, 9], [4, 8]], [1, 10]], [[22, [1, 5]], [25, [1, 3]], [[3, 9], [4, 8]], [1, 10]], [[22, [1, 5]], [25, [2, 3]], [[3, 9], [4, 8]], [1, 11]], [[22, [1, 6]], [25, [1, 2]], [[3, 9], [4, 8]], [1, 10]], [[22, [1, 6]], [25, [1, 3]], [[3, 9], [4, 8]], [1, 10]], [[22, [1, 6]], [25, [2, 3]], [[3, 9], [4, 8]], [1, 11]], [[22, [2, 3]], [25, [1, 2]], [[6, 9], [7, 8], [8, 7]], [2, 10]], [[22, [2, 3]], [25, [1, 3]], [[6, 9], [7, 8], [8, 7]], [2, 10]], [[22, [2, 3]], [25, [2, 3]], [[6, 9], [7, 8], [8, 7]], [2, 11]], [[22, [5, 6]], [25, [1, 2]], [[3, 9], [4, 8]], [5, 10]], [[22, [5, 6]], [25, [1, 3]], [[3, 9], [4, 8]], [5, 10]], [[22, [5, 6]], [25, [2, 3]], [[3, 9], [4, 8]], [5, 11]], [[25, [1, 2]], [18, [1, 3]], [[8, 7], [9, 6]], [1, 11]], [[25, [1, 2]], [22, [1, 2]], [[7, 8], [8, 7], [9, 6]], [1, 11]], [[25, [1, 2]], [22, [1, 3]], [[7, 8], [8, 7], [9, 6]], [1, 11]], [[25, [1, 2]], [22, [1, 5]], [[8, 4], [9, 3]], [1, 11]], [[25, [1, 2]], [22, [1, 6]], [[8, 4], [9, 3]], [1, 11]], [[25, [1, 2]], [22, [2, 3]], [[7, 8], [8, 7], [9, 6]], [1, 12]], [[25, [1, 2]], [22, [5, 6]], [[8, 4], [9, 3]], [1, 13]], [[25, [1, 2]], [25, [1, 2]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [1, 11]], [[25, [1, 2]], [25, [1, 3]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [1, 11]], [[25, [1, 2]], [25, [1, 4]], [[9, 3]], [1, 11]], [[25, [1, 2]], [25, [1, 5]], [[9, 3]], [1, 11]], [[25, [1, 2]], [25, [2, 3]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [1, 12]], [[25, [1, 2]], [25, [4, 5]], [[9, 3]], [1, 13]], [[25, [1, 3]], [18, [1, 3]], [[8, 7], [9, 6]], [1, 11]], [[25, [1, 3]], [22, [1, 2]], [[7, 8], [8, 7], [9, 6]], [1, 11]], [[25, [1, 3]], [22, [1, 3]], [[7, 8], [8, 7], [9, 6]], [1, 11]], [[25, [1, 3]], [22, [1, 5]], [[8, 4], [9, 3]], [1, 11]], [[25, [1, 3]], [22, [1, 6]], [[8, 4], [9, 3]], [1, 11]], [[25, [1, 3]], [22, [2, 3]], [[7, 8], [8, 7], [9, 6]], [1, 12]], [[25, [1, 3]], [22, [5, 6]], [[8, 4], [9, 3]], [1, 13]], [[25, [1, 3]], [25, [1, 2]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [1, 11]], [[25, [1, 3]], [25, [1, 3]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [1, 11]], [[25, [1, 3]], [25, [1, 4]], [[9, 3]], [1, 11]], [[25, [1, 3]], [25, [1, 5]], [[9, 3]], [1, 11]], [[25, [1, 3]], [25, [2, 3]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [1, 12]], [[25, [1, 3]], [25, [4, 5]], [[9, 3]], [1, 13]], [[25, [1, 4]], [25, [1, 2]], [[3, 9]], [1, 11]], [[25, [1, 4]], [25, [1, 3]], [[3, 9]], [1, 11]], [[25, [1, 4]], [25, [2, 3]], [[3, 9]], [1, 12]], [[25, [1, 5]], [25, [1, 2]], [[3, 9]], [1, 11]], [[25, [1, 5]], [25, [1, 3]], [[3, 9]], [1, 11]], [[25, [1, 5]], [25, [2, 3]], [[3, 9]], [1, 12]], [[25, [2, 3]], [18, [1, 3]], [[8, 7], [9, 6]], [2, 11]], [[25, [2, 3]], [22, [1, 2]], [[7, 8], [8, 7], [9, 6]], [2, 11]], [[25, [2, 3]], [22, [1, 3]], [[7, 8], [8, 7], [9, 6]], [2, 11]], [[25, [2, 3]], [22, [1, 5]], [[8, 4], [9, 3]], [2, 11]], [[25, [2, 3]], [22, [1, 6]], [[8, 4], [9, 3]], [2, 11]], [[25, [2, 3]], [22, [2, 3]], [[7, 8], [8, 7], [9, 6]], [2, 12]], [[25, [2, 3]], [22, [5, 6]], [[8, 4], [9, 3]], [2, 13]], [[25, [2, 3]], [25, [1, 2]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [2, 11]], [[25, [2, 3]], [25, [1, 3]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [2, 11]], [[25, [2, 3]], [25, [1, 4]], [[9, 3]], [2, 11]], [[25, [2, 3]], [25, [1, 5]], [[9, 3]], [2, 11]], [[25, [2, 3]], [25, [2, 3]], [[5, 9], [6, 8], [7, 7], [8, 6], [9, 5]], [2, 12]], [[25, [2, 3]], [25, [4, 5]], [[9, 3]], [2, 13]], [[25, [4, 5]], [25, [1, 2]], [[3, 9]], [4, 11]], [[25, [4, 5]], [25, [1, 3]], [[3, 9]], [4, 11]], [[25, [4, 5]], [25, [2, 3]], [[3, 9]], [4, 12]]]


def make(seed,ds):
    left,right,match,preferred=seed
    a,b=ds[left[0]],ds[right[0]];mp={j:i for i,j in match};n=a['order']
    for j in range(b['order']):
        if j not in mp:mp[j]=n;n+=1
    edges={}
    for i,j,m in a['edges']:edges[tuple(sorted((i,j)))]=m
    for i,j,m in b['edges']:
        ij=tuple(sorted((mp[i],mp[j])))
        assert ij not in edges or edges[ij]==m
        edges[ij]=m
    G=2*s.eye(n)
    for (i,j),m in edges.items():
        assert m==3 # every one of these exceptional full cores is simply laced
        G[i,j]=G[j,i]=-1
    return G,mp


def components(ids,G):
    remaining=set(ids);out=[]
    while remaining:
        c={min(remaining)};todo=list(c)
        while todo:
            i=todo.pop()
            for j in remaining-c:
                if G[i,j]:c.add(j);todo.append(j)
        remaining-=c;out.append(sorted(c))
    return out


def distances(G,marks):
    d={i:0 for i in marks};todo=deque(marks)
    while todo:
        i=todo.popleft()
        for j in range(G.rows):
            if i!=j and G[i,j] and j not in d:d[j]=d[i]+1;todo.append(j)
    assert len(d)==G.rows
    return d


def matrix_out(G):return [[str(x)for x in row]for row in G.tolist()]


def inertia(G):
    # Independent exact symmetric congruence; all full cores nonsingular.
    a=G.copy();p=n=z=0
    while a.rows:
        k=next((i for i in range(a.rows)if a[i,i]),None)
        if k is not None:
            ids=[k]+[i for i in range(a.rows)if i!=k];a=a.extract(ids,ids);v=a[0,0]
            p+=1 if v>0 else 0;n+=1 if v<0 else 0
            a=a[1:,1:]-a[1:,0]*a[0,1:]/v
        else:
            ij=next(((i,j)for i in range(a.rows)for j in range(i+1,a.rows)if a[i,j]),None)
            if ij is None:z+=a.rows;break
            i,j=ij;ids=[i,j]+[k for k in range(a.rows)if k not in ij];a=a.extract(ids,ids)
            p+=1;n+=1;a=a[2:,2:]-a[2:,:2]*a[:2,:2].inv()*a[:2,2:]
    return [p,n,z]


def run():
    cp=ROOT/f'results/{CENSUS}.json';ds={d['index']:d for d in json.loads(cp.read_text())['diagrams']}
    rows=[];cores={};types=Counter();residuals=Counter();bounds=Counter();gcount=0
    assert len(SEEDS)==69 and len({json.dumps(seed)for seed in SEEDS})==69
    for seed in SEEDS:
        left,right,match,pref=seed;G,mp=make(seed,ds);J=G.inv();n=G.rows
        corekey=json.dumps([left[0],right[0],match]);sig=inertia(G);assert sig==[n-1,1,0]
        if corekey not in cores:cores[corekey]=dict(gram=matrix_out(G),inverse=matrix_out(J),signature=sig)
        black=left[1]+[mp[i]for i in right[1]];omitted=sorted(set(black)-set(pref))
        ids=[i for i in range(n)if i not in omitted];allowed=list(omitted);marked=None;parts=[]
        for part in components(ids,G):
            H=G.extract(part,part)
            if H.det():
                assert H.is_positive_definite
                Q=H.inv();assert all(t>=0 for t in Q)
                vals=[Q[j,j]for j in range(len(part))]
                allowed.extend(i for j,i in enumerate(part)if vals[j]<=2)
                parts.append(dict(nodes=part,kind='finite',inverse_diagonals=list(map(str,vals))))
            else:
                ns=H.nullspace();assert len(ns)==1;z=ns[0]
                if all(t<0 for t in z):z=-z
                assert all(t>0 for t in z) and H.is_positive_semidefinite
                parts.append(dict(nodes=part,kind='affine',positive_kernel=list(map(str,z))))
            if set(pref)<=set(part):marked=part
        assert marked is not None;allowed=sorted(allowed)
        du,dv=distances(G,left[1]),distances(G,[mp[i]for i in right[1]])
        assert min(du[i]for i in [mp[t]for t in right[1]])>9
        choices=[]
        for radius in range(8):
            L=[i for i in allowed if du[i]<=radius];R=[i for i in allowed if dv[i]<=7-radius]
            if set(L+R)==set(allowed)and not set(L)&set(R)and all(J[i,j]>0 for end in(L,R)for i in end for j in end):choices.append((L,R))
        assert choices;L,R=choices[0]
        pairbounds=[dict(nodes=[i,j],minimum_projection=str(J[i,i]+2*J[i,j]+J[j,j]))for end in(L,R)for i,j in combinations(end,2)]
        assert all(s.Rational(r['minimum_projection'])>2 for r in pairbounds)
        candidates=[i for i in allowed if J[i,i]<=2]
        assert all(J[i,i]>=1 for i in candidates)
        relevant=[i for i in candidates if i in marked]
        assert all(J[i,i]>1 for i in relevant)
        extensions=[]
        for i in relevant:
            p=J[i,i];rho=2-p;assert 0<=rho<1
            alone=int(s.floor(rho/(1-rho)));coexisting=[];tests=[]
            for j in candidates:
                rr=2-J[j,j];P=J[i,j];dist=min(du[i]+dv[j],du[j]+dv[i])+3
                rec=dict(node=j,projection=str(P),residual_product=str(rho*rr),shortcut_bound=dist)
                if P>0 and P*P>rho*rr:rec['kind']='positive_actual_product'
                elif P<0 and P*P>rho*rr and dist<=9:rec['kind']='forced_short_path'
                elif P>0 and P*P<=rho*rr<(P+1)**2:
                    assert j not in relevant;rec['kind']='orthogonal_by_discrete_gap';coexisting.append(j)
                else:raise AssertionError((seed,i,j,P,rho,rr))
                tests.append(rec)
            assert len(coexisting)<=1
            doubled=[]
            for j in candidates:
                if J[j,j]!=1:continue
                P=J[i,j];dist=min(du[i]+dv[j],du[j]+dv[i])+3
                assert P>0 or(P<0 and dist<=9)
                doubled.append(dict(node=j,projection=str(P),kind='positive_actual_product'if P>0 else'forced_short_path',shortcut_bound=dist))
            joint=None
            if coexisting:
                j=coexisting[0];rr=2-J[j,j];P=J[i,j];S=s.Matrix([[rho,-P],[-P,rr]])
                assert 0<rr<=1 and P>0 and S.det()>=0
                if not S.det():
                    ns=S.nullspace();assert len(ns)==1;z=ns[0]
                    if all(t<0 for t in z):z=-z
                    assert all(t>0 for t in z)
                    joint=dict(outside_node=j,residual_gram=matrix_out(S),determinant='0',positive_kernel=list(map(str,z)),kind='positive_residual_relation',tail_bound=0)
                else:
                    SI=S.inv();assert all(t>=0 for t in SI) and SI[0,0]>=2
                    joint=dict(outside_node=j,residual_gram=matrix_out(S),determinant=str(S.det()),inverse=matrix_out(SI),kind='no_tail'if SI[0,0]>2 else'saturated_single_tail',tail_bound=0 if SI[0,0]>2 else 1)
                residuals[(str(p),str(J[j,j]),str(P),str(S.det()),None if not S.det()else str(S.inv()[0,0]))]+=1
            extensions.append(dict(node=i,projection=str(p),rho=str(rho),alone_tail_bound=alone,pair_tests=tests,doubled_contact_tests=doubled,coexisting_outside=joint))
        bound=len(marked)+max([1+max(e['alone_tail_bound'],(e['coexisting_outside']or{}).get('tail_bound',0))for e in extensions],default=0)
        assert bound<=16;gcount+=len(extensions);bounds[bound]+=1;types[corekey]+=1
        rows.append(dict(left=left,right=right,match=match,preferred_corner=pref,omitted_black=omitted,blocks=parts,possible_contacts=allowed,marked_component=marked,end_sets=[L,R],multicontact_tests=pairbounds,contact_domains=candidates,relevant_nodes=relevant,extensions=extensions,component_bound=bound))
    paths=[Path(__file__),cp]
    return dict(status='PASS_INDEPENDENT_DISTANCE9_SUPPORTING_EXTENSION_CONTEXTS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(contexts=len(rows),core_types=len(cores),relevant_root_cases=gcount,no_relevant_root_contexts=sum(not r['relevant_nodes']for r in rows),max_marked_component_order=max(bounds)),component_bounds={str(k):v for k,v in sorted(bounds.items())},residual_types=[dict(data=list(k),count=v)for k,v in sorted(residuals.items(),key=lambda kv:str(kv[0]))],cores=cores,contexts=rows,scope='Independent exact local Gram reconstruction and analytical supporting-extension bounds on the supplied69context domain; no complete high-overlap or distance9 dimension classification.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],d['counts'],flush=True)
