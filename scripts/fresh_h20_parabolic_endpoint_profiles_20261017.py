#!/usr/bin/env python3
"""Exact one-root parabolic endpoint profiles on the three old E8 partner edges.

Affine Coxeter classification is the completeness input. No wall heights are
chosen and no coherent chamber or finite-volume completion is claimed here.
"""
import argparse
import hashlib
import itertools
import json
import networkx as nx
import sympy as s
import fresh_h20_e8_partner_endpoints_20261017 as P

ROOT=P.ROOT
TAG='fresh_h20_parabolic_endpoint_profiles_20261017'
END='results/fresh_h20_e8_partner_endpoints_20261017.json'

def graph(n,edges):
    g=nx.Graph();g.add_nodes_from(range(n))
    g.add_edges_from((i,j,dict(label=m))for i,j,m in edges)
    return g

def star(arms):
    edges=[];v=1
    for arm in arms:
        last=0
        for _ in range(arm):edges.append((last,v,3));last=v;v+=1
    return graph(v,edges)

def affine_models():
    yield 'affine_A1',graph(2,[(0,1,'parallel')])
    for r in range(2,20):
        yield f'affine_A{r}',graph(r+1,[(i,(i+1)%(r+1),3)for i in range(r+1)])
        yield f'affine_C{r}',graph(r+1,[(i,i+1,4 if i in(0,r-1) else 3)for i in range(r)])
    for r in range(3,20):
        yield f'affine_B{r}',graph(r+1,[(0,2,3),(1,2,3)]+[(i,i+1,4 if i==r-1 else 3)for i in range(2,r)])
    for r in range(4,20):
        yield f'affine_D{r}',graph(r+1,[(0,2,3),(1,2,3)]+[(i,i+1,3)for i in range(2,r-2)]+[(r-2,r-1,3),(r-2,r,3)])
    for r,arms in [(6,(2,2,2)),(7,(1,3,3)),(8,(1,2,5))]:yield f'affine_E{r}',star(arms)
    yield 'affine_F4',graph(5,[(0,1,3),(1,2,3),(2,3,4),(3,4,3)])
    yield 'affine_G2',graph(3,[(0,1,3),(1,2,6)])

def magnitude(label):return s.Integer(2) if label=='parallel' else 2*s.cos(s.pi/label)

def run():
    D=json.loads((ROOT/END).read_text());out=[];models=list(affine_models())
    for ci,cen in enumerate(D['endpoint_censuses']):
        W=cen['support'];old=nx.Graph();old.add_nodes_from(W)
        G=s.Matrix(cen['gram']);J=s.Matrix(cen['inverse'])
        for i in range(19):
            for j in range(i+1,19):
                if G[i,j]==0:continue
                labels=[m for m in(3,4,5,6)if s.simplify(G[i,j]+magnitude(m))==0]
                assert len(labels)==1
                old.add_edge(W[i],W[j],label=labels[0])
        comps=sorted((sorted(c)for c in nx.connected_components(old)),key=lambda c:(-len(c),c))
        subsets={}
        for bits in itertools.product((False,True),repeat=len(comps)):
            v=sum((c for c,b in zip(comps,bits)if b),[])
            subsets.setdefault(len(v),[]).append(v)
        records={};raw=0
        for typ,g in models:
            if len(g)-1 not in subsets:continue
            for deleted in g:
                rest=g.copy();rest.remove_node(deleted)
                for nodes in subsets[len(rest)]:
                    gm=nx.algorithms.isomorphism.GraphMatcher(rest,old.subgraph(nodes),edge_match=lambda a,b:a['label']==b['label'])
                    for mapping in gm.isomorphisms_iter():
                        raw+=1
                        contacts=tuple(sorted(((mapping[j],g[deleted][j]['label'])for j in g.neighbors(deleted)),key=lambda x:x[0]))
                        key=json.dumps(contacts)
                        rec=records.setdefault(key,dict(contacts=[list(z)for z in contacts],models=[]))
                        rec['models'].append(dict(type=typ,deleted=deleted,mapping={str(k):v for k,v in mapping.items()}))
        rows=[];T=W[:17];JT=G[:17,:17].inv()
        for key,rec in sorted(records.items()):
            d=s.zeros(19,1)
            for j,m in rec['contacts']:d[W.index(j)]=magnitude(m)
            cost=s.simplify((d.T*J*d)[0]);assert cost==2
            coeff=(J*d).applyfunc(s.simplify);assert all(x.is_nonnegative for x in coeff)
            A=G.row_join(-d).col_join((-d.T).row_join(s.Matrix([[2]])))
            delta=coeff.col_join(s.Matrix([1]));assert A*delta==s.zeros(20,1)
            assert A.rank()==19
            dt=d[:17,0];res=s.simplify(2-(dt.T*JT*dt)[0]);assert res.is_nonnegative
            rec.update(id=len(rows),cost='2',projection_coefficients=list(map(str,-coeff)),
                affine_relation=list(map(str,delta)),T_residual_norm=str(res),
                can_meet_target_H3=bool(res.is_positive))
            rows.append(rec)
        out.append(dict(edge_index=ci,edge=cen['edge'],support=W,components=comps,
            raw_matches=raw,profiles=rows,meeting_profile_ids=[r['id']for r in rows if r['can_meet_target_H3']]))
        print('EDGE',ci,'parabolic_profiles',len(rows),'positive_T_norm',sum(r['can_meet_target_H3']for r in rows),flush=True)
    paths=[f'scripts/{TAG}.py',END,'results/fresh_h20_e8_cusp_completion_20261017.json',
           'results/fresh_h20_e8_exact_cluster_20261017.json',
           'results/fresh_h20_e8_f4_cusps_20261017.json',f'scripts/{P.TAG}.py']
    return dict(status='PASS_EXACT_PARABOLIC_ENDPOINT_PROFILE_CENSUS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        affine_model_names=[name for name,g in models],edges=out,
        profile_counts=[len(x['profiles'])for x in out],
        meeting_profile_counts=[len(x['meeting_profile_ids'])for x in out],
        scope='Complete one-new-root PSD corank1 endpoint profiles on the three spherical rank19 supports, using normalized affine Coxeter classification and whole-component deletion matching. Exact cost2, positive component relations and rank19 are checked. This does not choose a future null ray/height, certify old168 compatibility, or complete a cusp/target face.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
