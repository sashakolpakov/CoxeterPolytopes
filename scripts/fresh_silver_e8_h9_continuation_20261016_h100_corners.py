#!/usr/bin/env python3
"""Exact finite corners and literal sector-map proposals for the silver E8 host."""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_silver_e8_h9_20261015 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_silver_e8_h9_continuation_20261016_h100_corners'
SOURCE='results/fresh_silver_e8_h9_continuation_20261016.json'
S,H=A.S,A.H
Z,O=S(0),S(1)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dec(x):return S.from_json(x)
def enc(x):
    if isinstance(x,S):return x.to_json()
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    return x
def transpose(m):return tuple(zip(*m))
def mm(a,b):return tuple(tuple(H.dot(r,c)for c in zip(*b))for r in a)
def pd(m):
    m=[list(r)for r in m];p=[]
    for k in range(len(m)):
        q=m[k][k];p.append(q)
        if q.sign()<=0:return False,p
        for i in range(k+1,len(m)):
            for j in range(k+1,len(m)):m[i][j]-=m[i][k]*m[k][j]/q
    return True,p
def graph(ids,gram):
    g=nx.Graph()
    for i in ids:g.add_node(i,norm=gram[i][i])
    for i,j in combinations(ids,2):
        if gram[i][j]:g.add_edge(i,j,pair=gram[i][j])
    return g

def run():
    d=json.loads((ROOT/SOURCE).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    G=tuple(tuple(map(dec,r))for r in d['gram'])
    N=tuple(tuple(map(dec,r))for r in d['all_walls'])
    C=tuple(map(dec,d['strict_interior']));base=(O,)+(Z,)*9
    cov=[H.mv(G,n)for n in N]
    gram=[[H.dot(n,c)for c in cov]for n in N]
    assert G==A.GRAM and len(N)>=20
    assert H.pair(C,G,C).sign()<0 and all(H.dot(C,c).sign()<0 for c in cov)
    adjacent=[set()for _ in N]
    for i,j in combinations(range(len(N)),2):
        assert gram[i][j].sign()<=0
        if (gram[i][i]*gram[j][j]-gram[i][j]*gram[i][j]).sign()>0:
            adjacent[i].add(j);adjacent[j].add(i)
    supports=[];nodes=[0]*10
    def visit(chosen,remaining):
        k=len(chosen);nodes[k]+=1
        if k==9:supports.append(chosen);return
        for u,j in enumerate(remaining):
            tail=[x for x in remaining[u+1:]if x in adjacent[j]]
            if len(tail)>=8-k:visit(chosen+[j],tail)
    visit([],list(range(len(N))))
    tests=[];vertices=[]
    for ids in supports:
        yes,piv=pd([[gram[i][j]for j in ids]for i in ids])
        tests.append(dict(support=ids,positive=yes,pivots=piv))
        if not yes:continue
        spatial=H.mv(H.inverse(tuple(tuple(cov[i][j]for j in range(1,10))for i in ids)),tuple(-cov[i][0]for i in ids))
        V=(O,*spatial);vals=[H.dot(V,c)for c in cov];norm=H.pair(V,G,V)
        assert norm.sign()<0 and [i for i,x in enumerate(vals)if not x]==ids
        assert all(x.sign()<=0 for x in vals)
        vertices.append(dict(incident=ids,point=V,norm=norm,wall_pairings=vals))
    edges=[dict(vertices=[i,j],support=sorted(set(a['incident'])&set(b['incident'])))
           for (i,a),(j,b)in combinations(enumerate(vertices),2)
           if len(set(a['incident'])&set(b['incident']))==8]
    src=list(range(9));SG=[[gram[i][j]for j in src]for i in src]
    sf=transpose([N[i]for i in src]+[base]);si=H.inverse(sf)
    maps=[]
    for vi,v in enumerate(vertices):
        target=v['incident']
        matcher=nx.algorithms.isomorphism.GraphMatcher(graph(src,gram),graph(target,gram),
            node_match=lambda x,y:x['norm']==y['norm'],edge_match=lambda x,y:x['pair']==y['pair'])
        for mapping in matcher.isomorphisms_iter():
            dst=[mapping[i]for i in src];V=v['point']
            df=transpose([N[i]for i in dst]+[V])
            scale=A.silver_determinant(sf)/A.silver_determinant(df)
            if scale.sign()<0:scale=-scale
            assert H.pair(base,G,base)==scale*scale*v['norm']
            T=mm(transpose([N[i]for i in dst]+[tuple(scale*x for x in V)]),si)
            TI=H.inverse(T)
            assert mm(mm(transpose(T),G),T)==G
            assert all(H.mv(T,N[i])==N[j]for i,j in zip(src,dst))
            integral=all(x.is_integral for r in (*T,*TI)for x in r)
            trace=sum((T[i][i]for i in range(10)),Z)
            maps.append(dict(target_vertex=vi,source_indices=src,target_indices=dst,
                             scale=scale,matrix=T,inverse=TI,integral=integral,
                             determinant=A.silver_determinant(T),trace=trace,
                             trace_exceeds_dimension=(trace-10).sign()>0))
    paths=[f'scripts/{TAG}.py',SOURCE,'scripts/fresh_silver_e8_h9_20261015.py',
           'scripts/fresh_silver_f4_h9_20261013.py','src/coxeter_polytopes/exact/silver.py']
    return enc(dict(status='PASS_COMPLETE_FINITE_CORNERS_AND_EXACT_BASE_SECTOR_MAPS',
        sources={p:sha(p)for p in paths},form=G,basepoint=base,
        counts=dict(walls=len(N),finite_neighbor_supports=len(supports),finite_vertices=len(vertices),
                    finite_edges=len(edges),open_finite_germs=9*len(vertices)-2*len(edges),
                    base_sector_maps=len(maps),integral_sector_maps=sum(m['integral']for m in maps)),
        clique_recursion_nodes=nodes,support_tests=tests,finite_vertices=vertices,finite_edges=edges,
        sector_maps=maps,
        scope='Complete finite-corner census of the completeheight100 prefix and all exact norm-coloured maps from its initial E8+A1 sector to those corners. No complete ideal/hyperideal cone census. The primitive norm alphabet/chronological completeness is an explicit inherited input; host-wide conclusions require a separate persistence/isometry certificate and independent review.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved finite-corner replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
    print('integral traces',[(m['target_vertex'],m['trace'])for m in d['sector_maps']if m['integral']])
