#!/usr/bin/env python3
"""Bounded exact finite-corner map search on the norm131 prefix."""
import argparse,json
from hashlib import sha256
from pathlib import Path
from itertools import combinations
import networkx as nx
import fresh_h9_beta11_vinberg_h300_20261014_review as A
import fresh_h9_beta55_corner_symmetries_20261014 as M
K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta131_h100_20261016_maps'
PREFIX='results/fresh_h9_beta131_h100_20261016.json'
GEOMETRY='results/fresh_h9_beta131_h100_20261016_geometry.json'
LIMIT=10000
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def graph(ids,G):
    g=nx.Graph()
    for i in ids:
        g.add_node(i,q=G[i][i])
    for i,j in combinations(ids,2):
        if G[i][j]:g.add_edge(i,j,g=G[i][j])
    return g
def common_fixed_timelike(maps,G,H):
    rows=[]
    for m in maps:
        rows.extend([[x-K(i==j)for j,x in enumerate(r)]for i,r in enumerate(m['matrix'])])
    pivots=[]
    for j in range(10):
        k=len(pivots);t=next((i for i in range(k,len(rows))if rows[i][j]),None)
        if t is None:continue
        rows[t],rows[k]=rows[k],rows[t];q=rows[k][j];rows[k]=[x/q for x in rows[k]]
        for i in range(len(rows)):
            if i!=k:
                q=rows[i][j];rows[i]=[x-q*y for x,y in zip(rows[i],rows[k])]
        pivots.append(j)
    orthogonal=[]
    for j in range(10):
        if j in pivots:continue
        v=[K(i==j)for i in range(10)]
        for r,i in zip(rows,pivots):v[i]=-r[j]
        for w,q in orthogonal:
            coefficient=A.pair(v,w,G)/q;v=[x-coefficient*y for x,y in zip(v,w)]
        q=A.pair(v,v,G);assert q
        if q.sign()<0:
            if A.pair(v,H,G).sign()>0:v=[-x for x in v]
            assert all(A.mv(m['matrix'],v)==v for m in maps)
            return dict(point=v,norm=q,future_basepoint_pairing=A.pair(v,H,G),fixed_equation_rank=len(pivots))
        orthogonal.append((v,q))
    return None
def run():
    d,g=map(read,(PREFIX,GEOMETRY));G=[[A.dec(x)for x in r]for r in d['form']]
    N=[list(map(A.dec,r))for r in d['roots']];gram=[[A.pair(a,b,G)for b in N]for a in N]
    vertices=g['finite_vertices'];graphs=[graph(v['incident'],gram)for v in vertices]
    frames=[]
    for v in vertices:
        F=M.tr([N[i]for i in v['incident']]+[list(map(A.dec,v['point']))])
        FI,det=M.inverse_det(F);frames.append((F,FI,det))
    tested=0;integral=[];hit=None;complete=True
    for i in range(len(vertices)):
        for j in range(i,len(vertices)):
            gm=nx.algorithms.isomorphism.GraphMatcher(graphs[i],graphs[j],
                node_match=lambda a,b:a['q']==b['q'],edge_match=lambda a,b:a['g']==b['g'])
            for mapping in gm.isomorphisms_iter():
                if tested>=LIMIT:complete=False;break
                tested+=1;src=vertices[i]['incident'];dst=[mapping[v]for v in src]
                V=list(map(A.dec,vertices[i]['point']));P=list(map(A.dec,vertices[j]['point']))
                F,FI,detF=frames[i];D=M.tr([N[k]for k in dst]+[P]);_,detD=M.inverse_det(D)
                scale=detF/detD
                if scale.sign()<0:scale=-scale
                assert A.pair(V,V,G)==scale*scale*A.pair(P,P,G)
                for row in D:row[-1]*=scale
                T=M.mm(D,FI)
                if not all(x.integral()for r in T for x in r):continue
                TI,det=M.inverse_det(T);assert all(x.integral()for r in TI for x in r)
                assert M.mm(M.tr(T),M.mm(G,T))==G
                assert all(A.mv(T,N[a])==N[b]for a,b in zip(src,dst))
                assert A.mv(T,V)==[scale*x for x in P]
                trace=sum((T[k][k]for k in range(10)),K())
                row=dict(source_vertex=i,target_vertex=j,source_ids=src,target_ids=dst,
                    matrix=T,inverse=TI,scale=scale,determinant=det,trace=trace)
                integral.append(row)
                if (trace-10).sign()>0:hit=dict(kind='single',map=len(integral)-1,trace=trace);break
            if hit or not complete:break
        if hit or not complete:break
    if not hit:
        for i,a in enumerate(integral):
            for j,b in enumerate(integral[:i]):
                T=M.mm(a['matrix'],b['matrix']);trace=sum((T[k][k]for k in range(10)),K())
                if (trace-10).sign()>0:
                    hit=dict(kind='product',maps=[i,j],trace=trace,matrix=T);break
            if hit:break
    fixed=common_fixed_timelike(integral,G,list(map(A.dec,d['basepoint'])))if not hit else None
    paths=[f'scripts/{TAG}.py',PREFIX,GEOMETRY,
        'scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py',
        'scripts/fresh_h9_beta55_corner_symmetries_20261014.py']
    return A.enc(dict(status='PASS_EXACT_BOUNDED_BETA131_CORNER_MAP_SCREEN',
        sources={p:sha(p)for p in paths},finite_vertices=len(vertices),
        checked_maps=tested,limit=LIMIT,all_corner_pairs_exhausted=complete and not hit,
        integral_maps=integral,hit=hit,common_fixed_timelike=fixed,
        scope='Only exact norm-colored outward-sector maps between known finite vertices. A common timelike fixed vector proves the group generated by all saved integral maps finite, because the ambient integral isometry group is discrete. This does not exclude new maps at later corners or arbitrary additional symmetries. No global reflectivity claim.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved corner-map screen replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'maps',d['checked_maps'],'integral',len(d['integral_maps']),'hit',d['hit'])
