#!/usr/bin/env python3
"""Complete finite-vertex census of the 28-wall height100 prefix."""
import argparse,json
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import fresh_frontier_construction_20261015 as A
import fresh_frontier_construction_20261015_geometry as B
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_norm31_continuation_20261015_corners'
SOURCE='results/fresh_norm31_continuation_20261015.json'
HELPER='scripts/fresh_frontier_construction_20261015.py'
NULL='scripts/fresh_frontier_construction_20261015_geometry.py'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    d=json.loads((ROOT/SOURCE).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    G=[[A.dec(x)for x in r]for r in d['form']]
    N=[list(map(A.dec,r))for r in d['roots']]
    GG=[[A.dec(x)for x in r]for r in d['gram']]
    C=list(map(A.dec,d['strict_center']));H=list(map(A.dec,d['basepoint']))
    NG=[A.mv(G,n)for n in N];n=len(N);assert n==28
    adjacent=[set()for _ in N]
    for i in range(n):
        for j in range(i):
            if (GG[i][i]*GG[j][j]-GG[i][j]*GG[i][j]).sign()>0:
                adjacent[i].add(j);adjacent[j].add(i)
    supports=[];nodes=[0]*10
    def visit(chosen,candidates):
        k=len(chosen);nodes[k]+=1
        if k==9:supports.append(chosen);return
        if len(candidates)<9-k:return
        for u,j in enumerate(candidates):visit(chosen+[j],[x for x in candidates[u+1:]if x in adjacent[j]])
    visit([],list(range(n)))
    vertices=[];tests=[]
    for ids in supports:
        M=[[GG[i][j]for j in ids]for i in ids];piv=[];ok=True
        for k in range(9):
            q=M[k][k];piv.append(q)
            if q.sign()<=0:ok=False;break
            for i in range(k+1,9):
                for j in range(k+1,9):M[i][j]-=M[i][k]*M[k][j]/q
        tests.append(dict(support=ids,pivots=piv,positive=ok))
        if not ok:continue
        V=B.null([NG[i]for i in ids]);assert V
        if A.pair(V,C,G).sign()>0:V=[-x for x in V]
        norm=A.pair(V,V,G);assert norm.sign()<0
        vals=[A.pair(V,r,G)for r in N]
        assert [i for i,x in enumerate(vals)if not x]==ids and all(x.sign()<=0 for x in vals)
        vertices.append(dict(incident=ids,point=V,norm=norm,basepoint_pairing=A.pair(V,H,G),wall_pairings=vals))
    edges=[]
    for i,j in combinations(range(len(vertices)),2):
        common=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
        if len(common)==8:edges.append(dict(vertices=[i,j],support=common))
    assert len(supports)==35 and len(vertices)==11
    return A.enc(dict(status='PASS_COMPLETE_NORM31_HEIGHT100_FINITE_VERTEX_CENSUS',
        sources={f:sha(f)for f in [f'scripts/{TAG}.py',HELPER,NULL,SOURCE]},
        counts=dict(finite_neighbor_supports=len(supports),finite_vertices=len(vertices),
            finite_edges=len(edges),open_finite_germs=9*len(vertices)-2*len(edges)),
        clique_recursion_nodes=nodes,support_tests=tests,finite_vertices=vertices,finite_edges=edges,
        scope='Complete finite-vertex census only: every finite corner has nine pairwise finite-related walls. No complete ideal or hyperideal cone census is claimed.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS finite-corner saved replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'])
