#!/usr/bin/env python3
"""Complete finite-corner census via nondotted cliques and exact tree Schur tests."""
import argparse,json
from collections import Counter
from itertools import combinations
from hashlib import sha256
from pathlib import Path
import networkx as nx
import fresh_h9_beta11_vinberg_h300_20261014_review as A
import fresh_h9_beta55_corner_symmetries_20261014 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta131_h300_20261016_finite'
PREFIX='results/fresh_h9_beta131_h300_20261016.json'
OLD='results/fresh_h9_beta131_h100_20261016_geometry.json'
K=A.K
def sha(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def read(path):
    d=json.loads((ROOT/path).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    return d
def vector(row):return [A.dec(x) for x in row]

def run():
    d,old=read(PREFIX),read(OLD)
    G=[vector(row) for row in d['form']];N=[vector(row) for row in d['roots']]
    H=vector(d['strict_center']);gram=[vector(row) for row in d['gram']]
    graph=nx.Graph();graph.add_nodes_from(range(len(N)))
    labels=[[0]*len(N) for _ in N]
    for i in range(len(N)):
        for j in range(i):
            value=4*gram[i][j]*gram[i][j]/(gram[i][i]*gram[j][j])
            if (value-4).sign()>=0:continue
            graph.add_edge(i,j)
            if value:
                assert gram[i][i]==gram[j][j]==2
                assert gram[i][j] in (K(-1),K(0,-1))
                labels[i][j]=labels[j][i]=3 if gram[i][j]==-1 else 5
    maxima=[tuple(sorted(c)) for c in nx.find_cliques(graph)]
    candidates=set()
    for c in maxima:
        if len(c)>=9:candidates.update(combinations(c,9))
    cache={};statistics=Counter();positive=[]
    def schur(code):
        if code in cache:return cache[code][0]
        pivot=K(2);good=True
        for label,child in code:
            p=schur(child)
            if p is None:good=False;break
            pivot-=(K(1) if label==3 else K(1,1))/p
        result=pivot if good and pivot.sign()>0 else None
        cache[code]=(result,pivot if good else None)
        return result
    def definite(support):
        adjacency={i:[] for i in support}
        for i,j in combinations(support,2):
            if labels[i][j]:adjacency[i].append(j);adjacency[j].append(i)
        unseen=set(support)
        while unseen:
            root=min(unseen);stack=[root];component=set()
            while stack:
                i=stack.pop()
                if i in component:continue
                component.add(i);stack.extend(adjacency[i])
            unseen-=component
            if sum(len(adjacency[i]) for i in component)!=2*(len(component)-1):
                statistics['cycle_rejections']+=1;return False
            def code(i,parent):
                return tuple(sorted((labels[i][j],code(j,i)) for j in adjacency[i] if j!=parent))
            if schur(code(root,None)) is None:
                statistics['tree_schur_rejections']+=1;return False
        return True
    for k,support in enumerate(sorted(candidates)):
        if definite(support):positive.append(support)
        if (k+1)%50000==0:print('finite census',k+1,'of',len(candidates),'positive',len(positive),flush=True)
    covectors=[A.mv(G,n) for n in N]
    vertices=[]
    for support in positive:
        B=[[gram[i][j] for j in support] for i in support]
        inverse,det=M.inverse_det(B)
        assert all(x.sign()>=0 for row in inverse for x in row)
        coefficients=A.mv(inverse,[A.pair(H,N[i],G) for i in support])
        v=[H[k]-sum((c*N[i][k] for c,i in zip(coefficients,support)),K()) for k in range(10)]
        q=A.pair(v,v,G)
        assert q.sign()<0 and A.pair(v,H,G).sign()<0
        values=A.mv(covectors,v)
        assert all(x.sign()<=0 for x in values)
        assert tuple(i for i,x in enumerate(values) if not x)==support
        vertices.append(dict(incident=list(support),point=v,norm=q))
    assert {tuple(v['incident']) for v in old['finite_vertices']}<=set(positive)
    E=vector(old['strict_escape']['point']);ev=A.mv(covectors,E)
    assert A.pair(E,E,G).sign()>0 and A.pair(E,H,G).sign()<0 and all(x.sign()<0 for x in ev)
    paths=[PREFIX,OLD,f'scripts/{TAG}.py','scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py',
           'scripts/fresh_h9_beta55_corner_symmetries_20261014.py']
    return A.enc(dict(status='PASS_COMPLETE_BETA131_H300_FINITE_CORNER_CENSUS',sources={p:sha(p) for p in paths},
        counts=dict(walls=len(N),maximal_nondotted_cliques=len(maxima),
                    unique_nine_subsets=len(candidates),finite_vertices=len(vertices),
                    cached_rooted_trees=len(cache),**statistics),
        maximal_clique_size_counts={str(k):v for k,v in sorted(Counter(map(len,maxima)).items())},
        tree_schur_checks=[dict(code=c,positive=x is not None,pivot=p) for c,(x,p) in sorted(cache.items())],
        finite_vertices=vertices,strict_escape=dict(point=E,norm=A.pair(E,E,G),pairings=ev),
        completeness='A finite vertex of an acute Coxeter chamber is simple and supplies a positive-definite nine-clique. Every nine-clique occurs in a maximal nondotted clique. A nonright cycle is non-PD because all normalized edges have magnitude at least one. All remaining trees are tested exactly by leaf Schur pivots. The only non-norm-two root is orthogonal to every finite neighbor, so its positive isolated factor is harmless. Every positive support is reconstructed and tested against every actual wall.',
        scope='Complete finite-vertex census only; ideal and hyperideal rays are not enumerated. The same strict spacelike escape proves this finite prefix infinite-volume.'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text());print('PASS saved complete finite-corner replay')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
