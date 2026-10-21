#!/usr/bin/env python3
"""Separate algebraic-field census of the changed H8-based H9 seed."""
import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
from fresh_h9_h8_all_high_20261014_review import (
    K,Z,O,SQ5,decode,source_decode,ev,pair,mv,dot,eliminate,ldl,encode)
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_disjoint_cap_20261014_review'
PROOF='results/fresh_h9_h8_disjoint_cap_20261014.json'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def key(v):
    q=next(x for x in v if x);return tuple(x/q for x in v)
def dm(M):return [[decode(x)for x in row]for row in M]
def dv(v):return [decode(x)for x in v]
def neg(x):return bool(x)and ev(x).h<0
def pos(x):return bool(x)and ev(x).l>0
def pd(G):return all(pos(x)for x in ldl(G))
def kernel(M):
    R,piv=eliminate(M);n=len(M[0])
    if len(piv)!=n-1:return None
    free=next(i for i in range(n)if i not in piv);v=[Z]*n;v[free]=O
    for i,p in enumerate(piv):v[p]=-R[i][free]
    assert all(not x for x in mv(M,v));return v
def run():
    p=json.loads((ROOT/PROOF).read_text())
    for path,h in p['sources'].items():assert sha(path)==h,path
    src=json.loads((ROOT/SOURCE).read_text())
    Q,N=dm(p['form']),dm(p['normals']);H=dv(p['center'])
    norms=[pair(n,n,Q)for n in N]
    assert len(N)==15 and len(Q)==10
    signs=[ev(x)for x in ldl(Q)]
    assert sum(x.l>0 for x in signs)==9 and sum(x.h<0 for x in signs)==1
    assert neg(pair(H,H,Q))and all(neg(pair(H,n,Q))for n in N)
    Q8=[[source_decode(x)for x in row]for row in src['lattice_gram']]
    U=[[source_decode(x)for x in row]for row in src['simple_roots']]
    assert Q[:9]==[row+[Z]for row in Q8]
    assert N[:11]==[u+[Z]for u in U]
    d=K.from_sympy(SQ5/4);a=K.from_sympy((SQ5-1)/4)
    assert Q[9]==[Z]*9+[d]and N[11]==[Z]*9+[O]
    assert norms==[2*O]*11+[d]+[O/4]*3
    G=[[pair(n,m,Q)for m in N]for n in N]
    assert G==dm(p['gram'])
    cosine5sq=K.from_sympy((3+SQ5)/8)
    types={};counts=Counter()
    for i,j in combinations(range(15),2):
        g=G[i][j];assert not g or neg(g)
        if not g:typ='right'
        elif 4*g*g==norms[i]*norms[j]:typ='pi3'
        elif g*g==cosine5sq*norms[i]*norms[j]:typ='pi5'
        else:assert pos(g*g-norms[i]*norms[j]);typ='dotted'
        types[i,j]=typ;counts[typ]+=1
    assert dict(counts)==p['pair_counts']==dict(right=79,pi3=8,pi5=2,dotted=16)
    # Complete source census, independently solved over the SymPy number field.
    source_rows=[mv(Q8,u)for u in U]
    dual=[sum((u[i]for u in U),Z)for i in range(9)]
    assert neg(pair(dual,dual,Q8))
    source_vertices={};rank_deficient=0
    for ids in combinations(range(11),8):
        v=kernel([source_rows[i]for i in ids])
        if v is None:rank_deficient+=1;continue
        orient=pair(v,dual,Q8)
        if not orient:continue
        if pos(orient):v=[-x for x in v]
        vals=[dot(v,row)for row in source_rows]
        if any(pos(x)for x in vals):continue
        assert all(not x or neg(x)for x in vals)
        assert neg(pair(v,v,Q8))
        active=tuple(i for i,x in enumerate(vals)if not x);assert len(active)==8
        assert pd([[G[i][j]for j in active]for i in active])
        assert active not in source_vertices;source_vertices[active]=v
    assert len(source_vertices)==41
    assert set(source_vertices)=={tuple(x)for x in src['geometry']['vertices']}
    assert rank_deficient==p['source_support_census']['rank_deficient']
    vertices=[]
    for ids in src['geometry']['vertices']:
        vertices.append((source_vertices[tuple(ids)]+[Z],sorted(ids+[11])))
    for j,vi in zip((12,13,14),(0,32,33)):
        ids=src['geometry']['vertices'][vi]
        near=[i for i in range(15)if i!=j and types[tuple(sorted((i,j)))]!='dotted']
        assert near==ids
        cap=N[j];v=cap[:9]
        assert cap[9]==-O and pair(v,v,Q8)==-a
        assert key(v)==key(source_vertices[tuple(ids)])
        vertices.append((v+[-a/d],sorted(ids+[j])))
    assert len(vertices)==44
    for (v,ids),saved in zip(vertices,p['complete_finite_vertices']):
        assert key(v)==key(dv(saved['point']))and ids==sorted(saved['incident'])
        vals=[pair(v,n,Q)for n in N]
        assert all(not x or neg(x)for x in vals)
        assert [i for i,x in enumerate(vals)if not x]==ids
        assert neg(pair(v,v,Q))and neg(pair(v,H,Q))
        assert pd([[G[i][j]for j in ids]for i in ids])
    # Source-only positive support cannot have rank9 in signature(8,1).
    source_signs=[ev(x)for x in ldl(Q8)]
    assert sum(x.l>0 for x in source_signs)==8 and sum(x.h<0 for x in source_signs)==1
    edges=[]
    for i,j in combinations(range(44),2):
        shared=sorted(set(vertices[i][1])&set(vertices[j][1]))
        if len(shared)==8:
            assert pd([[G[k][l]for l in shared]for k in shared])
            edges.append(dict(vertices=[i,j],support=shared))
    assert edges==p['complete_finite_edges']and len(edges)==167
    assert [row['vertices']for row in edges if max(row['vertices'])>=41]==[[0,41],[32,42],[33,43]]
    assert len(p['facet_points'])==15
    for i,row in enumerate(p['facet_points']):
        v=dv(row['point']);vals=[pair(v,n,Q)for n in N]
        assert not vals[i]and all(neg(x)for j,x in enumerate(vals)if j!=i)
        assert neg(pair(v,v,Q))and neg(pair(v,H,Q))
    escape=dv(p['strict_escape']['point']);escape_norm=pair(escape,escape,Q)
    assert pos(escape_norm)and neg(pair(escape,H,Q))
    assert all(neg(pair(escape,n,Q))for n in N)
    cycle=4*G[11][12]*G[11][12]/(norms[11]*norms[12])
    assert cycle==K.from_sympy(4*SQ5)
    assert s.simplify(K.to_sympy(cycle).xreplace({SQ5:-SQ5}))==-4*SQ5
    phi=K.from_sympy((1+SQ5)/2)
    ratio_norms=[]
    for old in (2*phi,phi-O):
        expr=K.to_sympy(d/old);conj=expr.xreplace({SQ5:-SQ5})
        nn=s.simplify(expr*conj);assert nn.is_Rational and not s.sqrt(nn).is_Rational
        ratio_norms.append(str(nn))
    assert ratio_norms==['5/64','5/16']
    paths=[PROOF,SOURCE,f'scripts/{TAG}.py','scripts/fresh_h9_h8_all_high_20261014_review.py',
           'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_CHANGED_H8_CORE_CENSUS_AND_ESCAPE_REVIEW',
                sources={path:sha(path)for path in paths},pair_counts=dict(counts),
                source_supports=165,finite_vertices=44,finite_edges=167,actual_facets=15,
                remaining_finite_corner_germs=44*9-2*167,exact_escape_norm=encode(escape_norm),
                conjugate_two_cycle='-4*sqrt(5)',determinant_ratio_norms=ratio_norms,
                scope='Separate SymPy number-field arithmetic reconstructs the105pairs, all165source supports,44complete finite vertices,167finite edges,15facet witnesses and strict escape. Completeness uses each new cap having exactly8finite neighbors, bottom-face completeness and source signature(8,1). Does not enumerate ideal/spacelike vertices. The strict escape proves infinite volume; the negative conjugate two-cycle obstructs quasi-arithmetic finite-covolume completions retaining the displayed pair, not arbitrary nonarithmetic completions.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['finite_vertices'],result['finite_edges'])
