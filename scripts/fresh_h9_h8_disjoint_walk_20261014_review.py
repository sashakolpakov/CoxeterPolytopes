#!/usr/bin/env python3
"""Independent number-field/nullspace review of the four-cap extension."""
import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
import fresh_h9_h8_all_high_20261014_review as C
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_disjoint_walk_20261014_review'
PROOF='results/fresh_h9_h8_disjoint_walk_20261014_exact.json'
FRONT='results/fresh_h9_h8_disjoint_walk_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dm(M):return [[C.decode(x)for x in row]for row in M]
def dv(v):return [C.decode(x)for x in v]
def positive(x):return bool(x)and C.ev(x).l>0
def negative(x):return bool(x)and C.ev(x).h<0
def projective(v):
    a=next(x for x in v if x);return tuple(x/a for x in v)
def pd(M):
    W=[row[:]for row in M]
    for k in range(len(W)):
        p=W[k][k]
        if not positive(p):return False
        for i in range(k+1,len(W)):
            for j in range(k+1,len(W)):W[i][j]-=W[i][k]*W[k][j]/p
    return True
def run():
    data=json.loads((ROOT/PROOF).read_text());front=json.loads((ROOT/FRONT).read_text())
    pending,seen=[PROOF,FRONT],set()
    while pending:
        path=pending.pop()
        if path in seen:continue
        seen.add(path);doc=json.loads((ROOT/path).read_text())
        for p,h in doc.get('sources',{}).items():
            assert sha(p)==h,p
            if p.endswith('.json'):pending.append(p)
    rad=data['field']['new_positive_radical_square']
    R=s.Rational(rad[0])+s.Rational(rad[1])*s.sqrt(5)
    assert s.simplify(R-(6+s.Rational(16,5)*s.sqrt(5)))==0
    C.RR=s.sqrt(R);C.K=s.QQ.algebraic_field(s.sqrt(5),C.RR);C.Z,C.O=C.K.zero,C.K.one
    K,Z,O=C.K,C.Z,C.O
    Q,N=dm(data['form']),dm(data['normals']);H=dv(data['center'])
    assert len(N)==19 and len(Q)==10
    pivots=[C.ev(x)for x in C.ldl(Q)]
    assert sum(x.l>0 for x in pivots)==9 and sum(x.h<0 for x in pivots)==1
    assert negative(C.pair(H,H,Q))and all(negative(C.pair(H,n,Q))for n in N)
    norms=[C.pair(n,n,Q)for n in N];assert all(positive(x)for x in norms)
    assert norms[15:]==[2*O]*4
    G=[[C.pair(n,m,Q)for m in N]for n in N];assert G==dm(data['gram'])
    cosine5=K.from_sympy((3+s.sqrt(5))/8);adj=[set()for _ in N];new_counts=Counter()
    for i,j in combinations(range(19),2):
        g=G[i][j];assert not g or negative(g)
        if not g:typ='right'
        elif 4*g*g==norms[i]*norms[j]:typ='pi3'
        elif g*g==cosine5*norms[i]*norms[j]:typ='pi5'
        else:assert positive(g*g-norms[i]*norms[j]);typ='dotted'
        if typ!='dotted':adj[i].add(j);adj[j].add(i)
        if j>=15:new_counts[typ]+=1
    assert dict(new_counts)==data['new_pair_counts']and sum(new_counts.values())==66
    assert {str(i):sorted(a)for i,a in enumerate(adj)}==data['finite_neighbors']
    # Independently filter all C(19,9) supports, rather than use the
    # producer's recursive clique generator.
    candidates=[]
    for ids in combinations(range(19),9):
        if all(j in adj[i]for i,j in combinations(ids,2)):candidates.append(ids)
    vertices=[];bad=0
    for ids in candidates:
        M=[[G[i][j]for j in ids]for i in ids]
        if not pd(M):bad+=1;continue
        v=C.null_one([C.mv(Q,N[i])for i in ids])
        if positive(C.pair(v,H,Q)):v=[-x for x in v]
        assert negative(C.pair(v,H,Q))and negative(C.pair(v,v,Q))
        vals=[C.pair(v,n,Q)for n in N]
        assert [i for i,x in enumerate(vals)if not x]==list(ids)
        assert all(not x or negative(x)for x in vals)
        vertices.append((ids,v))
    saved=data['complete_finite_vertex_census']
    assert len(candidates)==saved['candidate_cliques']==799 and bad==saved['nonpositive']==746
    assert len(vertices)==saved['count']==53
    for (ids,v),row in zip(vertices,saved['vertices']):
        assert list(ids)==row['incident']and projective(v)==projective(dv(row['point']))
    assert len({projective(v)for _,v in vertices})==53
    edges=[]
    for i,j in combinations(range(53),2):
        support=sorted(set(vertices[i][0])&set(vertices[j][0]))
        if len(support)==8:
            assert pd([[G[k][l]for l in support]for k in support])
            edges.append(dict(vertices=[i,j],support=support))
    assert edges==data['complete_finite_edges']and len(edges)==175
    for row in data['closed_old_germs']:
        germ=front['queue'][row['germ']];ids,v=vertices[row['new_vertex']]
        start,direction=dv(germ['start_point']),dv(germ['direction'])
        alpha=C.pair(v,start,Q)/C.pair(start,start,Q)
        beta=C.pair(v,direction,Q)/C.pair(direction,direction,Q)
        assert positive(alpha)and positive(beta)
        assert all(x==alpha*y+beta*z for x,y,z in zip(v,start,direction))
        assert beta/alpha==C.decode(row['forward_parameter'])
        assert set(row['support'])<=set(ids)
    assert len(data['closed_old_germs'])==6
    for i,row in enumerate(data['facet_points']):
        v=dv(row['point']);vals=[C.pair(v,n,Q)for n in N]
        assert not vals[i]and all(negative(x)for j,x in enumerate(vals)if j!=i)
        assert negative(C.pair(v,v,Q))and negative(C.pair(v,H,Q))
    assert len(data['facet_points'])==19
    E=dv(data['strict_escape']['point'])
    assert positive(C.pair(E,E,Q))and negative(C.pair(E,H,Q))
    assert all(negative(C.pair(E,n,Q))for n in N)
    paths=[PROOF,FRONT,f'scripts/{TAG}.py','scripts/fresh_h9_h8_all_high_20261014_review.py',
           'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_FOUR_CAP_COMPLETE_FINITE_CENSUS_REVIEW',
        sources={p:sha(p)for p in paths},new_pairs=66,candidate_supports=799,
        nonpositive_supports=746,finite_vertices=53,finite_edges=175,closed_old_germs=6,
        actual_facets=19,scope='No quartet producer or Fraction extension-class imports. Separate SymPy number field and physical nullspaces reconstruct all171wall pairs, all92378nine-root subsets (799finite cliques),53complete finite vertices,175finite edges,six forward old-edge closures,19actualfacets and the surviving strict escape. Ideal/exterior vertices and the renewed full frontier are not classified here. Volume remains infinite.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['finite_vertices'],result['finite_edges'])
