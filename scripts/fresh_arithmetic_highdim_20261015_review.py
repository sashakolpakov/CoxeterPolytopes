#!/usr/bin/env python3
"""Independent graph, algebra and compact-face audit of the arithmetic bound.

No producer imports. Full face argument and arithmeticity criterion are source
proof obligations; this replay independently checks the finite combinatorics,
all-rank graph formulas and exact scalar obstruction they use.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb,gcd
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_arithmetic_highdim_20261015_review'
PRIMARY='fresh_arithmetic_highdim_20261015'


def path(n):
    return [(i,i+1)for i in range(n-1)]


def arms(lengths):
    edges=[];k=1
    for length in lengths:
        last=0
        for _ in range(length):
            edges.append((last,k));last=k;k+=1
    return k,edges


def pair_count(n,edges):
    adj=[set()for _ in range(n)]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    pairs=[]
    for a in range(n):
        dist={a:0};queue=deque([a])
        while queue:
            u=queue.popleft()
            for v in adj[u]:
                if v not in dist:dist[v]=dist[u]+1;queue.append(v)
        pairs.extend((a,b)for b in range(a+1,n)if dist.get(b,n+1)<=3)
    deg=list(map(len,adj))
    # A tree has exactly one path between endpoints. This independent count
    # partitions unordered pairs by their path length1,2 or3.
    by_length=[len(edges),sum(comb(d,2)for d in deg),
               sum((deg[a]-1)*(deg[b]-1)for a,b in edges)]
    assert sum(by_length)==len(pairs)
    return len(pairs),by_length


def connected(nodes,edges):
    reach={nodes[0]}
    while True:
        larger=reach|{b for a,b in edges if a in reach}|{a for a,b in edges if b in reach}
        if larger==reach:return reach==set(nodes)
        reach=larger


def connected_graphs(nodes):
    possible=list(combinations(nodes,2))
    for mask in range(1<<len(possible)):
        edges=[edge for j,edge in enumerate(possible)if mask&(1<<j)]
        if connected(nodes,edges):yield edges


def shortest_white_path(edges):
    # Quad black sets {0,2} and{1,3}; a path may use only white internal nodes.
    best=100
    for a in(0,2):
        for b in(1,3):
            nodes={a,b,4,5};adj={v:set()for v in nodes}
            for x,y in edges:
                if x in nodes and y in nodes:adj[x].add(y);adj[y].add(x)
            queue=deque([(a,0)]);seen={a}
            while queue:
                v,d=queue.popleft()
                if v==b:best=min(best,d);break
                for u in adj[v]-seen:seen.add(u);queue.append((u,d+1))
    return best


def run():
    source=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert source['status'].startswith('PASS')
    for p,h in source.get('sources',{}).items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    graph_rows=[]
    maxima={1:0,2:1}
    for n in range(1,81):
        w,counts=pair_count(n,path(n))
        expected=0 if n==1 else 1 if n==2 else 3*n-6
        assert w==expected
        graph_rows.append(dict(type='path',rank=n,weight=w,counts=counts))
        maxima[n]=max(maxima.get(n,0),w)
        if n>=4:
            N,edges=arms([1,1,n-3]);assert N==n
            w,counts=pair_count(n,edges)
            assert w==(6 if n==4 else 3*n-5)
            graph_rows.append(dict(type='D',rank=n,weight=w,counts=counts))
            maxima[n]=max(maxima[n],w)
    for n in(6,7,8):
        N,edges=arms([1,2,n-4]);assert N==n
        w,counts=pair_count(n,edges);assert w==3*n-4
        graph_rows.append(dict(type='E',rank=n,weight=w,counts=counts))
        maxima[n]=max(maxima[n],w)
    # Bind the producer's literal finite graph census to this separately
    # implemented BFS and degree-path formula.
    for row in source['graph_checks']:
        weight,counts=pair_count(row['rank'],row['edges'])
        assert weight==row['weight'] and counts==row['distance_counts']
        assert row['deficit']==3*row['rank']-weight
    assert all(maxima[n]<=3*n-4 for n in range(3,81))
    assert 3-maxima[1]==3 and 6-maxima[2]==5
    # Forest partition dynamic programming is an independent check at the
    # critical dimension, including arbitrary finite I2(m) via its graph shape.
    dp=[0]+[-10000]*29
    for total in range(1,30):
        dp[total]=max(dp[total-r]+maxima[r]for r in range(1,total+1))
    assert dp[29]==82<=3*29-4==83
    n=s.symbols('n',integer=True,positive=True)
    assert s.expand((n-1)+(n-2)+(n-3))==3*n-6
    assert s.expand((n-1)+(n-1)+(n-3))==3*n-5
    assert s.expand((n-1)+(n-1)+(n-2))==3*n-4

    # Every connected three-black graph has at least two marked pairs joined
    # directly; each is a different actual corner in a triangular2-face.
    triangles=[]
    for edges in connected_graphs([0,1,2]):
        assert len(edges)>=2
        triangles.append(dict(edges=edges,forced_corner_weights=len(edges)))
    assert len(triangles)==4
    # Independent exhaustive tiny topology check for the quadrilateral bridge.
    # Each minimal indefinite opposite-pair graph has at most one white node.
    # If disjoint, at least one cross-edge exists; if overlapping, the only
    # possible intersection is their white node. Extra edges only shorten paths.
    quadrilaterals=[]
    for left in([0,2],[0,2,4]):
        for right in([1,3],[1,3,5],[1,3,4]):
            if 4 in right and 4 not in left:continue
            for A in connected_graphs(left):
                for B in connected_graphs(right):
                    overlap=set(left)&set(right)
                    bridges=[None]if overlap else[(a,b)for a in left for b in right]
                    for bridge in bridges:
                        edges=A+B+([]if bridge is None else[bridge])
                        length=shortest_white_path(edges)
                        assert length<=3
                        quadrilaterals.append(dict(left_nodes=left,right_nodes=right,
                                                   left_edges=A,right_edges=B,bridge=bridge,
                                                   shortest_white_path=length))
    assert max(row['shortest_white_path']for row in quadrilaterals)==3
    assert len(quadrilaterals)==len(source['local_face_graphs']['quadrilateral_cases'])==212

    # Exact h-vector coefficient minimum at29 and symbolic continuation.
    coeff=[comb(29-i,2)+comb(i,2)for i in range(15)]
    assert min(coeff)==196 and coeff[0]>196
    assert [c-196 for c in coeff]==source['scalar']['h29_coefficients']
    assert F(5*196,2)-comb(29,2)==84
    t=s.symbols('t',nonnegative=True,integer=True)
    odd=s.expand(((29+2*t)-1)*((29+2*t)-5)/8-(3*(29+2*t)-4))
    even=s.expand((30+2*t)*((30+2*t)-6)/8-(3*(30+2*t)-4))
    assert odd==t*t/2+7*t+1 and even==t*t/2+s.Rational(15,2)*t+4

    # An order>=4 minimal-indefinite connected diagram cannot have a label>=6:
    # that edge has an adjacent nonright edge and the resulting proper triple
    # already has nonpositive determinant. Extra third edges only worsen it.
    a,b,c=s.symbols('a b c',nonnegative=True)
    det3=s.Matrix([[1,-a,-c],[-a,1,-b],[-c,-b,1]]).det()
    assert s.expand(det3-(1-a*a-b*b-c*c-2*a*b*c))==0
    x,y,z=s.symbols('x y z',nonnegative=True)
    obstruction=s.Poly(s.expand(det3.subs({a:s.sqrt(3)/2+x,b:s.Rational(1,2)+y,c:z})),x,y,z)
    assert obstruction.coeff_monomial(1)==0
    assert all(v<=0 for v in obstruction.coeffs())
    assert gcd(3,4)==1
    # The cubic-seven generator polynomial is irreducible and totally real.
    X=s.symbols('X');f=s.Poly(X**3-X**2-2*X+1,X)
    assert f.is_irreducible and s.discriminant(f.as_expr(),X)==49
    assert s.Poly(f).count_roots(-s.oo,s.oo)==3

    paths=[f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_COMPACT_ARITHMETIC_DISTANCE_THREE_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                producer_status=source['status'],graph_checks=graph_rows,
                connected_all_rank_formulas=dict(path='3r-6 for r>=3',D='3r-5 for r>=5; D4=6',
                                                 E='3r-4 at r6,7,8'),
                critical_rank29_forest_max=dp[29],uniform_rank29_upper=83,
                triangle_topologies=triangles,quadrilateral_topologies=quadrilaterals,
                h_vector=dict(paired_coefficients=coeff,average_lower_strictly_greater_than=84,
                              odd_margin=str(odd),even_margin=str(even)),
                arithmetic=dict(triple_determinant=str(det3),large_label_obstruction=str(obstruction.as_expr()),
                                cubic_seven_polynomial=str(f.as_expr()),discriminant=49,
                                determinant_field='Q(sqrt2,sqrt5)',
                                cycle_field_argument='Permutation monomials are products of Gram cycles; determinant belongs to the ground field.',
                                general_field_argument='Vinberg Proposition17, original printed p91, checked directly: k=Q(det L). Alternatively a proper extension k/Q(det L) has a nonidentity real embedding fixing the negative determinant, contradicting conjugate PSD.',
                                intersection_argument='Odd-degree k intersects the degree4 multiquadratic field only in Q.',
                                conjugate_argument='A negative rational principal determinant stays negative under a nonidentity ground-field embedding, contradicting conjugate PSD.'),
                source_proof_scope=dict(compactness='Required globally to exclude affine subdiagrams and make all proper missing-face sets spherical.',
                                        triangle='A minimal indefinite star subdiagram contains all three black roots, hence is exactly their connected triple.',
                                        quadrilateral='Opposite-pair subdiagrams contain minimal indefinite diagrams with <=1 white node each; index1 forces overlap or a cross-edge.',
                                        arithmetic='Any arithmetic ground field not contained in Q(sqrt2,sqrt5), including every nontrivial odd-degree field. Only the determinant-generation consequence of admissibility is used; no integrality census or fixed-lattice assumption.',
                                        omitted='No independent enumeration of the producer optional seven Lanner triangles; not needed for the dimension proof.'),
                conclusion='Compact Coxeter polytopes with every Lanner diagram of order<=3 have dimension<=28. The arithmetic corollary applies to every ground field not contained in Q(sqrt2,sqrt5), in particular cubic-seven. No unrestricted compact28 claim.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['graph_checks']),len(result['quadrilateral_topologies']))
