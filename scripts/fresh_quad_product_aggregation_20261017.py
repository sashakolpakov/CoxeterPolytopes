#!/usr/bin/env python3
"""Canonical full-star quadrilateral aggregation and overlapping product faces."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_quad_product_aggregation_20261017'
DOTTED='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017'
WITNESS='fresh_sqrt5_adjacent_discount_witness_20261016'
K,Z=A.K,A.Z
FINITE={C.code(g):name for name,g in C.models()}


def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def pd_forest(ids,edges):
    names=[]
    for g in C.pieces(C.graph(ids,edges)):
        if sum(map(len,g.values()))!=2*len(g)-2:return None
        if any(m==0 for row in g.values()for m in row.values()):return None
        name=FINITE.get(C.code(g))
        if name is None:return None
        names.append(name)
    return sorted(names)


def quads(P,N,E):
    return [(p,q)for p in combinations(sorted(P),2)for q in combinations(sorted(N),2)
            if all((i,j)in E for i in p for j in q)]


def maximal_bicliques(P,N,E):
    ps=sorted(P);ans=set()
    for mask in range(1<<len(ps)):
        I={p for i,p in enumerate(ps)if mask>>i&1}
        if len(I)<2:continue
        J={q for q in N if all((p,q)in E for p in I)}
        if len(J)<2:continue
        I={p for p in P if all((p,q)in E for q in J)}
        ans.add((tuple(sorted(I)),tuple(sorted(J))))
    return sorted(ans)


def dotted_census():
    prior=read(DOTTED);assert prior['counts']['cores']==1399
    records=[];counts=Counter();dimension_counts=Counter()
    for source in prior['records']:
        n=source['nodes'];labels=source['row']['labels']
        assert set(labels)<={0,2,3,4,5,6,10}
        edges=dict(zip(combinations(range(n),2),labels))
        k=[K(*x)for x in source['kernel']]
        signs=[x.sign()for x in k]
        P={i for i,s in enumerate(signs)if s>0};N={i for i,s in enumerate(signs)if s<0}
        zero={i for i,s in enumerate(signs)if not s}
        pd=[];E=set()
        for i,j in combinations(range(n),2):
            forest=pd_forest([a for a in range(n)if a not in(i,j)],edges)
            if forest is None:continue
            assert signs[i]*signs[j]==-1
            edge=(i,j)if i in P else(j,i)
            E.add(edge);pd.append(dict(omitted=[i,j],forest=forest))
        Q=quads(P,N,E);assert (tuple(sorted(P&set(range(4)))),tuple(sorted(N&set(range(4)))))in Q
        dots=[p for p,m in edges.items()if m==0]
        assert dots==list(map(tuple,source['dotted_pairs']))and dots
        assert all(set(d)in[{0,2},{1,3}]for d in dots)
        d=set(dots[0]);assert len({signs[i]for i in d})==1 and 0 not in{signs[i]for i in d}
        assert all(d&set(e)for e in E)
        if d<=P:
            common={q for q in N if all((p,q)in E for p in d)};prism=(tuple(sorted(d)),tuple(sorted(common)))
        else:
            common={p for p in P if all((p,q)in E for q in d)};prism=(tuple(sorted(common)),tuple(sorted(d)))
        ksize=len(common);assert ksize>=2
        maxes=maximal_bicliques(P,N,E);assert maxes==[prism]
        assert len(Q)==ksize*(ksize-1)//2
        for I,J in Q:assert set(I)<=set(prism[0])and set(J)<=set(prism[1])
        records.append(dict(index=source['index'],source_index=source['source_index'],nodes=n,
            kernel_signs=signs,positive=sorted(P),negative=sorted(N),zero=sorted(zero),
            positive_definite_two_deletions=pd,pd_deletion_edges=sorted(E),
            dotted_base=sorted(d),maximal_prism_factors=[list(x)for x in prism],
            prism_dimension=ksize,prism_vertex_count=2*ksize,
            quadrilateral_count=len(Q),quadrilateral_black_partitions=[[list(a),list(b)]for a,b in Q],
            aggregate_corner_count=4*len(Q)))
        counts['stars']+=1;counts['pd_deletions']+=len(E);counts['quadrilaterals']+=len(Q)
        counts['aggregate_corners']+=4*len(Q);dimension_counts[str(ksize)]+=1
    return dict(counts=dict(counts),prism_dimension_counts=dict(sorted(dimension_counts.items())),records=records,
        scope='All1399 saved marked dotted local-star records, with every two-root deletion recognized by the complete finite Coxeter catalogue. Quadrilateral totals are sums over records; no cross-record full-star isomorphism quotient is claimed. Exact field Gram/kernel/admissibility proofs are inherited from the frozen source. Orthogonal white spectators leave this deletion graph and prism unchanged.')


def exact_overlap():
    saved=read(WITNESS)['geometry'];phi=K(A.F(1,2),A.F(1,2))
    G=[[K(2*(i==j))for j in range(8)]for i in range(8)];G[6][6]=K(4)
    for i,j,m in[(0,5,5),(1,3,5),(2,3,3),(2,5,3),(3,7,3),(4,5,3),(6,7,4)]:
        G[i][j]=G[j][i]=-{3:K(1),4:K(2),5:phi}[m]
    assert A.enc(G)==saved['scaled_gram']
    assert A.inertia(G)==[6,1,1]and A.inertia([[x.conj()for x in r]for r in G])==[7,0,1]
    assert all((2*G[i][j]/G[j][j]).integral()for i in range(8)for j in range(8))
    kernel=A.kernel(G);assert len(kernel)==1;k=kernel[0]
    assert all(not x for x in A.mv(G,k))
    P={i for i,x in enumerate(k)if x.sign()>0};N={i for i,x in enumerate(k)if x.sign()<0}
    assert P=={1,3,6,7}and N=={0,2,4,5}
    basis=list(range(1,8));Q=A.sub(G,basis);Qi=A.inverse(Q)
    normals=[A.mv(Qi,[G[j][i]for j in basis])for i in range(8)]
    assert [[A.pair(Q,x,y)for y in normals]for x in normals]==G
    H=A.mv(Qi,[-K(1)]*7);hp=[A.pair(Q,H,n)for n in normals]
    assert A.pair(Q,H,H).sign()<0 and all(x.sign()<0 for x in hp)
    E=set();deletions=[];vertices={}
    for i,j in combinations(range(8),2):
        active=[a for a in range(8)if a not in(i,j)];B=A.sub(G,active);inertia=A.inertia(B)
        rec=dict(omitted=[i,j],active=active,inertia=inertia)
        if inertia==[6,0,0]:
            assert k[i].sign()*k[j].sign()==-1
            edge=(i,j)if i in P else(j,i);E.add(edge)
            coeff=A.mv(A.inverse(B),[hp[a]for a in active])
            v=[H[z]-sum((c*normals[a][z]for c,a in zip(coeff,active)),Z)for z in range(7)]
            ips=[A.pair(Q,v,n)for n in normals]
            assert all(not ips[a]for a in active)and ips[i].sign()<0 and ips[j].sign()<0
            norm,future=A.pair(Q,v,v),A.pair(Q,v,H);assert norm.sign()<0 and future.sign()<0
            rec.update(point=[x.enc()for x in v],norm=norm.enc(),future=future.enc(),products=[x.enc()for x in ips])
            vertices[edge]=v
        deletions.append(rec)
    maxes=maximal_bicliques(P,N,E)
    expected=[((1,3),(0,2,5)),((1,3,7),(2,5))];assert maxes==expected
    allquads=quads(P,N,E);assert len(allquads)==5
    products=[];quadsets=[]
    for I,J in maxes:
        white=sorted(set(range(8))-set(I)-set(J));assert A.inertia(A.sub(G,white))==[3,0,0]
        vv=[vertices[(i,j)]for i in I for j in J]
        assert all(A.pair(Q,x,y).sign()<0 for x,y in combinations(vv,2))
        qs=quads(I,J,E);quadsets.append(set(qs));assert len(qs)==3
        products.append(dict(factors=[list(I),list(J)],white=white,dimension=3,combinatorial_type='Delta1 x Delta2',
            vertex_deletion_pairs=[[i,j]for i in I for j in J],quadrilaterals=[[list(a),list(b)]for a,b in qs]))
    shared=quadsets[0]&quadsets[1];assert shared=={((1,3),(2,5))}
    assert A.inertia(A.sub(G,[0,4,6,7]))==[4,0,0]
    assert set.union(*quadsets)==set(allquads)
    assert next(r['inertia']for r in deletions if r['omitted']==[0,7])==[5,1,0]
    aggregate=[]
    for i,j in sorted(E):
        for a in sorted(P-{i}):
            for b in sorted(N-{j}):
                if all(e in E for e in[(i,b),(a,j),(a,b)]):
                    aggregate.append(dict(vertex_omitted=[i,j],marked_corner_roots=[a,b]))
    assert len(aggregate)==20
    return dict(core='4202',ambient_dimension=6,gram=A.enc(G),kernel=[x.enc()for x in k],
        basis_indices=basis,basis_gram=A.enc(Q),root_coordinates=A.enc(normals),
        strict_center=[x.enc()for x in H],center_norm=A.pair(Q,H,H).enc(),center_products=[x.enc()for x in hp],
        physical_inertia=[6,1,1],conjugate_inertia=[7,0,1],integral_cartan_entries=64,
        positive=sorted(P),negative=sorted(N),all_two_root_deletions=deletions,pd_deletion_edges=sorted(E),
        maximal_products=products,shared_quad=dict(positive=[1,3],negative=[2,5],white=[0,4,6,7]),
        missing_union_vertex=dict(omitted=[0,7],inertia=[5,1,0]),
        distinct_quadrilaterals=[[list(a),list(b)]for a,b in allquads],aggregate_corner_occurrences=aggregate,
        aggregate_required_payment=5,incorrect_sum_of_maximal_product_demands=6,
        scope='Two actual compact triangular-prism faces overlap in a quadrilateral in this local acute arithmetic H6 chamber. Compactness of the entire chamber, or a compact completion containing it, is not asserted. This disproves uniqueness from the local corank1/acute/finite-corner premises, not a stronger theorem requiring an independently established compact ambient realization.')


def run():
    dotted=dotted_census();overlap=exact_overlap()
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{C.TAG}.py',
        f'results/{DOTTED}.json',f'scripts/{DOTTED}.py',f'results/{WITNESS}.json',
        'results/fresh_singular_quad_local_realization_20261017.md',
        'results/fresh_sqrt5_rank13_triangle_simplex_flags_20261017.md']
    return dict(status='PASS_FULL_STAR_QUADRILATERAL_AGGREGATION_AND_LOCAL_PRODUCT_OVERLAP',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dotted_prism_census=dotted,
        exact_finite_overlap=overlap,
        theorem='For the full n+2-wall corank1 star of a compact quadrilateral in an acute hyperbolic chamber, PD two-root deletions form a bipartite graph on the nonzero kernel signs. Its K2,2 subgraphs are in bijection with all actual compact quadrilaterals having that same full wall set. Consequently aggregate payment at least the number of distinct K2,2, with each quadrilateral corner counted once, is globally sufficient when applied to every actual full star. A dotted pair forces a unique maximal prism containing all these quadrilaterals; all-finite stars can have overlapping maximal products.',
        scope='General geometric aggregation lemma plus exact local tests. This does not supply a complete arithmetic face catalogue, a feasible weighting rule, a compact ambient example, or a dimension improvement.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(path.read_text());print('PASS saved full-star quad replay')
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['dotted_prism_census']['counts'],d['dotted_prism_census']['prism_dimension_counts'])
