#!/usr/bin/env python3
"""Independent marked-component weights for the Q(sqrt2) quad classification.

The global implication is conditional on the separately reviewed exhaustive
star classification. This checker does not replace its attachment proof.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import networkx as nx

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt2_weight_bound_20261015'
BASE='fresh_sqrt2_quad_exhaustion_20261015'

def graph(n,edges):
    g=nx.Graph();g.add_nodes_from(range(n))
    for e in edges:
        a,b,*m=e;g.add_edge(a,b,label=m[0] if m else 3)
    return g

def star(arms):
    edges=[];paths=[];next_node=1
    for length in arms:
        path=list(range(next_node,next_node+length));next_node+=length
        edges+=list(zip([0]+path,path));paths.append(path)
    return graph(next_node,edges),paths

def extras():
    out=[]
    g=graph(4,[(0,1,3),(1,2,4),(2,3,3)])
    out.append(('F4',g,{(0,3):F(1)}))
    g,p=star([1,1,2]);out.append(('D5',g,{tuple(sorted((p[i][0],p[2][-1]))):F(1,2) for i in (0,1)}))
    g,p=star([1,2,2]);out.append(('E6',g,{tuple(sorted((p[i][0],p[j][-1]))):F(1,2) for i,j in [(1,2),(2,1)]}))
    for name,arms in [('E7',[1,2,3]),('E8',[1,2,4])]:
        g,p=star(arms);out.append((name,g,{tuple(sorted((p[2][0],p[1][-1]))):F(1,4)}))
    return out

EXTRAS=extras()
def weights(g):
    out={p:F(int(nx.has_path(g,*p) and nx.shortest_path_length(g,*p)<=2)) for p in combinations(sorted(g),2)}
    for nodes in nx.connected_components(g):
        h=g.subgraph(nodes)
        for name,model,add in EXTRAS:
            if len(h)!=len(model):continue
            matcher=nx.isomorphism.GraphMatcher(h,model,edge_match=lambda a,b:a==b)
            maps=list(matcher.isomorphisms_iter())
            if not maps:continue
            assignments=[]
            for mp in maps:
                inverse={v:k for k,v in mp.items()}
                assignments.append({tuple(sorted((inverse[a],inverse[b]))):v for (a,b),v in add.items()})
            assert all(a==assignments[0] for a in assignments)
            for pair,value in assignments[0].items():
                assert out[pair]==0
                out[pair]+=value
    return out

def total(g):return sum(weights(g).values(),F())

def run():
    base=json.loads((ROOT/f'results/{BASE}.json').read_text())
    cores=json.loads((ROOT/f'results/{BASE}_cores.json').read_text())
    records=[]
    for rec in cores['records']:
        if rec['status']!='PASS_PRINCIPAL_UNION_NECESSARY_TESTS' or rec['physical']['determinant']:continue
        row=base['labelled_union_supports'][rec['union_index']]
        edges=[(a,b,m) for (a,b),m in zip(combinations(range(row['nodes']),2),row['labels']) if m!=2]
        g=graph(row['nodes'],edges);angles=[]
        for a,b in [(0,1),(1,2),(2,3),(3,0)]:
            h=g.subgraph(row['white']+[a,b]);value=weights(h)[tuple(sorted((a,b)))]
            angles.append(str(value))
        assert sum(map(F,angles))>=1
        records.append(dict(class_index=rec['class_index'],corner_weights=angles,total=str(sum(map(F,angles)))))
    assert len(records)==7
    # Rebuild the tail survivors as labelled diagrams, independently of their
    # scaled Gram representation in the primary certificate.
    tails=[]
    for length,double in [(0,False),(1,False),(0,True)]:
        n=7+length;edges=[(2,3,4),(2,4,3),(4,5,3),(5,6,3),(6,3,3),(0,5,3),(0,1,8 if length else 0)]
        if double:edges.append((1,5,3))
        if length:edges.append((0,7,3))
        g=graph(n,edges);W=[4,5,6]+list(range(7,n));angles=[]
        for a,b in [(0,2),(2,1),(1,3),(3,0)]:
            angles.append(str(weights(g.subgraph(W+[a,b]))[tuple(sorted((a,b)))]))
        assert sum(map(F,angles))>=1
        tails.append(dict(length=length,double=double,corner_weights=angles,total=str(sum(map(F,angles)))))
    finite=[]
    for r in range(1,61):
        A=graph(r,[(i,i+1,3) for i in range(r-1)])
        assert total(A)==max(0,2*r-3)
        finite.append(dict(type=f'A{r}',rank=r,total=str(total(A))))
        if r>=2:
            B=graph(r,[(i,i+1,4 if i==0 else 3) for i in range(r-1)])
            assert total(B)==2*r-3
            finite.append(dict(type=f'B{r}',rank=r,total=str(total(B))))
        if r>=4:
            D,_=star([1,1,r-3]);assert total(D)==2*r-2+(r==5)
            finite.append(dict(type=f'D{r}',rank=r,total=str(total(D))))
    for name,g,add in EXTRAS:
        t=total(g);expected={'F4':F(6),'D5':F(9),'E6':F(11),'E7':F(49,4),'E8':F(57,4)}[name]
        assert t==expected
        finite.append(dict(type=name,rank=len(g),total=str(t)))
    for m in (6,8):
        g=graph(2,[(0,1,m)]);assert total(g)==1
        finite.append(dict(type=f'I2({m})',rank=2,total='1'))
    # Independent finite-rank partition check accompanies the all-rank proof.
    best=[F(0)]+[F(-10000)]*60
    for r in range(1,61):
        best[r]=max(best[r-z['rank']]+F(z['total']) for z in finite if z['rank']<=r)
        if r>=9:assert best[r]<=2*r-2
    # A minimal Lanner diagram for a triangular face contains all 3 marks.
    # For orders4/5 it is a cycle; enumerate every marking and count gaps<=2.
    triangles=[]
    for n in (4,5):
        for marked in combinations(range(n),3):
            gaps=[(marked[(i+1)%3]-marked[i])%n for i in range(3)]
            assert sum(g<=2 for g in gaps)>=2
            triangles.append(dict(order=n,marks=list(marked),gaps=gaps))
    dimension=[]
    for n in range(21,121):
        lower=F(n*(n-6),8) if n%2==0 else F((n-1)*(n-5),8)
        assert lower>=2*n-2
        dimension.append(dict(dimension=n,strict_lower=str(lower),upper=2*n-2))
    # Symbolic all-dimension gap coefficients in u>=0 for n=21+2u/22+2u.
    for u in range(4):
        odd=F((20+2*u)*(16+2*u),8)-(40+4*u)
        even=F((22+2*u)*(16+2*u),8)-(42+4*u)
        assert odd==F(u*u,2)+5*u and even==F(u*u,2)+F(11*u,2)+2
    paths=[f'scripts/{TAG}.py',f'results/{BASE}.json',f'results/{BASE}_cores.json',
           'results/fresh_sqrt2_dotted_tail_20261015.json','results/fresh_compact_global_bound_20261015_review.md']
    return dict(status='PASS_QSQRT2_WEIGHTS_CONDITIONAL_ON_COMPLETE_STAR_CLASSIFICATION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        finite_core_weights=records,tail_weights=tails,finite_diagrams=finite,
        rank_maxima=[str(x) for x in best],cycle_triangle_marks=triangles,dimension_inequalities=dimension,
        scope='Weights and all-rank dimension implication independently checked. Bound20 requires the separate exhaustive star and arbitrary external attachment proof; this script alone does not establish it.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
