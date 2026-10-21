#!/usr/bin/env python3
"""Exact finite checks for the all-white closure and Qsqrt2 weight theorem.

The unbounded attachment and face-counting arguments are proved in the report.
They are not inferred from a bounded attachment search.
"""
import argparse,importlib.util,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product,permutations
from pathlib import Path
import networkx as nx
import sympy as s
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt2_quad_exhaustion_20261015';TAG=BASE+'_closure'
sp=importlib.util.spec_from_file_location('union_supports',ROOT/'scripts'/f'{BASE}.py')
E=importlib.util.module_from_spec(sp);sp.loader.exec_module(E)

def graph(row):
    G=nx.Graph();G.add_nodes_from((i,dict(black=i<4)) for i in range(row['nodes']))
    for p,label in zip(combinations(range(row['nodes']),2),row['labels']):G.add_edge(*p,label=label,opposite=p in ((0,2),(1,3)))
    return G

def iso(left,right):
    gm=nx.algorithms.isomorphism.GraphMatcher(graph(left),graph(right),node_match=lambda a,b:a==b,edge_match=lambda a,b:a==b)
    assert gm.is_isomorphic()
    mp=gm.mapping;L=graph(left);R=graph(right)
    assert len(mp)==left['nodes'] and len(set(mp.values()))==left['nodes']
    assert all(L.edges[u,v]==R.edges[mp[u],mp[v]] for u,v in L.edges)
    return [mp[i] for i in range(left['nodes'])]

def add_path(row,length,ends):
    n=row['nodes'];g=dict(zip(combinations(range(n),2),row['labels']));N=n+length
    gg={p:g.get(p,2) for p in combinations(range(N),2)}
    path=[ends[0]]+list(range(n,N))+(list(ends[1:]))
    for a,b in zip(path,path[1:]):gg[E.pair(a,b)]=3
    return dict(nodes=N,white=row['white']+list(range(n,N)),labels=[gg[p] for p in combinations(range(N),2)])

def sparse_graph(nodes,g):
    G=nx.Graph();G.add_nodes_from(nodes)
    for a,b in combinations(sorted(nodes),2):
        if g[(a,b)]!=2:G.add_edge(a,b,label=g[(a,b)])
    return G

def extra_pairs(H):
    n=len(H);deg=dict(H.degree());extra={};kind=None
    if n==4 and sorted(deg.values())==[1,1,2,2] and sorted(nx.get_edge_attributes(H,'label').values())==[3,3,4]:
        ends=sorted(v for v in H if deg[v]==1)
        internal=next(E.pair(a,b) for a,b,d in H.edges(data=True) if d['label']==4)
        if all(deg[v]==2 for v in internal):extra[E.pair(*ends)]=F(1);kind='F4'
    branch=[v for v in H if deg[v]==3]
    if len(branch)==1 and max(deg.values())==3 and all(d['label']==3 for a,b,d in H.edges(data=True)):
        c=branch[0];arms=[]
        for v in H[c]:
            arm=[v];prev,u=c,v
            while deg[u]==2:
                nxt=next(x for x in H[u] if x!=prev);arm.append(nxt);prev,u=u,nxt
            arms.append(arm)
        lengths=sorted(map(len,arms))
        if lengths==[1,1,2]:
            tip=next(a[-1] for a in arms if len(a)==2)
            for a in arms:
                if len(a)==1:extra[E.pair(tip,a[0])]=F(1,2)
            kind='D5'
        elif lengths==[1,2,2]:
            a,b=[a for a in arms if len(a)==2]
            extra[E.pair(a[0],b[-1])]=extra[E.pair(b[0],a[-1])]=F(1,2);kind='E6'
        elif lengths in ([1,2,3],[1,2,4]):
            long=max(arms,key=len);mid=next(a for a in arms if len(a)==2)
            extra[E.pair(long[0],mid[-1])]=F(1,4);kind='E'+str(n)
    return extra,kind

def corner_weights(row):
    g=dict(zip(combinations(range(row['nodes']),2),row['labels']));out=[]
    for a,b in product((0,2),(1,3)):
        G=sparse_graph(row['white']+[a,b],g);weight=F(0);kind=None;distance=None
        for nodes in nx.connected_components(G):
            H=G.subgraph(nodes)
            if a in H and b in H:
                distance=nx.shortest_path_length(H,a,b)
                if distance<=2:weight+=1
                extras,kind=extra_pairs(H);weight+=extras.get(E.pair(a,b),0)
        out.append(dict(marked=[a,b],distance=distance,exceptional_type=kind,weight=str(weight)))
    assert sum(F(r['weight']) for r in out)>=1
    return out

def finite_graph(kind,n):
    if kind in ('A','B','I','F'):
        edges=[(i,i+1,4 if kind=='B' and i==0 or kind=='F' and i==1 else 3) for i in range(n-1)]
    else:
        arms=[1,1,n-3] if kind=='D' else [1,2,n-4];edges=[];nxt=1
        for length in arms:
            prev=0
            for _ in range(length):edges.append((prev,nxt,3));prev=nxt;nxt+=1
        assert nxt==n
    G=nx.Graph();G.add_nodes_from(range(n));G.add_edges_from((a,b,dict(label=l)) for a,b,l in edges);return G

def component_weights():
    out=[]
    for kind,n in [('A',n) for n in range(1,41)]+[('B',n) for n in range(2,41)]+[('D',n) for n in range(4,41)]+[('F',4)]+[('E',n) for n in (6,7,8)]:
        G=finite_graph(kind,n);base=sum(nx.shortest_path_length(G,a,b)<=2 for a,b in combinations(G,2))
        extras,_=extra_pairs(G);bonus=sum(extras.values(),F(0));weight=base+bonus
        expected=0 if n==1 else 2*n-3 if kind in ('A','B','F') else 2*n-2
        assert base==expected
        deficit=2*n-weight;assert deficit>=1
        if n>=9:assert deficit>=2
        out.append(dict(type=kind+str(n),rank=n,base=base,extra=str(bonus),deficit=str(deficit),
                        exceptional_pairs=[[a,b,str(v)] for (a,b),v in sorted(extras.items())]))
    return out

def run():
    base=json.loads((ROOT/'results'/f'{BASE}.json').read_text())
    core=json.loads((ROOT/'results'/f'{BASE}_cores.json').read_text())
    ext=json.loads((ROOT/'results'/f'{BASE}_extensions.json').read_text())
    tail=json.loads((ROOT/'results/fresh_sqrt2_dotted_tail_20261015.json').read_text())
    byclass={i:base['labelled_union_supports'][r['representative']] for i,r in enumerate(base['marked_isomorphism_classes'])}
    erows={r['class_index']:r for r in ext['records']}
    sealed=[r['class_index'] for r in core['records'] if r['status']=='PASS_PRINCIPAL_UNION_NECESSARY_TESTS' and not r['physical']['determinant']]
    assert sealed==[5,9,17,32,41,47,50]
    assert all(not erows[i]['accepted'] for i in sealed)
    assert erows[14]['accepted']==[[2,2,3,2,2,2,2]]
    assert set(map(tuple,erows[37]['accepted']))=={(2,2,3,2,2,2,2),(3,2,2,2,2,2,2),(3,2,3,2,2,2,2)}
    maps=[dict(source_core=14,path_length=1,ends=[2],target_core=17,
               bijection=iso(add_path(byclass[14],1,[2]),byclass[17]))]
    for length,target in [(1,41),(2,47),(3,50)]:
        maps.append(dict(source_core=37,path_length=length,ends=[0,2],target_core=target,
                         bijection=iso(add_path(byclass[37],length,[0,2]),byclass[target])))
    # The finite D5/E6/E7/E8 arm classification is all-size. These exact
    # checks anchor its four finite cases and first affine obstruction.
    assert [r['inertia'] for r in tail['corner_tails']]==[[5,0,0],[6,0,0],[7,0,0],[8,0,0],[8,0,1]]
    p,q=s.symbols('p q');z=s.sqrt(2);t=(1-z)/2
    core37_poly=s.expand(z*z+2*t*z-4*p*q+2*t*(p+q))
    assert s.expand(core37_poly-(p+q-4*p*q+(1-p-q)*z))==0
    tail_values=[s.Rational(r['reduced_unit_diagonal']) for r in tail['tail_diagonals'][:4]]
    rational_tails=[]
    for a,b in product(tail_values,repeat=2):
        val=s.expand(core37_poly.subs({p:a,q:b}));assert val!=0 and a+b>1
        rational_tails.append(dict(p=str(a),q=str(b),rank_polynomial=str(val)))
    # Connected exterior components lower BOTH endpoint diagonals by >=1/4
    # and increase the effective opposite contact. The upper bounds below
    # are strictly below the original dotted threshold2 (norm-two roots).
    single_bound=s.Rational(3,4)*(2+z)
    double_bound=s.Rational(1,2)+z
    assert single_bound<4 and double_bound<2
    # Direct full determinant check of the parent double-family formula.
    row=byclass[37];G=s.eye(row['nodes']);G[0,0]=p;G[2,2]=q
    values={2:0,3:s.Rational(1,2),4:z/2}
    for (a,b),lab in zip(combinations(range(row['nodes']),2),row['labels']):G[a,b]=G[b,a]=-values[lab]
    determinant=0
    for perm in permutations(range(row['nodes'])):
        terms=[G[i,perm[i]] for i in range(row['nodes'])]
        if any(v==0 for v in terms):continue
        parity=sum(perm[i]>perm[j] for i in range(row['nodes']) for j in range(i+1,row['nodes']))
        determinant+=(-1)**parity*s.prod(terms)
    factor=(1+z)/64
    assert s.expand(determinant-factor*core37_poly)==0
    weights={str(i):corner_weights(byclass[i]) for i in sealed}
    # The dotted parent has already certified entire D5/E6 marked components.
    parent_weights=[]
    for star in tail['allowed_local_stars']:
        n=len(star['gram']);G=[[s.Rational(a)+s.Rational(b)*z for a,b in row] for row in star['gram']]
        cweights=[]
        for corner in star['corners']:
            nodes=corner['support'];a,b=corner['marked'];H=nx.Graph();H.add_nodes_from(nodes)
            for i,j in combinations(nodes,2):
                if G[i][j]!=0:
                    sq=s.simplify(G[i][j]**2/(G[i][i]*G[j][j]));assert sq==s.Rational(1,4)
                    H.add_edge(i,j,label=3)
            w=F(0);kind=None
            for comp in nx.connected_components(H):
                J=H.subgraph(comp)
                if a in J and b in J:
                    d=nx.shortest_path_length(J,a,b);w+=int(d<=2);bonus,kind=extra_pairs(J);w+=bonus.get(E.pair(a,b),0)
            cweights.append(dict(marked=[a,b],type=kind,weight=str(w)))
        assert sum(F(x['weight']) for x in cweights)>=1
        parent_weights.append(dict(double=star['double'],tail_length=star['ell'],corners=cweights))
    triangle_checks=[]
    for labs in E.TRIANGLES:
        assert sum(v!=2 for v in labs)>=2
        triangle_checks.append(dict(kind='L3',labels=list(labs),short_paths=sum(v!=2 for v in labs)))
    for n in (4,5):
        for black in combinations(range(n),3):
            gaps=[(black[(i+1)%3]-black[i])%n for i in range(3)]
            assert sum(g<=2 for g in gaps)>=2
            triangle_checks.append(dict(kind='C'+str(n),black=list(black),gaps=gaps,short_paths=sum(g<=2 for g in gaps)))
    N=s.symbols('n',integer=True)
    even=s.expand(N*(N-6)/8-(2*N-2));odd=s.expand((N-1)*(N-5)/8-(2*N-2))
    assert s.expand(even-(N*(N-22)+16)/8)==0 and s.factor(odd)==(N-21)*(N-1)/8
    files=[f'scripts/{TAG}.py']+[f'{d}/{BASE}{sf}.{extn}' for sf in ('','_cores','_extensions') for d,extn in [('scripts','py'),('results','json')]]
    files += [f'{d}/fresh_sqrt2_dotted_tail_20261015{sf}.{extn}' for sf in ('','_review') for d,extn in [('scripts','py'),('results','json')]]
    files += ['results/fresh_compact_dimension_bound_20261015.json','results/fresh_compact_dimension_bound_20261015_review.json']
    return dict(status='PASS_EXACT_CHECKS_FOR_COMPLETE_SQRT2_QUADRILATERAL_CLOSURE_AND_BOUND20_PENDING_INDEPENDENT_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                sealed_finite_cores=sealed,exact_core_extension_maps=maps,
                remaining_core37_rank=dict(polynomial=str(core37_poly),determinant_factor=str(factor),tail_cases=rational_tails),
                connected_dotted_component_bounds=dict(single_squared=str(single_bound),double_contact=str(double_bound),
                    monotonicity='positive inverse of white forest; p,q<=3/4; effective opposite contact increases'),
                finite_core_corner_weights=weights,parent_tail_corner_weights=parent_weights,
                triangle_checks=triangle_checks,component_weight_checks=component_weights(),
                all_rank_identities=dict(even_lower_minus_upper=str(even),odd_lower_minus_upper=str(odd),
                    vertex_upper='2n-2 for n>=9',contradiction='strict face lower bound for every n>=21'),
                proof_scope='Compact arithmetic Coxeter polytopes over Qsqrt2 only. Full external-white closure and h-vector argument are in the report; exact census is not substituted for those unbounded proofs.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
