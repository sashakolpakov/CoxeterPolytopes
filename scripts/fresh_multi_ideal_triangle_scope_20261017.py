#!/usr/bin/env python3
"""Exact Q(sqrt2) Coxeter counterexamples and bounded-block identities.

No imported producer arithmetic. The geometric star argument and its published
rank-ten classification input are proved/scoped in the companion report.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_multi_ideal_triangle_scope_20261017'


@dataclass(frozen=True)
class Q:
    a:F=F(0)
    b:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a)); object.__setattr__(self,'b',F(self.b))
    def __add__(self,x):
        x=coerce(x); return Q(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self): return Q(-self.a,-self.b)
    def __sub__(self,x): return self+-coerce(x)
    def __rsub__(self,x): return coerce(x)+-self
    def __mul__(self,x):
        x=coerce(x); return Q(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=coerce(x); den=x.a*x.a-2*x.b*x.b
        assert den
        return self*Q(x.a/den,-x.b/den)
    def __bool__(self): return bool(self.a or self.b)
    def sign(self):
        if not self.a: return (self.b>0)-(self.b<0)
        if not self.b or self.a*self.b>0: return (self.a>0)-(self.a<0)
        den=self.a*self.a-2*self.b*self.b
        assert den
        return ((self.a>0)-(self.a<0))*((den>0)-(den<0))
    def data(self): return [str(self.a),str(self.b)]


def coerce(x): return x if isinstance(x,Q) else Q(x)


def dot(a,b): return -a[0]*b[0]+sum((x*y for x,y in zip(a[1:],b[1:])),Q())


def solve(a,b):
    m=[[coerce(x) for x in row]+[coerce(y)] for row,y in zip(a,b)]
    n=len(b)
    for i in range(n):
        j=next((j for j in range(i,n) if m[j][i]),None)
        if j is None:return None
        m[i],m[j]=m[j],m[i]
        den=m[i][i];m[i]=[x/den for x in m[i]]
        for j in range(n):
            if j!=i:
                c=m[j][i];m[j]=[x-c*y for x,y in zip(m[j],m[i])]
    return [row[-1] for row in m]


def det(a):
    n=len(a)
    if n==0:return Q(1)
    return sum(((-1)**j*a[0][j]*det([[a[i][k] for k in range(n) if k!=j]
                                   for i in range(1,n)]) for j in range(n)),Q())


def pyramid(labels):
    c={3:Q(F(1,2)),4:Q(0,F(1,2))}
    side=[]
    for m,(x,y) in zip(labels,[(1,0),(0,1),(-1,0),(0,-1)]):
        side.append((c[m],Q(x),Q(y),c[m]))
    N=[tuple(map(Q,(0,0,0,-1)))]+side
    G=[[dot(a,b) for b in N] for a in N]
    assert all(G[i][i]==Q(1) for i in range(5))
    cycles=[]
    for i,j in itertools.combinations(range(5),2):
        assert G[i][j].sign()<=0
        cyc=4*G[i][j]*G[i][j]
        assert cyc in list(map(Q,[0,1,2,4]))
        cycles.append([i,j,cyc.data()])
    H=tuple(map(Q,(2,0,0,1)))
    assert dot(H,H)==Q(-3) and all(dot(H,r).sign()<0 for r in N)
    vertices=set()
    for ids in itertools.combinations(range(5),3):
        p=solve([N[i][1:] for i in ids],[N[i][0] for i in ids])
        if p is None:continue
        v=tuple([Q(1)]+p)
        if all(dot(v,r).sign()<=0 for r in N):vertices.add(v)
    vertices=sorted(vertices,key=lambda v:tuple((x.a,x.b) for x in v))
    assert len(vertices)==5
    active=[{i for i,r in enumerate(N) if not dot(v,r)} for v in vertices]
    faces=[[k for k,a in enumerate(active) if i in a] for i in range(5)]
    edges=set()
    for i,j in itertools.combinations(range(5),2):
        inds=tuple(k for k,a in enumerate(active) if i in a and j in a)
        if len(inds)==2:edges.add(inds)
    assert len(edges)==8 and sorted(map(len,faces))==[3,3,3,3,4]
    assert all(dot(v,v).sign()<=0 for v in vertices)
    nu=next(i for i,v in enumerate(vertices) if v==tuple(map(Q,(1,0,0,1))))
    v=next(i for i,a in enumerate(active) if a=={0,1,2})
    assert not dot(vertices[v],vertices[v])
    assert active[nu]=={1,2,3,4}
    assert tuple(sorted((v,nu))) in edges
    # All ideal supports are PSD. Positive cusp relations below give bounded links.
    for p,a in zip(vertices,active):
        if not dot(p,p):
            ids=sorted(a)
            for k in range(1,len(ids)+1):
                for js in itertools.combinations(ids,k):
                    assert det([[G[i][j] for j in js] for i in js]).sign()>=0
    R2=Q(0,1)
    ray=[R2*N[0][i]+N[1][i]+N[2][i] for i in range(4)]
    assert ray==[R2*x for x in vertices[v]]
    assert all(N[1][i]+N[3][i]==(c[labels[0]]+c[labels[2]])*vertices[nu][i]
               for i in range(4))
    assert all(N[2][i]+N[4][i]==(c[labels[1]]+c[labels[3]])*vertices[nu][i]
               for i in range(4))
    idealcounts=[sum(not dot(vertices[k],vertices[k]) for k in face) for face in faces]
    expected=2 if labels==[4,4,3,3] else 3
    assert len(faces[1])==len(faces[2])==3
    assert idealcounts[1]==idealcounts[2]==expected
    assert len(faces[0])==4
    if expected==2:assert sum(not dot(p,p) for p in vertices)==2
    else:assert all(not dot(p,p) for p in vertices)
    return dict(base_side_labels=labels,normal_coordinates=[[x.data() for x in r] for r in N],
                gram=[[x.data() for x in r] for r in G],cartan_two_cycles=cycles,
                strict_center=[x.data() for x in H],
                vertices=[[x.data() for x in p] for p in vertices],
                vertex_norms=[dot(p,p).data() for p in vertices],
                active_sets=[sorted(a) for a in active],faces=faces,edges=[list(e) for e in sorted(edges)],
                ideal_counts_by_face=idealcounts,chosen_simple_cusp=v,nonsimple_apex=nu,
                two_triangles=[1,2],nontriangle_twoface=0,
                failed_triangle_scope=f'exactly_{expected}_ideal_vertices')


def block_maxima():
    # Exact DP checks the closed formula; the report proves it for arbitrary m.
    dp=[0]
    rows=[]
    for m in range(1,257):
        dp.append(max([dp[m-1]]+[dp[m-k]+math.comb(k,2) for k in range(2,min(9,m)+1)]))
        formula=36*(m//9)+math.comb(m%9,2)
        assert dp[m]==formula
        rows.append(dict(ports=m,maximum_triangles=formula,
                         two_triangle_payment_per_port=str(F(2*formula,m))))
    return rows


def run():
    paths=[f'scripts/{TAG}.py',
           'results/fresh_finite_volume_simplex_truncation_20261017.md',
           'results/fresh_cusp_one_ideal_simplex_capacity_20261017.md',
           'results/fresh_finite_ended_simplex_ports_20261017.md',
           'results/fresh_nikulin_narrow_20261015_compact_facet.md']
    return dict(status='PASS_EXACT_MULTI_IDEAL_COUNTEREXAMPLES_AND_SIMPLEX_BLOCK_FORMULA',
                scope='Two literal finite-volume Coxeter counterexamples; original-star simplex-face bound9; no ambient dimension claim.',
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                examples=[pyramid([4,4,3,3]),pyramid([4,4,4,4])],
                block_maxima=block_maxima(),
                classification_input=dict(maximum_minimal_indefinite_order=10,
                    primary_url='https://www.mathnet.ru/eng/im1484',
                    primary_pdf='https://www.mathnet.ru/links/1e830e3298f5df6e148b8dc52fd64e1c/im1484_eng.pdf',
                    location='Section2, Proposition1, printed p407',
                    independently_reclassified=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==data
    else:p.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(data['status'],'saved_replay=',args.verify)
