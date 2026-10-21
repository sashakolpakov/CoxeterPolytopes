#!/usr/bin/env python3
"""Independent exact symmetry and finite-corner checks; no producer imports."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import networkx as nx

ROOT=Path(__file__).resolve().parents[1];BASE='fresh_norm31_continuation_20261015';TAG=BASE+'_review'
def sg(x):return (x>0)-(x<0)
class Q:
    def __init__(self,a=0,b=0):
        if isinstance(a,Q):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,o):o=Q(o);return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,o):return self+-Q(o)
    def __rsub__(self,o):return Q(o)+-self
    def __mul__(self,o):o=Q(o);return Q(self.a*o.a+2*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=Q(o);d=o.a*o.a-2*o.b*o.b;assert d
        return self*Q(o.a/d,-o.b/d)
    def __eq__(self,o):o=Q(o);return self.a==o.a and self.b==o.b
    def __bool__(self):return bool(self.a or self.b)
    def conj(self):return Q(self.a,-self.b)
    def sign(self):
        if not self.b:return sg(self.a)
        if not self.a or sg(self.a)==sg(self.b):return sg(self.b)
        return sg(self.a)*sg(self.a*self.a-2*self.b*self.b)
    def integral(self):return self.a.denominator==self.b.denominator==1
    def enc(self):return [str(self.a),str(self.b)]
def dec(x):return Q(x['a'],x['b'])
def vec(v):return list(map(dec,v))
def mat(M):return list(map(vec,M))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q())
def mv(A,v):return [dot(row,v) for row in A]
def tr(A):return list(map(list,zip(*A)))
def mm(A,B):return [[dot(a,b) for b in tr(B)] for a in A]
def form(a,b,G):return dot(a,mv(G,b))
def sub(A,ids):return [[A[i][j] for j in ids] for i in ids]
def det(A):
    M=[r[:] for r in A];d=Q(1)
    for k in range(len(M)):
        j=next((j for j in range(k,len(M)) if M[j][k]),None)
        if j is None:return Q()
        if j!=k:M[k],M[j]=M[j],M[k];d=-d
        pivot=M[k][k];d*=pivot
        for i in range(k+1,len(M)):
            factor=M[i][k]/pivot
            M[i]=[x-factor*y for x,y in zip(M[i],M[k])]
    return d
def positive(A):
    M=[r[:] for r in A]
    while M:
        d=M[0][0]
        if d.sign()<=0:return False
        M=[[M[i][j]-M[i][0]*M[0][j]/d for j in range(1,len(M))] for i in range(1,len(M))]
    return True
def run():
    read=lambda suffix:json.loads((ROOT/f'results/{BASE}{suffix}.json').read_text())
    P,C,S=read(''),read('_corners'),read('_symmetry')
    seen={}
    def manifest(d):
        for p,h in d.get('sources',{}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
            if p in seen:assert seen[p]==h
            else:
                seen[p]=h
                if p.endswith('.json'):manifest(json.loads((ROOT/p).read_text()))
    for d in (P,C,S):manifest(d)
    G=mat(P['form']);N=mat(P['roots']);H=vec(P['basepoint']);center=vec(P['strict_center'])
    assert len(N)==28 and G[0][0]==Q(0,-2) and G[0][9]==G[9][0]==-1 and G[9][9]==2
    assert det(G)==Q(-1,-4) and positive([[x.conj() for x in row] for row in G])
    GG=[[form(a,b,G) for b in N] for a in N];assert GG==mat(P['gram'])
    assert form(center,center,G).sign()<0 and form(center,H,G).sign()<0
    for i,r in enumerate(N):
        q=GG[i][i];assert q in (Q(2),Q(18,10))
        assert all(x.integral() for x in r) and all((2*x/q).integral() for x in mv(G,r))
        assert all(GG[i][j].sign()<=0 for j in range(i))
        assert form(center,r,G).sign()<0
        point=vec(P['actual_facet_points'][i])
        assert form(point,point,G).sign()<0 and form(point,H,G).sign()<0
        assert all(form(point,n,G).sign()==(0 if i==j else -1) for j,n in enumerate(N))
    # An independent maximal-clique enumeration, followed by all 9-subsets,
    # covers every candidate finite corner (the producer uses a size recursion).
    graph=nx.Graph();graph.add_nodes_from(range(28))
    for i,j in combinations(range(28),2):
        if (GG[i][i]*GG[j][j]-GG[i][j]*GG[i][j]).sign()>0:graph.add_edge(i,j)
    supports=set()
    for clique in nx.find_cliques(graph):
        supports.update(combinations(sorted(clique),9))
    finite=sorted(ids for ids in supports if positive(sub(GG,ids)))
    assert len(supports)==35 and len(finite)==11
    assert finite==sorted(tuple(r['incident']) for r in C['finite_vertices'])
    for v in C['finite_vertices']:
        point=vec(v['point']);assert form(point,point,G).sign()<0 and form(point,H,G).sign()<0
        assert all(form(point,r,G).sign()==(0 if i in v['incident'] else -1) for i,r in enumerate(N))
    finite_edges=sum(len(set(a)&set(b))==8 for a,b in combinations(finite,2));assert finite_edges==10
    T,Ti=mat(S['matrix']),mat(S['inverse_matrix']);I=[[Q(i==j) for j in range(10)] for i in range(10)]
    assert mm(T,Ti)==mm(Ti,T)==I and mm(mm(tr(T),G),T)==G
    assert all(x.integral() for row in T+Ti for x in row)
    assert det(T)==-1
    trace=sum((T[i][i] for i in range(10)),Q());assert trace==Q(9,5) and (trace-10).sign()>0
    source,target=S['source_corner'],S['target_corner'];scale=dec(S['positive_perpendicular_scale'])
    assert scale==Q(18,14) and scale.sign()>0
    assert mv(T,vec(source['point']))==[scale*x for x in vec(target['point'])]
    for i,j in zip(source['ordered_indices'],target['ordered_indices']):assert mv(T,N[i])==N[j]
    for corner in (source,target):
        V=vec(corner['point']);ids=corner['ordered_indices'];assert positive(sub(GG,ids))
        assert form(V,V,G).sign()<0 and form(V,H,G).sign()<0
        assert all(form(V,n,G).sign()==(0 if j in ids else -1) for j,n in enumerate(N))
        D=corner['preservation_decomposition'];a=dec(D['vertex_coefficient']);b=vec(D['normal_coefficients'])
        assert a.sign()>0 and all(x.sign()<=0 for x in b)
        assert [a*V[k]+sum((b[i]*N[j][k] for i,j in enumerate(ids)),Q()) for k in range(10)]==H
    paths=[f'scripts/{TAG}.py']+[f'{folder}/{BASE}{suffix}.{ext}' for suffix in ('','_corners','_symmetry') for folder,ext in [('scripts','py'),('results','json')]]
    paths+=['results/fresh_frontier_construction_20261015.md',f'results/{BASE}.md']
    return dict(status='PASS_INDEPENDENT_NORM31_NONREFLECTIVITY_WITNESS_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        transitive_manifest_count=len(seen),walls=28,finite_candidate_supports=len(supports),
        finite_vertices=len(finite),finite_edges=finite_edges,trace=trace.enc(),determinant=['-1','0'],
        symmetry_order='infinite',scope='Independent literal field arithmetic, wall witnesses, complete finite corners, integral isometry/inverse, and both complete local sectors. Full Vinberg enumeration completeness is reviewed mathematically and replayed separately, not independently reimplemented.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
