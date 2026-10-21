#!/usr/bin/env python3
"""Independent exact checks for eight-simplex ideal count <=3 and charge12.

No project arithmetic or SymPy imports. The frozen census JSON supplies literal
edge lists; every required matrix/inertia/inverse is rebuilt over Q(sqrt2) using
Fraction. Classification completeness is explicitly an inherited proof input.
"""
import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite_edge_cusp_tradeoff12_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'

@dataclass(frozen=True)
class Q:
    a:F=F(0)
    b:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a));object.__setattr__(self,'b',F(self.b))
    def __add__(self,x):
        x=q(x);return Q(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,x):return self+-q(x)
    def __rsub__(self,x):return q(x)+-self
    def __mul__(self,x):
        x=q(x);return Q(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=q(x);d=x.a*x.a-2*x.b*x.b;assert d
        return self*Q(x.a/d,-x.b/d)
    def __rtruediv__(self,x):return q(x)/self
    def __bool__(self):return bool(self.a or self.b)
    def sign(self):
        if not self.a:return (self.b>0)-(self.b<0)
        if not self.b or self.a*self.b>0:return (self.a>0)-(self.a<0)
        d=self.a*self.a-2*self.b*self.b;assert d
        return ((self.a>0)-(self.a<0))*((d>0)-(d<0))
    def data(self):return [str(self.a),str(self.b)]
    def rational(self):assert not self.b;return self.a

def q(x):return x if isinstance(x,Q) else Q(x)
def matrix_data(A):return [[x.data() for x in r] for r in A]
def sub(A,ids):return [[A[i][j] for j in ids] for i in ids]
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def transpose(A):return [list(r) for r in zip(*A)]
def mm(A,B):return [[sum((x*y for x,y in zip(r,c)),Q()) for c in zip(*B)] for r in A]

def inverse(A):
    n=len(A);M=[list(r)+eye(n)[i] for i,r in enumerate(A)]
    for i in range(n):
        p=next(j for j in range(i,n) if M[j][i]);M[i],M[p]=M[p],M[i]
        d=M[i][i];M[i]=[x/d for x in M[i]]
        for j in range(n):
            if j!=i:
                c=M[j][i];M[j]=[x-c*y for x,y in zip(M[j],M[i])]
    J=[r[n:] for r in M];assert mm(A,J)==eye(n)
    return J

def inertia(A):
    A=[r[:] for r in A];p=m=z=0
    while A:
        n=len(A);i=next((i for i in range(n) if A[i][i]),None)
        if i is None:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if A[i][j]),None)
            if pair is None:z+=n;break
            i,j=pair;c=A[i][j];ids=[k for k in range(n) if k not in pair]
            A=[[A[a][b]-(A[a][i]*A[j][b]+A[a][j]*A[i][b])/c for b in ids] for a in ids]
            p+=1;m+=1
        else:
            d=A[i][i];p+=d.sign()>0;m+=d.sign()<0;ids=[j for j in range(n) if j!=i]
            A=[[A[a][b]-A[a][i]*A[i][b]/d for b in ids] for a in ids]
    return [p,m,z]

def kernel1(A):
    n=len(A);M=[r[:] for r in A];row=0;piv=[]
    for col in range(n):
        j=next((j for j in range(row,n) if M[j][col]),None)
        if j is None:continue
        M[row],M[j]=M[j],M[row];d=M[row][col];M[row]=[x/d for x in M[row]]
        for j in range(n):
            if j!=row:
                c=M[j][col];M[j]=[x-c*y for x,y in zip(M[j],M[row])]
        piv.append(col);row+=1
    free=[i for i in range(n) if i not in piv];assert len(free)==1
    v=[Q()]*n;v[free[0]]=Q(1)
    for i,j in enumerate(piv):v[j]=-M[i][free[0]]
    if all(x.sign()<0 for x in v):v=[-x for x in v]
    assert all(x.sign()>0 for x in v)
    assert mm(A,[[x] for x in v])==[[Q()] for _ in v]
    return v

def gram(n,edges):
    G=[[Q(2*(i==j)) for j in range(n)] for i in range(n)]
    for a,b,m in edges:
        assert m in (3,4);G[a][b]=G[b][a]=-Q(1) if m==3 else -Q(0,1)
    return G

def cases(diagrams):
    out=[];checked=[]
    for L in diagrams:
        n=L['order']
        if n<9:continue
        G=gram(n,L['edges']);assert inertia(G)==[n-1,1,0]
        aff=[];deletions=[]
        for j in range(n):
            ids=[x for x in range(n) if x!=j];H=sub(G,ids);I=inertia(H)
            assert I in ([n-1,0,0],[n-2,0,1])
            rec=dict(deleted=j,inertia=I)
            if I[-1]:
                aff.append(j);rec['positive_kernel']=[x.data() for x in kernel1(H)]
            deletions.append(rec)
        assert aff==L['ideal_deletions']
        checked.append(dict(diagram=L['index'],order=n,edges=L['edges'],
                            gram=matrix_data(G),deletions=deletions))
        for w in ([None] if n==9 else range(n)):
            black=[x for x in range(n) if x!=w];I=[x for x in aff if x in black]
            rec=dict(diagram=L['index'],internal_white=w,black=black,black_affine_deletions=I)
            if len(I)>=2:
                rec.update(method='external_whites_orthogonal_to_L',maximum_ideal_vertices=len(I))
            elif len(I)==1:
                a=I[0];thresholds=[]
                for j in black:
                    if j==a:continue
                    ids=[x for x in range(n) if x!=j];H=sub(G,ids);J=inverse(H)
                    t=(1/J[ids.index(a)][ids.index(a)]).rational();assert t>0
                    thresholds.append([j,str(t)])
                low=min(F(t) for _,t in thresholds);ties=[j for j,t in thresholds if F(t)==low]
                for j in black:
                    if j==a:continue
                    ids=[x for x in range(n) if x!=j];H=sub(G,ids)
                    H[ids.index(a)][ids.index(a)]-=Q(low)
                    assert inertia(H)==([n-2,0,1] if j in ties else [n-1,0,0])
                rec.update(method='rank_one_diagonal_decrease_at_black',decreased_black=a,
                           thresholds=thresholds,minimum=str(low),ties=ties,
                           maximum_ideal_vertices=1+len(ties))
            else:
                assert L['index']==25 and w==9
                rec.update(method='affine_E8_full_white_block',maximum_ideal_vertices=2)
            assert rec['maximum_ideal_vertices']<=3;out.append(rec)
    assert len(out)==34 and len(checked)==7
    assert Counter(x['method'] for x in out)==Counter(
        external_whites_orthogonal_to_L=19,rank_one_diagonal_decrease_at_black=14,
        affine_E8_full_white_block=1)
    return checked,out

def e8_case():
    # Independent literal affine E8, not assembled from the producer's affine model.
    edges=[[0,1,3],[0,2,3],[0,4,3],[2,3,3],[4,5,3],[5,6,3],[6,7,3],[7,8,3]]
    B=gram(9,edges);assert inertia(B)==[8,0,1]
    ker=kernel1(B);deletions=[]
    for j in range(9):
        ids=[i for i in range(9) if i!=j];H=sub(B,ids);assert inertia(H)==[8,0,0]
        J=inverse(H);assert all(x.sign()>=0 for r in J for x in r)
        deletions.append(dict(deleted=j,indices=ids,inverse=matrix_data(J)))
    ids8=[i for i in range(9) if i!=8];J8=inverse(sub(B,ids8))
    diag=[J8[i][i].rational() for i in range(8)]
    assert diag==list(map(F,[30,8,14,4,20,12,6,2]))
    ids3=[i for i in range(9) if i!=3];J3=inverse(sub(B,ids3))
    pair=[[J3[ids3.index(i)][ids3.index(j)] for j in (7,8)] for i in (7,8)]
    assert pair==[[Q(2),Q(1)],[Q(1),Q(1)]]
    profiles=[('A',8,Q(1)),('B',8,Q(0,1)),('C',7,Q(1))]
    cols=[]
    for name,node,value in profiles:
        cols.append([value if i==node else Q() for i in ids3])
    cost=[
        [sum((cols[i][a]*J3[a][b]*cols[j][b] for a in range(8) for b in range(8)),Q())
         for j in range(3)] for i in range(3)]
    assert cost==[[Q(1),Q(0,1),Q(1)],[Q(0,1),Q(2),Q(0,1)],[Q(1),Q(0,1),Q(2)]]
    # Nonnegative native contact makes residual offdiagonals no larger than
    # minus these strictly positive products. A-A is possible only at right angle.
    aa=[[Q(1),Q(-1)],[Q(-1),Q(1)]];aaa=[[Q(2*(i==j)-1) for j in range(3)] for i in range(3)]
    assert inertia(aa)==[1,0,1] and inertia(aaa)==[2,1,0]
    stars=[]
    for names in [['A'],['B'],['C'],['A','A']]:
        n=9+len(names);G=[[Q(2*(i==j)) for j in range(n)] for i in range(n)]
        for i in range(9):
            for j in range(9):G[i][j]=B[i][j]
        for k,name in enumerate(names):
            _,node,value=next(p for p in profiles if p[0]==name)
            G[9+k][node]=G[node][9+k]=-value
        assert inertia(G)==[n-1,1,0]
        ideals=[];vertices=[]
        for j in range(9):
            ids=[i for i in range(n) if i!=j];H=sub(G,ids);I=inertia(H)
            assert I in ([n-1,0,0],[n-2,0,1])
            if I[-1]:ideals.append(j)
            vertices.append(dict(deleted_black=j,inertia=I))
        expected={'A':[],'B':[3],'C':[3,8],'AA':[3]}[''.join(names)];assert ideals==expected
        rec=dict(direct_profiles=names,gram=matrix_data(G),vertex_inertias=vertices,
                 ideal_black_deletions=ideals)
        if names==['A']:
            thresholds=[]
            for j in range(9):
                ids=[i for i in range(n) if i!=j];H=sub(G,ids);J=inverse(H)
                t=(1/J[ids.index(9)][ids.index(9)]).rational()
                thresholds.append([j,str(t)])
                H[ids.index(9)][ids.index(9)]-=Q(1)
                assert inertia(H)==([n-2,0,1] if j==3 else [n-1,0,0])
            assert min(F(t) for _,t in thresholds)==1
            assert [j for j,t in thresholds if F(t)==1]==[3]
            rec.update(outside_diagonal_thresholds=thresholds,ideal_black_deletions_at_D1=[3])
        stars.append(rec)
    return dict(affine_gram=matrix_data(B),affine_positive_kernel=[x.data() for x in ker],
       black_deletion_inverses=deletions,delete8_inverse_diagonal=list(map(str,diag)),
       delete3_inverse_7_8=matrix_data(pair),profiles=[dict(name=a,node=b,magnitude=c.data()) for a,b,c in profiles],
       delete3_profile_cost_gram=matrix_data(cost),two_A_residual=matrix_data(aa),
       three_A_residual_inertia=inertia(aaa),all_direct_white_stars=stars,
       exact_profile_completeness='Deleteblack8 bounds every contact except7/8 tozero;7 has only magnitude1. Deleteblack3 rules out both and bounds8 bysqrt2. All nonzero Coxeter contacts have magnitude>=1, so the list A,B,C is complete without a label cutoff.',
       outside_white_reduction='Z is orthogonal toB; after eliminating PD G_Z, D=R^t G_Z^-1 R is PSD and entrywise nonnegative. SingleB/C residual0 forcesD0; twoA residual has positive null(1,1), forcingD0; singleA permits0<=D<=1 and only deletion3 becomesideal atD1.')

def run():
    source=json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    assert source['status']=='PASS_EXACT_HIGH_RANK_QUASILANNER_ATTACHMENT_CENSUS'
    assert source['orders']==[[6,12],[7,3],[8,4],[9,4],[10,3]]
    checked,rows=cases(source['diagrams']);special=e8_case();coeff=[]
    for vertices in range(2,11):
        for f in range(2,vertices+1):
            i=vertices-f
            if vertices==9 and i>3:continue
            if vertices==10 and (f,i) not in [(9,1),(8,2),(7,3)]:continue
            c=F(2)*(i+F(f-2,3));assert c<=12
            coeff.append(dict(dimension=vertices-1,f=f,i=i,coefficient=str(c)))
    assert [r for r in coeff if r['coefficient']=='12']==[dict(dimension=7,f=2,i=6,coefficient='12')]
    paths=[f'scripts/{TAG}.py',f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
           'scripts/fresh_finite_edge_cusp_tradeoff_20261017.py',
           'results/fresh_finite_edge_cusp_tradeoff_20261017.json',
           'results/fresh_finite_edge_cusp_tradeoff_20261017.md',
           'results/fresh_finite_ended_simplex_ports_20261017.md',
           'results/fresh_multi_ideal_triangle_scope_20261017.json',
           'results/fresh_multi_ideal_triangle_scope_20261017.md',
           'results/fresh_finite_volume_simplex_truncation_20261017.md']
    return dict(status='PASS_EXACT_EIGHT_SIMPLEX_AT_MOST_THREE_IDEAL_AND_FINITE_EDGE_COEFFICIENT12',
        scope='Complete internally reviewed geometric lemma, independent exact small arithmetic; no ambient dimension bound.',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        independent_arithmetic='Standalone Fraction Qsqrt2; no producer arithmetic imports; literal census graphs are input.',
        full_census_replayed_here=False,checked_high_rank_diagrams=checked,
        eight_simplex_cases=rows,eight_simplex_case_count=len(rows),exceptional_affine_E8=special,
        coefficient_cases=coeff,uniform_coefficient='12',unresolved_equality_type=dict(f=2,i=6,dimension=7),
        inequality='2(T0+T1) <= 12 E_ff = 6 n v_f - 6 E_fi <= 6 n v_f')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(out.read_text())==data
    else:out.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(data['status'],'saved_replay=',args.verify)
