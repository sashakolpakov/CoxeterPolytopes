#!/usr/bin/env python3
"""Exact local arithmetic D4+2E8 edge with D5+2E8 endpoints.

Standalone Fraction arithmetic in Q(phi), phi^2=phi+1. The displayed chamber
has an explicit strict null escape and infinite volume; no lattice or new
finite-volume Coxeter polytope is claimed.
"""
import argparse,json
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_d4_two_e8_self_edge_20261016'

@dataclass(frozen=True)
class Q:
    a:F=F(0)
    b:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a));object.__setattr__(self,'b',F(self.b))
    @staticmethod
    def of(x):return x if isinstance(x,Q) else Q(x)
    def __add__(self,y):
        y=Q.of(y);return Q(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,y):return self+-Q.of(y)
    def __rsub__(self,y):return Q.of(y)+-self
    def __mul__(self,y):
        y=Q.of(y)
        return Q(self.a*y.a+self.b*y.b,self.a*y.b+self.b*y.a+self.b*y.b)
    __rmul__=__mul__
    def conj(self):return Q(self.a+self.b,-self.b)
    def inv(self):
        n=self.a*self.a+self.a*self.b-self.b*self.b;assert n
        return Q((self.a+self.b)/n,-self.b/n)
    def __truediv__(self,y):return self*Q.of(y).inv()
    def __rtruediv__(self,y):return Q.of(y)*self.inv()
    def sign(self):
        A=2*self.a+self.b;B=self.b
        sg=lambda x:(x>0)-(x<0)
        if not A:return sg(B)
        if not B:return sg(A)
        if sg(A)==sg(B):return sg(A)
        return sg(A)*sg(A*A-5*B*B)
    def __bool__(self):return bool(self.a or self.b)
    def saved(self):return [str(self.a),str(self.b)]

phi=Q(0,1)
def dot(x,y):return sum((a*b for a,b in zip(x,y)),Q())
def mv(A,v):return [dot(row,v) for row in A]
def sub(A,ids):return [[A[i][j] for j in ids] for i in ids]
def inverse(A):
    n=len(A);M=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        p=next(i for i in range(k,n) if M[i][k]);M[p],M[k]=M[k],M[p]
        z=M[k][k];M[k]=[v/z for v in M[k]]
        for i in range(n):
            if i!=k:
                z=M[i][k];M[i]=[a-z*b for a,b in zip(M[i],M[k])]
    J=[row[n:] for row in M]
    assert all(dot(A[i],[J[k][j] for k in range(n)])==Q(i==j) for i in range(n) for j in range(n))
    return J

def pivots(A):
    A=[list(row) for row in A];out=[]
    while A:
        d=A[0][0];assert d;out.append(d)
        A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,len(A))] for i in range(1,len(A))]
    return out

def schur(A,eliminate,keep):
    J=inverse(sub(A,eliminate))
    return [[A[i][j]-dot([A[i][k] for k in eliminate],mv(J,[A[k][j] for k in eliminate])) for j in keep] for i in keep]

def graph_gram(n,edges):
    A=[[Q(2*(i==j)) for j in range(n)] for i in range(n)]
    for i,j,v in edges:A[i][j]=A[j][i]=Q.of(v)
    return A

def matrices():
    core=graph_gram(6,[(0,1,-1),(0,2,-1),(0,3,-1),(1,4,-1),(1,5,-1),(4,5,-phi)])
    # Finite E8 is the simply-laced star with arm lengths1,2,4.
    E=graph_gram(8,[(0,1,-1),(0,2,-1),(2,3,-1),(0,4,-1),(4,5,-1),(5,6,-1),(6,7,-1)])
    A=[[Q() for _ in range(22)] for _ in range(22)]
    for off,M in ((0,core),(6,E),(14,E)):
        for i,row in enumerate(M):
            for j,v in enumerate(row):A[off+i][off+j]=v
    return core,E,A

def save_matrix(A):return [[v.saved() for v in row] for row in A]
def save_vector(v):return [z.saved() for z in v]

def run():
    assert phi*phi==phi+1 and phi.sign()>0 and phi.conj().sign()<0
    core,E,A=matrices();J=inverse(core);EJ=inverse(E)
    assert all(v.a.denominator==v.b.denominator==1 for row in A for v in row)
    assert all(A[i][j].sign()<=0 for i in range(22) for j in range(i))
    assert all(A[i][j] in (Q(),Q(-1),-phi) for i in range(22) for j in range(i))
    D=sub(core,list(range(4)));dp=pivots(D);ep=pivots(E)
    assert all(v.sign()>0 for v in dp+ep)
    partial=schur(core,[0,2,3],[1,4,5])
    assert partial==[[Q(1),Q(-1),Q(-1)],[Q(-1),Q(2),-phi],[Q(-1),-phi,Q(2)]]
    black=schur(core,list(range(4)),[4,5])
    assert black==[[Q(1),-1-phi],[-1-phi,Q(1)]]
    bp=pivots(black);assert bp==[Q(1),-1-3*phi]
    assert [v.sign() for v in bp]==[1,-1]
    conjugate=[[v.conj() for v in row] for row in A]
    cp=pivots(conjugate);assert all(v.sign()>0 for v in cp)
    ap=pivots(A);assert CounterSigns(ap)==dict(positive=21,negative=1,zero=0)
    # r=4,s=5 are the two framing walls. v_r omits r and v_s omits s.
    vertices=[];rays=[]
    for omitted in (4,5):
        v=[-J[i][omitted] for i in range(6)]+[Q()]*16
        Av=mv(A,v);square=dot(v,Av)
        assert Av==[-Q(i==omitted) for i in range(22)] and square.sign()<0
        active=[i for i in range(22) if i!=omitted]
        vp=pivots(sub(A,active));assert all(z.sign()>0 for z in vp)
        core_ids=[i for i in range(6) if i!=omitted]
        degrees=sorted(sum(bool(core[i][j]) for j in core_ids if j!=i) for i in core_ids)
        assert degrees==[1,1,1,2,3]
        assert all(core[i][j] in (Q(),Q(-1)) for i in core_ids for j in core_ids if i!=j)
        assert square==Q(F(4,5),F(-3,5))
        rays.append(v);vertices.append(dict(omitted_framing_wall=omitted,
            finite_vertex_type=['D5','E8','E8'],coefficients=save_vector(v),
            wall_pairings=save_vector(Av),square=square.saved(),active_gram_positive_pivots=save_vector(vp)))
    xy=dot(rays[0],mv(A,rays[1]));assert xy.sign()<0
    cosh=xy/Q(F(4,5),F(-3,5));assert cosh==1+phi
    common=list(range(4))+list(range(6,22));assert len(common)==20
    assert all(v.sign()>0 for v in pivots(sub(A,common)))
    # The positive segment between the two timelike rays has exactly these
    # common walls active and the two framing inequalities strict inside.
    midpoint=[a+b for a,b in zip(*rays)];Am=mv(A,midpoint)
    assert dot(midpoint,Am).sign()<0 and [Am[i] for i in (4,5)]==[Q(-1),Q(-1)]
    assert all(not Am[i] for i in common)
    # Strict timelike core vector; append an orthogonal negative Weyl vector
    # in each E8 block. No nonzero K-rational null ray can exist because the
    # conjugate form is positive definite, so use a declared real radical t.
    eps=F(1,10);p=[Q(eps)]*4+[Q(1),Q(1)]
    u=[-z for z in mv(J,p)]+[Q()]*16;Au=mv(A,u);U=dot(u,Au)
    assert U.sign()<0 and Au[:6]==[-z for z in p] and not any(Au[6:])
    eh=[-z for z in mv(EJ,[Q(1)]*8)]
    h=[Q()]*6+eh+eh;Ah=mv(A,h);H=dot(h,Ah)
    assert H==Q(1240) and Ah[:6]==[Q()]*6 and Ah[6:]==[Q(-1)]*16
    assert not dot(u,Ah)
    kappa=-U/H;assert kappa.sign()>0
    # z=u+t*h, t=sqrt(kappa)>0. Verify the norm and causal pairing by exact
    # polynomial coefficients modulo t^2=kappa; every wall pairing is strict.
    norm_constant=U+kappa*H;norm_t=2*dot(u,Ah)
    assert not norm_constant and not norm_t
    assert dot(u,Au).sign()<0 and not dot(h,Au)
    interior_scale=F(1,1000)
    w=[a+interior_scale*b for a,b in zip(u,h)];Aw=mv(A,w)
    assert all(v.sign()<0 for v in Aw) and dot(w,Aw).sign()<0 and dot(w,Au).sign()<0
    # The cycle leaf1-r-s has cyclic Gram product -phi, so the field really
    # contains Q(phi), rather than merely admitting a presentation over it.
    cycle=core[1][4]*core[4][5]*core[5][1];assert cycle==-phi
    # Reflection matrices in the simple-root basis have coefficients in Z[phi]:
    # s_i(e_j)=e_j-G_ij e_i, since G_ii=2.
    assert all(A[i][j].a.denominator==A[i][j].b.denominator==1 for i in range(22) for j in range(22))
    path=f'scripts/{TAG}.py'
    return dict(status='PASS_EXACT_ADMISSIBLE_D4_TWO_E8_SELF_EDGE_WITH_STRICT_NULL_ESCAPE',
        sources={path:sha256((ROOT/path).read_bytes()).hexdigest()},
        field=dict(name='Q(sqrt(5))',basis=['1','phi'],relation='phi^2=phi+1',
            real_embedding='phi=(1+sqrt(5))/2',conjugation='phi -> 1-phi',
            cyclic_product_generating_field=cycle.saved()),
        root_gram=save_matrix(A),root_count=22,ambient_hyperbolic_dimension=21,
        framing_walls=[4,5],common_edge_walls=common,common_edge_type=['D4','E8','E8'],
        gram_signature=dict(positive=21,negative=1,zero=0),physical_ldl_pivots=save_vector(ap),
        conjugate_signature=dict(positive=22,negative=0,zero=0),conjugate_positive_pivots=save_vector(cp),
        partial_schur_on_leaf1_r_s=save_matrix(partial),black_schur=save_matrix(black),
        integral_reflection_matrices=True,nonright_coxeter_labels=[3,5],
        endpoints=vertices,endpoint_pairing=xy.saved(),edge_length_cosh=cosh.saved(),
        edge_midpoint=save_vector(midpoint),edge_midpoint_square=dot(midpoint,Am).saved(),
        strict_interior_timelike_point=dict(coefficients=save_vector(w),square=dot(w,Aw).saved(),wall_pairings=save_vector(Aw)),
        strict_null_escape=dict(definition='z=u+t*h with t=sqrt(kappa)>0',
            u=save_vector(u),h=save_vector(h),kappa=kappa.saved(),u_square=U.saved(),h_square=H.saved(),
            norm_constant=norm_constant.saved(),norm_t_coefficient=norm_t.saved(),
            wall_pairing_coefficients=[[a.saved(),b.saved()] for a,b in zip(Au,Ah)],
            causal_pairing_with_u=U.saved(),all_wall_pairings_strictly_negative=True,
            scalar_extension='Real radical extension; no nonzero null vector over Q(phi) is asserted.'),
        conclusion='The displayed chamber contains a genuine finite edge with two D5+2E8 vertices, has an interior timelike point, and has an open ideal-boundary patch from the strict null ray; hence its volume is infinite.',
        scope='Local endpoint compatibility in an integral first-type admissible quadratic form only. The integral reflection subgroup is not claimed to be an arithmetic lattice. No finite-volume polytope, completion by further walls, new-dimensional example, or improved dimension bound is claimed.')

def CounterSigns(vs):
    return dict(positive=sum(v.sign()>0 for v in vs),negative=sum(v.sign()<0 for v in vs),zero=sum(v.sign()==0 for v in vs))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print('endpoint square',out['endpoints'][0]['square'],'u square',out['strict_null_escape']['u_square'],'kappa',out['strict_null_escape']['kappa'])
