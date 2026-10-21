#!/usr/bin/env python3
"""Exact compact Makarov H5 face inside a new H9 seed.

This is a geometric seed, not a finite-volume claim.  The bounded cap test
uses ten face vertices, four H4 weight rays, and labels 3,4,5,6,10.
Arithmetic is independently implemented here over Q(sqrt(5)).
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_h9_makarov_face_extension_20261013.json'

@dataclass(frozen=True)
class K:
    a:F=F(0)
    b:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a));object.__setattr__(self,'b',F(self.b))
    def __add__(self,o):
        o=kk(o);return K(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,o):return self+-kk(o)
    def __rsub__(self,o):return kk(o)+-self
    def __mul__(self,o):
        o=kk(o);return K(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=kk(o);d=o.a*o.a-5*o.b*o.b;assert d
        return self*K(o.a/d,-o.b/d)
    def __rtruediv__(self,o):return kk(o)/self
    def __bool__(self):return bool(self.a or self.b)
    def sign(self):
        if not self:return 0
        if not self.a:return (self.b>0)-(self.b<0)
        if not self.b:return (self.a>0)-(self.a<0)
        if self.a*self.b>0:return (self.a>0)-(self.a<0)
        d=self.a*self.a-5*self.b*self.b
        return ((self.a>0)-(self.a<0))*((d>0)-(d<0))
    def enc(self):return [str(self.a),str(self.b)]
    def approx(self):return float(self.a)+float(self.b)*5**.5

def kk(x):return x if isinstance(x,K) else K(x)
Z=K();O=K(1);PHI=K(F(1,2),F(1,2))
def dot(x,y):return sum((a*b for a,b in zip(x,y)),Z)
def mv(A,x):return [dot(a,x) for a in A]
def pair(x,y,Q):return dot(x,mv(Q,y))
def ident(n):return [[K(int(i==j)) for j in range(n)]for i in range(n)]
def inv(A):
    n=len(A);M=[list(r)+s for r,s in zip(A,ident(n))]
    for j in range(n):
        k=next(i for i in range(j,n) if M[i][j]);M[j],M[k]=M[k],M[j]
        t=M[j][j];M[j]=[x/t for x in M[j]]
        for i in range(n):
            if i!=j:
                t=M[i][j];M[i]=[x-t*y for x,y in zip(M[i],M[j])]
    return [r[n:]for r in M]
def null1(A):
    M=[list(r)for r in A];n=len(M[0]);piv=[];k=0
    for j in range(n):
        found=next((i for i in range(k,len(M))if M[i][j]),None)
        if found is None:continue
        M[k],M[found]=M[found],M[k];t=M[k][j];M[k]=[x/t for x in M[k]]
        for i in range(len(M)):
            if i!=k:
                t=M[i][j];M[i]=[x-t*y for x,y in zip(M[i],M[k])]
        piv.append(j);k+=1
    if len(piv)!=n-1:return None
    free=next(j for j in range(n)if j not in piv);v=[Z]*n;v[free]=O
    for i,j in enumerate(piv):v[j]=-M[i][free]
    assert all(not dot(r,v)for r in A)
    return v
def inertia(A):
    A=[list(r)for r in A];p=n=0
    while A:
        assert A[0][0]
        s=A[0][0].sign();p+=s>0;n+=s<0
        A=[[A[i][j]-A[i][0]*A[0][j]/A[0][0]for j in range(1,len(A))]for i in range(1,len(A))]
    return [p,n]
def matrix_enc(A):return [[x.enc()for x in row]for row in A]

RAD=K(F(336,11),F(128,11))
RAD_SQRT11_MULTIPLIER=K(F(16,11),F(4,11))
assert 11*RAD_SQRT11_MULTIPLIER*RAD_SQRT11_MULTIPLIER==RAD
assert RAD_SQRT11_MULTIPLIER.sign()>0
@dataclass(frozen=True)
class L:
    a:K=Z
    b:K=Z
    def __post_init__(self):
        object.__setattr__(self,'a',kk(self.a));object.__setattr__(self,'b',kk(self.b))
    def __add__(self,o):
        o=ll(o);return L(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return L(-self.a,-self.b)
    def __sub__(self,o):return self+-ll(o)
    def __rsub__(self,o):return ll(o)+-self
    def __mul__(self,o):
        o=ll(o);return L(self.a*o.a+RAD*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=ll(o);d=o.a*o.a-RAD*o.b*o.b;assert d
        return self*L(o.a/d,-o.b/d)
    def sign(self):
        a,b=self.a.sign(),self.b.sign()
        if not b:return a
        if not a:return b
        if a==b:return a
        return a*(self.a*self.a-RAD*self.b*self.b).sign()
    def enc(self):return [self.a.enc(),self.b.enc()]
def ll(x):return x if isinstance(x,L)else L(kk(x))
def lp(x,y,Q):
    return sum((ll(x[i])*Q[i][j]*y[j]for i in range(len(x))for j in range(len(y))),L())

def build():
    G=ident(6)
    for i in range(6):G[i][i]=K(2)
    for i,a in enumerate([PHI,O,O,O,O]):G[i][i+1]=G[i+1][i]=-a
    I=inv(G);assert inertia(G)==[5,1]
    W=ident(6)+[[-I[i][5]for i in range(6)]]
    norms=[pair(w,w,G)for w in W]
    assert norms[-1]==K(F(7,11),F(-1,11))
    assert all(pair(W[-1],W[i],G)==K(-int(i==5))for i in range(6))
    # The sole dotted edge is the original last wall and its polar.
    assert (1/(2*norms[-1])-1).sign()>0
    V=[];active=[];ray_count=0
    for inds in combinations(range(7),5):
        v=null1([mv(G,W[i])for i in inds])
        if v is None:continue
        ps=[pair(v,w,G)for w in W]
        if all(x.sign()>=0 for x in ps):v=[-x for x in v];ps=[-x for x in ps]
        if any(x.sign()>0 for x in ps):continue
        ray_count+=1
        assert pair(v,v,G).sign()<0
        assert [i for i,x in enumerate(ps)if not x]==list(inds)
        assert inertia([[pair(W[i],W[j],G)for j in inds]for i in inds])==[5,0]
        V.append(v);active.append(list(inds))
    assert ray_count==len(V)==10
    assert all(pair(v,w,G).sign()<0 for v in V for w in V)
    h=[sum((v[i]for v in V),Z)for i in range(6)]
    assert pair(h,h,G).sign()<0 and all(pair(h,w,G).sign()<0 for w in W)
    # Positive dual certificate: h is in the interior of cone(W). It makes
    # the h=-1 slice bounded, so the exhaustive extreme-ray check above
    # really proves compactness, without borrowing the prism classification.
    # Choose a positive coefficient for W6 and solve the remaining basis ones.
    candidates=[K(F(i,100))for i in range(1,10001)]
    cert=None
    for t in candidates:
        q=[h[i]-t*W[6][i]for i in range(6)]
        if all(x.sign()>0 for x in q):cert=q+[t];break
    assert cert is not None
    assert all(sum((cert[j]*W[j][i]for j in range(7)),Z)==h[i]for i in range(6))
    H=ident(4)
    for i in range(4):H[i][i]=K(2)
    for i,a in enumerate([PHI,O,O]):H[i][i+1]=H[i+1][i]=-a
    J=inv(H);assert inertia(H)==[4,0]
    Q=[r+[Z]*4 for r in G]+[[Z]*6+r for r in H]
    assert inertia(Q)==[9,1]
    # This exact inward perturbation makes every seed inequality strict.
    c=h+[-sum(J[i],Z)/100 for i in range(4)]
    roots=[w+[Z]*4 for w in W]+[[Z]*6+e for e in ident(4)]
    assert pair(c,c,Q).sign()<0 and all(pair(c,w,Q).sign()<0 for w in roots)
    lifted_vertices=[v+[Z]*4 for v in V]
    for v,inds in zip(lifted_vertices,active):
        a=inds+list(range(7,11))
        assert [i for i,w in enumerate(roots)if not pair(v,w,Q)]==a
        assert all(pair(v,w,Q).sign()<=0 for w in roots)
        assert inertia([[pair(roots[i],roots[j],Q)for j in a]for i in a])==[9,0]
    # Cap r=(sqrt(k)*v,-sqrt(X)*H^-1 e_j). It has norm two,
    # angle pi/m to H4 wall j and right angles to the other three.
    labels={3:O,4:K(2),5:PHI+1,6:K(3),10:PHI+2}
    # Every tested squared cosine lies in Q(sqrt(5)).  Degree of cos(2pi/n)
    # is phi(n)/2, so only n=2,3,4,5,6,8,10,12 are possible; n=8,12
    # generate Q(sqrt(2)), Q(sqrt(3)), and hence cannot occur here.
    oldangles={2:Z,3:K(F(1,4)),4:K(F(1,2)),5:(PHI+1)/4,6:K(F(3,4)),10:(PHI+2)/4}
    caps=[];rejections=[]
    for vi,v in enumerate(V):
        for j in range(4):
            for m,X in labels.items():
                k=(X*J[j][j]-2)/(-pair(v,v,G));assert k.sign()>0
                vals=[k*pair(v,w,G)*pair(v,w,G)/(2*q)for w,q in zip(W,norms)]
                tags=[];bad=[]
                for i,y in enumerate(vals):
                    if y==O:tag='parallel'
                    elif (y-1).sign()>0:tag='dotted'
                    else:
                        lab=next((n for n,z in oldangles.items()if y==z),None)
                        tag='m'+str(lab)if lab else 'forbidden'
                    tags.append(tag)
                    if tag=='forbidden':bad.append(i)
                row={'vertex':vi,'h4_index':j,'label':m,'X':X.enc(),'old_scale_square':k.enc(),'old_pair_types':tags,'old_cosine_squares':[y.enc()for y in vals]}
                (rejections if bad else caps).append(row)
    gaps={f'{i},{j}':(J[i][j]*J[i][j]-(J[i][i]-F(1,2))*(J[j][j]-F(1,2))).enc()for i,j in [(0,1),(1,2),(2,3)]}
    assert all((J[i][j]*J[i][j]-(J[i][i]-F(1,2))*(J[j][j]-F(1,2))).sign()>0 for i,j in [(0,1),(1,2),(2,3)])
    # Exact strict escape for ALL forty vertex-centered weight-ray families,
    # even if every incompatible family is imposed simultaneously.  At m=3
    # each pairing is sqrt(k)*A-E_j. Here A<0,E_j<0, so one squared
    # inequality certifies strict negativity. For all m>=3, dividing the
    # cap by 2cos(pi/m) increases the coefficient of the negative A.
    E=[K(i)for i in [152191883198,192187169028,138122715479,84058286892,29993903461,12504740634,-257820105850,-312871271018,-208580847566,-104290423893]]
    assert pair(E,E,Q).sign()>0 and pair(E,c,Q).sign()<0
    assert all(pair(E,w,Q).sign()<0 for w in roots)
    proof=[]
    for vi,v in enumerate(V):
        a=pair(v,E[:6],G);assert a.sign()<0
        for j in range(4):
            k=(J[j][j]-2)/(-pair(v,v,G))
            gap=k*a*a-E[j+6]*E[j+6]
            assert E[j+6].sign()<0 and gap.sign()>0
            proof.append({'vertex':vi,'h4_index':j,'old_pair':a.enc(),'squared_strict_gap':gap.enc()})
    # Changed mechanism: release the vertex-center condition, using the
    # compact face's A3 face 1,2,3. Impose (u,w0)=-1,(u,w5)=0 and the
    # three A3 zero pairings; q=(u,w4) is determined by the norm equation.
    aa=I[4][4];bb=-2*I[0][4];cc=I[0][0]-2+J[0][0]
    assert bb*bb-4*aa*cc==RAD and RAD.sign()>0
    q=(L(-bb,O))/(2*aa)
    u=[L(-I[i][0])+q*I[i][4]for i in range(6)]
    r=u+[L(-J[i][0])for i in range(4)]
    assert lp(r,r,Q)==L(K(2)) and lp(r,c,Q).sign()<0
    types=[]
    for i,w in enumerate(roots):
        p=lp(r,w,Q);assert p.sign()<=0
        if i in [0,7]:assert p==L(K(-1));types.append('m3')
        elif i in [4,6]:assert (p*p-2*pair(w,w,Q)).sign()>0;types.append('dotted')
        else:assert p==L();types.append('m2')
    assert lp(r,E,Q).sign()>0
    assert all(lp(r,v,Q).sign()<0 for v in lifted_vertices)
    # Actual new facet: project the strict center to its supporting plane.
    cp=lp(c,r,Q);f=[L(x)-cp*y/2 for x,y in zip(c,r)]
    assert lp(f,f,Q).sign()<0 and lp(f,r,Q)==L()
    assert all(lp(f,w,Q).sign()<0 for w in roots)
    # The twelve-wall chamber still has infinite volume.
    EE=[K(i)for i in [209087944259,258446925556,178582467912,98718035229,18853647703,110402,-264878711768,-327408093401,-226233428045,-125058762468]]
    assert pair(EE,EE,Q).sign()>0 and pair(EE,c,Q).sign()<0
    assert all(pair(EE,w,Q).sign()<0 for w in roots) and lp(r,EE,Q).sign()<0
    return {'status':'EXACT_COMPACT_H5_FACE_H9_SEED_AND_CHANGED_CENTER_CUTTER','base_field':'a+b*sqrt(5); used by the form, seed normals and base sign certificates','changed_center_field':'Q(sqrt(5),sqrt(11)); encoded as A+B*sqrt(RAD) with A,B in the base field','positive_sqrt_RAD_equals_sqrt11_times':RAD_SQRT11_MULTIPLIER.enc(),'form':matrix_enc(Q),'seed_normals':matrix_enc(roots),'strict_center':[x.enc()for x in c],'strict_center_norm':pair(c,c,Q).enc(),'compact_H5_vertices':matrix_enc(V),'vertex_active_sets':active,'positive_dual_center_coefficients':[x.enc()for x in cert],'H4_inverse':matrix_enc(J),'same_center_adjacent_gaps':gaps,'tested_cap_labels':list(labels),'tested_cap_count':200,'accepted_caps':caps,'rejected_caps':rejections,
            'all_label_vertex_center_family_obstruction':{'strict_integer_escape':[int(x.a)for x in E],'norm':pair(E,E,Q).enc(),'center_pair':pair(E,c,Q).enc(),'forty_strict_m3_checks':proof,'scope':'Every real transverse pairing magnitude t>=1, hence all finite labels m>=3, in all forty vertex-center/fundamental-weight families. Arbitrary finite additions from these families preserve this strict escape.'},
            'changed_center_cap':{'extension_radicand':RAD.enc(),'normal':[x.enc()for x in r],'norm':lp(r,r,Q).enc(),'q_old4':q.enc(),'old_pair_types':types,'center_pair':lp(r,c,Q).enc(),'old_escape_cut':lp(r,E,Q).enc(),'actual_facet_point':[x.enc()for x in f],'retains_ten_old_vertices':True,'new_strict_integer_escape':[int(x.a)for x in EE],'new_escape_norm':pair(EE,EE,Q).enc(),'new_escape_cap_pair':lp(r,EE,Q).enc()},
            'scope':'The 11-wall H9 seed has a compact H5 face and ten actual finite H9 vertices. The 200-cap screen and all-label forty-family obstruction are scoped. A new exact cap centered away from those vertices cuts their common escape and gives a coherent 12-wall chamber, still infinite volume. No unrestricted H9 completion theorem.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=build()
    if args.verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS',len(out['accepted_caps']),'accepted caps;',len(out['rejected_caps']),'rejections; compact face/10 finite vertices')
if __name__=='__main__':main()
