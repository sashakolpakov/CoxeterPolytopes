#!/usr/bin/env python3
"""Mixed compact H5 prism + A4: exact seed and bounded global cap-family test.

This does NOT compactify a noncompact H9 chamber.  It materializes the
literature-led fallback after auditing the prior T1/T2 compactification lane.
No old construction/source is edited. Arithmetic K is imported with a bound
source hash; the new quadratic extension and all geometry are checked here.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import argparse,hashlib,json,math
from fresh_h9_makarov_face_extension_20261013 import K,kk,PHI

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_h9_compactification_20261013.json'
DEP=ROOT/'scripts/fresh_h9_makarov_face_extension_20261013.py'
R=K(F(16,11),F(4,11))
@dataclass(frozen=True)
class L:
    a:K=K()
    b:K=K()
    def __post_init__(self):
        object.__setattr__(self,'a',kk(self.a));object.__setattr__(self,'b',kk(self.b))
    def __add__(self,o):
        o=ll(o);return L(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return L(-self.a,-self.b)
    def __sub__(self,o):return self+-ll(o)
    def __rsub__(self,o):return ll(o)+-self
    def __mul__(self,o):
        o=ll(o);return L(self.a*o.a+R*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=ll(o);d=o.a*o.a-R*o.b*o.b;assert d
        return self*L(o.a/d,-o.b/d)
    def __rtruediv__(self,o):return ll(o)/self
    def __bool__(self):return bool(self.a or self.b)
    def sign(self):
        a,b=self.a.sign(),self.b.sign()
        if not b:return a
        if not a:return b
        if a==b:return a
        return a*(self.a*self.a-R*self.b*self.b).sign()
    def enc(self):return [self.a.enc(),self.b.enc()]
    def approx(self):return self.a.approx()+self.b.approx()*math.sqrt(R.approx())
def ll(x):return x if isinstance(x,L)else L(kk(x))
Z=L();O=L(1)
def sm(xs):return sum(xs,Z)
def dot(a,b):return sm(x*y for x,y in zip(a,b))
def mv(A,x):return [dot(a,x)for a in A]
def pair(a,b,Q):return dot(a,mv(Q,b))
def eye(n):return [[L(int(i==j))for j in range(n)]for i in range(n)]
def inv(A):
    n=len(A);M=[list(r)+s for r,s in zip(A,eye(n))]
    for j in range(n):
        k=next(i for i in range(j,n)if M[i][j]);M[j],M[k]=M[k],M[j]
        d=M[j][j];M[j]=[x/d for x in M[j]]
        for i in range(n):
            if i!=j:
                d=M[i][j];M[i]=[x-d*y for x,y in zip(M[i],M[j])]
    return [r[n:]for r in M]
def null1(A):
    n=len(A[0]);M=[list(r)for r in A];piv=[];k=0
    for j in range(n):
        z=next((i for i in range(k,len(M))if M[i][j]),None)
        if z is None:continue
        M[k],M[z]=M[z],M[k];d=M[k][j];M[k]=[x/d for x in M[k]]
        for i in range(len(M)):
            if i!=k:
                d=M[i][j];M[i]=[x-d*y for x,y in zip(M[i],M[k])]
        piv.append(j);k+=1
    if len(piv)!=n-1:return None
    free=next(j for j in range(n)if j not in piv);v=[Z]*n;v[free]=O
    for i,j in enumerate(piv):v[j]=-M[i][free]
    assert all(not dot(a,v)for a in A)
    return v
def inertia(A):
    A=[list(r)for r in A];p=n=0
    while A:
        assert A[0][0];s=A[0][0].sign();p+=s>0;n+=s<0
        A=[[A[i][j]-A[i][0]*A[0][j]/A[0][0]for j in range(1,len(A))]for i in range(1,len(A))]
    return [p,n]
def enc(v):return [x.enc()for x in v]
def mat(A):return [enc(r)for r in A]
def sqrat(q):
    return math.isqrt(q.numerator)**2==q.numerator and math.isqrt(q.denominator)**2==q.denominator
def geometry():
    G=eye(5)
    for i in range(5):G[i][i]=L(2)
    for i,a in enumerate([PHI,K(1),K(1),K(1)]):G[i][i+1]=G[i+1][i]=-L(a)
    I=inv(G);assert inertia(G)==[4,1]
    d=I[4][4];D0=2-d;assert D0==L(K(F(7,4),F(1,4)))
    assert (4-4*d)/D0==L(R)
    Q5=[r+[Z]for r in G]+[[Z]*5+[D0]]
    g=[-I[i][4]for i in range(5)]
    W=[e+[Z]for e in eye(5)]+[g+[O],[2*x for x in g]+[L(K(),K(-1))]]
    assert inertia(Q5)==[5,1]
    norms=[pair(w,w,Q5)for w in W];assert norms==[L(2)]*6+[L(4)]
    assert pair(W[5],W[6],Q5).sign()<0
    assert (pair(W[5],W[6],Q5)*pair(W[5],W[6],Q5)-8).sign()>0
    V=[];inds=[]
    for ids in combinations(range(7),5):
        v=null1([mv(Q5,W[i])for i in ids])
        if v is None:continue
        ps=[pair(v,w,Q5)for w in W]
        if all(x.sign()>=0 for x in ps):v=[-x for x in v];ps=[-x for x in ps]
        if any(x.sign()>0 for x in ps):continue
        assert pair(v,v,Q5).sign()<0
        assert [i for i,x in enumerate(ps)if not x]==list(ids)
        assert inertia([[pair(W[i],W[j],Q5)for j in ids]for i in ids])==[5,0]
        V.append(v);inds.append(list(ids))
    assert len(V)==10 and all(pair(v,w,Q5).sign()<0 for v in V for w in V)
    h=[-sm(row)for row in I]+[Z]
    assert pair(h,h,Q5).sign()<0 and all(pair(h,w,Q5).sign()<0 for w in W)
    # Positive dual-center certificate => bounded projective h=-1 slice.
    eps=L(F(1,100));a5=L(K(),K(1))*eps
    cert=[h[i]-(a5+2*eps)*g[i]for i in range(5)]+[a5,eps]
    assert all(x.sign()>0 for x in cert)
    assert all(sm(cert[j]*W[j][i]for j in range(7))==h[i]for i in range(6))
    H=eye(4)
    for i in range(4):H[i][i]=L(2)
    for i in range(3):H[i][i+1]=H[i+1][i]=L(-1)
    J=inv(H);assert inertia(H)==[4,0]
    Q=[r+[Z]*4 for r in Q5]+[[Z]*6+r for r in H]
    roots=[w+[Z]*4 for w in W]+[[Z]*6+e for e in eye(4)]
    c=h+[-sm(r)/100 for r in J]
    assert inertia(Q)==[9,1] and pair(c,c,Q).sign()<0
    assert all(pair(c,w,Q).sign()<0 for w in roots)
    for v,ids in zip(V,inds):
        vv=v+[Z]*4;active=ids+list(range(7,11))
        ps=[pair(vv,w,Q)for w in roots]
        assert [i for i,z in enumerate(ps)if not z]==active
        assert all(z.sign()<=0 for z in ps) and pair(vv,vv,Q).sign()<0
        assert inertia([[pair(roots[i],roots[j],Q)for j in active]for i in active])==[9,0]
    # Every seed inequality is a genuine facet by exact strict projection.
    for w in roots:
        f=[x-pair(c,w,Q)/pair(w,w,Q)*y for x,y in zip(c,w)]
        assert pair(f,f,Q).sign()<0 and not pair(f,w,Q)
        assert all(pair(f,v,Q).sign()<0 for v in roots if v!=w)
    return G,Q5,Q,W,V,inds,h,cert,H,J,roots,c

def build():
    G,Q5,Q,W,V,inds,h,cert,H,J,roots,c=geometry()
    # L is nonnormal quartic: N(R)=16/11, neither N nor N/5 is rational square.
    N=R.a*R.a-5*R.b*R.b;assert N==F(16,11)
    assert not sqrat(N) and not sqrat(N/5)
    # Thus the only cyclotomic subfield of L is K; squared finite angles
    # in L have labels2,3,4,5,6,10. This gate is complete for these profiles.
    ph=L(PHI)
    angles={2:Z,3:L(F(1,4)),4:L(F(1,2)),5:(ph+1)/4,6:L(F(3,4)),10:(ph+2)/4}
    vals={2:Z,3:O,5:ph}
    seed_types=[]
    for i,j in combinations(range(len(roots)),2):
        p=pair(roots[i],roots[j],Q);assert p.sign()<=0
        y=p*p/(pair(roots[i],roots[i],Q)*pair(roots[j],roots[j],Q))
        if (y-1).sign()>=0:tag='parallel'if y==O else'dotted'
        else:
            m=next((m for m,a in angles.items()if y==a),None);assert m is not None
            tag='m'+str(m)
        seed_types.append([i,j,tag])
    caps=[];bad=[];nonsp=[]
    for vi,v in enumerate(V):
        for labels in product(vals,repeat=4):
            b=[vals[m]for m in labels];w=mv(J,b);s=dot(b,w)-2
            if s.sign()<=0:nonsp.append([vi,list(labels)]);continue
            k=s/(-pair(v,v,Q5));types=[];forbidden=[]
            for i,r in enumerate(W):
                p=pair(v,r,Q5);assert p.sign()<=0
                cs=k*p*p/(2*pair(r,r,Q5))
                if (cs-1).sign()>=0:tag='parallel'if cs==O else'dotted'
                else:
                    m=next((m for m,a in angles.items()if cs==a),None)
                    tag='m'+str(m)if m else'forbidden'
                types.append(tag)
                if tag=='forbidden':forbidden.append(i)
            row={'vertex':vi,'transverse_labels':list(labels),'old_scale_square':k.enc(),'transverse_part':enc([-x for x in w]),'old_pair_types':types}
            if forbidden:row['forbidden_old_indices']=forbidden;bad.append(row)
            else:caps.append(row)
    # Integer strict escape chosen from whole-bank numerical proposals, then
    # checked exactly against every accepted profile. Normal sqrt(k)*v + b.
    EE=[30861, 39915, 29906, 19910, 9932, -705, -7565, -10211, -10432, -10624]
    # Freeze a robust inward-rounded proposal; discovery never supplies predicates.
    E=[L(x)for x in EE]
    assert pair(E,E,Q).sign()>0 and pair(E,c,Q).sign()<0
    assert all(pair(E,w,Q).sign()<0 for w in roots)
    gaps=[]
    for row in caps:
        vi=row['vertex'];v=V[vi]
        k=L(K(*row['old_scale_square'][0]),K(*row['old_scale_square'][1]))
        b=[L(K(*x[0]),K(*x[1]))for x in row['transverse_part']]
        a=pair(v,E[:6],Q5);z=pair(b,E[6:],H);assert a.sign()<0
        # a sqrt(k)+z < 0, by sign or one exact squared comparison.
        gap=k*a*a-z*z
        assert z.sign()<=0 or gap.sign()>0
        gaps.append({'old_pair':a.enc(),'transverse_pair':z.enc(),'squared_gap':gap.enc()})
    # Trace-field obstruction: end-pair squared has nonzero sqrt(R) part.
    p=pair(W[5],W[6],Q5);cycle=p*p/2
    assert cycle.b and inertia(G)==[4,1]
    # The involution sqrt(R)->-sqrt(R) fixes this indefinite principal G.
    return {'status':'EXACT_MIXED_COMPACT_H5_A4_SEED_BOUNDED_GLOBAL_FAMILY_FAILURE',
      'source_hashes':{str(DEP.relative_to(ROOT)):hashlib.sha256(DEP.read_bytes()).hexdigest()},
      'field':{'base':'Q(sqrt5)','extension':'sqrt(R)>0','R':R.enc(),'norm_R':str(N),'nonnormal_quartic':True},
      'H9_form':mat(Q),'seed_normals':mat(roots),'seed_gram':mat([[pair(a,b,Q)for b in roots]for a in roots]),
      'strict_center':enc(c),'strict_center_norm':pair(c,c,Q).enc(),
      'compact_H5_vertices':mat(V),'compact_H5_active_sets':inds,'positive_dual_center_coefficients':enc(cert),
      'H5_center':enc(h),'A4_inverse':mat(J),'end_cartan_two_cycle':cycle.enc(),
      'seed_pair_types':seed_types,'seed_signature':[9,1],'H5_signature':[5,1],'H5_extreme_ray_count':10,
      'bounded_profile_labels':[2,3,5],'profile_count':810,'nonpositive_old_scale':nonsp,
      'accepted_caps':caps,'rejected_caps':bad,'strict_escape':enc(E),'strict_escape_norm':pair(E,E,Q).enc(),
      'strict_escape_center_pair':pair(E,c,Q).enc(),'accepted_cap_escape_gaps':gaps,
      'scope':'Exact seed, compact H5 face, and finite 810-proposal vertex-centered profile census only. No exhaustive H9 cap search or new compact H9 claim.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=build()
    if a.verify:
        assert d==json.loads(OUT.read_text());print('PASS saved exact seed and bounded global-family certificate')
    else:
        OUT.write_text(json.dumps(d,indent=2)+'\n');print('WROTE',OUT)
    print('accepted',len(d['accepted_caps']),'rejected',len(d['rejected_caps']),'nonpositive',len(d['nonpositive_old_scale']))
if __name__=='__main__':main()
