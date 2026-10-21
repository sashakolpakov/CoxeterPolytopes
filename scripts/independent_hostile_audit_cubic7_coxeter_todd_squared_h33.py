#!/usr/bin/env python3
"""Clean-room hostile replay of the cubic-seven K12^2 H33 cutter.

No production search module is imported.  The script reconstructs the cubic
field, catalogue K12 Gram matrix, complete local norm classes for the
``a=2*theta`` host, all controls through the cutting height, necessary
coefficient Gram types, the explicit root, chronological cross maxima,
same-height acuteness, and the repaired 7/5 replacement.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
from pathlib import Path

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_decomp


ROOT = Path(__file__).resolve().parents[1]
OWNER_MODULE = ROOT / "src/coxeter_polytopes/search/cubic7_coxeter_todd_squared.py"
OWNER_SCRIPT = ROOT / "scripts/audit_cubic7_coxeter_todd_squared_h33_cutter.py"
OWNER_JSON = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h33_cutter.json"
AUDIT_JSON = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h33_independent_hostile_audit.json"
AUDIT_REPORT = ROOT / "results/cubic7_coxeter_todd_squared_h33_independent_hostile_audit.md"


def F(value=0, denominator=None) -> Fraction:
    if denominator is not None:
        return Fraction(value,denominator)
    return value if isinstance(value, Fraction) else Fraction(value)


@dataclass(frozen=True, slots=True)
class C7:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)
    c: Fraction = Fraction(0)

    def __post_init__(self):
        object.__setattr__(self, "a", F(self.a)); object.__setattr__(self, "b", F(self.b)); object.__setattr__(self, "c", F(self.c))

    @staticmethod
    def coerce(x): return x if isinstance(x, C7) else C7(x)
    def __add__(self, x):
        x=C7.coerce(x); return C7(self.a+x.a,self.b+x.b,self.c+x.c)
    __radd__=__add__
    def __neg__(self): return C7(-self.a,-self.b,-self.c)
    def __sub__(self,x): return self+(-C7.coerce(x))
    def __rsub__(self,x): return C7.coerce(x)-self
    def __mul__(self,x):
        x=C7.coerce(x);a,b,c=self.a,self.b,self.c;d,e,f=x.a,x.b,x.c
        return C7(a*d+b*f+c*e-c*f,a*e+b*d+2*b*f+2*c*e-c*f,a*f+b*e+c*d-b*f-c*e+3*c*f)
    __rmul__=__mul__
    def conjugate(self): return C7(self.a-2*self.b+3*self.c,-self.c,self.b-self.c)
    def norm(self):
        x=self*self.conjugate()*self.conjugate().conjugate(); assert not x.b and not x.c; return x.a
    def inverse(self): return self.conjugate()*self.conjugate().conjugate()/self.norm()
    def __truediv__(self,x):
        x=C7.coerce(x)
        return C7(self.a/x.a,self.b/x.a,self.c/x.a) if not x.b and not x.c else self*x.inverse()
    def __rtruediv__(self,x): return C7.coerce(x)/self
    def __pow__(self,n):
        if n<0:return self.inverse()**(-n)
        r=C7(1);b=self
        while n:
            if n&1:r=r*b
            b=b*b;n//=2
        return r
    def conjugate_power(self,n):
        r=self
        for _ in range(n%3):r=r.conjugate()
        return r
    @lru_cache(maxsize=None)
    def sign(self,place=0):
        x=self.conjugate_power(place); lo,hi=F(1),F(2)
        if not (x.a or x.b or x.c):return 0
        while True:
            vals=[x.a+x.b*t+x.c*t*t for t in (lo,hi)]
            vertex=-x.b/(2*x.c) if x.c else None
            if vertex is not None and lo<vertex<hi:vals.append(x.a+x.b*vertex+x.c*vertex*vertex)
            if min(vals)>0:return 1
            if max(vals)<0:return -1
            mid=(lo+hi)/2;p=mid**3+mid**2-2*mid-1
            if p<0:lo=mid
            else:hi=mid
    def __lt__(self,x):return (self-x).sign()<0
    def __le__(self,x):return (self-x).sign()<=0
    @property
    def integral(self):return self.a.denominator==self.b.denominator==self.c.denominator==1
    def json(self):return {"a":str(self.a),"b":str(self.b),"c":str(self.c)}
    def display(self):return f"{self.a}+{self.b}t+{self.c}t2"


T=C7(0,1); A=2*T; CUT_H=C7(33,94,42); CUT_X=C7(1,6,3)
K12=((4,0,0,-2,0,0,2,-1,-1,-1,2,-1),(0,4,0,0,-2,0,2,-1,-1,-1,-1,2),(0,0,4,0,0,-2,2,2,2,-1,-1,-1),(-2,0,0,4,0,0,-1,-1,2,2,-1,-1),(0,-2,0,0,4,0,-1,2,-1,2,-1,-1),(0,0,-2,0,0,4,-1,-1,-1,2,2,2),(2,2,2,-1,-1,-1,4,0,0,-2,0,0),(-1,-1,2,-1,2,-1,0,4,0,0,-2,0),(-1,-1,2,2,-1,-1,0,0,4,0,0,-2),(-1,-1,-1,2,2,2,-2,0,0,4,0,0),(2,-1,-1,-1,-1,2,0,-2,0,0,4,0),(-1,2,-1,-1,-1,2,0,0,-2,0,0,4))
H=(F(4,3),0,0,F(2,3),0,0,F(-2,3),F(1,3),F(1,3),F(-1,3),F(-1,3),F(2,3))
U=(4,0,0,2,0,0,-2,1,1,-1,-1,2);V=(10,0,0,2,0,0,-5,1,4,-1,-4,5);W=(4,0,0,0,0,0,-2,0,2,0,-2,2)

def pair(x,y):return sum(F(x[i])*K12[i][j]*F(y[j]) for i in range(12) for j in range(12))

def gp_matrix(M):return "["+";".join(",".join(map(str,row)) for row in M)+"]"

def gp_minimum(M,bound):
    executable=shutil.which("gp")
    if executable is None:raise RuntimeError("PARI/GP is required for the hostile replay")
    program=f"R=qfminim({gp_matrix(M)},{bound},1,0);print(R[1]);print(floor(R[2]+1/2));quit()\n"
    done=subprocess.run([executable,"-q","-s","256000000"],input=program,text=True,capture_output=True,timeout=60,check=True)
    lines=[line.strip() for line in done.stdout.splitlines() if line.strip()]
    if len(lines)!=2:raise AssertionError(f"unexpected qfminim output: {lines!r}")
    return int(lines[0]),int(lines[1])

def integer_upper(value,place):
    if value.sign(place)<0:raise ValueError("negative embedding cap")
    upper=0 if value==C7() else 1
    while (C7(upper)-value).sign(place)<0:upper*=2
    return upper

DUAL_ABS=((F(1,6),F(1),F(1,7)),(F(4,9),F(1,4),F(1,5)),(F(1,5),F(4,9),F(1,4)))

def ceil_fraction(x):return -((-x.numerator)//x.denominator)

def ceil_sqrt(n):
    f=isqrt(n)
    return f if f*f==n else f+1

def bounded(caps):
    # Trace-dual coefficient bounds, followed by exact embedding-cap signs.
    upper=tuple(integer_upper(cap,e) for e,cap in enumerate(caps))
    embedding=tuple(ceil_sqrt(n) for n in upper)
    bounds=tuple(ceil_fraction(sum(embedding[e]*DUAL_ABS[i][e] for e in range(3))) for i in range(3))
    return tuple(C7(a,b,c) for a in range(-bounds[0],bounds[0]+1) for b in range(-bounds[1],bounds[1]+1) for c in range(-bounds[2],bounds[2]+1) if all((caps[e]-C7(a,b,c)**2).sign(e)>=0 for e in range(3)))

@dataclass(frozen=True,slots=True,order=True)
class Level: height:C7; norm:C7; control:C7; target:C7

@lru_cache(None)
def controls():
    out=set()
    for q in (C7(2),C7(6)):
        for x in bounded((q*CUT_H/A,q/(-A),q/(-A))):
            if x.sign()<=0 or not (2*A*x/q).integral:continue
            target=q+A*x*x
            if any(target.sign(e)<0 for e in (1,2)):continue
            h=A*x*x/q
            if h<=CUT_H:out.add(Level(h,q,x,target))
    return tuple(sorted(out))

def det3(v):
    u,vv,w,uv,uw,vw=v;return u*vv*w+2*uv*uw*vw-u*vw*vw-vv*uw*uw-w*uv*uv

@lru_cache(maxsize=None)
def quotient(M):
    R=Matrix(M)
    if R.det():return M
    D,L,Rt=smith_normal_decomp(R,domain=ZZ); nz=[i for i in range(R.rows) if D[i,i]]
    assert L*R*Rt==D and abs(int(L.det()))==abs(int(Rt.det()))==1
    zero=[i for i in range(R.rows) if i not in nz]
    assert not zero or R*Rt[:,zero]==Matrix.zeros(R.rows,len(zero))
    C=Rt[:,nz];Q=C.T*R*C
    assert Q.det()>0
    return tuple(tuple(int(Q[i,j]) for j in range(len(nz))) for i in range(len(nz)))

@lru_cache(maxsize=None)
def min_at_least(M,m):
    Q=quotient(M)
    if not Q:return False
    R=Matrix(Q);I=R.inv();bs=[isqrt(int((I[i,i]*(m-1)).p)//int((I[i,i]*(m-1)).q)) for i in range(R.rows)]
    return all(not(0<sum(Q[i][j]*z[i]*z[j] for i in range(len(z)) for j in range(len(z)))<m) for z in product(*(range(-b,b+1) for b in bs)) if any(z))

@lru_cache(None)
def types(level):
    assert level.target.integral
    a,b,c=map(int,(level.target.a,level.target.b,level.target.c));tr=int(3*a-b+5*c);out=[];ternary=level.norm==C7(6);minimum=12 if ternary else 4;mod=3 if ternary else 1
    for vv in range(0,2*tr//7+1,2):
      for w in range(0,2*tr//7+1,2):
       lower=max(-isqrt(vv*w),(a+w-tr+1)//2);upper=min(isqrt(vv*w),(a+w)//2)
       for vw in range(lower,upper+1):
        u=a-2*vw+w; tuv=b+w-4*vw;tuw=c-vv-3*w+2*vw
        if u<0 or u>tr or u%2 or tuv%2 or tuw%2:continue
        uv,uw=tuv//2,tuw//2; vals=(u,vv,w,uv,uw,vw)
        if any(0<x<minimum for x in (u,vv,w)) or any(x%mod for x in vals):continue
        if uv*uv>u*vv or uw*uw>u*w or det3(vals)<0:continue
        M=((u,uv,uw),(uv,vv,vw),(uw,vw,w))
        if min_at_least(M,minimum):out.append(vals)
    return tuple(sorted(out))

@lru_cache(maxsize=1)
def local_arithmetic():
    # x^3+x^2+1 has no F_2 root, so (2) is inert and squarefree.  The
    # discriminant is 49; Minkowski's cubic bound is 14/9<2, proving h=1.
    inert_two=all((x**3+x**2+1)%2 for x in (0,1))
    minkowski=Fraction(14,9)
    units=(C7(-1),T,T.conjugate())
    assert all(abs(u.norm())==1 for u in units)
    sign_rows=tuple(tuple(int(u.sign(e)<0) for e in range(3)) for u in units)
    # Row-reduce the 3-by-3 signature matrix over F_2.
    rows=[sum(bit<<j for j,bit in enumerate(row)) for row in sign_rows];rank=0
    for bit in range(3):
        pivot=next((i for i in range(rank,3) if rows[i]&(1<<bit)),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        for i in range(3):
            if i!=rank and rows[i]&(1<<bit):rows[i]^=rows[rank]
        rank+=1
    raw=Matrix(K12);D,L,R=smith_normal_decomp(raw,domain=ZZ)
    assert L*raw*R==D
    smith=tuple(sorted(abs(int(D[i,i])) for i in range(12)))
    assert smith==(1,)*6+(3,)*6
    factors=tuple(3//gcd(3,abs(int(D[i,i]))) for i in range(12))
    radical=hermite_normal_form(R*Matrix.diag(*factors))
    radical_gram=radical.T*raw*radical
    radical_q=tuple(tuple(int(radical_gram[i,j]) for j in range(12)) for i in range(12))
    assert all(int(x)%3==0 for x in raw*radical)
    trace_dual=(C7(1,F(-1,7),F(-3,7)),C7(F(-1,7),F(2,7),F(1,7)),C7(F(-3,7),F(1,7),F(2,7)))
    assert all((C7(DUAL_ABS[i][e])-s*trace_dual[i]).sign(e)>0 for i in range(3) for e in range(3) for s in (-1,1))
    ordinary=gp_minimum(K12,4);radical_min=gp_minimum(radical_q,12)
    assert ordinary==(756,4) and radical_min==(756,12)
    # Normalize tuples to their JSON-array representation so an in-memory
    # replay compares exactly with the frozen JSON artifact.
    return json.loads(json.dumps({
        "minimal_polynomial_mod_2_values": [1,1],
        "two_is_inert": inert_two,
        "field_discriminant":49,
        "minkowski_bound":str(minkowski),
        "class_number_one":minkowski<2,
        "unit_signature_rows":sign_rows,
        "unit_signature_rank_mod_2":rank,
        "totally_positive_units_are_squares":rank==3,
        "k12_smith":smith,
        "k12_qfminim":{"ordinary":ordinary,"radical_three":radical_min},
        "canonical_norms":[2,6],
        "norm_gate":(
            "At inert (2), evenness plus primitivity forces valuation one. "
            "Away from (2),(3), local unimodularity forces valuation zero. "
            "At (3), the elementary Smith factor makes the valuation at most "
            "one; it is optional. Class number one and U^+=U^2 leave q=2,6."
        ),
    }))

def schur_pairs(square):
    # Cross columns (a,c) against the fixed cutter basis.  The congruence is
    # forced by W_c=2(V_c-U_c)/3 being a lattice vector.
    out=[]
    for a in range(-isqrt(12*square)-1,isqrt(12*square)+2):
      for c in range(-isqrt(84*square)-1,isqrt(84*square)+2):
        residual=18*square-(14*a*a-10*a*c+2*c*c)
        if residual>=0 and (c-a)%3==0:out.append((a,c,residual))
    return tuple(out)

@lru_cache(maxsize=None)
def cross_rank_one(control):
    scalar=C7(0,2,1);alpha=C7(1,0,F(-2,3));beta=C7(0,1,F(2,3));best=None;kept=0
    for a,c,_residual in schur_pairs(4):
        M=((12,30,a),(30,84,c),(a,c,4))
        if not min_at_least(M,4):continue
        kept+=1;value=-A*CUT_X*control+scalar*(alpha*a+beta*c)
        if best is None or best<value:best=value
    assert best is not None
    return best,kept

@lru_cache(maxsize=None)
def cross_rank_two(base,control,minimum_gate=True):
    # Earlier coefficient vector is U+d*V, d=theta+theta^2/2.  Its third
    # coefficient is V/2, so both cross pairings with V must be even.
    p,r,s=base;delta=C7(0,1,F(1,2));alpha=C7(1,0,F(-2,3));beta=C7(0,1,F(2,3));best=None;schur=kept=0
    left=schur_pairs(p);right=schur_pairs(s)
    for a,c,r11 in left:
      for b,d,r22 in right:
        if b%2 or d%2:continue
        r12=18*r-(14*a*b-5*(a*d+c*b)+2*c*d)
        if r11*r22<r12*r12:continue
        schur+=1
        M=((12,30,a,b),(30,84,c,d),(a,c,p,r),(b,d,r,s))
        if minimum_gate and not min_at_least(M,4):continue
        kept+=1;value=-A*CUT_X*control+alpha*a+alpha*delta*b+beta*c+beta*delta*d
        if best is None or best<value:best=value
    assert best is not None
    return best,schur,kept

@lru_cache(maxsize=1)
def chronological_cross():
    first=cross_rank_one(C7(1,2,1))
    second=cross_rank_two((4,10,40),C7(1,4,2))
    third=cross_rank_two((4,14,64),C7(2,5,2))
    other_same=cross_rank_two((8,24,80),CUT_X,False)
    return first,second,third,other_same

@lru_cache(maxsize=None)
def directional_max(base,relation):
    # h=G^-1 e_1 has h^2=4/3 and h.X=X_1.  Enumerating its integral
    # pairings with a rank-one/two coefficient Gram is therefore complete.
    if len(base)==1:
        square=base[0];best=None;kept=0
        for m in range(-isqrt((4*square)//3)-2,isqrt((4*square)//3)+3):
            if 3*m*m>4*square:continue
            kept+=1;value=C7(0,2,1)*m
            if best is None or best<value:best=value
        return best,kept
    p,r,s=base;best=None;kept=0
    for u in range(-isqrt((4*p)//3)-2,isqrt((4*p)//3)+3):
      for v in range(-isqrt((4*s)//3)-2,isqrt((4*s)//3)+3):
        # PSD of [[4/3,u,v],[u,p,r],[v,r,s]], cleared by 3.
        if 4*p<3*u*u or 4*s<3*v*v:continue
        if 4*(p*s-r*r)<3*(s*u*u-2*r*u*v+p*v*v):continue
        if relation=="half" and v%2:continue
        if relation=="thirds" and (v-u)%3:continue
        kept+=1
        if relation=="half":value=C7(u)+C7(0,1,F(1,2))*v
        else:value=C7(1,0,F(-2,3))*u+C7(0,1,F(2,3))*v
        if best is None or best<value:best=value
    assert best is not None
    return best,kept

@lru_cache(maxsize=1)
def same_height():
    alpha=C7(1,0,F(-2,3));beta=C7(0,1,F(2,3));best=None;schur=kept=0
    for a in range(-12,13):
     for c in range(-31,32):
      r11=216-(14*a*a-10*a*c+2*c*c)
      if r11<0:continue
      for b in range(-31,32):
       for d in range(-84,85):
        r22=1512-(14*b*b-10*b*d+2*d*d);r12=540-(14*a*b-5*(a*d+c*b)+2*c*d)
        if r22<0 or r11*r22<r12*r12 or (b-a)%3 or (d-c)%3:continue
        schur+=1
        if (a,b,c,d)==(12,30,30,84):continue
        M=((12,30,a,b),(30,84,c,d),(a,c,12,30),(b,d,30,84))
        if not min_at_least(M,4):continue
        kept+=1; val=-A*CUT_X*CUT_X+alpha*alpha*a+alpha*beta*(b+c)+beta*beta*d
        if best is None or best<val:best=val
    return schur,kept,best

@lru_cache(maxsize=1)
def build():
    arithmetic=local_arithmetic()
    assert Matrix(K12).det()==729 and pair(H,H)==F(4,3)
    lv=controls(); non=[(x,types(x)) for x in lv if types(x)]
    frontier=[(x.height,x.norm,len(t)) for x,t in non]
    assert len(lv)==11
    assert frontier==[(C7(5,14,6),C7(2),1),(C7(16,45,20),C7(2),1),(C7(25,70,31),C7(2),1),(CUT_H,C7(2),2)]
    matrices=[tuple(((v[0],v[3],v[4]),(v[3],v[1],v[5]),(v[4],v[5],v[2])) for v in t) for _x,t in non]
    assert matrices==[
        (((0,0,0),(0,16,8),(0,8,4)),),
        (((4,10,5),(10,40,20),(5,20,10)),),
        (((4,14,7),(14,64,32),(7,32,16)),),
        (((8,24,12),(24,80,40),(12,40,20)),((12,30,12),(30,84,36),(12,36,16))),
    ]
    G=tuple(tuple(int(pair(x,y)) for y in (U,V,W)) for x in (U,V,W));assert G==((12,30,12),(30,84,36),(12,36,16));assert all(2*u-2*v+3*w==0 for u,v,w in zip(U,V,W))
    spatial=sum((T**i)*(T**j)*int(pair((U,V,W)[i],(U,V,W)[j])) for i in range(3) for j in range(3));assert -A*CUT_X*CUT_X+spatial==C7(2)
    old=-A*CUT_X+C7(F(3,2))*C7(U[0],V[0],W[0]);new=-A*CUT_X+C7(F(7,5))*C7(U[0],V[0],W[0]);assert old==T and new==C7(F(-2,5),0,F(-2,5))
    cross=chronological_cross()
    assert cross[0]==(C7(0),43)
    assert cross[1]==(C7(0,-3,-2),1117,1081)
    assert cross[2]==(C7(-4,-8,-3),1117,1081)
    assert cross[3]==(C7(-2,-8,-4),1225,1225)
    schur,kept,same=same_height();assert (schur,kept,same)==(8173,6960,C7(0,-4,-2))
    chronological=(cross[0][0],cross[1][0],cross[2][0],cross[3][0],same)
    assert all(x.sign()<=0 for x in chronological)
    directions=(directional_max((4,),"rank1"),directional_max((4,10,40),"half"),directional_max((4,14,64),"half"),directional_max((8,24,80),"half"),directional_max((12,30,84),"thirds"))
    maxima=tuple(x[0] for x in directions)
    assert maxima==(C7(0,4,2),C7(2,6,3),C7(2,8,4),C7(3,10,5),C7(4,10,4))
    level_controls=(C7(1,2,1),C7(1,4,2),C7(2,5,2),CUT_X,CUT_X)
    replacement_bounds=tuple(-A*x+F(7,5)*d for x,d in zip(level_controls,maxima))
    assert replacement_bounds==(C7(-2,F(-2,5),F(4,5)),C7(F(-6,5),F(-8,5),F(1,5)),C7(F(-6,5),F(-4,5),F(-2,5)),C7(F(-9,5),0,1),C7(F(-2,5),0,F(-2,5)))
    assert all(x.sign()<0 for x in replacement_bounds)
    omitted_old_bound=-A*CUT_X+F(36,25)*maxima[3]
    assert omitted_old_bound==C7(F(-42,25),F(2,5),F(6,5)) and omitted_old_bound.sign()>0
    newnorm=C7(F(196,75))-A;assert newnorm.sign()>0
    owner=json.loads(OWNER_JSON.read_text())
    assert owner["replacement_hyperideal_ray"]["scale"]=="7/5"
    assert owner["first_cutting_wall"]["same_height"]["minimum_four_candidate_count"]==kept
    payload={
        "schema":"coxeter.cubic7_k12_h33.hostile.v1",
        "imports_production_search_code":False,
        "owner_sha256":hashlib.sha256(OWNER_JSON.read_bytes()).hexdigest(),
        "owner_module_sha256":hashlib.sha256(OWNER_MODULE.read_bytes()).hexdigest(),
        "owner_script_sha256":hashlib.sha256(OWNER_SCRIPT.read_bytes()).hexdigest(),
        "local_arithmetic":arithmetic,
        "control_count":len(lv),
        "frontier":[{"height":h.json(),"norm":q.json(),"types":n,"matrices":M} for (h,q,n),M in zip(frontier,matrices)],
        "explicit_root":{"gram":G,"relation":"2U-2V+3W=0","norm":C7(2).json(),"primitive":True,"reflection_integral":True,"old_pairing":old.json()},
        "chronological_cross":{
            "first":{"maximum":cross[0][0].json(),"minimum_four_count":cross[0][1]},
            "second":{"maximum":cross[1][0].json(),"schur_count":cross[1][1],"minimum_four_count":cross[1][2]},
            "third":{"maximum":cross[2][0].json(),"schur_count":cross[2][1],"minimum_four_count":cross[2][2]},
            "other_H33_type":{"maximum":cross[3][0].json(),"relaxed_schur_count":cross[3][1]},
            "same_cutter_type":{"maximum":same.json(),"schur_count":schur,"minimum_four_count":kept},
        },
        "hostile_finding":{"kind":"REPAIRED_PROOF_GAP","detail":"The original 36/25 ledger omitted the other necessary H33 type.","omitted_type_relaxed_pairing_bound_at_36_over_25":omitted_old_bound.json(),"positive":True},
        "replacement":{"scale":"7/5","pairing_with_cutter":new.json(),"norm":newnorm.json(),"directional_maxima":[x.json() for x in maxima],"directional_candidate_counts":[x[1] for x in directions],"strict_pairing_bounds":[x.json() for x in replacement_bounds]},
        "bounded_scope_only":True,
        "compact_or_finite_volume_claim":False,
        "terminal_verdict":"PASS",
    }
    # The frozen replay API is deliberately JSON-native so direct equality
    # with json.loads(artifact) cannot hide tuple/list drift.
    return json.loads(json.dumps(payload))

def report(p):return f"""# Independent hostile audit: cubic-seven K12^2 H33 cutter

Terminal verdict: **{p['terminal_verdict']}**.

This clean-room script imports no production search code.  It reconstructs
the cubic arithmetic and local norms `2,6`, the raw K12 Gram and ternary
radical gates, all 11 controls, and all five surviving necessary Gram types.
It independently reproduces the three earlier cross maxima.  At H33 it
checks both types: the relaxed cross maximum against the other type is
`-2-8theta-4theta^2`, while the cutter-type enumeration reproduces 8,173
Schur blocks, 6,960 minimum-four blocks, and maximum
`-4theta-2theta^2`.

The audit found a real proof gap in the first production certificate: its
`36/25` directional ledger omitted the other H33 type, whose relaxed bound
would be positive.  Production is repaired to scale `7/5`.  Its norm is
`196/75-2theta>0`, and exact necessary-Gram maxima give strict negative
pairings at all five types, including both H33 types.  This terminal PASS is
for the repaired bounded certificate only; it makes no compactness,
finite-volume, or Vinberg-termination claim.
"""

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--verify",action="store_true");a=ap.parse_args();p=build();j=json.dumps(p,indent=2,sort_keys=True)+"\n";r=report(p)
    if a.verify:
        if AUDIT_JSON.read_text()!=j or AUDIT_REPORT.read_text()!=r:raise SystemExit("hostile artifact drift")
        print("PASS");return
    AUDIT_JSON.write_text(j);AUDIT_REPORT.write_text(r);print("PASS")

if __name__=="__main__":main()
