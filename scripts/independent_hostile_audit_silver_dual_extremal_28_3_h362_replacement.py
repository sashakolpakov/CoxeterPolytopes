#!/usr/bin/env python3
"""Owner-import-free hostile replay of the silver rank-28 H362 continuation."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import sys
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OWNER = ROOT / "src/coxeter_polytopes/search/silver_dual_extremal_28_3.py"
PRODUCTION = ROOT / "data/enumerations/silver_dual_extremal_28_3_h362_replacement.json"
OUTPUT = ROOT / "data/enumerations/silver_dual_extremal_28_3_h362_replacement_independent_hostile_audit.json"
REPORT = ROOT / "results/silver_dual_extremal_28_3_h362_replacement_independent_hostile_audit.md"
EXPECTED_GRAM_SHA = "d61f9a887f3b31db364b6c38c9b0e59b615c4cd178201ae3e8d0036a87629c54"


@dataclass(frozen=True, slots=True)
class K2:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @classmethod
    def make(cls, a=0, b=0): return cls(Fraction(a), Fraction(b))
    @classmethod
    def coerce(cls, x): return x if isinstance(x, cls) else cls.make(x)
    def __add__(self, x): x=self.coerce(x); return K2(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self): return K2(-self.a,-self.b)
    def __sub__(self,x): return self+(-self.coerce(x))
    def __rsub__(self,x): return self.coerce(x)-self
    def __mul__(self,x): x=self.coerce(x); return K2(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self.coerce(x); n=x.a*x.a-2*x.b*x.b
        return self*K2(x.a/n,-x.b/n)
    def sign(self, place=0):
        a,b=(self.a,-self.b) if place else (self.a,self.b)
        if not a:return (b>0)-(b<0)
        if not b or (a>0)==(b>0):return (a>0)-(a<0)
        return ((a>0)-(a<0)) if a*a>2*b*b else ((b>0)-(b<0))
    def __lt__(self,x): return (self-self.coerce(x)).sign()<0
    def __le__(self,x): return (self-self.coerce(x)).sign()<=0
    @property
    def integral(self): return self.a.denominator==self.b.denominator==1
    def __str__(self):
        if not self.b:return str(self.a)
        if not self.a:return f"{self.b}*sqrt(2)"
        return f"{self.a} {'+' if self.b>0 else '-'} {abs(self.b)}*sqrt(2)"


ZERO=K2(); ONE=K2.make(1); SQRT2=K2.make(0,1); RHO=2*SQRT2-ONE; P7=K2.make(3,1)
BOUND=K2.make(362,256); THREAT=K2.make(362,182); SCALE=Fraction(12,25)
LOCAL_THREE = tuple(c for c in product(range(-3,4), repeat=3) if any(c))
LOCAL_FOUR = tuple(c for c in product(range(-4,5), repeat=4) if any(c))


def literal_gram():
    syntax=ast.parse(OWNER.read_text()); text=None
    for node in syntax.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="_GRAM_TEXT" for t in node.targets):
            text=ast.literal_eval(node.value)
    rows=tuple(tuple(map(int,line.split())) for line in text.splitlines() if line)
    stream="".join(" ".join(map(str,row))+"\n" for row in rows)
    if hashlib.sha256(stream.encode()).hexdigest()!=EXPECTED_GRAM_SHA: raise AssertionError("Gram changed")
    return rows


@dataclass(frozen=True,slots=True)
class Gate: norm:K2; dyadic:bool; p7:bool; ternary:bool
GATES=tuple(Gate(K2.make(2 if d else 1)*(P7 if p else 1)*(3 if t else 1),d,p,t) for d,p,t in product((False,True),repeat=3))


def ceil_sqrt(n):
    r=isqrt(n); return r if r*r==n else r+1


def embedding_upper(x,place):
    u=0 if x==ZERO else 1
    while (K2.make(u)-x).sign(place)<0:u*=2
    return u


def controls():
    out=set()
    for gate in GATES:
        cap0=gate.norm*BOUND/RHO;cap1=gate.norm/(-RHO)
        n=(ceil_sqrt(embedding_upper(cap0,0))+ceil_sqrt(embedding_upper(cap1,1))+1)//2
        for a,b in product(range(-n,n+1),repeat=2):
            x=K2.make(a,b);sq=x*x
            if (sq-cap0).sign(0)>0 or (sq-cap1).sign(1)>0 or x.sign(0)<=0:continue
            if not (2*RHO*x/gate.norm).integral:continue
            target=gate.norm+RHO*sq;height=RHO*sq/gate.norm
            if target.sign(1)<0 or (height-BOUND).sign(0)>0:continue
            out.add((height,gate.norm,x,target,gate))
    return tuple(sorted(out,key=lambda z:(z[0],z[1],z[2])))


def square(t,a,b):return a*a*t[0]+b*b*t[1]+2*a*b*t[2]


def short_binary(t,m):
    D=t[0]*t[1]-t[2]*t[2]
    if D<0:return True
    if D==0:return gcd(gcd(t[0],t[1]),abs(t[2]))<m
    A=isqrt(((m-1)*t[1])//D)+1;B=isqrt(((m-1)*t[0])//D)+1
    return any(0<square(t,a,b)<m for a in range(-A,A+1) for b in range(-B,B+1) if a or b)


def triples(row):
    height,norm,x,target,gate=row
    target=target/(RHO*RHO) if gate.p7 else target
    if not target.integral or (target.b/2).denominator!=1:return ()
    A=int(target.a);M=int(target.b/2);minimum=12 if gate.ternary else 8;out=[]
    for v in range(A//2+1):
        u=A-2*v;t=(u,v,M)
        if u<0 or u%2 or v%2 or (u and u<minimum) or (v and v<minimum) or M*M>u*v:continue
        if gate.ternary and any(z%3 for z in t):continue
        if not short_binary(t,minimum):out.append(t)
    return tuple(out)


def augmented(t):
    u,v,m=t;best=None;pair=None;detbest=0;count=0
    for p in range(-isqrt(8*u)-2,isqrt(8*u)+3):
        for q in range(-isqrt(8*v)-2,isqrt(8*v)+3):
            D=8*(u*v-m*m)-v*p*p-u*q*q+2*m*p*q
            if D<0:continue
            if any(
                0 < 8*a*a + u*b*b + v*c*c + 2*p*a*b + 2*q*a*c + 2*m*b*c < 8
                for a,b,c in LOCAL_THREE
            ):continue
            count+=1;z=K2.make(p,q)
            if best is None or best<z:best,pair,detbest=z,(p,q),D
    return best,pair,detbest,count


def dot(G,x,y):return sum(x[i]*G[i][j]*y[j] for i in range(28) for j in range(28) if x[i] and y[j] and G[i][j])


SPECS = (
    (K2.make(93,66),K2.make(3),K2.make(9,6),(138,72,99),((12,3),(3,12)),K2.make(-3,-2),K2.make(-1,-1),3),
    (K2.make(93,66),K2.make(6),K2.make(12,9),(300,132,198),((12,6),(6,36)),K2.make(-1,-1),K2.make(3,2),3),
    (K2.make(93,66),K2.make(6),K2.make(12,9),(288,138,198),((18,0),(0,30)),K2.make(-1,-1),K2.make(3,2),3),
    (K2.make(93,66),K2.make(6),K2.make(12,9),(276,144,198),((24,6),(6,24)),K2.make(-3,-2),K2.make(-1,-1),3),
    (K2.make(93,66),K2.make(6),K2.make(12,9),(264,150,198),((18,-6),(-6,24)),K2.make(-4,-3),K2.make(-1,-1),3),
    (K2.make(93,66),K2.make(6),K2.make(12,9),(252,156,198),((12,6),(6,12)),K2.make(-4,-3),K2.make(-1,-1),3),
    (K2.make(124,88),K2.make(2),K2.make(8,6),(122,64,88),((8,4),(4,10)),K2.make(-4,-3),K2.make(3,2),1),
    (K2.make(181,128),K2.make(2),K2.make(10,7),(188,88,128),((8,0),(0,20)),K2.make(-1,-1),K2.make(3,2),1),
    (K2.make(181,128),K2.make(2),K2.make(10,7),(184,90,128),((10,-2),(-2,18)),K2.make(-1,-1),K2.make(3,2),1),
    (K2.make(181,128),K2.make(2),K2.make(10,7),(180,92,128),((12,-4),(-4,16)),K2.make(-1,-1),K2.make(3,2),1),
    (K2.make(181,128),K2.make(2),K2.make(10,7),(176,94,128),((14,6),(6,14)),K2.make(-3,-2),K2.make(-1,-1),1),
    (K2.make(181,128),K2.make(2),K2.make(10,7),(172,96,128),((12,4),(4,12)),K2.make(-1,-1),K2.make(4,3),1),
    (K2.make(181,128),K2.make(2),K2.make(10,7),(168,98,128),((8,0),(0,10)),K2.make(-1,-1),K2.make(4,3),1),
    (K2.make(186,132),K2.make(6),K2.make(18,12),(594,264,396),((66,0),(0,0)),K2.make(3,2),ZERO,3),
    (K2.make(186,132),K2.make(6),K2.make(18,12),(570,276,396),((12,0),(0,42)),K2.make(-4,-3),K2.make(3,2),3),
    (K2.make(186,132),K2.make(6),K2.make(18,12),(558,282,396),((18,0),(0,30)),K2.make(-4,-3),K2.make(3,2),3),
    (K2.make(186,132),K2.make(6),K2.make(18,12),(546,288,396),((18,0),(0,24)),K2.make(3,2),K2.make(-4,-3),3),
)


def cross_bound(spec):
    """Independent Schur-complement census against the explicit cutter plane."""

    height,norm,control,triple,B,gamma,delta,modulus=spec;det=348
    def columns(n):
        pb=isqrt(14*n);rb=isqrt(26*n);ps=-(pb//modulus)*modulus;rs=-(rb//modulus)*modulus
        return tuple((p,r,z) for p in range(ps,pb+1,modulus) for r in range(rs,rb+1,modulus) if (z:=26*p*p-8*p*r+14*r*r)<=det*n)
    C1=columns(B[0][0]);C2=columns(B[1][1]);fa=K2.make(-4,-3);fb=K2.make(3,2)
    weights=(fa*gamma,fa*delta,fb*gamma,fb*delta);weights=tuple((int(z.a),int(z.b)) for z in weights);time=RHO*K2.make(14,10)*control;time=(int(time.a),int(time.b))
    relaxed=None;admitted=None;relaxed_entries=None;admitted_entries=None;positive=[];feasible=0
    def sign(pair):
        a,b=pair
        if not a:return (b>0)-(b<0)
        if not b or (a>0)==(b>0):return (a>0)-(a<0)
        return ((a>0)-(a<0)) if a*a>2*b*b else ((b>0)-(b<0))
    for p,r,n1 in C1:
        R1=det*B[0][0]-n1
        for q,s,n2 in C2:
            R2=det*B[1][1]-n2;mix=26*p*q-4*(p*s+r*q)+14*r*s;R12=det*B[0][1]-mix
            if R1*R2<R12*R12:continue
            feasible+=1;pair=(weights[0][0]*p+weights[1][0]*q+weights[2][0]*r+weights[3][0]*s-time[0],weights[0][1]*p+weights[1][1]*q+weights[2][1]*r+weights[3][1]*s-time[1]);entries=(p,q,r,s)
            if relaxed is None or sign((pair[0]-relaxed[0],pair[1]-relaxed[1]))>0:relaxed,relaxed_entries=pair,entries
            if sign(pair)>0:
                G=((14,4,p,q),(4,26,r,s),(p,r,B[0][0],B[0][1]),(q,s,B[1][0],B[1][1]))
                cert=next(((c,z) for c in LOCAL_FOUR if 0<(z:=sum(c[i]*G[i][j]*c[j] for i in range(4) for j in range(4)))<8),None)
                if cert is None:raise AssertionError("unresolved positive cross extension")
                positive.append((entries,pair,cert[0],cert[1]));continue
            if admitted is None or sign((pair[0]-admitted[0],pair[1]-admitted[1]))>0:admitted,admitted_entries=pair,entries
    return {"height":str(height),"norm":str(norm),"gram":list(triple),"reduced_gram":[list(z) for z in B],"feasible":feasible,"relaxed":str(K2.make(*relaxed)),"relaxed_entries":list(relaxed_entries),"positive":[{"entries":list(e),"pairing":str(K2.make(*p)),"short_combination":list(c),"short_square":z} for e,p,c,z in positive],"admitted":str(K2.make(*admitted)),"admitted_entries":list(admitted_entries)}


def audit():
    G=literal_gram();rows=controls()
    nonempty=tuple((r,triples(r)) for r in rows if triples(r))
    observed=tuple((r[0],r[1],len(ts)) for r,ts in nonempty)
    expected=((K2.make(93,66),K2.make(3),1),(K2.make(93,66),K2.make(6),5),(K2.make(124,88),K2.make(2),1),(K2.make(181,128),K2.make(2),6),(K2.make(186,132),K2.make(6),4),(BOUND,K2.make(2),8))
    if len(rows)!=103 or len({r[0] for r in rows})!=70 or observed!=expected:raise AssertionError("chronology changed")
    first=None;margins=[]
    for r,ts in nonempty:
        for t in ts:
            M,p,D,n=augmented(t);old=RHO*r[2]-M/2;new=RHO*r[2]-K2.make(SCALE)*M
            margins.append(new)
            if old.sign(0)<=0 and first is None:first=(r,t,M,p,D,n,old)
    if first[0][0]!=BOUND or first[1]!=(362,182,256) or first[2]!=K2.make(53,37) or first[3]!=(53,37) or first[4]!=0 or first[6]!=K2.make(Fraction(-1,2),Fraction(-1,2)):raise AssertionError("first threat changed")
    d=(1,)+(0,)*27
    v=(4,0,1,3,1,-1,-2,1,0,1,1,1,0,0,2,1,-1,1,0,-1,0,0,-1,0,-1,1,1,-2)
    u=tuple(v[i]+2*d[i] for i in range(28))
    if (dot(G,d,d),dot(G,d,u),dot(G,d,v),dot(G,u,u),dot(G,v,v),dot(G,u,v))!=(8,53,37,362,182,256):raise AssertionError("explicit realization changed")
    old_pair=-RHO*K2.make(14,10)+K2.make(53,37)/2
    replacement_norm=K2.make(8*SCALE*SCALE)-RHO
    if old_pair!=K2.make(Fraction(1,2),Fraction(1,2)) or replacement_norm.sign(0)<=0 or any(x.sign(0)<=0 for x in margins):raise AssertionError("witness conclusion changed")
    acceptance=tuple(cross_bound(spec) for spec in SPECS)
    if any(row["admitted"] not in ("0","-1","-1 - 1*sqrt(2)","-3 - 3*sqrt(2)","-6 - 3*sqrt(2)","-9 - 6*sqrt(2)") for row in acceptance):raise AssertionError("positive admitted earlier pairing")
    if sum(len(row["positive"]) for row in acceptance)!=2 or [row["positive"][0]["short_square"] for row in acceptance if row["positive"]]!=[6,2]:raise AssertionError("short-vector exclusion census changed")
    production=json.loads(PRODUCTION.read_text())
    if production["chronology"]["control_count"]!=103 or not production["first_cutting_wall"]["accepted"]:raise AssertionError("production disagrees")
    if not production["chronology"]["binary_minimum_gate"]["repaired_former_fixed_box_gate"] or len(production["chronology"]["binary_minimum_gate"]["removed_spurious_types"])!=6:raise AssertionError("production minimum-gate repair changed")
    prod_acceptance=production["first_cutting_wall"]["earlier_packet_cross_gram_bounds"]
    if [row["feasible_psd_cross_grams"] for row in prod_acceptance]!=[row["feasible"] for row in acceptance] or [row["admitted_maximum_root_pairing"]["display"] for row in prod_acceptance]!=[row["admitted"] for row in acceptance]:raise AssertionError("production acceptance table disagrees")
    return {
        "schema_version":"0.1","verdict":"PASS","scope":"Owner-import-free exact replay of chronology, first threat, explicit realization, all earlier cross-Gram acceptance bounds, and strict replacement through H362.",
        "independence":"Owner source parsed only as inert AST Gram data; no owner imports.",
        "literal_gram_sha256":EXPECTED_GRAM_SHA,"control_count":len(rows),"height_count":len({r[0] for r in rows}),
        "nonempty_levels":[{"height":str(r[0]),"norm":str(r[1]),"type_count":len(ts)} for r,ts in nonempty],
        "first_threat":{"height":str(first[0][0]),"norm":str(first[0][1]),"control":str(first[0][2]),"gram":list(first[1]),"maximum":str(first[2]),"pair":list(first[3]),"determinant":first[4],"old_witness_pairing":str(old_pair),"explicit_relation":"u=2d+v"},
        "independent_earlier_acceptance_bounds":acceptance,
        "replacement":{"scale":str(SCALE),"norm":str(replacement_norm),"minimum_margin":str(min(margins)),"strict":True},
        "production_sha256":hashlib.sha256(PRODUCTION.read_bytes()).hexdigest(),"conclusion":{"accepted_cutter_independently_replayed":True,"prefix_infinite_volume":True},
    }


def render(p):
    return f"""# Independent hostile replay: silver DualExtremal H362

**PASS.** This owner-import-free replay reconstructs all {p['control_count']} controls,
the six nonempty levels, the first nonpositive augmented type, and its explicit
lattice realization `u=2d+v`. The old ray is cut, while the `12/25` replacement
has norm `{p['replacement']['norm']}` and strict minimum margin
`{p['replacement']['minimum_margin']}` through the bound. All 17 earlier
cross-Gram bounds were independently re-enumerated; their only two positive
relaxations contain explicit forbidden vectors of norm 6 and 2.
"""


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--verify",action="store_true");args=ap.parse_args();p=audit();a=(json.dumps(p,indent=2,sort_keys=True)+"\n").encode();r=render(p).encode()
    if args.verify:
        if OUTPUT.read_bytes()!=a or REPORT.read_bytes()!=r:raise SystemExit("frozen H362 hostile replay differs")
    else:
        OUTPUT.parent.mkdir(parents=True,exist_ok=True);REPORT.parent.mkdir(parents=True,exist_ok=True);OUTPUT.write_bytes(a);REPORT.write_bytes(r)
    print("PASS",hashlib.sha256(a).hexdigest())


if __name__=="__main__":main()
