#!/usr/bin/env python3
"""Exact necessary endpoint clauses for the two orthogonal eight-wall edges.

No numerical or search imports. Every endpoint/cap sign is determined exactly
in a possibly dependent biquadratic extension of Q(phi).
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json"
OUT=ROOT/"data/enumerations/fresh_h9_split_leaf_edge_endpoints_20261013.json"
Z=(F(0),F(0));O=(F(1),F(0))
def k(n):return F(n),F(0)
def dec(x):return F(x["a"]),F(x["b"])
def enc(x):return {"a":str(x[0]),"b":str(x[1])}
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def scale(x,n):return x[0]*n,x[1]*n
def divide(x,y):
    den=y[0]*y[0]+y[0]*y[1]-y[1]*y[1]
    assert den
    return mul(x,((y[0]+y[1])/den,-y[1]/den))
def sgn(x):
    a,b=2*x[0]+x[1],x[1]
    sa=(a>0)-(a<0);sb=(b>0)-(b<0)
    if not sa:return sb
    if not sb or sa==sb:return sa
    t=a*a-5*b*b
    return sa*((t>0)-(t<0))
def dot(x,y):
    assert len(x)==len(y)
    ans=Z
    for a,b in zip(x,y):ans=add(ans,mul(a,b))
    return ans
def mv(g,x):return [dot(row,x)for row in g]
def pair(x,g,y):return dot(x,mv(g,y))

class E:
    def __init__(self,d):
        assert sgn(d)>=0
        self.d=d
    def lift(self,x):return x,Z
    def add(self,x,y):return add(x[0],y[0]),add(x[1],y[1])
    def neg(self,x):return neg(x[0]),neg(x[1])
    def sub(self,x,y):return self.add(x,self.neg(y))
    def mul(self,x,y):
        return add(mul(x[0],y[0]),mul(self.d,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))
    def sign(self,x):
        a,b=sgn(x[0]),sgn(x[1]) if self.d!=Z else 0
        if not a:return b
        if not b or a==b:return a
        return a*sgn(sub(mul(x[0],x[0]),mul(self.d,mul(x[1],x[1]))))
    def dot(self,x,y):
        assert len(x)==len(y)
        ans=self.lift(Z)
        for a,b in zip(x,y):ans=self.add(ans,self.mul(a,b))
        return ans
    def outer_sign(self,A,B,d):
        # Sign(A+B sqrt(d)); A,B lie in this E, while d is in Q(phi).
        assert sgn(d)>0
        a,b=self.sign(A),self.sign(B)
        if not a:return b
        if not b or a==b:return a
        return a*self.sign(self.sub(self.mul(A,A),self.mul(self.lift(d),self.mul(B,B))))

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def build():
    raw=SOURCE.read_bytes();data=json.loads(raw)
    upstream=ROOT/data["source"];assert digest(upstream)==data["source_sha256"]
    a=(F(1,2),F(1));b=(F(-5,22),F(1,11))
    assert b==neg(divide(a,add(scale(a,12),O)))
    assert sgn(a)>0 and sgn(b)<0 and sgn(add(b,k(F(1,2))))>0
    eps=[1]*4+[-1]*4
    g=[[dec(x)for x in row]for row in data["gram"]]
    assert len(g)==16 and all(len(row)==16 for row in g)
    for i in range(16):
        for j in range(16):
            expected=k(2)if i==j else Z
            if i//8!=j//8:
                expected=neg(add(add(k(2 if i%8==j%8 else 0),a),scale(b,eps[i%8]*eps[j%8])))
            assert g[i][j]==expected
    S=[O]*16
    center_square=pair(S,g,S)
    assert center_square==scale(a,-128) and sgn(center_square)<0
    assert all(sgn(z)<0 for z in mv(g,S))
    d=neg(divide(mul(a,add(O,scale(a,2))),mul(b,add(O,scale(b,2)))))
    assert sgn(d)>0
    edge_field=E(d)
    rows=[]
    for side in range(2):
        P=[add(O,scale(a,4))if i//8==side else O for i in range(16)]
        R=[scale(add(O,scale(b,4))if i//8==side else O,eps[i%8])for i in range(16)]
        p2=pair(P,g,P);r2=pair(R,g,R);pr=pair(P,g,R)
        assert p2==scale(mul(a,add(O,scale(a,2))),-128) and sgn(p2)<0
        assert r2==scale(mul(b,add(O,scale(b,2))),-128) and sgn(r2)>0
        assert pr==Z and add(p2,mul(d,r2))==Z
        assert pair(R,g,S)==Z and pair(P,g,S)==p2
        for sigma in(-1,1):
            point=[(x,scale(y,sigma))for x,y in zip(P,R)]
            cov=[edge_field.dot([edge_field.lift(z)for z in gr],point)for gr in g]
            signs=[edge_field.sign(v)for v in cov]
            assert signs==[0 if i//8==side else -1 for i in range(16)]
            rows.append({"seed_side":side,"radical_sign":sigma,
                         "endpoint_coefficients":[[enc(v)for v in x]for x in point],
                         "seed_pairings":[[enc(v)for v in x]for x in cov],
                         "seed_pairing_signs":signs,"endpoint_norm":"0",
                         "center_pairing":enc(p2),"projected_center_norm":enc(p2),
                         "spacelike_direction_norm":enc(r2),
                         "nonnegative_caps":[],"zero_caps":[],"positive_caps":[],
                         "cap_pairing_signs":[]})
    caps=data["caps"]
    for index,cap in enumerate(caps):
        assert cap["index"]==index
        f=E(dec(cap["radicand"]))
        normal=[tuple(map(dec,v))for v in cap["normal"]]
        cov=[tuple(map(dec,v))for v in cap["covector"]]
        assert len(normal)==len(cov)==16
        for i,gr in enumerate(g):
            actual=f.dot([f.lift(z)for z in gr],normal)
            assert all(x==y for x,y in zip(actual,cov[i]))
        assert f.sign(f.sub(f.dot(normal,cov),f.lift(k(2))))==0
        assert f.sign(f.dot([f.lift(O)]*16,cov))<0
        assert all(f.sign(z)<=0 for z in cov)
        all_signs=[]
        for row in rows:
            points=[tuple(map(dec,v))for v in row["endpoint_coefficients"]]
            A=f.dot([f.lift(x[0])for x in points],cov)
            B=f.dot([f.lift(x[1])for x in points],cov)
            value=f.outer_sign(A,B,d)
            assert f.sign(A)<0  # Projected center is strict for every cap.
            row["cap_pairing_signs"].append(value)
            if value>=0:row["nonnegative_caps"].append(index)
            if value==0:row["zero_caps"].append(index)
            if value>0:row["positive_caps"].append(index)
            all_signs.append(value)
        assert not(all_signs[0]>=0 and all_signs[1]>=0)
        assert not(all_signs[2]>=0 and all_signs[3]>=0)
    assert SOURCE.read_bytes()==raw
    assert E(O).sign((O,neg(O)))==0
    assert E(Z).sign((neg(O),k(100)))==-1
    assert E(k(2)).outer_sign((Z,O),(neg(O),Z),k(2))==0
    return {"status":"PASS_EXACT_NECESSARY_IDEAL_EDGE_ENDPOINT_CLAUSES",
            "scope":"For any finite-volume acute-compatible completion retaining the sixteen seed facets and selecting caps from this literal bank, each of the four nonnegative-endpoint lists must be hit. Zero pairings must be included. No compatibility or Boolean contradiction is asserted here.",
            "source":str(SOURCE.relative_to(ROOT)),"source_sha256":digest(SOURCE),
            "upstream_source":str(upstream.relative_to(ROOT)),"upstream_sha256":digest(upstream),
            "a":enc(a),"b":enc(b),"endpoint_radicand":enc(d),
            "center_square":enc(center_square),"entry_count":len(caps),
            "endpoint_count":4,"exact_cap_endpoint_signs":4*len(caps),
            "rows":rows,
            "checks":{"all_cap_norms_two":True,"all_literal_normal_covector_identities":True,
                      "strict_common_timelike_center":True,"all_four_endpoints_future_null":True,
                      "exactly_eight_zero_seed_pairings_at_each_endpoint":True,
                      "each_cap_nonnegative_at_at_most_one_end_of_each_edge":True,
                      "dependent_and_zero_radical_tests":True}}

def main():
    p=argparse.ArgumentParser();p.add_argument("--verify",action="store_true");args=p.parse_args()
    result=build()
    if args.verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+"\n")
    print({"status":result["status"],"signs":result["exact_cap_endpoint_signs"],
           "rows":[{"side":r["seed_side"],"sign":r["radical_sign"],"nonnegative":r["nonnegative_caps"],"zero":r["zero_caps"]}for r in result["rows"]],
           "verify":args.verify},flush=True)
if __name__=="__main__":main()
