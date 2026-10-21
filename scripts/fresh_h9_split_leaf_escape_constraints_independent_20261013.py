#!/usr/bin/env python3
"""Independent exact cutter-clause replay; no search or project imports.

Instead of interval signs, compare A^2 and B^2 D in Q(phi) to decide every
sign A+B sqrt(D) exactly. The archived clauses may include extra cutters,
but may not omit a single true cutter. No Boolean satisfiability claim.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NORMALS = ROOT / "data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json"
CLAUSES = ROOT / "data/enumerations/fresh_h9_split_leaf_escape_constraints_20261013.json"
OUT = ROOT / "data/enumerations/fresh_h9_split_leaf_escape_constraints_independent_20261013.json"
ZERO = (F(0), F(0))
ONE = (F(1), F(0))

def sign_int(n):
    return (n > 0) - (n < 0)

def sign(a):
    # a0+a1*phi = ((2a0+a1)+a1*sqrt(5))/2.
    r, s = 2*a[0]+a[1], a[1]
    if r == 0:
        return sign_int(s)
    if s == 0 or sign_int(r) == sign_int(s):
        return sign_int(r)
    return sign_int(r)*sign_int(r*r-5*s*s)

def parse(a):
    return F(a["a"]), F(a["b"])

def add(a,b):
    return a[0]+b[0], a[1]+b[1]

def neg(a):
    return -a[0], -a[1]

def sub(a,b):
    return add(a,neg(b))

def mul(a,b):
    return a[0]*b[0]+a[1]*b[1], a[0]*b[1]+a[1]*b[0]+a[1]*b[1]

def div(a,b):
    den=b[0]*b[0]+b[0]*b[1]-b[1]*b[1]
    assert den
    return mul(a,((b[0]+b[1])/den,-b[1]/den))

def inertia(matrix):
    m=[row[:] for row in matrix]
    out=[0,0,0]
    while m:
        n=len(m)
        pivot=next((i for i in range(n) if m[i][i]!=ZERO),None)
        if pivot is None:
            edge=next(((i,j) for i in range(n) for j in range(i) if m[i][j]!=ZERO),None)
            if edge is None:
                out[2]+=n
                break
            i,j=edge
            # Replace e_i by e_i+e_j by an invertible congruence.
            old=[row[:] for row in m]
            for k in range(n):
                m[i][k]=add(old[i][k],old[j][k])
                m[k][i]=add(old[k][i],old[k][j])
            m[i][i]=add(add(old[i][i],old[j][j]),add(old[i][j],old[j][i]))
            pivot=i
        p=m[pivot][pivot]
        out[0 if sign(p)>0 else 1]+=1
        keep=[i for i in range(n) if i!=pivot]
        m=[[sub(m[i][j],div(mul(m[i][pivot],m[pivot][j]),p)) for j in keep] for i in keep]
    return out

def prepared_sign(d0,d1,d):
    assert sign(d)>=0
    scale=math.lcm(d0[0].denominator,d0[1].denominator)
    aa,bb=int(d0[0]*scale),int(d0[1]*scale)
    # A=(scale*u+t*aa)+(scale*v+t*bb)*phi; B=scale*t*d1.
    c=mul(mul(d1,d1),d)
    c=(c[0]*scale*scale,c[1]*scale*scale)
    den=math.lcm(c[0].denominator,c[1].denominator)
    cc,dd=int(c[0]*den),int(c[1]*den)
    radical_sign=sign(d1) if sign(d)>0 else 0
    def evaluate(u,v,t):
        a,b=scale*u+t*aa,scale*v+t*bb
        sa=sign((a,b))
        sb=sign_int(t)*radical_sign
        if not sa:return sb
        if not sb or sa==sb:return sa
        # No field-independence assumption: equality of these squares gives zero.
        delta=(den*(a*a+b*b)-t*t*cc,den*(2*a*b+b*b)-t*t*dd)
        return sa*sign(delta)
    return evaluate

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build(clauses=CLAUSES):
    normal_bytes=NORMALS.read_bytes()
    clause_bytes=clauses.read_bytes()
    normal=json.loads(normal_bytes)
    data=json.loads(clause_bytes)
    assert data["normal_source"]==str(NORMALS.relative_to(ROOT))
    assert data["normal_sha256"]==digest(NORMALS)
    source=ROOT/normal["source"]
    assert digest(source)==normal["source_sha256"]
    G=[[parse(x) for x in row] for row in normal["gram"]]
    assert len(G)==16 and all(len(row)==16 for row in G)
    A=[[0]*16 for _ in range(16)]
    B=[[0]*16 for _ in range(16)]
    for i in range(16):
        for j in range(16):
            if i==j:
                a,b=22,0
            elif i//8==j//8:
                a,b=0,0
            else:
                same=(i%8)//4==(j%8)//4
                a=-(3 if same else 8)-(22 if i%8==j%8 else 0)
                b=-(12 if same else 10)
            assert G[i][j]==(F(a,11),F(b,11))
            A[i][j],B[i][j]=a,b
    signature=inertia(G)
    assert signature==[9,1,6]
    center=(F(sum(map(sum,A)),11),F(sum(map(sum,B)),11))
    assert center==parse(normal["center_square"])==(-64,-128)
    assert all(F(x)==1 for x in normal["center_coefficients"])
    assert sign(center)<0
    patterns=[]
    for index,cap in enumerate(normal["caps"]):
        assert cap["index"]==index
        d0,d1=map(parse,cap["dotted_parameter"])
        d=parse(cap["radicand"])
        evaluator=prepared_sign(d0,d1,d)
        p=[[parse(a) for a in v] for v in cap["covector"]]
        n=[[parse(a) for a in v] for v in cap["normal"]]
        assert len(p)==len(n)==16 and all(len(v)==2 for v in p+n)
        for i in range(16):
            for component in range(2):
                total=ZERO
                for j in range(16):
                    total=add(total,mul(G[i][j],n[j][component]))
                assert total==p[i][component],("normal/covector mismatch",index,i)
        norm0=ZERO
        norm1=ZERO
        center0=ZERO
        center1=ZERO
        for ni,pi in zip(n,p,strict=True):
            norm0=add(norm0,add(mul(ni[0],pi[0]),mul(d,mul(ni[1],pi[1]))))
            norm1=add(norm1,add(mul(ni[0],pi[1]),mul(ni[1],pi[0])))
            center0=add(center0,pi[0])
            center1=add(center1,pi[1])
        assert prepared_sign(sub(norm0,(F(2),F(0))),norm1,d)(0,0,1)==0
        assert prepared_sign(center0,center1,d)(0,0,1)<0
        masks=[[],[],[]]
        for j,(x,y) in enumerate(p):
            if x==ZERO and y==ZERO:continue
            if y==ZERO and x==(-1,0):masks[0].append(j)
            elif y==ZERO and x==(0,-1):masks[1].append(j)
            else:
                assert (x,y)==(neg(d0),neg(d1))
                masks[2].append(j)
        patterns.append((masks,evaluator))
    counters={"clauses":len(data["constraints"]),"caps":len(patterns),
              "exact_pairing_signs":0,"true_cutters":0,"allowed_extra_cutters":0,
              "strict_seed_checks":0}
    first_true=None
    for clause in data["constraints"]:
        w=clause["integer_coefficients"]
        assert len(w)==16 and all(type(x) is int for x in w)
        den=clause["denominator"]
        assert type(den) is int and den>0
        aw=[sum(x*y for x,y in zip(row,w,strict=True)) for row in A]
        bw=[sum(x*y for x,y in zip(row,w,strict=True)) for row in B]
        q=(sum(x*y for x,y in zip(w,aw,strict=True)),sum(x*y for x,y in zip(w,bw,strict=True)))
        c=(sum(aw),sum(bw))
        assert sign(q)>0 and sign(c)<0
        assert list(map(str,q))==clause["norm_numerator"]
        assert clause["norm_denominator"]==11*den*den
        assert list(map(str,c))==clause["center_pairing_numerator"]
        assert clause["center_pairing_denominator"]==11*den
        assert all(sign((a,b))<0 for a,b in zip(aw,bw,strict=True))
        counters["strict_seed_checks"]+=16
        true=0
        for k,(masks,evaluate) in enumerate(patterns):
            u,v,t=(sum(w[j] for j in mask) for mask in masks)
            if evaluate(u,v,t)<0:true|=1<<k
        allowed=int(clause["possible_cutters_hex"],16)
        assert 0<=allowed<1<<len(patterns)
        assert true & ~allowed == 0,"A necessary clause omits an actual cutter"
        if clause["interval_ambiguous_cap_pairs"]==0:
            assert true==allowed
        counters["exact_pairing_signs"]+=len(patterns)
        counters["true_cutters"]+=true.bit_count()
        counters["allowed_extra_cutters"]+=(allowed & ~true).bit_count()
        if true and first_true is None:first_true=(true,allowed)
    assert first_true is not None
    true,allowed=first_true
    forged=allowed & ~(true & -true)
    assert true & ~forged != 0
    assert sign((0,0))==0
    # Square/dependent-radical edge cases exercise exact cancellation.
    assert prepared_sign((F(0),F(0)),ONE,ONE)(-1,0,1)==0
    assert prepared_sign(ZERO,ONE,ZERO)(-1,0,1)==-1
    assert NORMALS.read_bytes()==normal_bytes and clauses.read_bytes()==clause_bytes,"Inputs changed during this audit"
    return {"status":"PASS_INDEPENDENT_EXACT_NECESSARY_CUTTERS",
            "scope":"Every archived escape is strict future spacelike, every cap has norm two and a common strict timelike center, every cap covector equals its literal normal's Gram pairing, and every actual cutter is included. This is not an UNSAT or finite-volume impossibility proof by itself.",
            "normal_source":str(NORMALS.relative_to(ROOT)),"normal_sha256":digest(NORMALS),
            "clause_source":str(clauses.relative_to(ROOT)),"clause_sha256":digest(clauses),
            "upstream_source":str(source.relative_to(ROOT)),"upstream_sha256":digest(source),
            "seed_gram_inertia":signature,"center_square":{"a":str(center[0]),"b":str(center[1])},
            "counts":counters,
            "mutations":{"deleted_true_cutter_rejected":True,"zero_escape_rejected":True,
                         "dependent_and_zero_radical_tests":True}}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--verify",action="store_true")
    parser.add_argument("--source",type=Path,default=CLAUSES)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    if not args.source.is_absolute():args.source=ROOT/args.source
    args.source=args.source.resolve()
    if args.output is None:
        if args.source!=CLAUSES.resolve():
            parser.error("An alternate --source requires --output to preserve the previous audit")
        args.output=OUT
    elif not args.output.is_absolute():args.output=ROOT/args.output
    result=build(args.source)
    if args.verify:
        assert result==json.loads(args.output.read_text()),"The audited sources or exact clauses changed; regenerate the audit."
    else:
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print({"status":result["status"],**result["counts"],"verify":args.verify},flush=True)

if __name__=="__main__":
    main()
