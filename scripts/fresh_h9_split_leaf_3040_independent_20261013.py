#!/usr/bin/env python3
"""Independent stdlib exact audit of the literal 3040-cap reconstruction.

Own Fraction arithmetic; no repository producer/arithmetic imports. Checks
seed inertia, every norm/covector/seed gate, and exact profile-root matching.
No cap-cap compatibility or finite-volume conclusion is made.
"""
import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_scout.json'
CERT=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
OUT=ROOT/'results/fresh_h9_split_leaf_3040_independent_20261013.json'
Z=(F(0),F(0));O=(F(1),F(0));PHI=(F(0),F(1))


def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def sub(x,y):return add(x,neg(y))
def sc(x,c):return(x[0]*c,x[1]*c)
def mul(x,y):
    if x==Z or y==Z:return Z
    a,b=x;c,d=y
    return(a*c+b*d,a*d+b*c+b*d)
def inv(x):
    y=(x[0]+x[1],-x[1]);n=mul(x,y)
    assert n[1]==0 and n[0]!=0
    return sc(y,1/n[0])
def sgn(x):return(x>0)-(x<0)
def sign(x):
    p,q=2*x[0]+x[1],x[1];sp,sq=sgn(p),sgn(q)
    if not sp:return sq
    if not sq or sp==sq:return sp
    return sp*sgn(p*p-5*q*q)
def total(xs):
    r=Z
    for x in xs:r=add(r,x)
    return r
def dot(x,y):return total(mul(a,b)for a,b in zip(x,y))
def mv(g,x):return tuple(dot(r,x)for r in g)
def parse(x):return(F(x['a']),F(x['b']))
def encode(x):return{'a':str(x[0]),'b':str(x[1])}


class Quadratic:
    def __init__(self,d):
        self.d=d;assert sign(d)>=0
    def mul(self,x,y):
        a,b=x;c,d=y
        return(add(mul(a,c),mul(self.d,mul(b,d))),add(mul(a,d),mul(b,c)))
    def add(self,x,y):return(add(x[0],y[0]),add(x[1],y[1]))
    def sub(self,x,y):return(sub(x[0],y[0]),sub(x[1],y[1]))
    def sign(self,x):
        a,b=x;sa,sb=sign(a),sign(b)
        if self.d==Z:return sa
        if not sa:return sb
        if not sb or sa==sb:return sa
        return sa*sign(sub(mul(a,a),mul(self.d,mul(b,b))))
    def total(self,xs):
        value=(Z,Z)
        for x in xs:value=self.add(value,x)
        return value
    def dot(self,x,y):return self.total(self.mul(a,b)for a,b in zip(x,y))


def source_pattern(profile):
    assert len(profile)==2 and all(len(r)==8 for r in profile)
    fixed=[];linear=[]
    for raw in sum(profile,[]):
        value=(F(str(raw)),F(0));eps=(F(1,10**8),F(0))
        close=[sign(add(sub(value,c),eps))>0 and sign(sub(sub(value,c),eps))<0 for c in(Z,O,PHI)]
        assert sum(close)<=1
        if any(close):fixed.append((Z,O,PHI)[close.index(True)]);linear.append(Z)
        else:
            assert sign(add(sub(value,sc(O,2)),eps))>=0
            fixed.append(Z);linear.append(O)
    return tuple(fixed),tuple(linear)


def invert_minor(g):
    # Independent exact elimination on a nonsingular ten-wall minor. The
    # complementary six coefficients are zero: a different gauge from the
    # producer's symmetric formula, with the same geometric normal.
    ids=list(range(8))+[8,12]
    matrix=[[g[i][j]for j in ids]+[O if i==j else Z for j in ids]for i in ids]
    for k in range(10):
        pivot=next(i for i in range(k,10)if matrix[i][k]!=Z)
        matrix[k],matrix[pivot]=matrix[pivot],matrix[k]
        c=inv(matrix[k][k]);matrix[k]=[mul(x,c)for x in matrix[k]]
        for i in range(10):
            if i==k:continue
            c=matrix[i][k]
            if c!=Z:matrix[i]=[sub(x,mul(c,y))for x,y in zip(matrix[i],matrix[k])]
    return ids,[row[10:]for row in matrix]


def check_cap(raw,profile,g,minor):
    fixed,linear=source_pattern(profile);ids,inverse=minor
    def solve(values):
        out=[Z]*16
        for i,j in enumerate(ids):out[j]=neg(dot(inverse[i],[values[k]for k in ids]))
        assert mv(g,out)==tuple(map(neg,values))
        return tuple(out)
    r0=solve(fixed);r1=solve(linear)
    g0=tuple(map(neg,fixed));g1=tuple(map(neg,linear))
    A=dot(r1,g1);B=sc(dot(r0,g1),2);C=sub(dot(r0,g0),sc(O,2))
    assert list(map(parse,raw['norm_polynomial']))==[A,B,C]
    d=tuple(map(parse,raw['dotted_parameter']));D=parse(raw['radicand']);field=Quadratic(D)
    if not any(v!=Z for v in linear):
        assert A==B==C==D==Z and d==(Z,Z)and raw['branch']==0
    elif A==Z:
        assert B!=Z and D==Z and d==(neg(mul(C,inv(B))),Z)and raw['branch']==0
    else:
        assert D==sub(mul(B,B),sc(mul(A,C),4))
        assert raw['branch']in(-1,1)
        assert d==(neg(mul(B,inv(sc(A,2)))),mul(sc(O,raw['branch']),inv(sc(A,2))))
    roots=field.add(field.add(field.mul((A,Z),field.mul(d,d)),field.mul((B,Z),d)),(C,Z))
    assert field.sign(roots)==0
    if any(v!=Z for v in linear):assert field.sign(field.sub(d,(sc(O,2),Z)))>=0
    normal=tuple(tuple(map(parse,x))for x in raw['normal'])
    cov=tuple(tuple(map(parse,x))for x in raw['covector'])
    assert len(normal)==len(cov)==16
    # Direct multiplication of the archived coefficients, not the reconstructed gauge.
    direct=tuple((dot(row,[x[0]for x in normal]),dot(row,[x[1]for x in normal]))for row in g)
    assert direct==cov
    square=field.dot(normal,cov)
    assert square==(sc(O,2),Z)
    center=field.total(cov);assert field.sign(center)<0
    counts=Counter()
    for i,v in enumerate(cov):
        expected=(neg(add(fixed[i],mul(linear[i],d[0]))),neg(mul(linear[i],d[1])))
        assert v==expected
        assert field.sign(v)<=0
        if linear[i]==O:
            gap=field.sign(field.add(v,(sc(O,2),Z)))
            assert gap<=0
            counts['parallel'if gap==0 else'strictly_dotted']+=1
        else:
            assert v[1]==Z and v[0]in(Z,neg(O),neg(PHI))
            counts[{Z:'right',neg(O):'pi_over_3',neg(PHI):'pi_over_5'}[v[0]]]+=1
        # Rational interval comparison, replacing the producer's floating root selection.
        target=F(str(sum(profile,[])[i]));tol=F(1,10**6)
        assert field.sign(field.add(v,((target+tol,F(0)),Z)))>0
        assert field.sign(field.add(v,((target-tol,F(0)),Z)))<0
    return fixed,linear,d,D,counts


def run():
    source=json.loads(SOURCE.read_text());data=json.loads(CERT.read_text())
    assert data['source']==str(SOURCE.relative_to(ROOT))
    assert data['source_sha256']==hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    a=(F(1,2),F(1));b=neg(mul(a,inv(add(O,sc(a,12)))))
    assert parse(data['a'])==a and parse(data['b'])==b==(-F(5,22),F(1,11))
    eps=[1]*4+[-1]*4
    g=tuple(tuple(sc(O,2)if i==j else Z if i//8==j//8 else neg(add(add(a,sc(b,eps[i%8]*eps[j%8])),sc(O,2 if i%8==j%8 else 0)))for j in range(16))for i in range(16))
    assert g==tuple(tuple(map(parse,row))for row in data['gram'])
    assert sign(a)>0 and sign(b)<0 and sign(add(b,sc(O,F(1,2))))>0
    assert sign(sub(add(a,b),sc(O,2)))>0
    assert data['center_coefficients']==['1']*16
    assert mv(g,(O,)*16)==(sc(a,-8),)*16
    square=dot((O,)*16,mv(g,(O,)*16));assert square==parse(data['center_square'])==sc(a,-128)
    # Explicit congruence decomposition, checked against every Gram entry.
    S=[1]*16;E=eps+eps;D=[1]*8+[-1]*8;Fv=eps+[-e for e in eps]
    kernel=[];positive=[]
    for start in(0,4):
        for j in range(start+1,start+4):
            k=[0]*16;p=[0]*16
            k[j]=k[j+8]=1;k[start]=k[start+8]=-1
            p[j]=p[start+8]=1;p[start]=p[j+8]=-1
            kernel.append(k);positive.append(p)
    basis=[S,E,D,Fv]+positive+kernel
    gram_basis=[[dot([sc(O,x)for x in v],mv(g,[sc(O,x)for x in w]))for w in basis]for v in basis]
    diagonal=[sc(a,-128),sc(b,-128),add(sc(O,64),sc(a,128)),add(sc(O,64),sc(b,128))]
    for i in range(16):
        for j in range(16):
            expected=diagonal[i]if i==j and i<4 else sc(O,16 if i==j else 8)if 4<=i<10 and 4<=j<10 and(i-4)//3==(j-4)//3 else Z
            assert gram_basis[i][j]==expected
    assert [sign(x)for x in diagonal]==[-1,1,1,1]
    # This displayed sum/difference basis is invertible: in each four-block,
    # three differences plus the block sum form a basis. Positive blocks
    # are 8(I3+J3), giving six further positive directions.
    minor=invert_minor(g)
    assert len(data['caps'])==len(source['all_profiles'])==3040
    totals=Counter();branches=Counter();radicands=set();patterns={};duplicates=[]
    for i,(raw,profile)in enumerate(zip(data['caps'],source['all_profiles'])):
        assert raw['index']==i
        fixed,linear,d,Dq,counts=check_cap(raw,profile,g,minor)
        totals.update(counts);branches[raw['branch']]+=1;radicands.add(Dq)
        key=(fixed,linear)
        # Identical finite/dotted positions give the same norm equation and
        # radicand; exact distinct roots certify mirror distinctness.
        for oldd,oldD,oldindex in patterns.get(key,[]):
            assert oldD==Dq
            if Quadratic(Dq).sign(Quadratic(Dq).sub(d,oldd))==0:
                duplicates.append([oldindex,i]);break
        else:
            patterns.setdefault(key,[]).append((d,Dq,i))
        if i%500==0:print('independently verified normals',i,flush=True)
    assert data['distinct_literal_radicands']==len(radicands)
    assert data['branches']=={str(k):branches[k]for k in(-1,0,1)}
    assert data['status']=='PASS_EXACT_3040_NORMALS'
    forged=deepcopy(data['caps'][192]);forged['normal'][0][0]['a']=str(F(forged['normal'][0][0]['a'])+2)
    try:check_cap(forged,source['all_profiles'][192],g,minor)
    except AssertionError:pass
    else:raise AssertionError('Forged normal accepted')
    forged=deepcopy(data['caps'][192]);forged['branch']=-forged['branch']
    try:check_cap(forged,source['all_profiles'][192],g,minor)
    except AssertionError:pass
    else:raise AssertionError('Forged branch accepted')
    return{'status':'PASS_INDEPENDENT_3040_NORMALS','standard_library_only':True,
           'certificate':str(CERT.relative_to(ROOT)),'source_sha256':data['source_sha256'],
           'seed_signature':[9,1,6],'seed_center_square':encode(square),
           'exact_norm_two_caps':3040,'exact_distinct_cap_mirrors':3040-len(duplicates),
           'duplicate_profile_pairs':duplicates,
           'cap_seed_pairs':3040*16,'cap_seed_pair_counts':dict(sorted(totals.items())),
           'literal_radicands':len(radicands),'branches':dict(sorted(branches.items())),
           'profile_matching':'Exact rational intervals of radius1/1000000 around every saved decimal pairing; finite-profile classification uses radius1/100000000.',
           'arithmetic':'Independent Fraction implementation; covector norm reconstructed in a different ten-wall gauge using exact Gaussian elimination.',
           'forged_normal_rejected':True,'forged_branch_rejected':True,
           'scope':'Exact literal3040 normal reconstruction, seed Coxeter relations, common strict timelike center, and distinct mirrors only. No cap-cap compatibility, completeness of all caps, compactness, or finite-volume claim.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run()
    if args.verify:assert json.loads(OUT.read_text())==json.loads(json.dumps(result))
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
