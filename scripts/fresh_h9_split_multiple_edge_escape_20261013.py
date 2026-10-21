#!/usr/bin/env python3
"""Rational interval reconstruction of every no-removal cap and a common escape.

Each cap is defined by its nine exact prescribed pairings and its norm-two
quadratic, with a unique branch matched to the numerical proposal. This does
not certify its other Coxeter pairings; those are unnecessary for an escape
that survives even their whole intersection. No floating value decides a sign.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
from pathlib import Path
import time

import fresh_h4_blocks_independent_verify as K

ROOT=Path(__file__).resolve().parents[1]
SCOUT='data/enumerations/fresh_h9_split_endpoint_walk_optional_anchor_20261013_scout.json'
OLD='data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
POINT='data/enumerations/fresh_h9_split_multiple_edge_caps_20261013_escape_scout.json'
OUT=ROOT/'results/fresh_h9_split_multiple_edge_escape_20261013.json'
BASIS=(0,1,2,3,4,5,6,7,8,12)
S=10**70
Z=(0,0)
O=(S,S)


def ceildiv(a,b):return -((-a)//b)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[1],-x[0]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
    p=(x[0]*y[0],x[0]*y[1],x[1]*y[0],x[1]*y[1])
    return min(p)//S,ceildiv(max(p),S)
def rat(x):
    x=Q(x);return x.numerator*S//x.denominator,ceildiv(x.numerator*S,x.denominator)
def div(x,y):
    assert y[1]<0 or y[0]>0,('zero-containing divisor',y)
    r=(S*S//y[1],ceildiv(S*S,y[0]))
    return mul(x,r)
def sqrt(x):
    assert x[0]>=0,('negative radicand',x)
    return isqrt(x[0]*S),isqrt(x[1]*S)+1
def total(xs):
    r=Z
    for x in xs:r=add(r,x)
    return r
def dot(a,b):return total(mul(x,y)for x,y in zip(a,b,strict=True))
def mv(a,b):return [dot(row,b)for row in a]
def pair(a,g,b):return dot(a,mv(g,b))
def enc(x):return [str(x[0]),str(x[1])]
PHI=div(add(O,sqrt(rat(5))),rat(2))
def golden(x):return add(rat(x[0]),mul(rat(x[1]),PHI))
def extension(x,d):return add(golden(x[0]),mul(golden(x[1]),sqrt(golden(d))))
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def solve_line(rows):
    """Enclose the exact RREF of a rank-nine 9x10 matrix with full pivoting."""
    a=[list(row)+[O if i==j else Z for j in range(9)]for i,row in enumerate(rows)]
    columns=list(range(10));pivots=[]
    for k in range(9):
        options=[(min(abs(a[i][j][0]),abs(a[i][j][1])),i,j)
                 for i in range(k,9)for j in columns if a[i][j][0]>0 or a[i][j][1]<0]
        assert options,'Cannot certify rank nine'
        _,i,j=max(options);a[k],a[i]=a[i],a[k];p=a[k][j]
        a[k]=[div(x,p)for x in a[k]];a[k][j]=O
        for i in range(9):
            if i==k:continue
            p=a[i][j];a[i]=[sub(x,mul(p,y))for x,y in zip(a[i],a[k],strict=True)];a[i][j]=Z
        columns.remove(j);pivots.append(j)
    free,=columns
    direction=[Z]*10;direction[free]=O
    transform=[[Z]*9 for _ in range(10)]
    for i,j in enumerate(pivots):
        direction[j]=neg(a[i][free]);transform[j]=a[i][10:]
    return free,transform,direction


def exact_old():
    scout=read(SCOUT);snap=scout['source_snapshot'];old=read(OLD)
    g16=K.matrix(old['gram'])
    a=(Q(1,2),Q(1));b=(Q(-5,22),Q(1,11));eps=(1,1,1,1,-1,-1,-1,-1)
    for i in range(16):
        for j in range(16):
            want=K.scale(K.O,2)if i==j else K.Z if i//8==j//8 else K.neg(
                K.add(K.add(a,K.scale(b,eps[i%8]*eps[j%8])),K.scale(K.O,2 if i%8==j%8 else 0)))
            assert g16[i][j]==want
    g=tuple(tuple(g16[i][j]for j in BASIS)for i in BASIS);inv=K.inverse(g)
    reduced=[list(row)for row in g];signs=[]
    for i in range(10):
        pivot=reduced[i][i];signs.append(K.sign(pivot));assert signs[-1]
        for j in range(i+1,10):
            for k in range(i+1,10):
                reduced[j][k]=K.sub(reduced[j][k],K.div(K.mul(reduced[j][i],reduced[i][k]),pivot))
    assert signs.count(1)==9 and signs.count(-1)==1
    gi=[[golden(x)for x in row]for row in inv]
    profiles=[];banks={OLD:old}
    for i in range(16):
        p=tuple(g16[j][i]for j in BASIS);c=K.mv(inv,p)
        profiles.append([golden(x)for x in c])
    for filename,index in snap['refs']:
        path='data/enumerations/'+filename
        if path not in banks:banks[path]=read(path)
        cap=banks[path]['caps'][index];d=K.parse(cap['radicand'])
        p=[extension(tuple(map(K.parse,cap['covector'][j])),d)for j in BASIS]
        profiles.append(mv(gi,p))
    assert len(profiles)==72
    h=tuple(K.neg(K.total(g16[j]))for j in BASIS)
    assert K.sign(K.pair(h,inv,h))<0
    assert all(dot(c,list(map(golden,h)))[0]>0 for c in profiles)
    return scout,snap,inv,gi,h,profiles,banks


def proposed_point(snapshot,inv,h,profiles):
    # Floats propose rational coordinates only. Every geometric sign below
    # is recomputed with exact Q(phi) or outward-rounded integer intervals.
    import numpy as np
    v=np.array(read(POINT)['selected_point'])*(1-1e-6)
    physical=np.r_[1.,v];N=np.array(snapshot['N'])[list(BASIS)]
    y=-(N*(np.r_[-1.,np.ones(9)]))@physical
    yy=tuple((Q(round(float(x)*10**12),10**12),Q(0))for x in y)
    q=K.pair(yy,inv,yy);future=K.pair(yy,inv,h)
    assert K.sign(q)>0 and K.sign(future)<0
    yi=list(map(golden,yy));old_pairs=[dot(c,yi)for c in profiles]
    assert min(x[0]for x in old_pairs)>0
    return yy,q,future,old_pairs


def run(saved_point=None):
    start=time.monotonic();scout,snap,inv,gi,h,old,banks=exact_old()
    if saved_point is None:yy,q,future,old_pairs=proposed_point(snap,inv,h,old)
    else:
        yy=tuple(map(K.parse,saved_point));q=K.pair(yy,inv,yy);future=K.pair(yy,inv,h)
        assert K.sign(q)>0 and K.sign(future)<0
        old_pairs=[dot(c,list(map(golden,yy)))for c in old]
        assert min(x[0]for x in old_pairs)>0
    Iy=list(map(golden,K.mv(inv,yy)));Ih=list(map(golden,K.mv(inv,h)));cache={};records=[]
    for index,proposal in enumerate(scout['candidates']):
        if proposal['conflicts']:continue
        key=tuple(proposal['edge']+[proposal['anchor']])
        if key not in cache:
            free,transform,v=solve_line([old[i]for i in key]);cache[key]=(free,transform,v,pair(v,gi,v))
        free,transform,v,A=cache[key]
        rhs=[Z]*8
        for j,m in proposal['finite_profile']:rhs[j]=neg(O if m==3 else PHI)
        rhs.append({2:Z,3:neg(O),5:neg(PHI)}[proposal['anchor_label']])
        u=mv(transform,rhs)
        B=mul(rat(2),pair(u,gi,v));C=sub(pair(u,gi,u),rat(2))
        # The exact quadratic leading coefficient may vanish; none of the
        # accepted fixed-bank entries is allowed to rely on guessing that.
        assert A[0]>0 or A[1]<0,('undecided leading coefficient',index,A)
        D=sub(mul(B,B),mul(rat(4),mul(A,C)));assert D[0]>0,('no positive discriminant',index)
        target=[Q(str(proposal['seed_covector'][j]))for j in BASIS]
        matches=[]
        for branch in(-1,1):
            t=div(add(neg(B),mul(rat(branch),sqrt(D))),mul(rat(2),A))
            p=[add(x,mul(t,z))for x,z in zip(u,v,strict=True)]
            tol=Q(1,10**6)
            # A full interval must lie inside the matching rational box.
            if all(Q(lo,S)>w-tol and Q(hi,S)<w+tol for (lo,hi),w in zip(p,target,strict=True)):
                matches.append((branch,p))
            else:
                assert any(Q(hi,S)<w-tol or Q(lo,S)>w+tol for (lo,hi),w in zip(p,target,strict=True)),(
                    'Other branch not separated from the proposal box',index,branch)
        assert len(matches)==1,('nonunique or unmatched exact branch',index,len(matches))
        branch,p=matches[0];escape=dot(p,Iy);assert escape[0]>0,('escape not certified',index,escape)
        interior=dot(p,Ih);assert interior[0]>0,('center not certified',index,interior)
        # Exact norm two follows from the quadratic definition, not a small
        # numerical residual. The interval is an additional consistency gate.
        norm=pair(p,gi,p);assert norm[0]<=2*S<=norm[1]
        records.append({'proposal':index,'free_coordinate':free,'branch':branch,
                        'discriminant':enc(D),'escape_pairing':enc(escape),'center_pairing':enc(interior),
                        'norm_enclosure':enc(norm),'seed_covector_enclosures':list(map(enc,p))})
        if len(records)%400==0:print({'reconstructed':len(records),'seconds':round(time.monotonic()-start,1)},flush=True)
    assert len(records)==3857
    return {'status':'PASS_RIGOROUS_COMMON_ESCAPE_FOR_3857_DEFINED_CAPS',
            'source_sha256':{p:sha(p)for p in sorted({SCOUT,POINT,*banks})},
            'interval_scale':str(S),'phi_interval':enc(PHI),'independent_basis':list(BASIS),
            'gram_inertia':[9,1,0],'strict_center':list(map(K.enc,h)),
            'point':list(map(K.enc,yy)),'point_norm':K.enc(q),'future_center_pairing':K.enc(future),
            'old_wall_pairings':list(map(enc,old_pairs)),
            'cap_count':len(records),'distinct_linear_systems':len(cache),
            'minimum_new_pairing_lower_bound':str(min(int(r['escape_pairing'][0])for r in records)),
            'caps':records,'seconds':round(time.monotonic()-start,3),
            'scope':'Each archived no-removal proposal is interpreted as the unique norm-two root branch in the certified rational neighborhood satisfying its nine specified exact pairing equations. One future strict spacelike point survives every such root and all72 old walls. Thus no subset of this defined3857-cap bank gives finite volume with the retained old walls, regardless of pair compatibility. Other roots, other branches, omitted anchors/labels, and further edge generations remain open. This does not certify every proposed Coxeter angle.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    saved=json.loads(OUT.read_text())if args.verify else None
    result=run(saved['point']if saved else None)
    if saved:
        saved.pop('seconds');result.pop('seconds');assert saved==result
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k in ['status','cap_count','distinct_linear_systems','seconds']}),flush=True)
