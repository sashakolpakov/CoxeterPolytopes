#!/usr/bin/env python3
"""Standalone Fraction replay and independent census for the prime-twisted host."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product,combinations,groupby
import json
from math import ceil,isqrt,sqrt
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
Z=(Q(0),)*4;O=(Q(1),Q(0),Q(0),Q(0));T=(Q(0),Q(1),Q(0),Q(0))
def add(x,y):return tuple(a+b for a,b in zip(x,y,strict=True))
def neg(x):return tuple(-a for a in x)
def sub(x,y):return add(x,neg(y))
def scale(x,n):return tuple(a*n for a in x)
def mul(x,y):
    z=[Q(0)]*7
    for i in range(4):
        for j in range(4):z[i+j]+=x[i]*y[j]
    for i in range(6,3,-1):z[i-2]+=4*z[i];z[i-4]-=2*z[i]
    return tuple(z[:4])
def silver_mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def silver_sign(x):
    a,b=x;sa=(a>0)-(a<0);sb=(b>0)-(b<0)
    if not sa:return sb
    if not sb or sa==sb:return sa
    c=a*a-2*b*b;return sa*((c>0)-(c<0))
def sign(x,e=0):
    a,b,c,d=x;sg=1 if e in(0,3)else-1
    u=(a+2*c,sg*c);v=(b+2*d,sg*d)
    su=silver_sign(u);sv=silver_sign(v)*(1 if e in(0,1)else-1)
    if not su:return sv
    if not sv or sv==su:return su
    uu=silver_mul(u,u);vv=silver_mul((Q(2),Q(sg)),silver_mul(v,v))
    return su*silver_sign((uu[0]-vv[0],uu[1]-vv[1]))
@lru_cache(None)
def inv(x):
    cols=[mul(x,tuple(Q(i==j)for i in range(4)))for j in range(4)]
    a=[[cols[j][i]for j in range(4)]+[Q(i==0)]for i in range(4)]
    for i in range(4):
        j=next(j for j in range(i,4)if a[j][i]);a[i],a[j]=a[j],a[i]
        q=a[i][i];a[i]=[t/q for t in a[i]]
        for j in range(4):
            if j!=i:
                q=a[j][i];a[j]=[s-q*t for s,t in zip(a[j],a[i])]
    return tuple(row[4]for row in a)
def divide(x,y):return mul(x,inv(y))
def sq(x):return mul(x,x)
def integral(x):return all(t.denominator==1 for t in x)
def real(x,e=0):
    t=(1 if e<2 else-1)*sqrt(2+(1 if e in(0,3)else-1)*sqrt(2))
    return sum(float(a)*t**i for i,a in enumerate(x))
def dec(x):return tuple(Q(t)for t in x)
def vec(x):return tuple(map(dec,x))
def enc(x):return list(map(str,x))
def dot(x,y):
    out=Z
    for a,b in zip(x,y,strict=True):out=add(out,mul(a,b))
    return out
def mv(g,v):return tuple(dot(row,v)for row in g)
def minverse(g):
    n=len(g);a=[list(row)+[O if i==j else Z for j in range(n)]for i,row in enumerate(g)]
    for i in range(n):
        j=next(j for j in range(i,n)if a[j][i]!=Z);a[i],a[j]=a[j],a[i]
        q=inv(a[i][i]);a[i]=[mul(x,q)for x in a[i]]
        for j in range(n):
            if i!=j:
                q=a[j][i];a[j]=[sub(x,mul(q,y))for x,y in zip(a[j],a[i])]
    return tuple(tuple(row[n:])for row in a)
ALPHA=mul(sub(T,O),add(scale(O,2),T));SQRT2=sub(sq(T),scale(O,2));PI=add(scale(O,2),T)
POOL=(O,PI,sq(T),divide(scale(O,2),PI),scale(O,2))
F4=tuple(tuple(scale(O,2)if i==j else neg(SQRT2)if{i,j}=={1,2}else neg(O)if abs(i-j)==1 else Z for j in range(4))for i in range(4))
FI=minverse(F4);ALLOWED={Z,O,scale(O,2),scale(O,3),add(scale(O,2),SQRT2),add(scale(O,2),T)}

def ceil_real(x,e):
    n=ceil(real(x,e))
    while sign(sub(x,scale(O,n)),e)>0:n+=1
    while sign(sub(x,scale(O,n-1)),e)<=0:n-=1
    return n
@lru_cache(None)
def alphabet(caps):
    if any(sign(c,e)<0 for e,c in enumerate(caps)):return()
    bound=ceil(4*caps[0][0]+8*caps[0][2])if len(set(caps))==1 else sum(ceil_real(c,e)for e,c in enumerate(caps))
    out=[]
    for c in range(-isqrt(bound//8),isqrt(bound//8)+1):
      for d in range(-isqrt((bound-8*c*c)//8),isqrt((bound-8*c*c)//8)+1):
        budget=bound-8*c*c-8*d*d
        for A in range(-isqrt(budget//4),isqrt(budget//4)+1):
          for B in range(-isqrt((budget-4*A*A)//8),isqrt((budget-4*A*A)//8)+1):
            x=(Q(A-2*c),Q(B-3*d),Q(c),Q(d))
            if sign(x)<0:continue
            xx=sq(x)
            if all(sign(sub(cap,xx),e)>=0 for e,cap in enumerate(caps)):out.append(x)
    return tuple(out)
@lru_cache(None)
def block_vectors(energy):
    # Different from the producer's recursive common-alphabet Schur bounds:
    # positive inverse entries bound each dominant covector separately.
    weights=[alphabet((divide(energy,FI[i][i]),*(scale(energy,2)for _ in range(3))))for i in range(4)]
    out=[]
    for w in product(*weights):
        y=mv(FI,w);en=dot(w,y)
        if all(sign(sub(energy,en),e)>=0 for e in range(4)):
            assert all(map(integral,y));out.append((en,y))
    return tuple(out)

class Host:
    def __init__(self,blocks):
        self.blocks=blocks;self.n=4*blocks+1;n=self.n
        self.g=tuple(tuple(neg(ALPHA)if i==j==0 else O if i==j==n else F4[(i-1)%4][(j-1)%4]
            if 1<=i<=4*blocks and 1<=j<=4*blocks and(i-1)//4==(j-1)//4 else Z for j in range(n+1))for i in range(n+1))
    def pair(self,x,y):return dot(x,mv(self.g,y))
    def controls(self,height):
        rows=[(divide(mul(ALPHA,sq(x)),q),q,x)for q in POOL
            for x in alphabet((divide(scale(q,height),ALPHA),*(divide(neg(q),ALPHA)for _ in range(3))))if sign(x)>0]
        from functools import cmp_to_key
        rows.sort(key=cmp_to_key(lambda a,b:next((sign(sub(x,y))for x,y in zip(a,b)if x!=y),0)))
        return rows
    def candidates(self,x,q):
        E=add(q,mul(ALPHA,sq(x)));bs=block_vectors(E)
        ss={sq(z):z for z in alphabet((E,)*4)};out=[];imprimitive=0
        for packet in product(bs,repeat=self.blocks):
            rem=E
            for en,y in packet:rem=sub(rem,en)
            if rem not in ss:continue
            r=(x,*(v for en,y in packet for v in y),ss[rem])
            assert self.pair(r,r)==q
            if all(integral(divide(v,T))for v in r):imprimitive+=1;continue
            if all(integral(divide(scale(v,2),q))for v in mv(self.g,r)):out.append(r)
        return out,imprimitive,len(bs),len(ss)
    def compatible(self,r,walls):
        q=self.pair(r,r)
        for s in walls:
            p=self.pair(r,s)
            if sign(p)>0:return False
            ratio=divide(scale(sq(p),4),mul(q,self.pair(s,s)))
            if sign(sub(ratio,scale(O,4)))<0 and ratio not in ALLOWED:return False
        return True
    def check_geometry(self,data):
        walls=tuple(vec(r)for r in data['all_walls']);inside=vec(data['strict_interior'])
        assert sign(self.pair(inside,inside))<0 and all(sign(self.pair(inside,r))<0 for r in walls)
        for i,r in enumerate(walls):
            assert all(map(integral,r))and not all(integral(divide(x,T))for x in r)
            q=self.pair(r,r);assert q in POOL and all(sign(q,e)>0 for e in range(4))
            assert all(integral(divide(scale(x,2),q))for x in mv(self.g,r))
            assert self.compatible(r,walls[:i])
        geom=data['geometry'];assert geom['status']=='EXACT_STRICT_SPACELIKE_ESCAPE'
        v=vec(geom['vector']);strict=vec(geom['strict_vector'])
        assert self.pair(v,v)==dec(geom['norm'])and sign(self.pair(v,v))>0 and sign(v[0])>0
        assert all(sign(self.pair(v,r))<=0 for r in walls)
        assert self.pair(strict,strict)==dec(geom['strict_norm'])and sign(self.pair(strict,strict))>0 and sign(strict[0])>0
        assert all(sign(self.pair(strict,r))<0 for r in walls)
        assert self.pair(strict,inside)==dec(geom['strict_center_pairing'])
        return {'walls':len(walls),'all_pair_checks':len(walls)*(len(walls)-1)//2,'exact_strict_escape':True,'escape_norm':enc(self.pair(v,v))}

def verify(source):
    data=json.loads(source.read_text());host=Host(data['blocks']);assert tuple(vec(r)for r in data['gram'])==host.g
    assert data['completed']and vec(data['norm_pool'])==POOL
    zero={(q,r)for q in POOL for r in host.candidates(Z,q)[0]}
    assert zero=={(dec(r['norm']),vec(r['vector']))for r in data['zero_dominant_roots']}
    assert len(zero)==2*host.blocks+1
    controls=host.controls(data['height_bound']);saved={(dec(r['height']),dec(r['norm']),dec(r['time'])):r for r in data['levels']}
    assert len(controls)==len(saved)==data['control_count']and set(controls)==set(saved)
    walls=[tuple(neg(O)if i==j else Z for i in range(host.n+1))for j in range(1,host.n+1)]
    assert walls==[vec(r)for r in data['initial_walls']];counts=[]
    for height,packet in groupby(controls,key=lambda row:row[0]):
        added=[]
        for _,q,x in packet:
            roots,imprimitive,bs,ss=host.candidates(x,q);row=saved[height,q,x]
            assert set(roots)=={vec(r)for r in row['candidate_vectors']}
            assert len(roots)==row['candidates']and imprimitive==row['imprimitive']and bs==row['block_vectors']and ss==row['scalar_vectors']
            acute=[r for r in roots if all(sign(host.pair(r,s))<=0 for s in walls)]
            keep=[r for r in acute if host.compatible(r,walls)]
            assert len(keep)==row['accepted']and len(acute)-len(keep)==row['noncoxeter']
            added.extend(keep);counts.append({'height':enc(height),'norm':enc(q),'candidates':len(roots),'accepted':len(keep)})
        assert all(host.compatible(r,added[:i])for i,r in enumerate(added));walls.extend(added)
    assert set(walls)=={vec(r)for r in data['all_walls']}and len(walls)==data['wall_count']
    return {'source':str(source.relative_to(ROOT)),'status':'PASS_COMPLETE_INDEPENDENT_CENSUS','controls':counts,
        'complete_zero_orbits':len(zero),'geometry':host.check_geometry(data)}



def positive_definite(g):
    g=[list(row)for row in g]
    while g:
        p=g[0][0];assert sign(p)>0
        g=[[sub(g[i][j],divide(mul(g[i][0],g[0][j]),p))for j in range(1,len(g))]for i in range(1,len(g))]


def verify_completion(path):
    import hashlib
    d=json.loads(path.read_text());host=Host(2);walls=[vec(r)for r in d['all_walls']];inside=vec(d['strict_interior'])
    for p,h in d['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert tuple(vec(row)for row in d['gram'])==host.g
    geometry=host.check_geometry(d);face=Host(1).check_geometry(d['forced_H5_face'])
    for proof in d['finite_vertices']:
        point=vec(proof['vector']);assert sign(host.pair(point,point))<0
        assert host.pair(point,point)==dec(proof['norm'])and sign(host.pair(point,inside))<0
        values=[host.pair(point,r)for r in walls];assert all(sign(v)<=0 for v in values)
        assert values==list(vec(proof['all_wall_pairings']))
        ids=proof['indices'];assert len(ids)==9 and all(values[j]==Z for j in ids)
        gram=tuple(tuple(host.pair(walls[i],walls[j])for j in ids)for i in ids)
        assert gram==tuple(vec(row)for row in proof['gram']);positive_definite(gram)
    p=json.loads((ROOT/d['pell_source']).read_text());y=dec(p['unit_y']);s=dec(p['unit_s']);D=scale(ALPHA,3)
    assert dec(p['discriminant'])==D and sub(sq(y),mul(D,sq(s)))==O and sign(sub(y,O))>0
    for row in p['integral_candidates']:
        r=vec(row['vector']);q=dec(row['norm']);assert host.pair(r,r)==q and q in POOL
        assert all(map(integral,r))and not all(integral(divide(z,T))for z in r)
        assert all(integral(divide(scale(v,2),q))for v in mv(host.g,r))
    for row in p['compatible_candidates']:
        r=vec(row['vector']);old=json.loads((ROOT/d['source']).read_text())
        assert host.compatible(r,[vec(w)for w in old['all_walls']])
        pt=vec(row['endpoint']);rh=(r[0],*r[1:5],r[9]);fh=Host(1)
        assert fh.pair(pt,rh)==Z and fh.pair(pt,pt)==dec(row['endpoint_norm'])and sign(fh.pair(pt,pt))<0
        assert fh.pair(pt,pt)==divide(neg(ALPHA),sq(y))
    return {'source':str(path.relative_to(ROOT)),'geometry':geometry,'H5_geometry':face,
            'actual_finite_vertices':len(d['finite_vertices']),'Pell_integral_candidates':len(p['integral_candidates']),
            'Pell_compatible_alternatives':len(p['compatible_candidates']),'relative_unit_identity':True}


def structure():
    from fresh_h9_glued_host_20261013_structure import norm
    assert norm(PI)==2 and norm(ALPHA)==-2
    assert tuple(sign(ALPHA,e)for e in range(4))==(1,-1,-1,-1)
    assert all(sign(PI,e)>0 for e in range(4))
    assert integral(divide(scale(O,2),mul(mul(PI,PI),mul(PI,PI))))
    assert abs(norm(divide(scale(O,2),mul(mul(PI,PI),mul(PI,PI)))))==1
    assert abs(norm(divide(ALPHA,PI)))==1
    h=Host(2);gi=minverse(h.g)
    assert all(integral(mul(PI,z))for row in gi for z in row)
    assert not integral(gi[0][0])
    assert [norm(q)for q in (*POOL,scale(PI,2))]==[1,2,4,8,16,32]
    assert all(integral(divide(scale(PI,2),q))for q in (*POOL,scale(PI,2)))
    return {'field_discriminant':2048,'alpha':enc(ALPHA),'alpha_field_norm':'-2',
            'alpha_signatures':[1,-1,-1,-1],'dual_exponent':enc(PI),
            'safe_pool':[enc(q)for q in (*POOL,scale(PI,2))],
            'primitive_pool':[enc(q)for q in POOL],
            'excluded_norm':enc(scale(PI,2)),
            'excluded_norm_proof':'For q=2*pi, reflection integrality forces all spatial coordinates into pi*O because the spatial Gram is unimodular. If the time coordinate were a pi-unit, its term would have valuation1, whereas the spatial square has valuation at least2 and q has valuation5. Thus time also lies in pi*O and the normal is imprimitive.',
            'norm_completeness_basis':'The prior exact cyclotomic field structure certificate proves O=Z[T], class number1, a unique prime pi over2, and totally positive units are squares. For this dual exponent pi, a primitive reflection norm divides2*pi. The six valuations0..5 exhaust the possible norms modulo square units; valuation5 is excluded above.',
            'rational_space_distinction':'The new/old determinant ratio in rank10 (and rank6) is pi. Its field norm2 is not a rational square, so pi is not a field square. Scalar similarity in even dimension preserves determinant squareclass; hence this form is not similar to the old alpha=T-1 form. This compares specified forms, not all lattices or literature.'}


def main():
    import hashlib
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    files=[ROOT/'results/fresh_h9_new_host_20261013_h5_h100.json',ROOT/'results/fresh_h9_new_host_20261013_h9_h100.json']
    completion=ROOT/'results/fresh_h9_new_host_20261013_completion.json'
    outpath=ROOT/'results/fresh_h9_new_host_20261013_independent.json'
    result={'status':'PASS_COMPLETE_INITIAL_CENSUS_AND_CONSTRUCTED_GEOMETRY','structure':structure(),
       'complete_bounded_cases':[verify(p)for p in files],'completion':verify_completion(completion),
       'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[*files,completion]},
       'arithmetic':'Separate stdlib Fraction multiplication/sign/inverse implementation. Complete initial root census uses direct trace-ellipsoid alphabets and Cartesian spatial enumeration, independent of the producer unit-balanced recursive enumeration. Endpoint alternatives, all selected pairs, and geometry are replayed exactly. A separate edge_independent checker repeats the bounded edge-proposal enumeration.'}
    if args.verify:assert json.loads(outpath.read_text())==result
    else:outpath.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],flush=True)

if __name__=='__main__':main()
