#!/usr/bin/env python3
"""Separate Fraction A9 census through H150, plus all saved root/geometry replay."""
import argparse,hashlib,json
from functools import cmp_to_key,lru_cache
from itertools import groupby
from pathlib import Path
import fresh_h9_new_host_20261013_independent as a
from fresh_h9_new_host_20261013_edge_independent import balanced
from fresh_h9_glued_host_20261013_structure import norm

ALPHA=a.sub(a.T,a.O)
A=tuple(tuple(a.scale(a.O,2 if i==j else -1 if abs(i-j)==1 else 0)for j in range(9))for i in range(9))
AI=a.minverse(A)
G=tuple(tuple(a.neg(ALPHA)if i==j==0 else A[i-1][j-1]if i and j else a.Z for j in range(10))for i in range(10))
DYAD=(a.O,a.PI,a.sq(a.T),a.divide(a.scale(a.O,2),a.PI),a.scale(a.O,2),a.scale(a.PI,2),a.scale(a.sq(a.T),2),a.divide(a.scale(a.O,4),a.PI),a.scale(a.O,4))
POOL=DYAD+tuple(a.scale(q,5)for q in DYAD)
P=lambda x,y:a.dot(x,a.mv(G,y))

def compatible(r,walls):
    q=P(r,r)
    for s in walls:
        p=P(r,s)
        if a.sign(p)>0:return False
        ratio=a.divide(a.scale(a.sq(p),4),a.mul(q,P(s,s)))
        if a.sign(a.sub(ratio,a.scale(a.O,4)))<0 and ratio not in a.ALLOWED:return False
    return True

@lru_cache(None)
def spatial(E,q):
    alph=[]
    for j in range(9):
        caps=(a.divide(a.scale(E,4),a.mul(a.sq(q),AI[j][j])),*(a.divide(a.scale(E,8),a.sq(q))for _ in range(3)))
        alph.append(tuple(a.scale(a.mul(q,z),a.Q(1,2))for z in balanced(caps)if a.integral(a.scale(a.mul(q,z),a.Q(1,2)))))
    # Reverse the coordinate order and form conditional lower bounds by a
    # direct Schur complement of AI (the producer uses inverse Cartan minors).
    order=list(reversed(range(9)));Q=tuple(tuple(AI[i][j]for j in order)for i in order);schur=[]
    for k in range(1,10):
        tail=a.minverse(tuple(row[k:]for row in Q[k:]))if k<9 else ()
        schur.append(tuple(tuple(a.sub(Q[i][j],a.dot(Q[i][k:],a.mv(tail,tuple(Q[l][j]for l in range(k,9)))))if k<9 else Q[i][j]for j in range(k))for i in range(k)))
    found=[]
    def walk(pre):
        k=len(pre)
        if k:
            lower=a.dot(pre,a.mv(schur[k-1],pre))
            if any(a.sign(a.sub(E,lower),e)<0 for e in range(4)):return
            partial=a.dot(pre,a.mv(tuple(row[:k]for row in Q[:k]),pre))
            if a.sign(a.sub(E,partial))<0:return
        if k==9:
            if lower!=E:return
            cov=tuple(reversed(pre));y=a.mv(AI,cov)
            if all(map(a.integral,y)):found.append(y)
            return
        for v in alph[order[k]]:walk((*pre,v))
    walk(())
    return found

def run(path,limit):
    d=json.loads(path.read_text());assert d['completed']and tuple(a.vec(row)for row in d['gram'])==G and a.vec(d['norm_pool'])==POOL
    assert norm(ALPHA)==-1 and tuple(a.sign(ALPHA,e)for e in range(4))==(1,-1,-1,-1)
    assert all(a.integral(a.scale(z,10))for row in a.minverse(G)for z in row)
    assert AI==tuple(tuple(a.scale(a.O,a.Q(min(i,j)*(10-max(i,j)),10))for j in range(1,10))for i in range(1,10))
    # Inertness independently via x^(5^d) gcd tests for d=1,2; f has degree4.
    import sympy as sp
    x=sp.Symbol('x');f=sp.Poly(x**4-4*x*x+2,x,modulus=5)
    assert sp.gcd(f,sp.Poly(x**5-x,x,modulus=5)).degree()==0
    assert sp.gcd(f,sp.Poly(x**25-x,x,modulus=5)).degree()==0
    height=min(limit,d['height_bound']);controls=[]
    for q in POOL:
        caps=(a.divide(a.scale(q,height),ALPHA),*(a.divide(a.neg(q),ALPHA)for _ in range(3)))
        for t in balanced(caps):
            if a.sign(t)>0 and a.integral(a.divide(a.scale(a.mul(ALPHA,t),2),q)):
                controls.append((a.divide(a.mul(ALPHA,a.sq(t)),q),q,t))
    controls.sort(key=cmp_to_key(lambda r,s:next((a.sign(a.sub(u,v))for u,v in zip(r,s)if u!=v),0)))
    saved={(a.dec(row['height']),a.dec(row['norm']),a.dec(row['time'])):row for row in d['levels']if a.sign(a.sub(a.scale(a.O,height),a.dec(row['height'])))>=0}
    assert set(controls)==set(saved)
    walls=[tuple(a.neg(a.O)if i==j else a.Z for i in range(10))for j in range(1,10)];count=0
    for h,group in groupby(controls,key=lambda row:row[0]):
        new=[]
        for _,q,t in group:
            E=a.add(q,a.mul(ALPHA,a.sq(t)));raw=[(t,*y)for y in spatial(E,q)]
            roots=[r for r in raw if not all(a.integral(a.divide(z,a.T))for z in r)and not all(a.integral(a.scale(z,a.Q(1,5)))for z in r)]
            row=saved[h,q,t];assert set(roots)==set(a.vec(r)for r in row['candidates'])
            assert len(raw)-len(roots)==row['imprimitive'];count+=len(roots)
            keep=[r for r in roots if compatible(r,walls)]
            assert set(keep)==set(a.vec(r)for r in row['accepted']);new+=keep
        assert all(compatible(r,new[:i])for i,r in enumerate(new));walls+=new
    allwalls=[a.vec(r)for r in d['all_walls']];inside=a.vec(d['strict_interior'])
    assert a.sign(P(inside,inside))<0 and all(a.sign(P(inside,r))<0 for r in allwalls)
    for i,r in enumerate(allwalls):
        q=P(r,r);assert q in POOL and all(a.integral(z)for z in r)
        assert not all(a.integral(a.divide(z,a.T))for z in r)and not all(a.integral(a.scale(z,a.Q(1,5)))for z in r)
        assert all(a.integral(a.divide(a.scale(z,2),q))for z in a.mv(G,r))and compatible(r,allwalls[:i])
    geom=d['geometry'];v=a.vec(geom['strict_vector']);assert a.sign(v[0])>0 and a.sign(P(v,v))>0
    assert P(v,v)==a.dec(geom['strict_norm'])and P(v,inside)==a.dec(geom['strict_center_pairing'])and a.sign(P(v,inside))<0
    assert all(a.sign(P(v,r))<0 for r in allwalls)
    return {'source':str(path.relative_to(a.ROOT)),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'status':'PASS_INDEPENDENT_A9_CENSUS_AND_GEOMETRY','complete_spatial_census_height':height,
      'controls':len(controls),'primitive_candidates':count,'walls_through_complete_bound':len(walls),
      'all_saved_walls':len(allwalls),'all_pair_checks':len(allwalls)*(len(allwalls)-1)//2,
      'strict_escape_norm':a.enc(P(v,v)),'scope':'Separate Fraction field arithmetic, independently balanced trace alphabets, reversed spatial recursion with Schur-complement lower bounds; complete census only through the stated independent height. Every saved wall and the final exact escape are replayed.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sources',nargs='+',type=Path);p.add_argument('--limit',type=int,default=150);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify',action='store_true');args=p.parse_args()
    result={'status':'PASS_A9_INDEPENDENT_AUDIT','cases':[run(path.resolve(),args.limit)for path in args.sources]}
    if args.verify:assert json.loads(args.output.read_text())==result
    else:args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['cases'],flush=True)
