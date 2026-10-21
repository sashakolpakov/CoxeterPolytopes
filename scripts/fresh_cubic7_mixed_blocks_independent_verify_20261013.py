#!/usr/bin/env python3
"""Stdlib Fraction audit, including a second bounded positive-root census.

Imports only the independent polynomial field implementation. The producer,
project field class, thin-box enumerator, NumPy and SciPy are not imported.
"""
from fractions import Fraction as F
import argparse
from functools import lru_cache,cmp_to_key
from itertools import product,combinations,groupby
from math import isqrt
import json
from pathlib import Path
import time
import fresh_cubic7_dihedral_independent_audit_20261013 as f

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_cubic7_mixed_blocks_h9_h30_20261013.json'
OUTPUT=ROOT/'results/fresh_cubic7_mixed_blocks_independent_audit_20261013.json'
Z=f.ZERO;O=f.ONE;T=f.THETA;BETA=f.BETA;C=f.C;D=f.D
A3=tuple(tuple(f.scale(O,2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(3)) for i in range(3))
AINV=tuple(tuple(f.scale(O,F(x,4)) for x in row) for row in ((3,2,1),(2,4,2),(1,2,3)))
POOL=(f.TWO,f.FOUR,f.scale(BETA,2),f.scale(BETA,4))
ALLOWED={Z,O,f.TWO,f.scale(O,3),f.add(f.TWO,T),f.add(f.TWO,C)}
DUAL=((F(1),F(-1,7),F(-3,7)),(F(-1,7),F(2,7),F(1,7)),(F(-3,7),F(1,7),F(2,7)))
BOUNDS=((F(1,6),F(1),F(1,7)),(F(4,9),F(1,4),F(1,5)),(F(1,5),F(4,9),F(1,4)))
SORT=cmp_to_key(lambda a,b:f.sign(f.sub(a,b)))

def sq(x):return f.mul(x,x)
def all_le(a,b):return all(f.sign(f.sub(b,a),e)>=0 for e in range(3))
def ceil(x):return -((-x.numerator)//x.denominator)

@lru_cache(None)
def alphabet(caps):
    """Brute force the certified coefficient box; no interval pruning."""
    if any(f.sign(c,e)<0 for e,c in enumerate(caps)):return ()
    ms=[]
    for e,c in enumerate(caps):
        m=0 if c==Z else 1
        while f.sign(f.sub(f.scale(O,m*m),c),e)<0:m*=2
        ms.append(m)
    bounds=[ceil(sum(b*m for b,m in zip(row,ms))) for row in BOUNDS]
    out=[]
    for raw in product(*(range(-b,b+1) for b in bounds)):
        x=tuple(map(F,raw));xx=sq(x)
        if all(f.sign(f.sub(c,xx),e)>=0 for e,c in enumerate(caps)):out.append(x)
    return tuple(sorted(out,key=SORT))

def modular_norm8():
    # q divisible by 8 implies time even, B coordinates divisible by 4,
    # and A3*z divisible by 4. The norm divided by 2 is 0 modulo 4.
    ring=tuple(product(range(4),repeat=3))
    def add(a,b):return tuple((x+y)%4 for x,y in zip(a,b))
    def neg(a):return tuple((-x)%4 for x in a)
    def mul(a,b):return tuple(int(x)%4 for x in f.mul(tuple(map(F,a)),tuple(map(F,b))))
    square={a:mul(a,a) for a in ring};zero=(0,0,0);count=0;residues=[]
    for x in ring:
        # First and last Cartan rows imply y=2x and y=2z modulo 4.
        y=add(x,x)
        for z in ring:
            if add(z,z)!=y or add(add(neg(x),add(y,y)),neg(z))!=zero:continue
            energy=add(add(add(square[x],square[y]),square[z]),neg(mul(y,add(x,z))))
            if energy!=zero:continue
            count+=1;assert all(c%2==0 for v in (x,y,z) for c in v)
            residues.append([list(x),list(y),list(z)])
    return {'ring':'Z[t]/(t^3+t^2-2t-1,4)','admissible_residue_count':count,
            'all_coordinates_even':True,'admissible_residues':residues}

def a3_zero_census():
    result={}
    for q in (f.TWO,f.FOUR):
        aa=[alphabet((f.mul(AINV[i][i],q),)*3) for i in range(3)]
        roots=set()
        for r in product(*aa):
            if f.pair(r,A3,r)!=q:continue
            if not all(f.integral(f.divide(f.scale(x,2),q)) for x in f.mv(A3,r)):continue
            if q==f.FOUR and all(f.integral(f.scale(x,F(1,2))) for x in r):continue
            roots.add(r)
        assert len(roots)==(12 if q==f.TWO else 6)
        assert all(x[1:]==(F(0),F(0)) for r in roots for x in r)
        for r in roots:
            coefficients=(f.neg(f.add(r[0],r[2])),f.neg(r[1]),f.neg(r[2]))
            assert all(f.sign(c)>=0 for c in coefficients) or all(f.sign(c)<=0 for c in coefficients)
        result[q]=roots
    return result

def expected_gram():
    return tuple(tuple(f.neg(f.scale(T,2)) if i==j==0 else
        f.B[(i-1)%2][(j-1)%2] if 1<=i<=6 and 1<=j<=6 and (i-1)//2==(j-1)//2 else
        A3[i-7][j-7] if i>=7 and j>=7 else Z for j in range(10)) for i in range(10))

def expected_initial():
    rows=[]
    for k in range(3):
        r=[Z]*10;r[1+2*k]=f.neg(O);rows.append(tuple(r))
        r=[Z]*10;r[1+2*k]=D;r[2+2*k]=f.neg(D);rows.append(tuple(r))
    rows.extend((Z,*(Z for _ in range(6)),*tail) for tail in ((f.neg(O),Z,Z),(Z,f.neg(O),Z),(O,Z,f.neg(O))))
    return tuple(rows)

@lru_cache(None)
def cov(r):return f.mv(expected_gram(),r)
def norm(r):return f.dot(r,cov(r))
def compatible(r,walls):
    for s in walls:
        p=f.dot(r,cov(s))
        if f.sign(p)>0:return False
        c=f.divide(f.scale(sq(p),4),f.mul(norm(r),norm(s)))
        if f.sign(f.sub(c,f.FOUR))<0 and c not in ALLOWED:return False
    return True

@lru_cache(None)
def candidates(x,q):
    target=f.add(q,f.mul(f.scale(T,2),sq(x)))
    weights=tuple(w for w in alphabet((f.scale(target,2),)*3)
        if f.sign(w)>=0 and f.integral(f.divide(f.scale(w,2),q)))
    bs=[];tails=[]
    for u,w in product(weights,repeat=2):
        if f.sign(f.sub(w,u))<0:continue
        y=(f.divide(f.add(f.scale(u,2),f.mul(C,w)),BETA),f.divide(f.add(f.mul(C,u),f.scale(w,2)),BETA))
        energy=f.dot(y,(u,w))
        if all(map(f.integral,y)) and all_le(energy,target):bs.append((energy,y))
        for v in weights:
            y=f.mv(AINV,(u,v,w));energy=f.dot(y,(u,v,w))
            if all(map(f.integral,y)) and all_le(energy,target):tails.append((energy,y))
    result=set();imprimitive=0
    # Independent Cartesian join, without the producer's sorted multiset
    # selection and permutation expansion.
    for trip in product(bs,repeat=3):
        energy=f.sum_field(row[0] for row in trip)
        for e,tail in tails:
            if f.add(energy,e)!=target:continue
            r=(x,*(v for row in trip for v in row[1]),*tail)
            assert norm(r)==q and all(f.integral(f.divide(f.scale(v,2),q)) for v in cov(r))
            if q in (f.FOUR,f.scale(BETA,4)) and all(f.integral(f.scale(v,F(1,2))) for v in r):
                imprimitive+=1;continue
            result.add(r)
    return result,{'weights':len(weights),'block_vectors':len(bs),'a3_vectors':len(tails),'imprimitive':imprimitive}

def comparisons():
    # A non-square rational field norm certifies a non-square field ratio.
    determinant=f.neg(f.mul(f.scale(T,8),f.mul(sq(BETA),BETA)))
    out=[]
    for i in range(3):
        path=f'results/fresh_cubic7_h9_seed{i}_h100_20261013.json';data=json.loads((ROOT/path).read_text())
        alpha=f.decode(data['timelike_coefficient']);prior=f.neg(alpha)
        ratio=f.divide(determinant,prior);n=f.field_norm(ratio)
        assert n>0 and (isqrt(n.numerator)**2!=n.numerator or isqrt(n.denominator)**2!=n.denominator)
        out.append({'source':path,'family':'E8+1','determinant_ratio':f.encoded(ratio),'ratio_field_norm':str(n),'nonsquare':True})
    for name in ('theta','2theta','theta_beta'):
        path=f'results/fresh_cubic7_dihedral_h9_{name}_h30_20261013.json';data=json.loads((ROOT/path).read_text())
        alpha=f.decode(data['alpha']);prior=f.neg(f.mul(alpha,sq(sq(BETA))))
        ratio=f.divide(determinant,prior);n=f.field_norm(ratio)
        assert n>0 and (isqrt(n.numerator)**2!=n.numerator or isqrt(n.denominator)**2!=n.denominator)
        out.append({'source':path,'family':'B^4+1','determinant_ratio':f.encoded(ratio),'ratio_field_norm':str(n),'nonsquare':True})
    return determinant,out

def odd_overclosure():
    # Columns identify A3 with the even-coordinate-sum sublattice D3 of O^3.
    transform=tuple(tuple(f.scale(O,x) for x in row) for row in ((0,1,0),(1,-1,1),(-1,0,1)))
    assert tuple(tuple(f.sum_field(f.mul(transform[k][i],transform[k][j]) for k in range(3))
                       for j in range(3)) for i in range(3))==A3
    assert f.determinant(transform) in (f.TWO,f.neg(f.TWO))
    omega2=(f.scale(O,F(1,2)),O,f.scale(O,F(1,2)))
    assert f.pair(omega2,A3,omega2)==O and f.mv(transform,omega2)==(O,Z,Z)
    # A dyadic discriminant class of order two is a*e0/2+b*omega2.
    # Its norm is -theta*a^2/2+b^2. Integrality forces a=0 modulo 2.
    residue_alphabet=tuple(product(range(2),repeat=3))
    allowed=[]
    for a in residue_alphabet:
        aa=tuple(map(F,a));n=f.scale(f.mul(T,sq(aa)),F(-1,2))
        if f.integral(n):allowed.append(a)
    assert allowed==[(0,0,0)]
    return {'target':'<-2theta>+B^3+I3','a3_embedding':[[f.encoded(x) for x in row] for row in transform],
        'index_ideal':'(2)','index_cardinality':8,'order_two_a3_dual_vector':f.vecencode(omega2),
        'integral_time_residue_classes_mod2':[list(a) for a in allowed],
        'full_isometry_group_equality_argument':'The dyadic dual quotient is O/(2)+O/(4). Norm integrality leaves exactly the A3 order-two line, so the integral index-(2) overclosure M is unique. Conversely L is precisely {v in M:q(v) in2O}, hence is invariant under O(M). Therefore O(L)=O(M), including the full reflection group.'}

def run(source=SOURCE,output=OUTPUT):
    start=time.monotonic();data=json.loads(source.read_text());g=expected_gram()
    assert data['completed'] and data['dimension']==9 and data['blocks']==3
    assert tuple(f.vecdecode(row) for row in data['gram'])==g
    assert f.vecdecode(data['primitive_norm_pool'])==POOL
    basis=(O,T,sq(T))
    assert [[f.trace(f.mul(d,x)) for x in basis] for d in DUAL]==[[1,0,0],[0,1,0],[0,0,1]]
    assert all(f.sign(f.sub(f.scale(O,b),d),e)>0 and f.sign(f.add(f.scale(O,b),d),e)>0
               for d,bs in zip(DUAL,BOUNDS) for e,b in enumerate(bs))
    assert f.determinant(A3)==f.FOUR and f.determinant(f.B)==BETA
    exponent=f.scale(BETA,4)
    assert f.integral(f.divide(exponent,f.scale(T,2)))
    assert all(f.integral(f.mul(exponent,x)) for row in AINV for x in row)
    binv=((f.divide(f.TWO,BETA),f.divide(C,BETA)),(f.divide(C,BETA),f.divide(f.TWO,BETA)))
    assert all(f.integral(f.mul(exponent,x)) for row in binv for x in row)
    modular=modular_norm8();dominant,binary=f.binary_census();a3=a3_zero_census()
    expected_zero=set()
    for q,tail in dominant.items():
        for k in range(3):
            r=[Z]*10;r[1+2*k:3+2*k]=tail;expected_zero.add((q,tuple(r)))
    for q,roots in a3.items():
        for tail in roots:
            u=f.mv(A3,tail)
            if f.sign(u[0])>=0 and f.sign(u[1])>=0 and f.sign(f.sub(u[2],u[0]))>=0:
                expected_zero.add((q,(Z,*(Z for _ in range(6)),*tail)))
    actual_zero={(f.decode(r['norm']),f.vecdecode(r['vector'])) for r in data['zero_height_dominant_census']}
    assert actual_zero==expected_zero and len(actual_zero)==8
    for q in POOL:
        roots,_=candidates(Z,q)
        assert {(q,r) for r in roots}=={row for row in expected_zero if row[0]==q}
    controls=[]
    for q in POOL:
        cap=f.divide(q,f.scale(T,2))
        for x in alphabet((f.scale(cap,data['height_bound']),f.neg(cap),f.neg(cap))):
            if f.sign(x)>0 and f.integral(f.divide(f.mul(f.scale(T,4),x),q)):
                controls.append((f.divide(f.mul(f.scale(T,2),sq(x)),q),q,x))
    controls.sort(key=lambda row:SORT(row[0]))
    recorded={(f.decode(r['height']),f.decode(r['norm']),f.decode(r['time'])):r for r in data['levels']}
    assert set(controls)==set(recorded) and len(controls)==data['control_count']
    walls=list(expected_initial());assert tuple(f.vecdecode(r) for r in data['initial_walls'])==tuple(walls)
    seen={f.ray_key(r) for r in walls};checks=[]
    for height,packet in groupby(controls,key=lambda row:row[0]):
        added=[]
        for _,q,x in packet:
            roots,counts=candidates(x,q);row=recorded[(height,q,x)]
            assert len(roots)==row['candidates']
            assert all(row[k]==v for k,v in counts.items())
            keep=[];angle_rejected=0
            for r in sorted(roots):
                if any(f.sign(f.dot(r,cov(s)))>0 for s in walls):continue
                if not compatible(r,walls):angle_rejected+=1;continue
                if f.ray_key(r) in seen:continue
                seen.add(f.ray_key(r));keep.append(r)
            assert len(keep)==row['accepted'] and angle_rejected==row['acute_noncoxeter_rejected']
            added.extend(keep);checks.append({'height':f.encoded(height),'norm':f.encoded(q),'candidates':len(roots),'accepted':len(keep)})
        assert all(compatible(r,added[:i]) for i,r in enumerate(added));walls.extend(added)
    saved=tuple(f.vecdecode(r) for r in data['all_walls'])
    assert {f.ray_key(r) for r in walls}=={f.ray_key(r) for r in saved} and len(walls)==data['wall_count']
    assert all(compatible(r,saved[:i]) for i,r in enumerate(saved))
    for r in saved:
        q=norm(r);assert q in POOL and all(f.sign(q,e)>0 for e in range(3))
        assert all(map(f.integral,r)) and all(f.integral(f.divide(f.scale(x,2),q)) for x in cov(r))
        if q in (f.FOUR,f.scale(BETA,4)):assert not all(f.integral(f.scale(x,F(1,2))) for x in r)
    interior=f.vecdecode(data['strict_interior']);v=f.vecdecode(data['geometry']['hyperideal_witness']['vector'])
    assert f.sign(interior[0])>0 and f.sign(norm(interior))<0 and all(f.sign(f.dot(interior,cov(r)))<0 for r in saved)
    assert f.sign(v[0])>0 and f.sign(norm(v))>0 and norm(v)==f.decode(data['geometry']['hyperideal_witness']['norm'])
    assert all(f.sign(f.dot(v,cov(r)))<=0 for r in saved)
    eps=O
    while f.sign(norm(tuple(f.add(x,f.mul(eps,y)) for x,y in zip(v,interior))))<=0:eps=f.scale(eps,F(1,2))
    strict=tuple(f.add(x,f.mul(eps,y)) for x,y in zip(v,interior))
    assert all(f.sign(f.dot(strict,cov(r)))<0 for r in saved)
    determinant,comp=comparisons()
    out={'status':'PASS','source':str(source.relative_to(ROOT)),'field':'Q(2cos(2pi/7))',
         'dual_exponent':f.encoded(exponent),'determinant':f.encoded(determinant),
         'primitive_norm_pool':f.vecencode(POOL),'norm8_imprimitivity_residue_certificate':modular,
         'odd_overclosure_identification':odd_overclosure(),
         'zero_census':{'binary_independent_census':binary,'a3_oriented_roots':[
             {'norm':f.encoded(q),'vectors':[f.vecencode(r) for r in sorted(roots)]} for q,roots in a3.items()],
             'full_oriented_roots':102,'dominant_orbits':8,'simple_wall_count':9,'chamber':'I2(14)^3 x B3'},
         'compared_host_determinants':comp,'bounded_levels_independently_reenumerated':checks,
         'wall_count':len(saved),'all_wall_pairs_checked':len(saved)*(len(saved)-1)//2,
         'escape':{'vector':f.vecencode(v),'norm':f.encoded(norm(v)),
             'strict_spacelike_vector':f.vecencode(strict),'strict_spacelike_norm':f.encoded(norm(strict))},
         'height_bound':data['height_bound'],
         'scope':f'Six compared cubic hosts have different determinant squareclasses. Complete zero-height census and independent height-{data["height_bound"]} primitive-root enumeration replay. The {len(saved)}-wall prefix has infinite volume; the full host remains open.',
         'seconds':round(time.monotonic()-start,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k in ('status','wall_count','all_wall_pairs_checked','seconds')}))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--output',type=Path,default=OUTPUT)
    a=p.parse_args();run(a.source.resolve(),a.output)
