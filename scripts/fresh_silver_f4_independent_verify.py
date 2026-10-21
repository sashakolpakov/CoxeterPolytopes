#!/usr/bin/env python3
"""Independent Fraction-only audit of the silver F4 root and cone archives.

No project imports; this checks certificates, plus two small finite-volume
calibrations. It does not independently rerun bounded root enumeration.
"""
from fractions import Fraction as F
from itertools import combinations, product
from functools import lru_cache
from pathlib import Path
import argparse
import json

Z=(F(0),F(0));O=(F(1),F(0));P=(F(0),F(1))
def parse(x):return F(x['a']),F(x['b'])
def enc(x):return {'a':str(x[0]),'b':str(x[1])}
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def scale(x,k):return x[0]*k,x[1]*k
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def inv(x):
    n=x[0]*x[0]-2*x[1]*x[1]
    return x[0]/n,-x[1]/n
def div(x,y):return mul(x,inv(y))
def conj(x):return x[0],-x[1]
def sign(x):
    a=x[0];b=x[1]
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    return ((a>0)-(a<0))*((a*a>2*b*b)-(a*a<2*b*b))
def integral(x):return all(c.denominator==1 for c in x)
def total(values):
    out=Z
    for x in values:out=add(out,x)
    return out
def dot(x,y):return total(mul(a,b) for a,b in zip(x,y,strict=True))
def mv(a,v):return tuple(dot(row,v) for row in a)
def tr(a):return tuple(zip(*a))
def mm(a,b):return tuple(tuple(dot(row,col) for col in tr(b)) for row in a)
def pair(x,a,y):return dot(x,mv(a,y))
def matrix(x):return tuple(tuple(map(parse,row)) for row in x)
def eye(n):return tuple(tuple(O if i==j else Z for j in range(n)) for i in range(n))

def rref(a):
    a=[list(row) for row in a];cols=len(a[0]);pivots=[]
    for j in range(cols):
        k=len(pivots)
        pivot=next((i for i in range(k,len(a)) if a[i][j]!=Z),None)
        if pivot is None:continue
        a[k],a[pivot]=a[pivot],a[k]
        q=a[k][j];a[k]=[div(x,q) for x in a[k]]
        for i in range(len(a)):
            if i==k:continue
            q=a[i][j];a[i]=[sub(x,mul(q,y)) for x,y in zip(a[i],a[k])]
        pivots.append(j)
        if len(pivots)==len(a):break
    return a,pivots

def inverse(a):
    n=len(a);rows,pivots=rref(tuple(tuple(row)+ident for row,ident in zip(a,eye(n))))
    assert pivots==list(range(n))
    return tuple(tuple(row[n:]) for row in rows)

def determinant(a):
    a=[list(row) for row in a];out=O
    for j in range(len(a)):
        k=next((i for i in range(j,len(a)) if a[i][j]!=Z),None)
        if k is None:return Z
        if k!=j:a[k],a[j]=a[j],a[k];out=neg(out)
        pivot=a[j][j];out=mul(out,pivot)
        for i in range(j+1,len(a)):
            q=div(a[i][j],pivot)
            a[i]=[sub(x,mul(q,y)) for x,y in zip(a[i],a[j])]
    return out

def kernel_one(a):
    rows,pivots=rref(a);free=[i for i in range(len(a[0])) if i not in pivots]
    if len(free)!=1:return None
    out=[Z]*len(a[0]);out[free[0]]=O
    for row,pivot in zip(rows,pivots):out[pivot]=neg(row[free[0]])
    return tuple(out)

ROOT=Path(__file__).resolve().parents[1]
UNIT=add(O,P);ETA=add(scale(O,2),P)
F4=tuple(tuple(scale(O,2) if i==j else neg(P) if {i,j}=={1,2}
    else neg(O) if abs(i-j)==1 else Z for j in range(4)) for i in range(4))
FI=inverse(F4)

def host(alpha,blocks,scalar):
    n=4*blocks+int(scalar)
    return tuple(tuple(neg(alpha) if i==j==0 else O if scalar and i==j==n
        else F4[(i-1)%4][(j-1)%4] if 1<=i<=4*blocks and 1<=j<=4*blocks
        and (i-1)//4==(j-1)//4 else Z for j in range(n+1)) for i in range(n+1))

def read(path):return json.loads((ROOT/path).read_text())
def oriented_key(v):
    c=next(x for x in v if x!=Z)
    return sign(c),tuple(div(x,c) for x in v)
def angle_ok(x):return sign(sub(x,scale(O,4)))>=0 or x in (Z,O,scale(O,2),scale(O,3),ETA)

def bounded(caps):
    """Independent coefficient box derived from |a±b sqrt2| bounds."""
    bounds=[]
    for cap,embedding in zip(caps,(0,1)):
        k=0
        while sign(sub(scale(O,k*k),cap if embedding==0 else conj(cap)))<0:k+=1
        bounds.append(k)
    bound=sum(bounds)+1
    return tuple((F(a),F(b)) for a in range(-bound,bound+1) for b in range(-bound,bound+1)
        if sign(sub(caps[0],mul((F(a),F(b)),(F(a),F(b)))))>=0
        and sign(conj(sub(caps[1],mul((F(a),F(b)),(F(a),F(b))))))>=0)

@lru_cache(None)
def zero_blocks(q):
    alphabet=[x for x in bounded((scale(q,2),scale(q,2))) if sign(x)>=0 and integral(div(scale(x,2),q))]
    out=[]
    # A direct fourfold product, independent of the producer's Schur recursion.
    for d in product(alphabet,repeat=4):
        v=mv(FI,d);energy=dot(d,v)
        if sign(sub(q,energy))>=0 and sign(conj(sub(q,energy)))>=0:
            assert all(map(integral,v));out.append((energy,v))
    return tuple(out)

def zero_census(pool,blocks,scalar):
    out=set()
    for q in pool:
        block=zero_blocks(q)
        scalars=[x for x in bounded((q,q)) if sign(x)>=0 and integral(div(scale(x,2),q))] if scalar else [Z]
        bysquare={mul(x,x):x for x in scalars}
        for entry in product(block,repeat=blocks):
            left=sub(q,total(e for e,_ in entry))
            if left in bysquare:
                v=(Z,*(x for _,r in entry for x in r),*((bysquare[left],) if scalar else ()))
                out.add((q,v))
    return out

def finite_f4():
    roots=set(eye(4));todo=list(roots)
    while todo:
        r=todo.pop();products=mv(F4,r)
        for j in range(4):
            s=tuple(sub(x,products[j]) if i==j else x for i,x in enumerate(r))
            assert pair(s,F4,s)==scale(O,2)
            if s not in roots:roots.add(s);todo.append(s)
        assert len(roots)<=48
    assert len(roots)==48
    dominant={r for r in roots if all(sign(x)>=0 for x in mv(F4,r))}
    assert len(dominant)==2
    assert dominant=={r for e,r in zero_blocks(scale(O,2)) if e==scale(O,2)}
    return {'finite_roots':48,'dominant_norm2_orbits':2}

def check_escape(w,roots,gram):
    v=tuple(map(parse,w['vector']));assert sign(v[0])>0
    q=pair(v,gram,v);assert q==parse(w['norm']) and sign(q)>0
    cov=tuple(mv(gram,r) for r in roots);values=tuple(dot(v,c) for c in cov)
    assert values==tuple(map(parse,w['wall_pairings'])) and all(sign(x)<=0 for x in values)
    ids=w['active_facets'];assert len(ids)==len(v)-1 and all(values[i]==Z for i in ids)
    assert len(rref([cov[i] for i in ids])[1])==len(v)-1
    return q

def verify_case(name):
    d=read('results/fresh_silver_f4_'+name+'_20261013.json');alpha=parse(d['alpha'])
    blocks=d['blocks'];scalar=d['scalar'];gram=host(alpha,blocks,scalar)
    assert gram==matrix(d['gram']) and d['completed'] is True
    assert sign(alpha)>0 and sign(conj(alpha))<0
    assert determinant(F4)==O and all(sign(determinant(tuple(row[:k] for row in F4[:k])))>0 for k in range(1,5))
    assert all(sign(conj(determinant(tuple(row[:k] for row in F4[:k]))))>0 for k in range(1,5))
    pools={UNIT:(O,ETA,scale(O,2)),P:(O,ETA,scale(O,2),scale(ETA,2)),scale(UNIT,2):(O,ETA,scale(O,2),scale(ETA,2),scale(O,4))}
    pool=tuple(map(parse,d['norm_pool']));assert pool==pools[alpha]
    actual={(parse(c['norm']),tuple(map(parse,c['vector']))) for c in d['zero_height_dominant_census']}
    assert len(actual)==len(d['zero_height_dominant_census']) and actual==zero_census(pool,blocks,scalar)
    walls=matrix(d['all_walls']);n=len(gram)-1;inside=tuple(map(parse,d['strict_interior']))
    assert n==d['dimension'] and len(walls)==d['wall_count']
    assert len({oriented_key(r) for r in walls})==len(walls)
    assert walls[:n]==tuple(tuple(neg(O) if i==j else Z for i in range(n+1)) for j in range(1,n+1))
    assert sign(inside[0])>0 and sign(pair(inside,gram,inside))<0
    cov=tuple(mv(gram,r) for r in walls);norms=tuple(dot(r,c) for r,c in zip(walls,cov))
    heights=[]
    for i,(r,c,q) in enumerate(zip(walls,cov,norms)):
        assert q in pool and all(map(integral,r)) and all(integral(div(scale(x,2),q)) for x in c)
        assert sign(dot(inside,c))<0
        h=div(mul(alpha,mul(r[0],r[0])),q);assert sign(h)>=0 and sign(sub(parse(d['height']),h))>=0
        heights.append(h)
        for j in range(i):
            v=dot(r,cov[j]);assert sign(v)<=0 and angle_ok(div(scale(mul(v,v),4),mul(q,norms[j])))
    assert all(sign(sub(b,a))>=0 for a,b in zip(heights,heights[1:]))
    # Recompute the complete time-coordinate control census in an independent box.
    expected={(div(mul(alpha,mul(x,x)),q),q,x) for q in pool
        for x in bounded((div(mul(parse(d['height']),q),alpha),neg(div(q,alpha))))
        if sign(x)>0 and integral(div(scale(mul(alpha,x),2),q))}
    levels=d['levels'];actual_controls={(parse(l['height']),parse(l['norm']),parse(l['control'])) for l in levels}
    assert actual_controls==expected and len(levels)==len(expected)
    for h,q,x in actual_controls:assert div(mul(alpha,mul(x,x)),q)==h
    result={'case':name,'walls':len(walls),'zero_census':len(actual),'controls':len(expected),'all_pairs_coxeter':True}
    if n==4:
        assert len(walls)==5
        vertices=[]
        for ids in combinations(range(5),4):
            v=kernel_one([cov[i] for i in ids]);assert v is not None and v[0]!=Z
            v=tuple(div(x,v[0]) for x in v)
            assert all(sign(dot(v,c))<=0 for c in cov) and sign(pair(v,gram,v))<0
            vertices.append(v)
        assert len(set(vertices))==5 and d['geometry']['finite_vertices']==5
        # All five affine vertices finite, plus positive barycentric interior,
        # certify the bounded simplex, not just a sample of potential vertices.
        affine=tr(vertices);bary=mv(inverse(affine),inside)
        assert all(sign(x)>0 for x in bary)
        edges=[div(scale(mul(dot(walls[i],cov[j]),dot(walls[i],cov[j])),4),mul(norms[i],norms[j]))
               for i,j in combinations(range(5),2) if dot(walls[i],cov[j])!=Z]
        assert sorted(edges)==sorted([O,O,O,O,scale(O,2)])
        result.update(compact_simplex=True,finite_vertices=5)
    else:result['escape_norm']=enc(check_escape(d['geometry']['hyperideal_witness'],walls,gram))
    return result

def identifications():
    d=read('results/fresh_silver_f4_identifications_20261013.json')
    calibration=read('results/fresh_silver_f4_h4_calibration_20261013.json')
    B=matrix(calibration['gram']);walls=matrix(calibration['all_walls']);perm=d['calibration_to_unsaturated_C5_permutation']
    old_c5=tuple(tuple(scale(O,2) if i==j else neg(P) if {i,j}=={1,2}
        else neg(O) if (i-j)%5 in (1,4) else Z for j in range(5)) for i in range(5))
    assert tuple(tuple(pair(walls[perm[i]],B,walls[perm[j]]) for j in range(5)) for i in range(5))==old_c5
    sat=tuple(tuple(O if i==j and j<4 else scale(P,F(1,2)) if j==4 and i in (1,2,4) else Z for j in range(5)) for i in range(5))
    old5=mm(tr(sat),mm(old_c5,sat));t5=matrix(d['saturated_C5_to_time_plus_F4'])
    assert mm(tr(t5),mm(B,t5))==old5
    assert all(integral(x) for row in t5 for x in row) and all(integral(x) for row in inverse(t5) for x in row)
    old=tuple(tuple(old5[i][j] if i<5 and j<5 else F4[i-5][j-5] if 5<=i<9 and 5<=j<9
        else scale(O,2) if i==j==9 else Z for j in range(10)) for i in range(10))
    M=host(UNIT,2,True);T=matrix(d['old_carrier_to_F4_squared_scalar1'])
    assert old==matrix(d['old_carrier_gram']) and M==matrix(d['target_gram'])
    assert mm(tr(T),mm(M,T))==old and determinant(T)==parse(d['transport_determinant'])
    assert mul(determinant(T),determinant(T))==scale(O,2)
    assert all(integral(x) for row in T for x in row)
    assert all(integral(mul(P,x)) for row in inverse(T) for x in row)
    # Its image is exactly {last coordinate divisible by sqrt2}.
    scal=tuple(tuple(P if i==j==9 else O if i==j else Z for j in range(10)) for i in range(10))
    quotient=mm(inverse(scal),T)
    assert all(integral(x) for row in quotient for x in row) and all(integral(x) for row in inverse(quotient) for x in row)
    field=matrix(d['F4_to_diagonal_field_isometry']);assert mm(tr(field),field)==F4
    unit=matrix(read('results/fresh_silver_f4_h9_unit_h30_20261013.json')['all_walls'])
    doubled=matrix(read('results/fresh_silver_f4_h9_2unit_h30_20261013.json')['all_walls'])
    assert {oriented_key(r) for r in unit}=={oriented_key((mul(P,r[0]),*r[1:])) for r in doubled}
    # Exhaust O/(2)'s four-element ring: its only nonzero proper ideal has
    # two elements {0,sqrt2}, proving uniqueness of an index-sqrt2 overclosure
    # when the discriminant module is cyclic O/(2).
    residues=tuple((F(a),F(b)) for a,b in product((0,1),repeat=2))
    mod2=lambda x:(int(x[0])%2,int(x[1])%2)
    ideals=[]
    for mask in range(1<<4):
        subset={mod2(residues[i]) for i in range(4) if mask>>i&1}
        if (0,0) not in subset:continue
        if all(mod2(add(x,y)) in subset and mod2(mul(x,r)) in subset for x in subset for y in subset for r in residues):ideals.append(subset)
    assert ideals==[{(0,0)},{(0,0),(0,1)},set(map(mod2,residues))]
    # Check the intrinsic norm formula on all time/scalar residue choices;
    # the spatial contribution is identically in2O by the displayed Gram.
    assert all(integral(div(x,scale(O,2))) for x in (F4[i][i] for i in range(4)))
    for x,z in product(residues,repeat=2):
        q=sub(mul(z,z),mul(UNIT,mul(x,x)));ax=int(x[0]);az=int(z[0])
        assert mod2(q)==((ax+az)%2,ax)
        assert (ax==0)==(mod2(q)[1]==0)
        assert (az==0)==((sum(mod2(q))%2)==0)
    # F4 is unimodular; each smaller lattice's only nonunit diagonal
    # block is <-2lambda> or <2>, hence cyclic discriminant O/(2).
    assert integral(inv(UNIT)) and determinant(F4)==O
    assert not integral(inv(scale(UNIT,-2))) and not integral(scale(O,F(1,2)))
    return {'integral_saturated_C5_isometry':True,'old_carrier_embedding_cardinality':2,
        'old_carrier_image_is_scalar_index_sqrt2_sublattice':True,'norm_mod2_intrinsic_sublattices':True,
        'unique_index_sqrt2_overclosure':True,'equal_full_integral_groups_argument_verified':True,
        'equal_height30_oriented_walls':len(unit),'F4_field_diagonal_isometry':True,
        'scope':'The smaller lattices are not claimed integrally isometric to M; group equality follows via the intrinsic sublattices and unique overclosure.'}

def run():
    result={'status':'INDEPENDENT_SILVER_F4_AUDIT_PASS','finite_f4':finite_f4(),
        'cases':[verify_case(n) for n in ('h4_calibration','h9_unit_h30','h9_sqrt2_h30','h9_2unit_h30')],
        'identifications':identifications(),
        'scope':'Exact certificates, zero-height and control enumeration independently repeated. Positive-height spatial enumeration source reviewed, not independently repeated in full.'}
    print(json.dumps(result,sort_keys=True));return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    a=p.parse_args();result=run()
    if a.output:a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
