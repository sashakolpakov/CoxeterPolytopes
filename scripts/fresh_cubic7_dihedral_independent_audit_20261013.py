#!/usr/bin/env python3
"""Independent Fraction-polynomial audit of the cubic-seven dihedral hosts.

No repository modules, search code, field arithmetic, or numeric libraries
are imported. This independently enumerates the complete binary zero-height root system,
then verifies saved roots, level records, and exact geometry certificates.
It does not independently repeat the positive-height candidate enumeration.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path

ZERO=(F(0),F(0),F(0)); ONE=(F(1),F(0),F(0)); THETA=(F(0),F(1),F(0))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def sub(x,y):return add(x,neg(y))
def scale(x,n):return tuple(n*a for a in x)
def mul(x,y):
    p=[F(0)]*5
    for i,a in enumerate(x):
        for j,b in enumerate(y):p[i+j]+=a*b
    # Independent generic reduction by theta^3 + theta^2 -2 theta -1.
    for d in (4,3):
        p[d-3]+=p[d];p[d-2]+=2*p[d];p[d-1]-=p[d]
    return tuple(p[:3])
@lru_cache(None)
def inverse(x):
    cols=[mul(x,b) for b in (ONE,THETA,mul(THETA,THETA))]
    matrix=[[cols[j][i] for j in range(3)]+[F(i==0)] for i in range(3)]
    for k in range(3):
        pivot=next(i for i in range(k,3) if matrix[i][k])
        matrix[k],matrix[pivot]=matrix[pivot],matrix[k]
        t=matrix[k][k];matrix[k]=[z/t for z in matrix[k]]
        for i in range(3):
            if i!=k:
                t=matrix[i][k];matrix[i]=[a-t*b for a,b in zip(matrix[i],matrix[k])]
    return tuple(matrix[i][3] for i in range(3))
def divide(x,y):return mul(x,inverse(y))
def decode(x):return tuple(F(x[k]) for k in ('a','b','c'))
def encoded(x):return dict(zip(('a','b','c'),map(str,x)))
def polynomial(t):return t**3+t*t-2*t-1
def interval_mul(x,y):
    p=[a*b for a in x for b in y];return min(p),max(p)
@lru_cache(None)
def sign(x,place=0):
    if x==ZERO:return 0
    low,high=((F(1),F(2)),(F(-1),F(0)),(F(-2),F(-1)))[place]
    low_polynomial_positive=polynomial(low)>0
    while True:
        p=interval_mul((x[2],x[2]),(low,high))
        p=(p[0]+x[1],p[1]+x[1]);p=interval_mul(p,(low,high))
        p=(p[0]+x[0],p[1]+x[0])
        if p[0]>0:return 1
        if p[1]<0:return -1
        mid=(low+high)/2
        if (polynomial(mid)>0)==low_polynomial_positive:low=mid
        else:high=mid

TWO=scale(ONE,2); FOUR=scale(ONE,4)
C=sub(add(mul(THETA,THETA),THETA),ONE)
D=sub(TWO,mul(THETA,THETA)); BETA=sub(TWO,THETA)
B=((TWO,neg(C)),(neg(C),TWO))

def integral(x):return all(a.denominator==1 for a in x)
def dot(x,y):return sum_field(mul(a,b) for a,b in zip(x,y))
def sum_field(values):
    out=ZERO
    for x in values:out=add(out,x)
    return out
def mv(a,v):return tuple(dot(row,v) for row in a)
def pair(x,g,y):return dot(x,mv(g,y))
def trace(x):return 3*x[0]-x[1]+5*x[2]
def determinant(a):
    if not a:return ONE
    total=ZERO
    for j,x in enumerate(a[0]):
        minor=tuple(tuple(row[k] for k in range(len(a)) if k!=j) for row in a[1:])
        total=add(total,scale(mul(x,determinant(minor)),(-1)**j))
    return total
def field_norm(x):
    cols=[mul(x,b) for b in (ONE,THETA,mul(THETA,THETA))]
    a=[[cols[j][i] for j in range(3)] for i in range(3)]
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
def rank(a):
    a=[list(row) for row in a];k=0
    for j in range(len(a[0]) if a else 0):
        p=next((i for i in range(k,len(a)) if a[i][j]!=ZERO),None)
        if p is None:continue
        a[k],a[p]=a[p],a[k];z=inverse(a[k][j]);a[k]=[mul(z,x) for x in a[k]]
        for i in range(k+1,len(a)):
            z=a[i][j]
            if z!=ZERO:a[i]=[sub(x,mul(z,y)) for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k
def ray_key(v):
    x=next(x for x in v if x!=ZERO)
    return sign(x),tuple(divide(y,x) for y in v)
def vecdecode(v):return tuple(map(decode,v))
def vecencode(v):return list(map(encoded,v))

def binary_census():
    # Trace-dual bounds independently enclose EVERY integral coordinate of
    # F(x,y)=x^2-Cxy+y^2=1, since x^2,y^2 <= 4/beta at every place.
    assert mul(C,C)==add(TWO,THETA)
    assert mul(D,D)==sub(TWO,C)
    assert field_norm(BETA)==7 and field_norm(C)==-1
    assert abs(field_norm(D))==1 and field_norm(add(TWO,C))==7
    assert field_norm(sub(TWO,C))==1
    assert [sign(THETA,e) for e in range(3)]==[1,-1,-1]
    assert [sign(C,e) for e in range(3)]==[1,-1,1]
    assert all(sign(BETA,e)>0 for e in range(3))
    assert all(sign(sub(TWO,C),e)>0 and sign(add(TWO,C),e)>0 for e in range(3))
    dual=((F(1),F(-1,7),F(-3,7)),(F(-1,7),F(2,7),F(1,7)),(F(-3,7),F(1,7),F(2,7)))
    basis=(ONE,THETA,mul(THETA,THETA))
    assert [[trace(mul(x,y)) for y in basis] for x in dual]==[[1,0,0],[0,1,0],[0,0,1]]
    for i,x in enumerate(dual):
        limit=scale(ONE,1 if i==0 else F(1,2))
        assert all(sign(sub(limit,x),e)>0 and sign(add(limit,x),e)>0 for e in range(3))
    cap=divide(FOUR,BETA);embedding_bounds=(3,2,2)
    assert all(sign(sub(scale(ONE,m*m),cap),e)>0 for e,m in enumerate(embedding_bounds))
    coordinates=[]
    for raw in product(range(-7,8),range(-4,5),range(-4,5)):
        x=tuple(map(F,raw))
        if all(sign(sub(cap,mul(x,x)),e)>=0 for e in range(3)):coordinates.append(x)
    solutions={(x,y) for x in coordinates for y in coordinates
               if sub(add(mul(x,x),mul(y,y)),mul(C,mul(x,y)))==ONE}
    def orbit(v):
        out=[]
        for _ in range(14):
            out.append(v);v=(sub(mul(C,v[0]),v[1]),v[0])
        assert v==out[0] and len(set(out))==14
        return set(out)
    norm2=orbit((ONE,ZERO));assert solutions==norm2 and len(solutions)==14
    norm2beta=orbit((neg(D),D))
    assert all(pair(v,B,v)==scale(BETA,2) for v in norm2beta)
    # Reflection-integral q=2beta vectors are adj(B)t, with F_plus(t)=1.
    from_dual={(add(scale(x,2),mul(C,neg(y))),add(mul(C,x),scale(neg(y),2))) for x,y in solutions}
    assert from_dual==norm2beta
    dominant={}
    for q,roots in ((TWO,norm2),(scale(BETA,2),norm2beta)):
        for v in roots:
            assert all(integral(divide(scale(x,2),q)) for x in mv(B,v))
        selected=[v for v in roots if sign(mv(B,v)[0])>=0 and sign(sub(mv(B,v)[1],mv(B,v)[0]))>=0]
        assert len(selected)==1;dominant[q]=selected[0]
    a=(ONE,ZERO);b=(neg(D),D)
    assert divide(scale(mul(pair(a,B,b),pair(a,B,b)),4),mul(pair(a,B,a),pair(b,B,b)))==add(TWO,C)
    # The norm-2beta reflection is literally the coordinate swap.
    for e in ((ONE,ZERO),(ZERO,ONE)):
        reflected=tuple(sub(x,mul(divide(scale(pair(e,B,b),2),pair(b,B,b)),y)) for x,y in zip(e,b))
        assert reflected==tuple(reversed(e))
    result={'coordinate_coefficient_bounds':[7,4,4],'coordinate_embedding_integer_bounds':[3,2,2],
            'coordinate_count':len(coordinates),'binary_norm2_oriented_roots':len(norm2),
            'binary_norm2beta_oriented_roots':len(norm2beta),
            'dominant_binary_roots':[{'norm':encoded(q),'vector':vecencode(v)} for q,v in dominant.items()],
            'field_discriminant':49,'beta_norm':7,'two_plus_C_norm':7,
            'pi14_four_cos_squared':encoded(add(TWO,C)),
            'full_H9_primitive_oriented_zero_roots':114,
            'complete_zero_census_argument':'Trace positivity forces each non-scalar root into exactly one block; see companion report.'}
    return dominant,result

def make_gram(alpha,blocks,scalar):
    n=2*blocks+int(scalar)
    return tuple(tuple(neg(alpha) if i==j==0 else ONE if scalar and i==j==n
                 else B[(i-1)%2][(j-1)%2] if 1<=i<=2*blocks and 1<=j<=2*blocks and (i-1)//2==(j-1)//2
                 else ZERO for j in range(n+1)) for i in range(n+1))

def verify(data,dominant):
    assert data['completed'] is True
    alpha=decode(data['alpha']);blocks=data['blocks'];scalar=data['scalar'];n=2*blocks+int(scalar)
    assert alpha in (THETA,scale(THETA,2),mul(THETA,BETA))
    assert blocks in (1,2,4) and type(scalar) is bool and data['dimension']==n
    g=make_gram(alpha,blocks,scalar)
    assert tuple(vecdecode(row) for row in data['gram'])==g
    assert [sign(alpha,e) for e in range(3)]==[1,-1,-1]
    norms=(ONE,TWO,FOUR,BETA,scale(BETA,2),scale(BETA,4)) if alpha==scale(THETA,2) else (ONE,TWO,BETA,scale(BETA,2))
    assert vecdecode(data['norm_pool'])==norms
    dual_exponent=scale(BETA,2) if alpha==scale(THETA,2) else BETA
    assert integral(divide(dual_exponent,alpha))
    assert all(integral(divide(mul(dual_exponent,x),BETA)) for row in ((TWO,C),(C,TWO)) for x in row)
    assert all(integral(divide(scale(dual_exponent,2),q)) for q in norms)
    expected_zero=set()
    for q,v in dominant.items():
        for b in range(blocks):
            root=[ZERO]*(n+1);root[1+2*b:3+2*b]=v;expected_zero.add((q,tuple(root)))
    if scalar:
        for q,z in ((ONE,ONE),(FOUR,TWO)):
            if q in norms:expected_zero.add((q,tuple([ZERO]*n+[z])))
    actual_zero=[(decode(row['norm']),vecdecode(row['vector'])) for row in data['zero_height_dominant_census']]
    assert len(actual_zero)==len(expected_zero) and set(actual_zero)==expected_zero
    walls=tuple(vecdecode(row) for row in data['all_walls'])
    assert len(walls)==data['wall_count'] and all(len(w)==n+1 for w in walls)
    initial=[]
    for b in range(blocks):
        w=[ZERO]*(n+1);w[1+2*b]=neg(ONE);initial.append(tuple(w))
        w=[ZERO]*(n+1);w[1+2*b]=D;w[2+2*b]=neg(D);initial.append(tuple(w))
    if scalar:initial.append(tuple([ZERO]*n+[neg(ONE)]))
    assert walls[:n]==tuple(initial)
    assert len({ray_key(w) for w in walls})==len(walls)
    cov=tuple(mv(g,w) for w in walls);squares=tuple(dot(w,c) for w,c in zip(walls,cov))
    for i,(w,q) in enumerate(zip(walls,squares)):
        assert all(integral(x) for x in w) and q in norms
        assert all(sign(q,e)>0 for e in range(3))
        assert all(integral(divide(scale(x,2),q)) for x in cov[i])
    allowed={ZERO:2,ONE:3,TWO:4,scale(ONE,3):6,add(TWO,THETA):7,add(TWO,C):14}
    angle_counts={}
    for i,j in combinations(range(len(walls)),2):
        p=dot(walls[i],cov[j]);assert sign(p)<=0
        v=divide(scale(mul(p,p),4),mul(squares[i],squares[j]))
        if sign(sub(v,FOUR))<0:
            assert v in allowed
            m=str(allowed[v]);angle_counts[m]=angle_counts.get(m,0)+1
    height=decode(data['height_bound']);assert sign(height)>0
    assert decode(data['completed_height'])==height
    assert len(data['levels'])==data['control_count']
    previous=ZERO;cursor=n;controls=set()
    for level in data['levels']:
        h=decode(level['height']);q=decode(level['norm']);x=decode(level['time'])
        assert q in norms and integral(x) and sign(x)>0 and (x,q) not in controls
        controls.add((x,q));assert integral(divide(scale(mul(alpha,x),2),q))
        assert h==divide(mul(alpha,mul(x,x)),q)
        assert sign(sub(h,previous))>=0 and sign(sub(height,h))>=0
        assert all(sign(add(q,mul(alpha,mul(x,x))),e)>=0 for e in range(3))
        previous=h;count=level['accepted'];assert type(count) is int and count>=0
        assert type(level['candidate_roots']) is int and level['candidate_roots']>=count
        assert type(level['block_vectors']) is int and level['block_vectors']>=0
        assert all(w[0]==x and s==q for w,s in zip(walls[cursor:cursor+count],squares[cursor:cursor+count]))
        cursor+=count
    assert cursor==len(walls)
    center=vecdecode(data['strict_interior']);assert len(center)==n+1 and sign(center[0])>0
    assert sign(pair(center,g,center))<0 and all(sign(dot(center,c))<0 for c in cov)
    geom=data['geometry']
    def check_vertex(row,desired_sign=None):
        v=vecdecode(row['vector']);assert len(v)==n+1 and v[0]==ONE
        q=pair(v,g,v);assert q==decode(row['norm'])
        ps=tuple(dot(v,c) for c in cov);assert ps==vecdecode(row['wall_pairings']) and all(sign(p)<=0 for p in ps)
        ids=row['active_facets']
        if ids:
            assert len(ids)==n and len(set(ids))==n and all(0<=i<len(walls) for i in ids)
            assert all(ps[i]==ZERO for i in ids) and rank([cov[i][1:] for i in ids])==n
        else:
            assert desired_sign==1 and all(sign(p)<0 for p in ps)
        if desired_sign is not None:assert sign(q)==desired_sign
        return v,q,ps
    out={'dimension':n,'alpha':encoded(alpha),'height_bound':encoded(height),'walls':len(walls),
         'verified_control_levels':len(controls),'zero_dominant_count':len(actual_zero),
         'finite_angle_counts_by_label':angle_counts,'center_norm':encoded(pair(center,g,center))}
    if geom['status']=='EXACT_ESCAPE':
        v,q,ps=check_vertex(geom['hyperideal_witness'],1)
        # A strict spacelike point proves an open ideal-boundary patch, hence
        # infinite volume.  Construct it without any floating-point step.
        eps=F(1)
        while True:
            strict=tuple(add(x,scale(y,eps)) for x,y in zip(v,center))
            strict_q=pair(strict,g,strict)
            if sign(strict_q)>0:break
            eps/=2
        assert sign(strict[0])>0 and all(sign(dot(strict,c))<0 for c in cov)
        out.update(verdict='EXACT_INFINITE_VOLUME_PREFIX',escape_norm=encoded(q),
                   strict_escape=vecencode(strict),strict_escape_norm=encoded(strict_q),
                   strict_escape_center_coefficient=str(eps),escape_nonzero_spatial_blocks=[b for b in range(blocks) if any(x!=ZERO for x in v[1+2*b:3+2*b])],
                   escape_scalar_zero=(v[-1]==ZERO) if scalar else None)
    else:
        assert geom['hyperideal_witness'] is None
        vertices=[check_vertex(row) for row in geom['vertices']]
        assert all(sign(q)<=0 for v,q,ps in vertices)
        if geom['status']=='EXACT_COMPACT_TRIANGLE_CALIBRATION':
            assert n==2 and len(walls)==3 and len(vertices)==3
            assert {tuple(v['active_facets']) for v in geom['vertices']}==set(combinations(range(3),2))
            assert all(sign(q)<0 for v,q,ps in vertices)
            coeffs=vecdecode(geom['positive_affine_normal_dependence'])
            assert len(coeffs)==3 and all(sign(x)>0 for x in coeffs)
            assert all(sum_field(mul(a,c[j]) for a,c in zip(coeffs,cov))==ZERO for j in range(1,n+1))
            assert rank([c[1:] for c in cov])==n
            assert sorted(int(m) for m,count in angle_counts.items() for _ in range(count))==geom['coxeter_labels']==[2,3,14]
            out.update(verdict='EXACT_COMPACT_KNOWN_2_3_14_TRIANGLE',positive_affine_normal_dependence=vecencode(coeffs))
        else:
            assert geom['status']=='EXHAUSTIVE_AFFINE_VERTEX_CENSUS'
            out.update(verdict='SAVED_AFFINE_VERTICES_VALIDATED_NO_VOLUME_CONCLUSION',saved_vertices=len(vertices))
    return out

def verify_pell_cap(path):
    data=json.loads(path.read_text());source=json.loads(Path(data['source']).read_text())
    assert data['source']=='results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
    assert data['status']=='EXACT_COMPATIBLE_PELL_CAP_WITH_REMAINING_ESCAPE'
    alpha=scale(THETA,2);g=make_gram(alpha,4,True)
    root=vecdecode(data['root']);x=(F(6),F(18),F(8));z=(F(11),F(28),F(12))
    assert root==(x,*([ZERO]*8),z) and sign(x)>0 and sign(z)>0
    assert sub(mul(z,z),mul(alpha,mul(x,x)))==ONE
    assert pair(root,g,root)==decode(data['norm'])==ONE
    assert all(integral(scale(a,2)) for a in mv(g,root))
    h=mul(alpha,mul(x,x));assert h==decode(data['height']) and sign(sub(h,scale(ONE,30)))>0
    walls=tuple(vecdecode(v) for v in source['all_walls']);assert len(walls)==25
    ps=tuple(pair(root,g,w) for w in walls)
    assert ps==vecdecode(data['wall_pairings']) and all(sign(p)<=0 for p in ps)
    new_finite_angles=[]
    allowed={ZERO:2,ONE:3,TWO:4,scale(ONE,3):6,add(TWO,THETA):7,add(TWO,C):14}
    for w,p in zip(walls,ps):
        value=divide(scale(mul(p,p),4),pair(w,g,w))
        if sign(sub(value,FOUR))<0:
            assert value in allowed;new_finite_angles.append(allowed[value])
    center=vecdecode(source['strict_interior']);escape=vecdecode(source['geometry']['hyperideal_witness']['vector'])
    assert sign(pair(center,g,center))<0 and all(sign(pair(center,g,w))<0 for w in walls)
    assert sign(pair(escape,g,escape))>0 and all(sign(pair(escape,g,w))<=0 for w in walls)
    pc=pair(root,g,center);pe=pair(root,g,escape)
    assert pc==decode(data['strict_interior_pairing']) and sign(pc)<0
    assert pe==decode(data['retained_escape_pairing']) and sign(pe)<0
    return {'source':str(path),'verdict':'EXACT_COMPATIBLE_PELL_CAP_WITH_RETAINED_SPACELIKE_ESCAPE',
            'root_norm':encoded(ONE),'height':encoded(h),'center_pairing':encoded(pc),
            'retained_escape_pairing':encoded(pe),'new_finite_angles':new_finite_angles}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files',type=Path,nargs='*')
    parser.add_argument('--output',type=Path,default=Path('results/fresh_cubic7_dihedral_independent_audit_20261013.json'))
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();dominant,census=binary_census()
    paths=args.files or sorted(Path('results').glob('fresh_cubic7_dihedral_h[2349]_*.json'))
    results=[]
    for path in paths:
        data=json.loads(path.read_text());row=verify(data,dominant)
        forged=copy.deepcopy(data);forged['levels'][0]['height']=encoded(add(decode(forged['levels'][0]['height']),ONE))
        try:verify(forged,dominant)
        except AssertionError:pass
        else:raise AssertionError('forged control height was accepted')
        row['forged_level_height_rejected']=True
        if data['geometry']['status']=='EXACT_ESCAPE':
            forged=copy.deepcopy(data);bad=forged['geometry']['hyperideal_witness']
            bad['vector']=copy.deepcopy(forged['strict_interior'])
            g=make_gram(decode(forged['alpha']),forged['blocks'],forged['scalar']);v=vecdecode(bad['vector'])
            bad['norm']=encoded(pair(v,g,v));bad['wall_pairings']=vecencode(pair(v,g,vecdecode(w)) for w in forged['all_walls'])
            try:verify(forged,dominant)
            except AssertionError:pass
            else:raise AssertionError('timelike fake escape was accepted')
            row['timelike_fake_escape_rejected']=True
        row['source']=str(path);results.append(row)
    result={'status':'PASS','arithmetic':'standard-library Fraction triples; exact isolating-interval sign decisions; no repository imports',
            'scope':'Independent complete binary zero-height census, exact saved-root/chronology/interior/geometry replay. Positive-height enumeration is reviewed mathematically but not independently rerun.',
            'zero_height_census':census,'cases':results,
            'pell_cap':verify_pell_cap(Path('results/fresh_cubic7_dihedral_pell_cap_20261013.json'))}
    if args.verify:assert json.loads(args.output.read_text())==result
    else:args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','cases':len(results),'output':str(args.output),'binary_roots':census['binary_norm2_oriented_roots']+census['binary_norm2beta_oriented_roots'],'verdicts':[(r['source'],r['verdict']) for r in results]},indent=2))

if __name__=='__main__':main()
