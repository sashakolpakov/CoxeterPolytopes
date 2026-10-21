#!/usr/bin/env python3
"""Independent Fraction-only audit of the H4-block root and cone archives.

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
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def inv(x):
    n=x[0]*x[0]+x[0]*x[1]-x[1]*x[1]
    return (x[0]+x[1])/n,-x[1]/n
def div(x,y):return mul(x,inv(y))
def conj(x):return x[0]+x[1],-x[1]
def sign(x):
    a=2*x[0]+x[1];b=x[1]
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    return ((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))
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

H4=tuple(tuple(scale(O,2) if i==j else neg(P) if {i,j}=={0,1}
    else neg(O) if abs(i-j)==1 else Z for j in range(4)) for i in range(4))

def host(alpha,blocks,scalar):
    n=4*blocks+int(scalar)
    return tuple(tuple(neg(alpha) if i==j==0 else O if scalar and i==j==n
        else H4[(i-1)%4][(j-1)%4] if 1<=i<=4*blocks and 1<=j<=4*blocks
        and (i-1)//4==(j-1)//4 else Z for j in range(n+1)) for i in range(n+1))

def e8_host(alpha,scalar):
    size=9+int(scalar);target=[[Z]*size for _ in range(size)]
    target[0][0]=neg(alpha)
    for i in range(1,9):target[i][i]=scale(O,2)
    if scalar:target[9][9]=O
    for i,j in ((1,3),(2,4),(3,4),(4,5),(5,6),(6,7),(7,8)):target[i][j]=target[j][i]=neg(O)
    return tuple(map(tuple,target))

def bounded_integers(primary_bound,conjugate_bound):
    # If |x|<=M and |x'|<=N, then |b|<=(M+N)/2 and |a|<=N+|b|.
    m=n=0
    while sign(sub(scale(O,m*m),primary_bound))<0:m+=1
    while sign(conj(sub(scale(O,n*n),conjugate_bound)))<0:n+=1
    bmax=(m+n+1)//2;amax=n+bmax
    return tuple((F(a),F(b)) for a in range(-amax,amax+1) for b in range(-bmax,bmax+1)
        if sign(sub(primary_bound,mul((F(a),F(b)),(F(a),F(b)))))>=0
        and sign(conj(sub(conjugate_bound,mul((F(a),F(b)),(F(a),F(b))))))>=0)

@lru_cache(None)
def brute_zero_blocks(q):
    alphabet=[w for w in bounded_integers(scale(q,2),scale(q,2))
        if sign(w)>=0 and integral(div(scale(w,2),q))]
    h_inverse=inverse(H4);out=[]
    for w in product(alphabet,repeat=4):
        y=mv(h_inverse,w);energy=dot(w,y)
        if sign(sub(q,energy))>=0 and sign(conj(sub(q,energy)))>=0:
            assert all(map(integral,y));out.append((energy,y))
    return tuple(out)

def brute_zero_census(norms,blocks,scalar):
    out=set()
    for q in norms:
        entries=brute_zero_blocks(q)
        scalars=[u for u in bounded_integers(q,q) if sign(u)>=0 and integral(div(scale(u,2),q))] if scalar else [Z]
        for choices in product(entries,repeat=blocks):
            remaining=sub(q,total(c[0] for c in choices))
            for u in scalars:
                if remaining==mul(u,u):
                    r=(Z,*(x for _,y in choices for x in y),*((u,) if scalar else ()))
                    out.add((q,r))
    return out

def verify(path):
    d=json.loads(Path(path).read_text());alpha=parse(d['alpha'])
    assert sign(alpha)>0 and sign(conj(alpha))<0
    n=d['dimension'];assert n==4*d['blocks']+int(d['scalar'])
    gram=host(alpha,d['blocks'],d['scalar']);assert gram==matrix(d['gram'])
    walls=matrix(d['all_walls']);assert len(walls)==d['wall_count']
    assert d['completed'] and parse(d['completed_height'])==parse(d['requested_height'])
    bound=parse(d['completed_height']);norms=set(map(parse,d['norm_pool']))
    if alpha in (P,scale(P,2)):expected_norms={O,scale(O,2)}
    else:
        assert alpha==sub(scale(P,2),O)
        expected_norms={O,scale(O,2),add(scale(O,2),P),add(scale(O,4),scale(P,2))}
    assert norms==expected_norms
    inside=tuple(map(parse,d['strict_interior']))
    assert sign(pair(inside,gram,inside))<0 and sign(inside[0])>0
    for i,r in enumerate(walls):
        assert len(r)==n+1 and all(map(integral,r))
        q=pair(r,gram,r);assert q in norms and sign(q)>0 and sign(conj(q))>0
        assert all(integral(div(scale(v,2),q)) for v in mv(gram,r))
        assert sign(pair(inside,gram,r))<0
        assert all(sign(pair(r,gram,s))<=0 for s in walls[:i])
    assert walls[:n]==tuple(tuple(neg(O) if i==j else Z for i in range(n+1)) for j in range(1,n+1))
    cursor=n;previous=Z
    actual_controls=[]
    for level in d['levels']:
        h=parse(level['height']);q=parse(level['norm']);x=parse(level['control'])
        assert q in norms and integral(x) and sign(x)>0
        assert h==div(mul(alpha,mul(x,x)),q)
        assert sign(sub(h,previous))>=0 and sign(sub(bound,h))>=0
        assert sign(add(conj(q),mul(conj(alpha),mul(conj(x),conj(x)))))>=0
        assert integral(div(scale(mul(alpha,x),2),q))
        actual_controls.append((q,x))
        assert isinstance(level['roots'],int) and 0<=level['roots']<=level['candidates']
        selected=walls[cursor:cursor+level['roots']];assert len(selected)==level['roots']
        for r in selected:assert r[0]==x and pair(r,gram,r)==q
        cursor+=level['roots'];previous=h
    assert cursor==len(walls)
    expected_controls={(q,x) for q in norms for x in bounded_integers(div(mul(bound,q),alpha),neg(div(q,alpha)))
        if sign(x)>0 and integral(div(scale(mul(alpha,x),2),q))}
    assert len(set(actual_controls))==len(actual_controls) and set(actual_controls)==expected_controls
    zero=d['zero_height_dominant_census'];assert len(zero)==d['blocks']+int(d['scalar'])
    assert len({tuple(map(parse,z['vector'])) for z in zero})==len(zero)
    support=[]
    for z in zero:
        r=tuple(map(parse,z['vector']));q=parse(z['norm'])
        assert r[0]==Z and all(map(integral,r)) and pair(r,gram,r)==q
        assert all(sign(v)>=0 for v in mv(gram,r)[1:])
        assert all(integral(div(scale(v,2),q)) for v in mv(gram,r))
        if q==O:
            assert d['scalar'] and r[-1]==O and all(v==Z for v in r[:-1]);support.append('scalar')
        else:
            assert q==scale(O,2) and (not d['scalar'] or r[-1]==Z)
            active=[b for b in range(d['blocks']) if any(v!=Z for v in r[1+4*b:5+4*b])]
            assert len(active)==1;support.append(active[0])
    assert set(support)==set(range(d['blocks']))|({'scalar'} if d['scalar'] else set())
    assert {(parse(z['norm']),tuple(map(parse,z['vector']))) for z in zero}==brute_zero_census(norms,d['blocks'],d['scalar'])
    witness=d.get('hyperideal_witness')
    result={'artifact':str(path),'dimension':n,'walls':len(walls),'status':'ROOT_CERTIFICATE_PASS',
        'independently_complete_controls':len(expected_controls),'independent_zero_height_census':len(zero)}
    if witness:
        v=tuple(map(parse,witness['vector']));q=pair(v,gram,v)
        assert q==parse(witness['norm']) and sign(q)>0 and sign(v[0])>0
        products=tuple(pair(v,gram,r) for r in walls)
        assert products==tuple(map(parse,witness['wall_pairings'])) and all(sign(x)<=0 for x in products)
        active=witness['active_facets'];assert len(set(active))==len(active)==n
        assert all(products[i]==Z for i in active)
        assert len(rref([mv(gram,walls[i]) for i in active])[1])==n
        result.update(status='EXACT_SPACELIKE_ESCAPE_PASS',escape_norm=enc(q))
    print(json.dumps(result,sort_keys=True))
    return result

def h4_root_calibration():
    assert determinant(H4)==sub(scale(O,5),scale(P,3))
    assert all(integral(v) for row in inverse(H4) for v in row)
    for k in range(1,5):
        minor=tuple(row[:k] for row in H4[:k]);det=determinant(minor)
        assert sign(det)>0 and sign(conj(det))>0
    roots=set(eye(4));pending=list(roots)
    while pending:
        r=pending.pop()
        for i in range(4):
            pairing=dot(H4[i],r)
            v=tuple(sub(x,pairing) if j==i else x for j,x in enumerate(r))
            if v not in roots:roots.add(v);pending.append(v)
        assert len(roots)<=120
    assert len(roots)==120
    assert all(pair(r,H4,r)==scale(O,2) and all(map(integral,r)) for r in roots)
    assert all(all(sign(x)>=0 for x in r) or all(sign(x)<=0 for x in r) for r in roots)
    assert sum(all(sign(x)>=0 for x in mv(H4,r)) for r in roots)==1
    return {'h4_exact_root_closure':120,'dominant_norm_two_roots':1,'both_embeddings_positive_definite':True}

def finite_cone_calibration(path,expected_vertices):
    d=json.loads(Path(path).read_text());gram=matrix(d['gram']);walls=matrix(d['all_walls']);n=d['dimension']
    normals=tuple(mv(gram,r) for r in walls);vertices=set()
    assert len(rref(normals)[1])==n+1
    for active in combinations(range(len(walls)),n):
        v=kernel_one([normals[i] for i in active])
        if v is None:continue
        values=mv(normals,v)
        if all(sign(x)>=0 for x in values):v=tuple(map(neg,v));values=tuple(map(neg,values))
        if not all(sign(x)<=0 for x in values):continue
        assert sign(pair(v,gram,v))<0 and sign(v[0])>0
        vertices.add(tuple(div(x,v[0]) for x in v))
    assert len(vertices)==expected_vertices
    return {'artifact':str(path),'exact_finite_vertices':len(vertices),'all_extreme_rays_future_timelike':True}

def simplex_calibration(path):
    d=json.loads(Path(path).read_text());g=matrix(d['gram']);r=matrix(d['all_walls'])
    expected=tuple(tuple(scale(O,2) if i==j and i<4 else O if i==j else
        neg(P) if {i,j}=={0,1} else neg(O) if abs(i-j)==1 else Z for j in range(5)) for i in range(5))
    assert tuple(tuple(pair(x,g,y) for y in r) for x in r)==expected
    result=finite_cone_calibration(path,5);result['coxeter_chain']=[5,3,3,4]
    return result

def bugaenko_calibration(root):
    d=json.loads((root/'results/fresh_h4_blocks_h8_2phi_h30_20261013.json').read_text())
    m=json.loads((root/'results/fresh_h4_blocks_h8_bugaenko_match_20261013.json').read_text())
    source=matrix(d['gram']);walls=matrix(d['all_walls']);t=matrix(m['transform']);ti=inverse(t)
    target=e8_host(scale(P,2),False)
    assert mm(tr(t),mm(target,t))==source
    assert all(integral(x) for row in t for x in row) and all(integral(x) for row in ti for x in row)
    assert determinant(t)==parse(m['determinant'])==sub(scale(O,5),scale(P,3))
    primitive=[tuple(neg(O) if i==j else Z for i in range(9)) for j in range(1,9)]
    primitive.extend((
        (sub(P,O),*(scale(O,k) for k in (2,3,4,6,5,4,3,2))),
        (O,*(scale(O,k) for k in (4,5,7,10,8,6,4,2))),
        tuple((F(a),F(b)) for a,b in ((1,2),(4,6),(7,9),(8,12),(12,18),(9,15),(6,12),(4,8),(2,4)))
    ))
    published=[tuple(mul(P,x) for x in primitive[9]),primitive[0],primitive[2],primitive[3],primitive[4],
        primitive[5],primitive[6],primitive[7],tuple(mul(P,x) for x in primitive[8]),primitive[1],tuple(mul(P,x) for x in primitive[10])]
    assert sorted(map(int,m['wall_map'].values()))==list(range(11))
    for i,r in enumerate(walls):assert mv(t,r)==published[m['wall_map'][str(i)]]
    for i,r in enumerate(walls):
        for j,s in enumerate(walls):
            assert pair(r,source,s)==pair(published[m['wall_map'][str(i)]],target,published[m['wall_map'][str(j)]])
    return {'bugaenko_exact_wall_gram_match':11,'integral_isometry_and_inverse':True,'determinant':enc(determinant(t))}

def h9_isometry_calibrations(root):
    output=[]
    for name,alpha in (('phi',P),('2phi',scale(P,2)),('sqrt5',sub(scale(P,2),O))):
        path=root/f'results/fresh_h4_blocks_h9_{name}_isometry_20261013.json'
        d=json.loads(path.read_text());assert parse(d['alpha'])==alpha
        t=matrix(d['transform']);ti=inverse(t);target=e8_host(alpha,True)
        assert target==matrix(d['target_gram'])
        assert mm(tr(t),mm(target,t))==host(alpha,2,True)
        assert all(integral(x) for row in t for x in row) and all(integral(x) for row in ti for x in row)
        assert determinant(t)==parse(d['determinant'])==sub(scale(O,5),scale(P,3))
        output.append({'alpha':name,'status':'EXACT_INTEGRAL_ISOMETRY_TO_PREVIOUS_E8_HOST_PASS',
            'determinant':enc(determinant(t)),'artifact':str(path)})
    return output

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('artifacts',type=Path,nargs='*')
    p.add_argument('--calibrations',action='store_true');p.add_argument('--output',type=Path)
    args=p.parse_args();results={'certificates':[verify(a) for a in args.artifacts]}
    if args.calibrations:
        root=Path(__file__).resolve().parents[1]
        results['calibrations']=[h4_root_calibration(),bugaenko_calibration(root),*h9_isometry_calibrations(root),
            simplex_calibration(root/'results/fresh_h4_blocks_independent_h4_phi_h30.json'),
            finite_cone_calibration(root/'results/fresh_h4_blocks_h8_2phi_h30_20261013.json',41)]
        print(json.dumps(results['calibrations'],sort_keys=True))
    if args.output:args.output.write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
