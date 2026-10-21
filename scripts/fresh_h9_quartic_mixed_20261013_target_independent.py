#!/usr/bin/env python3
"""Independent Fraction census of mixed/neighboring slices and exact cones."""
import argparse
import copy
from functools import lru_cache
from itertools import product
import json
from math import isqrt
import time
import fresh_h9_quartic_mixed_20261013_independent as old
import fresh_h9_glued_host_20261013_independent as a
from fresh_h9_glued_host_20261013_independent import *

ORIGINAL_MUL=a.mul


@lru_cache(maxsize=300000)
def fastmul(x,y):
    if x==Z or y==Z:return Z
    if all(z==0 for z in x[1:]):return scale(y,x[0])
    if all(z==0 for z in y[1:]):return scale(x,y[0])
    return ORIGINAL_MUL(x,y)


a.mul=fastmul
a.alphabet=old.balanced
mul=fastmul


def rational_sqrt(x):
    if x<0:return None
    n,d=isqrt(x.numerator),isqrt(x.denominator)
    return Q(n,d)if n*n==x.numerator and d*d==x.denominator else None


def sadd(x,y):return x[0]+y[0],x[1]+y[1]
def sscale(x,c):return x[0]*c,x[1]*c
def sdiv(x,y):
    norm=y[0]*y[0]-2*y[1]*y[1]
    return sscale(silver_mul(x,(y[0],-y[1])),1/norm)


def silver_roots(x):
    aa,bb=x;out=set()
    if not bb:
        r=rational_sqrt(aa)
        if r is not None:out.add((r,Q(0)));out.add((-r,Q(0)))
        r=rational_sqrt(aa/2)
        if r is not None:out.add((Q(0),r));out.add((Q(0),-r))
    else:
        n=rational_sqrt(aa*aa-2*bb*bb)
        if n is None:return ()
        for t in(n,-n):
            r=rational_sqrt((aa+t)/2)
            if r:
                s=bb/(2*r);out.add((r,s));out.add((-r,-s))
    assert all(silver_mul(y,y)==x for y in out)
    return tuple(out)


def field_sqrt(x):
    """Complete square test using the quadratic tower over Q(sqrt2)."""
    aa,bb,cc,dd=x;u=(aa+2*cc,cc);v=(bb+2*dd,dd);d=(Q(2),Q(1))
    candidates=[]
    if v==(0,0):
        candidates.extend((p,(Q(0),Q(0)))for p in silver_roots(u))
        candidates.extend(((Q(0),Q(0)),q)for q in silver_roots(sdiv(u,d)))
    else:
        norm=sadd(silver_mul(u,u),sscale(silver_mul(d,silver_mul(v,v)),-1))
        for delta in silver_roots(norm):
            for p in silver_roots(sscale(sadd(u,delta),Q(1,2))):
                if p!=(0,0):candidates.append((p,sdiv(v,sscale(p,2))))
    for p,q in candidates:
        r=(p[0]-2*p[1],q[0]-2*q[1],p[1],q[1])
        assert sq(r)==x
        if sign(r)>=0:return r
    return None


def compat(r,q,cov,norms):
    for c,n in zip(cov,norms):
        p=dot(r,c)
        if sign(p)>0:return False
        ratio=divide(scale(sq(p),4),mul(q,n))
        if sign(sub(ratio,scale(O,4)))<0 and ratio not in ALLOWED:return False
    return True


def geometry_check(data,blocks):
    host=Host(blocks);g=host.g;walls=tuple(vec(r)for r in data['all_walls'])
    cov=[mv(g,r)for r in walls];norms=[dot(r,c)for r,c in zip(walls,cov)]
    center=vec(data['strict_interior']);assert sign(host.pair(center,center))<0
    for i,(r,c,q)in enumerate(zip(walls,cov,norms)):
        assert all(map(integral,r))and not all(integral(divide(x,T))for x in r)
        assert q in POOL and all(sign(q,e)>0 for e in range(4))
        assert all(integral(divide(scale(x,2),q))for x in c)
        assert compat(r,q,cov[:i],norms[:i])and sign(dot(center,c))<0
    geo=data['geometry'];assert geo['status']=='EXACT_STRICT_SPACELIKE_ESCAPE'
    for name,normkey,strict in [('vector','norm',False),('strict_vector','strict_norm',True)]:
        v=vec(geo[name]);norm=host.pair(v,v)
        assert norm==dec(geo[normkey])and sign(norm)>0 and sign(v[0])>0
        assert all(sign(dot(v,c))<(0 if strict else 1)for c in cov)
    assert host.pair(vec(geo['strict_vector']),center)==dec(geo['strict_center_pairing'])
    assert sign(dec(geo['strict_center_pairing']))<0
    return {'walls':len(walls),'all_pairs':len(walls)*(len(walls)-1)//2,
            'strict_interior':True,'strict_escape':True,'escape_norm':geo['norm']}


def orbit(r):
    found=set()
    for i,j,k in product((0,1),repeat=3):
        aa=r[1:5][::-1]if i else r[1:5];bb=r[5:9][::-1]if j else r[5:9]
        if k:aa,bb=bb,aa
        found.add((r[0],*aa,*bb,r[9]))
    return found


def replay(source,face=False):
    started=time.monotonic();data=json.loads(source.read_text());host=Host(2)
    for i in range(20):
        z=(Q(i-7),Q(i%3-1),Q(i%5-2),Q(i%7-3))
        expected=z if sign(z)>=0 else neg(z)
        assert field_sqrt(sq(z))==expected
    assert field_sqrt(T)is None and field_sqrt(ALPHA)is None
    base=json.loads((ROOT/data['source']).read_text());walls=tuple(vec(r)for r in base['all_walls'])
    cov=[mv(host.g,r)for r in walls];norms=[dot(r,c)for r,c in zip(walls,cov)];oldkeys={old.key(r)for r in walls}
    controls=host.controls(data['height_bound']);assert len(controls)==data['control_count']
    if face:
        saved={(dec(r['height']),dec(r['norm']),dec(r['time']),tuple(r['covector_indices'])):r for r in data['levels']}
        assert len(saved)==3*len(controls)
    else:
        saved={(dec(r['height']),dec(r['norm']),dec(r['time'])):r for r in data['levels']}
        assert set(saved)==set(controls)
    dictionary={};root_count=0
    for index,(h,q,x)in enumerate(controls):
        energy=add(q,mul(ALPHA,sq(x)))
        for j,k in(((0,2),(1,2),(2,3))if face else((0,1),)):
            if face:
                aa,bb,cc=FI[j][j],FI[k][k],FI[j][k]
                determinant=sub(mul(aa,bb),sq(cc))
                # Reverse the producer: enumerate v, then complete the u square.
                vs=old.balanced((divide(energy,bb),*(divide(mul(energy,aa),determinant)for _ in range(3))))
            else:
                aa,bb,cc=scale(O,2),scale(O,6),Z
                vs=old.balanced((divide(energy,bb),)*4)
                square_roots={sq(z):z for z in old.balanced((divide(energy,aa),)*4)}
            roots=[]
            for v in vs:
                if face:
                    square=sub(mul(aa,energy),mul(determinant,sq(v)))
                    root=field_sqrt(square)
                    if root is None:continue
                    u=divide(sub(root,mul(cc,v)),aa)
                else:
                    square=divide(sub(energy,mul(bb,sq(v))),aa)
                    if square not in square_roots:continue
                    u=square_roots[square]
                if sign(u)<=0 or sign(v)<=0 or not integral(u):continue
                if face:r=(x,*(add(mul(row[j],u),mul(row[k],v))for row in FI),Z)
                else:r=(x,*(mul(row[0],u)for row in FI),*(mul(row[1],v)for row in FI),Z)
                if all(integral(divide(z,T))for z in r):continue
                actual=Host(1).pair(r,r)if face else host.pair(r,r)
                assert actual==q and all(map(integral,r))
                roots.append(r)
                rr=(r[0],*r[1:5],*(Z,)*5)if face else r
                for image in orbit(rr):
                    key=old.key(image)
                    if key in oldkeys or not compat(image,q,cov,norms):continue
                    dictionary[key]=(image,h,q)
            row=saved[(h,q,x,(j,k))if face else(h,q,x)]
            assert set(roots)=={vec(r)for r in row['roots']}
            root_count+=len(roots)
        if index%100==0:print({'control':index,'face':face,'elapsed':round(time.monotonic()-started,2)},flush=True)
    rows=data['new_compatible_dictionary'if face else'compatible_dictionary']
    archived={old.key(vec(r['vector'])):r for r in rows}
    assert archived.keys()==dictionary.keys()
    for key,(r,h,q)in dictionary.items():
        assert(vec(archived[key]['vector']),dec(archived[key]['height']),dec(archived[key]['norm']))==(r,h,q)
    merged=dict(dictionary)
    if face:
        prior=json.loads((ROOT/data['prior_source']).read_text())
        merged.update({old.key(vec(r['vector'])):(vec(r['vector']),dec(r['height']),dec(r['norm']))for r in prior['compatible_dictionary']})
        assert merged.keys()=={old.key(vec(r['vector']))for r in data['combined_dictionary']}
    bank=data['combined_dictionary'if face else'compatible_dictionary']
    chosen=[];rejected=[];selected=list(walls);sc=list(cov);sn=list(norms)
    for i,row in enumerate(bank):
        r=vec(row['vector']);q=dec(row['norm'])
        if compat(r,q,sc,sn):chosen.append(i);selected.append(r);sc.append(mv(host.g,r));sn.append(q)
        else:rejected.append(i)
    assert chosen==data['selected_indices']and rejected==data['rejected_indices']
    assert selected==list(map(vec,data['all_walls']))
    geometry=geometry_check(data,2)
    result={'status':'PASS_COMPLETE_INDEPENDENT_TARGETED_QUARTIC_CENSUS','source':str(source.relative_to(ROOT)),
            'controls':len(controls),'slice_controls':len(saved),'primitive_slice_roots':root_count,
            'new_compatible_dictionary':len(dictionary),'combined_dictionary':len(merged),
            'selected':len(chosen),'rejected':len(rejected),'H9_geometry':geometry,
            'arithmetic':'Standalone Fraction quartic arithmetic; independent greedy unit balancing and reversed completed-square enumeration. Face replay uses a complete quadratic-tower field square-root test instead of a square alphabet. Cached exact products only.',
            'square_root_regressions':22,
            'scope':'Complete bounded slice census with both spatial weights strictly positive, plus all selected root/angle/interior/escape replay. Boundary one-weight planes are not enumerated here. No full-host nonreflectivity claim.'}
    if face:
        fg=data['forced_H5_face'];expected={(*r[:5],r[9])for r in selected if all(z==Z for z in r[5:9])}
        assert expected=={vec(r)for r in fg['all_walls']}
        result['H5_geometry']=geometry_check(fg,1)
    result['seconds']=round(time.monotonic()-started,3)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--face',action='store_true');p.add_argument('--verify',action='store_true');p.add_argument('--geometry-only',action='store_true');args=p.parse_args()
    suffix='face_target'if args.face else'target'
    source=ROOT/f'results/fresh_h9_quartic_mixed_20261013_{suffix}.json'
    out=ROOT/f'results/fresh_h9_quartic_mixed_20261013_{suffix}_independent.json'
    if args.geometry_only:
        data=json.loads(source.read_text());checks={'H9':geometry_check(data,2)}
        if args.face:checks['H5']=geometry_check(data['forced_H5_face'],1)
        bad=copy.deepcopy(data);bad['all_walls'][0][0]='1' if isinstance(bad['all_walls'][0][0],str)else['1','0','0','0']
        try:geometry_check(bad,2)
        except AssertionError:pass
        else:raise AssertionError('Forged normal accepted')
        print({'status':'PASS_EXACT_GEOMETRY_REPLAY','checks':checks,'forged_normal_rejected':True})
        raise SystemExit(0)
    result=replay(source,args.face)
    if args.verify:
        saved=json.loads(out.read_text());assert{k:v for k,v in saved.items()if k!='seconds'}=={k:v for k,v in result.items()if k!='seconds'}
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result)
