#!/usr/bin/env python3
"""Independent edge-profile/time census and unrelated polynomial square test."""
import argparse,hashlib,json
from functools import lru_cache
from itertools import product
from math import isqrt
import sympy as sp
import fresh_h9_new_host_20261013_independent as a
from fresh_h9_glued_host_20261013_structure import norm

SOURCE=a.ROOT/'results/fresh_h9_new_host_20261013_edge.json'
OUT=a.ROOT/'results/fresh_h9_new_host_20261013_edge_independent.json'

@lru_cache(None)
def balanced(caps):
    if any(a.sign(c,e)<0 for e,c in enumerate(caps)):return ()
    if any(a.sign(c,e)==0 for e,c in enumerate(caps)):return (a.Z,)
    cub=a.mul(a.sq(a.T),a.T)
    units=[a.sub(a.T,a.O),a.sub(a.sub(cub,a.scale(a.T,3)),a.O),a.sub(a.sub(a.scale(a.T,3),cub),a.O)]
    units=units+[a.inv(u)for u in units]
    unit=a.O
    def score(u):return sum(a.ceil_real(a.divide(c,a.sq(u)),e)for e,c in enumerate(caps))
    while True:
        proposals=[a.mul(unit,u)for u in units];best=min(proposals,key=score)
        if score(best)>=score(unit):break
        unit=best
    if a.sign(unit)<0:unit=a.neg(unit)
    assert a.integral(unit)and a.integral(a.inv(unit))
    transformed=tuple(a.divide(c,a.sq(unit))for c in caps)
    return tuple(a.mul(unit,x)for x in a.alphabet(transformed))

def run(verify=False):
    data=json.loads(SOURCE.read_text());face=json.loads((a.ROOT/data['face_source']).read_text())
    for p,h in data['source_sha256'].items():assert hashlib.sha256((a.ROOT/p).read_bytes()).hexdigest()==h
    host=a.Host(1);walls=[a.vec(r)for r in face['all_walls']];edge=[walls[i]for i in data['edge_indices']]
    gi=a.minverse(tuple(tuple(host.pair(r,s)for s in edge)for r in edge));direction=a.vec(data['direction'])
    assert host.pair(direction,direction)==a.scale(a.O,6)
    saved={(a.dec(p['norm']),a.vec(p['fixed_magnitudes'])):p for p in data['profiles']}
    profiles=set();quadratics=0;norm_rejections=0;residual=[]
    for q in a.POOL:
        alph=[]
        for e in edge:
            values=set()
            for ratio in a.ALLOWED:
                square=a.scale(a.mul(a.mul(q,host.pair(e,e)),ratio),a.Q(1,4))
                values.update(x for x in a.alphabet((square,)*4)if a.sq(x)==square)
            alph.append(values)
        for fixed in product(*alph):
            cf=a.mv(gi,tuple(a.neg(z)for z in fixed))
            r0=tuple(a.dot(cf,tuple(e[j]for e in edge))for j in range(6))
            R=a.sub(q,host.pair(r0,r0))
            if any(a.sign(R,e)<=0 for e in range(4)):continue
            key=(q,fixed);profiles.add(key);row=saved[key]
            assert a.vec(row['fixed_component'])==r0 and a.dec(row['residual_norm'])==R
            caps=(a.divide(a.scale(q,data['height_bound']),a.ALPHA),*(a.divide(a.neg(R),a.ALPHA)for _ in range(3)))
            assert caps==a.vec(row['time_caps'])
            times=set(balanced(caps));assert times==set(a.vec(row['time_coordinates']))
            for x in times:
                if a.sign(x)<=0:continue
                quadratics+=1;value=a.scale(a.add(R,a.mul(a.ALPHA,a.sq(x))),a.Q(1,6));n=norm(value)
                if n<0 or isqrt(n.numerator)**2!=n.numerator or isqrt(n.denominator)**2!=n.denominator:
                    norm_rejections+=1
                else:residual.append(value)
    assert profiles==set(saved)and quadratics==data['counts']['quadratics']
    # Exact factorization in the same abstract number field is independent of
    # the producer's quadratic-tower square routine; embedding choice does not
    # affect existence of a field square root.
    y,x=sp.symbols('y x');th=sp.CRootOf(y**4-4*y**2+2,3);records=[]
    for value in residual:
        expression=sum(sp.Rational(c.numerator,c.denominator)*th**j for j,c in enumerate(value))
        coefficient,factors=sp.factor_list(x*x-expression,x,extension=th)
        assert len(factors)==1 and factors[0][1]==1 and sp.degree(factors[0][0],x)==2
        records.append({'square_candidate':a.enc(value),'irreducible_degree':2})
    assert not data['integral_candidates']and not data['compatible_candidates']
    result={'status':'PASS_COMPLETE_BOUNDED_EDGE_REPLAY','source':str(SOURCE.relative_to(a.ROOT)),
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'height_bound':data['height_bound'],
        'profiles':len(profiles),'quadratics':quadratics,'rational_norm_nonsquare_rejections':norm_rejections,
        'number_field_irreducibility_checks':records,'integral_roots':0,
        'scope':'Independent direct trace enumeration after exact greedy integral-unit balancing, all finite coupling profiles and time coordinates. 374 square candidates have nonsquare rational norm; the remaining17 quadratics are irreducible by exact SymPy number-field factorization. This proves only the stated bounded edge-proposal failure.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],len(profiles),quadratics,norm_rejections,len(residual),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
