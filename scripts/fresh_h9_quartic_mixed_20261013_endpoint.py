#!/usr/bin/env python3
"""Finite-label edge profiles and a constructive relative-Pell endpoint cap."""
from fresh_h9_quartic_mixed_20261013_target import *
from fresh_h9_quartic_mixed_20261013_target_independent import field_sqrt
from fractions import Fraction
from itertools import combinations

SOURCE=ROOT/'results/fresh_h9_quartic_mixed_20261013_face_target.json'
OUT=ROOT/'results/fresh_h9_quartic_mixed_20261013_endpoint.json'


def square_root(x):
    value=field_sqrt(tuple(map(Fraction,x.to_json())))
    return None if value is None else K(*value)


def finite_profiles():
    """All K-integral finite couplings to the fixed A2+A1+A1 edge."""
    rows=[];alphabets=[]
    for q in NORMS:
        ws=sorted({r for label in ALLOWED if (r:=square_root(q*label/2)) is not None and r.is_integral})
        zs=sorted({r for label in ALLOWED if (r:=square_root(q*label/4)) is not None and r.is_integral})
        alphabets.append({'norm':q.to_json(),'F4_covectors':enc(ws),'scalar_coordinates':enc(zs)})
        for a,b,c,z in product(ws,ws,ws,zs):
            qp=q-K(Fraction(2,3))*(a*a+a*b+b*b)-c*c/2-z*z
            if any(qp.sign(e)<=0 for e in range(4)):continue
            rows.append({'norm':q.to_json(),'fixed_covectors':enc((a,b,c,z)),
                         'orthogonal_norm':qp.to_json(),
                         'completed_square_shift':(2*SQRT2*a+4*SQRT2*b+3*c).to_json()})
    return alphabets,rows


def certificate(walls,g,inside,v):
    values=tuple(pair(r,g,v)for r in walls)
    assert pair(v,g,v)>0 and all(z<=0 for z in values)
    eps=O
    while True:
        strict=tuple(x+eps*y for x,y in zip(v,inside))
        if pair(strict,g,strict)>0 and pair(strict,g,inside)<0:break
        eps/=2
    assert all(pair(r,g,strict)<0 for r in walls)
    return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','vector':enc(v),'norm':pair(v,g,v).to_json(),
            'wall_pairings':enc(values),'strict_vector':enc(strict),'strict_norm':pair(strict,g,strict).to_json(),
            'strict_center_pairing':pair(strict,g,inside).to_json()}


def scaled_geometry(walls,g,inside,n):
    import numpy as np
    from scipy.optimize import linprog
    cov=tuple(mv(g,r)for r in walls);raw=np.array([[z.real()for z in c]for c in cov])
    matrix=raw/np.max(np.abs(raw),axis=1)[:,None]
    seen=set();rng=np.random.default_rng(20261013)
    for _ in range(500):
        lp=linprog(rng.normal(size=n),A_ub=matrix[:,1:],b_ub=-matrix[:,0],bounds=(None,None),method='highs')
        if not lp.success:continue
        active=tuple(map(int,np.flatnonzero(abs(matrix[:,1:]@lp.x+matrix[:,0])<1e-8)))
        if len(active)>n+2:continue
        for ids in combinations(active,n):
            if ids in seen:continue
            seen.add(ids)
            try:v=(O,*mv(inverse(tuple(tuple(cov[i][j]for j in range(1,n+1))for i in ids)),tuple(-cov[i][0]for i in ids)))
            except(StopIteration,ZeroDivisionError):continue
            if pair(v,g,v)>0 and all(dot(c,v)<=0 for c in cov):
                out=certificate(walls,g,inside,v);out['active_facets']=list(ids);return out
    return {'status':'UNDECIDED','exact_supports_tested':len(seen)}


def run():
    started=time.monotonic();source=json.loads(SOURCE.read_text());g=Host(2).g;fg=Host(1).g
    old=[tuple(map(K.from_json,r))for r in source['all_walls']]
    y=K(-17023554335,-9213081216,29061025044,15727713420)
    s=-K(7548113028,4085009334,-12885434936,-6973547136)
    assert y>0 and s>0 and y*y-6*ALPHA*s*s==O
    x=3*s;t=(y-1)/2;r=(x,*(row[2]*t+row[3]for row in FI),*(Z,)*5)
    roots=sorted(orbit(r));assert len(roots)==4
    wc=[mv(g,w)for w in old];wq=[dot(w,c)for w,c in zip(old,wc)]
    for w in roots:
        assert pair(w,g,w)==2 and all(z.is_integral for z in w)
        assert not all((z/T).is_integral for z in w)
        assert all((2*z/2).is_integral for z in mv(g,w))
        assert coxeter(w,K(2),wc,wq)
        wc.append(mv(g,w));wq.append(K(2))
    walls=old+roots;inside=tuple(map(K.from_json,source['strict_interior']))
    assert pair(inside,g,inside)<0 and all(pair(w,g,inside)<0 for w in walls)
    v=(O,*(row[2]*ALPHA*x/(3*y)for row in FI),*(Z,)*5)
    assert pair(v,g,v)==-ALPHA/y**2<0
    assert pair(r,g,v)==0 and all(pair(w,g,v)<=0 for w in walls)
    zero=[i for i,w in enumerate(walls)if pair(w,g,v)==0]
    assert len(zero)==9
    endpoint={'vector':enc(v),'norm':pair(v,g,v).to_json(),'active_facets':zero,
              'wall_pairings':enc(tuple(pair(w,g,v)for w in walls)),
              'center_pairing':pair(v,g,inside).to_json(),
              'diagram':'F4 + A2 + A2 + A1','forced_H5_diagram':'A2 + A2 + A1'}
    fw=[(*w[:5],w[9])for w in walls if not any(w[5:9])];fi=(*inside[:5],inside[9])
    oldv=tuple(map(K.from_json,source['geometry']['vector']))
    alphabets,profiles=finite_profiles()
    out={'status':'EXACT_RELATIVE_PELL_FINITE_ENDPOINT','source':str(SOURCE.relative_to(ROOT)),
         'relative_pell':{'y':y.to_json(),'s':s.to_json(),'identity':'y^2-6*(T-1)*s^2=1',
                          'time':x.to_json(),'inner_covector':t.to_json(),
                          'height':(ALPHA*x*x/2).to_json(),'discovery':'PARI unit calculation; replay needs no PARI'},
         'finite_endpoint_alphabets':alphabets,'finite_endpoint_profiles':profiles,
         'new_roots':menc(roots),'all_walls':menc(walls),'wall_count':len(walls),
         'strict_interior':enc(inside),'finite_endpoint':endpoint,
         'previous_H5_escape_cut':pair((*r[:5],r[9]),fg,tuple(map(K.from_json,source['forced_H5_face']['geometry']['vector']))).to_json(),
         'geometry':certificate(walls,g,inside,oldv),
         'forced_H5_face':{'all_walls':menc(fw),'wall_count':len(fw),'strict_interior':enc(fi),
                           'geometry':scaled_geometry(fw,fg,fi,5)},
         'seconds':round(time.monotonic()-started,3),
         'scope':'Four explicit reflective roots close the specified A2+A1+A1 edge at actual finite vertices. The finite-label profile census is complete only for K-integral norm-pool roots extending that fixed edge elliptically; it is not a complete Pell solution census. The full H9 host remains open.'}
    assert K.from_json(out['previous_H5_escape_cut'])>0
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'walls':len(walls),'profiles':len(profiles),'H5_walls':len(fw),
                      'H5_geometry':out['forced_H5_face']['geometry']['status'],'seconds':out['seconds']}),flush=True)


if __name__=='__main__':run()
