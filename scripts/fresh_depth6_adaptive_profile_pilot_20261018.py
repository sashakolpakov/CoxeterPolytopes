#!/usr/bin/env python3
"""Adaptive small-core diagnostic. Not a universal geometric theorem."""
import json
from itertools import product
from pathlib import Path
import sympy as S
import fresh_depth8_finite7_profile_transfer_20261018 as P

ROOT=P.ROOT
TAG='fresh_depth6_adaptive_profile_pilot_20261018'

def run():
    doc=json.loads((ROOT/'results/fresh_depth6_depth7_finite7_profile_pilot_20261018.json').read_text())
    ps=[{int(i):S.sympify(q) for i,q in x.items()} for x in doc['rows'][0]['domain']['accepted_profiles']]
    t,u,T,U=S.symbols('t u T U',nonnegative=True)
    rows=[];bad=[];discarded=0
    def inner(x,y,n,shift):
        return S.expand(sum(v*w*P.inv(n,i+shift,j+shift) for i,v in x.items() for j,w in y.items()))
    for ci,c in enumerate(ps):
        for di,rr in enumerate(ps):
            d={19-i:q for i,q in rr.items()}
            if any(abs(i-j)+2<=6 and not any(min(i,j)<k<max(i,j) for k in(6,13)) for i in c for j in d):continue
            x={**c,6:t};y={**d,13:u}
            oq=inner(x,x,20,0);orr=inner(y,y,20,0)
            for tc,ud in product((0,1),repeat=2):
                clow=P.field(oq.subs(t,tc)-2)['sign']<=0
                dlow=P.field(orr.subs(u,ud)-2)['sign']<=0
                if clow and (tc or c!={1:1}) or dlow and (ud or d!={18:1}):discarded+=1;continue
                for cp,dp in product(('extend','isolated') if clow else ('none',),('extend','isolated') if dlow else ('none',)):
                    shift=int(cp=='extend');n=20+shift+int(dp=='extend')
                    q=inner(x,x,n,shift)+S.Rational(int(cp=='isolated'),2)
                    r=inner(y,y,n,shift)+S.Rational(int(dp=='isolated'),2)
                    ss=inner(x,y,n,shift)
                    q0=P.field(q.subs(t,tc)-2);r0=P.field(r.subs(u,ud)-2)
                    assert q0['sign']>0 and r0['sign']>0
                    p=S.Poly(S.expand(((q-2)*(r-2)-ss**2).subs({t:1+T if tc else 0,u:1+U if ud else 0})),T,U)
                    coeff=[dict(powers=list(m),value=P.field(v)) for m,v in p.terms()]
                    const=P.field(p.eval({T:0,U:0}))
                    row=dict(c=ci,d=di,t=tc,u=ud,c_patch=cp,d_patch=dp,coefficients=coeff,constant=const)
                    rows.append(row)
                    if const['sign']<=0 or any(z['value']['sign']<0 for z in coeff):bad.append(row)
    out=dict(scope='Conditional adaptive extensions and isolated whites; geometric classification and contact seals require proof.',profiles=ps,branches=rows,unproved=bad,discarded_low_projection=discarded)
    out['profiles']=[{str(i):str(q)for i,q in x.items()}for x in ps]
    path=ROOT/f'results/{TAG}.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print('ADAPTIVE',len(rows),'unproved',len(bad),'discarded',discarded,flush=True)
    for r in bad:print(r,flush=True)

if __name__=='__main__':run()
