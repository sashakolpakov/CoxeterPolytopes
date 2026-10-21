#!/usr/bin/env python3
"""Bounded H3-neighbor rank-four slice, beyond the mod-3 closed slice."""
import argparse
from fractions import Fraction as F
from math import isqrt
import json
from pathlib import Path
import time
from scripts.fresh_h3_a6_target_search import *

def run(name,bound,output):
    start=time.monotonic();path,d,gram,walls,alpha,inside,target,h,a=setup(name)
    k=(Z,*(INVERSES[0][i][1] for i in range(3)),*((Z,)*6))
    assert all(c.is_integral for c in k)
    assert dot(k,mv(gram,k))==2*PHI**2 and dot(h,mv(gram,k))==2*PHI**3
    normals=tuple(mv(gram,r) for r in walls);norms=tuple(dot(r,n) for r,n in zip(walls,normals))
    rows=tuple((-alpha*r[0],dot(h,n),dot(k,n),dot(a,n)) for r,n in zip(walls,normals))
    pool=(1,2,4) if name=='phi' else (2,4);factor=1 if name=='phi' else 2
    cuts=[];censuses=[]
    for q in pool:
        ss=alphabet(bound,G(F(q,2)));ts=alphabet(bound,G(F(3*q,2)));us=alphabet(bound,G(F(q,6)))
        if q==4:
            ts=tuple(t for t in ts if (t/2).is_integral);us=tuple(u for u in us if (u/2).is_integral)
        uu=[(u,mul((int(u.a),int(u.b)),(int(u.a),int(u.b)))) for u in us]
        tt=[(t,(int(t.a),int(t.b)),mul((int(t.a),int(t.b)),(int(t.a),int(t.b)))) for t in ts if t.sign()>0]
        count={'norm':q,'s_alphabet':len(ss),'t_positive_alphabet':len(tt),'u_alphabet':len(us),'triples':len(ss)*len(tt)*len(us),
            'conjugate_admissible':0,'norm_square_candidates':0,'integral_norm_solutions':0,
            'primitive_reflections':0,'acute_compatible':0,'coxeter_compatible':0,'cut_old_target':0}
        for s in ss:
            si=(int(s.a),int(s.b));sa,sb=mul(si,si)
            for t,ti,(ta,tb) in tt:
                sta,stb=mul(si,ti)
                ha,hb=6*(sa+sb)+4*(sta+2*stb)+2*(ta+tb)-q,6*(sa+2*sb)+4*(2*sta+3*stb)+2*(ta+2*tb)
                for u,(ua,ub) in uu:
                    ea,eb=ha+6*ua,hb+6*ub
                    if ri_sign(ea,eb)<=0 or ri_sign(ea+eb,-eb)>0:continue
                    count['conjugate_admissible']+=1
                    if ea%factor or eb%factor:continue
                    xa,xb=(eb-ea)//factor,ea//factor
                    norm=xa*xa+xa*xb-xb*xb
                    if norm<0 or isqrt(norm)**2!=norm:continue
                    count['norm_square_candidates']+=1
                    x=golden_square_root(G(xa,xb))
                    if x is None or not x.is_integral:continue
                    if x.sign()<0:x=-x
                    count['integral_norm_solutions']+=1
                    r=(x,*(s*h[i]+t*k[i]+u*a[i] for i in range(1,10)))
                    assert dot(r,mv(gram,r))==q and all(c.is_integral for c in r)
                    if q==4 and all((c/2).is_integral for c in r):continue
                    if not all((2*c/q).is_integral for c in mv(gram,r)):continue
                    count['primitive_reflections']+=1
                    products=tuple(c*x+b*s+f*t+e*u for c,b,f,e in rows)
                    if any(p.sign()>0 for p in products):continue
                    count['acute_compatible']+=1
                    if not all(coxeter_pair_value(4*p*p/(q*n)) for p,n in zip(products,norms)):continue
                    if r in walls:continue
                    count['coxeter_compatible']+=1
                    target_pair=dot(r,mv(gram,target));count['cut_old_target']+=int(target_pair.sign()>0)
                    assert dot(r,mv(gram,inside)).sign()<0
                    cuts.append({'root':enc(r),'norm':q,'parameters':enc((x,s,t,u)),
                        'target_pairing':target_pair.to_json(),'height':(alpha*x*x/q).to_json()})
        censuses.append(count);print(json.dumps({'alpha':name,**count}),flush=True)
    payload={'source':path,'alpha':alpha.to_json(),'gram':d['gram'],'basis':encm(((O,*((Z,)*9)),h,k,a)),
        'target':enc(target),'target_norm':dot(target,mv(gram,target)).to_json(),'strict_interior':d['strict_interior'],
        'norm_pool':pool,'bounds':{'s_primary_max':bound,'t_primary_max':bound,'u_primary_max':bound},'censuses':censuses,
        'cuts':cuts,'seconds':round(time.monotonic()-start,3),
        'scope':'All new primitive old-compatible Coxeter normals in span(time,omega1(H3),omega2(H3),omega3(A6)+omega4(A6)), with 0<=s,u<=bound and 0<t<=bound. The t=0 slice is globally root-free by the separate mod-3 proof. Not a complete ambient height search.'}
    output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'alpha':name,'cuts':len(cuts),'seconds':payload['seconds']}),flush=True)
    return payload

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--alpha',choices=('phi','2phi'),required=True)
    p.add_argument('--bound',type=int,default=100);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.alpha,a.bound,a.output)
