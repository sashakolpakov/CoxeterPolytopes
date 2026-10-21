#!/usr/bin/env python3
"""Exact bounded symmetric ternary-slice completion of H3+A6 H9 cones."""
import argparse
from fractions import Fraction as F
from math import isqrt
import json
from pathlib import Path
import time
from scripts.fresh_h3_blocks_h9_a6_search import G,PHI,O,Z,INVERSES
from scripts.fresh_multibasepoint_h9_merge import ROOT,read,matrix,mv,dot,enc,encm,coxeter_pair_value
from scripts.fresh_rank3_h9_search import alphabet,ri_sign
from coxeter_polytopes.geometric.golden_fusion import golden_square_root

def mul(a,b):return (a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])

def setup(name):
    path='results/fresh_h3_blocks_h9_a6_'+name+'_h30_20261013.json'
    d=read(path);gram=matrix(d['gram']);walls=matrix(d['all_walls']);alpha=G.from_json(d['alpha'])
    inside=tuple(map(G.from_json,d['strict_interior']));v=tuple(map(G.from_json,d['hyperideal_witness']['vector']))
    reversal=(0,1,2,3,9,8,7,6,5,4)
    target=tuple((a+v[reversal[i]])/2 for i,a in enumerate(v))
    h=(Z,*(2*INVERSES[0][i][0] for i in range(3)),*((Z,)*6))
    a=(Z,Z,Z,Z,*(INVERSES[1][i][2]+INVERSES[1][i][3] for i in range(6)))
    assert all(c.is_integral for c in (*h,*a))
    assert dot(h,mv(gram,h))==6*PHI**2 and dot(a,mv(gram,a))==6 and dot(h,mv(gram,a))==0
    assert dot(target,mv(gram,target)).sign()>0 and all(dot(target,mv(gram,w)).sign()<=0 for w in walls)
    return path,d,gram,walls,alpha,inside,target,h,a

def run(name,bound,output):
    start=time.monotonic();path,d,gram,walls,alpha,inside,target,h,a=setup(name)
    normals=tuple(mv(gram,r) for r in walls);norms=tuple(dot(r,n) for r,n in zip(walls,normals))
    rows=tuple((-alpha*r[0],dot(h,n),dot(a,n)) for r,n in zip(walls,normals))
    pool=(1,2,4) if name=='phi' else (2,4);factor=1 if name=='phi' else 2
    cuts=[];censuses=[]
    for q in pool:
        ss=alphabet(bound,G(q)/(6*PHI**2));us=alphabet(bound,G(F(q,6)))
        uu=[(u,(int(u.a),int(u.b)),mul((int(u.a),int(u.b)),(int(u.a),int(u.b)))) for u in us]
        count={'norm':q,'s_alphabet':len(ss),'u_alphabet':len(us),'pairs':len(ss)*len(us),
            'conjugate_admissible':0,'integral_square_candidates':0,'integral_norm_solutions':0,
            'primitive_reflections':0,'acute_compatible':0,'coxeter_compatible':0,'cut_old_target':0}
        for s in ss:
            sa,sb=mul((int(s.a),int(s.b)),(int(s.a),int(s.b)))
            for u,ui,(ua,ub) in uu:
                ea,eb=6*(sa+sb+ua)-q,6*(sa+2*sb+ub)
                if ri_sign(ea,eb)<=0 or ri_sign(ea+eb,-eb)>0:continue
                count['conjugate_admissible']+=1
                if ea%factor or eb%factor:continue
                xa,xb=(eb-ea)//factor,ea//factor
                norm=xa*xa+xa*xb-xb*xb
                if norm<0 or isqrt(norm)**2!=norm:continue
                count['integral_square_candidates']+=1
                x=golden_square_root(G(xa,xb))
                if x is None or not x.is_integral:continue
                if x.sign()<0:x=-x
                count['integral_norm_solutions']+=1
                r=(x,*(s*h[i]+u*a[i] for i in range(1,10)))
                assert dot(r,mv(gram,r))==q and all(c.is_integral for c in r)
                if q==4 and all((c/2).is_integral for c in r):continue
                if not all((2*c/q).is_integral for c in mv(gram,r)):continue
                count['primitive_reflections']+=1
                products=tuple(c*x+b*s+e*u for c,b,e in rows)
                if any(p.sign()>0 for p in products):continue
                count['acute_compatible']+=1
                if not all(coxeter_pair_value(4*p*p/(q*n)) for p,n in zip(products,norms)):continue
                if r in walls:continue
                count['coxeter_compatible']+=1
                target_pair=dot(r,mv(gram,target));count['cut_old_target']+=int(target_pair.sign()>0)
                assert dot(r,mv(gram,inside)).sign()<0
                cuts.append({'root':enc(r),'norm':q,'parameters':enc((x,s,u)),
                    'target_pairing':target_pair.to_json(),'height':(alpha*x*x/q).to_json()})
        censuses.append(count);print(json.dumps({'alpha':name,**count}),flush=True)
    payload={'source':path,'alpha':alpha.to_json(),'gram':d['gram'],'basis':encm(((O,*((Z,)*9)),h,a)),
        'target':enc(target),'target_norm':dot(target,mv(gram,target)).to_json(),'strict_interior':d['strict_interior'],
        'norm_pool':pool,'bounds':{'s_primary_max':bound,'u_primary_max':bound},'censuses':censuses,
        'cuts':cuts,'seconds':round(time.monotonic()-start,3),
        'scope':'All new primitive old-compatible Coxeter normals in the specified symmetric ternary slice and parameter box, including those not cutting the old target. Not a complete height search or full-slice nonreflectivity claim.'}
    output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'alpha':name,'cuts':len(cuts),'seconds':payload['seconds']}),flush=True)
    return payload

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--alpha',choices=('phi','2phi'),required=True)
    p.add_argument('--bound',type=int,default=1000);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.alpha,a.bound,a.output)
