#!/usr/bin/env python3
"""Bounded ternary slice search freeing the scalar/E8 ratio of an H9 escape."""
from __future__ import annotations
import argparse
import json
from math import isqrt,lcm
import time
from coxeter_polytopes.exact.golden import Golden as G,PHI,ONE as O,ZERO as Z
from coxeter_polytopes.geometric.golden_fusion import golden_square_root
from scripts.fresh_multibasepoint_h9_merge import ROOT,read,matrix,mv,dot,enc,encm,coxeter_pair_value
from scripts.fresh_multibasepoint_h9_escape_plane import ri_sign,scale_bound_coefficients

E8_DIRECTION=(7,14,20,16,12,8,4,10)
def ri_mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def ri_sum(*xs):return sum(x[0] for x in xs),sum(x[1] for x in xs)
def ri_scale(x,n):return x[0]*n,x[1]*n

def alphabet(limit,conjugate_cap):
    k=1
    while (G(k*k)-conjugate_cap).sign(1)<0:k+=1
    ca,cb,den=scale_bound_coefficients(conjugate_cap);out=[]
    for b in range(-k-1,(limit+k)//2+3):
        ab=abs(b);floor=(isqrt(5*ab*ab)-ab)//2
        if b<0:floor=-floor-1
        for a in range(floor-k-1,floor+k+3):
            if ri_sign(a,b)<0 or ri_sign(limit-a,-b)<0:continue
            aa=a*a+b*b;bb=2*a*b+b*b
            if ri_sign(ca-den*(aa+bb),cb+den*bb)<0:continue
            out.append(G(a,b))
    return tuple(sorted(out))

def run(t_bound,u_bound):
    start=time.monotonic();source='results/fresh_multibasepoint_h9_2phi_targeted_rounds_20261013.json';d=read(source)
    gram=matrix(d['gram']);walls=matrix(d['all_walls']);normals=tuple(mv(gram,r) for r in walls)
    old_norms=tuple(dot(r,n) for r,n in zip(walls,normals));alpha=G.from_json(d['alpha'])
    v=tuple(map(G.from_json,d['hyperideal_witness']['vector']));inside=tuple(map(G.from_json,d['strict_interior']))
    z=(Z,*(G(k) for k in E8_DIRECTION),Z);scalar=(Z,*(Z for _ in range(8)),O)
    assert alpha==2*PHI and dot(z,mv(gram,z))==14 and dot(z,mv(gram,scalar))==0
    assert all(v[i]==v[1]*G(E8_DIRECTION[i-1])/7 for i in range(1,9))
    cone_coefficients=tuple((-alpha*r[0],dot(z,n),r[9]) for r,n in zip(walls,normals))
    av=dot(z,mv(gram,v));bv=v[9]
    assert av.sign()>0 and bv.sign()>0
    escape_coefficients=(av*av-14*alpha,2*av*bv,bv*bv-alpha)
    escape_den=lcm(*(c.denominator for k in escape_coefficients for c in (k.a,k.b)))
    c11,c12,c22=tuple((int(k.a*escape_den),int(k.b*escape_den)) for k in escape_coefficients)
    results=[];cuts=[]
    for q in (G(1),G(2)):
        ts=alphabet(t_bound,q/14);us=alphabet(u_bound,q)
        ts=tuple(t for t in ts if t.sign()>0)
        us=tuple(u for u in us if int(u.a)%2==int(q.a)%2 and int(u.b)%2==0)
        counters={'t_count':len(ts),'u_count_after_parity':len(us),'pairs':0,'conjugate_admissible':0,
            'potential_escape_cuts':0,'integral_square_candidates':0,'integral_norm_solutions_after_escape_prefilter':0,'cuts_old_escape':0,
            'cone_compatible':0,'coxeter_compatible':0,'outside_previous_plane':0}
        u_entries=[(u,(int(u.a),int(u.b)),ri_mul((int(u.a),int(u.b)),(int(u.a),int(u.b)))) for u in us]
        constant=(int((alpha*q).a*escape_den),int((alpha*q).b*escape_den))
        for t in ts:
            ti=(int(t.a),int(t.b));tt=ri_mul(ti,ti);energy_t=ri_scale(tt,14)
            first=ri_mul(c11,tt)
            for u,ui,uu in u_entries:
                counters['pairs']+=1;ea,eb=ri_sum(energy_t,uu)
                if ri_sign(int(q.a)-ea-eb,eb)<0 or ri_sign(ea-int(q.a),eb)<=0:continue
                counters['conjugate_admissible']+=1
                # L>0, so L-alpha*x>0 is exactly L^2>alpha*(energy-q).
                gate=ri_sum(first,ri_mul(c12,ri_mul(ti,ui)),ri_mul(c22,uu),constant)
                if ri_sign(*gate)<=0:continue
                counters['potential_escape_cuts']+=1;square=(G(ea,eb)-q)/alpha
                if not square.is_integral:continue
                counters['integral_square_candidates']+=1;x=golden_square_root(square)
                if x is None or not x.is_integral:continue
                if x.sign()<0:x=-x
                counters['integral_norm_solutions_after_escape_prefilter']+=1
                r=(x,*(t*G(k) for k in E8_DIRECTION),u)
                if dot(r,mv(gram,v)).sign()<=0:continue
                counters['cuts_old_escape']+=1
                products=tuple(a*x+b*t+c*u for a,b,c in cone_coefficients)
                if any(p.sign()>0 for p in products):continue
                counters['cone_compatible']+=1
                if not all(coxeter_pair_value(4*p*p/(q*old_q)) for p,old_q in zip(products,old_norms)):continue
                counters['coxeter_compatible']+=1
                # The old rank-two slice requires u/t = v_scalar/(v_E8/7).
                if u*v[1]==7*t*v[9]:continue
                counters['outside_previous_plane']+=1
                assert r not in walls and dot(r,mv(gram,r))==q
                assert all((2*y/q).is_integral for y in mv(gram,r))
                assert dot(r,mv(gram,inside)).sign()<0
                cuts.append({'root':enc(r),'norm':q.to_json(),'x':x.to_json(),'t':t.to_json(),'u':u.to_json(),
                    'old_basepoint_height':(alpha*x*x/q).to_json(),
                    'escape_pairing':dot(r,mv(gram,v)).to_json(),'wall_pairings':enc(products),
                    'outside_previous_plane_witness':(u*v[1]-7*t*v[9]).to_json()})
        results.append({'norm':q.to_json(),**counters})
        print(json.dumps({'norm':str(q),**counters}),flush=True)
    cuts.sort(key=lambda c:G.from_json(c['old_basepoint_height']))
    out={'source':source,'alpha':alpha.to_json(),'dimension':9,'rank3_gram':[[(-alpha).to_json(),Z.to_json(),Z.to_json()],
        [Z.to_json(),G(14).to_json(),Z.to_json()],[Z.to_json(),Z.to_json(),O.to_json()]],
        'e8_direction':list(E8_DIRECTION),'t_bound':t_bound,'u_bound':u_bound,
        'parameter_scope':'0<t<=t_bound, 0<=u<=u_bound in Z[phi]; exact conjugate norm bounds; q in {1,2}.',
        'completed_parameter_box':True,'counters':results,'cuts':cuts,
        'targeted_escape':d['hyperideal_witness'],'seconds':round(time.monotonic()-start,3),
        'scope':'Exact bounded saturated rank-three slice search, not ambient root exhaustion or a volume conclusion.'}
    output=ROOT/'results/fresh_rank3_h9_search_20261013.json'
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'cuts':len(cuts),'seconds':out['seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--t-bound',type=int,default=200)
    p.add_argument('--u-bound',type=int,default=400);a=p.parse_args();run(a.t_bound,a.u_bound)
