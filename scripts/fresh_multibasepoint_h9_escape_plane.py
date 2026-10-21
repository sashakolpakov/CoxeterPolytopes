#!/usr/bin/env python3
"""Target a surviving escape in the saturated plane spanned by it and p.

This is a bounded rank-two slice search, not exhaustive root enumeration
of the ambient lattice. A discovered compatible reflective cut is exact.
"""
from __future__ import annotations
import argparse
from math import isqrt,lcm
import time
import json
from sympy import factorint
from coxeter_polytopes.exact.golden import Golden as G,PHI,ONE as O,ZERO as Z
from coxeter_polytopes.geometric.golden_fusion import golden_square_root
from scripts.fresh_multibasepoint_h9_merge import ROOT,read,matrix,mv,dot,enc,coxeter_pair_value

def nearest(f):return (2*f.numerator+f.denominator)//(2*f.denominator)
def gcd(a,b):
    while b:
        z=a/b;q=G(nearest(z.a),nearest(z.b));r=a-q*b
        assert not r or abs(r.norm())<abs(b.norm())
        a,b=b,r
    return a

def ri_sign(a,b):
    x=2*a+b
    if not b:return (x>0)-(x<0)
    if not x or (x>0)==(b>0):return (b>0)-(b<0)
    return ((x>0)-(x<0))*((x*x>5*b*b)-(x*x<5*b*b))

def scale_bound_coefficients(cap):
    c=cap.conjugate();den=lcm(c.a.denominator,c.b.denominator)
    return int(c.a*den),int(c.b*den),den

def local_obstructions(beta,alpha,q):
    output=[]
    for p in factorint(abs(int(beta.norm()))):
        p=int(p)
        if p==2:continue
        for f in range(p):
            if (f*f-f-1)%p or (int(beta.a)+int(beta.b)*f)%p:continue
            a=(int(alpha.a)+int(alpha.b)*f)%p
            value=(int(q.a)+int(q.b)*f)%p
            if not a:continue
            residue=(-value*pow(a,-1,p))%p
            if residue and pow(residue,(p-1)//2,p)==p-1:
                output.append({'prime':p,'phi_residue':f,'required_square_residue':residue})
    return output

def plane_parameters(gram,v):
    denominator=lcm(*(c.denominator for x in v[1:] for c in (x.a,x.b)))
    raw=tuple(denominator*x for x in v[1:]);content=Z
    for x in raw:content=gcd(content,x)
    if content.sign()<0:content=-content
    z=tuple(x/content for x in raw);c=G(denominator)/content
    spatial=lambda y:dot((Z,*y),mv(gram,(Z,*y)))
    beta=spatial(z)
    # Balance the totally positive norm by a square of the fundamental unit.
    while (beta*beta*PHI*PHI-beta.norm()).sign()<0:
        z=tuple(PHI*x for x in z);c*=PHI;beta*=PHI*PHI
    while (beta*beta-beta.norm()*PHI*PHI).sign()>0:
        z=tuple(x/PHI for x in z);c/=PHI;beta/=PHI*PHI
    assert all(x.is_integral for x in z) and c.sign()>0 and beta==spatial(z)
    content_check=Z
    for x in z:content_check=gcd(content_check,x)
    assert abs(content_check.norm())==1
    assert tuple(x/c for x in z)==v[1:]
    return z,c,beta

def search_plane(d,name,limit,source):
    start=time.monotonic()
    gram=matrix(d['gram']);walls=matrix(d['all_walls']);normals=tuple(mv(gram,r) for r in walls)
    alpha=G.from_json(d['alpha']);v=tuple(map(G.from_json,d['hyperideal_witness']['vector']))
    z,c,beta=plane_parameters(gram,v)
    norms=(G(1),G(2)) if name=='2phi' else (G(1),G(2),G(2,1),G(4,2))
    records=[];found=[]
    for q in norms:
        obstructions=local_obstructions(beta,alpha,q)
        record={'norm':q.to_json(),'local_obstructions':obstructions}
        if obstructions:
            record['status']='NO_ROOT_IN_THIS_PLANE_BY_LOCAL_OBSTRUCTION';records.append(record);continue
        cap=q/beta;k=1
        while (G(k*k)-cap).sign(1)<0:k+=1
        ca,cb,den=scale_bound_coefficients(cap)
        tested=conjugate_admissible=integral_squares=0
        ratio=beta/alpha;constant=q/alpha
        stop=False
        for b in range(-k-1,(limit+k)//2+3):
            ab=abs(b);floor=(isqrt(5*ab*ab)-ab)//2
            if b<0:floor=-floor-1
            for a in range(floor-k-1,floor+k+3):
                if ri_sign(a,b)<=0 or ri_sign(limit-a,-b)<0:continue
                aa=a*a+b*b;bb=2*a*b+b*b
                if ri_sign(ca-den*(aa+bb),cb+den*bb)<0:continue
                conjugate_admissible+=1;t=G(a,b)
                if not (2*t/q).is_integral:continue
                tested+=1;square=ratio*G(aa,bb)-constant
                if not square.is_integral:continue
                integral_squares+=1;x=golden_square_root(square)
                if x is None or not x.is_integral:continue
                if x.sign()<0:x=-x
                if (x-c*t).sign()<0:continue
                r=(x,*(t*y for y in z));q_actual=dot(r,mv(gram,r))
                assert q_actual==q and all(y.is_integral for y in r)
                assert all((2*y/q).is_integral for y in mv(gram,r))
                assert all(dot(r,n).sign()<=0 for n in normals)
                assert dot(r,mv(gram,v)).sign()>0
                if not all(coxeter_pair_value(4*dot(r,n)**2/(q*dot(w,n))) for w,n in zip(walls,normals)):
                    continue
                inside=tuple(map(G.from_json,d['strict_interior']))
                assert dot(r,mv(gram,inside)).sign()<0
                found.append({'root':enc(r),'norm':q.to_json(),'x':x.to_json(),'t':t.to_json(),
                    'escape_pairing':dot(r,mv(gram,v)).to_json(),'wall_pairings':enc(tuple(dot(r,n) for n in normals)),
                    'old_basepoint_height':(alpha*x*x/q).to_json()})
                stop=True;break
            if stop:break
        record.update(status='EXACT_COMPATIBLE_ESCAPE_CUT_FOUND' if stop else 'NO_CUT_WITHIN_STATED_PLANE_BOUND',
            conjugate_admissible_t=conjugate_admissible,reflection_admissible_t=tested,integral_square_candidates=integral_squares,
            completed_bound=not stop)
        records.append(record)
        print(json.dumps({'host':name,**record}),flush=True)
        if found:break
    u=matrix(d['standard_block_swap_involution']);uv=mv(u,v);qv=dot(v,mv(gram,v))
    out={'source':source,'alpha':alpha.to_json(),'saturated_spatial_direction':enc(z),'spatial_norm':beta.to_json(),
        'escape_spatial_multiplier':c.to_json(),'plane_equation':'beta*t^2-alpha*x^2=q; t,x in Z[phi]',
        'positive_t_bound':limit,'norm_searches':records,'cuts':found,'escape_fixed_by_block_swap':uv==v,
        'escape_swap_pairing':dot(v,mv(gram,uv)).to_json(),
        'escape_polar_reflective':all((2*y/qv).is_integral for y in mv(gram,v)),
        'seconds':round(time.monotonic()-start,3),
        'scope':'Exact cuts or bounded saturated rank-two slice failure only. No ambient root-exhaustion or global nonreflectivity claim.'}
    print(json.dumps({'host':name,'cuts':len(found),'seconds':out['seconds']}),flush=True)
    return out

def run(name,limit):
    source=f'results/fresh_multibasepoint_h9_{name}_20261013.json'
    out=search_plane(read(source),name,limit,source)
    output=ROOT/f'results/fresh_multibasepoint_h9_{name}_escape_plane_20261013.json'
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--limit',type=int,default=1000000);a=p.parse_args()
    for name in ('2phi','sqrt5'):run(name,a.limit)
