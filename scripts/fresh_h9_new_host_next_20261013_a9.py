#!/usr/bin/env python3
"""Changed irreducible A9 host: exact complete bounded positive-root census."""
import argparse,hashlib,json,time
from functools import lru_cache
from itertools import groupby
from pathlib import Path
from fresh_h9_glued_host_20261013_arithmetic import *
from fresh_h9_quartic_mixed_20261013 import balanced_alphabet
from fresh_h9_new_host_20261013 import compatible,geometry

ROOT=Path(__file__).resolve().parents[1]
ALPHA=T-1
A=tuple(tuple(K(2 if i==j else -1 if abs(i-j)==1 else 0)for j in range(9))for i in range(9))
AI=inverse(A)
G=tuple(tuple(-ALPHA if i==j==0 else A[i-1][j-1]if i and j else Z for j in range(10))for i in range(10))
SCH=tuple(inverse(tuple(row[:k]for row in A[:k]))for k in range(1,10))
DYADIC=(O,PI,T*T,2/PI,K(2),2*PI,2*T*T,4/PI,K(4))
POOL=DYADIC+tuple(5*q for q in DYADIC)

@lru_cache(None)
def spatial(E,q):
    alph=[]
    for j in range(9):
        caps=(4*E/(q*q*AI[j][j]),*(8*E/(q*q)for _ in range(3)))
        alph.append(tuple(q*z/2 for z in balanced_alphabet(caps)if(q*z/2).is_integral))
    found=[];nodes=0
    def walk(prefix):
        nonlocal nodes
        nodes+=1;k=len(prefix)
        if k:
            lower=pair(prefix,SCH[k-1],prefix)
            if not le_all(lower,E):return
            if pair(prefix,tuple(row[:k]for row in AI[:k]),prefix)>E:return
        if k==9:
            if lower!=E:return
            y=mv(AI,prefix)
            if all(v.is_integral for v in y):found.append(y)
            return
        for v in alph[k]:walk((*prefix,v))
    walk(())
    return tuple(found),{'coordinate_alphabets':list(map(len,alph)),'nodes':nodes}

def structure():
    # p(T) is irreducible modulo5: no root and no monic quadratic divisor.
    def rem(poly,div):
        p=list(poly)
        while len(p)>=len(div):
            c=p[-1];d=len(p)-len(div)
            for i,z in enumerate(div):p[d+i]=(p[d+i]-c*z)%5
            while p and p[-1]==0:p.pop()
        return p
    poly=[2,0,1,0,1]
    assert all(sum(c*x**i for i,c in enumerate(poly))%5 for x in range(5))
    assert all(rem(poly,[a,b,1])for a in range(5)for b in range(5))
    # exact inverse and determinant are standard A9 formulas, checked here.
    assert AI==tuple(tuple(K(F(min(i,j)*(10-max(i,j)),10))for j in range(1,10))for i in range(1,10))
    assert all((10*z).is_integral for row in AI for z in row)
    assert ALPHA.norm()==-1 and all(ALPHA.sign(e)<0 for e in range(1,4))and ALPHA>0
    assert all(q.sign(e)>0 for q in POOL for e in range(4))
    for i,q in enumerate(DYADIC):
        unit=q/(PI**i);assert unit.is_integral and unit.inverse().is_integral and unit.norm()==1
    return {'determinant':'-10*(T-1)','dual_exponent':10,'safe_primitive_norm_divisor':20,
      'norm_pool':enc(POOL),'five_inert_mod_polynomial_verified':True,
      'zero_chamber':'A9','zero_mirrors':45,
      'zero_completeness_argument':'The positive-definite integral root group is finite and contains irreducible A9. An irreducible finite real reflection group of rank9 has type A9, B9, or D9. The order10! does not divide2^9*9! or2^8*9!, so the full zero-height reflection group is A9.',
      'norm_pool_scope':'Safe complete primitive pool from q dividing2 times the dual exponent, PID and totally-positive-units-are-squares certificate. The dyadic prime has ramification index4;5 is inert. Some displayed classes may contain no primitive roots.',
      'comparison':'Prime5 has odd determinant valuation, unlike the two compared F4²+1 hosts; rank10 similarity preserves determinant squareclass. No full-group or worldwide novelty claim.'}

def run(height,out):
    start=time.monotonic();s=structure();walls=[tuple(-O if i==j else Z for i in range(10))for j in range(1,10)]
    controls=[];discard=0
    for q in POOL:
        for x in balanced_alphabet((height*q/ALPHA,*(-q/ALPHA for _ in range(3)))):
            if x<=0:continue
            if not(2*ALPHA*x/q).is_integral:discard+=1;continue
            controls.append((ALPHA*x*x/q,q,x))
    controls.sort();data={'status':'RUNNING','completed':False,'structure':s,'height_bound':height,
      'gram':menc(G),'norm_pool':enc(POOL),'control_count':len(controls),'time_integrality_rejections':discard,
      'initial_walls':menc(walls),'levels':[]}
    print(json.dumps({'phase':'positive','controls':len(controls),'time_rejections':discard}),flush=True)
    keys={ray_key(r)for r in walls}
    for h,packet in groupby(controls,key=lambda row:row[0]):
        added=[]
        for _,q,x in packet:
            E=q+ALPHA*x*x;ys,stats=spatial(E,q);candidates=[];imprimitive=0
            for y in ys:
                r=(x,*y);assert pair(r,G,r)==q and all((2*z/q).is_integral for z in mv(G,r))
                if all((z/T).is_integral for z in r)or all((z/5).is_integral for z in r):imprimitive+=1;continue
                candidates.append(r)
            kept=[]
            for r in candidates:
                if ray_key(r)in keys or not compatible(r,q,walls,G):continue
                keys.add(ray_key(r));kept.append(r);added.append(r)
            data['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
              'candidates':menc(candidates),'accepted':menc(kept),'imprimitive':imprimitive,**stats})
            print(json.dumps({'height':h.real(),'q':str(q),'candidates':len(candidates),'new':len(kept),'elapsed':round(time.monotonic()-start,2),**stats}),flush=True)
        for i,r in enumerate(added):assert compatible(r,pair(r,G,r),added[:i],G)
        walls+=added;data.update(all_walls=menc(walls),wall_count=len(walls));out.write_text(json.dumps(data,indent=2)+'\n')
    direction=mv(AI,(O,)*9);eps=O
    while True:
        inside=(O,*(eps*y for y in direction))
        if pair(inside,G,inside)<0 and all(pair(inside,G,r)<0 for r in walls):break
        eps/=2
    assert all(compatible(r,pair(r,G,r),walls[:i],G)for i,r in enumerate(walls))
    data.update(status='COMPLETE_BOUNDED_A9_HOST',completed=True,all_walls=menc(walls),wall_count=len(walls),
      strict_interior=enc(inside),geometry=geometry(walls,G,inside,9),seconds=round(time.monotonic()-start,3),
      scope='Complete dominant primitive-root census in the irreducible A9 host through the stated height. Zero-height completeness uses finite reflection-group classification and group orders. Saved cone geometry is separate from the unresolved full host.')
    out.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'status':data['status'],'walls':len(walls),'geometry':data['geometry']['status'],'seconds':data['seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=30);p.add_argument('--output',type=Path,required=True);args=p.parse_args();run(args.height,args.output)
