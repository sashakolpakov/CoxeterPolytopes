#!/usr/bin/env python3
"""Solve complete local Coxeter-angle equations on two-direction cusp planes."""
from itertools import combinations
from collections import Counter
import json
import numpy as np
import sympy as sp
import fresh_h20_joint_cusp_20261013 as J
import fresh_h20_global_section_parameters_20261013 as algebra
P=J.old.P;t=J.old.t;ROOT=J.ROOT
OUT=ROOT/'data/enumerations/fresh_h20_joint_cusp_directions_20261013.json'

def local(Q,W,G,u,h,comp,e0,e1,drop):
 ids=[i for i in comp if i!=drop];q=P(e0@Q@e0+2*t*(e0@Q@e1)+t*t*(e1@Q@e1))
 v={i:P(W[i]@Q@e0+t*(W[i]@Q@e1)) for i in ids}
 gr=lambda i,j:q*G[i,j]-v[i]*v[j]
 gcd=P(0);constraints=[]
 for i,j in combinations(ids,2):
  if v[i].is_zero and v[j].is_zero:continue
  N=gr(i,j)**2;den=gr(i,i)*gr(j,j)
  poly=N*(4*N-den)*(2*N-den)*(4*N-3*den)*(N-den)
  if poly.is_zero:continue
  gcd=sp.gcd(gcd,poly).monic();constraints.append([i+1,j+1])
  if gcd.degree()==0:break
 return q,v,gr,gcd,constraints

def kernel_field(matrix,F):
 n=len(matrix);M=[[F.rem(x) for x in r] for r in matrix];piv=[];rank=0
 for c in range(n):
  p=next((r for r in range(rank,n) if not M[r][c].is_zero),None)
  if p is None:continue
  M[rank],M[p]=M[p],M[rank];inv=sp.invert(M[rank][c],F.p);M[rank]=[F.rem(x*inv) for x in M[rank]]
  for r in range(n):
   if r==rank:continue
   d=M[r][c];M[r]=[F.rem(a-d*b) for a,b in zip(M[r],M[rank])]
  piv.append(c);rank+=1
 free=[i for i in range(n) if i not in piv]
 if len(free)!=1:return None
 v=[P(0) for _ in range(n)];v[free[0]]=P(1)
 for r,c in enumerate(piv):v[c]=-M[r][free[0]]
 signs=[F.sign(x) for x in v]
 if max(signs)<0:v=[-x for x in v];signs=[-s for s in signs]
 return v if min(signs)>0 else None

def run():
 Q,W,G,h=J.old.base.load();source=json.loads(J.old.CUSP_SOURCE.read_text());counts=Counter();results=[];survivors=[]
 for ci,raw in enumerate(source['cusps']):
  u=np.array(raw['vector'],dtype=object)
  if max(W@Q@u)>0:continue
  ids=[i for i,p in enumerate(W@Q@u) if p==0]
  for comp in J.old.components(G,ids):
   dirs=J.two_node_directions(Q,W,G,h,u,comp);lookup={frozenset(d['bridge']):d for d in dirs}
   for b in comp:
    neighbors=[i for i in comp if i!=b and G[i,b]]
    for a,c in combinations(neighbors,2):
     A=lookup[frozenset((a,b))];B=lookup[frozenset((b,c))]
     e0,e1=A['e'],B['e']
     for drop in [None]+comp:
      q,v,gr,gcd,constraints=local(Q,W,G,u,h,comp,e0,e1,drop);counts['local_plane_facet_patterns']+=1
      record={'cusp_index':ci,'path':[a+1,b+1,c+1],'component':[i+1 for i in comp],
              'dropped':None if drop is None else drop+1,'e0':list(map(str,e0)),'e1':list(map(str,e1)),
              'gcd':algebra.coeffs(gcd),'constraints':constraints,'roots':[]}
      if gcd.is_zero:
       active=[i for i in comp if i!=drop]
       adj=np.asarray([[int(not gr(i,j).is_zero) for j in active] for i in active])
       fixed=next(([active[k]+1 for k in ks] for ks in J.old.components(adj) if all(v[active[k]].is_zero for k in ks)),None)
       if fixed is not None:
        record['continuous_status']='FIXED_POSITIVE_DEFINITE_COMPONENT';record['fixed_component']=fixed;counts['fixed_finite_block']+=1
       else:
        record['continuous_status']='CONTINUOUS_LOCAL_ANGLE_FAMILY';counts['continuous_local_angle_family']+=1
        print('CONTINUOUS',ci,[a+1,b+1,c+1],None if drop is None else drop+1,flush=True)
       results.append(record);continue
      active=[i for i in comp if i!=drop]
      for poly,power in sp.factor_list(gcd)[1]:
       for (lo,hi),mult in poly.intervals(eps=sp.Rational(1,10**14)):
        if lo==hi:lo-=sp.Rational(1,10**18);hi+=sp.Rational(1,10**18)
        F=algebra.Field(poly,lo,hi);r={'polynomial':algebra.coeffs(poly.monic()),'interval':[str(lo),str(hi)]};counts['necessary_real_roots']+=1
        H=[[F.rem(gr(i,j)) for j in active] for i in active]
        bad=None
        for k,i in enumerate(active):
         if F.sign(H[k][k])<=0:bad={'kind':'nonpositive_norm','wall':i+1};break
        if bad is None:
         for k,l in combinations(range(len(active)),2):
          N=F.rem(H[k][l]**2);den=F.rem(H[k][k]*H[l][l])
          if F.sign(H[k][l])>0 or all(not F.rem(N-z*den).is_zero for z in J.ALLOWED):
           bad={'kind':'non_Coxeter_pair','walls':[active[k]+1,active[l]+1]};break
        cs=[]
        if bad is None:
         ints=np.asarray([[int(not p.is_zero) for p in row] for row in H]);
         for ks in J.old.components(ints):
          marks=kernel_field([[H[k][l] for l in ks] for k in ks],F)
          if marks is None:bad={'kind':'no_positive_affine_kernel','walls':[active[k]+1 for k in ks]};break
          cs.append({'walls':[active[k]+1 for k in ks],'marks':[algebra.coeffs(p) for p in marks],'rank':len(ks)-1})
        if bad is None:
         assert sum(c['rank'] for c in cs)==len(comp)-2
         r['components']=cs;r['status']='LOCAL_COXETER_DIRECTION';counts['local_Coxeter_directions']+=1
         survivors.append({**record,'root':r,'u':list(map(int,u))})
         print('SURVIVOR',ci,[a+1,b+1,c+1],r['polynomial'],None if drop is None else drop+1,flush=True)
        else:r['status']=bad['kind'];r['rejection']=bad;counts[bad['kind']]+=1
        record['roots'].append(r)
      results.append(record)
 print(dict(counts),flush=True)
 OUT.write_text(json.dumps({'schema':'two-direction-cusp-plane-local-v1','counts':dict(counts),'cases':results,'survivors':survivors},indent=2)+'\n')
if __name__=='__main__':run()
