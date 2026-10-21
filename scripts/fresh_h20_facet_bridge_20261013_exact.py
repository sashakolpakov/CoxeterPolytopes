#!/usr/bin/env python3
"""Two exact obstructions to specified H20 facet-bridge constructions.

(1) No compatible positive-norm normal is perpendicular to facet175 and all
18 of its finite neighbors. The positive two-plane and all 14 possible
unit-normal branches are certified with outward-rounded integer intervals.
(2) The saved strict176 escape survives every polar reflection of the175
predecessor normals. This does not require screening reflected pair angles.
Neither statement excludes general intersecting additions or replacements.
"""
import argparse,json
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
from fresh_h20_coupled_completion_20261013_second_exact import setup,OUT as SOURCE,e,point,pairing,norm
from fresh_h20_coupled_completion_20261013_orbit_exact import OUT as POLAR
from fresh_h20_coupled_completion_20261013 import ROOT,read,sha
from fresh_h20_edge_family_completion_20261013_cusp_repairs import base_root
from fresh_dimension_bridge_20261013_verify import E,Z,O,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv
from fresh_h9_split_leaf_3034_safe_graph_20261013 import atan_bounds,cos_interval,TRIG_SCALE
OUT='data/enumerations/fresh_h20_facet_bridge_20261013_exact.json'
S=10**60

def ceildiv(a,b):return-((-a)//b)
def rat(q):
 q=F(q)*S;return(q.numerator//q.denominator,ceildiv(q.numerator,q.denominator))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 p=[x*y for x in a for y in b];return(min(p)//S,ceildiv(max(p),S))
def div(a,b):
 assert b[0]>0 or b[1]<0
 p=[F(x*S,y)for x in a for y in b];lo,hi=min(p),max(p)
 return(lo.numerator//lo.denominator,ceildiv(hi.numerator,hi.denominator))
def sqrt(a):
 assert a[0]>=0
 lo=isqrt(a[0]*S);hi=isqrt(a[1]*S)+1
 assert lo*lo<=a[0]*S and hi*hi>=a[1]*S
 return lo,hi
@lru_cache(None)
def interval(f,x):
 if f is B:return add(rat(x.a),mul(rat(x.b),sqrt(rat(6))))
 return add(interval(f.base,x[0]),mul(interval(f.base,x[1]),sqrt(interval(f.base,f.d))))
def isum(xs):
 s=rat(0)
 for x in xs:s=add(s,x)
 return s
def idot(a,b):return isum(mul(x,y)for x,y in zip(a,b,strict=True))
def enc(x):return list(map(str,x))
def nullspace(rows):
 a=[list(r)for r in rows];piv=[];n=len(a[0]);k=0
 for j in range(n):
  p=next((i for i in range(k,len(a))if a[i][j].sign()),None)
  if p is None:continue
  a[k],a[p]=a[p],a[k];q=a[k][j];a[k]=[x/q for x in a[k]]
  for i in range(len(a)):
   if i!=k:
    q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[k],strict=True)]
  piv.append(j);k+=1
  if k==len(a):break
 free=[j for j in range(n)if j not in piv];basis=[]
 for j in free:
  x=[Z]*n;x[j]=O
  for i,p in enumerate(piv):x[p]=-a[i][j]
  assert all(not dot(r,x).sign()for r in rows);basis.append(tuple(x))
 return basis,piv

def angles():
 al,ah=atan_bounds(5,60);bl,bh=atan_bounds(239,18);pl,ph=16*al-4*bh,16*ah-4*bl
 t2=2*F(1,5)/(1-F(1,25));t4=2*t2/(1-t2*t2);assert(t4-F(1,239))/(1+t4*F(1,239))==1 and 3<pl<ph<4
 ans={2:rat(0)}
 for m in range(3,10):
  lo=pl*TRIG_SCALE/m;hi=ph*TRIG_SCALE/m
  cl,ch=cos_interval(lo.numerator//lo.denominator,ceildiv(hi.numerator,hi.denominator))
  ans[m]=(cl*S//TRIG_SCALE,ceildiv(ch*S,TRIG_SCALE))
 return ans,[str(pl),str(ph)]

def run():
 Q,h,hc,roots,f=setup();d=read(SOURCE);g=Extension(f,f.parse(d['field']['radicand']));r=tuple(map(g.parse,d['new_cap']['normal']));roots.append((g,r));assert len(roots)==176
 ids=[p['old']for p in d['new_cap']['old_pair_checks']if p['kind']!='dotted'];assert ids==[0,2,3,4,5,6,8,9,11,12,13,14,15,16,17,20,22,40]
 rows=[mv(Q,base_root(*roots[i]))for i in ids];K,piv=nullspace(rows);assert len(piv)==18 and len(K)==3
 rc=fmv(g,tuple(tuple(e(g,x)for x in row)for row in Q),r);a=[fdot(g,rc,tuple(e(g,x)for x in v))for v in K];pivot=next(i for i,x in enumerate(a)if g.sign(x));rest=[i for i in range(3)if i!=pivot]
 UV=[tuple(g.sub(g.mul(a[pivot],e(g,x)),g.mul(a[j],e(g,y)))for x,y in zip(K[j],K[pivot],strict=True))for j in rest]
 assert all(g.sign(fdot(g,v,rc))==0 for v in UV)
 assert all(g.sign(fdot(g,v,tuple(e(g,x)for x in row)))==0 for v in UV for row in rows)
 FQ=tuple(tuple(e(g,x)for x in row)for row in Q);cov=[fmv(g,FQ,v)for v in UV];H=[[fdot(g,v,c)for c in cov]for v in UV];assert g.sign(H[0][0])>0
 det=g.sub(g.mul(H[0][0],H[1][1]),g.mul(H[0][1],H[1][0]));assert g.sign(det)>0
 # Orthonormal coordinates: e1=U/sqrt(H00), e2=(V-H01/H00*U)/sqrt(det/H00).
 hi=[[interval(g,x)for x in row]for row in H];sq0=sqrt(hi[0][0]);mix=div(hi[0][1],hi[0][0]);sq1=sqrt(div(interval(g,det),hi[0][0]));proj={}
 for i in [10,18,30,50]:
  rr=base_root(*roots[i]);p=[interval(g,fdot(g,v,tuple(e(g,x)for x in mv(Q,rr))))for v in UV];q=dot(rr,mv(Q,rr));assert q.sign()>0
  raw=[div(p[0],sq0),div(sub(p[1],mul(mix,p[0])),sq1)];qn=sqrt(interval(B,q));proj[i]=[div(x,qn)for x in raw]
 A=proj[18];aa=idot(A,A);cs,pi=angles();assert aa[1]<mul(cs[9],cs[9])[0] and aa[0]>mul(cs[8],cs[8])[1]
 results=[]
 for m in range(2,9):
  cm=cs[m];t=div(neg(cm),aa);scale=div(sqrt(sub(rat(1),div(mul(cm,cm),aa))),sqrt(aa));parallel=[mul(t,x)for x in A];perp=[mul(scale,neg(A[1])),mul(scale,A[0])]
  for branch in [-1,1]:
   u=[add(x,y if branch==1 else neg(y))for x,y in zip(parallel,perp,strict=True)]
   # Only one strictly forbidden old pair is needed for each possible branch.
   witnesses=[]
   for i,b in proj.items():
    z=neg(idot(u,b))
    if z[1]<0:witnesses.append({'wall':i,'cosine':enc(z),'reason':'positive_pairing'})
    else:
     for k in range(2,9):
      if z[0]>cs[k][1] and z[1]<cs[k+1][0]:
       witnesses.append({'wall':i,'cosine':enc(z),'reason':'strict_between_adjacent_labels','labels':[k,k+1]});break
   assert witnesses,(m,branch)
   results.append({'anchor_label':m,'branch':branch,'unit_coordinates':list(map(enc,u)),'rejection':witnesses[0]})
 # Existing exact escape excludes every reflected predecessor inequality.
 P=tuple(map(f.parse,read(POLAR)['final_hyperideal_point']['vector']));pq=norm(f,P,Q);assert f.sign(pq)>0
 pc=fmv(f,tuple(tuple(e(f,x)for x in row)for row in Q),P);pp=[]
 for ff,v in roots[:175]:
  H0,z=pairing(f,pc,ff,v);assert H0.sign(z)<=0;pp.append({'pairing':H0.encode(z),'sign':H0.sign(z)})
 X=vec(d['strict_escape']['vector']);escape=point(B,X,Q,hc,roots,True);px=fdot(f,P,tuple(e(f,x)for x in mv(Q,X)));assert f.sign(px)<0
 return {'status':'EXACT_TWO_SCOPED_H20_FACET_BRIDGE_OBSTRUCTIONS','source':SOURCE,'source_sha256':sha(SOURCE),'polar_source':POLAR,'polar_sha256':sha(POLAR),'interval_scale':str(S),'perpendicular_bridge':{'finite_neighbors':ids,'constraint_rank':19,'base_constraint_pivots':piv,'base_nullspace':[[x.enc()for x in v]for v in K],'eliminated_coordinate':pivot,'plane_gram':[[g.encode(x)for x in row]for row in H],'plane_gram_determinant':g.encode(det),'normalized_wall18_projection_square':enc(aa),'cosine_intervals':{str(k):enc(v)for k,v in cs.items()},'pi_interval':pi,'branches':results,'scope':'Every positive-norm normal perpendicular to facet175 and all18 finite-neighbor normals fails a necessary Coxeter pair relation. All real directions and all integer labels m>=2 are covered; no center assumption is needed.'},'polar_reflection_bank':{'polar_norm':f.encode(pq),'polar_old_pairings':pp,'strict176_escape':escape,'polar_escape_pairing':f.encode(px),'scope':'The strict escape survives all175 reflected predecessor inequalities r-2(P,r)P/P^2, simultaneously with all176 current inequalities. Thus no selection from this entire reflection bank gives finite volume, even if compatibility screening or deletion of current walls is allowed.'}}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print('14 branches and all175 reflected inequalities certified')
