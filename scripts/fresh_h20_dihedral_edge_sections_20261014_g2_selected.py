#!/usr/bin/env python3
"""Complete exact G2-face closure census on three specified source edges.

The floating-point packet proposes only defective-ridge witness indices.
Every closure root, branch, transverse interval, normal and actual ridge is
independently rebuilt here in rational/quadratic arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse, hashlib, json, math
from sympy import factorint
import fresh_h20_cusp_sections_20261013 as P
import fresh_h22_multiquadratic_20261013_review as R
ROOT=P.ROOT;TAG='fresh_h20_dihedral_edge_sections_20261014_g2'
RECIPES=f'results/{TAG}_selected_recipes.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def sqrtq(x):
 x=F(x)
 if not x:return{}
 assert x>0;n=x.numerator*x.denominator;fac=factorint(n);assert math.prod(int(p)**int(e)for p,e in fac.items())==n
 a=math.prod(int(p)**(int(e)//2)for p,e in fac.items());b=math.prod(int(p)for p,e in fac.items()if e%2)
 out={b:F(a,x.denominator)};assert R.mul(out,out)==R.C(x);return out
def numeric(a):return sum(float(c)*math.sqrt(d)for d,c in a.items())
def key(v):return tuple(tuple(sorted(x.items()))for x in v)
def run():
 recipes=json.loads((ROOT/RECIPES).read_text())
 assert sha(f'results/{TAG}_scout.json.gz')==recipes['archive_sha256']
 raw=json.loads((ROOT/P.RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']];C=[P.mv(Q,w)for w in W]
 G=[[P.dot(w,c)for c in C]for w in W];walk=json.loads((ROOT/P.WALK).read_text())['corner_orbits']
 h=[sum(w[j]for w in W)for j in range(22)];gg=math.gcd(*h);h=[F(x,gg)for x in h];hc=P.mv(Q,h)
 qe=F(3,4);gamma=F(1,12);families=[];out=[];counts=Counter()
 for fid in recipes['selected_family_ids']:
  f=recipes['families'][str(fid)];corner=walk[f['corner']];edge=corner['edge_orbits'][f['edge_index']];E=f['edge_walls'];a=f['anchor'];cc=f['component'];ai=E.index(a)
  assert E==edge['edge_walls']and len(cc)==2 and a in cc and G[cc[0]][cc[1]]==-1 and all(G[i][i]==2 for i in cc)
  assert all(G[i][j]==0 for i in cc for j in E if j not in cc)
  A=[[F(G[i][j])for j in E]for i in E];piv=P.pd(A);I=P.inv(A);assert I[ai][ai]==F(2,3)
  ww=[[sum(I[k][j]*W[l][d]for j,l in enumerate(E))for d in range(22)]for k in range(20)];w=ww[ai]
  left=list(map(F,corner['primitive_ray']));right=list(map(F,edge['endpoint_primitive_ray']))
  left=P.scale(-1/P.dot(left,hc),left);right=P.scale(-1/P.dot(right,hc),right)
  B=[[P.pair(x,y,Q)for y in[left,right]]for x in[left,right]];assert B[0][0]*B[1][1]-B[0][1]**2<0
  eligible=[l for l in range(210)if l not in E and all(G[l][i]==0 for i in cc if i!=a)]
  assert eligible==f['eligible_closure_walls'];records=[];unique={};tangencies=0
  for l in eligible:
   p=[P.dot(x,C[l])for x in[left,right]];pp=P.dot(p,p);assert pp>0
   k=gamma*G[l][a];z0=P.scale(k/pp,p);y=[p[1],-p[0]]
   aa=P.dot(y,P.mv(B,y));bb=2*P.dot(z0,P.mv(B,y));dd=P.dot(z0,P.mv(B,z0))-gamma;disc=bb*bb-4*aa*dd
   rec={'wall':l,'quadratic':[aa,bb,dd],'discriminant':disc};records.append(rec)
   if aa:
    if disc<0:rec['status']='NO_REAL_ROOT';continue
    roots=[R.add(R.C(-bb/(2*aa)),R.scale(sqrtq(disc),F(s,2)/aa))for s in([-1,1]if disc else[0])]
   elif bb:roots=[R.C(-dd/bb)]
   else:assert dd;rec['status']='NO_ROOT';continue
   rec['root_statuses']=[]
   for root in roots:
    zc=[R.add(R.C(z0[i]),R.scale(root,y[i]))for i in range(2)];vals=R.mv(B,zc)
    assert R.pair(zc,vals)==R.C(gamma)and R.pair(zc,R.ratvec(p))==R.C(k)
    ss=[R.sign(x)for x in vals]
    if ss[0]*ss[1]>=0:
     rec['root_statuses'].append('OUTSIDE_TRANSVERSE_INTERIOR');tangencies+=int(0 in ss);continue
    rec['root_statuses'].append('TRANSVERSE_INTERIOR');unique.setdefault(key(zc),(l,zc,vals))
  proposed=[r for r in recipes['candidate_witnesses']if r['family']==fid]
  assert len(unique)==len(proposed)
  for l,zc,vals in unique.values():
   nearby=sorted((max(abs(numeric(zc[i])-r['edge_plane_coefficients'][i])for i in range(2)),j)for j,r in enumerate(proposed))
   assert nearby[0][0]<1e-5 and (len(nearby)==1 or nearby[1][0]>1e-4)
   r=proposed[nearby[0][1]];witness=r['witness'];assert r['status']=='NUMERICAL_ACTUAL_RIDGE_DEFECT'
   z=R.vadd(R.vm(R.ratvec(left),zc[0]),R.vm(R.ratvec(right),zc[1]));normal=R.vadd(R.ratvec(w),z);ep=R.mv(C,normal);nc=R.mv(Q,normal)
   assert R.pair(normal,nc)==R.C(qe)and [ep[i]for i in E]==[R.C(int(i==a))for i in E]
   point=R.vadd(R.vm(R.ratvec(left),vals[1]),R.vm(R.ratvec(right),R.scale(vals[0],-1)))
   if R.sign(R.pair(point,R.ratvec(hc)))>0:point=R.vs(point,-1)
   qp=R.pair(point,R.mv(Q,point));ps=R.mv(C,point);assert R.sign(qp)<0 and R.pair(point,nc)=={}
   assert [i for i,x in enumerate(ps)if not x]==E and all(R.sign(ps[i])<0 for i in range(210)if i not in E)
   release=witness['release'];new=witness['new_wall'];ri=E.index(release)
   dual=R.vadd(R.ratvec(ww[ri]),R.vm(z,R.C(-I[ri][ai]/gamma)))
   assert R.pair(dual,nc)=={}and R.pair(dual,R.mv(Q,point))=={}
   assert [R.pair(dual,R.ratvec(C[i]))for i in E]==[R.C(int(i==release))for i in E]
   dp=R.pair(dual,R.ratvec(C[new]));assert R.sign(dp)<0 and R.sign(ps[new])<0
   X=R.vadd(R.vm(point,R.scale(dp,-1)),R.vm(dual,ps[new]));qX=R.pair(X,R.mv(Q,X));xp=R.mv(C,X);signs=[R.sign(x)for x in xp]
   assert R.sign(qX)<0 and R.sign(R.pair(X,R.ratvec(hc)))<0 and max(signs)==0 and R.pair(X,nc)=={}
   active=[i for i,x in enumerate(signs)if not x];assert active==sorted(set(E)-{release}|{new})==witness['active']
   i,j=witness['pair'];assert i in active and j in active
   p=R.add(R.C(G[i][j]),R.scale(R.mul(ep[i],ep[j]),-1/qe))
   ni=R.add(R.C(G[i][i]),R.scale(R.mul(ep[i],ep[i]),-1/qe));nj=R.add(R.C(G[j][j]),R.scale(R.mul(ep[j],ep[j]),-1/qe));prod=R.mul(ni,nj);sq=R.mul(p,p)
   determinant=R.add(prod,R.scale(sq,-1));assert R.sign(ni)>0 and R.sign(nj)>0 and R.sign(determinant)>0
   kind=witness['kind']
   if kind=='obtuse':assert R.sign(p)>0
   elif kind=='between_0_and_1/4':assert R.sign(p)<0 and R.sign(R.add(prod,R.scale(sq,-4)))>0
   else:raise AssertionError(kind)
   counts[kind]+=1;out.append({'family':fid,'closure_wall':l,'edge_plane_coefficients':zc,
    'normal':normal,'initial_point':point,'initial_point_norm':qp,'release':release,'new_wall':new,
    'actual_contact':X,'contact_norm':qX,'all210_pairing_signs':signs,'active_walls':active,
    'actual_ridge':[i,j],'projected_ridge_pairing':p,'ridge_norms':[ni,nj],
    'positive_ridge_determinant':determinant,'defect':kind})
  families.append({'family':fid,'corner':f['corner'],'edge_index':f['edge_index'],'edge_walls':E,
   'component':cc,'anchor':a,'old_edge_LDL':piv,'eligible_old_wall_quadratics':records,
   'distinct_transverse_normals':len(unique),'boundary_roots':tangencies,
   'source_endpoint_types':[corner['type'],walk[edge['endpoint_orbit']]['type']]})
 assert len(out)==131 and counts=={'obtuse':105,'between_0_and_1/4':26}
 sources=[P.RAW,P.WALK,RECIPES,f'results/{TAG}_scout.json.gz',f'scripts/{TAG}_selected.py',
  'scripts/fresh_h20_cusp_sections_20261013.py','scripts/fresh_h22_multiquadratic_20261013_review.py']
 return enc({'status':'PASS_COMPLETE_THREE_CONTINUOUS_G2_EDGE_FAMILIES_131_ACTUAL_RIDGE_DEFECTS',
  'sources':{p:sha(p)for p in sources},'counts':dict(counts),'families':families,'cases':out,
  'scope':'Three specified anchored source-edge families, covering ideal-finite, ideal-ideal and finite-finite endpoint types. Exact exhaustive necessary G2-face closure roots and all131 actual-ridge obstructions. Does not certify every other G2 family.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run();out=ROOT/f'results/{TAG}_selected.json'
 if args.verify:assert z==json.loads(out.read_text())
 else:out.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status']);print(z['counts'])
