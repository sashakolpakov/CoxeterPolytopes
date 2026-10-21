#!/usr/bin/env python3
"""Exact sparse-support/local-profile audit and bounded new-section replay.

No floating coordinate or geometry decision is trusted. The numerical file
only supplies rational construction recipes and proposed defect indices.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json
import fresh_h20_cusp_sections_20261013 as P
import fresh_h22_multiquadratic_20261013_review as R
ROOT=P.ROOT;TAG='fresh_h20_compact_section_coverage_20261014'
SCOUT=f'results/{TAG}_scout.json';OUT=ROOT/f'results/{TAG}.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def rat(a):assert not(set(a)-{1});return a.get(1,F(0))
def run():
 scout=json.loads((ROOT/SCOUT).read_text())
 for path,digest in scout['sources'].items():assert sha(path)==digest,path
 raw=json.loads((ROOT/P.RAW).read_text())['coxeter_input']
 Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[P.mv(Q,w)for w in W];G=[[P.dot(w,c)for c in C]for w in W]
 walk=json.loads((ROOT/P.WALK).read_text());corners=walk['corner_orbits']
 total=[sum(w[k]for w in W)for k in range(22)]
 import math
 gcd=math.gcd(*total);h=[F(x,gcd)for x in total];hc=P.mv(Q,h)
 assert P.pair(h,h,Q)==-14 and max(P.dot(h,c)for c in C)<0

 seeds=[];cache={};counts=Counter();out=[]
 for ci,row in enumerate(scout['cases']):
  assert row['status']=='NUMERICAL_ACTUAL_VERTEX_DEFECT'
  si=row['seed'];seed=scout['selected_edges'][si];contact_squares=row['contact_squares']
  key=(si,tuple(contact_squares))
  if key not in cache:
   c=corners[seed['corner']];edge=c['edge_orbits'][seed['edge_index']];E=seed['edge_walls']
   assert E==edge['edge_walls'] and len(E)==20 and len(set(E))==20
   A=[[F(G[i][j])for j in E]for i in E];piv=P.pd(A);I=P.inv(A)
   a=[{}for i in E]
   for j,rad in zip(seed['contact_walls'],contact_squares):
    coef='1'
    assert not a[E.index(j)];a[E.index(j)]=R.scale(R.sqrtq(F(rad)),F(coef))
   assert sum(bool(x)for x in a)==2
   v=R.mv(I,a);w=R.mv([[W[j][k]for j in E]for k in range(22)],v)
   spatial=rat(R.pair(a,v));gamma=1-spatial
   assert 0<spatial<1 and gamma>0
   assert R.pair(w,R.mv(Q,w))==R.C(spatial)
   assert [R.pair(w,R.ratvec(C[i]))for i in E]==a
   H=[[R.add(R.C(A[i][j]),R.scale(R.mul(a[i],a[j]),-1))for j in range(20)]for i in range(20)]
   # G is PD and a^T G^-1 a<1, so H is PD. Every pair is checked separately.
   diag=[rat(H[i][i])for i in range(20)];assert min(diag)>0;labels=Counter()
   for i,j in combinations(range(20),2):
    z=H[i][j];assert R.sign(z)<=0;sq=R.mul(z,z)
    matches=[m for m,c2 in[(2,F(0)),(3,F(1,4)),(4,F(1,2)),(6,F(3,4))]if sq==R.C(diag[i]*diag[j]*c2)]
    assert len(matches)==1;labels[str(matches[0])]+=1
   left=list(map(F,c['primitive_ray']));right=list(map(F,edge['endpoint_primitive_ray']))
   left=P.scale(-1/P.dot(left,hc),left);right=P.scale(-1/P.dot(right,hc),right)
   assert P.pair(left,left,Q)<=0 and P.pair(right,right,Q)<=0
   assert max(P.dot(left,r)for r in C)==max(P.dot(right,r)for r in C)==0

   aa,bb=map(F,contact_squares)
   assert aa*bb==F(1,4)*(2-aa)*(2-bb)
   assert spatial==(aa+bb)/2
   assert sum(labels.values())==190 and labels['3']>=1
   assert all(G[j][k]==0 for j in seed['contact_walls'] for k in E if j!=k)
   seeds.append({'seed':si,'contact_squares':contact_squares,'edge_walls':E,
                 'spatial_norm':spatial,'positive_residual_norm':gamma,
                 'old_edge_LDL_pivots':piv,'projected_pair_labels':dict(labels)})
   cache[key]=(E,I,a,v,w,gamma,left,right)
  E,I,a,v,w,gamma,left,right=cache[key]
  sgn=row['radical_sign'];assert sgn in[-1,1]
  t=F(row['parameter']);assert t in[F(1,5),F(1,2),F(4,5)]
  point=P.add(P.scale(1-t,left),P.scale(t,right))
  qp=P.pair(point,point,Q);assert qp<0
  pp=[P.dot(point,c)for c in C];assert [i for i,p in enumerate(pp)if p==0]==E
  assert max(pp[i]for i in range(210)if i not in E)<0
  z0=P.add(right,P.scale(-P.pair(right,point,Q)/qp,point));qz=P.pair(z0,z0,Q);assert qz>0
  assert all(P.dot(z0,C[i])==0 for i in E)and P.pair(z0,point,Q)==0
  z=R.vm(R.ratvec(z0),R.scale(R.sqrtq(gamma/qz),sgn));normal=R.vadd(w,z)
  nc=R.mv(Q,normal);ep=R.mv(C,normal)
  assert R.pair(normal,nc)==R.C(1)and R.pair(R.ratvec(point),nc)=={}
  assert [ep[i]for i in E]==a
  witness=row['witness'];release=witness['release'];new=witness['new_wall'];k=E.index(release)
  assert new not in E
  wj=[sum(I[k][l]*W[j][d]for l,j in enumerate(E))for d in range(22)]
  dual=R.vadd(R.ratvec(wj),R.vm(z,R.scale(v[k],-1/gamma)))
  dp=R.pair(dual,R.ratvec(C[new]));assert R.sign(dp)<0 and pp[new]<0
  assert R.pair(dual,nc)=={} and R.pair(dual,R.ratvec(P.mv(Q,point)))=={}
  assert [R.pair(dual,R.ratvec(C[j]))for j in E]==[R.C(int(j==release))for j in E]
  qt=R.pair(dual,R.mv(Q,dual));assert R.sign(qt)>0
  # Positive scaling avoids division in the multiquadratic coefficient ring.
  X=R.vadd(R.vm(R.ratvec(point),R.scale(dp,-1)),R.vs(dual,pp[new]))
  qX=R.pair(X,R.mv(Q,X));future=R.pair(X,R.ratvec(hc))
  assert qX==R.add(R.scale(R.mul(dp,dp),qp),R.scale(qt,pp[new]**2))
  assert R.sign(qX)<0 and R.sign(future)<0 and R.pair(X,nc)=={}
  xp=R.mv(C,X);signs=[R.sign(p)for p in xp];assert max(signs)==0
  active=[i for i,s in enumerate(signs)if s==0]
  assert active==sorted(set(E)-{release}|{new})==witness['active']
  assert len(active)==20 and all(signs[i]<0 for i in range(210)if i not in active)
  # The retained19 normals are PD. The new normal is independent because
  # (dual,new)=dp !=0, and all20 are orthogonal to timelike X, hence PD.
  i,j=witness['bad_pair'];assert i in active and j in active and i!=j
  p=R.add(R.C(G[i][j]),R.scale(R.mul(ep[i],ep[j]),-1))
  qi=R.add(R.C(G[i][i]),R.scale(R.mul(ep[i],ep[i]),-1))
  qj=R.add(R.C(G[j][j]),R.scale(R.mul(ep[j],ep[j]),-1));prod=R.mul(qi,qj);square=R.mul(p,p)
  assert R.sign(qi)>0 and R.sign(qj)>0 and R.sign(R.add(prod,R.scale(square,-1)))>0
  kind=witness['kind']
  if kind=='obtuse':assert R.sign(p)>0
  elif kind=='between_1/4_and_1/2':
   assert R.sign(p)<0 and R.sign(R.add(R.scale(square,4),R.scale(prod,-1)))>0
   assert R.sign(R.add(R.scale(square,2),R.scale(prod,-1)))<0
  elif kind=='between_0_and_1/4':
   assert R.sign(p)<0 and R.sign(square)>0 and R.sign(R.add(R.scale(square,4),R.scale(prod,-1)))<0
  else:raise AssertionError(kind)
  counts[kind]+=1
  out.append({'case':ci,'seed':si,'contact_squares':row['contact_squares'],'parameter':str(t),'radical_sign':sgn,
    'normal':normal,'initial_point':point,'initial_point_norm':qp,
    'release':release,'new_wall':new,'positive_scaled_first_contact':X,'contact_norm':qX,
    'contact_center_pairing':future,'active_walls':active,'all210_wall_pairing_signs':signs,
    'defect':kind,'actual_ridge': [i,j],'projected_pairing':p,
    'projected_norms':[qi,qj],'positive_two_by_two_determinant':R.add(prod,R.scale(square,-1))})
 assert len(out)==162 and len(seeds)==27
 sources=[P.RAW,P.WALK,SCOUT,f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',
          'scripts/fresh_h22_multiquadratic_20261013_review.py']
 result={'status':'PASS_EXACT_162_NEW_SAME_SIGN_A1_A1_TO_A2_SECTIONS',
         'sources':{p:sha(p) for p in sources},'seed_profile_count':len(seeds),
         'case_count':len(out),'counts':dict(counts),'seeds':seeds,'cases':out,
         'scope':'Nine source edges, three same-sign A1+A1 to A2 contact ratios, three rational positions and two residual signs. Every recipe has a certified actual finite-vertex defect. Neither full local-profile classification nor all-real-parameter exclusion.'}
 return enc(result)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run()
 if args.verify:assert z==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status']);print(z['counts'])
