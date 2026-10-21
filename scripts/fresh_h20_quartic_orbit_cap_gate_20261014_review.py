#!/usr/bin/env python3
"""Independent SymPy algebraic-field replay of the fixed quartic cap gate.

No producer or project arithmetic imports. Re-enumerates the necessary old
PSD triples, solves their endpoint equations, including singular rational
identities, and reconstructs the three positive affine seed components.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import sympy as S
from sympy.polys.domains import QQ

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_quartic_orbit_cap_gate_20261014'
SOURCE=f'results/{TAG}.json'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
t,T=S.symbols('t T')

def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def ff(a): return F(int(a.numerator),int(a.denominator))
def dot(a,b): return sum(x*y for x,y in zip(a,b,strict=True))

class Field:
 def __init__(self,p,iv):
  self.p=S.Poly.from_list(list(map(S.Rational,p)),t)
  assert self.p.degree()==4 and self.p.is_irreducible
  self.lo,self.hi=map(F,iv)
  assert self.p.count_roots(self.lo,self.hi)==1
  self.K=QQ.algebraic_field((self.p.as_expr(),t))
  self.ring,self.var=S.polys.rings.ring('T',self.K)
  self.zero,self.one=self.K.zero,self.K.one
  self.cache={self.zero:0}
 def scalar(self,a):
  a=F(a);return self.K.new([QQ(a.numerator,a.denominator)])
 def read(self,a):return self.K.new([QQ(F(x).numerator,F(x).denominator)for x in a])
 def asc(self,a):return self.read(list(reversed(a)))
 def coeffs(self,a):return [F(0)]*(4-len(a.to_list()))+list(map(ff,a.to_list()))
 def encode(self,a):return [str(x)for x in self.coeffs(a)]
 def sign(self,a):
  if a in self.cache:return self.cache[a]
  lo,hi=self.lo,self.hi
  for d in (20,40,80,160,320):
   l=h=F(0)
   for c in self.coeffs(a):
    v=(l*lo,l*hi,h*lo,h*hi);l,h=min(v)+c,max(v)+c
   if l>0:self.cache[a]=1;return 1
   if h<0:self.cache[a]=-1;return -1
   lo,hi=map(F,self.p.refine_root(lo,hi,eps=S.Rational(1,10**d)))
  raise AssertionError('unresolved exact sign')
 def poly(self,cs):return self.ring.from_dict({(i,):c for i,c in enumerate(cs)if c})

def run():
 saved=json.loads((ROOT/SOURCE).read_text())
 local=json.loads((ROOT/LOCAL).read_text())
 walk=json.loads((ROOT/WALK).read_text())
 for doc in (saved,local,walk):
  for p,h in doc['sources'].items():assert sha(p)==h,(p,'stale source')
 row=local['directions'][0];spec=local['systems'][row['system']]
 assert spec['retained_walls']==[23,50,200] and spec['omitted_wall']==7 and spec['cusp']==4
 K=Field(row['minimal_polynomial'],row['root_interval']);c=K.scalar;Z=K.zero
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=[[-F(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(F,r))for r in raw['wall_coordinates_as_columns_transposed']]
 pair=lambda a,b:dot(a,[dot(r,b)for r in Q])
 R=spec['retained_walls'];B=[[pair(W[i],W[j])for j in R]for i in R]
 assert B==[[F(x)for x in r]for r in saved['old_seed_Gram']]
 q=K.read(row['normal_square']);L=K.read(row['parameter_lower_magnitude'])
 a=list(map(K.read,row['retained_pairings']));k=list(map(K.read,row['positive_kernel']))
 delta=list(map(F,row['null_shift']))
 e0=[sum((k[i]*c(W[R[i]][j])for i in range(3)),Z)for j in range(22)]
 assert e0==list(map(K.read,row['normal_constant']))
 kp=lambda x,y:sum((x[i]*c(Q[i][j])*y[j]for i in range(22)for j in range(22)),Z)
 dc=list(map(c,delta));ec=[x-L*d/c(2)for x,d in zip(e0,dc)]
 assert kp(e0,e0)==q and kp(e0,dc)==Z and kp(dc,dc)==Z
 assert K.sign(q)>0 and K.sign(L)>0 and all(K.sign(x)>0 for x in k)
 assert a==[kp(e0,list(map(c,W[i])))for i in R]
 center=[sum(w[j]for w in W)for j in range(22)]
 assert pair(center,center)<0 and pair(delta,center)<0 and all(pair(center,w)<0 for w in W)

 # Exact positive affine kernels. Untouched ranks 11 and 6, target rank 2.
 corner=walk['corner_orbits'][4];seed=[i for i in corner['active_walls']if i!=7]
 assert len(seed)==22 and len(corner['affine_components'])==3
 projs={i:[c(w)-kp(list(map(c,W[i])),ec)*ee/q for w,ee in zip(W[i],ec)]for i in seed}
 components=[]
 for f in corner['affine_components']:
  ids=[i for i in f['walls']if i in seed]
  weights=k if ids==R else list(map(c,map(F,f['positive_kernel'])))
  v=[sum((ww*projs[i][j]for i,ww in zip(ids,weights)),Z)for j in range(22)]
  pivot=next(i for i,x in enumerate(delta)if x)
  mu=v[pivot]/c(delta[pivot]);assert K.sign(mu)>0 and v==[mu*d for d in dc]
  H=[[kp(projs[i],projs[j])for j in ids]for i in ids]
  # Delete any one wall: exact positive pivots prove rank n-1 and PSD.
  for omit in range(len(ids)):
   M=[[H[i][j]for j in range(len(ids))if j!=omit]for i in range(len(ids))if i!=omit]
   while M:
    d=M[0][0];assert K.sign(d)>0
    M=[[M[i][j]-M[i][0]*M[0][j]/d for j in range(1,len(M))]for i in range(1,len(M))]
  assert all(K.sign(x)>0 for x in weights)
  components.append({'walls':ids,'quotient_rank':len(ids)-1,'positive_null_multiple':K.encode(mu)})
 for f,g in combinations(components,2):assert all(kp(projs[i],projs[j])==Z for i in f['walls']for j in g['walls'])
 assert sum(f['quotient_rank']for f in components)==19
 edge_count=math.prod(len(f['walls'])for f in components);assert edge_count==252
 nn=[q*c(B[i][i])-a[i]**2 for i in range(3)]
 for i,j in combinations(range(3),2):
  p=q*c(B[i][j])-a[i]*a[j]
  assert K.sign(p)<0 and 4*p*p==nn[i]*nn[j]

 def label(p,den):
  if p==Z:return 2
  if K.sign(p)>0:return None
  n=p*p
  for m,v in ((3,F(1,4)),(4,F(1,2)),(6,F(3,4))):
   if n==c(v)*den:return m
  # Includes label five without assuming sqrt(5) absent from the field.
  if 16*n*n-12*n*den+den*den==Z and K.sign(2*n-den)>0:return 5
  return None

 counts=Counter();records=[];quad=[];old_cases=[]
 def identity(v,i,j,gi,gj):
  z=next(x for x in range(3)if x not in (i,j))
  A=S.Matrix.hstack(S.Matrix(K.coeffs(k[z])),S.Matrix(K.coeffs(-L/c(2))))
  b=S.Matrix(K.coeffs(v-k[i]*c(gi)-k[j]*c(gj)))
  rank=A.rank()
  assert rank in (1,2)
  if A.row_join(b).rank()!=rank:return None
  counts['T_identity_consistent']+=1
  if rank==2:
   sol=A.gauss_jordan_solve(b)[0];gz,h=map(F,sol)
   ok=gz.denominator==h.denominator==1 and h<0
   if ok:counts['integral_negative_h']+=1
   return {'third':str(gz),'h':str(h),'admissible':ok}
  p=next(i for i,x in enumerate(A[:,0])if x)
  rho=F(A[p,1]/A[p,0]);tau=F(b[p]/A[p,0]);assert rho
  assert A[:,1]==rho*A[:,0] and b==tau*A[:,0]
  den=math.lcm(rho.denominator,tau.denominator)
  ok=(tau*den).numerator%math.gcd(den,(rho*den).numerator)==0
  if ok:counts['integral_negative_h']+=1
  return {'third_plus_rho_h_equals_tau':[str(rho),str(tau)],'admissible':ok,'h':'negative integer in the indicated congruence class'}
 def evaluate(v,i,j,gi,gj,norm,mode):
  nl=q*c(norm)-v*v
  if K.sign(nl)<=0:return
  ps=[q*c(g)-a[x]*v for x,g in ((i,gi),(j,gj))]
  labels=[label(p,nn[x]*nl)for p,x in zip(ps,(i,j))]
  if None in labels:return
  if mode=='right':
   assert 2 in labels;counts['right_contact_angle_pass']+=1
  elif labels!=[3,3]:return
  pij=q*c(B[i][j])-a[i]*a[j]
  det=nn[i]*nn[j]*nl+2*pij*ps[0]*ps[1]-nn[i]*ps[1]**2-nn[j]*ps[0]**2-nl*pij**2
  if K.sign(det)<0:return
  records.append({'pair':[i,j],'old_g':[gi,gj],'ambient_norm':norm,'mode':mode,'T':K.encode(v),'labels':labels,'projected_determinant':K.encode(det),'identity_solution':identity(v,i,j,gi,gj)})
 for norm in (2,4):
  for i,j in combinations(range(3),2):
   di,dj,g=B[i][i],B[j][j],B[i][j]
   for gi,gj in product(range(-math.isqrt(int(di*norm)),math.isqrt(int(di*norm))+1),range(-math.isqrt(int(dj*norm)),math.isqrt(int(dj*norm))+1)):
    old=S.Matrix([[di,g,gi],[g,dj,gj],[gi,gj,norm]])
    if any(old.extract(ids,ids).det()<0 for n in (1,2,3)for ids in combinations(range(3),n)):continue
    counts['old_PSD3_integer_cases']+=1
    for ii,gg in ((i,gi),(j,gj)):
     counts['right_contact_trials']+=1;evaluate(q*c(gg)/a[ii],i,j,gi,gj,norm,'right')
    def polynomial(ii,gg):return K.poly([4*q*q*c(gg*gg)-nn[ii]*q*c(norm),-8*q*a[ii]*c(gg),4*a[ii]**2+nn[ii]])
    gcd=polynomial(i,gi).gcd(polynomial(j,gj)).monic()
    old_cases.append((i,j,norm,gi,gj,F(old.det()),gcd))
    if gcd.degree()==1:
     counts['both_label3_common_linear']+=1
     v=-gcd[(0,)];evaluate(v,i,j,gi,gj,norm,'both3')
    elif gcd.degree()==2:counts['both_label3_common_quadratic']+=1;quad.append((i,j,gi,gj,norm))
    else:assert gcd.degree()==0
 counts['both_label3_common_quadratic']+=0
 admissible=[r for r in records if r['identity_solution']and r['identity_solution']['admissible']]
 unique={(tuple(r['pair']),r['ambient_norm'],tuple(r['old_g']),tuple(r['T']))for r in admissible}
 counts['distinct_necessary_survivor_families']=len(unique)
 counts['seed_cusp_edge_germs']=edge_count;counts['unclosable_edge_germs']=edge_count*2//3
 assert quad==[] and dict(counts)==saved['counts'],(quad,dict(counts),saved['counts'])
 assert len(old_cases)==len(saved['old_PSD3_cases'])==120
 for (i,j,norm,gi,gj,d,g),r in zip(old_cases,saved['old_PSD3_cases']):
  assert r['pair']==[i,j] and r['ambient_norm']==norm and r['old_g']==[gi,gj] and F(r['old_Gram_determinant'])==d
  assert [g[(ii,)]for ii in range(g.degree()+1)]==[K.asc(x)for x in r['both3_quadratic_gcd']]
 assert len(records)==len(saved['angle_pass_records'])==15
 for r,s in zip(records,saved['angle_pass_records']):
  for key in ('pair','old_g','ambient_norm','mode','labels','identity_solution'):assert r[key]==s[key],(key,r,s)
  assert K.read(r['T'])==K.asc(s['T']) and K.read(r['projected_determinant'])==K.asc(s['projected_determinant'])
 assert len(admissible)==5 and len(unique)==3 and all(r['pair']==[0,1]for r in admissible)
 assert saved['seed_cusp']['seed_walls']==sorted(seed) and saved['seed_cusp']['projected_span_rank']==20 and saved['seed_cusp']['projected_Gram_rank']==19
 assert saved['unclosable_triangle_pairs']==[[23,200],[50,200]]
 paths=[SOURCE,LOCAL,RAW,WALK,f'scripts/{TAG}.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return {'status':'PASS_INDEPENDENT_FIXED_QUARTIC_CAP_ENDPOINT_GATE','sources':{p:sha(p)for p in paths},'counts':dict(counts),'components':components,'edge_germs':252,'unclosable_edge_germs':168,'unclosable_triangle_pairs':[[23,200],[50,200]],'angle_pass_records':records,'scope':'Independent necessary-condition enumeration for integral ambient norm2/4 normals at local record0,u=1/2, retaining the22 seed cusp walls. No producer or project arithmetic imports; uses SymPy algebraic fields and rational real-root isolation. Five admissible records are necessary conditions only, not constructed caps. No all-normal or other-section exclusion.'}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 out=run();path=ROOT/'results'/f'{TAG}_review.json'
 if args.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['counts'])
