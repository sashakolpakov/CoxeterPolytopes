#!/usr/bin/env python3
"""Separate Fraction/radical check of the global witness certificate.

Checks face geometry, polynomial positivity certificates and both exceptional
vertices.  Does not repeat the complete local direction enumeration.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_omitted_factor_classification_20261014'
SRC=f'results/{TAG}_global.json';LOCAL=f'results/{TAG}.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
D=(1,2,3,6,5,10,15,30);ZERO=(F(0),)*8;ONE=(F(1),)+(F(0),)*7
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dec(a):return tuple(map(F,a))
def add(a,b):return tuple(x+y for x,y in zip(a,b,strict=True))
def neg(a):return tuple(-x for x in a)
def sub(a,b):return add(a,neg(b))
def scale(a,c):return tuple(x*c for x in a)
def rat(a):return (F(a),)+(F(0),)*7
def mul(a,b):
 out=[F(0)]*8
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:
     g=math.gcd(D[i],D[j]);k=D.index(D[i]*D[j]//(g*g));out[k]+=x*y*g
 return tuple(out)
@lru_cache(None)
def sg(a):
 if a==ZERO:return 0
 for bits in (60,120,240,480):
  den=1<<bits;lo=hi=F(0)
  for c,d in zip(a,D):
   k=math.isqrt(d*den*den);l=F(k,den);h=l if k*k==d*den*den else F(k+1,den)
   lo+=c*(l if c>=0 else h);hi+=c*(h if c>=0 else l)
  if lo>0:return 1
  if hi<0:return -1
 raise AssertionError('unresolved radical sign')
def linear(v,c):
 out=ZERO
 for a,b in zip(v,c,strict=True):out=add(out,scale(a,b))
 return out
def scalar(v,w):
 out=ZERO
 for a,b in zip(v,w,strict=True):out=add(out,mul(a,b))
 return out
def dot(v,w):return sum(a*b for a,b in zip(v,w,strict=True))
def pd(A):
 A=[[F(x)for x in r]for r in A]
 while A:
  p=A[0][0];assert p>0
  A=[[A[i][j]-A[i][0]*A[0][j]/p for j in range(1,len(A))]for i in range(1,len(A))]
def solve(A,b):
 A=[[F(x)for x in r]for r in A];b=list(b);n=len(A)
 for k in range(n):
  pivot=next(i for i in range(k,n)if A[i][k]);A[k],A[pivot]=A[pivot],A[k];b[k],b[pivot]=b[pivot],b[k]
  d=A[k][k];A[k]=[v/d for v in A[k]];b[k]=scale(b[k],1/d)
  for i in range(n):
   if i!=k:
    c=A[i][k];A[i]=[u-c*v for u,v in zip(A[i],A[k])];b[i]=sub(b[i],scale(b[k],c))
 return b
def trim(p):
 p=list(p)
 while len(p)>1 and p[-1]==ZERO:p.pop()
 return p
def pmul(a,b):
 c=[ZERO]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=add(c[i+j],mul(x,y))
 return trim(c)
def padd(a,b):return trim([add(a[i]if i<len(a)else ZERO,b[i]if i<len(b)else ZERO)for i in range(max(len(a),len(b)))])
def psub(a,b):return padd(a,[neg(x)for x in b])
def peval(p,u):
 v=ZERO
 for c in reversed(p):v=add(scale(v,u),c)
 return v
def check_positive(p,proof):
 reduced=list(map(dec,proof['reduced_polynomial']));q=reduced
 for _ in range(proof['zero_order_at_0']):q=pmul(q,[ZERO,ONE])
 for _ in range(proof['zero_order_at_1']):q=pmul(q,[ONE,neg(ONE)])
 assert trim(q)==trim(p)
 previous=F(0)
 for row in proof['positive_Bernstein_intervals']:
  lo,hi=map(F,row['interval']);assert lo==previous and lo<hi;previous=hi
  coeff=[ZERO]*len(reduced)
  for i,c in enumerate(reduced):
   for j in range(i+1):coeff[j]=add(coeff[j],scale(c,F(math.comb(i,j))*lo**(i-j)*(hi-lo)**j))
  n=len(reduced)-1
  bern=[linear(coeff[:k+1],[F(math.comb(k,i),math.comb(n,i))for i in range(k+1)])for k in range(n+1)]
  assert bern==list(map(dec,row['Bernstein'])) and all(sg(a)>0 for a in bern)
 assert previous==1

def run():
 data=json.loads((ROOT/SRC).read_text());local=json.loads((ROOT/LOCAL).read_text())
 for p,h in data['sources'].items():assert sha(p)==h
 for p,h in local['sources'].items():assert sha(p)==h
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[[dot(q,w)for q in Q]for w in W];G=[[dot(w,c)for c in C]for w in W]
 h=[sum(w[j]for w in W)for j in range(22)];g=math.gcd(*h);h=[F(x,g)for x in h]
 pair=lambda a,b:dot(a,[dot(q,b)for q in Q])
 fmixed=lambda a,b:linear(a,[dot(q,b)for q in Q])
 fpair=lambda a,b:scalar(a,[linear(b,q)for q in Q])
 assert pair(h,h)<0 and all(dot(h,c)<0 for c in C)
 endpoint_checks=0;polynomial_checks=0
 for r in data['new_plane_global_certificates']:
  row=local['local_directions'][r['row']];e=list(map(dec,row['unit_normal_constant']));delta=list(map(F,row['null_shift']));L=dec(row['parameter_lower_magnitude'])
  assert fpair(e,e)==ONE and fmixed(e,delta)==ZERO and pair(delta,delta)==0 and sg(L)>0
  points=[];zeros=[]
  for j in (0,1):
   f=r[f'face_at_{j}'];p=list(map(F,f['point']));vals=[dot(p,c)for c in C]
   assert pair(p,p)==F(f['norm'])<0 and pair(p,h)==F(f['future'])<0
   ids=[i for i,a in enumerate(vals)if a==0];assert ids==f['walls'] and all(a<=0 for a in vals)
   assert max(a for a in vals if a)!=0 and max(a for a in vals if a)==F(f['max_other_pairing'])
   pd([[G[a][b]for b in ids]for a in ids])
   point=list(map(F,r[f'normalized_point_{j}']));d=-pair(p,delta);assert d>0 and point==[a/d for a in p]
   assert pair(point,delta)==-1 and fmixed(e,point)==(ZERO if j==0 else neg(L))
   points.append(point);zeros.append(set(ids));endpoint_checks+=1
  K=sorted(zeros[0]&zeros[1]);assert K==r['common_zero_set'] and set(r['ridge'])<=set(K)
  assert zeros[1]-set(K)=={row['omitted_wall']}
  assert pair(*points)==F(r['point_cross_pair'])<0
  pd([[G[a][b]for b in K]for a in K])
  i,j=r['ridge'];a=[fmixed(e,W[k])for k in (i,j)];b=[scale(L,-pair(delta,W[k]))for k in (i,j)]
  pp=psub([rat(G[i][j])],pmul([a[0],b[0]],[a[1],b[1]]))
  ni=psub([rat(G[i][i])],pmul([a[0],b[0]],[a[0],b[0]]));nj=psub([rat(G[j][j])],pmul([a[1],b[1]],[a[1],b[1]]))
  assert pp==list(map(dec,r['pairing_polynomial'])) and ni==list(map(dec,r['first_norm_polynomial'])) and nj==list(map(dec,r['second_norm_polynomial']))
  tests=r['sign_certificates'];check_positive(ni,tests['first_norm']);check_positive(nj,tests['second_norm']);polynomial_checks+=2
  if r['mode']=='obtuse':check_positive(pp,tests['positive_pairing']);polynomial_checks+=1
  elif r['mode']=='gap':
   lo,hi=map(F,r['angle_gap']);NN=pmul(pp,pp);DD=pmul(ni,nj)
   for p,name in [([neg(v)for v in pp],'negative_pairing'),(psub(NN,[scale(v,lo)for v in DD]),'above_lower'),(psub([scale(v,hi)for v in DD],NN),'below_upper')]:check_positive(p,tests[name]);polynomial_checks+=1
  else:assert r['row']==110 and pmul(pp,pp)==[ZERO,ZERO,ONE] and pmul(ni,nj)==[rat(4),ZERO,neg(ONE)]
  coeff=list(map(dec,r['u1_omitted_wall_redundancy_coefficients']));ret=row['retained_walls']
  assert len(coeff)==len(ret) and all(sg(c)>=0 for c in coeff) and all(c==ZERO or i in K for i,c in zip(ret,coeff))
  # Recover the omitted coefficient from the independently checked null relation.
  k=list(map(dec,row['positive_kernel']))
  fac=list(row['factor']);om=row['omitted_wall']
  from_local_e=[linear(k,[W[i][j]for i in ret])for j in range(22)];assert from_local_e==e
  difference=[sub(sub(e[j],scale(L,delta[j])),linear(coeff,[W[i][j]for i in ret]))for j in range(22)]
  pivot=next(j for j,a in enumerate(W[om])if a)
  lam=scale(difference[pivot],F(1,W[om][pivot]));assert sg(lam)<0
  assert difference==[scale(lam,a)for a in W[om]]
  assert r['u1_is_old_facet_hyperplane']==all(c==ZERO for c in coeff)
  if r['u1_bad_angle_verified'] and not r['u1_is_old_facet_hyperplane']:
   pp1,n1,n2=peval(pp,F(1)),peval(ni,F(1)),peval(nj,F(1));assert sg(n1)>0 and sg(n2)>0
   if r['mode']=='obtuse':assert sg(pp1)>0
   else:
    lo,hi=map(F,r['angle_gap']);assert sg(pp1)<0
    assert sg(sub(mul(pp1,pp1),scale(mul(n1,n2),lo)))>0 and sg(sub(scale(mul(n1,n2),hi),mul(pp1,pp1)))>0
 for r in data['isolated_parameter_actual_vertices']:
  point=list(map(dec,r['actual_vertex']));normal=list(map(dec,r['normal']));ids=r['old_edge_walls']
  assert fpair(normal,normal)==ONE and fpair(point,normal)==ZERO
  assert fpair(point,point)==dec(r['norm']) and sg(dec(r['norm']))<0
  assert fmixed(point,h)==dec(r['future_pair']) and sg(dec(r['future_pair']))<0
  vals=[linear(point,c)for c in C];assert vals==list(map(dec,r['all210_pairings']))
  assert all(sg(v)<=0 for v in vals) and [i for i,a in enumerate(vals)if a==ZERO]==sorted(ids)
  A=[[G[a][b]for b in ids]for a in ids];pd(A)
  a=[fmixed(normal,W[i])for i in ids];ce=solve(A,a)
  qz=sub(ONE,scalar(a,ce));assert qz==dec(r['orthogonal_normal_square']) and sg(qz)>0
  i,j=r['bad_pair'];assert i in ids and j in ids
  ai,aj=fmixed(normal,W[i]),fmixed(normal,W[j]);p=sub(rat(G[i][j]),mul(ai,aj));ni=sub(rat(G[i][i]),mul(ai,ai));nj=sub(rat(G[j][j]),mul(aj,aj))
  n=mul(p,p);d=mul(ni,nj)
  assert sg(ni)>0 and sg(nj)>0 and sg(p)<0 and sg(n)>0 and sg(sub(scale(d,F(1,4)),n))>0
  assert p==dec(r['pairing']) and n==dec(r['squared_cosine_numerator']) and d==dec(r['squared_cosine_denominator'])
 paths=[SRC,LOCAL,RAW,f'scripts/{TAG}_global.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return {'status':'PASS_SEPARATE_GLOBAL_WITNESS_REVIEW','sources':{p:sha(p)for p in paths},
         'rational_face_points':endpoint_checks,'continuous_families':63,'positive_polynomial_certificates':polynomial_checks,
         'isolated_actual_vertices':2,
         'scope':'Separate Fraction/radical reconstruction of63 interpolation families, all126 ambient face points, all saved interval-polynomial positivity certificates, endpoint redundant-wall identities and2 exceptional vertices. Source/proof review covers actual ridge reasoning and old-plane transfer. Does not repeat local classification or prior global angle theorem.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run()
 p=ROOT/'results'/f'{TAG}_global_review.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['positive_polynomial_certificates'])
