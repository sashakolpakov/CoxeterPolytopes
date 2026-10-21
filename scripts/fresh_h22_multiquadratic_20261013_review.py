#!/usr/bin/env python3
"""Independent sparse square-root replay of the pi/6 860-wall construction.

Only rational source data are read.  No producer arithmetic/projection or
interval implementation is imported.  Square-root expressions are collected
in the Q-linearly independent basis sqrt(squarefree positive integer).
"""
import argparse,hashlib,json,math,time
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
P=Path(__file__).resolve().parents[1];TAG='fresh_h22_multiquadratic_20261013';OUT=P/'results'/f'{TAG}_review.json'
def read(s):return json.loads((P/s).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def C(x):return {1:F(x)}if x else{}
def add(a,b):
 z=a.copy()
 for k,v in b.items():
  z[k]=z.get(k,F(0))+v
  if not z[k]:del z[k]
 return z
def scale(a,x):return {k:v*x for k,v in a.items()if v*x}
def mul(a,b):
 z={}
 for r,x in a.items():
  for s,y in b.items():
   g=math.gcd(r,s);k=(r//g)*(s//g);z[k]=z.get(k,F(0))+x*y*g
 return {k:v for k,v in z.items()if v}
def sqrtq(x):
 x=F(x);assert x>0;n=x.numerator*x.denominator;original=n;sq=sf=1;p=2
 while p*p<=n:
  e=0
  while n%p==0:n//=p;e+=1
  sq*=p**(e//2)
  if e%2:sf*=p
  p=3 if p==2 else p+2
 sf*=n;assert sq*sq*sf==original;return {sf:F(sq,x.denominator)}
@lru_cache(None)
def interval_root(n,digits):
 D=10**digits;r=math.isqrt(n*D*D);return F(r,D),F(r+(r*r<n*D*D),D)
def sign(a):
 if not a:return 0
 for digits in(50,100,200):
  lo=hi=F(0)
  for n,c in a.items():
   x,y=interval_root(n,digits)
   if c<0:x,y=y,x
   lo+=c*x;hi+=c*y
  if lo>0:return 1
  if hi<0:return -1
 raise AssertionError(('unresolved algebraic sign',a))
def psum(xs):
 z={}
 for x in xs:z=add(z,x)
 return z
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mv(M,a):return [psum(scale(x,F(c))for x,c in zip(a,row)if c)for row in M]
def pair(a,covb):return psum(mul(x,y)for x,y in zip(a,covb)if x and y)
def ratvec(a):return list(map(C,a))
def vadd(a,b):return [add(x,y)for x,y in zip(a,b)]
def vs(a,c):return [scale(x,c)for x in a]
def vm(a,c):return [mul(x,c)for x in a]
def flip(a):return a[:-1]+[scale(a[-1],-1)]
def enc(a):return {str(k):str(v)for k,v in sorted(a.items())}
def pd(A):
 A=[list(map(F,r))for r in A]
 for k in range(len(A)):
  d=A[k][k];assert d>0
  for i in range(k+1,len(A)):
   for j in range(i,len(A)):
    A[i][j]-=A[i][k]*A[j][k]/d;A[j][i]=A[i][j]
def replay():
 started=time.monotonic();d=read('results/'+TAG+'.json')
 for p,s in d['sources'].items():assert sha(p)==s,p
 rawpath='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';corepath='results/fresh_h22_incidence_core_20261013.json';bankpath='data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json';recpath='data/enumerations/fresh_h22_incidence_prism_20261013_recipes.json';coupledpath='results/fresh_h22_pi7_coupled_caps_20261013.json'
 raw=read(rawpath)['coxeter_input'];Q=[[-int(x)for x in row]for row in raw['basis_gram']];M=[r+[0]for r in Q]+[[0]*22+[2]];W=[list(map(int,row))+[0]for row in raw['wall_coordinates_as_columns_transposed']];core=read(corepath);points=core['incidence_parts'][0];vp=core['old_projections'][1];old=[]
 for r in read(bankpath)['candidates']:
  if r['twice_height_square']!=4:continue
  u=[(r['mask']>>i)&1 for i in range(21)];assert sum(u)==r['weight'];v=[r['A']*vp[j]-sum(u[i]*W[points[i]][j]for i in range(21))for j in range(22)];assert v==r['twice_old_projection']and all(x%2==0 for x in v);old.append([x//2 for x in v])
 assert len(old)==301;R=W+[r+[1]for r in old]+[r+[-1]for r in old];h=core['center']+[0];Rc=[ [dot(row,r)for row in M]for r in R];hc=[dot(row,h)for row in M];rn=[dot(r,c)for r,c in zip(R,Rc)];assert Counter(rn)=={2:644,4:168};assert dot(h,hc)==-14 and max(dot(r,hc)for r in R)<0
 recs=read(recpath)['recipes'];coupled=read(coupledpath);types=d['selected_original_types']+[23,24];norms={};normals={};rad={};T=sqrtq(3);original={}
 for k,r in enumerate(recs):
  u=list(map(F,r['u']));v=list(map(F,r['v']));q=F(r['norm']);uc=[dot(row,u)for row in M];vc=[dot(row,v)for row in M];qu=dot(u,uc);qv=dot(v,vc);assert q>0 and qv>0 and qu<=0 and dot(u,vc)==0
  rr=q/qv-3*q*q*qu/(4*qv);rad[k]=rr;normal=[add(scale(T,-q*x/2),scale(sqrtq(rr),r['radical_branch']*y))for x,y in zip(u,v)];assert pair(normal,mv(M,normal))==C(q);original[k]=normal;normals[k]=normal if sign(normal[-1])>0 else flip(normal);norms[k]=q
 a=list(map(F,recs[17]['u']));b=list(map(F,recs[17]['v']));assert norms[17]==2;Z=[add(scale(sqrtq(31),x/3),scale(T,-3*y/2))for x,y in zip(a,b)];Z=Z if sign(pair(Z,ratvec(hc)))<0 else vs(Z,-1);assert pair(Z,mv(M,Z))==C(-2)and pair(Z,mv(M,normals[17]))=={}
 u,v=[list(map(F,x))for x in coupled['rational_duals']];normals[23]=vs(vadd(vs(Z,209),ratvec([-241*x for x in u])),F(1,120));normals[24]=vm(vadd(vs(Z,241),ratvec([-418*x for x in v])),scale(sqrtq(29281),F(1,29281)));norms[23]=norms[24]=F(2)
 record={r['type']:r for r in d['cap_records']};caplist=[];qq=[];idx=[]
 parentpairs=Counter()
 for k in types:
  r=normals[k];q=norms[k];assert pair(r,mv(M,r))==C(q)and sign(pair(r,ratvec(hc)))<0 and sign(r[-1])>0
  ff={int(i):m for i,m in record[k]['finite_parent_labels'].items()}
  for j,c in enumerate(Rc):
   p=pair(r,ratvec(c));assert sign(p)<=0
   if j in ff:
    m=ff[j];expected={2:F(0),4:F(1,2),6:F(3,4)}[m];assert mul(p,p)==C(q*rn[j]*expected);assert (not p)==(m==2);parentpairs[str(m)]+=1
   else:assert sign(add(mul(p,p),C(-q*rn[j])))>0;parentpairs['dotted']+=1
  caplist.extend((r,flip(r)));qq.extend((q,q));idx.extend((k,k))
 assert len(caplist)==48;capcov=[mv(M,r)for r in caplist];counts=Counter()
 for i,r in enumerate(caplist):
  assert pair(r,capcov[i])==C(qq[i])and sign(pair(r,ratvec(hc)))<0
  for j in range(i):
   p=pair(r,capcov[j]);assert sign(p)<=0
   if not p:kind='right_by_exact_frame_identity'
   elif mul(p,p)==C(qq[i]*qq[j]):kind='parallel_by_exact_null_projection'
   else:assert sign(add(mul(p,p),C(-qq[i]*qq[j])))>0;kind='dotted'
   counts[kind]+=1
 assert dict(counts)==d['cap_pair_counts'];assert counts=={'dotted':1120,'right_by_exact_frame_identity':6,'parallel_by_exact_null_projection':2}
 x=list(map(F,d['strict_escape']['vector']));xc=[dot(row,x)for row in M];assert dot(x,xc)==1586 and dot(x,hc)==-814 and max(dot(r,xc)for r in R)<0
 assert all(sign(pair(r,ratvec(xc)))<0 for r in caplist)
 delta=list(map(F,d['untouched_cusp']['ray']));dc=[dot(row,delta)for row in M];assert dot(delta,dc)==0 and dot(delta,hc)<0 and max(dot(r,dc)for r in R)==0 and all(sign(pair(r,ratvec(dc)))<0 for r in caplist)
 # Four actual points in the new positive-height H2 face.
 S=coupled['old_incident'][:-1];pd([[dot(R[i],Rc[j])for j in S]for i in S]);common=[j for j in S if j not in(27,9)];vertices=[]
 assert len(S)==21 and len(common)==19
 assert all(pair(normals[17],ratvec(Rc[j]))=={}for j in S)
 for k,released in((23,27),(24,9)):
  assert all(pair(normals[k],ratvec(Rc[j]))=={}for j in S if j!=released)
  assert pair(normals[k],mv(M,normals[17]))=={}
 assert pair(normals[23],mv(M,normals[24]))=={}
 # Each of the four incident rank22 systems is a principal PD block of
 # S, plus one/two/three mutually perpendicular positive norm2 roots.
 triples=[(1,0,0),(241,-209,0),(50369,-29281,-28800),(209,0,-241)]
 for A,B,Cc in triples:
  p=vadd(vs(Z,A),ratvec([B*x+Cc*y for x,y in zip(u,v)]));qp=pair(p,mv(M,p));assert qp==C(-2*A*A+2*B*B+Cc*Cc+2*B*Cc)and sign(qp)<0 and sign(pair(p,ratvec(hc)))<0
  pp=[pair(p,ratvec(c))for c in Rc]+[pair(p,c)for c in capcov];assert max(sign(y)for y in pp)==0
  active=[i for i,y in enumerate(pp)if not y];assert len(active)==22 and all(j in active for j in common);vertices.append({'frame':[A,B,Cc],'norm':enc(qp),'active':active})
 for own,saved in zip(vertices,d['compact_H2_vertices']):
  assert own['frame']==saved['frame_coefficients']and own['active']==sorted(saved['incident'])and own['norm']==enc(C(saved['norm']))
 sf=sorted({next(iter(sqrtq(t)))for t in [F(3),F(31),F(29281)]+[rad[k]for k in types if k<23]}-{1})
 assert sf==d['field']['squarefree_generators']
 bad=pair(original[12],ratvec(Rc[40]));assert mul(bad,bad)==C(F(192,49));assert F(4*48,49)-2==F(94,49)
 cycle=d['field']['quadratic_irrational_cycle'];j=cycle['parent_wall'];cy=scale(mul(pair(original[1],ratvec(Rc[j])),pair(original[1],ratvec(Rc[j]))),F(4,norms[1]*rn[j]));assert cy=={1:F(cycle['rational_part']),3:F(cycle['sqrt3_coefficient'])}and cy[3]
 short=[i for i in range(210)if rn[i]==2];assert len(short)==42 and sum(dot(R[i],Rc[j])for i in short for j in short)==-126
 sources={p:sha(p)for p in(rawpath,corepath,bankpath,recpath,coupledpath,'scripts/'+TAG+'_review.py','scripts/'+TAG+'.py','results/'+TAG+'.json')}
 return {'status':'PASS_INDEPENDENT_MULTIQUADRATIC_EXACT_REVIEW','sources':sources,'normals_rebuilt':48,'facets':860,'positive_height_parent_pairs':dict(parentpairs),'all_signed_cap_pairs':dict(counts),'all_860_strict_escape_norm':1586,'all_860_strict_escape_center_pair':-814,'untouched_old_cusp':True,'exact_quadrilateral_vertices':vertices,'irrational_two_cycle':enc(cy),'rational_indefinite_block_ones_value':-126,'scope':'Source/theorem review plus separate squarefree-radical arithmetic replay of all new norms, center signs, parent contacts, new-new angles, four finite vertices, escape and untouched cusp. The existing 812 rational-wall pair certificate is inherited. This does not claim finite volume, exact full tracefield, pseudo-arithmeticity, or novelty.','elapsed_seconds':time.monotonic()-started}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=replay()
 if a.verify:
  old=json.loads(OUT.read_text());assert {k:v for k,v in old.items()if k!='elapsed_seconds'}=={k:v for k,v in d.items()if k!='elapsed_seconds'}
 else:OUT.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps(d,indent=2))
