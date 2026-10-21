#!/usr/bin/env python3
"""Independent sparse-radical reconstruction of42 replacement recipes.

No producer/project arithmetic imports. General squarefree monomials and
outward rational enclosures also handle the extra sqrt3 and norm radicals.
"""
import argparse,json,math
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from hashlib import sha256
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_replacement_escape_20261014'
PRIMARY=f'results/{TAG}.json';BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
PARENT='results/fresh_h20_e7_intrinsic_caps_20261014.json';PROP='results/fresh_h20_e7_causal_residuals_20261014_interior_scout.json'
CONT=f'results/{TAG}_continuous.json'
RAD=[1,2,79,158,11,22,869,1738]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for s,h in d.get('sources',{}).items():assert sha(s)==h,s
 return d
def c(x):return{1:F(x)}if x else{}
def add(a,b):
 z=a.copy()
 for k,v in b.items():
  z[k]=z.get(k,F(0))+v
  if not z[k]:del z[k]
 return z
def scale(a,s):return{k:v*s for k,v in a.items()if v*s}
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
 z={}
 for r,x in a.items():
  for s,y in b.items():
   g=math.gcd(r,s);k=(r//g)*(s//g);z[k]=z.get(k,F(0))+g*x*y
 return{k:v for k,v in z.items()if v}
def total(xs):
 z={}
 for x in xs:z=add(z,x)
 return z
def parse(a):
 assert len(a)==8
 return{r:F(v)for r,v in zip(RAD,a)if F(v)}
def ext(a):assert len(a)==2;return add(parse(a[0]),mul({3:F(1)},parse(a[1])))
def vector(a):assert len(a)==21;return list(map(parse,a))
def enc(a):return{str(k):str(v)for k,v in sorted(a.items())}
@lru_cache(None)
def root_bounds(n,bits):
 D=1<<bits;x=math.isqrt(n*D*D)
 return F(x,D),F(x+int(x*x!=n*D*D),D)
def bounds(a,bits):
 lo=hi=F(0)
 for n,q in a.items():
  x,y=root_bounds(n,bits)
  lo+=q*(x if q>=0 else y);hi+=q*(y if q>=0 else x)
 return lo,hi
def sign(a):
 if not a:return 0
 for b in [128,256,512,1024]:
  l,u=bounds(a,b)
  if l>0:return 1
  if u<0:return -1
 raise AssertionError(('unresolved field sign',a))
def sqrt_interval(lo,hi,bits):
 assert lo>0;D=1<<bits
 a=math.isqrt((lo.numerator*D*D)//lo.denominator)
 b=math.isqrt((hi.numerator*D*D)//hi.denominator)
 return F(a,D),F(b+int(F(b*b,D*D)<hi),D)
def nested_sign(a,b,q):
 # No assumption that sqrt(q) is independent of the old field.
 for bits in [128,256,512,1024]:
  l,u=bounds(q,bits)
  if l<=0:continue
  rl,ru=sqrt_interval(l,u,bits);al,au=bounds(a,bits);bl,bu=bounds(b,bits)
  v=[x*y for x in [bl,bu]for y in [rl,ru]];lo=al+min(v);hi=au+max(v)
  if lo>0:return 1,bits
  if hi<0:return -1,bits
 raise AssertionError(('unresolved nested sign',a,b,q))
def inverse(A):
 n=len(A);B=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for k in range(n):
  p=next(i for i in range(k,n)if B[i][k]);B[k],B[p]=B[p],B[k];q=B[k][k];B[k]=[x/q for x in B[k]]
  for i in range(n):
   if i!=k:
    q=B[i][k];B[i]=[x-q*y for x,y in zip(B[i],B[k])]
 return[r[n:]for r in B]
def rational(a):assert set(a)<={1};return a.get(1,F(0))
def pd(A):
 B=[r[:]for r in A]
 for k in range(len(B)):
  q=B[k][k];assert q>0
  for i in range(k+1,len(B)):
   for j in range(k+1,len(B)):B[i][j]-=B[i][k]*B[k][j]/q

def run():
 p,base,parent,proposals=map(read,[PRIMARY,BASE,PARENT,PROP]);Q=[list(map(F,r))for r in base['form']]
 N=list(map(vector,base['normals']));assert len(N)==87 and base['normals'][:85]==parent['normals']
 H=vector(base['strict_center']);oldE=vector(base['strict_escape']['point']);E=vector(p['point'])
 def cov(v):return[total(scale(x,a)for x,a in zip(v,row)if a)for row in Q]
 def dot(v,w):return total(mul(a,b)for a,b in zip(v,cov(w))if a and b)
 def linear(v,co):return[total(mul(a[k],q)for a,q in zip(v,co))for k in range(21)]
 def vs(v,s):return[mul(x,s)for x in v]
 S=base['E7_support'];A=parent['face_facet_walls'][:13]
 assert S==[12,13,14,15,18,31,65]and A==list(range(11))+[17,79]
 assert all(dot(N[i],N[i])==c(2)for i in S+A)
 GS=[[rational(dot(N[i],N[j]))for j in S]for i in S]
 GA=[[rational(dot(N[i],N[j]))for j in A]for i in A]
 pd(GS);pd(GA);I=inverse(GS);J=inverse(GA)
 W=[linear([N[i]for i in S],list(map(c,row)))for row in I]
 L=[linear([N[i]for i in A],list(map(c,row)))for row in J]
 assert all(dot(w,N[j])==c(int(i==k))for i,w in enumerate(W)for k,j in enumerate(S))
 assert all(dot(l,N[j])==c(int(i==k))for i,l in enumerate(L)for k,j in enumerate(A))
 assert all(not dot(w,l)for w in W for l in L)
 T=vs(vector(parent['face_basis'][-1]),{2:F(1,2)})
 assert dot(T,T)==c(-2)and sign(dot(T,H))<0
 assert all(not dot(T,z)for z in W+L)
 assert dot(E,E)==parse(p['norm'])and sign(dot(E,E))>0
 assert dot(E,H)==parse(p['center_pairing'])and sign(dot(E,H))<0 and sign(dot(H,H))<0
 retained=[i for i in range(87)if i not in[61,75]];assert retained==p['retained_base_walls']
 ep=[dot(E,N[i])for i in retained];assert ep==list(map(parse,p['retained_wall_pairings']))
 assert all(sign(x)<0 for x in ep)and all(sign(dot(H,N[i]))<0 for i in retained)
 values=[c(1),{2:F(1)},{3:F(1)},c(2),c(3),c(4)]
 rows=[];radicands=Counter();maxbits=0
 candidates=proposals['delete61_75_candidates'];assert len(candidates)==42 and len(p['records'])==42
 for index,(r,saved)in enumerate(zip(candidates,p['records'])):
  si=r['E7_support'];ii,cc=r['source_recipe']
  assert len(si)==3 and len(set(si))==3 and all(0<=j<7 for j in si)
  assert r['E7_labels']==[3]*3 and len(ii)==len(cc)and len(set(ii))==len(ii)
  assert all(0<=i<13 for i in ii)and all(0<=j<6 for j in cc)
  d=[{}for _ in range(13)]
  for i,j in zip(ii,cc):d[i]=values[j]
  B=linear(W+L,[c(-int(i in si))for i in range(7)]+[scale(x,-1)for x in d])
  assert not dot(B,T)
  q=scale(sub(dot(B,B),c(2)),F(1,2));assert sign(q)>0 and q==ext(saved['time_coefficient_squared'])
  assert add(dot(B,B),mul(q,dot(T,T)))==c(2)
  assert [dot(B,N[j])for j in S]==[c(-int(i in si))for i in range(7)]
  assert [dot(B,N[j])for j in A]==[scale(x,-1)for x in d]
  rec={'index':index,'time_squared':enc(q)}
  for key,v in [('escape_pairing',E),('center_pairing',H)]:
   a=dot(B,v);b=dot(T,v);sg,bits=nested_sign(a,b,q);maxbits=max(maxbits,bits)
   assert sg==-1 and sg==saved[key]['sign']
   assert a==ext(saved[key]['constant'])and b==ext(saved[key]['radical_coefficient'])
   rec[key]=dict(constant=enc(a),radical_coefficient=enc(b),sign=sg,bits=bits)
  assert saved['index']==index and saved['E7_support']==si and saved['source_indices']==ii and saved['source_value_indices']==cc
  radicands['rational'if set(q)<={1}else'multiquadratic']+=1
  if index==34:
   assert si==[2,4,6]and ii==[0,12]and cc==[1,3]and q==c(25)
   cap=[add(a,scale(b,5))for a,b in zip(B,T)]
  rows.append(rec)
 selected=p['selected_cap'];assert selected['proposal_index']==34 and cap==vector(selected['normal'])
 assert dot(cap,cap)==c(2)==parse(selected['norm'])
 assert dot(cap,oldE)==parse(selected['old_escape_pairing'])and sign(dot(cap,oldE))>0
 assert dot(cap,E)==parse(selected['new_escape_pairing'])and sign(dot(cap,E))<0
 continuous=read(CONT)
 allpairs=[dot(E,n)for n in N]
 assert all(sign(a)<0 for a in allpairs)
 assert allpairs==list(map(parse,continuous['all87_pairings']))
 V=linear([W[2],W[4],W[6],L[0]],[c(-1),c(-1),c(-1),{2:F(-1)}])
 assert V==vector(continuous['family']['constant_vector'])
 assert L[12]==vector(continuous['family']['source_weight'])
 assert T==vector(continuous['family']['time_vector'])
 assert dot(V,V)==c(F(148,3))and dot(L[12],L[12])==c(F(2,3))
 assert not dot(V,L[12])and not dot(V,T)and not dot(L[12],T)
 C=dot(V,E);aa=scale(dot(L[12],E),-1);BB=scale(dot(T,E),-1)
 DD=sub(mul(BB,BB),scale(mul(aa,aa),3))
 margin=sub(scale(DD,F(71,3)),mul(C,C))
 assert all(sign(z)>0 for z in [C,aa,BB,DD,margin])
 cb=continuous['bound']
 assert [C,aa,BB,DD,margin]==[parse(cb[k])for k in ['C','a','B','D=B2-3a2','positive_gap=71D/3-C2']]
 paths=[PRIMARY,BASE,PARENT,PROP,CONT,f'scripts/{TAG}.py',f'scripts/{TAG}_continuous.py',f'scripts/{TAG}_review.py']
 return dict(status='PASS_INDEPENDENT_42_REPLACEMENT_RECIPE_COMMON_ESCAPE',sources={x:sha(x)for x in paths},
  field_basis_squarefree=RAD,extra_coefficient_radical=3,physical_E7_roots=S,physical_A_roots=A,
  physical_A_norms=[2]*13,retained_wall_signs=85,recipes_rebuilt=42,cap_norm_identities=42,
  strict_new_pairings=84,radicand_classes=dict(radicands),maximum_interval_bits=maxbits,
  escape_norm=enc(dot(E,E)),escape_center_pairing=enc(dot(E,H)),
  selected34=dict(norm=enc(dot(cap,cap)),time_squared={'1':'25'},positive_time_coefficient=5,
   old_escape_pairing=enc(dot(cap,oldE)),new_escape_pairing=enc(dot(cap,E))),records=rows,
  continuous_family=dict(all87_strict_escape_signs=87,constant_norm='148/3',weight_norm='2/3',time_norm='-2',
   C=enc(C),a=enc(aa),B=enc(BB),D=enc(DD),positive_margin=enc(margin),
   theorem='For all real x, r=V-x*L12+sqrt((71+x*x)/3)*T has norm2 and E*r=C+a*x-B*sqrt((71+x*x)/3)<0; the latter follows from B>0,D=B*B-3*a*a>0 and 71*D/3>C*C with C>0.'),
  scope='Independent sparse squarefree Fraction reconstruction of both physical dual bases, all42 norm-two recipes and both signs for each,85 retained-wall escape signs, common center, and cap34. All radical branches use outward rational enclosures, not producer sign routines. The separate continuous profile is checked for all real x by a Lorentz-Cauchy bound. No mutual-Coxeter/facet-actuality or other-profile claim. Every subset of the42 recipes with the retained85 has a common strict future spacelike escape and strict timelike interior.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();p=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['recipes_rebuilt'],d['maximum_interval_bits'])
