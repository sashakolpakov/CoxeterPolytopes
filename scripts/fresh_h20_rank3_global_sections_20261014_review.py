#!/usr/bin/env python3
"""Separate algebraic-field audit of all91 rank3 global ridge certificates.

Uses SymPy ANP/QQ for quotient-field operations, instead of the producer's
Fraction polynomial/Sturm implementation. Imports no project arithmetic.
"""
import argparse,hashlib,json,math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as S
from sympy.polys.domains import QQ
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_rank3_global_sections_20261014'
SRC=f'results/{TAG}.json'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
LOCAL_REVIEW='results/fresh_h20_rank3_omitted_sections_20261014_review.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
t=S.symbols('t')
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def frac(q):return F(int(q.numerator),int(q.denominator))
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def pd(A):
 A=[[F(x)for x in r]for r in A]
 while A:
  d=A[0][0];assert d>0
  A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,len(A))]for i in range(1,len(A))]

class Field:
 def __init__(self,poly,iv):
  self.p=S.Poly.from_list(list(map(S.Rational,poly)),t)
  assert self.p.is_irreducible
  self.lo,self.hi=map(F,iv)
  assert self.p.count_roots(S.Rational(self.lo),S.Rational(self.hi))==1
  self.K=QQ.algebraic_field((self.p.as_expr(),t));self.zero=self.K.zero;self.one=self.K.one
  self.cache={self.zero:0}
 def scalar(self,a):
  a=F(a);return self.K.new([QQ(a.numerator,a.denominator)])
 def desc(self,p):return self.K.new([QQ(F(a).numerator,F(a).denominator)for a in p])
 def asc(self,p):return self.desc(list(reversed(p)))
 def sign(self,a):
  if a in self.cache:return self.cache[a]
  cs=list(map(frac,a.to_list()));lo,hi=self.lo,self.hi
  for digits in (20,40,80,160,320,640):
   l=h=F(0)
   for c in cs:
    v=(l*lo,l*hi,h*lo,h*hi);l,h=min(v)+c,max(v)+c
   if l>0:self.cache[a]=1;return 1
   if h<0:self.cache[a]=-1;return -1
   lo,hi=map(F,self.p.refine_root(S.Rational(lo),S.Rational(hi),eps=S.Rational(1,10**digits)))
  raise AssertionError('unresolved algebraic sign')
 def linear(self,v,c):return sum((x*self.scalar(y)for x,y in zip(v,c,strict=True)),self.zero)
 def trim(self,p):
  p=list(p)
  while len(p)>1 and p[-1]==self.zero:p.pop()
  return p
 def pmul(self,a,b):
  p=[self.zero]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):p[i+j]+=x*y
  return self.trim(p)
 def psub(self,a,b):return self.trim([(a[i]if i<len(a)else self.zero)-(b[i]if i<len(b)else self.zero)for i in range(max(len(a),len(b)))])
 def peval(self,p,v):
  a=self.zero;v=self.scalar(v)
  for c in reversed(p):a=a*v+c
  return a
 def positive(self,p,cert):
  reduced=list(map(self.asc,cert['reduced_polynomial']));q=reduced
  for _ in range(cert['zero_order_at_0']):q=self.pmul(q,[self.zero,self.one])
  for _ in range(cert['zero_order_at_1']):q=self.pmul(q,[self.one,-self.one])
  assert self.trim(q)==self.trim(p)
  end=F(0)
  for row in cert['positive_Bernstein_intervals']:
   lo,hi=map(F,row['interval']);assert lo==end and lo<hi;end=hi
   c=[self.zero]*len(reduced)
   for i,a in enumerate(reduced):
    for j in range(i+1):c[j]+=a*self.scalar(F(math.comb(i,j))*lo**(i-j)*(hi-lo)**j)
   n=len(c)-1
   bern=[self.linear(c[:k+1],[F(math.comb(k,i),math.comb(n,i))for i in range(k+1)])for k in range(n+1)]
   assert bern==list(map(self.asc,row['Bernstein'])) and all(self.sign(a)>0 for a in bern)
  assert end==1

def run():
 saved=json.loads((ROOT/SRC).read_text());local=json.loads((ROOT/LOCAL).read_text());walk=json.loads((ROOT/WALK).read_text());lr=json.loads((ROOT/LOCAL_REVIEW).read_text())
 for doc in (saved,local,lr,walk):
  for p,h in doc['sources'].items():assert sha(p)==h
 assert saved['status']=='PASS_ALL91_RANK_THREE_FULL_PARAMETER_GLOBAL_EXCLUSIONS' and saved['pending_endpoint_rows']==[]
 assert len(saved['families'])==91 and sorted(r['row']for r in saved['families'])==list(range(91))
 assert lr['positive_directions']=={'A2':16,'C2':19,'G2':56}
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[[dot(q,w)for q in Q]for w in W];G=[[dot(w,c)for c in C]for w in W]
 h=[sum(w[j]for w in W)for j in range(22)];gg=math.gcd(*h);h=[F(a,gg)for a in h]
 pair=lambda a,b:dot(a,[dot(q,b)for q in Q])
 assert pair(h,h)<0 and all(dot(h,c)<0 for c in C)
 endpoints={r['row']:r for r in saved['separate_endpoint_certificates']};assert len(endpoints)==22
 counts=Counter();needed=[];records=[]
 for r in saved['families']:
  ix=r['row'];row=local['directions'][ix];s=local['systems'][row['system']]
  assert r['minimal_polynomial']==row['minimal_polynomial'] and r['root_interval']==row['root_interval']
  K=Field(row['minimal_polynomial'],row['root_interval']);Z=K.zero
  e=list(map(K.desc,row['normal_constant']));q=K.desc(row['normal_square']);L=K.desc(row['parameter_lower_magnitude']);delta=list(map(F,row['null_shift']));k=list(map(K.desc,row['positive_kernel']));R=s['retained_walls'];om=s['omitted_wall']
  mixed=lambda v,p:K.linear(v,[dot(a,p)for a in Q])
  assert q==K.asc(r['normal_square']) and L==K.asc(r['parameter_lower_magnitude']) and K.sign(q)>0 and K.sign(L)>0
  assert all(K.sign(a)>0 for a in k)
  assert e==[K.linear(k,[W[i][j]for i in R])for j in range(22)]
  aa=[mixed(e,W[i])for i in R];assert sum((a*b for a,b in zip(aa,k)),Z)==q and mixed(e,delta)==Z
  assert pair(delta,delta)==0 and pair(delta,h)<0
  factor=next(f for c in walk['corner_orbits']if c['id']==s['cusp']for f in c['affine_components']if f['walls']==s['factor'])
  marks=dict(zip(factor['walls'],map(F,factor['positive_kernel']),strict=True))
  assert delta==[sum(marks[i]*W[i][j]for i in marks)for j in range(22)]
  coeff=[a-L*K.scalar(marks[i])for a,i in zip(k,R)]
  assert coeff==list(map(K.asc,r['u1_redundancy_coefficients'])) and all(K.sign(a)>=0 for a in coeff) and Z in coeff
  e1=[a-L*K.scalar(d)for a,d in zip(e,delta)]
  assert e1==[K.linear(coeff,[W[i][j]for i in R])-L*K.scalar(marks[om]*W[om][j])for j in range(22)]
  def face_pair(record):
   points=[];sets=[]
   for j in (0,1):
    fr=record[f'face_at_{j}'];p=list(map(F,fr['point']));v=[dot(p,c)for c in C];ids=[i for i,a in enumerate(v)if a==0]
    assert ids==fr['walls'] and all(a<=0 for a in v) and max(a for a in v if a)==F(fr['max_other_pairing'])<0
    assert pair(p,p)==F(fr['norm'])<0 and pair(p,h)==F(fr['future'])<0
    pd([[G[a][b]for b in ids]for a in ids])
    d=-pair(p,delta);assert d>0;pn=[a/d for a in p]
    assert pn==list(map(F,record[f'normalized_point_{j}'])) and pair(pn,delta)==-1
    assert mixed(e,pn)==(Z if j==0 else -L)
    points.append(pn);sets.append(set(ids));counts['rational_face_points']+=1
   common=sorted(sets[0]&sets[1]);assert common==record['common_zero_set'] and sets[1]-set(common)=={om}
   assert set(record['ridge'])<=set(common) and pair(*points)<0
   pd([[G[a][b]for b in common]for a in common])
   assert all(c==Z or i in common for i,c in zip(R,coeff))
   return points,common
  points,common=face_pair(r)
  assert pair(*points)==F(r['point_cross_pair'])
  i,j=r['ridge'];a=[mixed(e,W[ii])for ii in (i,j)];b=[-L*K.scalar(pair(delta,W[ii]))for ii in (i,j)]
  pp=K.psub([q*K.scalar(G[i][j])],K.pmul([a[0],b[0]],[a[1],b[1]]))
  ni=K.psub([q*K.scalar(G[i][i])],K.pmul([a[0],b[0]],[a[0],b[0]]));nj=K.psub([q*K.scalar(G[j][j])],K.pmul([a[1],b[1]],[a[1],b[1]]))
  assert pp==list(map(K.asc,r['scaled_pairing_polynomial'])) and ni==list(map(K.asc,r['scaled_first_norm_polynomial'])) and nj==list(map(K.asc,r['scaled_second_norm_polynomial']))
  tests=r['sign_certificates'];K.positive(ni,tests['first_norm']);K.positive(nj,tests['second_norm']);counts['positive_polynomials']+=2
  NN=K.pmul(pp,pp);DD=K.pmul(ni,nj)
  if r['mode']=='obtuse':K.positive(pp,tests['positive_pairing']);counts['positive_polynomials']+=1;good=K.sign(K.peval(pp,1))>0
  else:
   lo,hi=map(F,r['gap']);assert (lo,hi)in ((F(0),F(1,4)),(F(1,4),F(1,2)))
   low=K.psub(NN,[a*K.scalar(lo)for a in DD]);high=K.psub([a*K.scalar(hi)for a in DD],NN)
   for p,name in (([-a for a in pp],'negative_pairing'),(low,'above_lower'),(high,'below_upper')):K.positive(p,tests[name]);counts['positive_polynomials']+=1
   good=K.sign(K.peval(low,1))>0 and K.sign(K.peval(high,1))>0
  boundary=all(c==Z for c in coeff)
  assert boundary==r['u1_is_old_facet_hyperplane'] and good==r['u1_bad_angle_verified']
  if boundary:counts['old_facet_endpoint_degeneracies']+=1
  elif good:
   # Actuality follows from the same face pair plus redundant omitted wall.
   assert K.sign(K.peval(ni,1))>0 and K.sign(K.peval(nj,1))>0
   if r['mode']=='gap':assert K.sign(K.peval(pp,1))<0
   counts['same_pair_endpoint_exclusions']+=1
  else:
   needed.append(ix);er=endpoints[ix];epoints,ec=face_pair(er)
   assert mixed(e1,epoints[1])==Z and mixed(e1,epoints[0])==L
   ii,jj=er['ridge'];ae=[mixed(e1,W[v])for v in (ii,jj)]
   p=q*K.scalar(G[ii][jj])-ae[0]*ae[1];n0=q*K.scalar(G[ii][ii])-ae[0]**2;n1=q*K.scalar(G[jj][jj])-ae[1]**2
   assert p==K.asc(er['scaled_pairing']) and n0==K.asc(er['scaled_first_norm']) and n1==K.asc(er['scaled_second_norm'])
   nn=p*p;dd=n0*n1;assert nn==K.asc(er['squared_cosine_numerator']) and dd==K.asc(er['squared_cosine_denominator'])
   assert K.sign(n0)>0 and K.sign(n1)>0 and K.sign(dd-nn)>0
   if er['mode']=='obtuse':assert K.sign(p)>0
   else:
    lo,hi=map(F,er['gap']);assert (lo,hi)in ((F(0),F(1,4)),(F(1,4),F(1,2)))
    assert K.sign(p)<0 and K.sign(nn-dd*K.scalar(lo))>0 and K.sign(dd*K.scalar(hi)-nn)>0
   assert coeff==list(map(K.asc,er['redundancy_coefficients']))
   counts['endpoint_'+er['mode']]+=1
  counts['families']+=1;counts['degree'+str(K.p.degree())]+=1;counts[r['mode']]+=1
  records.append({'row':ix,'field_degree':K.p.degree(),'continuous_mode':r['mode'],'endpoint':'old_facet'if boundary else 'same_pair'if good else 'separate_pair'})
 assert sorted(needed)==sorted(endpoints) and counts['families']==91
 assert counts['obtuse']==75 and counts['gap']==16 and counts['endpoint_obtuse']==4 and counts['endpoint_gap']==18
 paths=[SRC,LOCAL,LOCAL_REVIEW,RAW,WALK,f'scripts/{TAG}.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return {'status':'PASS_INDEPENDENT_ALL91_GLOBAL_ALGEBRAIC_RIDGE_REVIEW','sources':{p:sha(p)for p in paths},'counts':dict(counts),'records':records,
         'scope':'Independent saved-witness reconstruction with SymPy rational algebraic fields and rational root isolation; no producer/project arithmetic imports. Checks all91 intervals, all22 endpoint followups and actual-ridge premises. Local91 completeness is the separately audited, source-bound dependency.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}_review.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['counts'])
