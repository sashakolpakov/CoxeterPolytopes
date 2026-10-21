#!/usr/bin/env python3
"""Exact global ridge witnesses for the frozen rank-three local census.

Standard-library rational polynomial arithmetic. The local enumeration is an
explicit source dependency; no producer or numerical modules are imported.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
import argparse,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_rank3_global_sections_20261014'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
REC=f'results/{TAG}_recipes.json'

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 return x

def trim(p):
 p=list(p)
 while len(p)>1 and not p[-1]:p.pop()
 return p

def rpadd(a,b):return trim([(a[i]if i<len(a)else F(0))+(b[i]if i<len(b)else F(0))for i in range(max(len(a),len(b)))])
def rpscale(a,c):return trim([x*c for x in a])
def rpmul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return trim(out)
def rpdiv(a,b):
 a=trim(a);b=trim(b);assert b!=[0];q=[F(0)]*max(1,len(a)-len(b)+1)
 while a!=[0]and len(a)>=len(b):
  i=len(a)-len(b);c=a[-1]/b[-1];q[i]+=c
  for j,x in enumerate(b):a[i+j]-=c*x
  a=trim(a)
 return trim(q),a

def rpeval(p,x):
 y=F(0)
 for c in reversed(p):y=y*x+c
 return y

def variations(seq,x):
 signs=[1 if z>0 else -1 for p in seq if(z:=rpeval(p,x))]
 return sum(a!=b for a,b in zip(signs,signs[1:]))

class Field:
 def __init__(self,poly,interval):
  self.p=list(reversed(list(map(F,poly))));assert self.p[-1]==1
  self.d=len(self.p)-1;assert self.d in[1,2,4]
  self.lo,self.hi=map(F,interval)
  self.zero=(F(0),)*self.d;self.one=self.f(1)
  if self.d==1:
   assert self.lo==self.hi==-self.p[0]
  else:
   assert self.lo<self.hi
   seq=[self.p,[i*x for i,x in enumerate(self.p)][1:]]
   while seq[-1]!=[0]:
    rem=rpdiv(seq[-2],seq[-1])[1]
    if rem==[0]:break
    seq.append(rpscale(rem,-1))
   assert len(seq[-1])==1 and variations(seq,self.lo)-variations(seq,self.hi)==1
   assert rpeval(self.p,self.lo)*rpeval(self.p,self.hi)<0
  self.start_interval=[self.lo,self.hi];self.cache={self.zero:0}
 def f(self,x):return (F(x),)+(F(0),)*(self.d-1)
 def read(self,p):return self.reduce(list(reversed(list(map(F,p)))))
 def reduce(self,p):
  p=rpdiv(p,self.p)[1];return tuple(p+[F(0)]*(self.d-len(p)))
 def add(self,a,b):return tuple(x+y for x,y in zip(a,b,strict=True))
 def sub(self,a,b):return tuple(x-y for x,y in zip(a,b,strict=True))
 def scale(self,a,b):return tuple(x*b for x in a)
 def mul(self,a,b):return self.reduce(rpmul(a,b))
 def linear(self,vs,cs):
  out=self.zero
  for v,c in zip(vs,cs,strict=True):out=self.add(out,self.scale(v,c))
  return out
 def sign(self,a):
  if a in self.cache:return self.cache[a]
  for _ in range(1000):
   lo=hi=F(0)
   for c in reversed(a):
    xx=[lo*self.lo,lo*self.hi,hi*self.lo,hi*self.hi];lo=min(xx)+c;hi=max(xx)+c
   if lo>0:s=1;break
   if hi<0:s=-1;break
   assert self.d>1
   m=(self.lo+self.hi)/2;z=rpeval(self.p,m);assert z
   if z*rpeval(self.p,self.lo)>0:self.lo=m
   else:self.hi=m
  else:raise AssertionError(('unresolved sign',a,self.p,self.lo,self.hi))
  self.cache[a]=s;return s
 def ptrim(self,p):
  p=list(p)
  while len(p)>1 and p[-1]==self.zero:p.pop()
  return p
 def padd(self,a,b):return self.ptrim([self.add(a[i]if i<len(a)else self.zero,b[i]if i<len(b)else self.zero)for i in range(max(len(a),len(b)))])
 def pscale(self,a,c):return self.ptrim([self.scale(x,c)for x in a])
 def psub(self,a,b):return self.padd(a,self.pscale(b,-1))
 def pmul(self,a,b):
  z=[self.zero]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):z[i+j]=self.add(z[i+j],self.mul(x,y))
  return self.ptrim(z)
 def peval(self,p,t):
  z=self.zero
  for a in reversed(p):z=self.add(self.scale(z,t),a)
  return z
 def affine(self,p,l,r):
  q=[self.zero]
  for i,c in enumerate(p):
   for j in range(i+1):
    while len(q)<=j:q.append(self.zero)
    q[j]=self.add(q[j],self.scale(c,F(math.comb(i,j))*l**(i-j)*(r-l)**j))
  return self.ptrim(q)
 def positive(self,p):
  p=self.ptrim(p);assert p!=[self.zero];left=right=0
  while p[0]==self.zero:p=p[1:];left+=1
  while self.peval(p,F(1))==self.zero:
   q=[p[0]]
   for x in p[1:-1]:q.append(self.add(x,q[-1]))
   assert p[-1]==self.scale(q[-1],-1);p=self.ptrim(q);right+=1
  def bern(q,lo=F(0),hi=F(1),depth=0):
   n=len(q)-1;b=[self.linear(q[:k+1],[F(math.comb(k,i),math.comb(n,i))for i in range(k+1)])for k in range(n+1)]
   if all(self.sign(x)>0 for x in b):return [{'interval':[lo,hi],'Bernstein':b}]
   assert depth<12,('Bernstein failure',q,lo,hi,b)
   mid=(lo+hi)/2
   return bern(self.affine(q,F(0),F(1,2)),lo,mid,depth+1)+bern(self.affine(q,F(1,2),F(1)),mid,hi,depth+1)
  return {'zero_order_at_0':left,'zero_order_at_1':right,'reduced_polynomial':p,'positive_Bernstein_intervals':bern(p)}

def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def inv(A):
 n=len(A);M=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n)if M[i][j]);M[j],M[k]=M[k],M[j];v=M[j][j];M[j]=[x/v for x in M[j]]
  for i in range(n):
   if i!=j:
    z=M[i][j];M[i]=[a-z*b for a,b in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 M=[list(map(F,row))for row in A];out=[]
 for j in range(len(M)):
  q=M[j][j];assert q>0;out.append(q)
  for i in range(j+1,len(M)):
   for k in range(i,len(M)):M[k][i]=M[i][k]=M[i][k]-M[i][j]*M[k][j]/q
 return out

def run():
 packet=json.loads((ROOT/LOCAL).read_text());recipes=json.loads((ROOT/REC).read_text());walk=json.loads((ROOT/WALK).read_text())
 assert packet['status']=='PASS_EXACT_RANK_THREE_LOCAL_DIRECTION_CENSUS'and packet['direction_count']==91
 for doc in[packet,walk]:
  for p,h in doc['sources'].items():assert sha(p)==h
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[[dot(q,w)for q in Q]for w in W];G=[[dot(w,c)for c in C]for w in W]
 h=[sum(w[d]for w in W)for d in range(22)];gg=math.gcd(*h);h=[F(x,gg)for x in h]
 def pair(a,b):return dot(a,[dot(q,b)for q in Q])
 assert pair(h,h)<0 and all(pair(h,w)<0 for w in W)
 @lru_cache(None)
 def face(ids):
  A=[[G[i][j]for j in ids]for i in ids];piv=pd(A);I=inv(A);c=[dot(row,[dot(h,C[i])for i in ids])for row in I]
  x=[h[d]-sum(z*W[i][d]for z,i in zip(c,ids,strict=True))for d in range(22)];ps=[dot(x,c)for c in C]
  assert pair(x,x)<0 and pair(x,h)<0 and max(ps)==0 and [i for i,v in enumerate(ps)if not v]==list(ids)
  return {'walls':ids,'point':x,'norm':pair(x,x),'future':pair(x,h),'LDL':piv,'max_other_pairing':max(v for v in ps if v)}
 assert sorted(r['row']for r in recipes['families'])==list(range(91))
 proofs=[];pending=[];endpoint_proofs=[];counts=Counter()
 endpoint_recipes={r['row']:r for r in recipes['endpoints']};assert len(endpoint_recipes)==len(recipes['endpoints'])
 for rec in recipes['families']:
  ix=rec['row'];row=packet['directions'][ix];s=packet['systems'][row['system']];K=Field(row['minimal_polynomial'],row['root_interval']);Z=K.zero
  R=s['retained_walls'];factor=s['factor'];omit=s['omitted_wall'];i,j=rec['ridge']
  e=[K.read(v)for v in row['normal_constant']];D=list(map(F,row['null_shift']));L=K.read(row['parameter_lower_magnitude']);q=K.read(row['normal_square']);k=[K.read(v)for v in row['positive_kernel']]
  mixed=lambda e,x:K.linear(e,[dot(v,x)for v in Q])
  assert omit not in[i,j]and K.sign(L)>0 and K.sign(q)>0 and all(K.sign(v)>0 for v in k)
  assert all(e[d]==K.linear(k,[W[r][d]for r in R])for d in range(22))
  a=[K.linear(k,[G[r][t]for r in R])for t in R];qq=Z
  for aa,kk in zip(a,k,strict=True):qq=K.add(qq,K.mul(aa,kk))
  assert qq==q and mixed(e,D)==Z and pair(D,D)==0 and pair(D,h)<0
  fac=next(f for c in walk['corner_orbits']if c['id']==s['cusp']for f in c['affine_components']if f['walls']==factor)
  marks=dict(zip(factor,map(F,fac['positive_kernel']),strict=True));assert D==[sum(marks[r]*W[r][d]for r in factor)for d in range(22)]
  coeff=[K.sub(v,K.scale(L,marks[r]))for r,v in zip(R,k,strict=True)];assert all(K.sign(v)>=0 for v in coeff)
  assert any(v==Z for v in coeff)
  A0=tuple(sorted(set(R+[i,j])));left=[omit]if rec['direct_omitted_endpoint']else[r for r in factor if r!=row['lower_bound_witness_wall']]
  A1=tuple(sorted(set(left+[i,j])));F0=face(A0);F1=face(A1);common=sorted(set(A0)&set(A1))
  assert set(A1)-set(common)=={omit}and {i,j}<=set(common)
  pd([[G[a][b]for b in common]for a in common]);assert all(v==Z or r in common for r,v in zip(R,coeff,strict=True))
  p0=F0['point'];p1=F1['point'];d0=-pair(p0,D);d1=-pair(p1,D);assert d0>0 and d1>0;p0=[v/d0 for v in p0];p1=[v/d1 for v in p1]
  assert pair(p0,D)==pair(p1,D)==-1 and mixed(e,p0)==Z and mixed(e,p1)==K.scale(L,-1)and pair(p0,p1)<0
  aa=[mixed(e,W[r])for r in[i,j]];bb=[K.scale(L,-pair(D,W[r]))for r in[i,j]]
  pp=K.psub([K.scale(q,G[i][j])],K.pmul([aa[0],bb[0]],[aa[1],bb[1]]))
  ni=K.psub([K.scale(q,G[i][i])],K.pmul([aa[0],bb[0]],[aa[0],bb[0]]));nj=K.psub([K.scale(q,G[j][j])],K.pmul([aa[1],bb[1]],[aa[1],bb[1]]))
  NN=K.pmul(pp,pp);DD=K.pmul(ni,nj);tests={'first_norm':K.positive(ni),'second_norm':K.positive(nj)}
  if rec['mode']=='obtuse':tests['positive_pairing']=K.positive(pp);endpoint_good=K.sign(K.peval(pp,F(1)))>0
  else:
   lo,hi=map(F,rec['gap']);assert (lo,hi)in[(F(0),F(1,4)),(F(1,4),F(1,2))]
   tests['negative_pairing']=K.positive(K.pscale(pp,-1));low=K.psub(NN,K.pscale(DD,lo));high=K.psub(K.pscale(DD,hi),NN)
   tests['above_lower']=K.positive(low);tests['below_upper']=K.positive(high);endpoint_good=K.sign(K.peval(low,F(1)))>0 and K.sign(K.peval(high,F(1)))>0
  e1=[K.sub(v,K.scale(L,d))for v,d in zip(e,D,strict=True)]
  assert e1==[K.sub(K.linear(coeff,[W[r][d]for r in R]),K.scale(L,marks[omit]*W[omit][d]))for d in range(22)]
  boundary=all(v==Z for v in coeff)
  if not boundary:
   assert K.sign(K.peval(ni,F(1)))>0 and K.sign(K.peval(nj,F(1)))>0 and K.sign(K.peval(K.psub(DD,NN),F(1)))>0
  if not boundary and not endpoint_good:
   pending.append(ix);er=endpoint_recipes[ix];ii,jj=er['ridge'];assert omit not in[ii,jj]
   EA0=tuple(sorted(set(R+[ii,jj])));EA1=tuple(sorted(set([r for r in factor if r!=row['lower_bound_witness_wall']]+[ii,jj])))
   EF0=face(EA0);EF1=face(EA1);EK=sorted(set(EA0)&set(EA1));assert set(EA1)-set(EK)=={omit}and {ii,jj}<=set(EK)
   pd([[G[a][b]for b in EK]for a in EK]);assert all(v==Z or r in EK for r,v in zip(R,coeff,strict=True))
   ex0=EF0['point'];ex1=EF1['point'];ed0=-pair(ex0,D);ed1=-pair(ex1,D);assert ed0>0 and ed1>0;ex0=[v/ed0 for v in ex0];ex1=[v/ed1 for v in ex1]
   assert mixed(e,ex0)==Z and mixed(e,ex1)==K.scale(L,-1)and mixed(e1,ex1)==Z and mixed(e1,ex0)==L
   aaa=[mixed(e1,W[t])for t in[ii,jj]];ep=K.sub(K.scale(q,G[ii][jj]),K.mul(*aaa));eni=K.sub(K.scale(q,G[ii][ii]),K.mul(aaa[0],aaa[0]));enj=K.sub(K.scale(q,G[jj][jj]),K.mul(aaa[1],aaa[1]));en=K.mul(ep,ep);ed=K.mul(eni,enj)
   assert K.sign(eni)>0 and K.sign(enj)>0 and K.sign(K.sub(ed,en))>0
   if er['mode']=='obtuse':assert K.sign(ep)>0
   else:
    el,eh=map(F,er['gap']);assert(el,eh)in[(F(0),F(1,4)),(F(1,4),F(1,2))]
    assert K.sign(ep)<0 and K.sign(K.sub(en,K.scale(ed,el)))>0 and K.sign(K.sub(K.scale(ed,eh),en))>0
   endpoint_proofs.append({'row':ix,'ridge':[ii,jj],'mode':er['mode'],'gap':er['gap'],'face_at_0':EF0,'face_at_1':EF1,'normalized_point_0':ex0,'normalized_point_1':ex1,'common_zero_set':EK,
    'scaled_pairing':ep,'scaled_first_norm':eni,'scaled_second_norm':enj,'squared_cosine_numerator':en,'squared_cosine_denominator':ed,'redundancy_coefficients':coeff})
  counts[rec['mode']]+=1;counts[f'degree{K.d}']+=1
  proofs.append({'row':ix,'system':row['system'],'minimal_polynomial':row['minimal_polynomial'],'root_interval':row['root_interval'],'refined_root_interval':[K.lo,K.hi],
   'ridge':[i,j],'mode':rec['mode'],'gap':rec['gap'],'normal_square':q,'parameter_lower_magnitude':L,'face_at_0':F0,'face_at_1':F1,
   'normalized_point_0':p0,'normalized_point_1':p1,'point_cross_pair':pair(p0,p1),'common_zero_set':common,
   'scaled_pairing_polynomial':pp,'scaled_first_norm_polynomial':ni,'scaled_second_norm_polynomial':nj,'sign_certificates':tests,
   'u1_redundancy_coefficients':coeff,'u1_is_old_facet_hyperplane':boundary,'u1_bad_angle_verified':endpoint_good})
  print('row',ix,'PASS',rec['mode'],'endpoint',('separate certificate'if ix in endpoint_recipes else 'same ridge or old facet'),flush=True)
 assert sorted(pending)==sorted(endpoint_recipes)
 counts['separate_endpoint_certificates']=len(endpoint_proofs)
 counts['old_facet_endpoint_degeneracies']=sum(r['u1_is_old_facet_hyperplane']for r in proofs)
 counts['endpoint_obtuse']=sum(r['mode']=='obtuse'for r in endpoint_proofs)
 counts['endpoint_gap']=sum(r['mode']=='gap'for r in endpoint_proofs)
 assert counts['obtuse']==75 and counts['gap']==16 and len(endpoint_proofs)==22
 sources=[LOCAL,RAW,WALK,REC,f'scripts/{TAG}.py']
 return enc({'status':'PASS_ALL91_RANK_THREE_FULL_PARAMETER_GLOBAL_EXCLUSIONS','sources':{p:sha(p)for p in sources},'counts':dict(counts),'families':proofs,'separate_endpoint_certificates':endpoint_proofs,'pending_endpoint_rows':[],
  'field_encoding':'Coefficients in ascending powers of the isolated primitive root; minimal polynomial descending as in source.',
  'scope':'All91 archived local records excluded for the full admissible interval0<u<=1. Old-facet degeneracies miss the source interior. Local completeness remains a frozen dependency; no assertion about sections avoiding source cusps or changed chambers.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['counts']);print('pending endpoints',out['pending_endpoint_rows'])
