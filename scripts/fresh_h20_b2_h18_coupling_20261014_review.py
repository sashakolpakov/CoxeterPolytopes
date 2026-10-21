#!/usr/bin/env python3
"""Independent Fraction Q(sqrt2) replay of direct B2 profiles and escape."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_b2_h18_coupling_20261014_review'
PRIMARY='results/fresh_h20_b2_h18_coupling_20261014.json'
PRODUCER='scripts/fresh_h20_b2_h18_coupling_20261014.py'
SOURCE='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
class K:
 def __init__(self,a=0,b=0):
  if isinstance(a,K):self.a,self.b=a.a,a.b
  else:self.a,self.b=F(a),F(b)
 def __add__(self,y):y=K(y);return K(self.a+y.a,self.b+y.b)
 __radd__=__add__
 def __neg__(self):return K(-self.a,-self.b)
 def __sub__(self,y):return self+-K(y)
 def __rsub__(self,y):return K(y)+-self
 def __mul__(self,y):y=K(y);return K(self.a*y.a+2*self.b*y.b,self.a*y.b+self.b*y.a)
 __rmul__=__mul__
 def __truediv__(self,y):
  y=K(y);q=y.a*y.a-2*y.b*y.b;assert q;return self*K(y.a/q,-y.b/q)
 def __rtruediv__(self,y):return K(y)/self
 def __eq__(self,y):y=K(y);return(self.a,self.b)==(y.a,y.b)
 def __bool__(self):return bool(self.a or self.b)
 def sign(self):
  a,b=self.a,self.b
  if not b:return (a>0)-(a<0)
  if not a:return (b>0)-(b<0)
  if(a>0)==(b>0):return(a>0)-(a<0)
  return ((a>0)-(a<0))*((a*a>2*b*b)-(a*a<2*b*b))
 def enc(self):return dict(a=str(self.a),b=str(self.b))
def dec(x):return K(x['a'],x['b'])
def enc(x):
 if isinstance(x,K):return x.enc()
 if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def dot(x,y):return sum((a*b for a,b in zip(x,y)),K())
def mv(M,v):return[dot(r,v)for r in M]
def plus(x,y):return[a+b for a,b in zip(x,y)]
def scale(x,c):return[a*c for a in x]
def inverse(M):
 n=len(M);a=[r[:]+[K(i==j)for j in range(n)]for i,r in enumerate(M)]
 for k in range(n):
  j=next(j for j in range(k,n)if a[j][k]);a[k],a[j]=a[j],a[k];q=a[k][k];a[k]=[x/q for x in a[k]]
  for j in range(n):
   if j!=k:
    q=a[j][k];a[j]=[x-q*y for x,y in zip(a[j],a[k])]
 return[r[n:]for r in a]
def inertia(M):
 a=[r[:]for r in M];p=n=z=0
 while a:
  k=next((i for i in range(len(a))if a[i][i]),None)
  if k is not None:
   ids=[k]+[i for i in range(len(a))if i!=k];a=[[a[i][j]for j in ids]for i in ids];q=a[0][0]
   if q.sign()>0:p+=1
   else:n+=1
   a=[[a[i][j]-a[i][0]*a[0][j]/q for j in range(1,len(a))]for i in range(1,len(a))];continue
  pair=next(((i,j)for i in range(len(a))for j in range(i)if a[i][j]),None)
  if pair is None:z+=len(a);break
  i,j=pair;ids=[i,j]+[k for k in range(len(a))if k not in pair];a=[[a[i][j]for j in ids]for i in ids];q=a[0][1];p+=1;n+=1
  a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/q for j in range(2,len(a))]for i in range(2,len(a))]
 return p,n,z
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/PRIMARY).read_text());raw=json.loads((ROOT/SOURCE).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 G=[[dec(x)for x in r]for r in raw['gram']];basis=d['source_basis'];assert len(set(basis))==19
 Q0=[[G[i][j]for j in basis]for i in basis];J=inverse(Q0)
 old=[mv(J,[G[j][i]for j in basis])for i in range(40)]
 oldcov=[mv(Q0,r)for r in old]
 assert all(dot(old[i],oldcov[j])==G[i][j]for i in range(40)for j in range(40))
 rt=K(0,1);Q=[r+[K(),K()]for r in Q0]+[[K()]*19+[K(2),-rt],[K()]*19+[-rt,K(2)]]
 assert Q==[[dec(x)for x in r]for r in d['form']]and inertia(Q)==(20,1,0)
 support=[[K()]*19+[K(1),K()],[K()]*19+[K(),K(1)]]
 weights=[[K()]*19+[K(1),rt/2],[K()]*19+[rt/2,K(1)]]
 assert all(dot(weights[i],mv(Q,support[j]))==K(i==j)for i in range(2)for j in range(2))
 assert (2+rt-2).sign()>0 and(-rt/2).sign()<0 and(1-rt/2).sign()>0 and(rt/2-F(1,2)).sign()>0
 forced=set()
 for i in range(40):
  for j in range(i):
   if G[i][j].sign()<0 and 4*G[i][j]*G[i][j]==G[i][i]*G[j][j]:forced|={i,j}
 eligible=sorted(set(range(40))-forced)
 assert len(forced)==24 and len(eligible)==16 and eligible==d['eligible_special']
 assert {r['source_wall']for r in d['forced_ordinary']}==forced
 for r in d['forced_ordinary']:
  i,j=r['source_wall'],r['source_neighbor'];assert i in forced and 4*G[i][j]*G[i][j]==G[i][i]*G[j][j]and G[i][j].sign()<0
 profiles=[]
 for i,r in enumerate(old):
  assert G[i][i]in(K(1),K(2))
  ordinary=scale(r,rt if G[i][i]==1 else K(1))+[K(),K()];rows=[ordinary]
  if i in eligible:
   y=scale(r,K(1)if G[i][i]==1 else rt/2)+[K(),K()]
   rows += [plus(y,scale(w,-1))for w in weights]
  profiles.append(rows)
 assert profiles==[[list(map(dec,r))for r in p]for p in d['inducing_profiles']]
 bank=support+[r for p in profiles for r in p];cov=[mv(Q,r)for r in bank]
 assert len(bank)==74 and all(dot(r,c)==2 for r,c in zip(bank,cov))
 assert bank==[list(map(dec,r))for r in d['bank_normals']]
 x=list(map(dec,d['source_point']));h=[sum((r[j]for r in old),K())for j in range(19)]
 assert dot(x,mv(Q0,x))==K(F(3857,529),F(-2800,529))and dot(x,mv(Q0,x)).sign()<0
 xp=[dot(x,c)for c in oldcov];hp=[dot(h,c)for c in oldcov]
 assert all(v.sign()<=0 for v in xp)and all(v.sign()<0 for v in hp)
 assert all((xp[i]+(K(1)if G[i][i]==1 else rt)).sign()<=0 for i in eligible)
 assert dot(h,mv(Q0,h))==K(-2952,-1280)and dot(x,mv(Q0,h))==K(-16,-40)
 w=scale(plus(weights[0],weights[1]),F(-1,2));C=plus(h+[K(),K()],w);E=plus(plus(x,scale(h,F(1,1000)))+[K(),K()],w)
 assert C==list(map(dec,d['strict_center']))and E==list(map(dec,d['strict_escape']))
 cp=[dot(C,c)for c in cov];ep=[dot(E,c)for c in cov]
 assert all(v.sign()<0 for v in cp+ep)
 qC=dot(C,mv(Q,C));qE=dot(E,mv(Q,E));EC=dot(E,mv(Q,C))
 assert qC.sign()<0 and EC.sign()<0 and qE.sign()>0
 assert qE==K(F(512876299,66125000),F(-33884339,6612500))
 chosen=[1 if G[i][i]==1 else 0 for i in range(40)];N=support+[profiles[i][j]for i,j in enumerate(chosen)];NC=[mv(Q,r)for r in N]
 counts=Counter()
 for i in range(42):
  for j in range(i):
   g=dot(N[i],NC[j]);assert g.sign()<=0
   if g==0:kind='right'
   elif g==-1:kind='pi3'
   elif g==-rt:kind='pi4'
   elif g==-2:kind='parallel'
   else:assert(g+2).sign()<0;kind='dotted'
   counts[kind]+=1
 assert dict(counts)==d['coherent_example']['pair_counts']
 for i,r in enumerate(N):
  p=plus(C,scale(r,-dot(C,NC[i])/2))
  assert p==list(map(dec,d['coherent_example']['actual_facet_points'][i]))
  assert dot(p,mv(Q,p)).sign()<0
  assert all(dot(p,NC[j]).sign()==(0 if i==j else -1)for j in range(42))
 sourcepoint=h+[K(),K()]
 assert all(dot(sourcepoint,NC[i]).sign()==(0 if i<2 else -1)for i in range(42))
 names=raw['labels'];i,j=map(names.index,['L1','R10']);cycle=G[i][j]*G[i][j];assert cycle==K(12,8)
 i,j=map(names.index,['R23','R24']);assert G[i][j]==-6 and G[i][i]*G[j][j]-G[i][j]*G[i][j]==-32
 return enc(dict(status='PASS_INDEPENDENT_DIRECT_B2_COUPLING_REVIEW',sources={p:sha(p)for p in [f'scripts/{TAG}.py',PRIMARY,PRODUCER,SOURCE]},
  arithmetic='Standalone Fraction arithmetic in 1,sqrt2; no producer/project imports.',
  source_gram_entries=1600,ambient_inertia=[20,1,0],forced_ordinary=sorted(forced),eligible_special=eligible,
  exact_profile_norms=74,all_bank_center_pairings=cp,all_bank_escape_pairings=ep,
  coherent_example_pair_counts=dict(counts),coherent_example_pairs=861,actual_facets=42,facet_signs=1764,
  escape_norm=qE,escape_center_pairing=EC,conditional_nonquasi_cycle=cycle,
  proof_review='Positive B2 inverse gives exactly ordinary or single-pi3 inducing profiles. Every source pi3 edge forces both ends ordinary. The strict escape for the full74-normal union excludes every direct coherent assignment. Additional noninducing walls and changed sources are excluded from this theorem.',
  dependencies='The source finite-volume hybrid is inherited; this review reconstructs its literal40-root Gram but does not repeat its original gluing or complete vertex census.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS independent B2 saved replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print('PASS independent B2 generation')
