#!/usr/bin/env python3
"""Narrow source review with independent quadratic Schur and coverage checks."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_disjoint_walk_20261014_endpoint_scope_review'
class Q:
 def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
 def __add__(self,o):
  if not isinstance(o,Q):o=Q(o)
  return Q(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return Q(-self.a,-self.b)
 def __sub__(self,o):return self+-o if isinstance(o,Q)else self+Q(-o)
 def __rsub__(self,o):return Q(o)+-self
 def __mul__(self,o):
  if not isinstance(o,Q):o=Q(o)
  return Q(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
 __rmul__=__mul__
 def __truediv__(self,o):
  if not isinstance(o,Q):o=Q(o)
  d=o.a*o.a-5*o.b*o.b;assert d
  return self*Q(o.a/d,-o.b/d)
 def __bool__(self):return bool(self.a or self.b)
 def __eq__(self,o):
  if not isinstance(o,Q):o=Q(o)
  return self.a==o.a and self.b==o.b
 def sign(self):
  a,b=self.a,self.b
  if not b:return (a>0)-(a<0)
  if not a:return (b>0)-(b<0)
  if a>0 and b>0:return 1
  if a<0 and b<0:return -1
  d=a*a-5*b*b;assert d
  return ((d>0)-(d<0))*(1 if a>0 else -1)
 def enc(self):return [str(self.a),str(self.b)]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def pd(M):
 A=[row[:]for row in M]
 for k in range(len(A)):
  assert A[k][k].sign()>0
  for i in range(k+1,len(A)):
   for j in range(k+1,len(A)):A[i][j]=A[i][j]-A[i][k]*A[k][j]/A[k][k]
def tree_counts(row):
 c=row['counts'];leaves=sum(n for k,n in c.items()if k!='boxes')
 assert c['boxes']==2*leaves-1 and row['h_bound']>=1
 assert len(row['tree_sha256'])==64
 return c
def run():
 seen={}
 def load(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for f,h in d.get('sources',{}).items():
   assert sha(f)==h,f
   if f.endswith('.json'):load(f)
  return d
 pre='results/fresh_h9_h8_disjoint_walk_20261014'
 f=load(pre+'.json');u=load(pre+'_universal_gate.json')
 aff=load(pre+'_affine.json');gate=load(pre+'_affine_gate.json')
 assert f['open_count']==62 and len(f['queue'])==len(u['edges'])==62
 finite=tail=0;total=Counter();matrices={}
 for germ,record in zip(f['queue'],u['edges']):
  i=germ['id'];assert i==record['edge_id']
  z=germ['profile_data'];S=germ['support'];assert len(S)==len(set(S))==8 and z['support']==S
  M=[[Q(*v)for v in row]for row in z['normalized_Gram']]
  J=[[Q(*v)for v in row]for row in z['normalized_inverse']]
  pd(M);assert all(J[a][b].sign()>=0 for a in range(8)for b in range(8))
  for a in range(8):
   for b in range(8):assert sum((M[a][k]*J[k][b]for k in range(8)),Q())==Q(a==b)
  matrices[i]=(S,M,J);wanted=[]
  for p in z['finite_profiles']:
   ids=[S.index(w)for w,m in p['contacts']]
   assert len(ids)==len(set(ids))and all(not J[a][b]for a in ids for b in ids if a!=b)
   cost=sum((J[S.index(w)][S.index(w)]*{3:Q(1),4:Q(2),5:Q(F(3,2),F(1,2))}[m]for w,m in p['contacts']),Q())
   assert cost==Q(*p['projection_square']) and 2-cost==Q(*p['residual_norm'])and(2-cost).sign()>0
   wanted.append(('finite',tuple(map(tuple,p['contacts']))));finite+=1
  for p in z['dihedral_tail_profiles']:
   j=S.index(p['wall']);assert J[j][j]==Q(F(1,2))and all(not J[k][j]for k in range(8)if k!=j)
   wanted.append(('tail',p['wall']));tail+=1
  got=[]
  for p in record['profiles']:
   got.append((p['kind'],tuple(map(tuple,p['contacts']))if p['kind']=='finite'else p['wall']))
   total.update(tree_counts(p))
  assert got==wanted
 assert (finite,tail)==(1076,64)
 assert total==Counter({k:v for k,v in u['counts'].items()if k not in('finite_profiles','tail_profiles')})and total['boxes']==4964
 expected={};kernels=[]
 for row in aff['supports']:
  S,M,J=matrices[row['germ']];assert row['support']==S
  for p in row['affine_profiles']:
   contacts=p['contacts'];aa=[];bb=[]
   for i in range(8):
    a=sum((J[i][S.index(w)]for w,m in contacts if m==3),Q())
    b=sum((J[i][S.index(w)]for w,m in contacts if m==4),Q())
    assert a.sign()>=0 and b.sign()>=0 and(a or b)
    assert [a.enc(),b.enc()]==p['positive_kernel_first8_as_a_plus_b_sqrt2'][i]
    aa.append(a);bb.append(b)
   cost0=sum((aa[S.index(w)]if m==3 else 2*bb[S.index(w)]for w,m in contacts),Q())
   cost1=sum((bb[S.index(w)]if m==3 else aa[S.index(w)]for w,m in contacts),Q())
   assert cost0==2 and not cost1
   for i in range(8):
    assert sum((M[i][j]*aa[j]for j in range(8)),Q())==Q([S[i],3]in contacts)
    assert sum((M[i][j]*bb[j]for j in range(8)),Q())==Q([S[i],4]in contacts)
   key=(row['germ'],tuple(map(tuple,contacts)));assert key not in expected;expected[key]=1
   kernels.append(dict(germ=row['germ'],contacts=contacts,cost=2,all_kernel_weights_strictly_positive=True))
 atotal=Counter();actual=set()
 for r in gate['records']:
  key=(r['germ'],tuple(map(tuple,r['contacts'])));assert key not in actual;actual.add(key)
  atotal.update(tree_counts(r))
  # The stored dyadic endpoints use one positive common scale; integer
  # comparisons independently verify the monotonicity inequalities.
  aa=[tuple(map(int,v))for v in r['a']];cc=[tuple(map(int,v))for v in r['c']]
  assert len(aa)==len(cc)==len(r['outside'])+1
  assert all(a[1]<=0 and a[1]+c[1]<0 for a,c in zip(aa,cc))
  assert aa[-1][1]<0<cc[-1][0]
 assert actual==set(expected)and len(actual)==14
 assert atotal==Counter({k:v for k,v in gate['counts'].items()if k!='profiles'})and atotal['boxes']==50
 paths=[pre+s for s in('.json','_universal_gate.json','_affine.json','_affine_gate.json')]
 paths +=['scripts/fresh_h9_h8_disjoint_walk_20261014'+s for s in('_universal_gate.py','_affine.py','_affine_gate.py')]
 paths +=[f'scripts/{TAG}.py']
 return dict(status='PASS_ENDPOINT_PARAMETER_SOURCE_AND_SCHUR_REVIEW',sources={p:sha(p)for p in paths},
  recursive_manifests=len(seen),positive_definite_supports=62,finite_profiles=finite,continuous_tails=tail,
  affine_positive_kernels=kernels,finite_gate_boxes=4964,affine_gate_boxes=50,
  scope='Independent Fraction Qsqrt5 inverse/PD/Schur/kernel and literal coverage checks, plus source/proof review. Physical coefficient reconstruction and interval-tree arithmetic are source-reviewed and inherited, not independently reimplemented or replayed here. Complete local profile classification is inherited. Original62 supports only; no later71-support review or arbitrary-completion theorem.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('SAVED_REPLAY_PASS')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['finite_profiles'],d['continuous_tails'],len(d['affine_positive_kernels']))
