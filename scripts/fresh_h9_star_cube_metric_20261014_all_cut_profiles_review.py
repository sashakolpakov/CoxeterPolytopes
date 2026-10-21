#!/usr/bin/env python3
"""Separate square-root-basis arithmetic for all eleven cut profiles.

Reuses only the independent Fraction field engine from the preceding review;
no producer/SymPy imports. Fresh domain guards are rebuilt with unused cuts
right, since (2,5) is not a monotone increase of (3,3).
"""
import argparse,json
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import isqrt
from pathlib import Path
from fresh_h9_star_cube_metric_20261014_cut_profiles_review import F,Z,O,CO,add,sub,neg,mul,scale,total,inverse,determinant,anp,enc
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014_all_cut_profiles'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
@lru_cache(None)
def sign(x):
 if x==Z:return 0
 digits=20
 while True:
  den=10**digits;lo=hi=x[0]
  for a,n in zip(x[1:],[2,5,10]):
   p=isqrt(n*den*den);assert p*p<n*den*den<(p+1)*(p+1)
   l,h=F(p,den),F(p+1,den)
   if a>=0:lo+=a*l;hi+=a*h
   else:lo+=a*h;hi+=a*l
  if lo>0:return 1
  if hi<0:return -1
  digits*=2
def det_gauss(A):
 a=[list(r)for r in A];d=O;n=len(a)
 for j in range(n):
  k=next((k for k in range(j,n)if a[k][j]!=Z),None)
  if k is None:return Z
  if k!=j:a[j],a[k]=a[k],a[j];d=neg(d)
  p=a[j][j];d=mul(d,p);ip=inverse([[p]])[0][0]
  for k in range(j+1,n):
   q=mul(a[k][j],ip)
   if q!=Z:
    for l in range(j+1,n):a[k][l]=sub(a[k][l],mul(q,a[j][l]))
   a[k][j]=Z
 return d
def run():
 src=f'results/{TAG}.json';data=json.loads((ROOT/src).read_text())
 assert data['status']=='PASS_ALL_CUT_CONTACT_PROFILES_EXCLUDED'
 for p,h in data['sources'].items():assert sha(p)==h,p
 gp='results/fresh_h9_tree_cube_wedges_20261014.json';geometry=json.loads((ROOT/gp).read_text())
 for p,h in geometry['sources'].items():assert sha(p)==h,p
 g=geometry['representatives'][0]['geometry'];V=[v['incident']for v in g['vertices']];M={tuple(m)for m in g['minimal_nonfaces']}
 S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
 assert all(tuple(sorted((v,*S[2*i:2*i+2])))in M for i in range(4)for v in[P[i],C[i]])
 base={}
 for i in range(4):
  base[S[2*i],S[2*i+1]]=4
  for j in range(2):base[tuple(sorted((P[i],S[2*i+j])))]=3
 def gram(nodes,labels):return [[O if a==b else neg(CO[labels.get(tuple(sorted((a,b))),2)])for b in nodes]for a in nodes]
 def check_failure(row,labels):
  v=V[row['actual_vertex']];pair=tuple(row['pair']);assert set(pair)<=set(v)
  lab=labels|{pair:row['test_label']};A=gram(v,lab);k=row['leading_order']
  d=det_gauss([r[:k]for r in A[:k]]);p=det_gauss([r[:k-1]for r in A[:k-1]])if k>1 else O
  assert sign(p)>0 and sign(d)<=0 and mul(anp(row['last_pivot']),p)==d
  return dict(pair=list(pair),actual_vertex=row['actual_vertex'],prefix=k,determinant=enc(d))
 failures=[check_failure(r,base)for r in data['forced_cross_right_witnesses']]
 cross={tuple(sorted((a,b)))for i in range(4)for j in range(i+1,4)for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]}
 assert {tuple(r['pair'])for r in failures}==cross and len(failures)==24
 profiles=[p for p in product(range(2,6),repeat=2)if F(1,4)+sum(F(1,m)for m in p)<1]
 assert list(map(list,profiles))==data['profiles']and len(profiles)==11
 own=[]
 for row in data['own_pair_right_witnesses']:
  i=row['group'];profile=tuple(row['profile']);assert profile in profiles
  j=next(j for j,m in enumerate(profile)if m>=3);tri=sorted([P[i],C[i],S[2*i+j]])
  assert tri==row['three_support']and set(tri)<=set(V[row['actual_vertex']])
  labels=base|{tuple(sorted((P[i],C[i]))):3,tuple(sorted((C[i],S[2*i+j]))):profile[j]}
  A=gram(tri,labels);q=total(x for r in A for x in r)
  assert q==anp(row['ones_quadratic'])and sign(q)<=0
  own.append((i,profile))
 assert set(own)==set(product(range(4),profiles))and len(own)==44
 for row in data['own_contact_guards']:
  i=C.index(row['cut']);root=row['increased_own_root'];assert root in S[2*i:2*i+2]
  other=next(x for x in S[2*i:2*i+2]if x!=root)
  assert set(row['eight_support'])==(set(S)-{other})|{C[i]}
  assert set(row['eight_support'])<=set(V[row['actual_vertex']])
  assert set(row['label6_path'])=={P[i],C[i],root}<=set(V[row['path_actual_vertex']])
 assert len(data['own_contact_guards'])==8
 # The preceding symbolic argument does not assume m<=5: first PiCi is
 # forced right by any own m>=3, then the 3-m path rules out all m>=6.
 assert 1-F(1,4)-F(3,4)==0
 B=gram(S,base);J=inverse(B)
 assert all(total(mul(B[i][k],J[k][j])for k in range(8))==(O if i==j else Z)for i in range(8)for j in range(8))
 for leaf in[P[3]]+C:
  i=P.index(leaf)if leaf in P else C.index(leaf)
  for j in range(2):
   nodes=[x for x in S if x!=S[2*i+1-j]];support=set(nodes)|{leaf}
   assert any(support<=set(v)for v in V)
   K=inverse(gram(nodes,base));ownpos=nodes.index(S[2*i+j])
   assert K[ownpos][ownpos]==O and all(sign(x)>=0 for r in K for x in r)
   foreign=[k for k,x in enumerate(nodes)if x not in S[2*i:2*i+2]]
   assert len(foreign)==6 and all(K[k][k]==scale(O,2)for k in foreign)
 assert F(1,4)+2*F(1,4)*2>1 and F(1,4)+2*F(1,2)>1
 assert F(1,2)+2*F(1,4)==1 and sign(sub(add(mul(CO[5],CO[5]),CO[3]),O))>0
 upper=0
 for row in data['anchor_label_upper_witnesses']:
  left,right=tuple(row['C0_profile']),tuple(row['C3_profile']);assert left in profiles and right in profiles
  labels=base.copy()
  for i,profile in[(0,left),(3,right)]:
   for j,m in enumerate(profile):labels[tuple(sorted((C[i],S[2*i+j])))]=m
  assert {tuple(r['pair'])for r in row['witnesses']}=={(11,14),(14,17)}
  for r in row['witnesses']:assert r['test_label']==4;failures.append(check_failure(r,labels));upper+=1
 assert upper==242 and len(data['anchor_label_upper_witnesses'])==121
 @lru_cache(None)
 def pdot(x,y):return total(mul(x[i],mul(J[i][j],y[j]))for i in range(8)for j in range(8)if x[i]!=Z and J[i][j]!=Z and y[j]!=Z)
 def patterns(group,profile):
  ans=[]
  for extra in([None]+[j for j in range(8)if j//2!=group]if profile==(3,3)else[None]):
   q=[Z]*8;q[2*group]=neg(CO[profile[0]]);q[2*group+1]=neg(CO[profile[1]])
   if extra is not None:q[extra]=neg(CO[3])
   ans.append((extra,tuple(q)))
  return ans
 saved={}
 for row in data['records']:
  for case in row['cases']:saved[tuple(row['C0_profile']),tuple(row['C3_profile']),tuple(case['foreign_contacts']),tuple(case['anchor_labels'])]=anp(case['residual_determinant'])
 assert len(saved)==8092
 n=0;zero=[];anchors=set()
 for left,right in product(profiles,repeat=2):
  for ia,x in patterns(3,(3,3)):
   for ib,y in patterns(3,right):
    qx=sub(O,pdot(x,x));qy=sub(O,pdot(y,y));xy=neg(pdot(x,y))
    assert determinant([[qx,xy],[xy,qy]])!=Z;anchors.add((right,ia,ib))
    for ic,z in patterns(0,left):
     qz=sub(O,pdot(z,z))
     for ma,mb in product([2,3],repeat=2):
      xz=sub(neg(CO[ma]),pdot(x,z));yz=sub(neg(CO[mb]),pdot(y,z))
      d=determinant([[qx,xy,xz],[xy,qy,yz],[xz,yz,qz]])
      assert d==saved[left,right,(ia,ib,ic),(ma,mb)];n+=1
      if d==Z:zero.append((left,right,(ia,ib,ic),(ma,mb)))
 assert n==8092 and len(zero)==4 and len(anchors)==119
 closed=[]
 for row in data['actual_vertex_exclusions']:
  z=row['rank_survivor'];key=(tuple(z['C0_profile']),tuple(z['C3_profile']),tuple(z['foreign_contacts']),tuple(z['anchor_labels']))
  assert key in zero
  left,right,(ia,ib,ic),(ma,mb)=key
  assert left==(3,3)and right==(4,4)and ia in[0,1]and ib is None and ic in[6,7]and(ma,mb)==(2,3)
  cycle=[14,S[ia],11,S[ic]];vi=next(i for i,v in enumerate(V)if set(cycle)<=set(v))
  A=[[O if i==j else neg(CO[3])if(i-j)%4 in[1,3]else Z for j in range(4)]for i in range(4)]
  assert all(total(r)==Z for r in A)
  labels=base.copy()
  for i,pr in[(0,left),(3,right)]:
   for j,m in enumerate(pr):labels[tuple(sorted((C[i],S[2*i+j])))]=m
  for leaf,extra in zip([P[3],C[3],C[0]],[ia,ib,ic]):
   if extra is not None:labels[tuple(sorted((leaf,S[extra])))]=3
  for leaf,m in zip([P[3],C[3]],[ma,mb]):labels[tuple(sorted((C[0],leaf)))]=m
  v=V[row['actual_vertex']];k=row['leading_order'];G=gram(v,labels)
  d=det_gauss([r[:k]for r in G[:k]]);p=det_gauss([r[:k-1]for r in G[:k-1]])if k>1 else O
  assert sign(p)>0 and sign(d)<=0 and mul(anp(row['last_pivot']),p)==d
  closed.append(dict(C0_profile=list(left),C3_profile=list(right),foreign_contacts=list(key[2]),cycle=cycle,actual_vertex=vi,positive_null_vector=[1]*4))
 assert len(closed)==4 and not data['surviving_cases']
 paths=[src,f'scripts/{TAG}.py',gp,'scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py',f'scripts/{TAG}_review.py']
 return dict(status='PASS_SEPARATE_FULL_CUT_PROFILE_ARITHMETIC_AND_PROOF_REVIEW',sources={p:sha(p)for p in paths},
  ordered_cut_profiles=[list(p)for p in profiles],fresh_cross_block_witnesses=24,fresh_own_pair_witnesses=44,
  fresh_anchor_upper_witnesses=upper,exact_failure_minors=failures,
  profile_pairs=121,residual_cases=n,rank_nonzero_cases=n-len(zero),distinct_nondegenerate_anchor_cases=len(anchors),
  four_remaining_cases_closed=closed,
  scope='No producer or SymPy arithmetic. Reuses the separate Fraction biquadratic engine, but rebuilds all fresh24cross,44own-pair,242anchor-upper and8path/support guards. Unused cuts are relaxed to right, so no invalid monotone transfer through(2,5) occurs. All8092determinants are independently reconstructed and every primary value compared. All4rank-zero cases have explicit actual affineA3 null minors. Excludes compact realizations only with fixed fourB2blocks and positive-leaf3/3; all cut owncontacts are covered. Changed B2/positive contacts and the full star incidence remain open.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['residual_cases'],len(d['four_remaining_cases_closed']),flush=True)
