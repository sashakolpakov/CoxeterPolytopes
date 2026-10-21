#!/usr/bin/env python3
"""Complete finite/PSD endpoint profiles for affine A6 x A6 x A7."""
import argparse,json,hashlib
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from math import isqrt,sqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_balanced_endpoint_20261014'
RAD=[1,2,3,6,5,10,15,30];ZERO=(F(0),)*8;SCALE=1<<120
def c(x):return (F(x),)+(F(0),)*7
def add(a,b):return tuple(x+y for x,y in zip(a,b,strict=True))
def scale(a,x):return tuple(v*x for v in a)
def mul(a,b):
 out=[F(0)]*8
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if x and y:
    k=i&j;out[i^j]+=x*y*((2 if k&1 else 1)*(3 if k&2 else 1)*(5 if k&4 else 1))
 return tuple(out)
def sign(a):
 if a==ZERO:return 0
 lo=hi=F(0)
 for n,x in zip(RAD,a,strict=True):
  k=isqrt(n*SCALE*SCALE);l=F(k,SCALE);h=l if k*k==n*SCALE*SCALE else F(k+1,SCALE)
  lo+=x*(l if x>=0 else h);hi+=x*(h if x>=0 else l)
 assert lo>0 or hi<0
 return 1 if lo>0 else -1
def approx(a):return sum(float(x)*sqrt(n)for x,n in zip(a,RAD,strict=True))
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(tuple,list)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def run():
 factors=[list(range(7)),list(range(7,14)),list(range(14,22))];ranks=[6,6,7]
 G=[[2*int(i==j)for j in range(22)]for i in range(22)]
 for f in factors:
  for i,w in enumerate(f):v=f[(i+1)%len(f)];G[w][v]=G[v][w]=-1
 omitted=[f[0]for f in factors];ids=[i for i in range(22)if i not in omitted]
 B=[[G[i][j]for j in ids]for i in ids];I=[[F(0)for _ in ids]for _ in ids];off=0
 for n in ranks:
  for i in range(n):
   for j in range(n):I[off+i][off+j]=F(min(i+1,j+1)*(n-max(i,j)),n+1)
  off+=n
 assert all(sum(B[i][k]*I[k][j]for k in range(19))==int(i==j)for i in range(19)for j in range(19))
 assert all(x>=0 for row in I for x in row)and min(I[i][i]for i in range(19))==F(6,7)
 assert F(18,7)>2
 phi2=(F(3,2),F(0),F(0),F(0),F(1,2),F(0),F(0),F(0))
 assert sign(add(scale(phi2,F(6,7)),c(-2)))>0
 mags={3:c(1),4:(F(0),F(1),F(0),F(0),F(0),F(0),F(0),F(0))}
 profiles=[];trials=Counter();counts=Counter()
 for k in range(3):
  for inds in combinations(range(19),k):
   for labels in product((3,4),repeat=k):
    trials[k]+=1;cc=[ZERO]*19
    for i,m in zip(inds,labels):cc[i]=mags[m]
    cost=ZERO
    for i in inds:
     for j in inds:cost=add(cost,scale(mul(cc[i],cc[j]),I[i][j]))
    residual=add(c(2),scale(cost,-1));sg=sign(residual)
    if sg<0:continue
    cos=[scale(x,F(1,2))for x in cc];ratios=[]
    for f,n,o in zip(factors,ranks,omitted):
     b=ZERO
     for i,w in enumerate(ids):
      if w in f:b=add(b,cos[i])
     ratios.append(dict(factor_index=len(ratios),omitted=o,numerator_over_radius=n+1,
                        numerator_scale='factor_inradius',subtract=b,subtract_approx=approx(b)))
    kind='finite'if sg>0 else'ideal_PSD';counts[kind]+=1
    profiles.append(dict(id=len(profiles),contacts=[[ids[i],m]for i,m in zip(inds,labels)],
       cosine_vector=cos,cosine_vector_approx=list(map(approx,cos)),
       norm_two_projection_square=cost,norm_two_schur_residual=residual,
       normalized_spatial_cost=scale(cost,F(1,2)),normalized_spatial_cost_approx=approx(cost)/2,
       endpoint=kind,omitted_contact_ratios=ratios))
 assert dict(counts)=={'finite':37,'ideal_PSD':4}and sum(trials.values())==723
 assert [p['contacts']for p in profiles if p['endpoint']=='ideal_PSD']==[
  [[18,3]],[[1,3],[6,3]],[[8,3],[13,3]],[[15,3],[21,3]]]
 transports=[]
 for om in product(*factors):
  perm={w:f[(i+f.index(o))%len(f)]for f,o in zip(factors,om)for i,w in enumerate(f)}
  assert [perm[w]for w in omitted]==list(om)
  assert all(G[perm[i]][perm[j]]==G[i][j]for i in range(22)for j in range(22))
  mapped=[perm[w]for w in ids]
  assert [[G[i][j]for j in mapped]for i in mapped]==B
  transports.append(dict(omitted=list(om),edge_support_in_profile_order=mapped,wall_permutation=[perm[i]for i in range(22)]))
 assert len(transports)==392
 helper='scripts/fresh_h20_balanced_20261014.py'
 paths=[str(Path(__file__).resolve().relative_to(ROOT)),helper]
 typ=dict(id=0,omitted=omitted,edge_support=ids,literal_edge_multiplicity=392,
          PD19_Gram=B,inverse_Gram=I,counts=dict(counts),profiles=profiles,dihedral_tails=[])
 return enc(dict(status='PASS_COMPLETE_BALANCED_INTRINSIC_ENDPOINT_PROFILE_CENSUS',
  sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  field_basis=RAD,seed_order=list(range(22)),normalized_seed_Gram=G,factors=factors,
  factor_ranks=ranks,positive_affine_marks=[1]*22,types=[typ],literal_edge_transports=transports,
  finite_pool_trials=dict(trials),total_representative_profiles=41,
  literal_profile_counts=dict(finite=392*37,ideal_PSD=392*4,total=392*41),
  completeness=dict(inverse_diagonal_lower_bound='6/7',support_at_most=2,labels=[3,4],
   excluded_all_labels_at_least5=True,parallel_contact_excluded=True,no_dihedral_tail=True),
  scope='Every cap giving a finite or PSD ideal endpoint on any of the392 original cusp edges, retaining all22 seed walls, has one of these profiles after positive norm-two normalization. Arbitrary real cap coordinates allowed. Omitted-wall angles remain unrestricted and require the three saved radius equations. PSD does not assert a bounded cusp link; no finite-volume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['types'][0]['counts'])
