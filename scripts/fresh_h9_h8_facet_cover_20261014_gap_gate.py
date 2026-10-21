#!/usr/bin/env python3
"""Exact all-label projected-normal gap in the fixed 11-wall H8 source."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_facet_cover_20261014_gap_gate'
SOURCE='data/known_examples/bugaenko_compact_h8.json';BASE='results/fresh_h9_mixed_d4_coupling_20261014.json';CAP='results/fresh_h9_mixed_d4_global_caps_20261014.json'
Q=A.Q;Z=Q();ONE=Q(1);PHI=Q(F(1,2),F(1,2));SMALL=Q(F(-1,4),F(1,4));BOUND=4*SMALL
# BOUND=2(phi-1)=sqrt5-1, not 2(sqrt5-1).
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def sf(x):return Q(F(x['a'])+F(x['b'])/2,F(x['b'])/2)
def dot(x,y):return sum((a*b for a,b in zip(x,y)),Z)
def mv(M,x):return [dot(row,x)for row in M]
def inv(M):
 n=len(M);R=[row[:]+[Q(i==j)for j in range(n)]for i,row in enumerate(M)]
 for j in range(n):
  p=next((i for i in range(j,n)if R[i][j]),None)
  if p is None:return None
  R[j],R[p]=R[p],R[j];d=R[j][j];R[j]=[x/d for x in R[j]]
  for i in range(n):
   if i!=j and R[i][j]:
    d=R[i][j];R[i]=[x-d*y for x,y in zip(R[i],R[j])]
 return [row[n:]for row in R]
def enc(x):
 if isinstance(x,Q):return x.enc()
 if isinstance(x,(tuple,list)):return [enc(y)for y in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def run():
 src=json.loads((ROOT/SOURCE).read_text());assert src['dimension']==8 and len(src['simple_roots'])==11 and len(src['geometry']['vertices'])==41 and src['geometry']['cusps']==[]
 G=[[sf(x)for x in row]for row in src['wall_gram']];N=[[sf(x)for x in row]for row in src['simple_roots']];form=[[sf(x)for x in row]for row in src['lattice_gram']]
 assert all(dot(x,mv(form,y))==G[i][j]for i,x in enumerate(N)for j,y in enumerate(N))
 assert all(G[i][i]==2 for i in range(11))
 total=singular=feasible=0;unique={};bases=[]
 for ids in combinations(range(11),9):
  total+=1;M=inv([[G[i][j]for j in ids]for i in ids])
  if M is None:singular+=1;continue
  other=[i for i in range(11)if i not in ids]
  C=[[dot([G[i][k]for k in ids],[M[k][j]for k in range(9)])for j in range(9)]for i in other]
  norm=[Z]*512;outs=[[Z,Z]for _ in range(512)];accepted=0
  for mask in range(1,512):
   bit=mask&-mask;i=bit.bit_length()-1;rest=mask-bit
   norm[mask]=norm[rest]-M[i][i]-2*sum((M[i][j]for j in range(9)if(rest>>j)&1),Z)
   outs[mask]=[outs[rest][k]+C[k][i]for k in range(2)]
   if any(x and (x-ONE).sign()<0 for x in outs[mask]):continue
   h=[Z]*11
   for j,k in enumerate(ids):h[k]=Q((mask>>j)&1)
   for j,k in enumerate(other):h[k]=outs[mask][j]
   assert norm[mask].sign()>0
   zeros=tuple(j for j,x in enumerate(h)if not x);key=tuple(h);feasible+=1;accepted+=1
   if key in unique:assert unique[key][0]==norm[mask]
   unique[key]=(norm[mask],zeros)
  bases.append(dict(ids=list(ids),accepted=accepted))
 assert total==55 and singular==12 and len(bases)==43
 exceptional=[];ordinary_min=None;ordinary_count=0;patterns={}
 for h,(n,z)in unique.items():
  patterns.setdefault(z,[]).append(n)
  if (n-BOUND).sign()<0:
   assert n==SMALL and len(z)==8
   exceptional.append(dict(zero_set=z,contacts=h,norm=n,source_vertex=src['geometry']['vertices'].index(list(z))))
  else:
   ordinary_count+=1
   if ordinary_min is None or(n-ordinary_min).sign()<0:ordinary_min=n
 assert len(exceptional)==3 and sorted(x['source_vertex']for x in exceptional)==[0,32,33]
 assert ordinary_min==BOUND
 special={x['zero_set']for x in exceptional}
 assert all(z in special or all((v-BOUND).sign()>=0 for v in vals)for z,vals in patterns.items())
 # Every exceptional pattern has8 independent facet equations, hence one ray.
 for row in exceptional:
  z=row['zero_set'];assert inv([[G[i][j]for j in z]for i in z])is not None
  assert any(x==ONE for x in row['contacts'])
 b,d=[json.loads((ROOT/p).read_text())for p in[BASE,CAP]]
 Fm=A.dm(b['form']);old=A.dm(b['final_normals']);r=A.dv(d['raw_cap']);q=A.dec(d['raw_cap_norm'])
 gg=[A.pair(r,n,Fm)for n in old];nn=[A.pair(n,n,Fm)for n in old]
 kapp=[gg[i]*gg[i]/q-nn[i]for i in range(12)]
 assert all(nn[i]==2 for i in[0,1,2,3,7])
 assert all((A.E(BOUND)-kapp[i]).sign()>0 for i in[0,1,2,3])
 distinct=[]
 for i,j in combinations(range(4),2):
  p=A.pair(old[i],old[j],Fm)-gg[i]*gg[j]/q;det=p*p-kapp[i]*kapp[j]
  assert det.sign()>0;distinct.append(dict(pair=[i,j],positive_squared_distance_margin=A.enc(det)))
 # Wall7 gives an independent consecutive-cosine gap (1<t<sqrt2).
 assert (kapp[7]-A.E(SMALL)).sign()>0 and (2*A.E(SMALL)-kapp[7]).sign()>0
 paths=[SOURCE,BASE,CAP,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
 rows=sorted([(enc(h),enc(n),list(z))for h,(n,z)in unique.items()],key=str)
 return dict(status='PASS_FIXED_LITERAL_CAP_H8_SOURCE_ATTACHMENT_EXCLUDED',sources={p:sha(p)for p in paths},
  source_dimension=8,source_facets=11,source_finite_vertices=41,total_nine_subsets=total,singular_subsets=singular,independent_subsets=43,binary_rhs_per_subset=512,
  feasible_vertex_occurrences=feasible,unique_gap_polyhedron_vertices=len(unique),unique_zero_patterns=len(patterns),ordinary_vertex_count=ordinary_count,
  ordinary_minimum_norm=enc(ordinary_min),exceptional_minimum_norm=enc(SMALL),exceptional_rays=enc(sorted(exceptional,key=lambda x:x['source_vertex'])),
  bases=bases,canonical_vertex_digest=sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),
  literal_projected_norms=A.enc(kapp),four_forced_source_vertex_centers=[0,1,2,3],distinctness_margins=distinct,
  independent_wall7_gap=dict(lower=A.enc(kapp[7]-A.E(SMALL)),upper=A.enc(2*A.E(SMALL)-kapp[7])),
  scope='All isometric placements of this fixed compact11-facetH8 source perpendicular to the literal cap, retaining all12 old normals with acute Coxeter pairs. Norm-gap enumeration is complete over all real projected normals and all labels. Does not exclude different cap parameters, different H8 sources, nonperpendicular attachment or changed old walls.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['unique_gap_polyhedron_vertices'],'vertices',d['unique_zero_patterns'],'patterns',flush=True)
