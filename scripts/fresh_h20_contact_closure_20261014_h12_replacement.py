#!/usr/bin/env python3
"""Delete old wall0 and transplant the reversed full H12 A12 source.
Exact77-wall construction; its E8 face is finite volume, ambient volume is not.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
import fresh_h20_contact_closure_20261014_replacement as b
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014_h12_replacement'
BASE='results/fresh_h20_contact_closure_20261014_replacement.json'
FACE='results/fresh_h20_contact_closure_20261014_e8_face.json'
TRANSFER='results/fresh_h20_contact_closure_20261014_h12_transfer.json'
SOURCE='results/fresh_h12_a12_source_20261014.json';Z=b.ZZ
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d,f,tr,src=[json.loads((ROOT/p).read_text())for p in [BASE,FACE,TRANSFER,SOURCE]]
 for data in [d,f,tr,src]:
  for p,h in data['sources'].items():assert sha(p)==h,p
 Q=[[F(x)for x in row]for row in d['form']];old=[[tuple(map(F,x))for x in row]for row in d['normals']]
 caps=[[tuple(map(F,x))for x in row['normal']]for row in tr['orientations'][1]['caps']]
 N=old[1:]+caps;H=[b.c(F(x))for x in d['strict_center']];escape=[b.c(F(x))for x in d['strict_escape']['point']]
 cache={}
 def dot(x,y):
  key=tuple(y)
  if key not in cache:cache[key]=[tuple(sum(Q[i][j]*y[j][k]for j in range(21))for k in range(4))for i in range(21)]
  out=Z
  for a,v in zip(x,cache[key]):out=b.add(out,b.mul(a,v))
  return out
 def linear(vectors,coeff):
  out=[Z]*21
  for c,v in zip(coeff,vectors):out=[b.add(x,b.scale(y,c))for x,y in zip(out,v)]
  return out
 norms=[]
 for r in N:
  q=dot(r,r);assert q[1:]==(0,0,0)and q[0]>0;norms.append(q[0])
 G=[[dot(r,s)for s in N]for r in N];allcounts=Counter();newcounts=Counter()
 for i in range(77):
  for j in range(i):
   g=G[i][j];assert b.sign(g)<=0;sq=b.mul(g,g);prod=norms[i]*norms[j]
   if g==Z:kind='right'
   elif sq==b.c(prod/4):kind='pi3'
   elif sq==b.c(prod/2):kind='pi4'
   elif sq==b.c(prod):kind='parallel'
   else:assert b.sign(b.add(sq,b.c(-prod)))>0;kind='dotted'
   allcounts[kind]+=1
   if i>=75:newcounts[kind]+=1
 assert sum(allcounts.values())==2926 and sum(newcounts.values())==151
 assert b.sign(dot(H,H))<0 and all(b.sign(dot(H,r))<0 for r in N)
 facets=[]
 for i,r in enumerate(N):
  k=b.scale(dot(H,r),-1/norms[i]);p=[b.add(x,b.mul(k,y))for x,y in zip(H,r)]
  assert b.sign(dot(p,p))<0 and b.sign(dot(p,H))<0
  assert all(dot(p,s)==Z if j==i else b.sign(dot(p,s))<0 for j,s in enumerate(N));facets.append(p)
 assert b.sign(dot(escape,escape))>0 and b.sign(dot(escape,H))<0 and all(b.sign(dot(escape,r))<0 for r in N)
 # The oldE8 complement form/basis is unchanged. Reverse A12 chain orientation.
 E=[i-1 for i in f['E8_walls']];A=[i-1 for i in reversed(f['orthogonal_face_walls'])]
 T=[tuple(map(F,x))for x in f['normalized_timelike_basis']]
 source_chain=list(range(11))+[13];mapping=dict(zip(source_chain,A));mapping.update({12:75,11:76})
 assert all(G[mapping[i]][mapping[j]]==b.c(F(src['gram'][i][j]))for i in range(14)for j in range(14))
 C=[[F(src['gram'][i][j])for j in source_chain]for i in source_chain];CI=b.rat.inv(C)
 srcnormals=[[F(x)for x in row]for row in src['normals']];srcT=list(map(F,src['normalized_tangent_time']))
 def srcdot(x,y):return -x[0]*y[0]+sum(a*v for a,v in zip(x[1:],y[1:]))
 def lift(x):
  pairs=[srcdot(srcnormals[i],x)for i in source_chain];coefs=[sum(a*v for a,v in zip(row,pairs))for row in CI]+[-13*srcdot(srcT,x)]
  return linear([N[i]for i in A]+[T],coefs)
 assert all(lift(srcnormals[i])==N[mapping[i]]for i in range(14))
 finite=[];ideal=[]
 for row in src['projective_vertices']:
  p=lift(list(map(F,row['point'])));inc=sorted(E+[mapping[i]for i in row['incident']]);q=dot(p,p)
  assert q==b.c(F(row['norm']))and b.sign(dot(p,H))<0
  assert all(dot(p,r)==Z if i in inc else b.sign(dot(p,r))<0 for i,r in enumerate(N))
  record=dict(point=p,norm=q,incident=inc,source_incident=row['incident'])
  if b.sign(q)<0:
   assert len(inc)==20
   assert all(G[i][j][1:]==(0,0,0)for i in inc for j in inc)
   b.rat.pd([[G[i][j][0]for j in inc]for i in inc]);finite.append(record)
  else:assert q==Z;ideal.append(record)
 assert len(finite)==35 and len(ideal)==2
 EI=[[F(x)for x in row]for row in f['inverse']]
 negative_weight=linear([N[i]for i in E],[-x for x in EI[0]])
 assert dot(negative_weight,negative_weight)==b.c(2)
 assert all(dot(negative_weight,N[i])==b.c(-int(k==0))for k,i in enumerate(E))
 for row in ideal:
  assert dot(negative_weight,row['point'])==Z
  assert all(b.sign(dot(negative_weight,N[i]))<=0 for i in row['incident'])
  row['E8_link_recession']=negative_weight
 # Fresh witnesses survive deletion; do not inherit the old wall0 cycle.
 irrpair=[17,21];ratpair=[23,24]
 i,j=irrpair;cycle=b.scale(b.mul(G[i][j],G[i][j]),4/(norms[i]*norms[j]));assert cycle==(F(12),F(8),F(0),F(0))
 i,j=ratpair;assert [[G[a][v]for v in ratpair]for a in ratpair]==[[b.c(2),b.c(-6)],[b.c(-6),b.c(2)]]
 paths=[BASE,FACE,TRANSFER,SOURCE,f'scripts/{TAG}.py','scripts/fresh_h20_contact_closure_20261014_replacement.py','scripts/fresh_h20_coupled_lift_20261014.py']
 return b.enc(dict(status='PASS_77_WALL_CHANGED_CORE_WITH_FINITE_VOLUME_H12_FACE_AND_STRICT_ESCAPE',sources={p:sha(p)for p in paths},field_basis=d['field_basis'],form=Q,normals=N,norms=norms,deleted_predecessor_wall=0,retained_predecessor_walls=list(range(1,76)),new_cap_indices=[75,76],all_pair_counts=dict(allcounts),new_pair_counts=dict(newcounts),strict_center=H,facet_points=facets,E8_support=E,E8_inverse=EI,E8_minimum_weight_node=E[0],source_wall_map=mapping,H12_face_finite_vertices=finite,H12_face_ideal_vertices=ideal,source_face_volume='Finite, by exact isometry to the separately certified14-wall I12,1 source. All other ambient projections are future timelike.',strict_escape=dict(point=escape,norm=dot(escape,escape),center_pairing=dot(escape,H),new_pairings=[dot(escape,r)for r in caps]),conditional_nonquasi=dict(irrational_cycle_pair=irrpair,cycle=cycle,rational_indefinite_pair=ratpair,gram=[[2,-6],[-6,2]],determinant=-32),scope='All77 actual facets and all2926 Coxeter pairs. One E8 face is completed to the finite-volume H12 source, giving35 certified ambient finite vertices. Its2 ideal vertices still have an E8 recession direction in their H20 links. No complete ambient corner census. The saved76 strictescape survives, so ambient volume remains infinite.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'])
