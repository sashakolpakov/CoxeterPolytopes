#!/usr/bin/env python3
"""Three coherent null-source caps with bounded full H20 cusp links."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h20_new_projection_20261015_exact as A
Z=A.Z;P=A.P;K=A.K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_new_projection_20261015_cusps'
SOURCE='results/fresh_h20_new_projection_20261015_exact.json';HELPER='scripts/fresh_h20_new_projection_20261015_exact.py'
def dec(x):return Z(P(K(F(x[0]),F(x[1])),K(F(x[2]),F(x[3]))),P(K(F(x[4]),F(x[5])),K(F(x[6]),F(x[7]))))
dot=A.dot;mv=A.mv;enc=A.enc
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rref(M):
 a=[r[:]for r in M];piv=[];i=0
 for j in range(len(a[0])):
  k=next((k for k in range(i,len(a))if a[k][j]),None)
  if k is None:continue
  a[i],a[k]=a[k],a[i];v=a[i][j];a[i]=[x/v for x in a[i]]
  for k in range(len(a)):
   if k!=i and a[k][j]:
    v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[i])]
  piv.append(j);i+=1
  if i==len(a):break
 return a,piv
def pd(M):
 a=[r[:]for r in M];p=[]
 for i in range(len(a)):
  d=a[i][i];assert d.sign()>0;p.append(d)
  for j in range(i+1,len(a)):
   for k in range(i+1,len(a)):a[j][k]-=a[j][i]*a[i][k]/d
 return p
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 Q=[[dec(x)for x in r]for r in d['form']];N=[[dec(x)for x in r]for r in d['normals']]
 C=list(map(dec,d['strict_center']));E=list(map(dec,d['strict_escape']))
 pairs=[(17,23),(17,24),(23,24)]
 rays=[[x+y for x,y in zip(N[i],N[j])]for i,j in pairs]
 assert all(not dot(c,mv(Q,c)) and dot(c,mv(Q,C)).sign()<0 for c in rays)
 assert all(dot(c,mv(Q,e))==-2 for c,e in combinations(rays,2))
 W=[Z()]*19+[Z(P(K(0,F(-1,2))))]*2
 assert dot(W,mv(Q,W))==2
 for c in rays:N.append([x+y for x,y in zip(c,W)])
 cov=[mv(Q,n)for n in N];G=[[dot(n,c)for c in cov]for n in N];norm=[G[i][i]for i in range(47)]
 assert all(q==2 for q in norm[44:])
 counts=Counter();allcounts=Counter();labels={}
 for i,j in combinations(range(47),2):
  g=G[i][j];assert g.sign()<=0;t=4*g*g/(norm[i]*norm[j])
  if not g:lab='right'
  elif t==1:lab='pi3'
  elif t==2:lab='pi4'
  elif t==3:lab='pi6'
  elif t==4:lab='parallel'
  else:assert(t-4).sign()>0;lab='dotted'
  labels[i,j]=lab;allcounts[lab]+=1
  if j>=44:counts[lab]+=1
 assert sum(counts.values())==135
 assert all(not G[i][j]for i,j in combinations(range(44,47),2))
 assert all(labels[i,j]=='pi4'for i in (40,41)for j in (44,45,46))
 cp=[dot(C,c)for c in cov];qC=dot(C,mv(Q,C));assert qC.sign()<0 and all(p.sign()<0 for p in cp)
 facets=[]
 for i,n in enumerate(N):
  v=[x-cp[i]/norm[i]*y for x,y in zip(C,n)]
  q=qC-cp[i]*cp[i]/norm[i];assert q.sign()<0
  products=[cp[j]-cp[i]*G[i][j]/norm[i]for j in range(47)]
  assert [x.sign()for x in products]==[-1]*i+[0]+[-1]*(46-i)
  facets.append(dict(point=v,norm=q,pairings=products))
 cusp_rows=[]
 for k,c in enumerate(rays):
  products=[dot(c,n)for n in cov];assert all(x.sign()<=0 for x in products)
  incident=[i for i,x in enumerate(products)if not x];assert incident==[i for i in range(40)if not products[i]]+[40,41,44+k]
  assert len(incident)==23
  gram=[[G[i][j]for j in incident]for i in incident]
  _,piv=rref(gram);assert len(piv)==19
  _,span=rref([N[i]for i in incident]);assert len(span)==20
  pending=set(incident);components=[]
  while pending:
   start=min(pending);ids=[start];pending.remove(start)
   for i in ids:
    for j in sorted(pending.copy()):
     if G[i][j]:ids.append(j);pending.remove(j)
   ids.sort();g=[[G[i][j]for j in ids]for i in ids];r,p=rref(g)
   assert len(p)==len(ids)-1
   free=next(j for j in range(len(ids))if j not in p);weights=[Z()]*len(ids);weights[free]=Z(1)
   for a,j in enumerate(p):weights[j]=-r[a][free]
   if weights[0].sign()<0:weights=[-x for x in weights]
   assert all(x.sign()>0 for x in weights)
   principal=[row[:-1]for row in g[:-1]];pivots=pd(principal)
   v=[sum((a*N[i][j]for a,i in zip(weights,ids)),Z())for j in range(21)]
   ix=next(j for j,x in enumerate(c)if x);mu=v[ix]/c[ix]
   assert mu.sign()>0 and all(x==mu*y for x,y in zip(v,c))
   components.append(dict(walls=ids,rank=len(ids)-1,positive_kernel=weights,principal_pd_pivots=pivots,physical_kernel_multiple=mu))
  assert sorted(v['rank']for v in components)==[1,2,8,8]
  assert any(v['walls']==[40,41,44+k]for v in components)
  # A norm2 source wall with product -1 pins every positive scale t>=1,
  # regardless of the chosen norm2 transverse profile.
  anchor=next(i for i in range(40)if norm[i]==2 and products[i]==-1)
  finite_anchors=[i for i in range(40)if labels.get((i,44+k))=='pi3']
  assert len(finite_anchors)==4
  source_E=dot(c,mv(Q,E));assert source_E.sign()<0 and (source_E*source_E-2).sign()>0
  cusp_rows.append(dict(source_pair=pairs[k],ray=c,pairings=products,incident=incident,gram_rank=19,normal_span=20,positive_affine_components=components,cap_index=44+k,scale=1,minimum_scale_anchor=anchor,pi3_source_anchors=finite_anchors,source_escape_pairing=source_E,all_profile_escape_squared_margin=source_E*source_E-2))
 old_points=d['finite_edge']['endpoints']+d['retained_previous_edge_vertices'];retained=[]
 for r in old_points:
  v=list(map(dec,r['point']));p=[dot(v,c)for c in cov]
  assert dot(v,mv(Q,v)).sign()<0 and all(x.sign()==(0 if i in r['support']else -1)for i,x in enumerate(p))
  retained.append(dict(support=r['support'],point=v,pairings=p))
 ep=[dot(E,c)for c in cov];qE=dot(E,mv(Q,E));assert qE.sign()>0 and all(x.sign()<0 for x in ep)
 assert dot(E,mv(Q,C)).sign()<0
 paths=[f'scripts/{TAG}.py',SOURCE,HELPER]
 return enc(dict(status='PASS_EXACT47_WALL_THREE_BOUNDED_CUSPS_WITH_RETAINED_ESCAPE',sources={p:sha(p)for p in paths},field_basis=d['field_basis'],dimension=20,facets=47,
  form=Q,normals=N,gram=G,new_pair_counts=dict(counts),all_pair_counts=dict(allcounts),strict_center=C,center_norm=qC,center_pairings=cp,facet_witnesses=facets,
  cap_recipe=dict(source_pairs=pairs,common_transverse=W,scales=[1,1,1],mutual_new_cap_products=[0,0,0]),repaired_cusps=cusp_rows,retained_finite_vertices=retained,
  strict_escape=E,escape_norm=qE,escape_pairings=ep,escape_center_pairing=dot(E,mv(Q,C)),
  all_profile_three_ray_gate=dict(source_ray_count=3,transverse_escape_norm=1,inequality='For r=t*c+w, r²=2, w in the transverse plane, acuteness/Coxeter at the saved norm2 anchor forces t>=1. Since E_trans²=1, (w,E_trans)<=sqrt2. The exact c·E<-sqrt2 implies (r,E)<0 for every such profile, including every finite-angle alphabet. This gate is only for the three fixed source rays.'),
  scope='Three specific null-source caps at their least possible positive scale t1 retain all44 walls and give three actual bounded rank19 cusp links, each a product of affine ranks8,8,1,2. All47 facets and four retained finite vertices are checked. No complete source or ambient cusp/vertex census. The same exact strict spacelike escape proves infinite volume. The separate three-ray gate covers arbitrary transverse profiles on these fixed rays, not new rays or general null/interior source projections.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],[r['incident']for r in d['repaired_cusps']])
