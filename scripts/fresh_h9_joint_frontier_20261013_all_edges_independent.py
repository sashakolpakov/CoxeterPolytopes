#!/usr/bin/env python3
"""Separate replay of the exhaustive PD8 census and ten remote open ends."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations,product
import argparse,json
import fresh_h9_joint_frontier_20261013_independent as I
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_frontier_20261013_all_edges_independent'
SOURCE=ROOT/'results/fresh_h9_joint_frontier_20261013_all_edges.json';OUT=ROOT/'results'/f'{TAG}.json'
def run():
 prior=I.run();d=json.loads(SOURCE.read_text())
 for p,h in d['sources'].items():assert I.sha(ROOT/p)==h
 base=json.loads(I.SOURCE.read_text());src=json.loads(I.WALLS.read_text())
 Q=[[I.K(x)for x in row]for row in src['form']];W=[[I.dec(x)for x in row]for row in base['normalized_normals']]
 C=[I.cov(w,Q)for w in W];c=[I.con(I.K(x))for x in src['strict_center']]
 two=(F(2),F(0));phi=(F(1,2),F(1,2));G=[[None]*27 for _ in W];adj=[set()for _ in W]
 for i in range(27):G[i][i]=two
 for row in src['all_pair_gates']:
  if row['type']=='dotted':continue
  i,j=row['pair'];G[i][j]=G[j][i]={'m2':I.Z,'m3':I.kn(I.ONE),'m5':I.kn(phi)}[row['type']];adj[i].add(j);adj[j].add(i)
 stack=[((),set(range(27)))];count=0;pdsets=[]
 while stack:
  ids,rest=stack.pop()
  if len(ids)==8:
   count+=1
   if I.positive([[G[i][j]for j in ids]for i in ids]):pdsets.append(ids)
   continue
  if len(ids)+len(rest)<8:continue
  for i in rest:stack.append((ids+(i,),{j for j in rest&adj[i]if j>i}))
 assert count==1105 and len(pdsets)==446
 saved={tuple(e['support']):e for e in d['all_edge_supports']};assert set(pdsets)==set(saved)
 vmap={tuple(v['active']):v['id']for v in base['finite_vertices']};classes=Counter()
 for ids in pdsets:
  B=I.inverse([[G[i][j]for j in ids]for i in ids]);assert all(I.ks(x)>=0 for row in B for x in row)
  finite=[]
  for j in range(27):
   if j in ids or any(j not in adj[i]for i in ids):continue
   b=[G[i][j]for i in ids];q=I.ka(two,I.kn(I.ksum(I.km(I.km(x,B[i][k]),y)for i,x in enumerate(b)for k,y in enumerate(b))))
   assert q!=I.Z
   if I.ks(q)>0:finite.append(vmap[tuple(sorted(ids+(j,)))])
  row=saved[ids];assert sorted(finite)==sorted(row['finite_vertices']);assert len(finite)in[0,1,2]
  tag=['no_finite_ends','one_finite_end','two_finite_ends'][len(finite)];assert row['kind']==tag;classes[tag]+=1
 assert classes==d['edge_classification']
 new=d['new_edges_without_finite_endpoint'];assert len(new)==5 and {tuple(r['support'])for r in new}=={ids for ids,row in saved.items()if row['kind']=='no_finite_ends'}
 profile_count=tail_count=0
 for row in new:
  ids=row['support'];p=[I.dec(x)for x in row['interior_point']];v=[I.dec(x)for x in row['spacelike_direction']]
  pnorm=I.pair(p,p,Q);vnorm=I.pair(v,v,Q)
  assert pnorm==I.dec(row['interior_norm']) and I.sign(pnorm)<0 and I.sign(I.pair(p,c,Q))<0
  assert vnorm==I.dec(row['direction_norm']) and I.sign(vnorm)>0 and not I.pair(p,v,Q)
  ps=[I.sign(I.dot(p,z))for z in C];assert ps==row['interior_wall_signs'] and max(ps)<=0 and [i for i,x in enumerate(ps)if not x]==ids
  assert all(not I.dot(v,C[i])for i in ids)
  assert sorted(e['side']for e in row['spacelike_continuations'])==[-1,1]
  for end in row['spacelike_continuations']:
   t=F(end['parameter']);assert t>0
   y=[I.dec(x)for x in end['point']];assert y==[I.add(a,I.scale((end['side']*t,F(0)),b))for a,b in zip(p,v,strict=True)]
   q=I.pair(y,y,Q);assert q==I.dec(end['norm']) and I.sign(q)>0 and I.sign(I.pair(y,c,Q))<0
   signs=[I.sign(I.dot(y,z))for z in C];assert signs==end['all_wall_signs'] and max(signs)<=0 and [i for i,x in enumerate(signs)if not x]==ids
  # Recompute the spherical profiles using the independently replayed
  # normalized Gram and its inverse; exact labels need no coordinate field.
  M=[[G[i][j]for j in ids]for i in ids];B=I.inverse(M);left=set(range(8));components=[]
  while left:
   seen={min(left)};queue=list(seen)
   for i in queue:
    for j in sorted(left-seen):
     if M[i][j]!=I.Z:seen.add(j);queue.append(j)
   components.append(sorted(seen));left-=seen
  which={i:k for k,cc in enumerate(components)for i in cc};squares={3:I.ONE,4:two,5:(F(3,2),F(1,2))};found={():I.Z};tails=[]
  for n in [1,2,3]:
   for indices in combinations(range(8),n):
    if len({which[i]for i in indices})<n:continue
    if n==1 and len(components[which[indices[0]]])==1:tails.append(ids[indices[0]]);continue
    for labs in product([3,4,5],repeat=n)if n<3 else [(3,3,3)]:
     proj=I.ksum(I.km(squares[m],B[i][i])for i,m in zip(indices,labs,strict=True))
     if I.ks(I.ka(two,I.kn(proj)))>0:found[tuple((ids[i],m)for i,m in zip(indices,labs,strict=True))]=proj
  profile=row['spherical_profiles'];actual={tuple(map(tuple,r['contacts'])):I.K(r['projection_square'])for r in profile['finite_profiles']}
  assert found==actual and len(found)==len(profile['finite_profiles'])
  assert sorted(tails)==sorted(x['wall']for x in profile['dihedral_tail_profiles'])
  for r in profile['finite_profiles']:assert I.K(r['residual_norm'])==I.ka(two,I.kn(I.K(r['projection_square'])))
  profile_count+=len(found);tail_count+=len(tails)
 assert [profile_count,tail_count]==[143,14]
 return {'status':'PASS_INDEPENDENT_COMPLETE_PD8_AND_REMOTE_END_REPLAY','sources':{str(p.relative_to(ROOT)):I.sha(p)for p in [SOURCE,I.SOURCE,I.WALLS,Path(I.__file__),Path(__file__)]},'prior_frontier_status':prior['status'],'finite_pair_eight_cliques':1105,'actual_geodesic_edges':446,'edge_classification':dict(classes),'new_vertex_free_lines':5,'exact_spacelike_continuations':10,'new_finite_profiles':profile_count,'new_symbolic_dihedral_tails':tail_count,'scope':'Complete PD8 incidence and existing endpoint classification; no numerical geometry or producer arithmetic imports. Every new line has two exact open ends. This does not assert existence of a compatible completing cap.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
 if a.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['edge_classification'],flush=True)
