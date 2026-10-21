#!/usr/bin/env python3
"""Complete PD8 edge census, including five edges without finite endpoints."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import argparse,json
import fresh_h9_joint_frontier_20261013 as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_frontier_20261013_all_edges'
BASE=ROOT/'results/fresh_h9_joint_frontier_20261013.json';OUT=ROOT/'results'/f'{TAG}.json'
def run():
 base=json.loads(BASE.read_text())
 for p,h in base['sources'].items():assert A.sha(ROOT/p)==h
 s,Q,W,C,c,G,adj=A.setup();assert [A.encv(w)for w in W]==base['normalized_normals']
 cliques=A.clique_list(adj,8);assert len(cliques)==1105
 pd=[ids for ids in cliques if A.pd([[G[i][j]for j in ids]for i in ids])];assert len(pd)==446
 old={tuple(e['support']):e for e in base['edges']};vmap={tuple(v['active']):v['id']for v in base['finite_vertices']}
 new=[ids for ids in pd if ids not in old];assert len(new)==5
 rows=[];classification=Counter();all_supports=[]
 for ids in pd:
  M=[[G[i][j]for j in ids]for i in ids];B=A.ki(M);finite=[];affine=[]
  for j in range(27):
   if j in ids or any(j not in adj[i]for i in ids):continue
   b=[G[i][j]for i in ids];q=A.I.sub(A.K(2),A.sumk(A.I.mul(A.I.mul(x,B[i][k]),y)for i,x in enumerate(b)for k,y in enumerate(b)))
   if A.I.sign(q)>0:
    key=tuple(sorted(ids+(j,)));assert key in vmap;finite.append(vmap[key])
   elif q==A.Z:affine.append(j)
  assert not affine
  if ids in old:
   assert sorted(finite)==sorted(old[ids]['finite_vertices']);kind='two_finite_ends'if len(finite)==2 else'one_finite_end'
   all_supports.append({'support':list(ids),'kind':kind,'finite_vertices':finite,'finite_frontier_edge_id':old[ids]['id']})
   classification[kind]+=1;continue
  assert not finite;classification['no_finite_ends']+=1
  # Projection of the strict center to a PD face is timelike and strict
  # on every other old inequality because the acute inverse is nonnegative.
  assert all(A.I.sign(x)>=0 for row in B for x in row)
  cp=[A.dot(c,C[i])for i in ids];coeff=[A.total(A.scale(q,x)for q,x in zip(row,cp,strict=True))for row in B]
  p=[A.sub(c[k],A.total(A.mul(t,W[j][k])for t,j in zip(coeff,ids,strict=True)))for k in range(10)]
  pnorm=A.pair(p,p,Q);assert A.sgn(pnorm)<0 and A.sgn(A.pair(p,c,Q))<0
  signs=[A.sgn(A.dot(p,x))for x in C];assert max(signs)<=0 and [i for i,x in enumerate(signs)if not x]==list(ids)
  # Pick a coordinate direction, project to the Lorentz plane, then remove
  # its p-component without a division in the radical algebra.
  for axis in range(10):
   e=[A.C(A.O)if i==axis else{}for i in range(10)];b=[A.dot(e,C[i])for i in ids]
   coeff=[A.total(A.scale(q,x)for q,x in zip(row,b,strict=True))for row in B]
   z=[A.sub(e[k],A.total(A.mul(t,W[j][k])for t,j in zip(coeff,ids,strict=True)))for k in range(10)]
   pp=A.pair(z,p,Q);direction=[A.add(A.mul(A.neg(pnorm),x),A.mul(pp,y))for x,y in zip(z,p,strict=True)]
   dnorm=A.pair(direction,direction,Q)
   if A.sgn(dnorm)>0:break
  else:raise AssertionError('no spacelike line direction')
  assert not A.pair(p,direction,Q) and all(not A.dot(direction,C[i])for i in ids)
  critical=(-A.approx(pnorm)/A.approx(dnorm))**.5;ends=[]
  for side in [-1,1]:
   found=None
   for digits in range(2,16):
    den=10**digits;t=F(int(critical*den)+1,den)
    y=[A.add(x,A.scale(A.K(side*t),z))for x,z in zip(p,direction,strict=True)]
    ss=[A.sgn(A.dot(y,x))for x in C]
    if A.sgn(A.pair(y,y,Q))>0 and max(ss)<=0 and [i for i,x in enumerate(ss)if not x]==list(ids):found=t,y,ss;break
   assert found is not None
   t,y,ss=found;assert A.sgn(A.pair(y,c,Q))<0
   ends.append({'side':side,'parameter':str(t),'point':A.encv(y),'norm':A.enc(A.pair(y,y,Q)),'all_wall_signs':ss,'center_pairing':A.enc(A.pair(y,c,Q))})
  index=441+len(rows);profile=A.spherical_profiles(M,ids)
  rows.append({'id':index,'support':list(ids),'interior_point':A.encv(p),'interior_norm':A.enc(pnorm),'interior_wall_signs':signs,'direction_coordinate_proposal':axis,'spacelike_direction':A.encv(direction),'direction_norm':A.enc(dnorm),'spacelike_continuations':ends,'spherical_profiles':profile})
  all_supports.append({'support':list(ids),'kind':'no_finite_ends','finite_vertices':[],'new_edge_id':index})
 counts={'finite_profiles':sum(len(r['spherical_profiles']['finite_profiles'])for r in rows),'dihedral_tails':sum(len(r['spherical_profiles']['dihedral_tail_profiles'])for r in rows)}
 return {'status':'PASS_COMPLETE_ALL_GEODESIC_EDGES_27_WALL_CHAMBER','sources':{str(p.relative_to(ROOT)):A.sha(p)for p in [BASE,A.SOURCE,Path(A.__file__),Path(__file__)]},'finite_pair_eight_cliques':len(cliques),'not_positive_definite':len(cliques)-len(pd),'actual_PD8_edges':len(pd),'edge_classification':dict(classification),'ideal_endpoint_extension_count':0,'all_edge_supports':all_supports,'new_edges_without_finite_endpoint':rows,'new_profile_counts':counts,'scope':'Complete census of all actual geodesic edges of the frozen 27-wall H9 chamber: a PD8 acute support is actual by strict-center projection, and every interior edge has eight independent incident walls. Each new line has exact spacelike continuations in both directions. Necessary spherical profiles are complete locally but do not establish globally compatible new caps; affine cap profiles are not enumerated.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
 if a.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['edge_classification'],out['new_profile_counts'],flush=True)
