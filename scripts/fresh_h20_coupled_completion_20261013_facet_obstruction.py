#!/usr/bin/env python3
"""Verify the18-neighbor hypothesis for the dotted-only completion obstruction.

In r-perp, dimension20 and signature(19,1), project any strict completion center
onto r-perp to obtain c. A nonzero w orthogonal to c and the18 intersecting-wall
projections exists and is spacelike. The future null c+sqrt(-c^2/w^2)*w is strict
on all18 inequalities. A wall dotted to r projects to a future timelike normal,
so is strict on that null point too. A small outward perturbation within r-perp,
then an inward perturbation across r, gives a strict spacelike escape. This
works for every finite proposed addition consisting only of walls dotted to r;
the strict center used in the proof is that proposed completion's center.
"""
import argparse,json
from collections import Counter
from fresh_h20_coupled_completion_20261013_second_exact import setup,OUT as SOURCE,ROOT,read,sha,e,emb,pairing,norm
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,label,same
from fresh_dimension_bridge_20261013_verify import dot
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_facet_obstruction.json'
def run():
 Q,h,hc,roots,f=setup();d=read(SOURCE);g=Extension(f,f.parse(d['field']['radicand']));r=tuple(map(g.parse,d['new_cap']['normal']));qr=norm(g,r,Q);assert same(g,qr,g.rational(1))
 rc=fmv(g,tuple(tuple(e(g,x)for x in row)for row in Q),r);hp=fdot(g,r,tuple(e(g,x)for x in hc));assert g.sign(hp)<0
 cq=g.sub(e(g,dot(h,hc)),g.mul(hp,hp));assert g.sign(cq)<0
 rows=[];neighbors=[];dotted=[]
 for i,(f0,v)in enumerate(roots):
  H,p=pairing(g,rc,f0,v);qv=norm(f0,v,Q)
  if H is g:qold=emb(g,f0,qv)
  else:qold=(e(g,qv[0]),e(g,qv[1]))
  kind=label(H,p,H.rational(1),qold);assert kind==d['new_cap']['old_pair_checks'][i]['kind']
  qp=H.sub(qold,H.mul(p,p));sg=H.sign(qp)
  if kind=='dotted':assert sg<0;dotted.append(i)
  else:assert kind in('2','4','6')and sg>0;neighbors.append(i)
  rows.append({'wall':i,'kind':kind,'projection_norm':H.encode(qp),'projection_norm_sign':sg})
 assert len(roots)==175 and len(neighbors)==18 and len(dotted)==157
 return {'status':'EXACT_NO_FINITE_VOLUME_COMPLETION_BY_WALLS_DOTTED_TO_FACET175','source':SOURCE,'source_sha256':sha(SOURCE),'facet_wall':175,'ambient_hyperbolic_dimension':20,'facet_hyperbolic_dimension':19,'facet_vector_space_dimension':20,'facet_norm':g.encode(qr),'projected_center_norm':g.encode(cq),'finite_neighbors':neighbors,'finite_neighbor_count':18,'dotted_walls':dotted,'dotted_count':157,'all175_projection_checks':rows,'homogeneous_equations_at_most':19,'positive_recession_dimension_at_least':1,'scope':'Every full-dimensional chamber obtained by retaining these176 walls and adding finitely many walls all strictly dotted to wall175 has infinite volume. A direct null-direction perturbation proof uses a strict point of the proposed completion, not necessarily the archived center. Added walls intersecting or parallel to wall175, wall replacements, and other constructions remain unrestricted. The numerical second orbit is not needed for this theorem.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print('neighbors',d['finite_neighbors'])
