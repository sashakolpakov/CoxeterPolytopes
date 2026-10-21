#!/usr/bin/env python3
"""Explicit uncapped ideal directions in the13 new H22 facet representatives.

Local arithmetic is separate Fraction computation. The fact that every
nonlocal wall is strictly dotted is a declared primary-certificate dependency.
"""
import argparse,json
from fresh_h22_transported_cusp_completion_20261013_local_independent import (
 P,F,read,sha,dot,mv,RAW,PROPOSAL,SOURCE,SHIFTS)
PREFIX='fresh_h22_transported_cusp_completion_20261013'
LOCAL='results/'+PREFIX+'_local_independent.json';OUT='results/'+PREFIX+'_facet_obstruction.json'
def run():
 source=read(SOURCE);audit=read(LOCAL)
 assert source['status']=='PASS_EXACT_COMPLETE_TRANSPORTED_H22_CUSP_REPAIR'
 assert audit['status']=='PASS_INDEPENDENT_STDLIB_FRACTION_TRANSPORTED_H22_LOCAL_LINK'
 assert all(sha(p)==h for p,h in source['sources'].items())and all(sha(p)==h for p,h in audit['sources'].items())
 raw=read(RAW)['coxeter_input'];prop=read(PROPOSAL)
 Q=[[-F(x)for x in row]+[F(0)]for row in raw['basis_gram']]+[[F(0)]*22+[F(11)]]
 W=[[F(x)for x in row]+[F(0)]for row in raw['wall_coordinates_as_columns_transposed']]
 from functools import reduce
 from math import gcd
 h=[sum(int(r[j])for r in W)for j in range(23)];g=reduce(gcd,h);h=[F(x,g)for x in h];hc=mv(Q,h)
 delta=list(map(F,prop['target']['delta']));dc=mv(Q,delta);assert dot(delta,dc)==0 and dot(delta,hc)<0
 oldcap=list(map(F,prop['target']['cap']));incident=prop['target']['incident'];bases=[list(map(F,r))for r in prop['bases']]
 caps=[[a+t*b for a,b in zip(r,delta,strict=True)]for r,t in zip(bases,SHIFTS,strict=True)];local=[W[i]for i in incident]+[oldcap]+caps
 records=[]
 for k,r in enumerate(caps):
  rc=mv(Q,r);q=dot(r,rc);hr=dot(h,rc);assert q>0 and hr<0 and dot(r,dc)==0
  x0=[a-hr*b/q for a,b in zip(h,r,strict=True)];xq=dot(x0,mv(Q,x0));xd=dot(x0,dc);assert xq<0 and xd<0
  coefficient=-xq/(2*xd);assert coefficient<0
  eta=[a+coefficient*b for a,b in zip(x0,delta,strict=True)];ec=mv(Q,eta)
  assert dot(eta,ec)==0 and dot(eta,hc)<0 and dot(eta,dc)<0 and dot(eta,rc)==0
  pairs=[dot(w,ec)for w in local];assert [i for i,p in enumerate(pairs)if p==0]==[18+k]and all(p<=0 for p in pairs)
  oldpairs=[dot(w,ec)for w in W];assert max(oldpairs)<0
  records.append({'profile':k,'normal':list(map(str,r)),'normal_norm':str(q),'projected_center':list(map(str,x0)),'projected_center_norm':str(xq),'projected_center_delta_pairing':str(xd),'eta_delta_coefficient':str(coefficient),'uncapped_null_ray':list(map(str,eta)),'null_ray_center_pairing':str(dot(eta,hc)),'null_ray_delta_pairing':str(dot(eta,dc)),'local31_pairings':list(map(str,pairs)),'old210_max_pairing':str(max(oldpairs))})
 return {'status':'PASS_EXACT_THIRTEEN_UNCAPPED_FACET_DIRECTIONS','sources':{p:sha(p)for p in(SOURCE,LOCAL,RAW,PROPOSAL)},'field_coordinate':'last coordinate is coefficient ofsqrt11;formdiag(Q,11)','representatives':records,'preferred_target_profile':7,
 'nonlocal_geometry':'Every remaining retained wall has strictly dotted relation to the chosen newwall r, by the primary fullorbit certificate. Its projection to r-perp is future timelike (orientation follows from the strict projected center), hence pairs strictly negatively with each future null eta. Together with the explicit local checks, eta is strict on every other global wall.',
 'escape_argument':'eta-epsilon*delta stays in r-perp and has square -2*epsilon*(eta,delta)>0. All local inequalities remain strict because their normals annihilate delta; every nonlocal inequality remains strict for sufficiently small epsilon since the finite list was strict at eta. A further small perturbation toward the common strict center makes the r inequality strict while preserving positive norm. Thus the chamber is infinite-volume independently of the remote rational escape.',
 'completion_scope':'Retain r and the existing acute-compatible facets. Any finite additions for which every hyperplane meeting r also contains delta leave the intrinsic facet as a cone with no transverse closing wall, and cannot produce finite volume. Dotted or parallel hyperplanes do not bound the interior of the hyperbolic facet. If necessary choose a generic projected interior point to avoid the finitely many added parallel ideal tangencies. Hence a finite-volume completion must include a new wall meeting r at finite hyperbolic points and not containing delta, or must replace r. This is not a no-go for arbitrary additions or changed incidence patterns.',
 'no_finite_vertices_on_new_facets':'Every normal at a finite vertex containing r would annihilate delta, but delta-perp has positive index21, whereas an H22 finite vertex needs22 independent positive directions. Therefore none of these newfacets contains a finite vertex in this fixedchamber.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
