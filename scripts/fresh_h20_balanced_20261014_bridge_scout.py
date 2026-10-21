#!/usr/bin/env python3
"""Bounded product-face bridge scout; exact geometry data, numerical gates."""
import argparse,json,math,time
from itertools import product
import numpy as np
from fresh_h20_balanced_20261014 import ROOT,data,sha
TAG='fresh_h20_balanced_20261014_bridge_scout'
def reps(h):
 out={}
 for mask in range(1,1<<h):
  S=tuple(i for i in range(h)if mask>>i&1)
  orbit={tuple(sorted((j+eps*i)%h for i in S))for eps in(-1,1)for j in range(h)}
  canon=min(orbit)
  out[canon]=sorted(orbit)
 return sorted(out.items())
def allowed(co,cs):return (co>=1-1e-9)|np.any(abs(np.asarray(co)[...,None]-cs)<1e-9,axis=-1)
def run(golden):
 start=time.time();phi=(1+math.sqrt(5))/2;D=data((1,1,phi if golden else 1))
 factors=[reps(h)for h in[7,7,8]];fv=D['factor_vertices'];VV=np.array([v['point']for v in D['vertices']]);v2=np.sum(VV*VV,axis=1)
 cs=np.array([math.cos(math.pi/m)for m in range(2,13)]+[1.]);labels=list(range(2,13))+['parallel']
 Rold=math.sqrt(28);counts=dict(factor_orbit_representatives=list(map(len,factors)),product_representatives=0,seed_and_central=0,full_vertex_layer=0)
 accepted=[];allseed=[];seen=set()
 for ia,ib,ic in product(*[range(len(f))for f in factors]):
  if ia>ib:continue
  choice=[factors[k][a]for k,a in enumerate([ia,ib,ic])];subsets=[list(t[0])for t in choice]
  center=sum(np.mean(fv[k][s],axis=0)for k,s in enumerate(subsets));N=float(center@center);counts['product_representatives']+=1
  if N<=4+1e-10:continue
  gaps=D['offsets']-D['vectors']@center;assert min(gaps)>-1e-8
  positive=gaps[gaps>1e-8]
  candidate=[]
  # Every tested central label fixes a radius; old finite anchors add events
  # where the central angle can be tested independently.
  for m,c in zip(labels,cs):
   rr=-2*c+math.sqrt(N-4+4*c*c)
   if rr>1e-9:candidate.append((rr,{'central_label':m}))
  for a in positive:
   for m,c in zip(labels[1:],cs[1:]):
    rr=a/c;candidate.append((rr,{'seed_anchor_label':m,'seed_anchor_gap':float(a)}))
  radii=[]
  for rr,event in candidate:
   if any(abs(rr-r)<1e-8 for r in radii):continue
   radii.append(rr)
   central=(N-4-rr*rr)/(4*rr)
   if central>1+1e-9 or not allowed(central,cs):continue
   if not allowed(gaps/rr,cs).all():continue
   row=dict(factor_subsets=subsets,center=center.tolist(),center_square=N,radius=rr,central_cosine=central,event=event)
   allseed.append(row);counts['seed_and_central']+=1
   co=(N+v2-2*VV@center-rr*rr-Rold*Rold)/(2*rr*Rold)
   bad=np.flatnonzero(~allowed(co,cs))
   row['vertex_layer_conflicts']=bad.tolist();row['minimum_vertex_layer_cosine']=float(co.min())
   if len(bad):continue
   counts['full_vertex_layer']+=1;accepted.append(row)
 # Whole-bank relaxation, expanded through all diagram symmetries. No pair
 # compatibility or global finite-volume conclusion follows from this step.
 caps=[]
 for row in accepted:
  options=[dict(factors[k])[tuple(row['factor_subsets'][k])]for k in range(3)]
  for ss in product(*options):
   permutations=[ss]if ss[0]==ss[1]else[ss,(ss[1],ss[0],ss[2])]
   for subsets in permutations:
    center=sum(np.mean(fv[k][list(s)],axis=0)for k,s in enumerate(subsets))
    key=(tuple(np.round(center,9)),round(row['radius'],9))
    if key not in seen:seen.add(key);caps.append((center,row['radius'],[list(s)for s in subsets]))
 fractions=[.15,3/14,.25,.3,.4,.5,.6,.7,.75,.8,.9,1.]
 X=np.concatenate([a*VV for a in fractions]);norm=np.sum(X*X,axis=1)
 C=np.array([x[0]for x in caps])if caps else np.empty((0,19));R=np.array([x[1]for x in caps])
 # Include the central sphere and complete coherent sqrt28 vertex layer.
 minpower=np.minimum(norm-4,np.min(norm[:,None]+v2[None,:]-2*X@VV.T-28,axis=1))
 for i in range(0,len(C),256):
  CC=C[i:i+256];rr=R[i:i+256]
  pp=norm[:,None]+np.sum(CC*CC,axis=1)[None,:]-2*X@CC.T-rr[None,:]**2
  minpower=np.minimum(minpower,np.min(pp,axis=1))
 worst=int(np.argmax(minpower));idx=worst%392;frac=fractions[worst//392]
 paths=['scripts/fresh_h20_balanced_20261014.py',f'scripts/{TAG}.py','results/fresh_h20_intrinsic_cusp_benchmark_20261014.json']
 out=dict(status='BOUNDED_NUMERICAL_PRODUCT_FACE_CENTRAL_CONTACT_SCOUT',sources={p:sha(p)for p in paths},scales=D['scales'],labels_tested=labels,vertex_layer_radius_squared=28,counts=counts,seed_central_candidates=allseed,retained_layer_candidates=accepted,expanded_cap_count=len(caps),expansion_recipe='Apply each affine cycle dihedral group independently and swap the two A6 factors; remove identical (center,radius) rows.',probe_fraction_list=fractions,probe_count=len(X),uncovered_probes=int(np.sum(minpower>1e-8)),worst_probe=dict(vertex_index=idx,fraction=frac,minimum_ball_power=float(minpower[worst]),point=X[worst].tolist()),elapsed_seconds=time.time()-start,scope='Numerical bounded face-center radii from labels2..12 and parallel events, reduced by exact affine-cycle diagram symmetries. Initial core is centralR2 plus all392 vertexcapsR=sqrt28. Accepted bridge candidates pass seed/central/vertex pair proposals only; new-new compatibility and full volume are not certified. Probe gaps are numerical only. No all-label or unrestricted obstruction.')
 out['radial_probe_summary']=[dict(fraction=f,uncovered=int(np.sum(minpower[k*392:(k+1)*392]>1e-8)),worst_power=float(np.max(minpower[k*392:(k+1)*392])))for k,f in enumerate(fractions)]
 suffix='_golden'if golden else'';target=ROOT/'results'/f'{TAG}{suffix}.json';target.write_text(json.dumps(out,indent=2)+'\n')
 np.savez_compressed(ROOT/'data'/'enumerations'/f'{TAG}{suffix}.npz',centers=C,radii=R)
 print(counts,'expanded',len(caps),'radial',out['radial_probe_summary'])
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--golden',action='store_true');a=ap.parse_args();run(a.golden)
