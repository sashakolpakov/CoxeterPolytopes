#!/usr/bin/env python3
"""Reproducible numerical ideal-boundary scores in one fixed H3 chart.

These are low-discrepancy sample counts, not exact areas or volume proofs.
"""
import argparse,hashlib,json
import numpy as np
from scipy.stats import qmc
import fresh_h20_opposite_prism_cap_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_replacement_hexagon_scores_20261017'
FILES=['results/fresh_h20_prism_replacement_20261017.json','results/fresh_h20_replaced_old_edge_profiles_20261017.json','results/fresh_h20_h4_edge_continuation_20261017.json','results/fresh_h20_two_wall_replacement_cluster_20261017.json','results/fresh_h20_h4_hexagon_closure_20261017.json']
def run():
 G,Y,CT,J,q,V,C=P.geometry();rep,A,B,D,F=[json.loads((ROOT/p).read_text())for p in FILES];seed=A['candidates'][83];tri=B['candidates'][0];partner=next(z for z in D['candidates']if z['index']==62)['source_event']
 YY=np.r_[Y,[rep['residual4'][1],seed['residual4'],tri['residual4'],partner['residual4']]]
 u=qmc.Sobol(2,scramble=False).random_base2(19);z=2*u[:,0]-1;theta=2*np.pi*u[:,1];points=np.c_[np.sqrt(1-z*z)*np.cos(theta),np.sqrt(1-z*z)*np.sin(theta),z]
 stages=[('original168',Y),('one_wall_replacement168',np.r_[Y[:158],Y[159:],[rep['residual4'][1]]]),('nonright_T_seed168',YY[:170][[j for j in range(170)if j not in(31,158)]]),('four_new_roots170',YY[[j for j in range(172)if j not in(31,158)]])]
 for f in F['candidates']:stages.append((f'hexagon_final_label{f["label"]}',np.r_[YY[[j for j in range(172)if j not in(31,158)]],[f['residual4']]]))
 scores=[]
 for name,N in stages:
  # Zero projected normals are support equalities. A future timelike
  # projected normal is automatically nonpositive on the ideal sphere.
  active=[(j,n)for j,n in enumerate(N)if np.linalg.norm(n[1:])>1e-8 and np.linalg.norm(n[1:])-n[0]>1e-10];mask=np.ones(len(points),bool)
  for j,n in active:mask&=points@n[1:]<=n[0]
  scores.append(dict(stage=name,actual_wall_count=len(N),tested_spacelike_halfspaces=len(active),surviving_ideal_samples=int(sum(mask)),sample_count=len(points),fraction=float(sum(mask)/len(points)),mask_sha256=hashlib.sha256(mask.tobytes()).hexdigest()))
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py']+FILES
 return dict(status='REPRODUCIBLE_NUMERICAL_FIXED_CHART_IDEAL_SCORES',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},sample_rule='Unsmoothed scipy Sobol(2,scramble=False),2^19 points mapped to the uniform unit2-sphere by z=2u0-1 and longitude2pi*u1.',scores=scores,scope='One deterministic numerical angular-boundary score in the same fixed target-H3 Lorentz chart. No certified area, probability error bar, hyperbolic volume, exact negative or continuous-family completeness claim. All three final labels are alternatives, not simultaneous additions.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],[(z['stage'],z['surviving_ideal_samples'])for z in d['scores']],flush=True)
