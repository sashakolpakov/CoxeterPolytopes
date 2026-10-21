#!/usr/bin/env python3
"""Numerical complete-boundary cost of releasing old150. No exact promotion."""
import argparse,hashlib,json
import numpy as np
import fresh_h20_replacement_outer_patch_20261017 as S
O,L,ROOT=S.O,S.L,S.ROOT
TAG='fresh_h20_replacement150_boundary_20261017'
PAIR='results/fresh_h20_replacement_outer_patch_pairs_20261017.json'
BELT='results/fresh_h20_replacement_belt_step_20261017.json'
def geometry():
 D=O.P.P.load(S.INPUT);E=O.P.P.load(f'results/{S.TAG}.json');F=O.P.P.load(PAIR);g=D['geometry'];Y=np.array(g['residual_normals']);CT=np.array(g['T_products']);J=np.array(g['T_inverse']);q=np.array(g['native_norms']);m=np.array(g['residual_form']);T=D['support'];ret=D['stages'][-1]['retained_ids']+[173,174]
 for z in(E['candidates'][0],F['candidates'][77]):Y=np.r_[Y,[z['residual4']]];CT=np.r_[CT,[z['T_pairings']]];q=np.r_[q,2]
 return D,F,Y,CT,J,q,m,T,ret

def run():
 D,F,Y,CT,J,q,m,T,ret=geometry();B=O.P.P.load(BELT);rows=[z for z in B['at_most_three_bad_gates']if [u['wall']for u in z['failed_gates']]==[150]];assert len(rows)==18
 released=[j for j in ret if j!=150];delete=L.boundary(Y,m,released,T);out=[]
 for i,z in enumerate(rows):
  r=np.array(z['residual4']);dt=np.array(z['T_pairings']);co=-(r@m@Y.T+CT@J@dt)/np.sqrt(2*q);labels=[O.R.classify(float(x))for x in co];assert all(labels[j]not in('obtuse','angle_gap')for j in released)
  YY=np.r_[Y,[r]];b=L.boundary(YY,m,released+[175],T);rec=dict(index=i,recipe=z,retained_native_labels=[labels[j]for j in released],boundary=b);out.append(rec)
  print('CASE',i,z['contacts'],z['label'],b['poly_vertex_counts'],'open',len(b['open_edges']),'max',b['maximum_squared_radius'],flush=True)
 paths=[f'scripts/{TAG}.py',S.INPUT,f'results/{S.TAG}.json',PAIR,BELT,f'scripts/{O.TAG}.py',f'scripts/{L.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_OLD150_RELEASE_BOUNDARY_COST',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},base_retained_ids=ret,released_retained_ids=released,old150_norm=float(q[150]),old150_T_contacts=[[j,float(CT[150,k]),float(-CT[150,k]/np.sqrt(q[150]*q[j]))]for k,j in enumerate(T)if abs(CT[150,k])>1e-7],base_boundary=F['candidates'][77]['boundary'],deletion_only_boundary=delete,candidates=out,scope='Numerical complete target-H3 boundary after deleting old150 from the173-wall bounded two-cap branch, and adding each of18 previously recorded candidates whose only native Coxeter gate failure was150. Full retained native gates recounted and all positive-side vertices/ball-clipped edges recomputed. This is not an exact chamber certificate or an all-parameter search.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');z=a.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if z.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
