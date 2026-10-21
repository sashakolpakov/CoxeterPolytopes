#!/usr/bin/env python3
"""Bounded complete sampled clique cover over the 310 endpoint candidates."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fresh_h20_continuous_sample_joint_cover_20261016 import Cover
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_open_edge_closure_20261016_joint'
BASE='results/fresh_h20_open_edge_closure_20261016.json'
def run():
 D=json.loads((ROOT/BASE).read_text());G=geometry158();M=G['metric'];base,V,A,b,ids=face_data(G)
 seed=np.array(D['seed']['roots']);sa=-seed[:,1:]/seed[:,0,None];A=np.r_[A,sa@V];b=np.r_[b,1-sa@base];ids+=list(range(158,163))
 R=np.array([r['root']for r in D['candidates']]);AR=-R[:,1:]/R[:,0,None];search=Cover(R,M)
 Z=[np.array(v['z'])for v in D['seed']['face']['vertices']if v['radius2']>1+1e-6]
 center=np.mean([v['z']for v in D['seed']['face']['vertices']if v['radius2']<1-1e-6],axis=0)
 steps=[];final=None;reason='iteration_guard'
 for k in range(50):
  cov=np.array(Z)@AR.T>1+1e-7
  if not cov.any(axis=1).all():reason='uncovered_sample';steps.append(dict(iteration=k,uncovered=np.where(~cov.any(axis=1))[0].tolist()));break
  selected,meta=search.search(cov)
  if selected is None:reason='infeasible_finite_bank_cover'if meta['complete']else'node_guard';steps.append(dict(iteration=k,search=meta));break
  co=-R[selected]@M@R[selected].T/2;np.fill_diagonal(co,0);assert gate_mask(co).all()
  aa=AR[selected];face=poly(base,V,np.r_[A,aa@V],np.r_[b,1-aa@base],ids+list(range(163,163+len(selected))))
  row=dict(iteration=k,selected=selected,search=meta,face=face);steps.append(row)
  print('STEP',k,'caps',len(selected),'max',face['maximum_squared_radius'],flush=True)
  if not face['exterior']:reason='compact_candidate'if not face['near_sphere']else'near_sphere';final=row;break
  worst=max(face['vertices'],key=lambda x:x['radius2']);v=np.array(worst['z'])-center
  t=next(float(t)for t in np.roots([v@v,2*center@v,center@center-1])if np.isreal(t)and 0<float(t)<1)
  z=center+t*v;assert abs(z@z-1)<1e-7 and max(aa@z-1)<1e-6;Z.append(z)
 src=[f'scripts/{TAG}.py',BASE,'scripts/fresh_h20_continuous_sample_joint_cover_20261016.py','scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py']
 return dict(status='PASS_BOUNDED_OPEN_EDGE_CANDIDATE_JOINT_COVER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},candidate_count=len(R),termination=reason,steps=steps,samples=[z.tolist()for z in Z],final=final,scope='Numerical finite310-candidate bank only; inherited complete clique search unless explicit2million-node guard, and50iteration guard. No continuous-family or arbitrary-label conclusion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['termination'],len(d['steps']),flush=True)
