#!/usr/bin/env python3
"""Solve finite contacts from endpoint conics to the 310 three-side caps."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h20_two_anchor_families_20261016 import real_roots,point
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_seed_edge_families_20261016_pairs'
BASE='results/fresh_h20_seed_edge_families_20261016.json';SEED='results/fresh_h20_open_edge_closure_20261016.json'
def run():
 D=json.loads((ROOT/BASE).read_text());S=json.loads((ROOT/SEED).read_text());G=geometry158();M=G['metric'];H=np.array(S['seed']['inherited_H'])
 X=np.r_[G['X'],S['seed']['roots']];q=np.r_[G['norms'],[2.]*5];Q=np.array([x['root']for x in S['candidates']]);Z=np.array(D['exterior_points'])
 labels=[2,3,4,5,6,8,10];stats=Counter();pool=[];seen={};maps=[]
 for f in D['families']:
  num=np.array(f['numerator']);den=np.array(f['denominator']);ip=num@M@Q.T
  for j,r0 in enumerate(Q):
   for m in labels:
    c=0 if m==2 else np.cos(np.pi/m)
    for t in real_roots(ip[:,j]+2*c*den):
     stats['quadratic_roots']+=1
     if not any((a is None or t>a+1e-9*max(1,abs(a)))and(b is None or t<b-1e-9*max(1,abs(b)))for a,b in f['intervals']):continue
     r=point(f,t)
     if abs(r@M@r-2)>2e-4 or r@M@H>=0 or r[0]<=0:continue
     if abs(r@M@r0+2*c)>1e-5:continue
     if not gate_mask((-(r@M@X.T)/np.sqrt(2*q))[None,:])[0]:continue
     ar=-r[1:]/r[0];cuts=np.where(Z@ar>1+1e-7)[0].tolist()
     if not cuts:continue
     key=tuple(np.round(r,7));stats['accepted_relations']+=1
     if key not in seen:
      seen[key]=len(pool);pool.append(dict(id=len(pool),family=f['id'],parameter=t,root=r.tolist(),constraint=ar.tolist(),cuts=cuts))
     maps.append(dict(candidate=seen[key],partner=j,mutual_label=m))
 print('SOLVED',len(pool),len(maps),flush=True)
 src=[f'scripts/{TAG}.py',BASE,SEED,'scripts/fresh_h20_two_anchor_families_20261016.py','scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py']
 return dict(status='PASS_BOUNDED_ENDPOINT_FAMILY_PARTNER_SOLVER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},stats=dict(stats),mutual_labels=labels,candidates=pool,partner_maps=maps,scope='Numerical quadratic finite-contact recovery for151 saved endpoint families and310 fixed caps, stated seven labels. Every accepted root retains163 old gates. No complete arbitrary-label/family classification or exact construction.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['candidates']),flush=True)
