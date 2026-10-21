#!/usr/bin/env python3
"""Replay numerical selected-cap/probe geometry; not exact certification."""
import argparse,json
from pathlib import Path
import numpy as np
from fresh_h20_intrinsic_caps_20261014_quarter_face_scout import build,allowed,ROOT,TAG,sha

SOURCE=f'results/{TAG}_quarter_face_scout.json'
OUT=ROOT/f'results/{TAG}_quarter_face_geometry_scout.json'
def run():
 data=json.loads((ROOT/SOURCE).read_text())
 for p,h in data['sources'].items():assert sha(p)==h,p
 D,V,X,C,cover,dist,recipes,counts=build('contiguous');S=data['selected'];assert len(S)==60
 assert [recipes[i]for i in S]==data['selected_recipes']
 C=C[S];oldC=np.concatenate([np.zeros((1,19)),V,V/2,C]);oldR=np.concatenate([[2],np.ones(504),np.full(len(C),1.5)])
 P=V[0]/8;values=np.sum((P-oldC)**2,axis=1)-oldR**2
 qdist=dist[S]-2.25;dd=np.maximum(np.sum(C*C,axis=1)[:,None]+np.sum(C*C,axis=1)[None,:]-2*C@C.T,0)
 cp=(dd/4.5-1)[np.triu_indices(60,1)];assert allowed(cp,(1+5**.5)/2).all()
 assert values.min()>.8 and (D['offsets']-D['vectors']@P).min()>.87
 assert qdist.min(axis=0).max()<-.09
 paths=[SOURCE,f'scripts/{TAG}_quarter_face_scout.py',f'scripts/{TAG}_quarter_face_geometry_scout.py']
 return dict(status='NUMERICAL_587_WALL_PROBE_CHECKPOINT',sources={p:sha(p)for p in paths},total_walls=587,new_caps=60,
 new_pair_counts={'right':int((abs(cp)<1e-8).sum()),'strictly_dotted':int((cp>1+1e-8).sum())},
 maximum_quarter_probe_minimum_power=float(qdist.min(axis=0).max()),remaining_escape=dict(recipe='first old product vertex divided by8',old_vertex_index=0,old_vertex_omitted_local=D['vertices'][0]['omitted'],point=P.tolist(),minimum_spherical_power=float(values.min()),closest_spherical_cap=int(values.argmin()),minimum_seed_distance=float((D['offsets']-D['vectors']@P).min())),
 scope='Floating strict ideal-boundary gap and positive probe-cover diagnostics only. Exactification required; not a finite-volume claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['new_pair_counts'],out['remaining_escape']['minimum_spherical_power'])
