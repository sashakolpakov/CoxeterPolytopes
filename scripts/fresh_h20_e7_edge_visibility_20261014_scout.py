#!/usr/bin/env python3
"""Numerical visibility diagnostic on the116 currently verified finite vertices."""
import json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_edge_visibility_20261014_scout'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json';OLD='results/fresh_h20_e7_finite_edgewalk_20261014_exact.json';WALK='results/fresh_h20_e7_repaired_cusp_walk_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 b,o,w=[json.loads((ROOT/p).read_text())for p in[BASE,OLD,WALK]]
 for doc in[b,o,w]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 rr=np.sqrt(b['field']['base_square_root_basis']);dd=sum(float(F(x))*r for x,r in zip(b['field']['discriminant'],rr))**.5
 def bv(x):return sum(float(F(q))*r for q,r in zip(x,rr))
 def ev(x):return bv(x[0])+dd*bv(x[1])
 def vec(x):return np.array([ev(a)for a in x])
 def vecb(x):return np.array([bv(a)for a in x])
 Q=np.array([[float(F(x))for x in row]for row in b['form']]);N=np.array([vec(v)for v in b['normals']]);NQ=N@Q;H=vec(b['strict_center']);E=vec(b['strict_escape']['point']);norms=np.array(list(map(float,b['norms'])));gram=NQ@N.T
 vertices=[dict(incident=r['incident'],point=vecb(r['point']),source=['old90',i])for i,r in enumerate(o['preserved_old_vertices']+o['new_cap_finite_vertices'])]
 vertices +=[dict(incident=r['incident'],point=vecb(r['point']),source=['old95endpoint',i])for i,r in enumerate(b['new_edge_endpoints'])]
 vertices +=[dict(incident=r['incident'],point=vec(r['point']),source=['cuspwalk',i])for i,r in enumerate(w['finite_vertices'])if r['known_vertex']is None]
 assert len(vertices)==116 and len({tuple(sorted(r['incident']))for r in vertices})==116
 counts=Counter();rows=[];peak=0.
 for vi,row in enumerate(vertices):
  ids=row['incident'];V=row['point'];v=V/np.sqrt(-V@Q@V);GI=np.linalg.inv(gram[np.ix_(ids,ids)]);dual=GI@N[ids]
  assert np.max(abs(dual@Q@N[ids].T-np.eye(20)))<2e-6
  for j,oldwall in enumerate(ids):
   support=[i for i in ids if i!=oldwall];other=[i for i in range(95)if i not in support]
   z=v-dual[j]/np.sqrt(GI[j,j]);pairs=NQ@z;err=max(abs(pairs[support]));peak=max(peak,float(err));assert err<2e-6
   mx=float(max(pairs[other]));cut=float(E@Q@z);zn=float(z@Q@z)
   if mx>2e-7:kind='another_wall_before_null_numeric'
   elif mx< -2e-7:kind='strict_open_null_numeric'
   else:kind='null_contact_numeric'
   counts[kind]+=1
   if kind=='strict_open_null_numeric':
    counts['strict_open_with_positive_E_numeric']+=cut>1e-7
    rows.append(dict(vertex=vi,source=row['source'],released=oldwall,support=support,null_point=z.tolist(),escape_pairing=cut,max_nonshared_pairing=mx,null_norm_error=zn))
 rows.sort(key=lambda r:r['escape_pairing'],reverse=True)
 paths=[BASE,OLD,WALK,f'scripts/{TAG}.py']
 return dict(status='NUMERICAL_VISIBILITY_SCREEN_NOT_AN_EXACT_CENSUS',sources={p:sha(p)for p in paths},known_vertices=116,directed_germs=2320,counts=dict(counts),maximum_incident_roundoff=peak,best_strict_open_directions=rows[:40],scope='Numerical diagnostic of20directions at each116verified finitevertex, with repeated edges retained. No exact signs or full95vertex/edge completeness claimed. A positive strict feasible null direction would motivate an exact far-cap construction; negative values do not exclude nonright or finite-position cutters.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n');print(d['counts']);print(d['best_strict_open_directions'][:1])
