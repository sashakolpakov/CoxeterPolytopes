import json,numpy as np
from pathlib import Path
from fractions import Fraction as F
from math import lcm
P=Path(__file__).resolve().parents[1];d=json.loads((P/'results/fresh_h22_gap_face_family_20261013.json').read_text());s=json.loads((P/'data/enumerations/fresh_h22_new_cusp_completion_20261013_scout.json').read_text());raw=json.loads((P/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64)
gens=[np.array(g['matrix'],dtype=np.int64)for g in json.loads((P/'data/enumerations/fresh_h20_global_interior_sections_20261013.json').read_text())['ambient_symmetry_generators']]
bases=[[F(x)for x in s['caps'][i]['base']]for i in(7,9)];ds=[lcm(*(x.denominator for x in b))for b in bases];V,W=[np.array([int(x*n)for x in b],dtype=np.int64)for b,n in zip(bases,ds)]
delta=np.array(d['six_old_projection'])+np.array(d['new_cusp']['six_partner']);rows=[delta];vs=[V];ws=[W];lookup={delta.tobytes():0};hol={};vn=int(V@Q@V);wn=int(W@Q@W)
for i, row in enumerate(rows):
 for k,g in enumerate(gens):
  y=g@row;v=g@vs[i];w=g@ws[i];key=y.tobytes()
  if key not in lookup:
   lookup[key]=len(rows);rows.append(y);vs.append(v);ws.append(w)
  else:
   j=lookup[key]
   if not(np.array_equal(v,vs[j])and np.array_equal(w,ws[j])):
    a=F(int(v@Q@vs[j]),vn);b=F(int(v@Q@ws[j]),wn);c=F(int(w@Q@vs[j]),vn);e=F(int(w@Q@ws[j]),wn)
    assert all(a*int(x)+b*int(y)==int(z)for x,y,z in zip(vs[j],ws[j],v))
    assert all(c*int(x)+e*int(y)==int(z)for x,y,z in zip(vs[j],ws[j],w))
    t=tuple(map(str,(a,b,c,e)))
    if t not in hol:hol[t]={'i':i,'j':j,'generator':k};print('holonomy',t,hol[t],flush=True)
 if i%10000==0:print('orbit',i,len(rows),len(hol),flush=True)
 assert len(rows)<=241920
print('FINAL',len(rows),hol,flush=True)
np.savez_compressed(P/'data/enumerations/fresh_h22_new_cusp_completion_20261013_frames.npz',six_cusp_old=np.array(rows),first_tangent_numerator=np.array(vs),second_tangent_numerator=np.array(ws))
(P/'data/enumerations/fresh_h22_new_cusp_completion_20261013_holonomy.json').write_text(json.dumps({'cusp_orbit':len(rows),'tangent_denominators':ds,'tangent_numerator_norms':[vn,wn],'transition_matrices':[{'matrix':list(t),'edge':v}for t,v in hol.items()]},indent=2)+'\n')
