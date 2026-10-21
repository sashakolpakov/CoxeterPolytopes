#!/usr/bin/env python3
"""Bounded intrinsic endpoint caps and global escape probes; floating proposals only."""
import json,math
from pathlib import Path
from itertools import product
import numpy as np
from fresh_h20_intrinsic_caps_20261014 import ROOT,TAG,data,sha,compatible_ratios
PROFILE=f'results/{TAG}_endpoint.json'

def main():
 s=(1+math.sqrt(5))/2;D=data(s);P=json.loads((ROOT/PROFILE).read_text());labels=list(range(2,13))+[16,20,24,30];cosines=[math.cos(math.pi/m)for m in labels]+[1.]
 types=[]
 for typ in P['types']:
  profiles=[]
  for p in typ['profiles']:
   profiles.append({'id':p['id'],'cos':p['cosine_vector_approx'],'endpoint':p['endpoint'],'contact_labels':p['contacts'],'B':[v['subtract_approx']for v in p['omitted_contact_ratios']]})
  for tail in typ['dihedral_tails']:
   # Explicitly bounded representatives of the otherwise unbounded tail.
   node=tail.get('isolated_wall',tail.get('wall'))
   if node is None:raise ValueError(tail)
   for m in labels:
    if m<7:continue
    cc=np.zeros(19);cc[typ['edge_support'].index(node)]=math.cos(math.pi/m)
    marks={str(w):v for w,v in zip(P['seed_order'],P['positive_affine_marks'])};bb=[]
    for f,o in zip(D['factors'],typ['omitted']):
     bb.append(sum(float(marks[str(w)])*cc[i]for i,w in enumerate(typ['edge_support'])if w in f)/float(marks[str(o)]))
    profiles.append({'id':f'tail_{m}','cos':cc.tolist(),'endpoint':'finite','contact_labels':[[node,m]],'B':bb})
  A=np.array([float(v['numerator_over_radius'])*(s if v['numerator_has_triangle_inradius_s']else 1.)for v in typ['profiles'][0]['omitted_contact_ratios']])
  rows=[]
  for p in profiles:
   B=np.array(p['B']);radii=[]
   for a,b,c in zip(np.repeat(A,len(cosines)),np.repeat(B,len(cosines)),cosines*3):
    if b+c<1e-10:continue
    R=a/(b+c)
    if compatible_ratios(A/R-B,cosines)and not any(abs(R-z)<1e-8 for z in radii):radii.append(R)
   if radii:rows.append((p,sorted(radii)))
  types.append(rows)
 lookup={tuple(v['omitted_walls']):v for v in D['vertices']};wallindex={w:i for i,w in enumerate(D['walls'])}
 centers=[];radii=[];recipes=[];seen=set()
 for ti,tr in enumerate(P['literal_edge_transports']):
  ids=[wallindex[w]for w in tr['edge_support_in_profile_order']];V=D['vectors'][ids];z=lookup[tuple(tr['omitted'])]['point']
  for p,rs in types[tr['type']]:
   y=np.linalg.solve(V,np.array(p['cos']));assert y@y<=1+1e-8
   for R in rs:
    w=z-R*y;d=D['offsets']-D['vectors']@w
    assert min(d)>-1e-7 and compatible_ratios(d/R,cosines,1e-7)
    key=tuple(np.rint(np.append(w,R)*1e7).astype(np.int64))
    if key in seen:continue
    seen.add(key);centers.append(w);radii.append(R);recipes.append([ti,p['id'],R,p['endpoint']])
 C=np.array(centers);rr=np.array(radii);print('bank',len(rr),'bytype profile/radius counts',[(len(t),sum(len(rs)for p,rs in t))for t in types],flush=True)
 # Mandatory golden central cap, with radius2, is checked separately from the
 # seed. Endpoint caps not compatible with it are unavailable to this branch.
 c2=np.sum(C*C,axis=1);co=(c2-rr*rr-4)/(4*rr);compat=co>=1-1e-8
 for c in cosines:compat|=abs(co-c)<1e-8
 print('central compatible',int(compat.sum()),flush=True)
 points=[np.zeros(19)];names=['incenter'];VV=np.array([v['point']for v in D['vertices']])
 for scale in[.25,.5,.75,1.]:
  points.extend(scale*VV);names.extend([f'vertex_{i}_scaled_{scale}'for i in range(len(VV))])
 X=np.array(points);margins=[];central_margins=[];best=[];bestc=[]
 for p in X:
  values=np.sum((C-p)**2,axis=1)-rr*rr
  j=int(np.argmin(values));margins.append(float(values[j]));best.append(j)
  vc=np.where(compat,values,np.inf);j=int(np.argmin(vc));central_margins.append(min(float(vc[j]),float(p@p-4)));bestc.append(j)
 worst=int(np.argmax(central_margins));bestall=int(np.argmax(margins))
 np.savez_compressed(ROOT/f'data/enumerations/{TAG}_endpoint_bank.npz',centers=C,radii=rr,central_compatible=compat)
 out={'status':'NUMERICAL_INTRINSIC_ENDPOINT_BANK_AND_ESCAPE_PROBES','sources':{p:sha(p)for p in[PROFILE,f'scripts/{TAG}.py',f'scripts/{TAG}_completion_scout.py']},'triangle_inradius':s,'finite_labels_tested':labels,'tail_scope':'Only the listed finite labels; arbitrary-I2 tail remains unbounded.','bank_size':len(rr),'central_compatible_caps':int(compat.sum()),'profile_radius_counts':[sum(len(rs)for p,rs in t)for t in types],'probe_count':len(X),'bank_uncovered_probes':int(sum(v>1e-7 for v in margins)),'central_branch_uncovered_probes':int(sum(v>1e-7 for v in central_margins)),'worst_full_bank_probe':{'name':names[bestall],'point':X[bestall].tolist(),'minimum_margin':margins[bestall],'closest_cap':best[bestall]},'worst_central_branch_probe':{'name':names[worst],'point':X[worst].tolist(),'minimum_margin':central_margins[worst],'closest_cap':bestc[worst]},'recipes':recipes,'scope':'Numerical bounded cap pool. Compatibility with the seed and central cap is necessary only; endpoint caps are not tested mutually. Positive escape margins rule out the bank numerically, not arbitrary caps. No finite-volume conclusion.'}
 (ROOT/f'results/{TAG}_completion_scout.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items()if k not in['recipes','sources']},flush=True)
if __name__=='__main__':main()
