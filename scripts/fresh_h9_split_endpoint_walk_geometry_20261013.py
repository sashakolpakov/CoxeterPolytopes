#!/usr/bin/env python3
"""Numerical coherent replacement geometry for the spherical-edge walk."""
import hashlib,json,os
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection,QhullError
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_endpoint_walk_20261013_scout.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_endpoint_walk_geometry_20261013_scout.json'

def bits(m):
 while m:
  b=m&-m;yield b.bit_length()-1;m-=b

def run():
 Path('/tmp/fresh_h9_split_endpoint_walk_geometry_20261013.pid').write_text(str(os.getpid()))
 raw=SOURCE.read_bytes();src=json.loads(raw);old=src['source_snapshot'];base=np.array(old['physical_normals']);oldn=len(base);W=np.vstack([base,np.array([r['normal_physical']for r in src['candidates']])]);eta=np.diag([-1]+[1]*9);pair=W@eta@W.T
 C=pair< -2-1e-7
 for m in range(2,121):C|=abs(pair+2*np.cos(np.pi/m))<1e-7
 np.fill_diagonal(C,False);adj=[sum(1<<int(i)for i in np.flatnonzero(row))for row in C]
 mandatory=list(range(16))+sorted({i for v in old['finite_vertices']for i in v['active']if i>=16});allbits=(1<<len(W))-1;baseallowed=allbits
 for i in mandatory:baseallowed&=adj[i]
 rng=np.random.default_rng(61013023);seen=set();rows=[];directions=rng.normal(size=(20,9))
 def pick(start):
  chosen=mandatory[:]+start;rem=baseallowed
  for i in start:rem&=adj[i]
  while rem:
   ids=list(bits(rem));degree=np.array([(adj[i]&rem).bit_count()for i in ids]);new=np.array([i>=oldn for i in ids]);scores=degree+rng.exponential(30,size=len(ids))+new*80
   j=ids[int(np.argmax(scores))];chosen.append(j);rem&=adj[j]
  return sorted(chosen)
 starts=[[],[oldn+788],[oldn+777]]
 for trial in range(160):
  start=starts[trial%len(starts)] if trial<80 else[oldn+int(rng.integers(len(src['candidates'])))];chosen=pick(start);key=tuple(chosen)
  if key in seen:continue
  seen.add(key);norms=[];H=W[chosen]
  for d in directions:
   ans=linprog(-d,A_ub=H[:,1:],b_ub=H[:,0],bounds=[(None,None)]*9,method='highs');assert ans.success;norms.append(float(ans.x@ans.x))
  rows.append({'trial':trial,'selected':chosen,'new_cap_count':sum(i>=oldn for i in chosen),'facets':len(chosen),'sampled_max_square':max(norms),'sampled_mean_square':float(np.mean(norms))})
 ranked=sorted(rows,key=lambda r:(r['sampled_max_square'],r['sampled_mean_square']))[:8]
 for row in ranked:
  chosen=row['selected'];H=W[chosen]
  try:h=HalfspaceIntersection(np.c_[H[:,1:],-H[:,0]],np.zeros(9))
  except QhullError:h=HalfspaceIntersection(np.c_[H[:,1:],-H[:,0]],np.zeros(9),qhull_options='Qs')
  q=np.sum(h.intersections*h.intersections,axis=1);j=int(np.argmax(q));ff=np.flatnonzero(q<1-1e-7);ii=np.flatnonzero(abs(q-1)<=1e-7)
  row.update({'vertex_count':len(q),'finite_vertices':len(ff),'near_ideal_vertices':len(ii),'outside_vertices':sum(q>1+1e-7).item(),'max_square':float(max(q)),'leading_escape':{'point':h.intersections[j].tolist(),'active':[chosen[i]for i in h.dual_facets[j]]},'finite_vertex_data':[{'point':h.intersections[i].tolist(),'active':[chosen[k]for k in h.dual_facets[i]]}for i in ff]})
  print({k:v for k,v in row.items()if k not in['selected','leading_escape','finite_vertex_data']},flush=True)
  OUT.write_text(json.dumps({'status':'NUMERICAL_REPLACEMENT_GEOMETRY','source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(raw).hexdigest(),'scope':'Numerical finite-angle graph and actual affine Klein vertices for selected coherent replacements. No exact pair, all-affine, or compactness certificate.','mandatory_wall_ids':mandatory,'old_wall_count':oldn,'scout_clique_count':len(rows),'screened':rows,'geometry_cases':ranked},indent=2)+'\n')
if __name__=='__main__':run()
