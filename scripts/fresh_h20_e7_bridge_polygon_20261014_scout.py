#!/usr/bin/env python3
"""Independent numerical same-face compatibility diagnostic; not a proof."""
import argparse,json,math
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_bridge_polygon_20261014_scout'
SOURCE='results/fresh_h20_e7_bridge_face_20261014_scout.json'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text());b=json.loads((ROOT/BASE).read_text())
 for doc in[d,b]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 cs=d['candidates'];N=np.array([c['normal']for c in cs])
 Q=np.array([[float(Fraction(x))for x in row]for row in b['form']])
 G=N@Q@N.T;assert np.max(abs(np.diag(G)-2))<1e-5
 cutters={i for i,c in enumerate(cs)if c['escape_pairing']>1e-7}
 compatible=set();finite=[];counts=Counter();closest=[]
 for i in range(len(cs)):
  for j in range(i):
   if cs[i]['support']!=cs[j]['support']:continue
   counts['same_face_pairs']+=1;g=float(G[i,j])
   if g>1e-7:continue
   if g < -2-1e-7:kind='dotted'
   elif abs(g+2)<=1e-7:kind='parallel'
   elif abs(g)<=1e-7:kind='right'
   elif g>-1+1e-7:continue
   else:
    m=math.pi/math.acos(min(1.,-g/2));label=round(m)
    if abs(m-label)>2e-5:
     if i in cutters or j in cutters:closest.append(dict(candidates=[j,i],pairing=g,inferred_label=m,integer_gap=abs(m-label)))
     continue
    kind='pi'+str(label)
   compatible.add((j,i));counts[kind]+=1
   if kind!='dotted':finite.append(dict(candidates=[j,i],pairing=g,relation=kind))
 stored={tuple(r['candidates'])for r in d['compatible_pairs']if cs[r['candidates'][0]]['support']==cs[r['candidates'][1]]['support']}
 assert compatible==stored
 assert len(cs)==668 and len(cutters)==32 and len(compatible)==40265 and not finite
 paths=[SOURCE,BASE,f'scripts/{TAG}.py']
 return dict(status='NUMERICAL_SAME_FACE_POLYGON_DIAGNOSTIC_NO_FINITE_OR_IDEAL_CORNERS',
  sources={p:sha(p)for p in paths},candidate_normals=len(cs),escape_cutters=len(cutters),
  counts=dict(counts),same_face_compatible_pairs=len(compatible),finite_or_parallel_pairs=finite,
  closest_failed_finite_labels=sorted(closest,key=lambda r:r['integer_gap'])[:8],
  scope='Independent floating-point recomputation of pairings in the frozen bounded668-normal pool. Every accepted same-face pair is dotted, so this pool supplies no finite or ideal H2 corner and no finite-area polygon on either specified face. No exact pair-sign theorem, all-label exclusion, or obstruction to other caps is claimed.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'])
