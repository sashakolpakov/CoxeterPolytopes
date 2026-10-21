#!/usr/bin/env python3
"""Bounded transverse-swap diagnostic after the frozen43-wall construction."""
import argparse,json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_new_projection_20261015_scout'
BASE='results/fresh_h20_joint_route_20261015_exact.json'
BANK='results/fresh_h20_joint_route_20261015_two_anchor.json'
def run():
 d=json.loads((ROOT/BASE).read_text());bank=json.loads((ROOT/BANK).read_text())
 dec=lambda r:sum(float(Fraction(a))*np.sqrt(b)for a,b in zip(r,[1,2,3,6]))
 Q=np.array([[dec(x)for x in row]for row in d['form']]);N=np.array([[dec(x)for x in row]for row in d['normals']]);E=np.array(list(map(dec,d['strict_escape'])))
 norm=np.einsum('ij,jk,ik->i',N,Q,N);angles=np.cos(np.pi/np.arange(2,61))
 accepted=[];cutters=[];seen=set();counts=dict(tested=0,acute=0,compatible=0,cutting_recipes=0)
 for i,row in enumerate(bank['accepted']):
  a,b=row['support_contacts'];v=np.r_[row['source_projection'],-b/2,-a/2];g=N@Q@v;counts['tested']+=1
  if max(g)>1e-7:continue
  counts['acute']+=1;relations=[]
  for h in g/np.sqrt(2*norm):
   if h<-1-1e-7:label='dotted'
   elif abs(h+1)<1e-7:label='parallel'
   else:
    j=int(np.argmin(abs(angles+h)))
    if abs(angles[j]+h)>1e-7:break
    label=j+2
   relations.append(label)
  if len(relations)!=43:continue
  counts['compatible']+=1;accepted.append(v);e=float(v@Q@E)
  if e<=1e-7:continue
  counts['cutting_recipes']+=1;key=tuple(np.round(v,7))
  if key in seen:continue
  seen.add(key);cutters.append(dict(original_bank_index=i,source_rights=row['shared_source_support'],source_anchors=row['source_anchors'],source_anchor_labels=row['source_anchor_labels'],normal=v.tolist(),old_escape_pairing=e,relations=relations))
 V=np.array([c['normal']for c in cutters]);assert len(cutters)==6
 pairs=V@Q@V.T
 mutual=pairs[np.triu_indices(6,1)]
 pair_counts=dict(positive=int(sum(mutual>1e-7)),parallel=int(sum(abs(mutual+2)<1e-7)),dotted=int(sum(mutual<-2-1e-7)))
 assert pair_counts==dict(positive=6,parallel=3,dotted=6)
 E2=E.copy();E2[-2:]=[-.7,-.1];q2=float(E2@Q@E2);g2=np.array(accepted)@Q@E2
 paths=[BASE,BANK,f'scripts/{TAG}.py']
 return dict(status='NUMERICAL_BOUNDED_TRANSVERSE_SWAP_WITH_SIX_CUTTERS',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=counts,
  deduplicated_cutters=cutters,cutter_pair_matrix=pairs.tolist(),cutter_mutual_counts=pair_counts,renewed_escape_transverse=['-7/10','-1/10'],renewed_escape_norm=q2,
  renewed_escape_max_accepted_pairing=float(max(g2)),renewed_escape_remaining_cutters=int(sum(g2>1e-7)),
  scope='Only swaps of the1856 saved numerical recipes; finite labels2..60 with tolerance1e-7. No source projection census or all-label exclusion. The6cutters form two triples with positive products inside each triple and parallel/dotted products across triples. One simplest cutter is separately exactified without importing this scout. No finite-volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],'remaining cutters',d['renewed_escape_remaining_cutters'])
