#!/usr/bin/env python3
"""Numerical rescore of archived H9 nonpolar banks against all49 walls/E."""
import argparse,gzip,json
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_joint_regions_20261016_next_polar as P
import fresh_nonarithmetic_ridge_pair_20261016_polar as VP
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_ridge_pair_20261016_rescore'
BASIS=[0,1,2,3,4,5,6,7,8,11]
TAGS=['fresh_h9_after_third_vertex_free_scout_20261015','fresh_h9_after_third_endpoint_scout_20261015','fresh_h9_octic_base_scout_20261015','fresh_h9_octic_germ_scout_20261015','fresh_h9_nine_contact_line_scout_20261015','fresh_h9_nine_contact_germ_scout_20261015','fresh_h9_joint_pair_20261015_scout_escape','fresh_h9_joint_pair_20261015_scout_escape_right','fresh_h9_joint_corner_20261015_scout']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def smooth(x):return float(format(float(x),'.12g'))
def run():
 Q,N,H,E=VP.geometry()
 M=mp.matrix([[x for x in Q*N[i]]for i in BASIS]);inv=M**-1;assert mp.norm(M*inv-mp.eye(10))<mp.mpf('1e-60')
 transform=np.array([[float(x)for x in n.T*Q*inv]for n in N+[E,H]])
 paths=[f'scripts/{TAG}.py','results/fresh_nonarithmetic_common_escape_20261016_cap.json',VP.CO,VP.CL,f'scripts/{VP.TAG}.py',f'scripts/{P.TAG}.py'];records=[]
 for tag in TAGS:
  meta_path=f'results/{tag}.json';bp=f'data/enumerations/{tag}_bank.json.gz'
  meta=json.loads((ROOT/meta_path).read_text());saved=meta.get('candidate_bank',{})
  if saved:assert saved['path']==bp and saved['sha256']==sha(bp)
  bank=json.loads(gzip.decompress((ROOT/bp).read_bytes()));paths.extend([meta_path,bp]);V=np.array([r['normalized_wall_pairings']for r in bank]);old=V.shape[1]
  assert 12<=old<=33
  vals=V[:,BASIS]@transform.T
  # Retain inherited pairings for the old walls: the banks' original gates
  # remain numerical premises. The transform is used for later walls only.
  later=vals[:,old:52];acute=np.max(later,axis=1)<1e-6
  finite=(later>-2+1e-6)&(later<-1e-6);m=np.pi/np.arccos(np.clip(-later/2,-1+1e-15,1-1e-15));err=np.abs(m-np.rint(m))
  bad=np.sum(finite&(err>1e-5),axis=1);near=np.sum(abs(later+2)<=1e-6,axis=1)
  center=vals[:,-1]<-1e-8;cut=vals[:,-2]>1e-7
  good=acute&center&(bad==0)&(near==0);cutters=np.flatnonzero(good&cut).tolist()
  close=np.flatnonzero(acute&center&cut&(bad<=2)).tolist()
  details=[]
  for i in sorted(cutters,key=lambda i:-vals[i,-2]):
   recipe={k:v for k,v in bank[i].items()if k not in['normalized_wall_pairings','endpoint_other_pairings','physical_normal']}
   details.append(dict(index=i,recipe=recipe,escape_pairing=smooth(vals[i,-2]),new_pairings=[smooth(x)for x in vals[i,old:52]],normalized_basis_pairings=V[i,BASIS].tolist()))
  record=dict(bank=bp,count=len(bank),old_wall_count=old,compatible_proposals=int(good.sum()),escape_cutters=len(cutters),cutters=details,acute_cutters_with_at_most_two_unresolved=len(close),close_indices=close,maximum_compatible_escape_pairing=smooth(max(vals[good,-2]))if any(good)else None,accepted_mask_sha256=sha256(np.packbits(good).tobytes()).hexdigest())
  records.append(record);print(tag,{k:v for k,v in record.items()if k not in['cutters','close_indices']},flush=True)
 return dict(status='PASS_BOUNDED_NUMERICAL_52_WALL_NONPOLAR_BANK_RESCORE',sources={p:sha(p)for p in paths},pairing_basis=BASIS,escape_normalization='E/(-E.H)',banks=records,total_candidates=sum(r['count']for r in records),escape_cutter_proposals=sum(r['escape_cutters']for r in records),scope='Nine frozen bounded numerical banks only. Their old contact gates are inherited numerically; new contacts use binary64 rescore with1e-6acuteness/parallel margin and1e-5integer-label tolerance. No exact cap, complete nonpolar family or all-label exclusion follows. Every survivor must be exactified against all52walls.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved bank rescore replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['total_candidates'],d['escape_cutter_proposals'])
