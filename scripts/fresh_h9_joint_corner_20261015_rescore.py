#!/usr/bin/env python3
"""Rescore the old two-anchor bank against both exactjoint caps and their newescape."""
import argparse,gzip,json,math
from hashlib import sha256
from pathlib import Path
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_corner_20261015_rescore'
BASE='results/fresh_h9_joint_pair_20261015_two_extension.json';BANK='results/fresh_h9_joint_pair_20261015_scout_escape.json'
def scalar(x):return mp.mpf(x[0][0])+mp.mpf(x[0][1])*mp.sqrt(5)+(mp.mpf(x[1][0])+mp.mpf(x[1][1])*mp.sqrt(5))*mp.sqrt(6+16*mp.sqrt(5)/5)
def run():
 mp.mp.dps=100;d=json.loads((ROOT/BASE).read_text());b=json.loads((ROOT/BANK).read_text());bp=b['candidate_bank']['path'];blob=(ROOT/bp).read_bytes();assert sha256(blob).hexdigest()==b['candidate_bank']['sha256'];rows=json.loads(gzip.decompress(blob))
 Q=np.array([[float(scalar(x))for x in row]for row in d['form']]);T=np.zeros((10,10));T[0,0]=math.sqrt(-Q[0,0]);T[1:,1:]=np.linalg.cholesky(Q[1:,1:]);metric=np.diag([-1.]+[1.]*9)
 N=np.array([[float(scalar(x))for x in row]for row in d['normals']])@T;N*=np.array([math.sqrt(2/float(scalar(x)))for x in d['norms']])[:,None]
 E=np.array([float(scalar(x))for x in d['strict_escape']['point']])@T;normals=np.array([r['physical_normal']for r in rows]);pairs=normals@metric@N.T
 values=-2*np.cos(np.pi/np.arange(2,257));err=np.min(abs(pairs[:,:,None]-values),axis=2);ok=np.all((pairs<-2-1e-7)|(err<1e-7),axis=1);esc=normals@metric@E
 compatible=np.flatnonzero(ok);cutters=sorted([int(i)for i in compatible if esc[i]>1e-7],key=lambda i:-esc[i])
 records=[dict(bank_index=i,new_escape_pairing=float(esc[i]),pairings_to_joint_caps=pairs[i,30:].tolist(),recipe=rows[i])for i in cutters]
 paths=[BASE,BANK,bp,f'scripts/{TAG}.py']
 return dict(status='BOUNDED_32WALL_RESCORE_OF_RANK7_BANK',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},bank_count=len(rows),compatible_count=len(compatible),escape_cutter_count=len(cutters),cutters=records,
   best_compatible_escape_pairing=float(max(esc[compatible]))if len(compatible)else None,
   scope='Numerical rescore of4566 frozen proposals; same finiteangle tolerance1e-7 and labels2..256. No exact newwall or allparameter exclusion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['compatible_count'],d['escape_cutter_count'],d['best_compatible_escape_pairing'])
