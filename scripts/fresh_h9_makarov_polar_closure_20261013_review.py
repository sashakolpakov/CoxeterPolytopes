#!/usr/bin/env python3
"""Separate certificate-consistency and Boolean-cover replay.

This does not reconstruct algebraic polars or their pairing intervals. The
primary saved replay and accompanying source review cover that arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_makarov_polar_closure_20261013'
SOURCE=ROOT/f'results/{TAG}_verify.json'
OUT=ROOT/f'results/{TAG}_review.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bits(m):
 while m:
  b=m&-m;yield b.bit_length()-1;m-=b
def main():
 d=json.loads(SOURCE.read_text());n=d['polar_count'];nw=d['witness_count'];S=int(d['interval_scale'])
 assert n==1564 and nw==21
 for p,h in d['sources'].items():assert sha(ROOT/p)==h
 cuts=[];masks=[0]*n
 for k,w in enumerate(d['witnesses']):
  assert int(w['square_interval'][0])>0 and int(w['center_pair_interval'][1])<0
  assert len(w['strict_old_pair_intervals'])==27 and all(int(x[1])<0 for x in w['strict_old_pair_intervals'])
  ps=w['all_polar_pair_intervals'];assert len(ps)==n and all(int(a)<=int(b)for a,b in ps)
  ids=[i for i,p in enumerate(ps)if int(p[1])>=0];assert ids==w['possible_cutters']
  cuts.append(sum(1<<i for i in ids))
  for i in ids:masks[i]|=1<<k
 union=0
 for m in cuts:union|=m
 assert union.bit_count()==d['candidates_cutting_some_witness']==1515
 A=d['angle_intervals'];scale=int(A['pair_scale']);tail=F(int(A['tail_upper']),scale)
 angles=[(F(int(a['lower']),scale),F(int(a['upper']),scale))for a in A['finite_angles']]
 assert [a['label']for a in A['finite_angles']]==list(range(2,121))
 neighbors=[union&~(1<<i)for i in range(n)];seen=set();reasons={}
 for c in d['exact_conflicts']:
  i,j=c['polar_indices'];assert 0<=i<j<n and (i,j)not in seen;seen.add((i,j))
  lo,hi=(F(int(v),S)for v in c['pair_interval']);assert lo<=hi
  if c['reason']=='positive':assert lo>0
  else:
   assert c['reason']=='finite_angle_gap' and lo>tail
   assert all(hi<a or lo>b for a,b in angles)
  neighbors[i]&=~(1<<j);neighbors[j]&=~(1<<i)
  reasons[c['reason']]=reasons.get(c['reason'],0)+1
 assert len(seen)==31308
 states=0
 @lru_cache(None)
 def cover(available,missing):
  nonlocal states
  states+=1
  if not missing:return True
  options=min((available&cuts[k]for k in bits(missing)),key=int.bit_count)
  if not options:return False
  # Deliberately use literal index order rather than the producer's
  # descending witness-coverage order.
  for i in bits(options):
   if cover(available&neighbors[i],missing&~masks[i]):return True
   available&=~(1<<i)
  return False
 assert not cover(union,(1<<nw)-1)
 return {'status':'PASS_SEPARATE_POLAR_CERTIFICATE_CONSISTENCY_AND_BOOLEAN_COVER',
  'sources':{str(p.relative_to(ROOT)):sha(p)for p in[SOURCE,Path(__file__)]},
  'polar_count':n,'witness_count':nw,'pair_interval_count':n*nw,
  'rechecked_conflicts':len(seen),'conflict_reasons':reasons,'separate_DFS_states':states,
  'scope':'Rechecks saved sign/cutter consistency, every conflict interval against the saved all-label enclosures, and UNSAT of the symmetric compatibility supergraph by separate DFS. Does not independently reconstruct normals, pairing intervals, or trigonometric enclosures; those were source-reviewed and primary-replayed.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();r=main()
 if a.verify:assert r==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(r,indent=2)+'\n')
 print(r['status'],r['separate_DFS_states'])
