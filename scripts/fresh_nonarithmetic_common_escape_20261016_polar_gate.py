#!/usr/bin/env python3
"""Exact all-label test of every48-cone polar that could cut common E."""
import argparse,json,math
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_joint_regions_20261016_batch as B
import fresh_nonarithmetic_joint_regions_20261016_polar_gate as T
D,A,R,C,K,Z,O=B.D,B.A,B.R,B.C,B.K,B.Z,B.O
I=R.L.I;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_common_escape_20261016_polar_gate'
BATCH=f'results/{B.TAG}.json';CONE='results/fresh_nonarithmetic_joint_regions_20261016_cone.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def decode(a):return tuple(A.B.db(x)for x in a)
def ii(a):
 z=I()
 for x,b in zip(a,D.BIV):
  if x:z=z+D.biv(x)*b
 return B.rnd(z,256)
def dot(v,w):return B.rnd(sum((a*b for a,b in zip(v,w)),I()),256)
def mv(Q,v):return[dot(row,v)for row in Q]
def run():
 batch=B.read(BATCH);cone=B.read(CONE);_,old,Q,H,_,_,base_ray=B.geometry()
 N=[tuple(decode(x)for x in n)for n in batch['batch_normals']];cov=[tuple(D.mv(Q,n))for n in N]
 covI=[[ii(x)for x in c]for c in cov];nnI=[ii(decode(x))for x in batch['batch_norms']]
 QI=[[ii(D.lift(x))for x in row]for row in Q]
 E=tuple(decode(x)for x in batch['common_strict_escape']['point']);eI=[ii(x)for x in D.mv(Q,E)]
 crossings={s['wall']:{tuple(x['new_support']):x for x in s['crossings']}for s in cone['clipping_stages']}
 scales={};fallback=Counter()
 @lru_cache(maxsize=None)
 def exact(stage,ids):
  if stage==35:return base_ray(35,ids)
  if stage not in ids:return exact(stage-1,ids)
  rec=crossings[stage].get(ids)
  if rec is None:
   out=exact(stage-1,tuple(i for i in ids if i!=stage));assert D.sign(D.dot(out,cov[stage]))==0;return out
  a=exact(stage-1,tuple(rec['negative_support']));b=exact(stage-1,tuple(rec['positive_support']))
  ga,gb=D.dot(a,cov[stage]),D.dot(b,cov[stage]);assert D.sign(ga)<0<D.sign(gb)
  out=B.vecsub(B.vecscale(a,gb),B.vecscale(b,ga))
  return tuple(D.scale(x,K.convert(scales[(stage,ids)]))for x in out)
 @lru_cache(maxsize=None)
 def ray(stage,ids):
  if stage==35:return[ii(x)for x in base_ray(35,ids)]
  if stage not in ids:return ray(stage-1,ids)
  rec=crossings[stage].get(ids)
  if rec is None:
   out=ray(stage-1,tuple(i for i in ids if i!=stage))
   assert D.sign(D.dot(exact(stage-1,tuple(i for i in ids if i!=stage)),cov[stage]))==0
   return out
  a=ray(stage-1,tuple(rec['negative_support']));b=ray(stage-1,tuple(rec['positive_support']));ga,gb=dot(a,covI[stage]),dot(b,covI[stage])
  if not ga.h<0:ga=ii(D.dot(exact(stage-1,tuple(rec['negative_support'])),cov[stage]));fallback['crossing_sign']+=1
  if not gb.l>0:gb=ii(D.dot(exact(stage-1,tuple(rec['positive_support'])),cov[stage]));fallback['crossing_sign']+=1
  assert ga.h<0<gb.l
  out=[B.rnd(gb*x-ga*y,256)for x,y in zip(a,b)];bound=max(max(abs(x.l),abs(x.h))for x in out);p=bound.numerator.bit_length()-bound.denominator.bit_length();s=F(1,1<<p)if p>=0 else F(1<<-p)
  scales[(stage,ids)]=s
  return[B.rnd(s*x,256)for x in out]
 records=[];noncut=[];unresolved=[];hist=Counter()
 for index,active in enumerate(cone['complete_census']['hyperideal_rays']):
  ids=tuple(active);v=ray(47,ids);ep=dot(v,eI)
  if ep.h<0:noncut.append(dict(ray_index=index,escape_pairing=ep.enc()));continue
  if not ep.l>0:
   sg=D.sign(D.pair(exact(47,ids),E,Q));fallback['escape_sign']+=1
   if sg<=0:noncut.append(dict(ray_index=index,escape_pairing=ep.enc(),exact_sign=sg));continue
  q=dot(v,mv(QI,v))
  if not q.l>0:q=ii(D.pair(exact(47,ids),exact(47,ids),Q));fallback['norm']+=1
  assert q.l>0
  proof=None
  for wall in range(48):
   if wall in ids:continue
   g=dot(v,covI[wall])
   if not g.h<0:g=ii(D.dot(exact(47,ids),cov[wall]));fallback['old_pairing']+=1
   assert g.h<0
   cyc=B.rnd(4*g*g/(q*nnI[wall]),192)
   if cyc.l>=4:continue
   if cyc.h>=4 or cyc.l<=0:continue
   mid=float((cyc.l+cyc.h)/2);m=math.pi/math.acos(math.sqrt(mid/4))
   for lower in sorted({max(2,math.floor(m)+s)for s in[-1,0,1]}):
    lo,hi=T.threshold(lower),T.threshold(lower+1)
    if lo.h<cyc.l and cyc.h<hi.l:
     proof=dict(ray_index=index,active=active,escape_pairing=ep.enc(),wall=wall,neighbor_labels=[lower,lower+1],cycle_interval=cyc.enc(),lower_threshold_interval=lo.enc(),upper_threshold_interval=hi.enc());hist[lower]+=1;break
   if proof:break
  if proof:records.append(proof)
  else:unresolved.append(dict(ray_index=index,active=active,escape_pairing=ep.enc()));print('UNRESOLVED CUTTER',index,active,flush=True)
  if index%200==0:print('processed',index,'cut-gap',len(records),'noncut',len(noncut),flush=True)
 assert len(records)+len(noncut)+len(unresolved)==2109
 paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{T.TAG}.py',BATCH,CONE,B.CONE]
 return dict(status='PASS_NO_ESCAPE_CUTTING_POLAR_IN_COMPLETE_48_CONE'if not unresolved else'EXACT_POLAR_ESCAPE_FRONTIER_UNRESOLVED',sources={p:sha(p)for p in paths},complete_hyperideal_count=2109,strictly_noncutting_count=len(noncut),forbidden_cutting_count=len(records),unresolved_cutters=unresolved,noncutting=noncut,cutting_gap_records=records,gap_histogram={str(k):v for k,v in sorted(hist.items())},exact_fallbacks=dict(fallback),
  proof='All2109 complete-cone hyperideal rays are reconstructed through exact crossing expressions. Every noncutting sign has an outward rational certificate or exact fallback. Every positive escape pairing is accompanied by a squared normalized product strictly between consecutive Coxeter-angle thresholds. No finite-label cutoff or numerical rejection is used.',
  scope='No polar wall of this fixed48-wall cone can cut the current common strict escape. This does not exclude nonpolar walls, changed retained walls, or polars of later enlarged chambers. The separate46-proposal numerical screen is not used as a completeness premise.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact escape-polar gate replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['strictly_noncutting_count'],d['forbidden_cutting_count'],len(d['unresolved_cutters']))
