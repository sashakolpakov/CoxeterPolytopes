#!/usr/bin/env python3
"""All-label exact gap rejection of the remaining1053 exterior polars."""
import argparse,json,math
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_joint_regions_20261016_batch as B
D,A,R,C,K,Z,O=B.D,B.A,B.R,B.C,B.K,B.Z,B.O
I=R.L.I;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_joint_regions_20261016_polar_gate'
BATCH=f'results/{B.TAG}.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def atan_recip(n,terms=80):
 x=F(1,n);s=sum(((-1)**k*x**(2*k+1)/F(2*k+1)for k in range(terms)),F());nextterm=x**(2*terms+1)/F(2*terms+1)
 return I(s,s+nextterm)if terms%2==0 else I(s-nextterm,s)
PI=B.rnd(16*atan_recip(5)-4*atan_recip(239),256)
@lru_cache(maxsize=None)
def threshold(m):
 assert m>=2
 if m==2:return I()
 x=PI/m;t=I(1);s=I(1)
 for k in range(1,51):
  t=B.rnd(t*x*x/F((2*k-1)*(2*k)),256)
  s=B.rnd(s+((-1)**k)*t,256)
 rem=F(x.h)**102/F(math.factorial(102));s=I(s.l-rem,s.h+rem)
 assert s.l>0
 return B.rnd(4*s*s,192)
def run():
 batch=B.read(BATCH);scout=B.read(B.SCOUT);local,cone,Q,H,N,cov,ray=B.geometry()
 screens={r['ray_index']:r for r in scout['all_ray_screens']};survivors=set(B.RAY_IDS)
 assert len(screens)==1068 and len(survivors)==15
 assert batch['status']=='PASS_EXACT_FIFTEEN_POLARS_AND_COMPATIBLE_TWELVE_CAP_BATCH'
 norms=[D.dot(n,c)for n,c in zip(N,cov)];records=[];hist=Counter();additional=[]
 E=tuple(tuple(A.B.db(x)for x in a)for a in batch['common_strict_escape']['point'])
 chosen=[tuple(tuple(A.B.db(x)for x in a)for a in r['normal'])for r in batch['candidates']]
 chosen_norms=[D.pair(v,v,Q)for v in chosen]
 for i,ids in enumerate(cone['complete_census']['hyperideal_rays']):
  if i in survivors:continue
  V=ray(35,tuple(ids));q=D.pair(V,V,Q);qiv=B.iv(q);assert qiv.l>0
  hints=[x for x in screens[i]['unresolved_finite']if isinstance(x[-1],float)]
  hints.sort(key=lambda x:-min(x[-1]%1,1-x[-1]%1))
  accepted=None
  for wall,_,m in hints:
   g=D.dot(V,cov[wall]);giv=B.iv(g);assert giv.h<0
   cyc=B.rnd(4*giv*giv/(qiv*B.iv(norms[wall])),192)
   for lower in sorted(set([max(2,math.floor(m)+a)for a in[0,-1,1]])):
    lo,hi=threshold(lower),threshold(lower+1)
    if lo.h<cyc.l and cyc.h<hi.l:
     accepted=dict(ray_index=i,active=ids,wall=wall,neighbor_labels=[lower,lower+1],cycle_interval=cyc.enc(),lower_threshold_interval=lo.enc(),upper_threshold_interval=hi.enc());hist[lower]+=1;break
   if accepted:break
  if accepted is None:
   vals=[D.dot(V,c)for c in cov];rels=[B.relation(g,D.mul(q,nn))for g,nn in zip(vals,norms)]
   assert all(r['kind']in['right','finite','dotted']for r in rels),('No exact gap or valid polar',i,hints,rels)
   assert D.sign(D.pair(V,H,Q))<0 and D.sign(D.pair(V,E,Q))<0
   mutual=[B.relation(D.pair(V,w,Q),D.mul(q,nn))for w,nn in zip(chosen,chosen_norms)]
   assert all(r['kind']=='dotted'for r in mutual)
   additional.append(B.encode(dict(ray_index=i,active=ids,normal=V,norm=q,pairings=vals,relations=rels,center_pairing=D.pair(V,H,Q),common_escape_pairing=D.pair(V,E,Q),mutual_with_selected_fifteen=mutual)))
   print('additional exact polar',i,flush=True);continue
  records.append(accepted)
  if len(records)%100==0:print('exact polar gaps',len(records),flush=True)
 assert len(records)+len(additional)==1053
 paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',BATCH,B.CONE,B.PRIMARY,B.SCOUT]
 return dict(status='PASS_COMPLETE_ALL_LABEL_1068_EXTERIOR_POLAR_SCREEN',sources={p:sha(p)for p in paths},complete_polar_count=1068,compatible_polar_count=15+len(additional),compatible_ray_indices=sorted(survivors|{r['ray_index']for r in additional}),additional_exact_polars=additional,rejected_polar_count=len(records),gap_records=records,gap_histogram={str(k):v for k,v in sorted(hist.items())},pi_interval=PI.enc(),threshold_intervals={str(m):threshold(m).enc()for m in sorted({k+a for k in hist for a in[0,1]})},
  proof='4cos²(pi/m) is strictly increasing for integer m>=2 and tends to4. Every rejected polar has one squared normalized product strictly between consecutive thresholds, certified by rational intervals. Pi is enclosed by Machin identity and alternating arctan series; cosine uses50 Taylor terms with an explicit102nd-order remainder. No finite-angle alphabet or cutoff is assumed.',
  scope='Complete all-label classification of polar hyperplanes of the fixed36-wall cone only. All16 compatible polars retain the same exact strict escape, so no subset of this one-generation bank completes finite volume. The additionally recovered wall is dotted to the initial15 and right to recentwall33. New polars of enlarged cones and nonpolar walls remain open.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved all-label polar screen replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
