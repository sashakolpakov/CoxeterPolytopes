#!/usr/bin/env python3
"""Bounded coordinated H3 repair trials in four-dimensional projected coordinates."""
import argparse,itertools,json,math
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_h3_closure_20261016_scout as H3
P=H3.P;ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_closure_20261016_repair_trials';S=H3.S
TRIAL='results/fresh_nonarithmetic_h3_closure_20261016_trial.json';SCOUT='results/fresh_nonarithmetic_h3_closure_20261016_repair_scout.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def geometry():
 Q,N,H,E=P.geometry();cap=P.P.B.read(H3.CAP);N.append(mp.matrix([P.mid(x)for x in cap['normal_interval']]));tr=P.P.B.read(TRIAL);N.extend(mp.matrix([mp.mpf(x)for x in row])for row in tr['normal_intervals_or_proposals'])
 W=mp.matrix([[N[i][j]for i in S]for j in range(10)]);J=(W.T*Q*W)**-1
 proj=lambda n:n-W*J*W.T*Q*n
 F=mp.matrix([[proj(N[i])[j]for i in[7,8,30,33]]for j in range(10)]);G=F.T*Q*F;Gi=G**-1;lift=lambda n:Gi*F.T*Q*n
 Ns=np.array([list(lift(n))for n in N],dtype=float);g=np.array(G.tolist(),dtype=float);h=np.array(list(lift(H)),dtype=float);hc=g@h
 rays=[dict(mask=sum(1<<i for i in row['active']),v=np.array(list(lift(mp.matrix([mp.mpf(x)for x in row['point']]))),dtype=float))for row in tr['rays']]
 return tr,Ns,g,h,hc,rays
def null(rows):
 a=np.array(rows);_,s,vt=np.linalg.svd(a);assert s[-1]>1e-7;return vt[-1]
def clip(rays,n,wall,g,hc):
 cov=g@n;gs=[v['v']@cov for v in rays];neg=[i for i,x in enumerate(gs)if x<-1e-8];pos=[i for i,x in enumerate(gs)if x>1e-8];zero=[i for i,x in enumerate(gs)if abs(x)<=1e-8];new={rays[i]['mask']:rays[i]for i in neg}
 for i in zero:new[rays[i]['mask']|(1<<wall)]=dict(mask=rays[i]['mask']|(1<<wall),v=rays[i]['v'])
 for i,j in itertools.product(neg,pos):
  m=rays[i]['mask']&rays[j]['mask']
  if m.bit_count()<8:continue
  if sum(z['mask']&m==m for z in rays)!=2:continue
  v=gs[j]*rays[i]['v']-gs[i]*rays[j]['v'];v/=-(v@hc);key=m|(1<<wall)
  if key in new:assert np.max(abs(new[key]['v']-v))<1e-6
  new[key]=dict(mask=key,v=v)
 return list(new.values())
def badpairs(vals):
 if max(vals)>1e-6:return[['obtuse',int(np.argmax(vals)),float(max(vals))]]
 out=[]
 for j,x in enumerate(vals):
  if abs(x)<1e-6 or x<-2-1e-6:continue
  if abs(x+2)<1e-6:out.append(['near_parallel',j]);continue
  m=math.pi/math.acos(np.clip(-x/2,-1+1e-15,1-1e-15))
  if abs(m-round(m))>1e-5:out.append(['gap',j,float(m)])
 return out
def run(limit=8):
 tr,N,g,h,hc,rays=geometry();sc=P.P.B.read(SCOUT);base=len(N)
 counts=Counter();records=[];best=[];candidates=[(i,r)for i,r in enumerate(sc['candidates'])if r['target']in[1,3,4,5]and max(r['labels'])<=limit]
 for ordinal,(ix,r)in enumerate(candidates):
  frame=r['anchors'];W=N[frame].T;J=np.linalg.inv(W.T@g@W);c=np.array([-2*math.cos(math.pi/m)if m!=2 else 0 for m in r['labels']]);p=W@J@c;k=h-W@J@W.T@g@h;qk=k@g@k;rho=(2-p@g@p)/qk
  assert rho>0;n=p+r['branch']*math.sqrt(rho)*k;bad=badpairs(N@g@n)
  if bad:counts['incompatible_with_trial']+=1;continue
  rr=clip(rays,n,len(N),g,hc);NN=np.vstack([N,n]);newpol=[];bads=[];finite=ideal=0
  for v in rr:
   q=v['v']@g@v['v']
   qt=1e-11*np.sum(abs(np.outer(v['v'],v['v'])*g))
   if q<-qt:finite+=1;continue
   if abs(q)<=qt:ideal+=1;continue
   local=bool(v['mask']&(1<<base))or bool(v['mask']&(1<<(57 if r['target']in[1,3]else 58)))
   if not local:continue
   z=v['v']*math.sqrt(2/q);b=badpairs(NN@g@z);newpol.append(z)
   if b:bads.append(dict(active=[i for i in range(len(NN))if v['mask']>>i&1],bad=b))
  mutual=[]
  if not bads and not ideal:
   for i,j in itertools.combinations(range(len(newpol)),2):
    b=badpairs([newpol[i]@g@newpol[j]])
    if b:mutual.append(dict(pair=[i,j],bad=b))
  counts['tested']+=1;rec=dict(candidate=ix,region='B'if r['target']in[1,3]else'C',anchors=r['anchors'],labels=r['labels'],branch=r['branch'],finite=finite,ideal=ideal,local_exterior=len(newpol),bad_polars=bads,mutual_bad=mutual)
  if len(bads)<=1 or not mutual and not bads:records.append(rec)
  if not bads and not mutual and not ideal:best.append(rec);print('NUMERICAL_LOCAL_TRUNCATION',rec,flush=True)
  if ordinal%100==0:print('progress',ordinal,'/',len(candidates),dict(counts),'best',len(best),flush=True)
 paths=[f'scripts/{TAG}.py',TRIAL,SCOUT,H3.CAP,H3.CONE,f'scripts/{H3.TAG}.py',P.CL,P.CO,P.P.BATCH]
 return dict(status='BOUNDED_NUMERICAL_COORDINATED_H3_REPAIR_TRIALS',sources={p:sha(p)for p in paths},max_anchor_label=limit,stats=dict(counts),retained_records=records,complete_truncation_proposals=best,scope='Double-precision discovery only, using high-precision projected geometry. Tests the saved six-cap trial plus one label-bounded repair, then all replacement polars in its B or C region and their mutual gates. Other regions remain open; combining regional proposals requires further checks. No exact closure, affine-filter, actual-facet or ambient finite-volume claim is made.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--max-label',type=int,default=8);a=ap.parse_args();d=json.loads(json.dumps(run(a.max_label)));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved coordinated repair trials')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['complete_truncation_proposals']))
