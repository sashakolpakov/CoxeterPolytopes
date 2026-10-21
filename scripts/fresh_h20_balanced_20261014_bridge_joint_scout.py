#!/usr/bin/env python3
"""Numerical complementary-bank probe and bounded non-endpoint boundary caps."""
import json,math,time
from itertools import product,combinations
import numpy as np
from fresh_h20_balanced_20261014 import ROOT,data,sha
TAG='fresh_h20_balanced_20261014_bridge_joint_scout'
EP='data/enumerations/fresh_h20_balanced_20261014_unit_endpoint_bank.npz'
def run():
 st=time.time();D=data();V=np.array([v['point']for v in D['vertices']]);cs=np.array([math.cos(math.pi/m)for m in range(2,61)])
 rec=[]
 for n,m in product([3,4,5,6],[2,3,4,5,6,0]):
  a=math.cos(math.pi/n);b=1 if m==0 else math.cos(math.pi/m)
  A=196-1/a**2;B=2/a**2+4*b/a;C=-4-1/a**2-4*b/a
  g=(-B+math.sqrt(B*B-4*A*C))/(2*A);r=(1-g)/a
  gaps=(1-g+g*np.array([7,7,8]))/r;assert min(gaps)>1
  rec.append(dict(seed_label=n,central_label=m,gamma=g,radius=r))
 C=np.concatenate([r['gamma']*V for r in rec]);R=np.concatenate([np.full(392,r['radius'])for r in rec]);ds=np.unique(np.round(392-2*V@V.T,8));finitehits=[]
 for i,a in enumerate(rec):
  for j,b in enumerate(rec[:i+1]):
   d2=196*(a['gamma']-b['gamma'])**2+a['gamma']*b['gamma']*ds
   co=(d2-a['radius']**2-b['radius']**2)/(2*a['radius']*b['radius'])
   for k,l in np.argwhere(abs(co[:,None]-cs)<1e-8):finitehits.append([i,j,float(ds[k]),int(l+2)])
 fs=[.15,3/14,.25,.3,.35,.4,.5,.6,.75,1.];X=np.concatenate([f*V for f in fs]);q=np.sum(X*X,axis=1);power=q-4
 e=np.load(ROOT/EP)
 for CC,RR in[(C,R),(e['centers'],e['radii'])]:
  for i in range(0,len(CC),500):
   cc=CC[i:i+500];rr=RR[i:i+500];pp=q[:,None]+np.sum(cc*cc,axis=1)-2*X@cc.T-rr**2
   power=np.minimum(power,np.min(pp,axis=1))
 radial=[dict(fraction=f,uncovered=int(np.sum(power[k*392:(k+1)*392]>1e-8)),worst_power=float(np.max(power[k*392:(k+1)*392])))for k,f in enumerate(fs)]
 # Boundary-center caps: nineteen prescribed right/pi3 contacts, one omitted
 # dotted anchor in each factor. No PSD endpoint residual restriction is used.
 I=D['vertices'][0]['retained'];W=D['vectors'][I];dual=np.linalg.inv(W@W.T)@W;z=V[0];target=.4*V
 blocks=[range(6),range(6,12),range(12,19)];sets=[[s for k in range(4)for s in combinations(b,k)]for b in blocks]
 accepted=[];best=None;cutters=[];trials=0
 for parts in product(*sets):
  inds=sum((list(p)for p in parts),[]);y=.5*dual[inds].sum(axis=0)if inds else np.zeros(19);y2=y@y;zy=z@y;count=np.array(list(map(len,parts)))
  for m in list(range(2,17))+[0]:
   trials+=1;co=1 if m==0 else math.cos(math.pi/m);A=y2-1;B=-2*zy-4*co;cc=192;dis=B*B-4*A*cc
   if abs(A)<1e-10:rr=[-cc/B]
   elif dis>=0:rr=[(-B+math.sqrt(dis))/(2*A),(-B-math.sqrt(dis))/(2*A)]
   else:continue
   for r in rr:
    if r<=0:continue
    gaps=(np.array([7,7,8])-r*count/2)/r
    if min(gaps)<1-1e-8:continue
    c=z-r*y;pp=np.sum((target-c)**2,axis=1)-r*r;j=int(np.argmin(pp));row=dict(parts=[list(p)for p in parts],central_label=m,radius=r,schur_cost=float(y2),target_index=j,minimum_power=float(pp[j]))
    accepted.append(row)
    if best is None or row['minimum_power']<best['minimum_power']:best=row
    if pp[j]<-1e-8:cutters.append(row)
 out=dict(status='BOUNDED_NUMERICAL_JOINT_BRIDGE_ENDPOINT_PROBE',sources={p:sha(p)for p in['scripts/fresh_h20_balanced_20261014.py',f'scripts/{TAG}.py',EP,'results/fresh_h20_balanced_20261014_unit_endpoint_scout.json']},radial_profiles=rec,radial_bank_size=len(C),radial_new_new_finite_contact_hits=finitehits,finite_pair_labels_tested=list(range(2,61)),complementary_endpoint_bank_size=len(e['centers']),combined_radial_probe_summary=radial,boundary_search=dict(subsets_per_factor_maximum_size=3,positive_contact_label=3,central_labels=list(range(2,17))+['parallel'],profile_label_trials=trials,accepted_count=len(accepted),all_orbit_target_cutters=len(cutters),best=best,cutters=cutters),elapsed_seconds=time.time()-st,scope='Numerical bounded scouts only. Pair labels2..60 do not settle arbitrary finite labels. Combined bank coverage is a necessary probe test, not a volume certificate. Boundary recipe gates permit zero/negative endpoint residual and test all392 translates ofthe .4vertex target, but do not exhaust arbitrary seed-contact profiles.')
 (ROOT/'results'/f'{TAG}.json').write_text(json.dumps(out,indent=2)+'\n');print('radial',radial,'boundary',out['boundary_search'])
if __name__=='__main__':run()
