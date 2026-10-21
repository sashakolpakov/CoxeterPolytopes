#!/usr/bin/env python3
"""Numerical proposal of continuous two-wall dual certificates, no cap scan."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json,numpy as np
from scipy.linalg import null_space
import fresh_h9_joint_frontier_20261013 as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_frontier_20261013_pd8_relax'
BASE=ROOT/'results/fresh_h9_joint_frontier_20261013.json';EXTRA=ROOT/'results/fresh_h9_joint_frontier_20261013_all_edges.json'
GAP=ROOT/'results/fresh_h9_makarov_joint_completion_20261013_old_bank_gate.json';OUT=ROOT/'results'/f'{TAG}.json'
def run():
 src,Qe,We,_,ce,G,adj=A.setup();Q=np.array([[float(a)+np.sqrt(5)*float(b)for a,b in row]for row in Qe]);W=np.array([[A.approx(x)for x in row]for row in We]);c=np.array([A.approx(x)for x in ce]);E=np.array(json.loads(GAP.read_text())['escape']['point'],float);E/=np.linalg.norm(E)
 base=json.loads(BASE.read_text());extra=json.loads(EXTRA.read_text());rows=[];counts=Counter()
 queue=[(r['id'],r['support'],base['open_edge_spherical_profiles'][str(r['id'])])for r in base['edges']if r['kind']=='open_ideal_boundary']
 queue +=[(r['id'],r['support'],r['spherical_profiles'])for r in extra['new_edges_without_finite_endpoint']]
 for eid,ids,profiles in queue:
  WS=W[ids];M=np.array([[float(a)+np.sqrt(5)*float(b)for a,b in row]for row in [[G[i][j]for j in ids]for i in ids]]);B=np.linalg.inv(M)
  p=c-(B@(WS@Q@c))@WS;assert p@Q@p<0;u=p/np.sqrt(-p@Q@p)
  e=null_space(np.vstack([WS@Q,u@Q]))[:,0];assert e@Q@e>0;e/=np.sqrt(e@Q@e)
  V=W@Q@np.column_stack([u,e]);target=E@Q@np.column_stack([u,e]);outside=[i for i in range(27)if i not in ids];duals=[]
  for i,j in combinations(outside,2):
   C=V[[i,j]].T
   if abs(np.linalg.det(C))<1e-11:continue
   lam=np.linalg.solve(C,target)
   if min(lam)<-1e-10:continue
   remainder=E-lam@W[[i,j]];mu=B@(WS@Q@remainder)
   if np.max(abs(remainder-mu@WS))>1e-7:continue
   duals.append((i,j,lam,mu))
  tests=[]
  for kind,rr in [('finite',profiles['finite_profiles']),('tail',profiles['dihedral_tail_profiles'])]:
   for index,row in enumerate(rr):
    h=np.zeros(8)
    if kind=='tail':h[ids.index(row['wall'])]=1
    else:
     for wall,m in row['contacts']:h[ids.index(wall)]=2*np.cos(np.pi/m)
    if not duals:
     tests.append({'kind':kind,'index':index,'status':'NO_TWO_WALL_DUAL'});counts['no_dual']+=1;continue
    scores=[float(x[3]@h)for x in duals];owner=int(np.argmax(scores));i,j,lam,mu=duals[owner];value=scores[owner]
    tag='STRICT_DUAL_PROPOSAL'if value>1e-9 else'ZERO_DUAL_PROPOSAL'if abs(value)<=1e-9 else'CONTINUOUS_NORM_GATE_NEEDED'
    counts[tag]+=1;tests.append({'kind':kind,'index':index,'status':tag,'walls':[i,j],'lambda':lam.tolist(),'support_coefficients':mu.tolist(),'margin':value})
  rows.append({'edge':eid,'support':ids,'nonnegative_dual_pairs':len(duals),'profile_tests':tests})
 return {'status':'NUMERICAL_CONTINUOUS_DUAL_PROPOSALS_ONLY','sources':{str(p.relative_to(ROOT)):A.sha(p)for p in[BASE,EXTRA,GAP,Path(A.__file__),Path(__file__)]},'edge_supports':len(rows),'counts':dict(counts),'rows':rows,'scope':'Proposal search for E=sum(lambda_i*w_i)+sum(mu_j*support_j), lambda>=0. If mu dot h>=0, a fixed-profile acute normal cannot cut E even without its norm equation. None of these floating-point signs is an exact exclusion. Unresolved rows require a continuous norm-hyperbola gate, not a larger finite-anchor sample.'}
if __name__=='__main__':
 d=run();OUT.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'],flush=True)
