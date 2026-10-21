#!/usr/bin/env python3
"""Bounded release of one source-right contact, fixing the offending anchor."""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_joint_route_20261015_release'
SRC='results/fresh_h20_a1a1_h18_coupling_20261015.json';BANK='results/fresh_h20_joint_route_20261015_scout.json'
def run():
 d=json.loads((ROOT/SRC).read_text());bank=json.loads((ROOT/BANK).read_text())
 dec=lambda x:float(Fraction(x['a']))+np.sqrt(2)*float(Fraction(x['b']))
 Q=np.array([[dec(x)for x in row]for row in d['source_form']]);N=np.array([[dec(x)for x in row]for row in d['source_normals']]);H=np.array([dec(x)for x in d['source_center']]);E=np.array([dec(x)for x in d['strict_escape']]);q=np.einsum('ij,jk,ik->i',N,Q,N)
 cov=N@Q;seeds=sorted({(tuple(r['source_support']),r['bad_source_walls'][0]) for r in bank['single_deletion_cutters']})
 alphabet=np.cos(np.pi/np.arange(2,61));hits=[];cutters=[];counts=dict(seeds=len(seeds),linear_supports=0,quadratic_cases=0,near_null_linear_cases=0,degenerate_cases_skipped=0,small_negative_discriminants_clamped=0,real_branches=0,norm_residual_rejections=0,future_branches=0,acute_branches=0)
 best=[]
 for S,bad in seeds:
  for released in S:
   shared=tuple(i for i in S if i!=released);C=cov[list(shared)+[bad]];v=null_space(C);assert v.shape==(19,1);v=v[:,0]
   Dlin=C.T@np.linalg.inv(C@C.T);counts['linear_supports']+=1
   for m,n in [(4,6),(6,6)]:
    a,b=2*np.cos(np.pi/m),2*np.cos(np.pi/n);D=(a*a+b*b)/2-2
    for anchor_label in [2,3,4,5,6]:
     rhs=np.zeros(18);rhs[-1]=-np.sqrt(2*q[bad])*np.cos(np.pi/anchor_label);p=Dlin@rhs
     aa=v@Q@v;bb=2*p@Q@v;cc=p@Q@p+D;disc=bb*bb-4*aa*cc;counts['quadratic_cases']+=1
     if abs(aa)<1e-10:
      counts['near_null_linear_cases']+=1
      if abs(bb)<1e-10:counts['degenerate_cases_skipped']+=1;continue
      roots=[(0,-cc/bb)]
     else:
      if disc<-1e-9:continue
      if disc<0:counts['small_negative_discriminants_clamped']+=1;disc=0
      roots=[(branch,(-bb+branch*np.sqrt(disc))/(2*aa)) for branch in [-1,1]]
     for branch,t in roots:
      y=p+t*v;counts['real_branches']+=1
      if abs(y@Q@y+D)>1e-7:counts['norm_residual_rejections']+=1;continue
      if y@Q@H>=0:continue
      counts['future_branches']+=1;g=cov@y
      if max(g)>1e-7:continue
      counts['acute_branches']+=1;normalized=g/np.sqrt(2*q);relations=[];fails=[];maxerr=0
      for i,z in enumerate(normalized):
       if i in shared:relations.append(2);continue
       if z<-1-2e-8:relations.append('dotted')
       elif abs(z+1)<2e-8:relations.append('parallel')
       elif abs(z)<2e-8:relations.append(2)
       else:
        j=int(np.argmin(abs(alphabet+z)));err=abs(alphabet[j]+z);relations.append(j+2)
        if err>2e-8:fails.append(i);maxerr=max(maxerr,err)
      esc=y@Q@E[:19]+(a+b)/2
      row=dict(seed_support=list(S),shared_source_support=list(shared),released=released,anchor=bad,anchor_label=anchor_label,support_labels=[m,n],branch=branch,
               source_projection=y.tolist(),projection_square=float(y@Q@y),escape_pairing=float(esc),cuts_escape=bool(esc>1e-7),failed_source_walls=fails,maximum_angle_error=float(maxerr),source_relations=relations)
      if esc>1e-7:best.append(row)
      if not fails:
       hits.append(row)
       if esc>1e-7:cutters.append(row)
 best.sort(key=lambda r:(len(r['failed_source_walls']),r['maximum_angle_error'],-r['escape_pairing']))
 paths=[SRC,BANK,f'scripts/{TAG}.py'];return dict(status='NUMERICAL_BOUNDED_SOURCE_EDGE_NONINDUCING_CAP_REPAIR',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=counts,accepted_count=len(hits),cutting_count=len(cutters),accepted=hits,cutters=cutters,best_cutting_near_misses=best[:12],scope='Numerical proposals only. Five anchor labels and two transverse label pairs over one-release17-wall supports derived from the bounded corner scout. Near-null leading coefficients are treated as linear, and bad norm-residual branches are rejected; counters record these numerical limitations. No exhaustive quadratic-branch, all-profile exclusion, or exact construction claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],d['accepted_count'],d['cutting_count'])
