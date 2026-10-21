#!/usr/bin/env python3
"""Bounded noninducing two-source-anchor search; one support may be dotted."""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_joint_route_20261015_two_anchor';SRC='results/fresh_h20_a1a1_h18_coupling_20261015.json';BANK='results/fresh_h20_joint_route_20261015_scout.json'
def run():
 d=json.loads((ROOT/SRC).read_text());bank=json.loads((ROOT/BANK).read_text());dec=lambda x:float(Fraction(x['a']))+np.sqrt(2)*float(Fraction(x['b']))
 Q=np.array([[dec(x)for x in row]for row in d['source_form']]);N=np.array([[dec(x)for x in row]for row in d['source_normals']]);H=np.array([dec(x)for x in d['source_center']]);E=np.array([dec(x)for x in d['strict_escape']]);q=np.einsum('ij,jk,ik->i',N,Q,N);cov=N@Q
 seeds=sorted({(tuple(r['source_support']),r['bad_source_walls'][0])for r in bank['single_deletion_cutters']});alph=np.cos(np.pi/np.arange(2,61));hits=[];best=[];counts=dict(supports=0,source_anchor_tests=0,future_timelike=0,source_acute=0,source_coxeter=0,support_tests=0)
 def gate(z):
  if z<-1-2e-8:return'dotted'
  if abs(z+1)<2e-8:return'parallel'
  if abs(z)<2e-8:return 2
  if z>0:return None
  k=int(np.argmin(abs(alph+z)));return k+2 if abs(alph[k]+z)<2e-8 else None
 for S,bad in seeds:
  for released in S:
   shared=tuple(i for i in S if i!=released);C=cov[list(shared)+[bad,released]];J=np.linalg.inv(C);counts['supports']+=1
   for ma in [3,4,6]:
    for mb in [3,4,6]:
     rhs=np.zeros(19);rhs[-2]=-np.sqrt(2*q[bad])*np.cos(np.pi/ma);rhs[-1]=-np.sqrt(2*q[released])*np.cos(np.pi/mb);y=J@rhs;D=-(y@Q@y);counts['source_anchor_tests']+=1
     if D<=1e-7 or y@Q@H>=0:continue
     counts['future_timelike']+=1;g=cov@y
     if max(g)>1e-7:continue
     counts['source_acute']+=1
     rel=[2 if i in shared else gate(z) for i,z in enumerate(g/np.sqrt(2*q))];badrows=[i for i,x in enumerate(rel)if x is None]
     if not badrows:counts['source_coxeter']+=1
     for m in [2,3,4,6]:
      a=2*np.cos(np.pi/m);b2=2*D+4-a*a
      if b2<=0:continue
      b=np.sqrt(b2);support_relation=gate(-b/2);counts['support_tests']+=1
      if support_relation is None:continue
      esc=y@Q@E[:19]+(a+b)/2
      row=dict(shared_source_support=list(shared),source_anchors=[bad,released],source_anchor_labels=[ma,mb],source_projection=y.tolist(),projection_norm=-float(D),support_contacts=[float(a),float(b)],support_relations=[m,support_relation],source_relations=rel,failed_source_walls=badrows,escape_pairing=float(esc),cuts_escape=bool(esc>1e-7))
      if not badrows:hits.append(row)
      if esc>1e-7:best.append(row)
 best.sort(key=lambda r:(len(r['failed_source_walls']),-r['escape_pairing']))
 paths=[SRC,BANK,f'scripts/{TAG}.py'];return dict(status='NUMERICAL_BOUNDED_NONINDUCING_TWO_ANCHOR_TEST',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=counts,accepted_count=len(hits),cutters=[r for r in hits if r['cuts_escape']],accepted=hits,best_cutting_near_misses=best[:12],scope='Numerical source17-right/two-finite-anchor search. One support contact2/3/4/6 and the other unrestricted dotted or tested finite2..60. No all-label assertion; no exact construction or volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],'accepted',d['accepted_count'],'cutters',len(d['cutters']))
