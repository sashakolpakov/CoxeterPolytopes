#!/usr/bin/env python3
"""Bounded noninducing transverse-triangle scout over H18 finite corners.
Numerical proposal generation only; finite contacts are not exact certificates.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_joint_route_20261015_scout'
SOURCE='results/fresh_h20_a1a1_h18_coupling_20261015.json'
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 def dec(x):return float(x['a'])+np.sqrt(2)*float(x['b']) if '/' not in x['a']+x['b'] else float(__import__('fractions').Fraction(x['a']))+np.sqrt(2)*float(__import__('fractions').Fraction(x['b']))
 Q=np.array([[dec(x) for x in row] for row in d['source_form']]);N=np.array([[dec(x) for x in row] for row in d['source_normals']]);H=np.array([dec(x) for x in d['source_center']]);E=np.array([dec(x) for x in d['strict_escape']])
 val,vec=np.linalg.eigh(Q);assert sum(val<0)==1
 transform=np.diag(np.sqrt(abs(val)))@vec.T
 metric=np.diag([-1]+[1]*18)
 Z=N@transform.T;h=transform@H
 if h[0]<0:transform[0]*=-1;Z=N@transform.T;h=transform@H
 h/=np.sqrt(-h@metric@h)
 raw=np.eye(19)[:,1:]+np.outer(h,h[1:])
 J=raw.T@metric@raw;L=np.linalg.cholesky(J);T=raw@np.linalg.inv(L.T)
 assert np.max(abs(T.T@metric@T-np.eye(18)))<1e-8
 A=Z@metric@T;b=-Z@metric@h;assert min(b)>0
 rng=np.random.default_rng(20261015);directions=list(rng.normal(size=(512,18)))
 directions.extend(A);directions.extend(-A)
 e=transform@E[:19];directions.extend([T.T@metric@e,-T.T@metric@e])
 vertices={};fail=0
 for c in directions:
  sol=linprog(c,A_ub=A,b_ub=b,bounds=[(None,None)]*18,method='highs')
  assert sol.success,sol.message
  active=tuple(np.flatnonzero(abs(A@sol.x-b)<1e-7).tolist())
  y=h+T@sol.x;q=y@metric@y
  assert q<1e-7,(q,active)
  if q>=-1e-7:continue
  assert len(active)==18,(active,q)
  vertices[active]=y/np.sqrt(-q)
 norm=np.einsum('ij,jk,ik->i',N,Q,N)
 alphabet=np.cos(np.pi/np.arange(2,61));proposals=[];best=[];one_bad=[];tested=0
 for active,v in sorted(vertices.items()):
  vp=Z@metric@v;assert max(vp)<1e-7
  for m in range(3,25):
   for n in range(m,25):
    if 1/m+1/n>=.5-1e-13:continue
    a,bv=2*np.cos(np.pi/m),2*np.cos(np.pi/n);D=(a*a+bv*bv)/2-2
    assert D>0
    normalized=np.sqrt(D)*vp/np.sqrt(2*norm);tested+=1
    labels=[];error=0;bad=[]
    for j,z in enumerate(normalized):
     if z < -1-2e-8:labels.append('dotted')
     elif abs(z+1)<2e-8:labels.append('parallel')
     elif abs(z)<2e-8:labels.append(2)
     elif z>0:bad.append(j);error=max(error,z);labels.append(None)
     else:
      k=int(np.argmin(abs(alphabet+z)));gap=abs(alphabet[k]+z);labels.append(k+2)
      if gap>2e-8:bad.append(j);error=max(error,gap)
    esc=np.sqrt(D)*(v@metric@e)+(a+bv)/2
    row=dict(source_support=list(active),support_labels=[m,n],residual_norm=-D,
             escape_pairing=float(esc),cuts_saved_escape=bool(esc>1e-7),maximum_angle_error=float(error),bad_source_walls=bad,
             proposed_source_relations=labels,unit_source_corner_in_physical_basis=np.linalg.solve(transform,v).tolist())
    if not bad:proposals.append(row)
    if esc>0:best.append(row)
    if esc>0 and len(bad)==1:one_bad.append(row)
 best.sort(key=lambda x:(len(x['bad_source_walls']),x['maximum_angle_error'],-x['escape_pairing']))
 paths=[SOURCE,f'scripts/{TAG}.py']
 return dict(status='NUMERICAL_BOUNDED_NONINDUCING_TRIANGLE_SCOUT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
             lp_queries=len(directions),distinct_finite_source_corners=len(vertices),label_range=[3,24],target_angle_alphabet_tested=[2,60],
             corner_label_tests=tested,accepted_proposals=proposals,accepted_count=len(proposals),accepted_cutters=sum(p['cuts_saved_escape'] for p in proposals),
             best_cutting_near_misses=best[:12],single_deletion_cutters=one_bad,
             scope='Numerical bounded source-corner and finite-label scout. No source vertex census, all-label exclusion, exact angle acceptance, or finite-volume claim. Proposed caps are noninducing: source projection future timelike and transverse section a compact hyperbolic triangle.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print({k:d[k] for k in ['status','lp_queries','distinct_finite_source_corners','corner_label_tests','accepted_count','accepted_cutters']})
