#!/usr/bin/env python3
"""Bounded numerical discovery of caps joining recent H9 boundary regions.

Every accepted cap still requires a separate exact certificate. This program
does not use a finite-label cutoff for a nonexistence or acceptance theorem.
"""
import argparse,json,math
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_escape35_20261016 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_joint_regions_20261016'
PRIMARY='results/fresh_nonarithmetic_escape35_20261016.json'
CONE='results/fresh_nonarithmetic_escape35_20261016_cone.json'
COMMON=[0,1,2,5,6]
RECENT=[33,34,35]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def data():
 d=json.loads((ROOT/PRIMARY).read_text());c=json.loads((ROOT/CONE).read_text())
 for obj in[d,c]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 I=A.R.L.I
 def old(x):
  z=I()
  for a in x:z=z*A.R.THETA+I(F(str(a)))
  return z
 roots=[old(A.R.enc(q)).sqrt()for q in[A.B.P.D,A.B.D,A.D]]
 def tower(x,level=3):
  if not level:return old(x)
  return tower(x[0],level-1)+roots[level-1]*tower(x[1],level-1)
 def val(z):return float((z.l+z.h)/2)
 Q=np.array([[val(old(x))for x in row]for row in d['form']])
 raw=np.array([[val(tower(x))for x in row]for row in d['normals']])
 norms=np.array([val(tower(x))for x in d['norms']])
 h=np.array([val(old(x))for x in d['center']]);e=np.array([val(old(x))for x in d['strict_escape']['point']])
 ev,U=np.linalg.eigh(Q);L=np.sqrt(abs(ev))[:,None]*U.T;J=np.diag(np.sign(ev))
 N=raw@L.T;N*=np.sqrt(2/norms)[:,None]
 H=L@h;H/=math.sqrt(-H@J@H);E=L@e;E/=-E@J@H
 cov=N@J;G=cov@N.T;rays=[];residual=0.
 for ids in c['complete_census']['hyperideal_rays']:
  _,s,vh=np.linalg.svd(cov[ids]);v=vh[-1];v/=-v@J@H
  residual=max(residual,float(max(cov@v)),float(np.max(abs(cov[ids]@v))))
  assert v@J@v>1e-10 and max(cov@v)<1e-6
  rays.append(v)
 assert len(rays)==1068 and max(abs(np.diag(G)-2))<1e-6
 return N,cov,J,H,E,np.array(rays),dict(ray_reconstruction_max_residual=residual,norm_error=float(max(abs(np.diag(G)-2))))
def roots_line(cov,J,ids,targets):
 M=cov[ids];u,s,vh=np.linalg.svd(M,full_matrices=True)
 if min(s)<1e-9:return [],'rank_deficient'
 p=np.linalg.lstsq(M,np.array(targets),rcond=None)[0];k=vh[-1]
 aa=k@J@k;bb=2*p@J@k;cc=p@J@p-2
 if abs(aa)<1e-10:
  if abs(bb)<1e-10:return [],'degenerate'
  rr=[-cc/bb]
 else:
  dd=bb*bb-4*aa*cc
  if dd<-1e-8:return [],'no_real'
  if abs(dd)<1e-8:return [p-bb/(2*aa)*k],'double'
  rr=[(-bb+sig*math.sqrt(dd))/(2*aa)for sig in[-1,1]]
 return [p+x*k for x in rr],'two_or_linear'
def run(mode):
 N,cov,J,H,E,rays,diagnostic=data();stats=Counter();acute=[];good=[]
 if mode=='polar':
  for i,v in enumerate(rays):
   r=v*math.sqrt(2/(v@J@v));ps=cov@r
   if max(ps)>1e-6:stats['numerically_obtuse']+=1;continue
   finite=[];unresolved=[]
   for j,g in enumerate(ps):
    if g<-2-1e-7:continue
    if abs(g+2)<1e-7:unresolved.append([j,float(g),'near_parallel']);continue
    if abs(g)<1e-7:finite.append([j,2]);continue
    m=math.pi/math.acos(max(-1.,min(1.,-g/2)));mm=round(m)
    if mm>=2 and abs(m-mm)<1e-6:finite.append([j,mm])
    else:unresolved.append([j,float(g),float(m)])
   stats['screened']+=1;stats['unresolved_finite_contacts']+=len(unresolved)
   cuts=np.flatnonzero(rays@J@r>1e-7).tolist()
   row=dict(ray_index=i,finite=finite,unresolved_finite=unresolved,escape_pairing=float(r@J@E),cut_rays=cuts,cut_count=len(cuts),point=r.tolist())
   acute.append(row)
   if not unresolved:good.append(row);print('POLAR_PROPOSAL',{k:v for k,v in row.items()if k not in['point','cut_rays']},flush=True)
  good.sort(key=lambda x:-x['cut_count']);compat=[];mutual=[]
  for i,j in combinations(range(len(good)),2):
   g=float(np.array(good[i]['point'])@J@np.array(good[j]['point']));mutual.append(g)
   if g>1e-7:continue
   if g<-2-1e-7:compat.append([i,j,'dotted']);continue
   if abs(g)<1e-7:compat.append([i,j,2]);continue
   if abs(g+2)<1e-7:continue
   m=math.pi/math.acos(max(-1.,min(1.,-g/2)));mm=round(m)
   if mm>=2 and abs(m-mm)<1e-6:compat.append([i,j,mm])
  return dict(status='BOUNDED_NUMERICAL_POLAR_SCREEN_ONLY',sources={p:sha(p)for p in[f'scripts/{TAG}.py',PRIMARY,CONE]},mode=mode,stats=dict(stats),diagnostic=diagnostic,proposals=good,all_ray_screens=acute,compatible_proposal_pairs=compat,mutual_range=[min(mutual),max(mutual)]if mutual else None,scope='Screens all1068 saved hyperideal directions numerically. Unknown finite contacts remain unresolved; proposals require exactification. No all-cap or finite-volume conclusion.')
 # Stage0 fixes all three recent contacts right. Stage1 allows labels2..5
 # on those same three contacts. The one further anchor has labels2..12.
 profiles=[(2,2,2)]if mode=='right'else product(range(2,6),repeat=3)
 for labels in profiles:
  base=COMMON+RECENT;tg=[0.]*5+[-2*math.cos(math.pi/m)if m!=2 else 0. for m in labels]
  for anchor in range(36):
   if anchor in base:continue
   ids=base+[anchor]
   for label in range(2,13):
    target=-2*math.cos(math.pi/label)if label!=2 else 0.
    rr,kind=roots_line(cov,J,ids,tg+[target]);stats[kind]+=1
    for branch,r in enumerate(rr):
     stats['isolated_numeric_branches']+=1
     ps=cov@r
     if r@J@H>=-1e-8:stats['wrong_center']+=1;continue
     if max(ps)>1e-7:stats['obtuse']+=1;continue
     stats['acute']+=1
     finite=[];unresolved=[]
     for i,g in enumerate(ps):
      if g<-2-1e-7:continue
      if abs(g+2)<1e-7:unresolved.append([i,float(g),'near_parallel']);continue
      if abs(g)<1e-7:finite.append([i,2]);continue
      m=math.pi/math.acos(max(-1.,min(1.,-g/2)))
      mm=round(m)
      if mm>=2 and abs(m-mm)<1e-6:finite.append([i,mm])
      else:unresolved.append([i,float(g),float(m)])
     cuts=rays@J@r;cut=np.flatnonzero(cuts>1e-7).tolist()
     row=dict(recent_labels=list(labels),anchor=anchor,label=label,branch=branch,
      finite=finite,unresolved_finite=unresolved,escape_pairing=float(r@J@E),
      cut_rays=cut,cut_count=len(cut),maximum_old_pairing=float(max(ps)),point=r.tolist())
     acute.append(row)
     if not unresolved:
      good.append(row);print('NUMERIC_PROPOSAL',{k:v for k,v in row.items()if k not in['point','cut_rays']},flush=True)
  print('profile',list(labels),dict(stats),'proposals',len(good),flush=True)
 good.sort(key=lambda x:(-x['cut_count'],-x['escape_pairing']))
 paths=[f'scripts/{TAG}.py',PRIMARY,CONE]
 return dict(status='BOUNDED_NUMERICAL_DISCOVERY_ONLY',sources={p:sha(p)for p in paths},
  mode=mode,right_support=COMMON,recent_contacts=RECENT,recent_label_domain=[2]if mode=='right'else list(range(2,6)),anchor_label_domain=list(range(2,13)),stats=dict(stats),diagnostic=diagnostic,proposals=good,acute_branches=acute,
  scope='Numerical discovery with bounded prescribed labels. Every proposal needs exact all-wall Coxeter/facet/escape checks. Nonmatching unspecified finite contacts stay unresolved. No all-label exclusion or new exact chamber claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['right','finite','polar'],default='right');ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run(a.mode);p=ROOT/f'results/{TAG}_{a.mode}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved numerical replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['stats'],'proposals',len(d['proposals']))
