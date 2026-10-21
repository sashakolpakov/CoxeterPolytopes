#!/usr/bin/env python3
"""Numerical all-height screen on the 16 fresh open spherical edges.

All norm-two parameters on each of two hyperbola branches are intersected
with retained-wall acuteness and the current target-cut condition. This is
a scout, not an exact interval/no-go certificate. The infinitely many
isolated I2(m), m>=6, profiles are explicitly recorded as not screened.
"""
import itertools,json,math,time
from fractions import Fraction
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE='data/enumerations/fresh_h9_component_bridge_20261013';OUT=BASE+'_continuous_scout.json'

def negative_intervals(a,b,c):
 scale=max(abs(a),abs(b),abs(c),1e-100);a/=scale;b/=scale;c/=scale
 if abs(a)<1e-11:
  if abs(b)<1e-11:return[(0,math.inf)]if c<=1e-11 else[]
  root=-c/b
  return[(max(0,root),math.inf)]if b<0 else([(0,root)]if root>=0 else[])
 disc=b*b-4*a*c
 if disc< -1e-11:return[(0,math.inf)]if a<0 else[]
 disc=max(0,disc);x,y=sorted(((-b-math.sqrt(disc))/(2*a),(-b+math.sqrt(disc))/(2*a)))
 pieces=[(x,y)]if a>0 else[(-math.inf,x),(y,math.inf)]
 return[(max(0,x),y)for x,y in pieces if y>=0 and y>=max(0,x)]

def intersect(xs,ys):
 return[(max(a,c),min(b,d))for a,b in xs for c,d in ys if min(b,d)>=max(a,c)]

def roots(a,b,c):
 scale=max(abs(a),abs(b),abs(c),1e-100);a/=scale;b/=scale;c/=scale
 if abs(a)<1e-12:return[-c/b]if abs(b)>1e-12 else[]
 d=b*b-4*a*c
 return[]if d<0 else[(-b-math.sqrt(d))/(2*a),(-b+math.sqrt(d))/(2*a)]

def run():
 d=json.loads((ROOT/(BASE+'_faces_scout.json')).read_text());cert=json.loads((ROOT/(BASE+'_verify.json')).read_text());cap=json.loads((ROOT/(BASE+'_caps.json')).read_text())
 G=np.array(d['gram']);R=np.vstack([d['normals'],[[sum(map(int,x))/(2*10**80)for x in r]for r in cert['normal_enclosures']]])
 names=d['names']+[r['name']for r in cert['recipes']];center=np.array(d['center']);h=np.array(d['h']);phi=(1+np.sqrt(5))/2
 target=np.array([float(Fraction(x['a']))+phi*float(Fraction(x['b']))for x in cap['escape']]);w=-np.linalg.solve(G,target)
 edges=[r['edge']for r in cert['new_vertex_frontier']if len(r['end_walls'])==1];assert len(edges)==16
 fixed_labels=np.array([2,3,4,5]);alphabet=np.array([0.,-1.,-np.sqrt(2),-phi]);profiles=np.array(list(itertools.product(range(4),repeat=8)),dtype=int);P=alphabet[profiles]
 anchor_labels=np.array([2,3,4,5,6,8,10,12]);anchor_values=-2*np.cos(np.pi/anchor_labels);cov=R@G
 stats={'edges':16,'finite_profile_vectors':0,'branches':0,'acute_target_intervals':0,'anchor_proposals':0,'fully_compatible':0};families=[];proposals=[];seen=set();tails=[];start=time.time()
 for ei,ss in enumerate(edges):
  ids=[names.index(n)for n in ss];A=R[ids];M=A@G@A.T;Mi=np.linalg.inv(M);rem=2-np.einsum('bi,ij,bj->b',P,Mi,P);valid=np.flatnonzero(rem>1e-8)
  isol=[i for i in range(8)if np.max(abs(np.delete(M[i],i)))<1e-7];tails.extend({'edge':ei,'isolated_wall':ss[i],'uncovered_labels':'all integer m>=6','coupling_magnitude_range':'[sqrt(3),2)'}for i in isol)
  proj=lambda v:v-(A@G@v)@Mi@A
  tv=proj(center);qt=tv@G@tv;assert qt<0;tv/=np.sqrt(-qt);ev=proj(w);ev+=(ev@G@tv)*tv;qe=ev@G@ev;assert qe>0;ev/=np.sqrt(qe)
  exterior=[j for j in range(len(R))if j not in ids]
  for pi in valid:
   p=P[pi];base=p@Mi@A;rho=np.sqrt(rem[pi]);stats['finite_profile_vectors']+=1
   for branch in(-1,1):
    aa=rho*(tv+branch*ev)/2;bb=rho*(-tv+branch*ev)/2;stats['branches']+=1
    # Conditions r.y<=0 (cuts target), -r.h<=0 (correct orientation),
    # then all old pairings <=0, with the exact prescribed edge skipped.
    ac=np.r_[aa@target,-aa@h,cov[exterior]@aa];bc=np.r_[base@target,-base@h,cov[exterior]@base];cc=np.r_[bb@target,-bb@h,cov[exterior]@bb]
    ranges=[(0,math.inf)]
    for a,b,c in zip(ac,bc,cc):
     ranges=intersect(ranges,negative_intervals(a,b,c))
     if not ranges:break
    ranges=[(lo,hi)for lo,hi in ranges if hi>max(lo,0)+1e-9]
    if not ranges:continue
    family={'edge_index':ei,'edge':ss,'labels':fixed_labels[profiles[pi]].tolist(),'branch':branch,'intervals':[[float(lo),float(hi)if np.isfinite(hi)else None]for lo,hi in ranges]}
    fi=len(families);families.append(family);stats['acute_target_intervals']+=len(ranges)
    pa,pc,pb=cov@aa,cov@base,cov@bb
    for j in exterior:
     for lab,val in zip(anchor_labels,anchor_values):
      for z in roots(pa[j],pc[j]-val,pb[j]):
       if not z>0 or not any(lo-1e-9<=z<=hi+1e-9 for lo,hi in ranges):continue
       r=base+z*aa+bb/z
       if r@target>=-1e-8 or r@h<=1e-8:continue
       stats['anchor_proposals']+=1;pp=cov@r
       ok=(pp<-2-1e-7)|(np.min(abs(pp[:,None]-anchor_values),axis=1)<1e-6)
       if not all(ok):continue
       stats['fully_compatible']+=1;key=tuple(np.round(r,7))
       if key in seen:continue
       seen.add(key);proposals.append({'family':fi,'anchor':names[j],'anchor_label':int(lab),'parameter':float(z),'normal':r.tolist(),'target_cut':-float(r@target),'finite_neighbors':[names[k]for k in np.flatnonzero(pp>-2+1e-6)]})
  print('EDGE',ei,'profiles',len(valid),'families',len(families),'caps',len(proposals),'seconds',round(time.time()-start,1),flush=True)
 out={'status':'NUMERICAL_ALL_HEIGHT_FINITE_PROFILE_SCREEN','source_caps':BASE+'_caps.json','stats':stats,'finite_profile_labels':fixed_labels.tolist(),'anchor_labels':anchor_labels.tolist(),'families':families,'compatible_proposals':proposals,'uncovered_isolated_I2_tails':tails,
      'scope':'Continuous real-height acuteness/cutting intervals are numerically intersected for every PD rank9 extension profile with labels2,3,4,5 of16 archived open spherical8 edges. No interval or rejection is an exact certificate. Higher labels can occur only in an isolated I2(m) factor; all such m>=6 tails are explicitly left open. Fixed anchor probes are finite and not exhaustive. Positive proposals still require exact angle/affine/volume checks.'}
 (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print('DONE',stats,'unique',len(proposals),'tails',len(tails),flush=True)

if __name__=='__main__':run()
