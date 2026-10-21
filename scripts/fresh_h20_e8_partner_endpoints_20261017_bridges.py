#!/usr/bin/env python3
"""Numerical old-facet coincidences of the complete fixed endpoint profiles."""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import fresh_h20_e8_partner_endpoints_20261017 as P
ROOT=P.ROOT;TAG=P.TAG+'_bridges';INPUT='results/'+P.TAG+'.json'
def roots(a,b,c):
 scale=max(abs(a),abs(b),abs(c),1.);a/=scale;b/=scale;c/=scale
 if abs(a)<1e-12:return [-c/b]if abs(b)>1e-12 and -c/b>0 else[]
 disc=b*b-4*a*c
 if disc<0:return []
 return sorted({float(r)for r in [(-b-math.sqrt(disc))/(2*a),(-b+math.sqrt(disc))/(2*a)]if r>0})
def interval_range(a,b,c,lo,hi,scale):
 ks=[lo,hi]
 if a*b>0:
  k=math.sqrt(b/a)
  if lo<k<hi:ks.append(k)
 vals=[]
 for k in ks:
  if k==0:val=math.copysign(math.inf,-b)if b else-c/scale
  elif math.isinf(k):val=math.copysign(math.inf,-a)if a else-c/scale
  else:val=-(a*k+b/k+c)/scale
  vals.append(val)
 return min(vals),max(vals)
def run():
 D=P.load(INPUT);G,Y,CT,JT,q,H,live,verts=P.geometry();out=[];domains=[];tails=[];acute_domains=[]
 for ci,cen in enumerate(D['endpoint_censuses']):
  W=cen['support']
  for ff in cen['families']:
   row=cen['profiles'][ff['profile']]
   for br in ff['branches']:
    va=np.array(br['positive_k_coefficient']);vb=np.array(br['inverse_k_coefficient']);p=np.array(br['particular']);dt=-np.array(br['T_pairings'])
    offset=-CT@JT@dt;a=va@G['m']@Y.T;b=vb@G['m']@Y.T;c=p@G['m']@Y.T+offset
    hT=2-dt@JT@dt
    acute=[(0.,math.inf)]
    for root in range(162):
     if root not in W:acute=P.intersections(acute,P.quadratic_intervals(a[root],c[root],b[root]))
    acute=P.intersections(acute,P.quadratic_intervals(*br['center_quadratic']))
    for j in G['sides']:
     # A finite intrinsic intersection additionally requires the projected
     # two-root Gram to be PD, not merely a finite ambient pair.
     threshold=offset[j]-math.sqrt(hT*G['qt'][j])
     aj=P.intersections(acute,P.quadratic_intervals(-a[j],-c[j]+threshold,-b[j]))
     aj=P.intersections(aj,P.quadratic_intervals(-a[j],-c[j]-math.sqrt(2*q[j]),-b[j]))
     if aj:acute_domains.append(dict(edge=cen['edge'],profile=ff['profile'],branch=br['branch'],old_facet=j,intervals=[[lo,None if math.isinf(hi)else hi]for lo,hi in aj],scope='All other old gates relaxed to acuteness; projected pair positive definite in interval interiors. Not Coxeter feasibility.'))
     ints=[(0.,math.inf)]
     for cc in br['constraints']:
      if cc['root']!=j:ints=P.intersections(ints,P.quadratic_intervals(*cc['quadratic']))
     ints=P.intersections(ints,P.quadratic_intervals(*br['center_quadratic']))
     scale=math.sqrt(2*q[j])
     # The designated old contact may be any nondotted magnitude [0,scale].
     ints=P.intersections(ints,P.quadratic_intervals(a[j],c[j],b[j]))
     ints=P.intersections(ints,P.quadratic_intervals(-a[j],-c[j]-scale,-b[j]))
     if not ints:continue
     isolated=not np.any(abs(dt)>1e-10) and not np.any(abs(CT[j])>1e-7)
     ran=[interval_range(a[j],b[j],c[j],lo,hi,scale)for lo,hi in ints]
     upper=min(1.,max(z[1]for z in ran));lower=max(0.,min(z[0]for z in ran))
     if isolated:
      maxlabel=int(math.ceil(math.pi/math.acos(upper)))+1 if upper<1-1e-9 else 100
      if upper>=1-1e-9:tails.append(dict(edge=cen['edge'],profile=ff['profile'],branch=br['branch'],old_facet=j,unbounded_label_tail='m>100 not enumerated',intervals=[[lo,None if math.isinf(hi)else hi]for lo,hi in ints]))
     else:maxlabel=5
     dom=dict(edge=cen['edge'],profile=ff['profile'],branch=br['branch'],old_facet=j,isolated_I2=isolated,cosine_range=[lower,upper],intervals=[[lo,None if math.isinf(hi)else hi]for lo,hi in ints],tested_maximum_label=maxlabel)
     domains.append(dom)
     for label in range(2,maxlabel+1):
      co=0. if label==2 else math.cos(math.pi/label)
      if co<lower-1e-8 or co>upper+1e-8:continue
      target=-scale*co;projected=target-offset[j]
      if hT*G['qt'][j]-projected**2<=1e-8:continue
      for k in roots(a[j],c[j]+scale*co,b[j]):
       if not any(lo-1e-7*(1+lo)<=k and (math.isinf(hi)or k<=hi+1e-7*(1+hi))for lo,hi in ints):continue
       r=p+k*va+vb/k;cos=-(r@G['m']@Y.T+offset)/np.sqrt(2*q)
       other=[i for i in range(162)if i not in W and i!=j]
       if min(cos[other])<1-1e-6:continue
       if r@G['m']@H>=0:continue
       out.append(dict(edge=cen['edge'],profile=ff['profile'],branch=br['branch'],old_facet=j,label=label,k=k,residual4=r.tolist(),T_pairings=(-dt).tolist(),minimum_other_cosine=float(min(cos[other])),finite_edge_determinant=float(hT*G['qt'][j]-projected**2),face=P.poly(Y,live,H,r)))
 print('DOMAINS',len(domains),'CANDIDATES',len(out),'TAILS',len(tails),'ACUTE_RELAXATIONS',len(acute_domains),flush=True)
 src=[f'scripts/{TAG}.py',INPUT,f'scripts/{P.TAG}.py']
 return dict(status='NUMERICAL_ENDPOINT_TO_OLD_FACET_COINCIDENCE_DIAGNOSTIC',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},domains=domains,unbounded_tails=tails,candidates=out,acute_relaxations=acute_domains,scope='Complete numerical quadratic interval checks for fixed endpoint profiles with one designated nondotted old facet and every other unprescribed contact dotted. Nonisolated finite bonds have labels<=5 by finite-diagram classification. Isolated I2 tails are explicitly marked if not bounded by interval range. The three parametric endpoint151 profiles and multiple simultaneous old finite contacts are not exhausted. All-acute relaxations remain open and explicitly prevent a global obstruction claim. No exact-wall or global exclusion claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,P.load(str(out.relative_to(ROOT))))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
