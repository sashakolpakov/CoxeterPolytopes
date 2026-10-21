#!/usr/bin/env python3
"""Free the ambient158 bond in fixed endpoint conics; bounded angle scout."""
import argparse,hashlib,json,math,itertools
import numpy as np
import sympy as sp
import fresh_h20_e8_partner_endpoints_20261017 as P
import fresh_h20_e8_partner_endpoints_20261017_bridges as B
ROOT=P.ROOT;TAG='fresh_h20_e8_endpoint_belt_20261017_finite158';INPUT='results/'+P.TAG+'.json'
def run():
 D=P.load(INPUT);G,Y,CT,J,q,H,live,verts=P.geometry();m=G['m'];domains=[];out=[];near=[]
 C=P.load('results/fresh_h20_e8_cusp_completion_20261017.json');K=P.load('results/fresh_h20_e8_f4_cusps_20261017.json');ids=C['core_ids'];added=[]
 for cu in C['cusps']:
  nu=np.array([float(sp.sympify(x))for x in cu['null_coefficients']])
  for row in cu['base_new_root_coefficients']:added.append(np.array([float(sp.sympify(x))for x in row])+K['selected_common_height']*nu)
 AA=np.array(added);Y=np.r_[Y,AA@Y[ids]];CT=np.r_[CT,AA@CT[ids]];q=np.r_[q,[2.]*6];live=live+list(range(162,168));completion_rejections=[]
 for ci,c in enumerate(D['endpoint_censuses']):
  for row in c['profiles']:
   for br in c['families'][row['id']]['branches']:
    a=np.array(br['positive_k_coefficient']);b=np.array(br['inverse_k_coefficient']);p=np.array(br['particular']);dt=np.array(br['T_pairings']);offset=CT@J@dt
    aa=a@m@Y.T;bb=b@m@Y.T;cc=p@m@Y.T+offset
    ints=[(0.,math.inf)]
    for co in br['constraints']:
     if co['root']!=158:ints=P.intersections(ints,P.quadratic_intervals(*co['quadratic']))
    ints=P.intersections(ints,P.quadratic_intervals(*br['center_quadratic']))
    ints=P.intersections(ints,P.quadratic_intervals(aa[158],cc[158],bb[158]))
    ints=P.intersections(ints,P.quadratic_intervals(-aa[158],-cc[158]-2,-bb[158]))
    if not ints:continue
    ranges=[B.interval_range(aa[158],bb[158],cc[158],lo,hi,2)for lo,hi in ints]
    domains.append(dict(edge=c['edge'],profile=row['id'],branch=br['branch'],intervals=[[lo,None if math.isinf(hi)else hi]for lo,hi in ints],cosine_ranges=[list(z)for z in ranges]))
    for label in list(range(2,65))+['parallel']:
     co=1 if label=='parallel'else 0 if label==2 else math.cos(math.pi/label)
     for k in B.roots(aa[158],cc[158]+2*co,bb[158]):
      if not any(lo-1e-8*(1+lo)<=k and(math.isinf(hi)or k<=hi+1e-8*(1+hi))for lo,hi in ints):continue
      r=p+k*a+b/k;ps=r@m@Y.T+offset;other=[j for j in range(162)if j not in c['support']and j!=158]
      minimum=float(min(-ps[other]/np.sqrt(2*q[other])))
      if minimum<1-1e-7:continue
      last=-ps[162:]/2;labs=[]
      for v in last:
       if v>1+1e-7:labs.append('dotted')
       elif abs(v)<1e-7:labs.append(2)
       elif abs(v-1)<1e-7:labs.append('parallel')
       elif 0<v<1:
        z=math.pi/math.acos(v);ell=round(z)
        labs.append(ell if ell>=3 and abs(z-ell)<1e-5 else 'gap')
       else:labs.append('obtuse')
      if any(z in('gap','obtuse')for z in labs):
       completion_rejections.append(dict(edge_index=ci,profile=row['id'],label158=label,completion_cosines=last.tolist(),classifications=labs));continue
      rec=dict(completion_labels=labs,completion_cosines=last.tolist(),edge_index=ci,edge=c['edge'],profile=row['id'],contacts=row['contacts'],branch=br['branch'],label158=label,k=k,residual4=r.tolist(),T_pairings=dt.tolist(),minimum_other_cosine=minimum,face=P.poly(Y,live,H,r))
      out.append(rec)
 print('FINITE158 domains',len(domains),'caps',len(out),flush=True)
 pairrows=[];dotted_pairs=0
 if out:
  RR=np.array([r['residual4']for r in out]);DD=np.array([r['T_pairings']for r in out]);ci=np.array([r['edge_index']for r in out]);ar=np.arange(len(out))
  for lo in range(0,len(out),128):
   hi=min(lo+128,len(out));co=-(RR[lo:hi]@m@RR.T+DD[lo:hi]@J@DD.T)/2
   mask=(ar[None,:]<ar[lo:hi,None])&(ci[None,:]!=ci[lo:hi,None]);dotted_pairs+=int(np.sum(mask&(co>1+1e-7)))
   rows,cols=np.where(mask&(co>=-1e-7)&(co<=1+1e-7))
   for ri,j in zip(rows,cols):
    i=lo+int(ri);v=float(co[ri,j])
    if abs(v)<1e-7:lab=2
    elif abs(v-1)<1e-7:lab='parallel'
    else:
     z=math.pi/math.acos(max(-1,min(1,v)));lab=round(z)
     if lab<3 or abs(z-lab)>1e-5:continue
    pairrows.append(dict(i=i,j=int(j),label=lab,cosine=v))
 src=[f'scripts/{TAG}.py',INPUT,f'scripts/{P.TAG}.py',f'scripts/{B.TAG}.py','results/fresh_h20_e8_cusp_completion_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json']
 return dict(status='BOUNDED_AMBIENT158_FINITE_ENDPOINT_DIAGNOSTIC',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},domains=domains,completion_rejections=completion_rejections,candidates=out,finite_compatible_different_edge_pairs=pairrows,dotted_different_edge_pairs=dotted_pairs,scope='All159 fixed endpoint profiles, both conic branches; every unprescribed old pair except158 required dotted. All six completion-wall gates are separately classified; full168 replacement polyhedra use all retained constraints. Interval endpoints numerical; ambient158 labels2..64 plusparallel sampled exactly in formula. This relaxes the prior all-unprescribed-dotted bank, but is not complete over unbounded labels or multiple finite coincidences. No exact cap/finite-volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
