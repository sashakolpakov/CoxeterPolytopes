#!/usr/bin/env python3
"""Solve isolated old-facet finite events; classify every other old168 gate."""
import argparse,hashlib,json,math,itertools
import numpy as np
import sympy as sp
import fresh_h20_e8_partner_endpoints_20261017 as P
import fresh_h20_e8_partner_endpoints_20261017_bridges as B
import fresh_h20_e8_partner_endpoints_20261017_coordinated as K
ROOT=P.ROOT;TAG='fresh_h20_e8_endpoint_events_20261017';INPUT='results/'+P.TAG+'.json'
def classify(co):
 labels=[]
 for v in co:
  if v>1+1e-7:labels.append('dotted')
  elif abs(v)<1e-7:labels.append(2)
  elif abs(v-1)<1e-7:labels.append('parallel')
  elif 0<v<1:
   z=math.pi/math.acos(v);lab=round(z)
   if lab<3 or abs(v-math.cos(math.pi/lab))>1e-7:return None
   labels.append(lab)
  else:return None
 return labels

def run():
 D=P.load(INPUT);G,Y,CT,J,q,H,live,verts=P.geometry();C=P.load('results/fresh_h20_e8_cusp_completion_20261017.json');F=P.load('results/fresh_h20_e8_f4_cusps_20261017.json');ids=C['core_ids'];M=np.array([[float(sp.sympify(x))for x in row]for row in C['core_gram']]);MI=np.linalg.inv(M);added=[]
 for cu in C['cusps']:
  nu=np.array([float(sp.sympify(x))for x in cu['null_coefficients']])
  for row in cu['base_new_root_coefficients']:added.append(np.array([float(sp.sympify(x))for x in row])+F['selected_common_height']*nu)
 A=np.array(added);Y=np.r_[Y,A@Y[ids]];CT=np.r_[CT,A@CT[ids]];q=np.r_[q,[2.]*6];oldfacets=G['sides']+[158];m=G['m'];stats=dict(equations=0,positive_roots=0,physical_coherent=0,all_acute=0,all_Coxeter=0);out=[];near=[]
 for ci,c in enumerate(D['endpoint_censuses']):
  missing=next(j for j in(159,160,161)if j not in c['edge']);v=-MI[:,ids.index(missing)];V=v@Y[ids];assert V@m@V<0
  for row in c['profiles']:
   for br in c['families'][row['id']]['branches']:
    a=np.array(br['positive_k_coefficient']);b=np.array(br['inverse_k_coefficient']);p=np.array(br['particular']);dt=np.array(br['T_pairings']);offset=CT@J@dt;aa=a@m@Y.T;bb=b@m@Y.T;cc=p@m@Y.T+offset
    for j in G['sides']:
     isolated=max(abs(dt))<1e-9 and max(abs(CT[j]))<1e-8
     labels=range(2,65)if isolated else range(2,6)
     for lab in labels:
      target=-math.sqrt(2*q[j])*(0 if lab==2 else math.cos(math.pi/lab));stats['equations']+=1
      for k in B.roots(aa[j],cc[j]-target,bb[j]):
       stats['positive_roots']+=1;r=p+a*k+b/k
       if r@m@V>=-1e-9:continue
       stats['physical_coherent']+=1;ps=r@m@Y.T+offset;co=-ps/np.sqrt(2*q)
       if min(co)<-1e-6:continue
       stats['all_acute']+=1
       native=classify(co)
       if native is None:
        errors=[]
        for z in co:
         if 0<z<1:
          ell=max(2,round(math.pi/math.acos(z)));errors.append(abs(z-math.cos(math.pi/ell)))
        near.append(dict(edge=ci,profile=row['id'],anchor=j,label=lab,maximum_angle_error=max(errors,default=0.)));continue
       if j!=158:
        h=2-dt@J@dt;g=ps[j]-offset[j]
        if h*(Y[j]@m@Y[j])-g*g<=1e-8:continue
       if any(max(abs(np.array(z['residual4'])-r))<1e-6 and max(abs(np.array(z['T_pairings'])-dt))<1e-6 for z in out):continue
       stats['all_Coxeter']+=1
       finiteold=[u for u in oldfacets if u not in c['edge'] and (2-dt@J@dt)*(Y[u]@m@Y[u])-(ps[u]-offset[u])**2>1e-7]
       face=K.poly(np.r_[Y,[r]],G['T'])
       out.append(dict(edge_index=ci,edge=c['edge'],profile=row['id'],support_contacts=row['contacts'],branch=br['branch'],event_anchor=j,event_label=lab,k=k,residual4=r.tolist(),T_pairings=dt.tolist(),all168_labels=native,extra_intrinsic_old_neighbors=finiteold,face=face))
  print('EVENT_EDGE',ci,stats,'retained',len(out),flush=True)
 near.sort(key=lambda z:z['maximum_angle_error']);paths=[f'scripts/{TAG}.py',INPUT,P.CUSP if hasattr(P,'CUSP')else 'results/fresh_h20_e8_cusp_completion_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json',f'scripts/{P.TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{K.TAG}.py']
 return dict(status='BOUNDED_FINITE_OLD_FACET_ENDPOINT_EVENT_SOLVER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},statistics=stats,actual_old_H3_facets=oldfacets,candidates=out,best_rejected_acute_events=near[:100],scope='All159 fixed profiles, both conic branches. Each event sets a finite contact to an actual old H3 facet, then every168 old pair is allowed finite/right/parallel/dotted and explicitly classified numerically. No fixed old center is required; the opposite branch is removed by its positive product with the future old finite corner. Nonisolated finite corner bonds have labels<=5; isolated I2 events use bounded labels<=64; ambient158 is freely classified after every other event. The singleton151 profiles and unbounded isolated labels remain open; tolerance-based angle equalities need exact reconstruction.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
