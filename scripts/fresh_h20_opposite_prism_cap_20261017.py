#!/usr/bin/env python3
"""Complete fixed-profile intersection, bounded singleton scout for one prism cap.

One wall is required to make all three partner-edge endpoints spherical.
The profile intersection is exact; current ambient gate discovery is numerical.
"""
import argparse,hashlib,itertools,json,math
from collections import Counter
import numpy as np
import sympy as s
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection
import fresh_h20_e8_partner_endpoints_20261017 as P
import fresh_h20_e8_partner_endpoints_20261017_coordinated as K
import fresh_h20_e8_endpoint_events_20261017 as EV
ROOT=P.ROOT;TAG='fresh_h20_opposite_prism_cap_20261017';END='results/'+P.TAG+'.json'

def profiles(D):
 rows=[z['profiles']for z in D['endpoint_censuses']];Ws=[z['support']for z in D['endpoint_censuses']]
 by3={tuple(tuple(z)for z in r['contacts']):r for r in rows[2]};out=[]
 for a in rows[0]:
  da=dict(a['contacts'])
  for b in rows[1]:
   db=dict(b['contacts'])
   if any(da.get(j,2)!=db.get(j,2)for j in set(Ws[0])&set(Ws[1])):continue
   d={**da,**db};key=tuple(sorted((j,m)for j,m in d.items()if j in Ws[2]))
   if key not in by3:continue
   c=by3[key];out.append(dict(id=len(out),endpoint_profiles=[a['id'],b['id'],c['id']],contacts=sorted([j,m]for j,m in d.items())))
 return out

def geometry():
 G,Y,CT,J,q,H,live,verts=P.geometry();C=P.load('results/fresh_h20_e8_cusp_completion_20261017.json');F=P.load('results/fresh_h20_e8_f4_cusps_20261017.json');ids=C['core_ids'];M=np.array(s.Matrix(C['core_gram']),float);MI=np.linalg.inv(M)
 for cu in C['cusps']:
  nu=np.array(s.Matrix(cu['null_coefficients']),float).ravel()
  for b in cu['base_new_root_coefficients']:
   v=np.array(s.Matrix(b),float).ravel()+F['selected_common_height']*nu
   Y=np.r_[Y,[v@Y[ids]]];CT=np.r_[CT,[v@CT[ids]]];q=np.r_[q,2]
 V=[]
 for k in (159,160,161):
  v=-MI[:,ids.index(k)];vv=v@Y[ids];vv/=math.sqrt(-(vv@G['m']@vv));assert vv[0]>0;V.append(vv)
 return G,Y,CT,J,q,np.array(V),C

def run():
 D=P.load(END);rows=profiles(D);G,Y,CT,J,q,Vs,C=geometry();m=G['m'];T=G['T'];U=Y[[159,160,161]];gram=U@m@U.T;inv=np.linalg.inv(gram)
 ker=null_space(U@m)[:,0];kn=ker@m@ker;assert kn>0
 projections=q-np.einsum('ij,jk,ik->i',CT,J,CT);positive=[i for i,z in enumerate(projections)if z>1e-7 and i not in T]
 # In the coherent acute chamber, positivity of T+j is sufficient for an
 # actual target facet by projection of the common center.
 assert positive==[95,142,145,150,152,156,157,158,159,160,161],positive
 labels=list(range(3,65))+[96,128,256,1024]
 inputs=[dict(r,kind='fixed')for r in rows]+[dict(id=lab,kind='singleton151',contacts=[[151,lab]],endpoint_profiles=None)for lab in labels]
 stats=Counter();candidates=[];records=[]
 for row in inputs:
  dt=np.zeros(17);dp=np.zeros(3)
  for j,lab in row['contacts']:
   if j in T:dt[T.index(j)]=2*math.cos(math.pi/lab)
   else:dp[[159,160,161].index(j)]=2*math.cos(math.pi/lab)
  offset=-CT@J@dt;f=-dp-offset[[159,160,161]];p=(inv@f)@U;h=2-dt@J@dt;rad=(h-p@m@p)/kn
  stats['patterns']+=1;rec=dict(profile=row,projected_square=float(h),kernel_square=float(kn),norm_radicand=float(rad),branches=[])
  if rad<-1e-9:stats['negative_radicand']+=1;records.append(rec);continue
  for sign in(-1,1):
   r=p+sign*math.sqrt(max(0.,rad))*ker;co=-(r@m@Y.T+offset)/np.sqrt(2*q);cv=r@m@Vs.T
   br=dict(sign=sign,residual4=r.tolist(),old_vertex_products=cv.tolist(),minimum_cosine=float(min(co)))
   stats['real_branches']+=1
   if max(cv)>=-1e-9:br['rejection']='old finite vertex not strict';stats['wrong_orientation']+=1
   elif min(co)<-1e-7:br['rejection']='obtuse old wall';br['wall']=int(np.argmin(co));stats['old_obtuse']+=1
   else:
    stats['acute']+=1;native=EV.classify(co)
    if native is None:
     bad=[]
     for j,x in enumerate(co):
      if 1e-7<x<1-1e-7:
       z=math.pi/math.acos(x);lab=round(z)
       if lab<3 or abs(x-math.cos(math.pi/lab))>1e-7:bad.append(dict(wall=j,cosine=float(x),approx_label=float(z)))
     br.update(rejection='non-Coxeter old contact',bad_contacts=bad);stats['non_Coxeter']+=1
    else:
     stats['Coxeter']+=1;face=K.poly(np.r_[Y,[r]],T)
     br.update(all168_labels=native,face=face)
     candidates.append(dict(profile=row,branch=sign,residual4=r.tolist(),T_pairings=(-dt).tolist(),all168_labels=native,face=face))
   rec['branches'].append(br)
  records.append(rec)
  if row['kind']=='fixed':print('FIXED',row['id'],'accepted',len(candidates),flush=True)
 paths=[f'scripts/{TAG}.py',END,'results/fresh_h20_e8_cusp_completion_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json',f'scripts/{P.TAG}.py',f'scripts/{K.TAG}.py',f'scripts/{EV.TAG}.py']
 return dict(status='COMPLETE_FIXED_PRISM_PROFILE_INTERSECTION_WITH_NUMERICAL_AMBIENT_SCOUT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  fixed_profile_intersection=rows,isolated_parametric_profile=dict(contacts=[[151,'m']],domain='all integers m>=3; only selected m evaluated numerically here'),
  actual_old_T_facet_ids=positive,projected_norms=projections.tolist(),statistics=dict(stats),singleton_labels_tested=labels,records=records,candidates=candidates,
  scope='Exact complete intersection of all three fixed spherical endpoint catalogues, plus the shared singleton151 family. Both norm-completion branches tested numerically against all168 native products; singleton labels sampled, no all-label negative claim. Actual old facet incidence screen uses positive projected square and acute-coherence projection lemma. No exact new cap or finite-volume claim without separate promotion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],'fixed',len(d['fixed_profile_intersection']),'candidates',len(d['candidates']),flush=True)
