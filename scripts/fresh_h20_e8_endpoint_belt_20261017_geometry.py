#!/usr/bin/env python3
"""Numerical all-halfline polynomial screen of endpoint belt patterns."""
import argparse,hashlib,json,itertools,math
import numpy as np
from scipy.optimize import least_squares
import fresh_h20_e8_endpoint_belt_20261017 as B
P=B.P;ROOT=B.ROOT;TAG=B.TAG+'_geometry';INPUT='results/'+B.TAG+'_profiles.json'
def run():
 D=P.load(INPUT);E=P.load(B.INPUT);cc=E['endpoint_censuses'];G,Y,CT,J,q,H,live,verts=P.geometry();m=G['m'];families=[]
 for c in cc:
  ff={}
  for p in c['profiles']:
   bs=[b for b in c['families'][p['id']]['branches']if b['all_unprescribed_dotted_intervals']];assert len(bs)==1
   br=bs[0];assert len(br['all_unprescribed_dotted_intervals'])==1
   lo,hi=br['all_unprescribed_dotted_intervals'][0];assert hi is None
   a=np.array(br['positive_k_coefficient'])*lo;b=np.array(br['inverse_k_coefficient'])/lo;p0=np.array(br['particular']);dt=np.array(br['T_pairings'])
   ff[p['id']]=dict(a=a,b=b,p=p0,dt=dt,lo=lo,num=np.array([a+b+p0,2*a+p0,a]))
  families.append(ff)
 pairs={};survive=[];reject={};near=[]
 for pat in D['finite_patterns']:
  keys=[]
  for (i,j),lab in zip([(0,1),(0,2),(1,2)],pat['mutual_labels']):
   ia,ib=pat['profiles'][i],pat['profiles'][j];key=(i,j,ia,ib,lab);keys.append(key)
   if key in pairs:continue
   f,g=families[i][ia],families[j][ib];offset=f['dt']@J@g['dt'];target=0 if lab==2 else -2*math.cos(math.pi/lab)
   co=f['num']@m@g['num'].T;co[:2,:2]+=offset-target
   excluded=bool(co.max()<-1e-8)
   pairs[key]=dict(family_pair=[i,j],profiles=[ia,ib],label=lab,shifted_numerator_coefficients=co.tolist(),all_coefficients_strict_negative=excluded,max_coefficient=float(co.max()))
  bad=next((k for k in keys if pairs[k]['all_coefficients_strict_negative']),None)
  if bad:reject[str(bad)]=reject.get(str(bad),0)+1
  else:survive.append(pat)
 print('GEOMETRY',len(pairs),'pair cases',sum(r['all_coefficients_strict_negative']for r in pairs.values()),'excluded','surviving patterns',len(survive),flush=True)
 solutions=[];best=None;attempts=0
 for pi,pat in enumerate(survive):
  ff=[families[i][v]for i,v in enumerate(pat['profiles'])];off=np.array([[x['dt']@J@y['dt']for y in ff]for x in ff]);target=np.array([-2*np.cos(np.pi/l)if l>2 else 0 for l in pat['mutual_labels']])
  def rs(x):return np.array([f['p']+f['a']*np.exp(t)+f['b']/np.exp(t)for f,t in zip(ff,x)])
  def fun(x):
   rr=rs(x);gg=rr@m@rr.T+off;return np.array([gg[0,1],gg[0,2],gg[1,2]])-target
  for start in[(0,0,0),(.2,.2,.2),(1,1,1),(2,2,2),(.1,1,2),(2,.1,1),(1,2,.1)]:
   fit=least_squares(fun,np.array(start,float),bounds=(np.zeros(3),np.full(3,5.)),max_nfev=120,ftol=1e-10,xtol=1e-10,gtol=1e-10);attempts+=1;err=float(max(abs(fun(fit.x))))
   if best is None or err<best['error']:best=dict(pattern=pi,error=err,log_parameters=fit.x.tolist())
   if err<1e-6:
    rr=rs(fit.x);record=dict(pattern=pi,profiles=pat['profiles'],mutual_labels=pat['mutual_labels'],log_parameters=fit.x.tolist(),parameters=[float(f['lo']*np.exp(t))for f,t in zip(ff,fit.x)],residual4=rr.tolist(),T_pairings=[f['dt'].tolist()for f in ff],error=err)
    if not any(s['pattern']==pi and max(abs(np.array(s['log_parameters'])-fit.x))<1e-5 for s in solutions):solutions.append(record);print('SOLUTION',record,flush=True)
    break
  if pi%20==0:print('SOLVED',pi,'of',len(survive),'best',best,'solutions',len(solutions),flush=True)
 paths=[f'scripts/{TAG}.py',INPUT,B.INPUT,f'scripts/{B.TAG}.py',f'scripts/{P.TAG}.py']
 return dict(status='NUMERICAL_ENDPOINT_BELT_HALF_LINE_SCREEN_AND_SOLVE',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},input_patterns=len(D['finite_patterns']),pair_cases=list(pairs.values()),surviving_patterns=survive,rejection_counts=reject,solver=dict(log_bounds=[0,5],starts=7,max_nfev=120,attempts=attempts,best=best),solutions=solutions,scope='All coefficient signs and inherited conic endpoints are numerical. Strict negativity of every shifted polynomial coefficient would mathematically exclude the corresponding mutual angle on the whole parameter quadrant, but these saved doubles are not rigorous interval certificates. Numerical solves on survivors use explicit bounded log parameters and finitely many starts. No complete geometric exclusion or finite-volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
