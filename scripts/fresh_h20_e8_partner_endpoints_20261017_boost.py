#!/usr/bin/env python3
"""Solve E8 boost and partner rotation against old-compatible endpoint caps.

For a fixed old-compatible norm2 cap r, its contacts to two partner roots
give a quadratic in r.h; the E8 boost then satisfies a second quadratic.
No boost or rotation samples are used. The input cap list is bounded.
"""
import argparse,hashlib,json,math
import numpy as np
from scipy.linalg import null_space
import fresh_h20_e8_partner_endpoints_20261017 as P
from fresh_h20_e8_partner_endpoints_20261017_coordinated import poly
ROOT=P.ROOT;TAG=P.TAG+'_boost';INPUT='results/'+P.TAG+'.json'
BANK='results/fresh_h20_spherical_attachment_profiles_20261016_scout.json'
def quadratic(a,b,c):
 if abs(a)<1e-12:return [-c/b]if abs(b)>1e-12 else[]
 disc=b*b-4*a*c
 if disc<-1e-9:return []
 z=math.sqrt(max(0.,disc));return [(-b-z)/(2*a),(-b+z)/(2*a)]
def run():
 D=P.load(INPUT);B=P.load(BANK);G,Y,CT,J,q,H,live,verts=P.geometry();m=G['m']
 direction=Y[158,1:]/np.linalg.norm(Y[158,1:]);E=np.c_[np.zeros(2),null_space(direction[None,:]).T]
 M=Y[159:162]@m@Y[159:162].T;aa=-np.linalg.solve(M,np.ones(3));t=1/math.sqrt(-(aa@M@aa))
 spatialGram=M+t*t;ev,V=np.linalg.eigh(spatialGram);spatial=V[:,1:]*np.sqrt(ev[1:])[None,:]
 off0=CT[158]@J@CT[:158].T;off=CT[159:162]@J@CT[:158].T
 notT=[j for j in range(158)if j not in G['T']];out=[];allrows=[];stats=dict(input_caps=0,pair_profiles=0,quadratic_solutions=0,accepted=0)
 a0=math.sqrt(1/30)/2*(np.r_[1.,direction]@m@Y[:158].T)
 b0=math.sqrt(1/30)/2*(np.r_[-1.,direction]@m@Y[:158].T)
 acute=[(0.,math.inf)];dotted=[(0.,math.inf)]
 for j in notT:
  acute=P.intersections(acute,P.quadratic_intervals(a0[j],off0[j],b0[j]))
  dotted=P.intersections(dotted,P.quadratic_intervals(a0[j],off0[j]+math.sqrt(2*q[j]),b0[j]))
 for bi,cap in enumerate(B['candidates']):
  stats['input_caps']+=1;r=np.array(cap['residual4']);v=r@m@E.T
  cT=np.zeros(17)
  for j,lab in cap['contacts']:cT[G['T'].index(j)]=-2*math.cos(math.pi/lab)
  for cen in D['endpoint_censuses']:
   pair=[j-159 for j in cen['edge']];si=np.linalg.inv(spatial[pair])
   for row in cen['profiles']:
    if sorted([z for z in row['contacts']if z[0]in G['T']])!=sorted(cap['contacts']):continue
    target=-(CT[cen['edge']]@J@cT)
    for j,lab in row['contacts']:
     if j not in G['T']:target[cen['edge'].index(j)]-=2*math.cos(math.pi/lab)
    z0=si@target;z1=-t*(si@np.ones(2));stats['pair_profiles']+=1
    for rh in quadratic(z1@z1,2*z0@z1,z0@z0-v@v):
     ra=(-r[0]+r[1:]@direction)/2;rb=(-r[0]-r[1:]@direction)/2
     for k in quadratic(ra,-rh,rb):
      if k<=0:continue
      z=z0+z1*rh
      if np.linalg.norm(v)<1e-10:continue
      theta=math.atan2(v[0]*z[1]-v[1]*z[0],v@z);rot=np.array([[math.cos(theta),-math.sin(theta)],[math.sin(theta),math.cos(theta)]])
      hh=np.r_[(k+1/k)/2,(k-1/k)/2*direction];r0=math.sqrt(1/30)/2*np.r_[k-1/k,(k+1/k)*direction]
      WW=t*hh+spatial@rot@E
      assert max(abs(WW[pair]@m@r-target))<1e-6
      co0=-(r0@m@Y[:158].T+off0)/np.sqrt(2*q[:158]);cop=-(WW@m@Y[:158].T+off)/np.sqrt(2*q[:158])[None,:]
      remaining=next(j for j in range(3)if j not in pair);c0=-(r0@m@r+CT[158]@J@cT)/2;cx=-(WW[remaining]@m@r+CT[159+remaining]@J@cT)/2
      mm=dict(E8_old=float(min(co0[notT])),partners_old=float(min(cop[:,notT].ravel())),cap_E8=float(c0),unpaired_partner=float(cx));margin=min(mm.values())-1
      stats['quadratic_solutions']+=1
      rec=dict(bank_index=bi,T_contacts=cap['contacts'],old_contact_sides=cap['extra_old_sides'],old_contact_labels=cap['extra_labels'],new_edge=cen['edge'],endpoint_profile=row['id'],endpoint_contacts=row['contacts'],boost=k,rotation=theta,triangle_frame_time_product=float(rh),E8_residual4=r0.tolist(),partner_residual4=WW.tolist(),new_residual4=r.tolist(),gate_minima=mm,minimum_dotted_margin=float(margin))
      allrows.append(rec)
      if margin>=-1e-6:
       rec['face']=poly(np.r_[Y[:158],[r0],WW,[r]],G['T'])
       if rec['face']is not None:out.append(rec);stats['accepted']+=1
 print('BOOST',stats,flush=True);allrows.sort(key=lambda r:-r['minimum_dotted_margin'])
 src=[f'scripts/{TAG}.py',INPUT,BANK,P.RECIPE,f'scripts/{P.TAG}.py','scripts/'+P.TAG+'_coordinated.py']
 return dict(status='BOUNDED_OLD_CAP_E8_BOOST_ROTATION_QUADRATIC_SOLVER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},statistics=stats,direction=direction.tolist(),E8_necessary_acute_boost_intervals=[[lo,None if math.isinf(hi)else hi]for lo,hi in acute],E8_sufficient_dotted_boost_intervals=[[lo,None if math.isinf(hi)else hi]for lo,hi in dotted],largest_solved_boost=max(r['boost']for r in allrows),largest_E8_old_minimum_cosine=max(r['gate_minima']['E8_old']for r in allrows),candidates=out,all_quadratic_solutions=allrows,scope='Fixed ideal direction of the E8 root, complete two-quadratic solutions for551 bounded old-compatible caps and every matching fixed endpoint profile among159. Unprescribed contacts are required dotted. The whole E8 triangle cluster is rebuilt; this is an alternative162/163 branch, not an extension retaining the frozen162 recipe. Numerical discovery only.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,P.load(str(out.relative_to(ROOT))))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
