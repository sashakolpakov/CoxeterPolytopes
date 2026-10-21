#!/usr/bin/env python3
"""Solve coupled endpoint caps aiming to close the entire157 target facet."""
import argparse,hashlib,json,math
from collections import Counter
import numpy as np
import fresh_h20_replacement_outer_patch_20261017 as S
O,L,ROOT=S.O,S.L,S.ROOT;TAG='fresh_h20_replacement_outer_patch_pairs_20261017'
def run():
 D=O.P.P.load(S.INPUT);E=O.P.P.load(f'results/{S.TAG}.json');GG=D['geometry'];Y0=np.array(GG['residual_normals']);CT0=np.array(GG['T_products']);J=np.array(GG['T_inverse']);q0=np.array(GG['native_norms']);m=np.array(GG['residual_form']);T=D['support'];G=dict(T=T,m=m);ret0=D['stages'][-1]['retained_ids'];cen=next(z for z in E['censuses']if z['edge']==[145,157]);rows=list(cen['profiles'])
 for j in cen['parametric_singletons']:rows.extend(dict(id=f'single{j}_{lab}',contacts=[[j,lab]])for lab in range(3,13))
 stats=Counter();out=[];near=[];trials=[]
 for first in E['candidates'][:4]:
  Y=np.r_[Y0,[first['residual4']]];CT=np.r_[CT0,[first['T_pairings']]];q=np.r_[q0,2];ret=ret0+[173];base=L.boundary(Y,m,ret,T);V=np.array([[1.]+v['point']for v in base['poly_vertices']if v['kind']=='finite']);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
  for row in rows:
   for sg,p,a,b,dt,offset in O.conic(cen,row,G,Y,CT,J,q):
    aa=float(a@m@Y[173]);bb=float(b@m@Y[173]);cc=float(p@m@Y[173]+offset[173])
    for lab in range(2,7):
     stats['equations']+=1;target=-2*(0 if lab==2 else math.cos(math.pi/lab))
     for k in O.B.roots(aa,cc-target,bb):
      stats['real_branches']+=1;r=p+k*a+b/k;YY=np.r_[Y,[r]];CC=np.r_[CT,[dt]];whole=YY@m@YY.T+CC@J@CC.T;corner=T+[157,173,174]
      if min(np.linalg.eigvalsh(whole[np.ix_(corner,corner)]))<=1e-8:stats['joint_corner_not_PD']+=1;continue
      if max(r@m@V.T)>=-1e-8:stats['old_corner_violation']+=1;continue
      co=-(r@m@Y.T+offset)/np.sqrt(2*q);labels=[O.R.classify(float(x))for x in co];bad=[j for j in ret if labels[j]in('obtuse','angle_gap')];trials.append([first['id'],row['id'],sg,lab,float(k),bad]);rec=dict(first_candidate=first['id'],contacts=row['contacts'],profile=row['id'],branch=sg,label_to_first=lab,k=float(k),residual4=r.tolist(),T_pairings=dt.tolist(),failed_gates=[dict(wall=j,kind=labels[j],cosine=float(co[j]))for j in bad])
      if bad:
       if len(bad)<=3:near.append(rec)
       continue
      rec['all_retained_labels']=[labels[j]for j in ret];rec['boundary']=L.boundary(YY,m,ret+[174],T);rec['closed157']=next(z for z in rec['boundary']['facets']if z['wall']==157)['compact'];out.append(rec);stats['accepted']+=1
  print('FIRST',first['id'],'accepted',len(out),'near',len(near),dict(stats),flush=True)
 paths=[f'scripts/{TAG}.py',S.INPUT,f'results/{S.TAG}.json',f'scripts/{S.TAG}.py',f'scripts/{L.TAG}.py',f'scripts/{O.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_TWO_CAP157_CLOSURE_TEST',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},first_candidate_ids=[0,1,2,3],second_edge=[145,157],statistics=dict(stats),candidates=out,at_most_three_bad_gates=near,all_joint_corner_trials=trials,scope='Second caps from all fixed finite145/157 endpoint profiles plus singleton labels3..12, both conic branches, exact equation for first-cap angle2..6 evaluated numerically. Actual joint20-root corner must bePD, all retained171+first native gates and finite target vertices tested. Complete replacement T-H3 boundary recomputed for each survivor. Bounded numerical discovery only.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
