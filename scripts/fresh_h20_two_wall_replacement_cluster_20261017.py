#!/usr/bin/env python3
"""A small two-old-wall replacement with explicit finite partner tests.

Numerical discovery: delete exactly31 and158, keep166 old roots, and test a
specific nonright-T seed, its perpendicular triangular cap, and one partner.
"""
import argparse,hashlib,json,math
from collections import Counter
import numpy as np
from scipy.optimize import linprog
import fresh_h20_h4_edge_continuation_20261017 as H
P,R,ROOT=H.P,H.R,H.ROOT
TAG='fresh_h20_two_wall_replacement_cluster_20261017'
def run():
 G,Y,CT,J,q,V,C=P.geometry();m=G['m'];T=G['T'];rep=P.P.load('results/fresh_h20_prism_replacement_20261017.json');seed=P.P.load(H.INPUT)['candidates'][H.SEED];D=P.P.load(f'results/{H.TAG}.json');tri=D['candidates'][0]
 assert tri['contacts']==[] and tri['edge']==[142,169]and tri['anchor']==152 and tri['label']==2
 Y=np.r_[Y,[rep['residual4'][1],seed['residual4'],tri['residual4']]];CT=np.r_[CT,np.zeros((1,17)),[seed['T_pairings'],tri['T_pairings']]];q=np.r_[q,2,2,2];ret=[j for j in range(171)if j not in(31,158)]
 base=R.poly(Y,ret);rows=[]
 for z in D['at_most_three_bad_gates']:
  if z['census']!=0 or z['bad_wall_ids']!=[31]:continue
  if any(max(abs(np.array(w['residual4'])-z['residual4']))<1e-6 and max(abs(np.array(w['T_pairings'])-z['T_pairings']))<1e-6 for w in rows):continue
  rows.append(z)
 stats=Counter();out=[];trials=[]
 for ix,z in enumerate(rows):
  r=np.array(z['residual4']);dt=np.array(z['T_pairings']);co=-(r@m@Y.T+CT@J@dt)/np.sqrt(2*q);typ=R.classify(co[170])
  if isinstance(typ,int)and typ>=3 and abs(co[170]-math.cos(math.pi/typ))>1e-10:typ='angle_gap'
  stats[str(typ)]+=1;trials.append(dict(index=ix,edge=z['edge'],profile=z['profile'],anchor=z['anchor'],label=z['label'],triangle_pair_cosine=float(co[170]),triangle_pair_type=typ))
  if typ in('obtuse','angle_gap'):continue
  assert all(R.classify(co[j])not in('obtuse','angle_gap')for j in ret)
  YY=np.r_[Y,[r]];face=R.poly(YY,ret+[171]);sel=[j for j in ret+[171]if np.linalg.norm(YY[j,1:])>1e-9];A=YY[sel,1:];b=YY[sel,0]
  bounds=[]
  for i in range(3):
   for sg in(-1,1):
    f=np.zeros(3);f[i]=sg;lp=linprog(f,A_ub=A,b_ub=b,bounds=[(None,None)]*3,method='highs');assert lp.success;bounds.append(float(lp.fun))
  patch=[v for v in face['vertices']if 169 in v['active']];pd=[v['radius_squared']for v in patch]
  out.append(dict(index=ix,source_event=z,mutual_triangle_label=typ,face=face,chart_linear_bounds=bounds,seed_facet_vertices=patch,seed_facet_compact=bool(pd)and max(pd)<1-1e-6,seed_facet_maximum_radius_squared=max(pd)if pd else None))
 out.sort(key=lambda z:(not z['seed_facet_compact'],z['face']['maximum_squared_radius'],z['face']['exterior']))
 paths=[f'scripts/{TAG}.py',f'results/{H.TAG}.json',H.INPUT,'results/fresh_h20_prism_replacement_20261017.json',f'scripts/{H.TAG}.py',f'scripts/{P.TAG}.py',f'scripts/{R.TAG}.py']
 return dict(status='NUMERICAL_TWO_OLD_WALL_REPLACEMENT_CLUSTER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},removed_old_ids=[31,158],retained_old_count=166,base_new_roots=dict(replacement=168,nonright_T_seed=169,perpendicular_triangle=170),base_face=base,eligible_distinct_partner_events=len(rows),mutual_pair_counts=dict(stats),all_partner_trials=trials,candidates=out,scope='Separate branch deleting exactly31 and158, preserving the other166 original roots. Three fixed new roots plus one finite-event partner are tested against every retained native gate. Charts are numerically proved bounded by six coordinate LPs before compactness of the seed facet is inferred from all its vertices. All signs, finite coincidences, and cones remain numerical; no exact ambient promotion or wholeH20 volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'partners',d['eligible_distinct_partner_events'],'accepted',len(d['candidates']),'compact_seed',sum(z['seed_facet_compact']for z in d['candidates']),flush=True)
