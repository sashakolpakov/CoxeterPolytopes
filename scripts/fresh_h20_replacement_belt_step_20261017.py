#!/usr/bin/env python3
"""Next finite-adjacency belt step after the two-cap95/157 closure."""
import argparse,hashlib,json,math
from collections import Counter
import numpy as np
import fresh_h20_replacement_outer_patch_20261017 as S
O,L,ROOT=S.O,S.L,S.ROOT;TAG='fresh_h20_replacement_belt_step_20261017';PAIR='results/fresh_h20_replacement_outer_patch_pairs_20261017.json'
CASES=[dict(edge=[145,152],anchor=174,shared=145),dict(edge=[171,172],anchor=173,shared=171),dict(edge=[173,174],anchor=171,shared=173)]
def run():
 D=O.P.P.load(S.INPUT);E=O.P.P.load(f'results/{S.TAG}.json');F=O.P.P.load(PAIR);GG=D['geometry'];Y=np.array(GG['residual_normals']);CT=np.array(GG['T_products']);J=np.array(GG['T_inverse']);q=np.array(GG['native_norms']);m=np.array(GG['residual_form']);T=D['support'];G=dict(T=T,m=m);ret=D['stages'][-1]['retained_ids']+[173,174];first=E['candidates'][0];second=F['candidates'][77];assert second['first_candidate']==0 and second['contacts']==[[157,3]]and second['label_to_first']==3
 for z in(first,second):Y=np.r_[Y,[z['residual4']]];CT=np.r_[CT,[z['T_pairings']]];q=np.r_[q,2]
 base=second['boundary'];V=np.array([[1.]+v['point']for v in base['poly_vertices']if v['kind']=='finite']);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None];whole=Y@m@Y.T+CT@J@CT.T
 censuses=[];out=[];near=[];trials=[];stats=Counter()
 for ci,case in enumerate(CASES):
  W=T+case['edge'];cen=O.census(W,whole[np.ix_(W,W)]);cen.update(case);censuses.append(cen);rows=list(cen['profiles'])
  for j in cen['parametric_singletons']:rows.extend(dict(id=f'single{j}_{lab}',contacts=[[j,lab]])for lab in range(3,13))
  print('CENSUS',case,len(cen['profiles']),len(rows),flush=True)
  for row in rows:
   for sg,p,a,b,dt,offset in O.conic(cen,row,G,Y,CT,J,q):
    anchor=case['anchor'];aa=float(a@m@Y[anchor]);bb=float(b@m@Y[anchor]);cc=float(p@m@Y[anchor]+offset[anchor])
    for lab in range(2,7):
     stats['equations']+=1;target=-2*(0 if lab==2 else math.cos(math.pi/lab))
     for k in O.B.roots(aa,cc-target,bb):
      stats['real_branches']+=1;r=p+k*a+b/k;YY=np.r_[Y,[r]];CC=np.r_[CT,[dt]];full=YY@m@YY.T+CC@J@CC.T;corner=T+[case['shared'],anchor,175]
      if min(np.linalg.eigvalsh(full[np.ix_(corner,corner)]))<=1e-8:stats['joint_corner_not_PD']+=1;continue
      if max(r@m@V.T)>=-1e-8:stats['old_corner_violation']+=1;continue
      co=-(r@m@Y.T+offset)/np.sqrt(2*q);labs=[O.R.classify(float(x))for x in co];bad=[j for j in ret if labs[j]in('obtuse','angle_gap')];trials.append([ci,row['id'],sg,lab,float(k),bad]);rec=dict(case=ci,edge=case['edge'],contacts=row['contacts'],profile=row['id'],branch=sg,anchor=anchor,label=lab,k=float(k),residual4=r.tolist(),T_pairings=dt.tolist(),failed_gates=[dict(wall=j,kind=labs[j],cosine=float(co[j]))for j in bad])
      if bad:
       if len(bad)<=3:near.append(rec)
       continue
      if any(np.linalg.norm(r-z['residual4'])<1e-6 and np.linalg.norm(dt-z['T_pairings'])<1e-6 for z in out):stats['duplicate']+=1;continue
      rec['all_retained_labels']=[labs[j]for j in ret];bbound=L.boundary(YY,m,ret+[175],T);rec['boundary_summary']={k:bbound[k]for k in('poly_vertex_counts','actual_edge_count','maximum_squared_radius','open_facet_incidence_components','facets')};rec['open_pairs']=[z['pair']for z in bbound['open_edges']];rec['finite_vertices']=[z for z in bbound['poly_vertices']if z['kind']=='finite'];rec['exterior_vertices']=[z for z in bbound['poly_vertices']if z['kind']=='exterior'];rec['closed_target_facets']=[z['wall']for z in bbound['facets']if z['compact']];out.append(rec);stats['accepted_unique']+=1
  print('CASE',ci,dict(stats),flush=True)
 paths=[f'scripts/{TAG}.py',S.INPUT,f'results/{S.TAG}.json',PAIR,f'scripts/{O.TAG}.py',f'scripts/{L.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_COORDINATED_OUTER_BELT_STEP',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},base_first=0,base_second=77,retained_ids=ret,cases=CASES,censuses=censuses,statistics=dict(stats),candidates=out,at_most_three_bad_gates=near,trials=trials,scope='Three prescribed finite-adjacency continuations of a numerical two-cap branch, full fixed endpoint profiles plus singleton3..12, both branches and anchor2..6. All173 retained native gates, finite target vertices, and actual closing corner checked numerically. Full positive-side/ball-clipped boundary recomputed for each survivor. No exact or all-parameter conclusion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
