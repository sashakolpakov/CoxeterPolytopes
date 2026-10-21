#!/usr/bin/env python3
"""Jointly vary the first belt pair and a cap closing its new173 facet."""
import argparse,hashlib,json,math
from collections import Counter
import numpy as np
import fresh_h20_replacement_outer_patch_20261017 as S
O,L,ROOT=S.O,S.L,S.ROOT;TAG='fresh_h20_replacement_belt_variants_20261017';PAIR='results/fresh_h20_replacement_outer_patch_pairs_20261017.json'
def run():
 D=O.P.P.load(S.INPUT);E=O.P.P.load(f'results/{S.TAG}.json');F=O.P.P.load(PAIR);GG=D['geometry'];Y0=np.array(GG['residual_normals']);CT0=np.array(GG['T_products']);J=np.array(GG['T_inverse']);q0=np.array(GG['native_norms']);m=np.array(GG['residual_form']);T=D['support'];G=dict(T=T,m=m);ret=D['stages'][-1]['retained_ids']+[173,174];first=E['candidates'][0];choices=[(i,z)for i,z in enumerate(F['candidates'])if z['first_candidate']==0 and not any(j in T for j,lab in z['contacts'])];cache={};out=[];trials=[];near=[];stats=Counter()
 for idx,second in choices:
  Y=np.r_[Y0,[first['residual4'],second['residual4']]];CT=np.r_[CT0,[first['T_pairings'],second['T_pairings']]];q=np.r_[q0,2,2];whole=Y@m@Y.T+CT@J@CT.T;bond=O.R.classify(float(-whole[173,174]/2));assert isinstance(bond,int);W=T+[173,174]
  if bond not in cache:
   cache[bond]=O.census(W,whole[np.ix_(W,W)]);cache[bond]['edge']=[173,174]
  cen=cache[bond];rows=list(cen['profiles'])
  for j in cen['parametric_singletons']:rows.extend(dict(id=f'single{j}_{lab}',contacts=[[j,lab]])for lab in range(3,13))
  V=np.array([[1.]+v['point']for v in second['boundary']['poly_vertices']if v['kind']=='finite']);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
  for row in rows:
   for sg,p,a,b,dt,offset in O.conic(cen,row,G,Y,CT,J,q):
    aa=float(a@m@Y[171]);bb=float(b@m@Y[171]);cc=float(p@m@Y[171]+offset[171])
    for lab in range(2,7):
     stats['equations']+=1;target=-2*(0 if lab==2 else math.cos(math.pi/lab))
     for k in O.B.roots(aa,cc-target,bb):
      stats['real_branches']+=1;r=p+k*a+b/k;YY=np.r_[Y,[r]];CC=np.r_[CT,[dt]];full=YY@m@YY.T+CC@J@CC.T;corner=T+[173,171,175]
      if min(np.linalg.eigvalsh(full[np.ix_(corner,corner)]))<=1e-8:stats['closing_corner_not_PD']+=1;continue
      if max(r@m@V.T)>=-1e-8:stats['old_corner_violation']+=1;continue
      co=-(r@m@Y.T+offset)/np.sqrt(2*q);labs=[O.R.classify(float(x))for x in co];bad=[j for j in ret if labs[j]in('obtuse','angle_gap')];trials.append([idx,row['id'],sg,lab,float(k),bad]);rec=dict(second_candidate=idx,contacts=row['contacts'],profile=row['id'],branch=sg,label_to171=lab,k=float(k),residual4=r.tolist(),T_pairings=dt.tolist(),failed_gates=[dict(wall=j,kind=labs[j],cosine=float(co[j]))for j in bad])
      if bad:
       if len(bad)<=3:near.append(rec)
       continue
      rec['all_retained_labels']=[labs[j]for j in ret];bbound=L.boundary(YY,m,ret+[175],T);rec['boundary']=bbound;rec['closed173']=next(z for z in bbound['facets']if z['wall']==173)['compact'];out.append(rec);stats['accepted']+=1
  print('VARIANT',idx,'bond',bond,'accepted',len(out),dict(stats),flush=True)
 paths=[f'scripts/{TAG}.py',S.INPUT,f'results/{S.TAG}.json',PAIR,f'scripts/{O.TAG}.py',f'scripts/{L.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_THREE_CAP_BELT_VARIANT_TEST',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},base_first=0,second_candidate_ids=[i for i,z in choices],endpoint_censuses={str(k):v for k,v in cache.items()},statistics=dict(stats),candidates=out,at_most_three_bad_gates=near,trials=trials,scope='All right-T second-cap variants with first candidate0 in the frozen272-pair bank, followed by fixed173/174 profiles, singleton3..12, both conic branches, and finite171 angle2..6. All retained173 gates, previous finite corners and actual closing20-root corner checked numerically. Whole target boundary recomputed for survivors. Bounded numerical domain; no exact or all-family negative conclusion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
