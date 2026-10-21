#!/usr/bin/env python3
"""Exact-promotable rational free-contact recipe, numerically checked only."""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.linalg import null_space
import fresh_h20_spherical_attachment_profiles_20261016_scout as S
from fresh_h20_e8_triangle_cluster_20261016 import census
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_triangle_cluster_20261016_recipe'
INPUT='results/fresh_h20_e8_triangle_cluster_20261016.json'
def load(p):return json.loads((ROOT/p).read_text())
def finite_deletions():
 # The non-E8 component is a13-cycle with one4-edge. Each deletion is
 # a path, so exact rational leaf elimination uses only squared edge weights.
 vertices=list(range(1,11))+[159,160,161]
 edges=[(i,i+1,1)for i in range(1,10)]+[(1,159,1),(10,160,1),(159,161,2),(160,161,1)]
 out=[]
 for omit in vertices:
  diag={i:Fraction(2)for i in vertices if i!=omit};ee={(i,j):Fraction(q)for i,j,q in edges if omit not in(i,j)};pivots=[];pd=True
  while diag:
   leaf=min(i for i in diag if sum(i in e for e in ee)<=1)
   v=diag.pop(leaf);pivots.append([leaf,str(v)])
   if v<=0:pd=False;break
   pair=next((e for e in ee if leaf in e),None)
   if pair is not None:
    other=next(j for j in pair if j!=leaf);diag[other]-=ee.pop(pair)/v
  out.append(dict(omitted=omit,positive_definite=pd,leaf_pivots=pivots))
 assert [r['omitted']for r in out if r['positive_definite']]==[1,159,160,161]
 return out
def run():
 D,G=S.geometry();B=load(INPUT);rec=next(r for r in B['accepted']if(r['E8_profile'],r['direction'],r['k_factor'])==(0,5,1.25))
 desired=np.array(rec['residual4']);profiles=[0,38,42,56];roots=[];steps=[];m=G['m'];notT=[j for j in range(158)if j not in G['T']]
 # Full Gram of the positive T projections of the four new roots.
 P=np.array([[59/30,0,0,0],[0,10/11,1/11,0],[0,1/11,10/11,0],[0,0,0,0.]])
 target=np.array([[2,0,0,0],[0,2,0,-np.sqrt(2)],[0,0,2,-1],[0,-np.sqrt(2),-1,2.]])
 free=[[95,142,145],[95,142],[95],[]]
 for i,prof in enumerate(profiles):
  row=D['census']['profiles'][prof];offset=S.offsets(G,row);U=[];f=[];constraints=[]
  for j,r in enumerate(roots):
   U.append(r);f.append(target[i,j]-P[i,j]);constraints.append(dict(kind='earlier_new_root',index=158+j,full_pair=('0'if target[i,j]==0 else'-sqrt(2)'if target[i,j]==-np.sqrt(2)else'-1')))
  for j in free[i]:
   pair=desired[i]@m@G['Y'][j]+offset[j];rat=Fraction(round(pair*10000),10000)
   U.append(G['Y'][j]);f.append(float(rat)-offset[j]);constraints.append(dict(kind='old_root',index=j,full_pair=str(rat)))
  U=np.array(U);f=np.array(f);J=U@m@U.T;assert abs(np.linalg.det(J))>1e-10
  p=np.linalg.solve(J,f)@U;k=null_space(U@m)[:,0];kn=k@m@k;h=row['numeric_residual_norm'];rad=(h-p@m@p)/kn;assert rad>0
  branches=np.array([p-np.sqrt(rad)*k,p+np.sqrt(rad)*k]);which=int(np.argmin(np.linalg.norm(branches-desired[i],axis=1)));r=branches[which]
  mask,co=S.gates(G,r[None,:],offset);assert mask[0]
  roots.append(r);steps.append(dict(root_index=158+i,support_contacts=row['contacts'],other_T_contacts='right',free_and_previous_constraints=constraints,residual_target_norm=h,frame_gram=J.tolist(),quadratic_particular=p.tolist(),quadratic_kernel=k.tolist(),kernel_norm=float(kn),quadratic_radical=float(rad),branch=which,residual4=r.tolist(),distance_from_discovery=float(np.linalg.norm(r-desired[i])),minimum_old_nonsupport_cosine=float(min(co[0,notT]))))
 R=np.array(roots);got=R@m@R.T+P;assert np.max(abs(got-target))<1e-7
 # Numeric active sets/corners in the target H3; full Gram types are explicit.
 result=census(G,R)
 src=[f'scripts/{TAG}.py',INPUT,S.INPUT,f'scripts/{S.TAG}.py','scripts/fresh_h20_e8_triangle_cluster_20261016.py']
 return dict(status='NUMERICALLY_VALIDATED_EXACT_PROMOTABLE_E8_TRIANGLE_RECIPE',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},base_walls=158,discovery=dict(profile=0,direction=5,k_factor=1.25),steps=steps,full_new_gram=got.tolist(),target_gram=[['2','0','0','0'],['0','2','0','-sqrt(2)'],['0','0','2','-1'],['0','-sqrt(2)','-1','2']],positive_T_projection_gram=P.tolist(),new_finite_corner_supports=[G['T']+[158,159,160],G['T']+[158,159,161],G['T']+[158,160,161],[j for j in G['T']if j!=1]+[158,159,160,161]],new_finite_corner_types=['E8+A12','E8+B12','E8+A12','E8+B12'],exact_cycle_deletions=finite_deletions(),face=result,scope='All coefficients defining free old pairings are exact rationals; all finite contacts have exact Coxeter targets. The branch, inequalities, and geometry are numerically validated here only. Exact rational leaf pivots classify all13 cycle deletions. Exact radical reconstruction and interval signs remain required before promotion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,load(str(out.relative_to(ROOT))))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['face'],flush=True)
