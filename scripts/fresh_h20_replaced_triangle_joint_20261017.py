#!/usr/bin/env python3
"""Bounded numerical joint-contact and whole-face test of a replacement bank."""
import argparse,hashlib,json,math
from collections import Counter
import numpy as np
import fresh_h20_opposite_prism_cap_20261017 as P
import fresh_h20_prism_replacement_20261017 as R
ROOT=P.ROOT;TAG='fresh_h20_replaced_triangle_joint_20261017';INPUT='results/fresh_h20_replaced_triangle_events_20261017.json'
def run():
 D=P.P.load(INPUT);G,Y,CT,J,q,V,C=P.geometry();m=G['m'];rr=P.P.load('results/fresh_h20_prism_replacement_20261017.json');r0=np.array(rr['residual4'][1]);Y=np.r_[Y,[r0]];CT=np.r_[CT,np.zeros((1,17))];q=np.r_[q,2];ret=[i for i in range(169)if i!=158]
 rows=D['candidates'];A=np.array([z['residual4']for z in rows]);dt=np.array([z['T_pairings']for z in rows]);g=A@m@A.T+dt@J@dt.T
 h=2-np.einsum('ij,jk,ik->i',dt,J,dt);compat=[];finite=[];near=[];stats=Counter();adj=[[]for _ in rows];fadj=[[]for _ in rows]
 for i in range(len(rows)):
  for j in range(i):
   co=-g[i,j]/2;typ=R.classify(co)
   if isinstance(typ,int) and typ>=3:
    err=abs(co-math.cos(math.pi/typ))
    if err>1e-10:near.append(dict(candidates=[j,i],nearest_label=typ,cosine_error=err));typ='angle_gap'
   stats[str(typ)]+=1
   if typ not in('obtuse','angle_gap'):
    compat.append([j,i,typ]);adj[i].append(j);adj[j].append(i)
    z=A[i]@m@A[j]
    if h[i]*h[j]-z*z>1e-7:finite.append([j,i,typ]);fadj[i].append(j);fadj[j].append(i)
 base=rr['individual_replacements'][1]['face'];ext=[v for v in base['vertices']if v['radius_squared']>1+1e-6];XX=np.array([[1.]+v['point']for v in ext]);cuts=A@m@XX.T>1e-7
 missing=[i for i in range(len(ext))if not any(cuts[:,i])];newY=np.r_[Y,A]
 relaxation=R.poly(newY,ret+list(range(169,169+len(A))))
 # Small coordinated pairs have exact prescribed native pair labels only
 # after independent algebraic promotion; this table is numerical.
 pair_scores=[]
 for i,j,lab in finite:
  if not any(cuts[i]|cuts[j]):continue
  z=R.poly(newY,ret+[169+i,169+j]);pair_scores.append(dict(candidates=[i,j],mutual_label=lab,base_exterior_cuts=int(sum(cuts[i]|cuts[j])),face=z))
 pair_scores.sort(key=lambda z:(z['face']['maximum_squared_radius'],z['face']['exterior'],-z['base_exterior_cuts']))
 paths=[f'scripts/{TAG}.py',INPUT,'results/fresh_h20_prism_replacement_20261017.json',f'scripts/{P.TAG}.py',f'scripts/{R.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_REPLACEMENT_BANK_JOINT_DIAGNOSTIC',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},bank_count=len(rows),pair_type_counts=dict(stats),compatible_pairs=compat,intrinsic_finite_pairs=finite,rejected_loose_angle_coincidences=near,compatible_degrees=[len(z)for z in adj],finite_degrees=[len(z)for z in fadj],base_exterior_vertices=ext,base_exterior_cut_counts=cuts.sum(axis=0).tolist(),uncut_base_exterior_indices=missing,all_bank_relaxation=relaxation,finite_pair_scores=pair_scores,scope='All pairs in the finite542-root event bank, numerical native Coxeter classification and intrinsic finite-intersection tests. Nonright finite coincidences use the stricter1e-10 cosine tolerance; loose1e-7 matches are separately recorded. The all-bank halfspace intersection ignores mutual incompatibility and is only an optimistic geometric relaxation. No continuous family or all-label completion claim; exactification is required for any chosen pair.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['bank_count'],'finite_pairs',len(d['intrinsic_finite_pairs']),'uncut',d['uncut_base_exterior_indices'],'relaxation',d['all_bank_relaxation']['finite'],d['all_bank_relaxation']['exterior'],flush=True)
