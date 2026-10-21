#!/usr/bin/env python3
"""Close the one remaining seed-facet gap using a new finite-contact partner.

All fixed endpoint profiles of T+152+169 are tested; the further contact to171
is constrained by an actual finite T+169+171+new corner. Numerical discovery.
"""
import argparse,hashlib,json,math
from collections import Counter
import numpy as np
import fresh_h20_h4_edge_continuation_20261017 as H
import fresh_h20_two_wall_replacement_cluster_20261017 as Q
P,R,B,O,ROOT=H.P,H.R,H.B,H.O,H.ROOT
TAG='fresh_h20_h4_hexagon_closure_20261017'
def run():
 G,Y,CT,J,q,V,C=P.geometry();m=G['m'];T=G['T'];rep=P.P.load('results/fresh_h20_prism_replacement_20261017.json');seed=P.P.load(H.INPUT)['candidates'][H.SEED];D=P.P.load(f'results/{H.TAG}.json');tri=D['candidates'][0];BAS=P.P.load(f'results/{Q.TAG}.json');z=next(z for z in BAS['candidates']if z['index']==62)['source_event']
 assert z['contacts']==[[95,4],[169,3]]and z['anchor']==161 and z['label']==3
 Y=np.r_[Y,[rep['residual4'][1],seed['residual4'],tri['residual4'],z['residual4']]];CT=np.r_[CT,np.zeros((1,17)),[seed['T_pairings'],tri['T_pairings'],z['T_pairings']]];q=np.r_[q,2,2,2,2];ret=[i for i in range(172)if i not in(31,158)]
 face=R.poly(Y,ret);V=np.array([[1.]+v['point']for v in face['vertices']if v['radius_squared']<1-1e-6]);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
 cen=D['endpoint_censuses'][2];assert cen['edge']==[152,169];stats=Counter();out=[];bad=[]
 for row in cen['profiles']:
  for sg,p,a,b,dt,offset in O.conic(cen,row,G,Y,CT,J,q):
   aa=a@m@Y.T;bb=b@m@Y.T;cc=p@m@Y.T+offset
   for lab in range(2,6):
    target=-2*(0 if lab==2 else math.cos(math.pi/lab));stats['equations']+=1
    for k in B.roots(aa[171],cc[171]-target,bb[171]):
     stats['real_branches']+=1;r=p+k*a+b/k;ps=r@m@Y.T+offset;co=-ps/np.sqrt(2*q)
     # This third finite corner is what turns two individual endpoint caps
     # into a closed replacement of the final exterior seed-facet corner.
     YY=np.r_[Y,[r]];CC=np.r_[CT,[dt]];whole=YY@m@YY.T+CC@J@CC.T;corner=T+[169,171,172]
     if min(np.linalg.eigvalsh(whole[np.ix_(corner,corner)]))<=1e-8:stats['closing_corner_not_PD']+=1;continue
     if max(r@m@V.T)>=-1e-8:stats['old_vertex_not_strict']+=1;continue
     labels=[R.classify(x)for x in co];fails=[j for j in ret if labels[j]in('obtuse','angle_gap')]
     rec=dict(profile=row['id'],contacts=row['contacts'],branch=sg,anchor=171,label=lab,k=k,residual4=r.tolist(),T_pairings=dt.tolist(),failed_gates=[dict(wall=j,kind=labels[j],cosine=float(co[j]))for j in fails])
     if fails:bad.append(rec);continue
     rec['all_retained_labels']=[labels[j]for j in ret];rec['face']=R.poly(YY,ret+[172]);patch=[v for v in rec['face']['vertices']if 169 in v['active']];rec['seed_facet_vertices']=patch;rec['seed_facet_compact']=bool(patch)and max(v['radius_squared']for v in patch)<1-1e-6;out.append(rec);stats['accepted']+=1
 paths=[f'scripts/{TAG}.py',f'results/{Q.TAG}.json',f'results/{H.TAG}.json',H.INPUT,'results/fresh_h20_prism_replacement_20261017.json',f'scripts/{H.TAG}.py',f'scripts/{Q.TAG}.py',f'scripts/{O.TAG}.py']
 return dict(status='NUMERICAL_NEW_NEW_FINITE_CONTACT_HEXAGON_TEST',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},removed_old_ids=[31,158],base_new_roots=[168,169,170,171],partner_source_index=62,base_face=face,statistics=dict(stats),candidates=out,rejected_closing_corners=bad,scope='All27 fixed finite endpoint profiles on T+152+169 and both norm branches, with new171 finite labels2..5 and the actual closing20-root Gram PD. This is numerical discovery on a previously numerical branch; singleton151 profiles and ideal endpoints are not exhausted. No exact ambient or wholeH20 finite-volume assertion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],'compact',sum(z['seed_facet_compact']for z in d['candidates']),flush=True)
