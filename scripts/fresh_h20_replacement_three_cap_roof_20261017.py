#!/usr/bin/env python3
"""Bounded three-cap roof discovery closing the new173 target facet.

Numerical only. The complete replacement boundary remains open.
"""
import argparse,hashlib,json,itertools,math,collections
import numpy as np
from scipy.linalg import null_space
import fresh_h20_replacement_outer_patch_20261017 as S
O,L,ROOT=S.O,S.L,S.ROOT
TAG='fresh_h20_replacement_three_cap_roof_20261017'
def run():
 D=O.P.P.load(S.INPUT);E=O.P.P.load(f'results/{S.TAG}.json');F=O.P.P.load('results/fresh_h20_replacement_outer_patch_pairs_20261017.json');g=D['geometry'];Y0=np.array(g['residual_normals']);CT0=np.array(g['T_products']);J=np.array(g['T_inverse']);q0=np.array(g['native_norms']);m=np.array(g['residual_form']);T=D['support'];G=dict(T=T,m=m);ret0=D['stages'][-1]['retained_ids']+[173,174]
 for z in(E['candidates'][0],F['candidates'][77]):Y0=np.r_[Y0,[z['residual4']]];CT0=np.r_[CT0,[z['T_pairings']]];q0=np.r_[q0,2]
 base=F['candidates'][77]['boundary'];V=np.array([[1.]+z['point']for z in base['poly_vertices']if z['kind']=='finite']);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
 def rows(edge):return [dict(id=k,contacts=[[edge[i],x]for i,x in enumerate(z)if x>2])for k,z in enumerate([(2,2),(3,2),(4,2),(5,2),(2,3),(2,4),(2,5)])]
 def seeds(edge,anchors):
  out=[]
  for row in rows(edge):
   for sg,p,a,b,dt,offset in O.conic(dict(edge=edge),row,G,Y0,CT0,J,q0):
    aa=a@m@Y0.T;bb=b@m@Y0.T;cc=p@m@Y0.T+offset
    for anchor in anchors:
     for lab in range(2,7):
      for k in O.B.roots(aa[anchor],cc[anchor]+2*(0 if lab==2 else math.cos(math.pi/lab)),bb[anchor]):
       r=p+k*a+b/k
       if max(r@m@V.T)>=-1e-8:continue
       co=-(r@m@Y0.T+offset)/np.sqrt(2*q0);labs=[O.R.classify(float(x))for x in co]
       if any(labs[j]in('obtuse','angle_gap')for j in ret0):continue
       if any(np.linalg.norm(r-z['residual4'])<1e-6 for z in out):continue
       out.append(dict(edge=edge,contacts=row['contacts'],branch=sg,anchor=anchor,label=lab,k=float(k),residual4=r.tolist()))
  return out
 AA=seeds([173,174],[145,152,159,160,161,171,172]);BB=seeds([171,173],[145,152,159,160,161,172,174]);print('BANKS',len(AA),len(BB),flush=True)
 co=[(i,j)for i,a in enumerate(AA)for j,b in enumerate(BB)if O.R.classify(float(-np.array(a['residual4'])@m@np.array(b['residual4'])/2))not in('obtuse','angle_gap')];print('COMPATIBLE',len(co),flush=True)
 stat=collections.Counter();out=[]
 for pos,(i,j)in enumerate(co):
  ra=np.array(AA[i]['residual4']);rb=np.array(BB[j]['residual4']);Y=np.r_[Y0,[ra,rb]];q=np.r_[q0,2,2];ret=ret0+[175,176];U=Y[[173,175,176]];gram=U@m@U.T
  if abs(np.linalg.det(gram))<1e-10:stat['singular_frame']+=1;continue
  inv=np.linalg.inv(gram);ker=null_space(U@m)[:,0];kn=ker@m@ker
  for la,lb in itertools.product(range(2,7),repeat=2):
   f=np.array([0,-2*(0 if la==2 else math.cos(math.pi/la)),-2*(0 if lb==2 else math.cos(math.pi/lb))]);p=(inv@f)@U;rad=(2-p@m@p)/kn
   if rad<1e-10:stat['no_spacelike_branch']+=1;continue
   for sign in(-1,1):
    r=p+sign*math.sqrt(rad)*ker;stat['real_branches']+=1
    if max(r@m@V.T)>=-1e-8:stat['old_vertex_violation']+=1;continue
    coo=-(r@m@Y.T)/np.sqrt(2*q);labs=[O.R.classify(float(x))for x in coo];bad=[k for k in ret if labs[k]in('obtuse','angle_gap')]
    if bad:stat['gate_failure']+=1;continue
    YY=np.r_[Y,[r]];CT=np.r_[CT0,np.zeros((3,17))];whole=YY@m@YY.T+CT@J@CT.T
    if any(min(np.linalg.eigvalsh(whole[np.ix_(T+[173,k,177],T+[173,k,177])]))<=1e-8 for k in(175,176)):stat['closing_corner_not_PD']+=1;continue
    bnd=L.boundary(YY,m,ret+[177],T);closed=next(z for z in bnd['facets']if z['wall']==173)['compact'];z=dict(first_seed=i,second_seed=j,middle_labels=[la,lb],sign=sign,residual4=r.tolist(),middle_frame_gram=gram.tolist(),particular=p.tolist(),kernel=ker.tolist(),norm_radicand=float(rad),all_retained_native_labels=[labs[t]for t in ret],closed173=closed,boundary=bnd);out.append(z);print('SUCCESS',i,j,la,lb,'closed',closed,'open',len(bnd['open_edges']),'max',bnd['maximum_squared_radius'],flush=True)
  if pos%50==0:print('PROGRESS',pos,len(co),dict(stat),flush=True)
 paths=[f'scripts/{TAG}.py',S.INPUT,f'results/{S.TAG}.json','results/fresh_h20_replacement_outer_patch_pairs_20261017.json',f'scripts/{O.TAG}.py',f'scripts/{L.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_THREE_CAP_ROOF_WITH_COMPACT_HEPTAGON',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},base_first=0,base_second=77,retained_ids=ret0,first_seeds=AA,second_seeds=BB,compatible_pairs=[list(z)for z in co],statistics=dict(stat),candidates=out,scope='Numerical right-T three-cap roof on the173-wall two-cap branch. Seven finite right-T endpoint profiles, old anchor labels2..6, all35x32 endpoint pairs and25 middle-label pairs2..6. All175 prior native gates checked for each accepted middle root and complete target-H3 boundary recomputed. Eight compact173 heptagons found, but13/14 open edges remain. No exact promotion, all-label completeness, or finite-volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],'candidates',len(d['candidates']),flush=True)
