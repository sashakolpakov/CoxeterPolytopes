#!/usr/bin/env python3
"""Bounded singleton151 endpoint angles, including cross-bank finite pairs."""
import argparse,hashlib,json,math,itertools
import numpy as np
import sympy as sp
import fresh_h20_e8_partner_endpoints_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_e8_endpoint_belt_20261017_param151';FIXED='results/fresh_h20_e8_endpoint_belt_20261017_finite158.json'
def run():
 G,Y,CT,J,q,H,live,verts=P.geometry();C=P.load('results/fresh_h20_e8_cusp_completion_20261017.json');K=P.load('results/fresh_h20_e8_f4_cusps_20261017.json');ids=C['core_ids'];M=np.array([[float(sp.sympify(x))for x in r]for r in C['core_gram']]);MI=np.linalg.inv(M);added=[]
 for cu in C['cusps']:
  nu=np.array([float(sp.sympify(x))for x in cu['null_coefficients']])
  for row in cu['base_new_root_coefficients']:added.append(np.array([float(sp.sympify(x))for x in row])+K['selected_common_height']*nu)
 A=np.array(added);Y=np.r_[Y,A@Y[ids]];CT=np.r_[CT,A@CT[ids]];q=np.r_[q,[2.]*6];B=Y[ids]@G['m']@Y.T+CT[ids]@J@CT.T
 out=[];rejects=[]
 for ci,edge in enumerate(P.EDGES):
  missing=next(j for j in(159,160,161)if j not in edge);mi=ids.index(missing);ix=[i for i in range(21)if i!=mi];inv=np.linalg.inv(M[np.ix_(ix,ix)]);V=-MI[:,mi]/math.sqrt(-MI[mi,mi]);fixed=set(G['T']+list(edge)+[158])
  for label151,label158 in itertools.product(range(3,65),list(range(2,33))+['parallel']):
   f=np.zeros(20);f[ix.index(ids.index(151))]=-2*math.cos(math.pi/label151);f[ix.index(ids.index(158))]=-2*(1 if label158=='parallel'else 0 if label158==2 else math.cos(math.pi/label158));p=np.zeros(21);p[ix]=inv@f;rad=p@M@p-2
   if rad<=1e-10:continue
   r=p+math.sqrt(rad)*V;ps=r@B;others=[j for j in range(168)if j not in fixed];cs=-ps[others]/np.sqrt(2*q[others]);minimum=float(min(cs))
   if minimum<1-1e-7:rejects.append(dict(edge_index=ci,label151=label151,label158=label158,blocker=others[int(np.argmin(cs))],cosine=minimum));continue
   out.append(dict(edge_index=ci,edge=list(edge),label151=label151,label158=label158,core_coefficients=r.tolist(),residual4=(r@Y[ids]).tolist(),T_pairings=(r@CT[ids]).tolist(),minimum_other_cosine=minimum))
 print('PARAM151',len(out),'accepted',len(rejects),'rejected',flush=True)
 F=P.load(FIXED);bank=F['candidates']+out;RR=np.array([r['residual4']for r in bank]);DD=np.array([r['T_pairings']for r in bank]);ci=np.array([r['edge_index']for r in bank]);ar=np.arange(len(bank));pairs=[];dotted=0;finitegap=0
 for lo in range(len(F['candidates']),len(bank),128):
  hi=min(lo+128,len(bank));co=-(RR[lo:hi]@G['m']@RR.T+DD[lo:hi]@J@DD.T)/2;mask=(ar[None,:]<ar[lo:hi,None])&(ci[None,:]!=ci[lo:hi,None]);dotted+=int(np.sum(mask&(co>1+1e-7)))
  ii,jj=np.where(mask&(co>=-1e-7)&(co<=1+1e-7))
  for a,j in zip(ii,jj):
   i=lo+int(a);v=float(co[a,j])
   if abs(v)<1e-7:lab=2
   elif abs(v-1)<1e-7:lab='parallel'
   else:
    z=math.pi/math.acos(max(-1,min(1,v)));lab=round(z)
    if lab<3 or abs(z-lab)>1e-5:finitegap+=1;continue
   pairs.append(dict(i=i,j=int(j),label=lab,cosine=v))
 paths=[f'scripts/{TAG}.py',FIXED,'results/fresh_h20_e8_cusp_completion_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json',f'scripts/{P.TAG}.py']
 return dict(status='BOUNDED_SINGLETON151_ENDPOINT_CONTACT_DIAGNOSTIC',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},label151_domain=[3,64],label158_domain=[2,32,'parallel'],candidates=out,rejections=rejects,prior_fixed_count=len(F['candidates']),new_finite_or_parallel_pairs=pairs,new_dotted_pairs=dotted,new_noncoxeter_gaps=finitegap,scope='Numerical unique future norm-completion from each actual rank20 old corner. All168 unprescribed contacts are required dotted. Labels151=3..64 and158=2..32/parallel are bounded; no unbounded-family exclusion. Cross-pair check includes the prior fixed-profile finite158 bank. Finite-label numerical coincidences need exact promotion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
