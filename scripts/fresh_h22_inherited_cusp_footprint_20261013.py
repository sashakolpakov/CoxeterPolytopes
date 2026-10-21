#!/usr/bin/env python3
"""Exact inherited-cusp footprint and retained-core facet lower bound.

Uses the independently checked complete H21 census, without expanding its
1,095,990 literal cusp rays.  New-coordinate radicals use the separate sparse
squarefree arithmetic from the multiquadratic review, not producer intervals.
"""
import argparse,hashlib,json,math,time
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
import numpy as np
from sympy import Matrix
import fresh_h22_multiquadratic_20261013_review as A
P=Path(__file__).resolve().parents[1];TAG='fresh_h22_inherited_cusp_footprint_20261013';OUT=P/'results'/f'{TAG}.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json';WALK_AUDIT='results/fresh_h21_source_corner_walk_20261013_independent.json';MQ='results/fresh_h22_multiquadratic_20261013.json';MQ_AUDIT='results/fresh_h22_multiquadratic_20261013_review.json';REC='data/enumerations/fresh_h22_incidence_prism_20261013_recipes.json';COUPLED='results/fresh_h22_pi7_coupled_caps_20261013.json';RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';CORE='results/fresh_h22_incidence_core_20261013.json';BANK='data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def primitive(v):
 den=math.lcm(*(F(x).denominator for x in v));v=[int(F(x)*den)for x in v];g=math.gcd(*v);assert g;return tuple(x//g for x in v)
def run():
 started=time.monotonic();walk=read(WALK);audit=read(WALK_AUDIT);mq=read(MQ);mqa=read(MQ_AUDIT)
 assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS'and not walk['pending_orbit_ids']and audit['queue_closed'];assert audit['producer_sha256']==sha(WALK)
 for obj in(walk,mq,mqa):
  for p,h in obj['sources'].items():assert sha(p)==h,p
 assert mqa['sources'][MQ]==sha(MQ);assert sha(walk['group_array'])==walk['group_sha256'];perms=np.load(P/walk['group_array'])['permutations'];assert perms.shape==(241920,210)
 oldcusps=[r for r in walk['corner_orbits']if r['type']=='ideal'];assert len(oldcusps)==18;total=sum(r['orbit_size']for r in oldcusps);assert total==audit['literal_cusps']==1095990;orbit_lookup={tuple(r['canonical_active_walls']):r for r in oldcusps}
 raw=read(RAW)['coxeter_input'];Q=[[-int(x)for x in row]for row in raw['basis_gram']];M=[r+[0]for r in Q]+[[0]*22+[2]];W=[list(map(int,row))for row in raw['wall_coordinates_as_columns_transposed']];COV=[[dot(row,w)for row in Q]for w in W];core=read(CORE);h=core['center'];hc=[dot(row,h)for row in Q];assert dot(h,hc)==-14
 # h is a positive multiple of the sum of all source normals.  Therefore
 # it is strictly negative on every nonzero vector of the source cone.
 sw=[sum(w[j]for w in W)for j in range(22)];ii=next(i for i,x in enumerate(h)if x);center_multiple=F(sw[ii],h[ii]);assert center_multiple>0 and sw==[center_multiple*x for x in h]and Matrix(W).rank()==22
 points=core['incidence_parts'][0];vp=core['old_projections'][1];parent=[]
 for r in read(BANK)['candidates']:
  if r['twice_height_square']!=4:continue
  mask=[(r['mask']>>i)&1 for i in range(21)];assert sum(mask)==r['weight'];u=[r['A']*vp[j]-sum(mask[i]*W[points[i]][j]for i in range(21))for j in range(22)];assert u==r['twice_old_projection']and all(x%2==0 for x in u);parent.append([x//2 for x in u])
 assert len(parent)==301
 # Every one of these301 rational null directions has both height signs.
 nullrows=[];seen=set();covered_counts=Counter()
 def add_null(proj,origin,full_plus,full_minus):
  supports={j for x in proj for j in x};assert len(supports)==1;rad=next(iter(supports));coeff=[x.get(rad,F(0))for x in proj];ray=primitive(coeff);q=dot(ray,[dot(row,ray)for row in Q]);ps=[dot(ray,c)for c in COV];assert q==0 and dot(ray,hc)<0 and max(ps)==0 and ray not in seen;seen.add(ray)
  active=[i for i,x in enumerate(ps)if not x];assert Matrix([W[i]for i in active]).rank()==21
  aa=np.array(active,dtype=np.uint8);images=np.sort(perms[:,aa],axis=1);indices=np.arange(len(perms))
  for j in range(len(active)):
   vals=images[indices,j];indices=indices[vals==vals.min()]
  gi=int(indices[0]);key=tuple(map(int,images[gi]));record=orbit_lookup[key];assert sum(c['rank']for c in record['affine_components'])==20;covered_counts[record['id']]+=1
  # These two actual normals are parallel at this cusp, orthogonal to
  # all old cusp walls, and have opposite nonzero vertical projections.
  covp=A.mv(M,full_plus);covm=A.mv(M,full_minus);qn=A.pair(full_plus,covp);assert A.sign(qn)>0 and A.pair(full_minus,covm)==qn and A.pair(full_plus,covm)==A.scale(qn,-1)
  assert A.sign(full_plus[-1])>0 and A.sign(full_minus[-1])<0
  assert full_plus[:-1]==full_minus[:-1]and full_minus[-1]==A.scale(full_plus[-1],-1)
  for i in active:
   assert A.pair(full_plus,A.ratvec(COV[i]+[0]))=={}and A.pair(full_minus,A.ratvec(COV[i]+[0]))=={}
  factor=F(coeff[ii],ray[ii])if ray[ii]else next(F(c,r)for c,r in zip(coeff,ray)if r);assert factor>0 and coeff==[factor*x for x in ray]
  nullrows.append({'origin':origin,'primitive_old_cusp':list(ray),'old_active_walls':active,'old_cusp_orbit_id':record['id'],'canonical_group_index':gi,'projection_positive_multiple':A.enc({rad:factor}),'cap_squared_norm':A.enc(qn),'opposite_height_parallel_pair':True,'old_quotient_rank':20,'new_link_rank':21,'new_A1_positive_kernel':[1,1]})
 for i,u in enumerate(parent):add_null(A.ratvec(u),{'parent_null_type':i,'positive_wall':210+i,'negative_wall':511+i},A.ratvec(u+[1]),A.ratvec(u+[-1]))
 recs=read(REC)['recipes'];normals={};norms={};timelike=[]
 for k in mq['selected_original_types']:
  rec=recs[k];u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=dot(u,[dot(row,u)for row in M]);qv=dot(v,[dot(row,v)for row in M]);assert dot(u,[dot(row,v)for row in M])==0 and qv>0;rad=q/qv-3*q*q*qu/(4*qv)
  r=[A.add(A.scale(A.sqrtq(3),-q*x/2),A.scale(A.sqrtq(rad),rec['radical_branch']*y))for x,y in zip(u,v)];normals[k]=r if A.sign(r[-1])>0 else A.flip(r);norms[k]=q
 a=list(map(F,recs[17]['u']));b=list(map(F,recs[17]['v']));Z=[A.add(A.scale(A.sqrtq(31),x/3),A.scale(A.sqrtq(3),-3*y/2))for x,y in zip(a,b)];Z=Z if A.sign(A.pair(Z,A.ratvec(hc+[0])))<0 else A.vs(Z,-1)
 u,v=[list(map(F,x))for x in read(COUPLED)['rational_duals']];normals[23]=A.vs(A.vadd(A.vs(Z,209),A.ratvec([-241*x for x in u])),F(1,120));normals[24]=A.vm(A.vadd(A.vs(Z,241),A.ratvec([-418*x for x in v])),A.scale(A.sqrtq(29281),F(1,29281)));norms[23]=norms[24]=F(2)
 ids=mq['selected_original_types']+[23,24];new_null=[]
 for k in ids:
  r=normals[k];assert A.pair(r,A.mv(M,r))==A.C(norms[k])and A.sign(A.pair(r,A.ratvec(hc+[0])))<0
  projection=r[:-1];pq=A.pair(projection,A.mv(Q,projection));assert A.sign(pq)<=0
  assert all(A.sign(A.pair(projection,A.ratvec(c)))<=0 for c in COV)
  if not pq:
   new_null.append(k);add_null(projection,{'added_type':k,'positive_wall':812+2*ids.index(k),'negative_wall':813+2*ids.index(k)},r,A.flip(r))
  else:
   assert A.sign(pq)<0;timelike.append({'type':k,'old_projection_squared_norm':A.enc(pq),'old_center_pairing':A.enc(A.pair(projection,A.ratvec(hc)))})
 assert new_null==[4,5]and len(nullrows)==303 and len(timelike)==22
 # The two retained signs exclude every nonzero pure vertical added wall:
 # its pairing with a retained cap of the same sign is strictly positive.
 assert A.sign(normals[0][-1])>0
 missing=total-len(nullrows);lower=2*missing;assert missing==1095687 and lower==2191374
 counts=[{'old_cusp_orbit_id':r['id'],'orbit_size':r['orbit_size'],'bounded_literal_cusps':covered_counts[r['id']],'unbounded_literal_cusps':r['orbit_size']-covered_counts[r['id']]}for r in oldcusps];assert min(r['unbounded_literal_cusps']for r in counts)>=0
 sources={p:sha(p)for p in(WALK,WALK_AUDIT,MQ,MQ_AUDIT,REC,COUPLED,RAW,CORE,BANK,'scripts/'+TAG+'.py','scripts/fresh_h22_multiquadratic_20261013_review.py',walk['group_array'])}
 return {'status':'PASS_EXACT_INHERITED_CUSP_FOOTPRINT_AND_RETAINED_CORE_LOWER_BOUND','sources':sources,'old_cusp_orbits':18,'total_old_literal_cusps':total,'bounded_old_literal_cusps':len(nullrows),'unbounded_old_literal_cusps':missing,'old_vertical_facets':210,'rational_null_projection_caps':602,'new_null_projection_caps':4,'new_timelike_projection_caps':44,'pure_vertical_caps':0,'covered_null_directions':nullrows,'timelike_positive_height_types':timelike,'per_old_cusp_orbit':counts,'minimum_additional_facets_in_finite_completion':lower,'minimum_total_facets_in_finite_completion':860+lower,'source_center_positive_sum_multiple':str(center_multiple),'scope':'Exactly303 inherited cusp links are bounded in this frozen860-facet chamber. Every other inherited cusp has an unbounded new-coordinate link direction. Any finite-facet acute-compatible finite-volume completion retaining all860 actual facets needs at least two additional facets per missing cusp, because a nonzero causal old projection can be incident to only one old null ray and pure vertical caps conflict with retained height signs. This is a lower bound for retaining this core, not a nonexistence theorem and not a bound on other constructions.','elapsed_seconds':time.monotonic()-started}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run()
 if a.verify:
  old=read(str(OUT.relative_to(P)));assert {k:v for k,v in d.items()if k!='elapsed_seconds'}=={k:v for k,v in old.items()if k!='elapsed_seconds'}
 else:OUT.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'bounded',d['bounded_old_literal_cusps'],'missing',d['unbounded_old_literal_cusps'],'additional',d['minimum_additional_facets_in_finite_completion'],'seconds',d['elapsed_seconds'])
