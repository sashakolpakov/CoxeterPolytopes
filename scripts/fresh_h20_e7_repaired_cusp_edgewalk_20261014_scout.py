#!/usr/bin/env python3
"""Four-support finite-template/continuous-acute screen and bounded contacts."""
import json,math,itertools
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from hashlib import sha256
import numpy as np
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_repaired_cusp_edgewalk_20261014_scout'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json';WALK='results/fresh_h20_e7_repaired_cusp_walk_20261014.json'
SELECTED=[1675,1673,1677,1723];TAIL=[7,8,10,12,16,24,30,60];ANCHORS=[2,3,4,5,6,8,10,12,16,24,30,60,'parallel']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def q5sign(a,b):
 if not b:return(a>0)-(a<0)
 if not a or(a>0)==(b>0):return(b>0)-(b<0)
 n=a*a-5*b*b;assert n
 return((a>0)-(a<0))*((n>0)-(n<0))
def run():
 base,walk=[json.loads((ROOT/p).read_text())for p in[BASE,WALK]]
 for doc in[base,walk]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 rad=np.sqrt(base['field']['base_square_root_basis']);dd=sum(float(F(q))*r for q,r in zip(base['field']['discriminant'],rad))**.5
 def val(a):return sum(float(F(q))*r for q,r in zip(a[0],rad))+dd*sum(float(F(q))*r for q,r in zip(a[1],rad))
 def vec(v):return np.array([val(a)for a in v])
 Q=np.array([[float(F(x))for x in row]for row in base['form']]);N=np.array([vec(v)for v in base['normals']]);H=vec(base['strict_center']);E=vec(base['strict_escape']['point']);nu=vec(walk['cusp_point']);norms=np.array(list(map(float,base['norms'])))
 NQ=N@Q;nn=NQ@nu;h=H@Q@nu;en=E@Q@nu;assert h<0 and en<0
 sqrt=np.sqrt(2*norms);squares={3:(F(1),F(0)),4:(F(2),F(0)),5:(F(3,2),F(1,2)),6:(F(3),F(0))}
 components=base['completed_cusp']['components'];fullgram={}
 for c in components:
  for i,row in zip(c['walls'],c['gram']):
   for j,g in zip(c['walls'],row):fullgram[i,j]=F(g)
 def g(i,j):return fullgram.get((i,j),F(0))
 projections=[{r['omitted']:r for r in cc}for cc in walk['component_projections']]
 edges=[];hits=[];stats=Counter();raw_trials=0;feasible=[];rejected=[]
 for ei in SELECTED:
  edge=walk['edge_germs'][ei];assert edge['type']=='open_opposite_end';ids=edge['support'];omitted=edge['omitted'];parts=[projections[i][j]for i,j in enumerate(omitted)]
  W=H-sum((vec(r['projection'])for r in parts));w2=W@Q@W;z=W-w2/(2*h)*nu;zn=NQ@z;ze=E@Q@z
  assert w2<0 and ze<0 and max(zn[[i for i in range(95)if i not in ids]])<-1e-8
  M=[[g(i,j)for j in ids]for i in ids];MI=B.inv(M);scales2=[F(2)if i in[17,79]else F(1)for i in ids];diag=[MI[k][k]*s for k,s in enumerate(scales2)];assert min(diag)>=F(1,2)
  scales=np.sqrt(list(map(float,scales2)));GI=np.array(MI,float)*np.outer(scales,scales);assert np.max(abs(GI@(N[ids]@Q@N[ids].T)-np.eye(19)))<1e-7
  rem=set(ids);comps=[]
  while rem:
   cc={min(rem)};todo=list(cc)
   for i in todo:
    for j in sorted(rem-cc):
     if g(i,j):cc.add(j);todo.append(j)
   rem-=cc;comps.append(sorted(cc))
  which={i:k for k,c in enumerate(comps)for i in c};tails=[c[0]for c in comps if len(c)==1]
  prof=[dict(contacts=[],cost=['0','0'])]
  for size in [1,2,3]:
   for inds in itertools.combinations(range(19),size):
    if len({which[ids[i]]for i in inds})<size or sum(diag[i]for i in inds)>=2:continue
    for labels in itertools.product([3,4,5,6],repeat=size)if size<3 else[(3,3,3)]:
     aa=sum(diag[i]*squares[m][0]for i,m in zip(inds,labels));bb=sum(diag[i]*squares[m][1]for i,m in zip(inds,labels))
     if q5sign(2-aa,-bb)>0:prof.append(dict(contacts=[[ids[i],m]for i,m in zip(inds,labels)],cost=list(map(str,[aa,bb]))))
  sampled=prof+[dict(contacts=[[i,m]],I2_tail_sample=True)for i in tails for m in TAIL]
  erec=dict(cusp_germ=ei,support=ids,omitted=omitted,components=comps,finite_pool_profiles=prof,symbolic_I2_tails=[dict(wall=i,labels='all integers m>=7')for i in tails],sampled_profile_count=len(sampled));edges.append(erec)
  other=[i for i in range(95)if i not in ids]
  for pi,pf in enumerate(sampled):
   b=np.zeros(19)
   for i,m in pf['contacts']:b[ids.index(i)]=-2*math.cos(math.pi/m)
   coeff=GI@b;p=coeff@N[ids];cost=b@coeff;q=2-cost
   assert q>1e-10
   pn=NQ@p;pe=E@Q@p;C=q*nn/(2*h);CE=q*en/(2*h)
   # Every nonincident constraint is a downward quadratic with positive
   # constant; every omitted cusp constraint is linear. Float bounds only.
   lower=0.;blocking=None
   for j in other:
    if j in base['completed_cusp']['incident']:
     threshold=max(0.,-pn[j]/zn[j])
    else:
     assert C[j]>0;disc=pn[j]**2-4*zn[j]*C[j];threshold=(pn[j]+math.sqrt(disc))/(-2*zn[j])
    if threshold>lower:lower=threshold;blocking=j
   upper=(pe+math.sqrt(pe*pe-4*ze*CE))/(-2*ze)
   stats['sampled_profiles']+=1
   if lower>=upper-1e-9:
    stats['acute_cut_interval_empty_numeric']+=1;rejected.append(dict(edge=ei,profile=pf,lower=lower,upper=upper,blocking_wall=blocking));continue
   stats['acute_cut_interval_survives_numeric']+=1;feasible.append(dict(edge=ei,profile=pf,lower=lower,upper=upper,blocking_wall=blocking))
   for j in other:
    for m in ANCHORS:
     contact=sqrt[j]if m=='parallel'else sqrt[j]*math.cos(math.pi/m)
     aa=zn[j];bb=pn[j]+contact;cc=C[j]
     if abs(cc)<1e-8 and j in base['completed_cusp']['incident']:
      a=-bb/aa
     else:
      disc=bb*bb-4*aa*cc
      if disc<0:continue
      a=(bb+math.sqrt(disc))/(-2*aa)
     if a<=0 or a<lower-1e-8 or a>=upper-1e-9:continue
     raw_trials+=1;r=p+a*z+(q/(2*a*h))*nu;pair=NQ@r
     if max(pair)>2e-7 or r@Q@E<=1e-8:continue
     labels=[];fail=False
     for cosine in -pair/sqrt:
      if abs(cosine)<2e-7:labels.append(2)
      elif cosine>=1-2e-7:labels.append('parallel_or_dotted')
      elif cosine<.5-2e-7:fail=True;break
      else:
       label=math.pi/math.acos(min(cosine,1));ii=round(label)
       if abs(label-ii)>2e-5:fail=True;break
       labels.append(ii)
     if fail:continue
     point=z-q/(2*a*a*h)*nu;xp=NQ@point
     if point@Q@point>=-1e-8 or point@Q@H>=-1e-8 or max(abs(xp[ids]))>2e-7 or max(xp[other])>=-1e-8 or abs(point@Q@r)>2e-7:continue
     assert abs(r@Q@r-2)<2e-6
     hits.append(dict(edge=ei,profile=pf,anchor=j,label=m,a=a,normal=r.tolist(),point=point.tolist(),point_norm=float(point@Q@point),escape_pairing=float(r@Q@E),labels=labels))
     print('NUMERIC_ACTUAL_FINITE_ENDPOINT',ei,pf['contacts'],j,m,hits[-1]['escape_pairing'],flush=True)
  print('edge',ei,'profiles',len(prof),'tails',len(tails),'cumulative',dict(stats),'trials',raw_trials,'hits',len(hits),flush=True)
 paths=[BASE,WALK,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py']
 return dict(status='BOUNDED_FOUR_SUPPORT_NUMERICAL_ENDPOINT_SCOUT',sources={p:sha(p)for p in paths},selected_germs=SELECTED,counts=dict(stats,numerical_contact_trials=raw_trials,candidates=len(hits)),edge_templates=edges,continuous_survivors=feasible,continuous_rejections=rejected,candidates=hits,scope='Exact finite-profile Schur acceptance over Qsqrt5, symbolic I2 tails retained; sampled tail labels and anchor labels bounded. Acuteness/cut intervals and all candidate metric gates are numerical, not rigorous all-parameter exclusions. Candidate acceptance requires future timelike endpoint on the specified outgoing cusp edge and all95 inequalities, plus cutting the frozen95 strict escape.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n');print(d['counts'])
