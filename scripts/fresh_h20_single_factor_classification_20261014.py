#!/usr/bin/env python3
"""Complete all-retained single-affine-factor cusp-section classification.

Standard-library exact arithmetic only. Completeness of the star-cut reduction
and continuous interpolation is proved in the accompanying report.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from pathlib import Path
import argparse, hashlib, json, math

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_single_factor_classification_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
SUM='results/fresh_h20_single_factor_links_20261013.json'
DIFF='results/fresh_h20_cusp_sections_20261013.json'
RECIPES=f'results/{TAG}_recipes.json'
ALLOWED={F(0),F(1,4),F(1,2),F(3,4),F(1)}
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,x):return [dot(r,x)for r in A]
def pair(a,b,Q):return dot(a,mv(Q,b))
def inv(A):
 n=len(A);M=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  i=next(i for i in range(j,n)if M[i][j]);M[i],M[j]=M[j],M[i]
  q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 A=[[F(x)for x in r]for r in A];out=[]
 while A:
  q=A[0][0];assert q>0;out.append(q)
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))]for i in range(1,len(A))]
 return out
def canon(v):
 q=next(x for x in v if x);return tuple(F(x)/q for x in v)
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def star_partitions(A):
 n=len(A);N=[[j for j in range(n)if j!=i and A[i][j]]for i in range(n)]
 assert all(len(set(N[i])&set(N[j]))<2 for i,j in combinations(range(n),2))
 todo=[0]
 for i in todo:
  todo +=[j for j in N[i]if j not in todo]
 assert len(todo)==n
 parts=set();trials=0
 for center in range(n):
  for k in range(1,len(N[center])+1):
   for neighbors in combinations(N[center],k):
    trials+=1;cut={tuple(sorted((center,j)))for j in neighbors};col={0:0};queue=[0];ok=True
    for i in queue:
     for j in N[i]:
      value=col[i]^int(tuple(sorted((i,j)))in cut)
      if j in col:
       if col[j]!=value:ok=False;break
      else:col[j]=value;queue.append(j)
     if not ok:break
    if not ok:continue
    assert len(col)==n
    S=tuple(i for i in range(n)if col[i]==0);T=tuple(i for i in range(n)if col[i])
    if min(len(S),len(T))>=2:parts.add(S)
 return sorted(parts),trials
def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=[[-int(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[mv(Q,w)for w in W];G=[[dot(w,c)for c in C]for w in W]
 assert len(W)==210 and len(Q)==22 and all(G[i][j]<=0 for i,j in combinations(range(210),2))
 h=[sum(w[d]for w in W)for d in range(22)];g=math.gcd(*h);h=[x//g for x in h]
 assert pair(h,h,Q)==-14 and max(dot(h,c)for c in C)<0
 walk=json.loads((ROOT/WALK).read_text())
 assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS'and not walk['pending_orbit_ids']
 assert walk['sources'][RAW]==sha(RAW)
 cusps=[r for r in walk['corner_orbits']if r['type']=='ideal'];assert len(cusps)==18
 oldsum=json.loads((ROOT/SUM).read_text());olddiff=json.loads((ROOT/DIFF).read_text())
 for packet in[oldsum,olddiff]:
  assert packet['sources'][RAW]==sha(RAW)
 def key(e,D):
  assert pair(D,D,Q)==pair(e,D,Q)==0 and pair(e,e,Q)>0
  active=[i for i in range(210)if dot(D,C[i])==0]
  return canon(D),canon([dot(e,C[i])for i in active])
 sumkeys={key(list(map(F,r['normal_constant'])),list(map(F,r['normal_shift'])))for r in oldsum['sum_planes']}
 diffkeys={key(list(map(F,r['base_normal'])),list(map(F,r['cusp'])))for r in olddiff['pencils']}
 assert len(sumkeys)==79 and len(diffkeys)==59
 recipes=json.loads((ROOT/RECIPES).read_text());assert len(recipes)==28
 recipe_map={(r['cusp'],tuple(r['factor']),tuple(r['partition_S']),tuple(r['active_part'])):r for r in recipes}
 assert len(recipe_map)==28
 used_recipes=set();faces={};failed_faces=set()
 def face(ids):
  ids=tuple(sorted(set(ids)))
  if ids in faces:return faces[ids]
  if ids in failed_faces:return None
  A=[[G[i][j]for j in ids]for i in ids]
  try:piv=pd(A)
  except AssertionError:failed_faces.add(ids);return None
  B=inv(A);assert min(x for r in B for x in r)>=0
  cc=mv(B,[dot(h,C[i])for i in ids])
  p=[F(h[d])-sum(c*W[i][d]for c,i in zip(cc,ids,strict=True))for d in range(22)]
  ps=[dot(p,c)for c in C];norm=pair(p,p,Q);future=pair(p,h,Q)
  assert norm<0 and future<0 and max(ps)==0 and [i for i,x in enumerate(ps)if not x]==list(ids)
  out={'walls':list(ids),'point':p,'norm':norm,'future_pair':future,'maximum_other_pairing':max(x for x in ps if x),'LDL':piv}
  faces[ids]=out;return out
 def global_witness(act,rest,marks,D,half=None):
  e=[sum(marks[i]*W[i][d]for i in act)for d in range(22)];q=pair(e,e,Q);ep=[dot(e,c)for c in C]
  if half is None:
   choices=((i,j)for j in range(210)if j not in act+rest and all(G[j][k]==0 for k in rest)and any(G[j][k]for k in act)for i in act if ep[i]>0 and G[i][j]==0)
  else:choices=[tuple(half['ridge'])]
  for i,j in choices:
   assert i in act and j not in act+rest and ep[i]>0 and G[i][j]==0
   assert all(G[j][k]==0 for k in rest)and ep[j]==dot(D,C[j])<0
   L=face(rest+[i,j]);assert L is not None
   if half is None:
    R=face(act+[j])
    if R is None:continue
    if set(L['walls'])&set(R['walls'])!={i,j}:continue
    rp=R['point'];endpoint=F(0)
   else:
    rp=list(map(F,half['midpoint']));rr=[dot(rp,c)for c in C]
    assert max(rr)==0 and [k for k,x in enumerate(rr)if not x]==sorted([i,j])
    assert pair(rp,rp,Q)==F(half['norm'])<0 and pair(rp,h,Q)<0
    assert pair(rp,e,Q)==pair(rp,D,Q)/2
    endpoint=F(-1,2);R={'point':rp,'norm':pair(rp,rp,Q),'future_pair':pair(rp,h,Q),'zero_walls':[i,j],'maximum_other_pairing':max(x for x in rr if x)}
   lp=L['point'];lc=-pair(lp,D,Q);rc=-pair(rp,D,Q);assert lc>0 and rc>0
   lp=[x/lc for x in lp];rp=[x/rc for x in rp]
   assert pair(lp,D,Q)==pair(rp,D,Q)==-1 and pair(lp,e,Q)==-1 and pair(rp,e,Q)==endpoint
   assert pair(lp,rp,Q)<0
   u=1+endpoint;a=-ep[j];pcoef=ep[i]*a/q;ni=F(G[i][i])-ep[i]**2/q
   detlower=ni*G[j][j]-a*a*G[i][i]*u*u/q
   assert pcoef>0 and ni>0 and detlower>0
   return {'active_part':act,'other_part':rest,'normal_constant':e,'normal_norm':q,'ridge':[i,j],
     'parameter_interval':[-1,endpoint],'upper_endpoint_included':half is not None,
     'left_ambient_face':L,'right_ambient_face_or_midpoint':R,
     'normalized_left':lp,'normalized_right':rp,
     'projected_pairing_equals_coefficient_times_1_plus_t':pcoef,
     'projected_first_norm':ni,'ridge_determinant_lower_bound':detlower}
  return None
 factors=[];accepted=[];rejected=[];seen=set();counts=Counter()
 for cusp in cusps:
  delta=list(map(F,cusp['primitive_ray']));assert pair(delta,delta,Q)==0 and pair(delta,h,Q)<0
  vals=[dot(delta,c)for c in C];assert max(vals)==0 and [i for i,x in enumerate(vals)if not x]==cusp['active_walls']
  for comp in cusp['affine_components']:
   ids=comp['walls'];n=len(ids);m=list(map(F,comp['positive_kernel']));mu=F(comp['ray_multiple']);A=[[F(G[i][j])for j in ids]for i in ids]
   assert comp['rank']==n-1 and min(m)>0 and mu>0 and all(dot(r,m)==0 for r in A)
   D=[sum(m[k]*W[i][d]for k,i in enumerate(ids))for d in range(22)]
   assert D==[mu*x for x in delta]
   for omit in range(n):pd([[A[i][j]for j in range(n)if j!=omit]for i in range(n)if i!=omit])
   parts,trials=star_partitions(A);counts['factors']+=1;counts['star_trials']+=trials
   factors.append({'cusp':cusp['id'],'walls':ids,'marks':m,'null_sum':D,'star_trials':trials,'rank_one_partitions':len(parts)})
   for S in parts:
    T=tuple(i for i in range(n)if i not in S);ss=[ids[i]for i in S];tt=[ids[i]for i in T];counts['rank_one_partitions']+=1
    e=[sum(m[i]*W[ids[i]][d]for i in S)for d in range(22)];q=pair(e,e,Q);assert q>0 and pair(e,D,Q)==0
    ep=[dot(e,C[i])for i in ids];B=[[A[i][j]-ep[i]*ep[j]/q for j in range(n)]for i in range(n)]
    assert all(B[i][j]==0 for i in S for j in T)
    assert min(B[i][i]for i in range(n))>0 and all(B[i][j]<=0 for i,j in combinations(range(n),2))
    for block in[S,T]:
     assert all(sum(B[i][j]*m[j]for j in block)==0 for i in block)
     pd([[B[i][j]for j in block[:-1]]for i in block[:-1]])
    angle=[(i,j,B[i][j]**2/(B[i][i]*B[j][j]))for i,j in combinations(range(n),2)]
    bad=[(ids[i],ids[j],v)for i,j,v in angle if v not in ALLOWED]
    base={'cusp':cusp['id'],'factor':ids,'S':ss,'T':tt,'normal_constant':e,'null_shift':D,'normal_norm':q}
    if bad:
     counts['local_forbidden_angle']+=1;rejected.append(base|{'bad_pair':list(bad[0][:2]),'rational_squared_cosine':bad[0][2]});continue
    counts['locally_Coxeter']+=1;k=key(e,D);assert k not in seen and k not in diffkeys;seen.add(k)
    old=k in sumkeys;counts['old_leaf_sum_planes'if old else'additional_planes']+=1
    marks=dict(zip(ids,m,strict=True));witness=None
    for act,rest in[(ss,tt),(tt,ss)]:
     witness=global_witness(act,rest,marks,D)
     if witness:break
    if witness:
     counts['one_full_interval_witness']+=1;proof=[witness]
    else:
     proof=[]
     for act,rest in[(ss,tt),(tt,ss)]:
      rk=(cusp['id'],tuple(ids),tuple(ss),tuple(act));r=recipe_map[rk];used_recipes.add(rk)
      assert r['other_part']==rest
      proof.append(global_witness(act,rest,marks,D,half=r))
     assert len(proof)==2 and all(proof);counts['two_half_interval_witnesses']+=1
    accepted.append(base|{'old_leaf_sum_plane':old,'projected_local_Gram':B,'Coxeter_squared_cosine_counts':dict(Counter(str(v)for i,j,v in angle)),'global_actual_obtuse_ridge_certificates':proof})
 assert counts['factors']==53 and counts['rank_one_partitions']==287 and counts['local_forbidden_angle']==119
 assert counts['locally_Coxeter']==168 and counts['old_leaf_sum_planes']==77 and counts['additional_planes']==91
 assert counts['one_full_interval_witness']==154 and counts['two_half_interval_witnesses']==14
 assert used_recipes==set(recipe_map)
 return enc({'status':'PASS_COMPLETE_ALL_RETAINED_SINGLE_FACTOR_CLASSIFICATION_AND_GLOBAL_EXCLUSION',
  'sources':{p:sha(p)for p in[RAW,WALK,SUM,DIFF,RECIPES,f'scripts/{TAG}.py']},
  'counts':dict(counts)|{'cusps':18,'exact_finite_ambient_faces':len(faces),'exact_middle_ridge_points':len(recipes)},
  'factor_census':factors,'locally_rejected_partitions':rejected,'locally_Coxeter_planes':accepted,
  'scope':'Every real single-factor cusp section retaining all n+1 factor inequalities as distinct actual Euclidean-link facets is excluded, for all53 factors at the18 source cusp representatives. The one-omitted-facet case and other chambers are not classified.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();path=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert json.loads(path.read_text())==out
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['counts'])
