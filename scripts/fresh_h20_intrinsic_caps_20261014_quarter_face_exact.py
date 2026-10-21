#!/usr/bin/env python3
"""Exact selected quarter-face layer:587 intrinsic H20 facets, still infinite."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from fresh_h20_intrinsic_caps_20261014_golden_seed import ROOT,ZERO,ONE,PHI,q,add,sub,mul,scale,sign,enc,sha
TAG='fresh_h20_intrinsic_caps_20261014_quarter_face_exact'
SEED='results/fresh_h20_intrinsic_caps_20261014_golden_seed.json'
EXT='results/fresh_h20_intrinsic_caps_20261014_golden_extension.json'
SCOUT='results/fresh_h20_intrinsic_caps_20261014_quarter_face_scout.json'
HELPER='scripts/fresh_h20_intrinsic_caps_20261014_golden_seed.py'
def dec(x):return tuple(map(F,x))
def run():
 seed=json.loads((ROOT/SEED).read_text());ext=json.loads((ROOT/EXT).read_text());scout=json.loads((ROOT/SCOUT).read_text())
 for src in(seed,ext,scout):
  for p,h in src['sources'].items():assert sha(p)==h,p
 assert seed['status']=='PASS_EXACT_275_FACET_INTRINSIC_H20_GOLDEN_SEED_WITH_STRICT_ESCAPE'
 assert ext['status']=='PASS_EXACT_527_FACET_GOLDEN_EXTENSION_WITH_STRICT_ESCAPE'
 B=[[F(x)for x in r]for r in seed['spatial_Gram']]
 def sd(a,b):
  out=ZERO
  for i,row in enumerate(B):
   for j,x in enumerate(row):
    if x:out=add(out,scale(mul(a[i],b[j]),x))
  return out
 def ld(a,b):return sub(sub(sd(a,b),mul(a[19],b[20])),mul(a[20],b[19]))
 old=[[dec(x)for x in r]for r in seed['walls']+ext['new_walls']];assert len(old)==527
 oldnorm=[sd(w,w)if w[20]==ZERO else ld(w,w)for w in old]
 assert oldnorm==[ONE]*22+[q(4)]+[ONE]*504
 vertices=[w[:19]for w in old[23:275]]
 ids=seed['source_wall_order'];factors=[ids[:12],ids[12:15],ids[15:]];marks=list(map(F,seed['positive_affine_marks']))
 factor_marks=[marks[:12],marks[12:15],marks[15:]];basis=seed['spatial_basis_source_walls']
 # The scout uses factor-local vertex indices in this same literal wall order.
 raw=json.loads((ROOT/'results/fresh_h20_quartic_orbit_cap_gate_20261014.json').read_text())['seed_cusp']
 assert factors==[raw['unchanged_affine_factors'][0]['walls'],[23,50,200],raw['unchanged_affine_factors'][1]['walls']]
 assert factor_marks==[[F(1)]*12,[F(1)]*3,list(map(F,[3,1,1,2,1,2,2]))]
 fvertex=[]
 for k,f in enumerate(factors):
  inds=[j for j,w in enumerate(basis)if w in f]
  assert all(B[i][j]==0 for i in inds for j in range(19)if j not in inds)
  vv=[]
  for wall in f:
   matches=[v for v,row in zip(vertices,seed['finite_endpoints'])if row['omitted_source_walls'][k]==wall]
   z=[matches[0][j]if j in inds else ZERO for j in range(19)]
   assert all(v[j]==z[j]for v in matches for j in inds);vv.append(z)
  assert all(sum(factor_marks[k][i]*vv[i][j][a]for i in range(len(f)))==0 for j in range(19)for a in range(2))
  fvertex.append(vv)
 recipes=scout['selected_recipes'];assert len(recipes)==60 and len({json.dumps(r,sort_keys=True)for r in recipes})==60
 centers=[]
 for rec in recipes:
  z=[ZERO]*19
  for k,key in enumerate(['A11','triangle','E6']):
   subverts=rec[key];assert subverts and sorted(set(subverts))==subverts and all(type(i)is int and 0<=i<len(factors[k])for i in subverts)
   den=4*sum(factor_marks[k][i]for i in subverts)
   for i in subverts:
    z=[add(a,scale(b,factor_marks[k][i]/den))for a,b in zip(z,fvertex[k][i])]
  centers.append(z)
 assert len(set(tuple(z)for z in centers))==60
 R=F(3,2);R2=q(R*R)
 new=[z+[scale(sub(sd(z,z),R2),F(1,2)),ONE]for z in centers]
 assert all(ld(w,w)==R2 for w in new)
 def gate(p,den):
  if p==ZERO:return 'right'
  assert sign(p)<0
  if p==q(-den/F(2)):return 'label3'
  if p==scale(PHI,-den/F(2)):return 'label5'
  s=sign(add(p,q(den)))
  assert s<=0,(p,den)
  return 'parallel'if s==0 else'dotted'
 counts=Counter();seed_counts=Counter();cross_counts=Counter();new_counts=Counter()
 for r in new:
  for j,w in enumerate(old):
   label=gate(ld(r,w),R*(2 if j==22 else 1));counts[label]+=1
   if j<22:seed_counts[label]+=1
   else:assert label=='dotted';cross_counts[label]+=1
 for a,b in combinations(new,2):
  label=gate(ld(a,b),R*R);counts[label]+=1;new_counts[label]+=1
 assert sum(counts.values())==60*527+60*59//2
 center=[dec(x)for x in seed['common_strict_center']['point']]
 cps=[ld(center,r)for r in new];assert all(sign(p)<0 for p in cps)
 # Every previously supplied actual-facet certificate survives.
 oldpoints=[row['point']for row in seed['finite_endpoints']]+[seed['central_facet_point']['point']]+[row['point']for row in ext['new_actual_facet_points']]
 for row in oldpoints:
  x=list(map(dec,row));assert all(sign(ld(x,r))<0 for r in new)
 facet_points=[]
 for i,r in enumerate(new):
  x=r[:19]+[add(r[19],R2),ONE]
  assert ld(x,x)==scale(R2,-1)and sign(ld(x,center))<0
  ps=[ld(x,w)for w in old+new]
  assert [j for j,p in enumerate(ps)if p==ZERO]==[527+i]and all(sign(p)<=0 for p in ps)
  facet_points.append({'point':x,'norm':scale(R2,-1),'only_incident_wall':527+i})
 # Exact ball coverage of all252 sampled quarter-vertex positions.
 coverage=[]
 for v in vertices:
  z=[scale(x,F(1,4))for x in v];n=sd(z,z);x=z+[scale(sub(n,ONE),F(1,2)),ONE]
  covered=[]
  for i,c in enumerate(centers):
   d=sd([sub(a,b)for a,b in zip(z,c)],[sub(a,b)for a,b in zip(z,c)])
   if sign(sub(d,R2))<=0:covered.append(i);assert sign(ld(x,new[i]))>0
  assert covered;coverage.append(covered)
 oldescape=[dec(x)for x in ext['strict_remaining_escape']['point']]
 oldcuts=[ld(oldescape,r)for r in new];assert any(sign(p)>0 for p in oldcuts)
 z=[scale(x,F(1,8))for x in vertices[0]];n=sd(z,z);escape_norm=q(F(1,2));escape=z+[scale(sub(n,escape_norm),F(1,2)),ONE]
 assert ld(escape,escape)==escape_norm and sign(ld(escape,center))<0
 ep=[ld(escape,r)for r in old+new];assert all(sign(p)<0 for p in ep)
 # All spherical caps remain strictly dotted to the central hemisphere.
 assert all(gate(ld(old[22],w),2*(R if i>=527 else 1))=='dotted'for i,w in enumerate(old+new)if i>=23)
 firstnorm=sd(vertices[0],vertices[0]);assert sign(sub(firstnorm,q(4)))>0
 paths=[SEED,EXT,SCOUT,HELPER,'results/fresh_h20_quartic_orbit_cap_gate_20261014.json',str(Path(__file__).resolve().relative_to(ROOT))]
 total=Counter(ext['all527_pair_counts']);total.update(counts)
 return enc({'status':'PASS_EXACT_587_FACET_QUARTER_FACE_LAYER_WITH_STRICT_ESCAPE','sources':{p:sha(p)for p in paths},'field':seed['field'],'coordinates':seed['coordinates'],'spatial_Gram':B,'Lorentz_signature':[20,1],'retained_wall_count':527,'recipes':recipes,'radius':R,'center_scale':F(1,4),'new_centers':centers,'new_walls':new,'new_norms':[R2]*60,'new_pair_counts':dict(counts),'new_seed_pairs':dict(seed_counts),'new_old_spherical_pairs':dict(cross_counts),'new_new_pairs':dict(new_counts),'all587_pair_counts':dict(total),'strict_center':{'point':center,'norm':q(-9),'new_pairings':cps},'old252_actual_finite_endpoints_retained':True,'all527_old_facets_retained':True,'new_actual_facet_points':facet_points,'quarter_vertex_ball_cover':coverage,'previous_escape_cut':{'point':oldescape,'new_pairings':oldcuts},'strict_remaining_escape':{'point':escape,'norm':escape_norm,'center_pairing':ld(escape,center),'all587_pairings':ep},'central_cover_component_obstruction':{'central_radius':2,'first_vertex_square':firstnorm,'all_other_spherical_caps_strictly_dotted':True,'statement':'The central covering ball has a boundary point strictly inside the product alcove, along the first-vertex ray. Every other retained cap ball is strictly disjoint from it. Any finite collection of further caps strictly dotted to the central cap leaves an uncovered open set near that boundary point. Thus such additions cannot give finite volume; finite-contact bridges or replacement are needed.'},'scope':'Exact60 selected face-center recipes only, not the entire numerical bank. All587 inequalities are actual facets;252 prior finite vertices survive and252 quarter probes are covered. A new strict spacelike escape proves infinite volume. No full corner census or finite-volume claim.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['new_pair_counts']);print(out['new_new_pairs'])
