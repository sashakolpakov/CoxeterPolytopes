#!/usr/bin/env python3
"""Replace one radius-one vertex cap by an exact radius-four endpoint cap."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from fresh_h20_intrinsic_caps_20261014_golden_seed import ROOT,ZERO,ONE,PHI,q,add,sub,mul,scale,sign,enc,sha,inv,pd
TAG='fresh_h20_intrinsic_caps_20261014_radius4_replacement'
BASE='results/fresh_h20_intrinsic_caps_20261014_'
def dec(x):return tuple(map(F,x))
def run():
 paths=[BASE+s+'.json'for s in ['golden_seed','golden_extension','quarter_face_exact','endpoint','annular_gate','annulus_endpoint_scout']]
 seed,ext,quarter,endpoint,annulus,scout=[json.loads((ROOT/p).read_text())for p in paths]
 for src in(seed,ext,quarter,endpoint,annulus,scout):
  for p,h in src['sources'].items():assert sha(p)==h,p
 B=[[F(x)for x in r]for r in seed['spatial_Gram']]
 def sd(a,b):
  out=ZERO
  for i,row in enumerate(B):
   for j,x in enumerate(row):
    if x:out=add(out,scale(mul(a[i],b[j]),x))
  return out
 def ld(a,b):return sub(sub(sd(a,b),mul(a[19],b[20])),mul(a[20],b[19]))
 old=[[dec(x)for x in r]for r in seed['walls']+ext['new_walls']+quarter['new_walls']];assert len(old)==587
 tr=endpoint['literal_edge_transports'][0];assert tr['omitted']==[4,23,22]and tr['type']==2
 ty=endpoint['types'][2];pr=ty['profiles'][1];assert pr['contacts']==[[12,3]]
 ids=seed['source_wall_order'];edge=[ids.index(w)for w in tr['edge_support_in_profile_order']]
 assert len(edge)==19
 G=[[ld(old[i],old[j])[0]for j in edge]for i in edge]
 assert all(ld(old[i],old[j])[1]==0 for i in edge for j in edge);piv=pd(G);I=inv(G)
 cos=[F(1,2)if w==12 else F(0)for w in tr['edge_support_in_profile_order']]
 coeff=[sum(a*b for a,b in zip(row,cos))for row in I]
 y=[q(sum(a*old[i][j][0]for a,i in zip(coeff,edge)))for j in range(19)]
 assert all(old[i][j][1]==0 for i in edge for j in range(19))
 assert all(sd(y,old[i])==q(c)for i,c in zip(edge,cos))and sd(y,y)==q(F(5,6))
 z=old[23][:19];R=F(4);center=[sub(a,scale(b,R))for a,b in zip(z,y)]
 r=center+[scale(sub(sd(center,center),q(16)),F(1,2)),ONE]
 assert ld(r,r)==q(16)
 def gate(p,den):
  if p==ZERO:return 'right'
  assert sign(p)<0
  if p==q(-den/2):return 'label3'
  if p==scale(PHI,-den/2):return 'label5'
  sg=sign(add(p,q(den)));assert sg<=0,(p,den)
  return 'parallel'if sg==0 else'dotted'
 newpairs=[];counts=Counter()
 for i,w in enumerate(old):
  p=ld(r,w)
  if i==23:assert p==q(F(11,6));continue
  den=R*(2 if i==22 else(F(3,2)if i>=527 else 1))
  label=gate(p,den);counts[label]+=1;newpairs.append({'retained_index':i,'pairing':p,'relation':label})
 assert sum(counts.values())==586
 walls=old[:];walls[23]=r
 strictcenter=[dec(x)for x in seed['common_strict_center']['point']]
 assert sign(ld(strictcenter,r))<0
 # New actual finite endpoint on the same nineteen-wall cusp edge.
 height2=q(F(8,3));x=z+[scale(add(sd(z,z),height2),F(1,2)),ONE]
 assert ld(x,x)==scale(height2,-1)and sign(ld(x,strictcenter))<0
 xp=[ld(x,w)for w in walls];zeros=[i for i,p in enumerate(xp)if p==ZERO]
 assert zeros==sorted(edge+[23])and all(sign(p)<=0 for p in xp)
 incidence=[[ld(walls[i],walls[j])for j in zeros]for i in zeros]
 assert all(v[1]==0 for row in incidence for v in row)
 ipiv=pd([[v[0]for v in row]for row in incidence])
 kept=[]
 for j,row in enumerate(seed['finite_endpoints']):
  if j==0:continue
  pnt=list(map(dec,row['point']));pp=ld(pnt,r);assert sign(pp)<0;kept.append({'old_vertex_index':j,'new_pairing':pp})
 assert len(kept)==251
 # Existing explicit points establish the other spherical facets; seed
 # facets occur at the retained251 finite vertices (checked directly).
 assert all(any(i in row['incident_wall_indices']for row in seed['finite_endpoints'][1:])for i in range(22))
 facetpoints=[seed['central_facet_point']['point']]+[row['point']for row in ext['new_actual_facet_points']]+[row['point']for row in quarter['new_actual_facet_points']]
 assert all(sign(ld(list(map(dec,p)),r))<0 for p in facetpoints)
 # The explicit annular escape is cut, but the eighth-vertex escape remains.
 target=list(map(dec,annulus['escape_point']));assert ld(target,target)==q(F(annulus['escape_norm']))
 targetold=[ld(target,w)for w in old];assert all(sign(p)<0 for p in targetold)
 targetcut=ld(target,r);assert sign(targetcut)>0
 previous=list(map(dec,quarter['strict_remaining_escape']['point']));assert ld(previous,previous)==q(F(1,2))
 remaining=[ld(previous,w)for w in walls];assert all(sign(p)<0 for p in remaining)
 paths+=['scripts/fresh_h20_intrinsic_caps_20261014_golden_seed.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return enc({'status':'PASS_EXACT_RADIUS4_REPLACEMENT_587_FACETS_WITH_STRICT_ESCAPE','sources':{p:sha(p)for p in paths},'field':seed['field'],'coordinates':seed['coordinates'],'spatial_Gram':B,'Lorentz_signature':[20,1],'retained_prior_indices':[i for i in range(587)if i!=23],'replacement_index':23,'recipe':{'endpoint_transport':0,'endpoint_type':2,'profile':1,'radius':4,'only_nonright_old_edge_contact':[12,3],'support_source_walls':tr['edge_support_in_profile_order']},'projection_direction_y':y,'projection_square':q(F(5,6)),'new_center':center,'new_normal':r,'new_norm':q(16),'removed_cap_pairing':q(F(11,6)),'new_pair_counts':dict(counts),'all586_new_pairings':newpairs,'common_strict_center':{'point':strictcenter,'norm':q(-9),'new_pairing':ld(strictcenter,r)},'new_actual_finite_endpoint':{'point':x,'norm':scale(height2,-1),'incident_indices':zeros,'all587_pairings':xp,'Gram_LDL_pivots':ipiv},'old251_finite_endpoints_retained':kept,'all_other586_facets_retained':True,'annular_escape_cut':{'point':target,'norm':ld(target,target),'pairing':targetcut},'strict_remaining_escape':{'point':previous,'norm':q(F(1,2)),'center_pairing':ld(previous,strictcenter),'all587_pairings':remaining},'scope':'One exact changed-core replacement: remove radius-one vertex cap0 (index23), insert radius-four endpoint profile1. All other586 facets and251 old finite vertices survive; a new actual finite endpoint replaces vertex0. The annular escape is cut, but the earlier eighth-vertex escape survives. No finite-volume, complete vertex census or bounded-new-cusp claim.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['new_pair_counts'])
