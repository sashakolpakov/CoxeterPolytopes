#!/usr/bin/env python3
"""Two exact parallel mates close two unbounded H22 cusp links."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
import fresh_h22_incidence_prism_20261013_verify as V
P=V.P;TAG='fresh_h22_pi7_cusp_pair_20261013';OUT=P/'results'/f'{TAG}.json'
def enc(v):return list(map(str,v))
def run(verify=False):
 M,R,h=V.setup();cert=V.read(V.OUT);geo=V.read('results/'+V.PRE+'_geometry.json');ideal=V.read('results/fresh_h22_pi7_edgewalk_20261013_ideal.json');recipes=V.read(V.REC)['recipes']
 for obj in(cert,geo,ideal):
  for p,z in obj['sources'].items():assert V.sha(p)==z
 index={tuple(r):i for i,r in enumerate(R)};perm=[index[tuple(r[:-1]+[-r[-1]])]for r in R];assert all(perm[perm[i]]==i for i in range(812))and h[-1]==0
 assert all(M[-1][j]==M[j][-1]==0 for j in range(22))and M[-1][-1]==2
 caps=[]
 for rec in recipes:
  u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=V.pair(M,u,u);qv=V.pair(M,v,v);rad=(q/qv,F(0),-q*q*qu/(4*qv));root=V.scale(V.sqrt(V.ki(rad)),rec['radical_branch']);p=V.scale(V.T,-q/2);caps.append([V.add(V.scale(p,a),V.scale(root,b))for a,b in zip(u,v)])
 added=[];records=[];mates=[4,5]
 for k in mates:
  rec=recipes[k];u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);assert q==2 and u[-1]==0 and v[-1]==1 and rec['radical_branch']==-1
  assert V.pair(M,u,u)==0 and V.pair(M,v,v)==2 and V.pair(M,u,v)==0 and caps[k][-1]==V.rat(-1)
  normal=caps[k][:-1]+[V.rat(1)];center=V.isum(V.scale(x,y)for x,y in zip(normal,V.mv(M,h)));assert center[1]<0
  cp=[]
  for j,other in enumerate(caps):
   if j==k:value=V.rat(-2);kind='parallel'
   else:
    cov=[V.isum(V.scale(x,y)for x,y in zip(other,row))for row in M];value=V.isum(V.mul(x,y)for x,y in zip(normal,cov));V.dotted(value,q,F(recipes[j]['norm']));kind='dotted'
   cp.append({'old_cap':812+j,'kind':kind,'pairing_interval':enc(value)})
  # Pairings to all812 rational walls are inherited under an exact form
  # involution which permutes the parent normals and fixes the center.
  oldpairs=[cert['new_caps'][k]['old_pair_intervals'][perm[j]]for j in range(812)];finite={str(j):cert['new_caps'][k]['finite_contacts'][str(perm[j])]for j in range(812)if str(perm[j])in cert['new_caps'][k]['finite_contacts']}
  records.append({'new_wall':835+len(added),'reflected_cap_index':k,'candidate':rec['candidate_index'],'norm':2,'center_pairing':enc(center),'rational_parent_permutation':perm,'finite_parent_contacts':finite,'parent_pairing_intervals':oldpairs,'old_cap_pairings':cp,'normal_intervals':[enc(x)for x in normal]});added.append(normal)
 # Both new normals are images under the same isometry, hence their pairing
 # equals the certified old cap4-cap5 dotted pairing.
 mutual=tuple(map(int,cert['new_caps'][5]['previous_cap_pair_intervals'][4]));V.dotted(mutual,2,2)
 links=[]
 for row,normal,newidx in zip(ideal['rows'],added,(835,836)):
  cap=row['active_new_caps'][0];assert cap==812+mates[newidx-835];null=list(map(F,row['primitive_integer_null_ray']));assert null[-1]==0 and V.pair(M,null,null)==0
  active=row['active_base_walls'];assert all(R[i][-1]==0 for i in active)
  # New mate has exactly the same pairing as its original at any point with
  # last coordinate0. The old audit proves all remaining walls strict there.
  assert row['all_cap_pairings'][cap-812]['tau_coefficient']==row['all_cap_pairings'][cap-812]['radical_coefficient']=='0'
  other=mates[1-(newidx-835)];assert int(row['all_cap_pairings'][other]['pairing_interval'][1])<0
  assert all(cert['new_caps'][cap-812]['finite_contacts'][str(i)]==2 for i in active)
  # Active Gram is the orthogonal sum of the rank8/rank12 affine blocks
  # and [[2,-2],[-2,2]]. Each has a strictly positive null vector.
  assert sum(c['rank']for c in row['affine_components'])==20
  assert all(min(map(F,c['positive_null_coefficients']))>0 for c in row['affine_components'])
  links.append({'vertex':row['vertex'],'released':row['released'],'primitive_integer_null_ray':row['primitive_integer_null_ray'],'active_walls':active+[cap,newidx],'affine_component_ranks':[c['rank']for c in row['affine_components']]+[1],'positive_kernel_blocks':row['affine_components']+[{'walls':[cap,newidx],'rank':1,'positive_null_coefficients_order':[cap,newidx],'positive_null_coefficients':['1','1']}],'complete_active_Gram_inertia':[21,0,3],'link_dimension':21,'bounded_link':True,'proof':'Every active component has a positive kernel vector. Their quotient normals span Euclidean dimension21. A recession vector pairs nonpositively with every normal, so the positive relations force all pairings zero; spanning then forces that vector to vanish. The ideal endpoint is unchanged and every nonactive wall remains strict.'})
 # Reuse the previous rational strict escape only if it passes both additions.
 esc=geo['strict_escape'];x=[F(n,esc['denominator'])for n in esc['numerator']];assert V.pair(M,x,x)>0 and V.pair(M,x,h)<0;cov=V.mv(M,x);newpairs=[V.isum(V.scale(a,b)for a,b in zip(r,cov))for r in added];assert all(p[1]<0 for p in newpairs),'PRIOR_ESCAPE_CUT: new volume proposal required'
 out={'status':'PASS_EXACT_837_FACETS_AND_TWO_BOUNDED_CUSP_LINKS','sources':{s:V.sha(s)for s in(V.OUT,V.REC,'results/'+V.PRE+'_geometry.json','results/fresh_h22_pi7_edgewalk_20261013_ideal.json','scripts/'+Path(__file__).name,'scripts/'+V.PRE+'_verify.py')},'facets':837,'new_caps':records,'new_mutual_pairing':enc(mutual),'bounded_cusp_links':links,'strict_escape':{'numerator':esc['numerator'],'denominator':esc['denominator'],'norm':str(V.pair(M,x,x)),'center_pairing':str(V.pair(M,x,h)),'new_cap_pairings':[enc(z)for z in newpairs]},'scope':'Two exact last-coordinate reflections preserve all835 old walls and create parallel mates at the two previously unbounded ideal links. Both links are now bounded in all21 Euclidean directions. The same strict spacelike witness remains feasible for all837 walls, so the full chamber still has infinite volume.'}
 if verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'escape',out['strict_escape']['norm'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
