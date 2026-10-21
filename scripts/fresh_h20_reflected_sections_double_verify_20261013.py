#!/usr/bin/env python3
"""Independent exact replay of changed sections of an ambient two-chamber union."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
import json
import fresh_h20_joint_cusp_irrational_verify_20261013 as prior
A=prior.A;B=prior.B;K=prior.K;V=prior.V;U=prior.U;ROOT=prior.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_reflected_sections_double_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_reflected_sections_double_verify_20261013.json'

def boundary(Q,W,G,wall):
 rows=[];lookup={};q=G[wall][wall]
 for tile in (0,1):
  for j,w in enumerate(W):
   if j==wall:continue
   coefficient=F(-2*G[wall][j],q);assert coefficient>=0
   root=list(w) if tile==0 else [x+coefficient*y for x,y in zip(w,W[wall])]
   assert all(x.denominator==1 for x in map(F,root));key=tuple(map(int,root))
   if key not in lookup:lookup[key]=len(rows);rows.append({'normal':list(key),'owners':[]})
   rows[lookup[key]]['owners'].append([tile,j+1])
 return rows

def cusp(row,case,Q,W,G,e,field):
 u=list(map(F,case['u']));ids=[i-1 for i in case['component']];drop=case['dropped']-1
 co=[A.readp(p) for p in row['component_coefficients']];assert len(co)==len(ids)
 assert e==[field.rem(V.sum_polys(A.scale(c,W[i][d]) for c,i in zip(co,ids))) for d in range(22)]
 sg=field.sign(co[ids.index(drop)]);assert sg
 assert all(field.sign(c)*sg<0 for i,c in zip(ids,co) if i!=drop)
 _,q,v,gr=prior.products(e,Q,W,G,field);active=[i for i,w in enumerate(W) if A.dot(w,u,Q)==0 and i!=drop]
 H=[[gr(i,j) for j in active] for i in active]
 assert all(field.sign(H[i][i])>0 for i in range(len(active)))
 for i,j in combinations(range(len(active)),2):
  assert field.sign(H[i][j])<=0
  N=field.rem(A.mul(H[i][j],H[i][j]));den=field.rem(A.mul(H[i][i],H[j][j]))
  assert any(field.rem(A.sub(N,A.scale(den,z)))==A.Z for z in (0,F(1,4),F(1,2),F(3,4),1))
 comps=K.groups([[int(p!=A.Z) for p in r] for r in H],range(len(active)));proofs=row['cusp_link_in_original_wall_labels'];assert len(comps)==len(proofs)
 ranks=0
 for ks,proof in zip(comps,proofs):
  walls=[active[k] for k in ks];assert proof['walls']==[i+1 for i in walls];marks=[A.readp(p) for p in proof['marks']]
  assert len(marks)==len(ks) and all(field.sign(c)>0 for c in marks)
  for k in ks:assert field.rem(V.sum_polys(A.mul(c,H[k][l]) for c,l in zip(marks,ks)))==A.Z
  assert proof['rank']==len(ks)-1;ranks+=proof['rank'];mu=A.readp(proof['cusp_multiple']);assert field.sign(mu)>0
  for d in range(22):
   null=field.rem(V.sum_polys(A.mul(c,A.sub(A.scale(q,W[i][d]),A.mul(v[i],e[d]))) for c,i in zip(marks,walls)))
   assert null==field.rem(A.scale(mu,u[d]))
 assert ranks==19

def verify_row(row,case,Q,W,G,R,H,h,field,wall,control=False):
 base=[A.readp(p) for p in case['base_normal']];u=list(map(F,case['u']));tau=A.readp(row['translation'])
 e=[field.rem(A.add(p,A.scale(tau,z))) for p,z in zip(base,u)];assert e==[A.readp(p) for p in row['normal']]
 sd=prior.products(e,Q,R,H,field);_,q,v,gr=sd;assert field.sign(q)>0 and q==A.readp(row['normal_square'])
 center=[A.readp(p) for p in row['strict_center']];dual,cp,square=U.field_products(center,Q,R,field)
 assert field.sign(square)<0 and all(field.sign(p)<0 for p in cp)
 assert field.sign(V.sum_polys(A.scale(dual[d],h[d]) for d in range(22)))<0
 assert field.rem(V.sum_polys(A.mul(dual[d],e[d]) for d in range(22)))==A.Z
 cusp(row,case,Q,W,G,e,field)
 if not control:
  a=row['orthogonal_new_boundary_index']-1;assert R[a] not in W
  base_sd=prior.products(base,Q,R,H,field);up=A.dot(R[a],u,Q);assert up!=0
  assert field.rem(A.add(base_sd[2][a],A.scale(tau,up)))==A.Z and v[a]==A.Z
  assert tuple(row['translation']) not in {tuple(old['translation']) for old in case['candidates']}
  pairing=field.rem(V.sum_polys(A.scale(e[d],sum(Q[d][j]*W[wall][j] for j in range(22))) for d in range(22)))
  assert pairing==A.readp(row['section_not_fixed_by_reflection_pairing']) and pairing!=A.Z
  target=row['target_ridge'];assert target is not None and target['ridge_walls_one_based']==row['target_pair'] and a+1 in row['target_pair']
  aa,ab,bb=B.actual(target,Q,R,H,h,sd,field);assert ab==A.Z
 if row['doubled_old_ridge_right_angle'] is not None:
  aa,ab,bb=B.actual(row['doubled_old_ridge_right_angle'],Q,R,H,h,sd,field);assert ab==A.Z
 elif control:raise AssertionError('control failed to repair original ridge')
 proof=row['actual_forbidden_ridge'];aa,ab,bb=B.actual(proof,Q,R,H,h,sd,field)
 N=field.rem(A.mul(ab,ab));den=field.rem(A.mul(aa,bb))
 if proof['angle_rejection']=='obtuse':assert field.sign(ab)>0
 else:
  assert field.sign(ab)<0 and field.sign(N)>0 and field.sign(A.sub(A.scale(den,3),A.scale(N,4)))>0
  for p in (A.sub(A.scale(N,4),den),A.sub(A.scale(N,2),den),A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den))):assert field.sign(p)!=0
 assert row['status']=='NON_COXETER_ACTUAL_RIDGE'

def main():
 Q,W,G,h=A.source();data=json.loads(SOURCE.read_text());source=json.loads((ROOT/data['source']).read_text());case=source['cases'][data['source_case']]
 wall=data['doubling_wall_one_based']-1;assert wall==168
 expected=boundary(Q,W,G,wall);assert expected==data['ambient_boundary_normals'];R=[row['normal'] for row in expected];assert len(R)==data['ambient_boundary_count']==388
 H=[[A.dot(r,s,Q) for s in R] for r in R]
 u=list(map(F,case['u']));assert A.dot(u,u,Q)==0 and A.dot(u,h,Q)<0 and A.dot(W[wall],u,Q)<0
 assert {tuple(r) for r in R if A.dot(r,u,Q)==0}=={tuple(w) for w in W if A.dot(w,u,Q)==0}
 assert all(A.dot(r,u,Q)<=0 for r in R)
 assert u==[sum(W[i][d] for i in (33,45,203,209)) for d in range(22)]
 expected_base=[A.trim((F(W[45][d],2)-u[d]/36,F(W[33][d]-W[45][d]-W[203][d]+W[209][d],4))) for d in range(22)]
 assert expected_base==[A.readp(p) for p in case['base_normal']]
 field=K.Field(A.readp(data['field_polynomial']),*data['field_interval']);assert field.p==A.trim((-2,0,1)) and field.sign(A.T)>0
 assert len(data['candidates'])==6
 for row in data['candidates']:verify_row(row,case,Q,W,G,R,H,h,field,wall)
 control=data['control'];oldrow=case['candidates'][control['source_candidate']];assert control['normal']==oldrow['normal'];verify_row(control,case,Q,W,G,R,H,h,field,wall,True)
 oldsd=prior.products([A.readp(p) for p in control['normal']],Q,W,G,field);oldproof=oldrow['actual_forbidden_ridge'];assert oldproof['ridge_walls_one_based']==[1,169]
 aa,ab,bb=B.actual(oldproof,Q,W,G,h,oldsd,field);assert field.rem(A.sub(A.scale(A.mul(ab,ab),52),A.scale(A.mul(aa,bb),27)))==A.Z
 newids=control['doubled_old_ridge_right_angle']['ridge_walls_one_based'];assert [0,1] in expected[newids[0]-1]['owners'] and [1,1] in expected[newids[1]-1]['owners']
 forged=deepcopy(data['candidates'][0]);forged['target_ridge']['point_coordinate_polynomials'][0][-1]=str(F(forged['target_ridge']['point_coordinate_polynomials'][0][-1])+1)
 try:verify_row(forged,case,Q,W,G,R,H,h,field,wall)
 except AssertionError:pass
 else:raise AssertionError('forged reflected-mirror ridge accepted')
 result={'status':'PASS','standard_library_only':True,'ambient_chambers':2,'ambient_boundary_planes':388,'new_section_planes':6,'actual_new_mirror_right_ridges':6,
  'actual_forbidden_ridge_rejections':6,'control_old_forbidden_ridge_repaired':True,'control_still_non_Coxeter':True,'whole_cusp_rank':19,'forged_reflected_ridge_rejected':True,
  'scope':data['scope']}
 OUTPUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
