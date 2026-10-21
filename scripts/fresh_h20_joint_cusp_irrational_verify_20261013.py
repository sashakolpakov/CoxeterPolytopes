#!/usr/bin/env python3
"""Independent stdlib/Fraction replay of twelve irrational cusp-preserving cuts."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
import json
import fresh_h20_new_route_20261013_verify as K
A=K.A;V=K.V;U=K.U;B=K.B;ROOT=K.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_joint_cusp_irrational_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_joint_cusp_irrational_verify_20261013.json'

def products(e,Q,W,G,field):
 dual,vals,q=U.field_products(e,Q,W,field)
 def gr(a,b):return field.rem(A.sub(A.scale(q,G[a][b]),A.mul(vals[a],vals[b])))
 return e,q,vals,gr

def verify_candidate(row,case,Q,W,G,h,field):
 u=list(map(F,case['u']));base=[A.readp(p) for p in case['base_normal']];tau=A.readp(row['translation'])
 e=[field.rem(A.add(p,A.scale(tau,z))) for p,z in zip(base,u)]
 assert e==[A.readp(p) for p in row['normal']]
 # Irrational sqrt(2) forces both rational coefficient vectors into any
 # rational wall span containing e. Exact quotient-ray grouping tests all
 # spans of at most three original walls, including zero/repeated rays.
 E=[[p[0] for p in e],[p[1] if len(p)>1 else F(0) for p in e]]
 assert len(K.C.rref(E+[u])[1])==3
 assert K.C.plane_small_support(E,W) is None
 e,q,vals,gr=products(e,Q,W,G,field);assert field.sign(q)>0
 assert q==A.readp(row['normal_square'])
 assert field.rem(V.sum_polys(A.scale(e[d],sum(Q[d][j]*u[j] for j in range(22))) for d in range(22)))==A.Z
 base_e,base_q,bv,bgr=products(base,Q,W,G,field)
 assert q==base_q
 a,b=[x-1 for x in row['target_pair']];up=[A.dot(w,u,Q) for w in W]
 assert up[a]*up[b]==0 and (up[a]!=0 or up[b]!=0)
 slope=field.rem(A.add(A.scale(bv[a],up[b]),A.scale(bv[b],up[a])))
 assert slope!=A.Z and field.rem(A.sub(bgr(a,b),A.mul(tau,slope)))==A.Z
 assert gr(a,b)==A.Z # Unique exact solution of the declared linear target.
 coefficients=[A.readp(p) for p in row['component_coefficients']];ids=[x-1 for x in case['component']]
 assert len(coefficients)==len(ids)
 assert e==[field.rem(V.sum_polys(A.scale(c,W[i][d]) for c,i in zip(coefficients,ids))) for d in range(22)]
 drop=case['dropped']-1;di=ids.index(drop);sg=field.sign(coefficients[di]);assert sg
 assert all(field.sign(c)*sg<0 for k,c in enumerate(coefficients) if k!=di)
 # Hence the omitted wall is a strictly positive combination of the
 # retained component walls after restricting to e-perp.
 center=[A.readp(p) for p in row['strict_center']];dual,cp,square=U.field_products(center,Q,W,field)
 assert field.sign(square)<0 and all(field.sign(p)<0 for p in cp)
 assert field.sign(V.sum_polys(A.scale(dual[d],h[d]) for d in range(22)))<0
 assert field.rem(V.sum_polys(A.mul(dual[d],e[d]) for d in range(22)))==A.Z
 active=[i for i,p in enumerate(up) if p==0 and i!=drop]
 H=[[gr(i,j) for j in active] for i in active]
 assert all(field.sign(H[i][i])>0 for i in range(len(active)))
 for i,j in combinations(range(len(active)),2):
  assert field.sign(H[i][j])<=0
  N=field.rem(A.mul(H[i][j],H[i][j]));den=field.rem(A.mul(H[i][i],H[j][j]))
  assert any(field.rem(A.sub(N,A.scale(den,z)))==A.Z for z in (0,F(1,4),F(1,2),F(3,4),1))
 adjacency=[[int(x!=A.Z) for x in r] for r in H];components=K.groups(adjacency,range(len(active)))
 assert len(components)==len(row['cusp_link']);rank=0
 for ks,proof in zip(components,row['cusp_link']):
  walls=[active[k] for k in ks];assert proof['walls']==[i+1 for i in walls]
  marks=[A.readp(p) for p in proof['marks']];assert len(marks)==len(ks) and all(field.sign(c)>0 for c in marks)
  for k in ks:assert field.rem(V.sum_polys(A.mul(c,H[k][l]) for c,l in zip(marks,ks)))==A.Z
  # A connected Gram with nonpositive off-diagonals and positive null
  # vector is congruent to a weighted graph Laplacian: corank exactly one.
  assert proof['rank']==len(ks)-1;rank+=proof['rank']
  mu=A.readp(proof['cusp_multiple']);assert field.sign(mu)>0
  for d in range(22):
   null=field.rem(V.sum_polys(A.mul(c,A.sub(A.scale(q,W[i][d]),A.mul(vals[i],e[d]))) for c,i in zip(marks,walls)))
   assert null==field.rem(A.scale(mu,u[d]))
 assert rank==19
 sd=(e,q,vals,gr)
 if row['repaired_target_ridge'] is not None:
  assert row['repaired_target_ridge']['ridge_walls_one_based']==row['target_pair']
  aa,ab,bb=B.actual(row['repaired_target_ridge'],Q,W,G,h,sd,field);assert ab==A.Z
 for proof in row['additional_right_ridges']:
  aa,ab,bb=B.actual(proof,Q,W,G,h,sd,field);assert ab==A.Z
 assert row['status']=='NON_COXETER_ACTUAL_RIDGE'
 proof=row['actual_forbidden_ridge'];aa,ab,bb=B.actual(proof,Q,W,G,h,sd,field)
 N=field.rem(A.mul(ab,ab));den=field.rem(A.mul(aa,bb))
 if proof['angle_rejection']=='obtuse':assert field.sign(ab)>0
 else:
  assert field.sign(ab)<0 and field.sign(N)>0 and field.sign(A.sub(A.scale(den,3),A.scale(N,4)))>0
  for p in (A.sub(A.scale(N,4),den),A.sub(A.scale(N,2),den),A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den))):assert field.sign(p)!=0
 return row['repaired_target_ridge'] is not None

def main():
 Q,W,G,h=A.source();data=json.loads(SOURCE.read_text());count=targets=0;example=None
 cusp_source=json.loads((ROOT/data['cusp_source']).read_text())
 for case in data['cases']:
  assert case['cusp_index']==4 and case['u']==cusp_source['cusps'][4]['vector'];u=list(map(F,case['u']))
  assert A.dot(u,u,Q)==0 and A.dot(u,h,Q)<0 and all(A.dot(w,u,Q)<=0 for w in W)
  assert case['component']==[34,46,204,210]
  field=K.Field(A.readp(case['field_polynomial']),*case['field_interval']);assert field.p==A.trim((-2,0,1)) and field.sign(A.T)>0
  base=[A.readp(p) for p in case['base_normal']];sign=case['source_direction_sign'];assert sign in(-1,1)
  assert base==[A.trim((F(a),sign*F(b))) for a,b in zip(case['source_e0'],case['source_e1'])]
  assert field.rem(V.sum_polys(A.scale(p,sum(Q[d][j]*h[j] for j in range(22))) for d,p in enumerate(base)))==A.Z
  for row in case['candidates']:
   targets+=verify_candidate(row,case,Q,W,G,h,field);count+=1;example=(row,case,field)
 assert len(data['cases'])==2 and count==12 and targets==8
 row,case,field=example;forged=deepcopy(row);forged['component_coefficients'][0][-1]=str(F(forged['component_coefficients'][0][-1])+1)
 try:verify_candidate(forged,case,Q,W,G,h,field)
 except AssertionError:pass
 else:raise AssertionError('forged redundancy relation accepted')
 result={'status':'PASS','standard_library_only':True,'irrational_cusp_quotient_directions':2,'finite_volume_noncompact_sections':count,
  'minimum_original_wall_support_exactly_four':count,'actual_targeted_right_ridges':targets,'actual_forbidden_ridge_rejections':count,'whole_cusp_rank':19,'forged_redundancy_relation_rejected':True,
  'scope':'The twelve literal linear global-angle equations in two irrational quotient directions; no completeness over translations or local direction scouting.'}
 OUTPUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
