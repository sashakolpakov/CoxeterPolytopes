#!/usr/bin/env python3
"""All-residual-vector contact restrictions for the frozen79 H20 chamber.

Rational nonnegative wall multipliers propose exact Lorentz-Cauchy bounds.
No optimizer is imported or trusted by this certificate.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_contact_closure_20261014_replacement_review as K
import fresh_h12_a12_source_20261014 as Qlin

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_completion_gate_20261014'
BASE='results/fresh_h20_e8_cusp_completion_20261014.json'
PRE='results/fresh_h20_contact_closure_20261014_h12_replacement.json'
SOURCE='results/fresh_h12_a12_source_20261014.json'
RECIPE=f'results/{TAG}_recipes.json'
Z=K.ZERO
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def vec(r):return[tuple(map(F,x))for x in r]
def run():
 d,pre,source,recipes=[json.loads((ROOT/p).read_text())for p in[BASE,PRE,SOURCE,RECIPE]]
 for doc in[d,pre,source]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 assert sha(f'results/{TAG}_scout.json')==recipes['scout_sha256']
 Q=[[F(x)for x in row]for row in d['form']];N=list(map(vec,d['normals']));E=vec(d['strict_escape']['point']);H=vec(d['strict_center']);J=d['E8_support'];assert J==[12,13,14,15,18,26,33,67]
 cache={}
 def dot(v,w):
  key=tuple(w)
  if key not in cache:cache[key]=[tuple(sum(Q[i][j]*w[j][k]for j in range(21))for k in range(4))for i in range(21)]
  out=Z
  for a,b in zip(v,cache[key]):out=K.add(out,K.mul(a,b))
  return out
 def linear(vs,cs):
  out=[Z]*21
  for v,c in zip(vs,cs):out=[K.add(a,K.mul(b,c))for a,b in zip(out,v)]
  return out
 B=[[dot(N[i],N[j])for j in J]for i in J];assert all(x[1:]==(0,0,0)for r in B for x in r)
 B=[[x[0]for x in r]for r in B];Qlin.positive_definite(B)
 I=list(map(list,zip(*[Qlin.solve(B,[int(i==j)for i in range(8)])for j in range(8)])))
 assert I==[[F(x)for x in r]for r in pre['E8_inverse']]
 weights=[linear([N[j]for j in J],list(map(K.scalar,row)))for row in I]
 assert all(dot(w,N[j])==K.scalar(int(i==k))for i,w in enumerate(weights)for k,j in enumerate(J))
 def project(v):
  ew=linear(weights,[dot(v,N[j])for j in J])
  return[K.add(a,K.neg(b))for a,b in zip(v,ew)]
 EH=project(E);HP=project(H)
 assert K.sign(dot(EH,EH))<0 and K.sign(dot(EH,HP))<0
 assert K.sign(dot(HP,HP))<0
 # Identify the complete intrinsic source wall system by exact Gram data.
 mapping={int(k):v for k,v in pre['source_wall_map'].items()}
 assert set(mapping)==set(range(14))
 assert all(dot(N[i],N[j])==Z for i in mapping.values()for j in J)
 assert all(dot(N[mapping[i]],N[mapping[j]])==K.scalar(F(source['gram'][i][j]))for i in range(14)for j in range(14))
 assert all(K.sign(dot(EH,N[i]))<0 for i in mapping.values())
 assert source['projective_active_sets']==91 and source['finite_vertex_count']==35 and source['ideal_vertex_count']==2
 assert all(F(row['norm'])<=0 for row in source['projective_vertices'])
 rows=[]
 for recipe in recipes['recipes']:
  support=recipe['support'];lambdas=recipe['nonnegative_wall_multipliers']
  assert len(set(i for i,_ in lambdas))==len(lambdas)
  X=E[:]
  for i,value in lambdas:
   t=F(value);assert t>0 and 0<=i<79
   X=[K.add(a,K.scale(b,-t))for a,b in zip(X,N[i])]
  Y=project(X);kappa=K.neg(dot(Y,Y));assert K.sign(kappa)>0 and K.sign(dot(Y,HP))<0
  a=[K.neg(dot(w,X))for w in weights]
  M=[[K.add(K.scale(kappa,I[i][j]),K.neg(K.mul(a[i],a[j])))for j in range(8)]for i in range(8)]
  kind=recipe['kind'];margin=None
  if kind=='single_negative_linear':
   assert len(support)==1 and K.sign(a[support[0]])<0
   margin=K.neg(a[support[0]])
  elif kind=='single_threshold':
   assert len(support)==1;i=support[0];lower=F(recipe['minimum_squared_magnitude'])
   assert lower in[1,2]and K.sign(M[i][i])>0
   margin=K.add(K.scale(M[i][i],lower),K.scale(kappa,-2));assert K.sign(margin)>0
  elif kind=='independent_magnitudes_at_least_one':
   assert len(support)==2
   assert all(K.sign(M[i][j])>0 for i in support for j in support)
   margin=K.scale(kappa,-2)
   for i in support:
    for j in support:margin=K.add(margin,M[i][j])
   assert K.sign(margin)>0
  else:
   assert kind=='equal_magnitudes_at_least_one'and support==[0,7]
   total=Z
   for i in support:
    for j in support:total=K.add(total,M[i][j])
   assert K.sign(total)>0
   margin=K.add(total,K.scale(kappa,-2));assert K.sign(margin)>0
  rows.append(dict(kind=kind,support=support,minimum_squared_magnitude=recipe.get('minimum_squared_magnitude'),
                   nonnegative_wall_multipliers=lambdas,adjusted_escape=X,future_timelike_projection=Y,
                   kappa=kappa,linear_coefficients=a,quadratic_matrix=M,strict_algebraic_margin=margin))
 assert len(rows)==11
 assert{tuple(r['support'])for r in rows if r['kind']=='independent_magnitudes_at_least_one'}=={(0,1),(0,2),(1,2),(1,3),(1,7),(2,3)}
 # The minimal E8 weight has norm2. At its pi3 magnitude the residual is
 # null, hence lies on one of the only two source null extreme rays.
 rays=list(map(vec,d['cusp_rays']));minus_weight=[K.neg(x)for x in weights[0]]
 assert dot(minus_weight,minus_weight)==K.scalar(2)
 C0=dot(minus_weight,E);assert K.sign(C0)>0
 ray_rows=[]
 for nu in rays:
  assert dot(nu,nu)==Z and K.sign(dot(nu,HP))<0 and all(dot(nu,N[j])==Z for j in J)
  slope=dot(nu,E);assert K.sign(slope)<0
  gates=[]
  for i,r in enumerate(N):
   denominator=K.neg(dot(nu,r));offset=dot(minus_weight,r)
   if K.sign(denominator)<=0 or K.sign(offset)<=0:continue
   numerator=K.add(K.mul(denominator,C0),K.mul(offset,slope))
   if K.sign(numerator)<0:gates.append(dict(wall=i,positive_denominator=denominator,positive_offset=offset,strict_upper_numerator=numerator))
  assert gates
  ray_rows.append(dict(point=nu,escape_slope=slope,witness=gates[0]))
 assert len(ray_rows)==2
 # Verify that these two rays are the frozen source's complete null set.
 assert all(vec(row['point'])==nu for row,nu in zip(pre['H12_face_ideal_vertices'],rays))
 paths=[BASE,PRE,SOURCE,RECIPE,f'results/{TAG}_scout.json',f'scripts/{TAG}_scout.py',f'scripts/{TAG}.py',
        'scripts/fresh_h20_contact_closure_20261014_replacement_review.py','scripts/fresh_h12_a12_source_20261014.py']
 return Qlin.encode(dict(status='PASS_ALL_RESIDUAL_VECTOR_E8_CONTACT_RESTRICTIONS',sources={p:sha(p)for p in paths},
   E8_positions_to_walls=J,E8_gram=B,E8_inverse=I,escape=E,escape_E8_orthogonal_projection=EH,
   projected_escape_norm=dot(EH,EH),E8_fundamental_weights=weights,source_wall_map={str(k):v for k,v in mapping.items()},
   dual_certificates=rows,minimal_weight_pi3_null_ray_bounds=ray_rows,
   excluded_all_label_one_node_positions=[0,1,7],one_node_position2_requires_label3=True,
   excluded_arbitrary_magnitude_two_node_supports=[[0,1],[0,2],[1,2],[1,3],[1,7],[2,3]],
   excluded_equal_magnitude_two_node_support=[0,7],
   formula='For r=-w+y, w=sum c_i omega_i, y²=2-c^T I c and y in the future H12 cone. Nonnegative wall multipliers give r.E<=a^T c-sqrt(kappa*(c^T I c-2)). Each certificate proves this is strictlynegative on its stated whole contact domain.',
   scope='All real residual y and all positive-norm caps after normalization to norm2; retains all79 walls and their acute Coxeter pairings. Listed support restrictions only. No blanket obstruction, no complete contact-profile classification, and no assertion that surviving profiles give Coxeter caps. Numerical SOCP scouts are proposals only; exact acceptance uses their frozen rational multipliers.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['dual_certificates']))
