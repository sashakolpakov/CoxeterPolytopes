#!/usr/bin/env python3
"""Independent SymPy audit of all saved local Grams and global ridge proofs.

No producer arithmetic is imported. Graph-cut completeness is reviewed
mathematically, rather than enumerated a second time here.
"""
from pathlib import Path
import argparse, hashlib, json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_single_factor_classification_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT=ROOT/f'results/{TAG}_root_review.json'
def mat(a):return s.Matrix(a).applyfunc(s.Rational)
def run():
 z=json.loads((ROOT/f'results/{TAG}.json').read_text())
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=-mat(raw['basis_gram']);W=mat(raw['wall_coordinates_as_columns_transposed']);C=W*Q;G=C*W.T
 h=s.Matrix([sum(W[i,j]for i in range(210))for j in range(22)]);h=h/s.gcd_list(list(h))
 assert (h.T*Q*h)[0]==-14 and max(C*h)<0
 allowed={0,s.Rational(1,4),s.Rational(1,2),s.Rational(3,4),1}
 pointspec={};count=0;t=s.symbols('t')
 for row in z['locally_Coxeter_planes']:
  e=mat(row['normal_constant']);D=mat(row['null_shift']);q=(e.T*Q*e)[0]
  assert q==s.Rational(row['normal_norm'])>0 and (e.T*Q*D)[0]==(D.T*Q*D)[0]==0
  ids=row['factor'];p=(C*e).extract(ids,[0]);B=G.extract(ids,ids)-p*p.T/q
  assert B==mat(row['projected_local_Gram']) and min(B.diagonal())>0
  for i in range(len(ids)):
   for j in range(i):assert B[i,j]<=0 and B[i,j]**2/(B[i,i]*B[j,j])in allowed
  proofs=row['global_actual_obtuse_ridge_certificates'];assert len(proofs)in[1,2]
  if len(proofs)==2:assert {tuple(p['active_part'])for p in proofs}=={tuple(row['S']),tuple(row['T'])}
  for proof in proofs:
   count+=1;a=mat(proof['normal_constant']);q=(a.T*Q*a)[0]
   L=mat(proof['normalized_left']);R=mat(proof['normalized_right'])
   upper=s.Rational(proof['parameter_interval'][1]);i,j=proof['ridge']
   assert proof['parameter_interval'][0]==-1 and upper==(0 if len(proofs)==1 else s.Rational(-1,2))
   assert proof['upper_endpoint_included']==(len(proofs)==2)
   assert (a.T*Q*D)[0]==0 and (D.T*Q*L)[0]==(D.T*Q*R)[0]==-1
   assert (a.T*Q*L)[0]==-1 and (a.T*Q*R)[0]==upper
   x=-t*L+(1+t)*R if upper==0 else(-2*t-1)*L+2*(1+t)*R
   assert s.simplify(((a+t*D).T*Q*x)[0])==0 and (L.T*Q*R)[0]<0
   lz=tuple(proof['left_ambient_face']['walls'])
   rz=tuple(proof['right_ambient_face_or_midpoint']['walls'])if upper==0 else tuple(sorted([i,j]))
   assert set(lz)&set(rz)=={i,j}
   for P,zeros in[(L,lz),(R,rz)]:
    key=tuple(P)
    if key in pointspec:assert pointspec[key]==zeros
    pointspec[key]=zeros
   ap=C*a;dp=C*D;ai=ap[i];aj=ap[j]
   assert dp[i]==0 and aj==dp[j]<0 and G[i,j]==0 and ai>0
   ni=G[i,i]-ai*ai/q;coefficient=-ai*aj/q;bound=ni*G[j,j]-aj*aj*G[i,i]*(1+upper)**2/q
   assert coefficient==s.Rational(proof['projected_pairing_equals_coefficient_times_1_plus_t'])>0
   assert bound==s.Rational(proof['ridge_determinant_lower_bound'])>0 and ni>0
 P=s.Matrix.hstack(*(s.Matrix(p)for p in pointspec));signs=C*P;norms=Q*P;future=h.T*Q*P
 for k,(point,zeros)in enumerate(pointspec.items()):
  assert max(signs[:,k])==0 and tuple(i for i in range(210)if signs[i,k]==0)==tuple(sorted(zeros))
  assert sum(point[i]*norms[i,k]for i in range(22))<0 and future[k]<0
 for row in z['locally_rejected_partitions']:
  e=mat(row['normal_constant']);q=(e.T*Q*e)[0];i,j=row['bad_pair'];p=C*e
  value=(G[i,j]-p[i]*p[j]/q)**2/((G[i,i]-p[i]**2/q)*(G[j,j]-p[j]**2/q))
  assert value==s.Rational(row['rational_squared_cosine'])and value not in allowed
 files=[RAW,f'results/{TAG}.json',f'scripts/{TAG}.py',f'scripts/{TAG}_root_review.py']
 return {'status':'PASS_INDEPENDENT_SYMPY_LOCAL_GRAMS_AND_GLOBAL_RIDGE_CERTIFICATES',
  'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in files},
  'local_acceptances':len(z['locally_Coxeter_planes']),'local_rejections':len(z['locally_rejected_partitions']),
  'continuous_ridge_proofs':count,'distinct_exact_endpoint_points':len(pointspec),
  'scope':'Independent saved-certificate arithmetic; star-cut exhaustion and Euclidean classification are mathematical/source review.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run()
 if args.verify:assert z==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status']);print({k:v for k,v in z.items()if k not in['sources','scope']})
