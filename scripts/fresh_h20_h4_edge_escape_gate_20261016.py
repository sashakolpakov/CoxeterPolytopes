#!/usr/bin/env python3
"""All-parameter E4 obstruction for the four orthogonal residual edge axes."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from fresh_h20_s61_multicap_20261014_exact import I
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_h4_edge_escape_gate_20261016'
CLUSTER='results/fresh_h20_residual_h4_corner_20261016.json';MIXED='results/fresh_h20_mixed_isolated_cap_20261016.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def saved(x):
 d=1<<192
 return[str(F((x.l*d).__floor__(),d)),str(F((x.h*d).__ceil__(),d))]
def fd(v,w):return sum(g*a*b for g,a,b in zip([F(-1),F(1),F(1),F(1,2),F(1,2)],v,w))
def run():
 c=json.loads((ROOT/CLUSTER).read_text());m=json.loads((ROOT/MIXED).read_text())
 for doc in [c,m]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 assert c['status']=='PASS_EXACT_145_WALL_RESIDUAL_H4_FINITE_CORNER_CLUSTER'
 R=[list(map(F,row))for row in c['definition']['four_recipes']];q=list(map(F,c['definition']['four_norms']))
 C=list(map(F,c['new_finite_corner']['chart']));Cq=F(c['new_finite_corner']['norm']);assert fd(C,C)==Cq<0
 e=list(map(F,c['strict_escape']['frame_coefficients']))+[F(-3,10),F(-1,20)]
 ce=fd(C,e);assert ce<0
 L=[fd(r,e)/ce for r in R];assert all(x>0 for x in L)
 def prod(i,j):
  if i==141:
   row=m['pairs'][j];assert row['wall']==j
  else:row=next(row for row in c['pairs']if row['walls']==[j,i])
  return I(*row['pair_bounds'])
 def cp(j):return I(*c['new_finite_corner']['all145_pair_bounds'][j])
 ordinary=[]
 for i in range(3):
  a=prod(141+i,133);b=cp(133);assert a.h<0 and b.h<0
  lower=a/b;assert lower.l>L[i]
  ordinary.append(dict(axis=141+i,blocking_wall=133,axis_pair_bounds=saved(a),corner_pair_bounds=saved(b),
      acute_lambda_lower_bounds=saved(lower),E4_cut_lambda_upper=str(L[i]),strict_incompatibility_bounds=saved(lower-L[i])))
 # The w axis leaves only a forbidden-angle interval at chain132.
 a=prod(144,132);b=cp(132);assert a.h<0 and b.h<0
 A=a/b;U=L[3];assert 0<A.l<=A.h<U
 qU=q[3]+Cq*U*U;assert qU>0
 gU=-a+U*b;assert gU.h<0
 uppergap=I(2*qU)-4*gU*gU;assert uppergap.l>0
 # lambda=A is right to132, but then the isolated root31 is forbidden.
 assert prod(144,31).l<=F(-3,2)<=prod(144,31).h
 assert cp(31).l<=C[4]<=cp(31).h
 g31=I(F(3,2))+A*C[4];qA=I(q[3])+Cq*A*A
 assert qA.l>0 and g31.h<0
 endpointgap=2*qA-4*g31*g31;assert endpointgap.l>0
 tail=dict(axis=144,blocking_wall=132,acute_lambda_lower_bounds=saved(A),E4_cut_lambda_upper=str(U),
     axis_pair_bounds=saved(a),corner_pair_bounds=saved(b),norm_at_cut_upper=str(qU),
     blocking_pair_at_cut_upper_bounds=saved(gU),one_minus_cycle_numerator_bounds=saved(uppergap),
     right_endpoint_fallback_wall=31,fallback_pair_bounds=saved(g31),endpoint_norm_bounds=saved(qA),
     endpoint_one_minus_cycle_numerator_bounds=saved(endpointgap),
     monotonicity='For0<A<=lambda<U, -g132 increases and positiveq(lambda)=q144+Cq*lambda^2 decreases; hence4g132^2/(2q) is increasing and bounded above by its strictlysubunit value atU. AtA, root31 hascycle strictlybetween0and1.')
 sources=[CLUSTER,MIXED,f'scripts/{TAG}.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
 return dict(status='PASS_ALL_PARAMETER_FOUR_ORTHOGONAL_EDGE_E4_GATES',sources={p:sha(p)for p in sources},
    orthogonal_support='S16 plus anythreeofwalls141..144; positive rank19',
    general_normal='alpha*(-n_i+lambda*C), alpha>0, lambda real; positive norm required',
    normalization_proof='The orthogonal complement is exactlyspan(C,n_i). Acuteness to n_i makesits coefficientnegative; zero coefficientwouldgiveatimelikevector.',
    corner_norm=str(Cq),corner_E4_pair=str(ce),E4_cut_lambda_uppers=list(map(str,L)),
    direct_acute_obstructions=ordinary,forbidden_angle_axis=tail,
    scope='Every acuteCoxeter newwall retainingthe145walls andrighttoS16andanythreeofthefour residualA1caps missesE4. No norm,angle,heightcutoff. This doesnotexclude bent finitecontacts orwallsnotrighttoS16.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();data=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert data==json.loads(p.read_text())
 else:p.write_text(json.dumps(data,indent=2)+'\n')
 print(data['status'])
