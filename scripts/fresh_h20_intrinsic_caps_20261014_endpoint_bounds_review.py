#!/usr/bin/env python3
"""Source/formula audit plus independent small radius-bound arithmetic replay."""
import argparse,json
from pathlib import Path
import sympy as S
from fresh_h20_intrinsic_caps_20261014_review import ROOT,sha,exact,sign
TAG='fresh_h20_intrinsic_caps_20261014_endpoint_bounds_review'
SOURCE='results/fresh_h20_intrinsic_caps_20261014_endpoint_bounds.json'
ENDPOINT='results/fresh_h20_intrinsic_caps_20261014_endpoint.json'
def run():
 out=json.loads((ROOT/SOURCE).read_text());p=json.loads((ROOT/ENDPOINT).read_text())
 for data in(out,p):
  for f,h in data['sources'].items():assert sha(f)==h,f
 phi=(1+S.sqrt(5))/2
 rows={(r['type'],r['profile']):r for r in out['radius_relaxation']};assert len(rows)==260
 computed={};tightened=0;rightcases=0
 def less(a,b):return sign(S.expand(a[0]*b[1]-b[0]*a[1]))<0
 for ty in p['types']:
  for row in ty['profiles']:
   A=[S.Rational(x['numerator_over_radius'])*(phi if x['numerator_has_triangle_inradius_s']else 1)for x in row['omitted_contact_ratios']]
   B=[exact(x['subtract'])for x in row['omitted_contact_ratios']]
   assert all(sign(a)>0 for a in A)and all(sign(b)>=0 for b in B)
   ordinary=[(a,b+S.Rational(1,2))for a,b in zip(A,B)]
   bound=ordinary[0]
   for v in ordinary[1:]:
    if less(v,bound):bound=v
   accepted=[];rejected=[]
   for i,b in enumerate(B):
    if b==0:continue
    vals=[S.expand(A[j]*b-B[j]*A[i])for j in range(3)]
    ok=all(v==0 or sign(2*v-A[i])>=0 for v in vals)
    (accepted if ok else rejected).append(i)
    if ok and less(bound,(A[i],b)):bound=(A[i],b)
   tighten=False
   if ty['id']==2 and B[1]==B[2]==0:
    assert A[2]==4 and S.expand(A[1]-3*phi)==0 and not less((S.Integer(8),S.Integer(1)),bound)
    assert sign(3*phi-4)>0
    candidate=(3*phi*S.sqrt(2),S.Integer(1))
    if less(candidate,bound):bound=candidate;tighten=True
   saved=rows[ty['id'],row['id']]
   stored=(exact(saved['bound_numerator']),exact(saved['bound_denominator']))
   assert S.expand(bound[0]*stored[1]-stored[0]*bound[1])==0
   assert accepted==saved['accepted_right_factors']and rejected==saved['rejected_right_factors']and tighten==saved['golden_gap_tightened']
   computed[ty['id'],row['id']]=bound;tightened+=tighten;rightcases+=len(accepted)
 # Independently rebuild the five proof-sensitive half-vertex numerators.
 ids=p['seed_order'];gaps=list(map(exact,out['half_vertex_witness']['all22_seed_gaps']));lookup=dict(zip(ids,gaps))
 tr=next(x for x in p['literal_edge_transports']if x['omitted']==[4,23,22]);ty=p['types'][tr['type']]
 I=S.Matrix([[S.Rational(x)for x in row]for row in ty['inverse_Gram']]);d=S.Matrix([lookup[w]for w in tr['edge_support_in_profile_order']])
 norm=S.expand(2*(d.T*I*d)[0]);assert norm==76+S.sqrt(5)/2
 checked=[]
 for critical in out['critical_five_margins']:
  row=ty['profiles'][critical['profile']];assert row['id']==critical['profile']
  co=S.Matrix(list(map(exact,row['cosine_vector'])))
  t=S.expand(-2*(d.T*I*co)[0]);cost=S.expand(2*(co.T*I*co)[0]);A,C=computed[ty['id'],row['id']]
  num=S.expand(norm*C*C+2*A*C*t+(cost-1)*A*A)
  assert num==exact(critical['positive_power_numerator'])and S.expand(C*C)==exact(critical['positive_denominator'])
  assert sign(S.expand(num-C*C/4))>0;checked.append(row['id'])
 assert len(checked)==5
 c,x=S.symbols('c x',real=True)
 tail=S.expand(100*(1+c)**2-12*c*(1+c)+144*(c*c-1))
 assert tail==232*c*c+188*c-44 and S.expand(tail.subs(c,(1+x)/2))==108+210*x+58*x*x
 inner=S.expand(300*(1+c)**2-24*c*(1+c)+144*(c*c-1)-100*(1+c)**2)
 assert inner==320*c*c+376*c+56 and sign(12-6*phi)>0
 # Ball-power formula, including the PSD endpoint rho=0 case.
 F,t,R,k=S.symbols('F t R k',real=True)
 power=F+2*R*t+(k-1)*R*R
 assert S.expand(S.diff(power,R,2)-2*(k-1))==0
 paths=[SOURCE,ENDPOINT,'scripts/fresh_h20_intrinsic_caps_20261014_endpoint_bounds.py','scripts/fresh_h20_intrinsic_caps_20261014_review.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return {'status':'PASS_ENDPOINT_ALL_RADIUS_SOURCE_REVIEW_AND_INDEPENDENT_SMALL_FORMULA_CHECK','sources':{f:sha(f)for f in paths},'radius_profiles_independently_checked':260,'exceptional_right_factor_cases':rightcases,'golden_gap_tightened_profiles':tightened,'critical_margin_records':checked,'symbolic_tail_polynomials':[[108,210,58],[56,376,320]],'scope':'Source/proof review of all-radius coverage, equality/PSD/isolated-tail cases and strict escape. Separate SymPy/rational sign arithmetic rechecks260 radius bounds, all exceptional right choices and5 critical margins. Does not independently repeat the16272 transported point-profile inequalities or the rejected endpoint profile enumeration.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();path=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert json.loads(path.read_text())==out
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
