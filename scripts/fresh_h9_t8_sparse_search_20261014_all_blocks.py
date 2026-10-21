#!/usr/bin/env python3
"""Free all three remaining internal block labels in the sparse T8 gate.

The literal domains persist by Gram monotonicity. For p>=6 every foreign
contact to that block is excluded by a positive-vertex Schur budget.
Only p=5 and a foreign contact there change the frozen polynomial bank.
"""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product
from math import comb
import sympy as sy
import fresh_h9_star_cube_metric_20261014 as G
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_sparse_search_20261014_all_blocks'
BASE='results/fresh_h9_t8_asymmetric_20261014.json';STATE='results/fresh_h9_t8_sparse_search_20261014_case_08.json'
SPARSE='results/fresh_h9_t8_sparse_frontier_20261014.json';REVIEW='results/fresh_h9_t8_sparse_frontier_20261014_review.json';MIDDLE='results/fresh_h9_t8_middle_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def sign(z):
 a=(z[0],z[1]);b=(z[2],z[3]);sa=G.sign(a);sb=G.sign(b)
 if not sa:return sb
 if not sb or sa==sb:return sa
 diff=G.sub(G.mul(a,a),G.scale(G.mul(b,b),5));ss=G.sign(diff)
 assert ss!=0
 return sa if ss>0 else sb
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,s):return tuple(x*s for x in a)
def run():
 base,state,sparse,review,middle=[json.loads((ROOT/p).read_text())for p in[BASE,STATE,SPARSE,REVIEW,MIDDLE]]
 for obj in[base,state,sparse,review,middle]:
  for pp,hh in obj['sources'].items():assert sha(pp)==hh,pp
 assert base['status']=='PASS_T8_ALL_OWN_CONTACT_PATTERNS_N_GE_4_EXCLUDED'
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 assert sparse['status']=='PASS_SPARSE_LITERAL_DEPENDENCY_AND_SAVED_STATE_EXCLUSION'
 assert review['status']=='PASS_INDEPENDENT_SPARSE_T8_DEPENDENCY_REDUCTION'
 assert middle['status']=='PASS_T8_TWO_MIDDLE_AND_LAST_BLOCK_LABELS_GE_4_ALL_OWN_PATTERNS_EXCLUDED'
 fixed=[row for row in sparse['fixed_labels']if row[:2]!=[4,5]]
 bounds=sparse['finite_label_lower_bounds']+[[4,5,4]]
 labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
 assert all(labels[a,b]==m for a,b,m in fixed)
 assert all(labels[a,b]>=m for a,b,m in bounds)
 # Increasing either middle internal label preserves all sparse literal
 # non-PD witnesses and raises every lower Schur-cost bound. At label>=6
 # any foreign contact to that block costs>=1 before positive owncost.
 assert F(1,4)/(1-F(3,4))==1 and F(1,4)+1>1
 # The exact four-dimensional field is represented independently of the
 # old Q(sqrt2) field encoding. theta=sqrt2+sqrt5 has x^4-14x^2+9=0.
 rt=sy.sqrt(2);rv=sy.sqrt(5);K=sy.QQ.algebraic_field(rt+rv)
 assert K.mod.to_list()==[1,0,-14,0,9]
 def q4(xx):
  z=[F(str(vv))for vv in K.from_sympy(sy.sympify(xx)).to_list()];z=[F(0)]*(4-len(z))+z
  c3,c2,c1,c0=z
  return c0+7*c2,c1+17*c3,c1+11*c3,2*c2
 assert q4(rt)==(0,1,0,0)and q4(rv)==(0,0,1,0)and q4(rt*rv)==(0,0,0,1)
 golden=(1+rv)/4
 assert sy.simplify(1/(4*(1-golden**2))-(5+rv)/10)==0
 assert sy.simplify(golden/(4*(1-golden**2))-(5+3*rv)/20)==0
 assert sy.simplify(1/(4*(1-golden))-(3+rv)/4)==0
 assert sy.simplify(1/(2*(1-golden))-1-(1+rv)/2)==0
 e,a,b,d,c,u,v,p,delta,ac,bc,Q=sy.symbols('e a b d c u v p delta ac bc Q',real=True)
 B=1-c*c;P=e*e+2*c*e*a+a*a;C=b*b+2*c*b*d+d*d;R=e*b+a*d+c*(e*d+a*b);h=delta+p
 A=((1-u)*(1-v)-h*h)*B-(1-u)*C-(1-v)*P-2*h*R+(e*d-a*b)**2
 N=-Q*A-B*((1-u)*bc*bc+(1-v)*ac*ac+2*h*ac*bc)+C*ac*ac+P*bc*bc-2*R*ac*bc
 X=(1-u)*B-P;Y=(1-v)*B-C;Z=-h*B-R
 assert sy.expand(X*Y-Z*Z-B*A)==0
 assert sy.expand(-Q*B*(X*Y-Z*Z)-X*(B*bc)**2-Y*(B*ac)**2+2*Z*(B*ac)*(B*bc)-B*B*N)==0
 def bernstein(poly,variables,intervals):
  xs=sy.symbols('x0:'+str(len(variables)));sub={vv:lo+(hi-lo)*xx for vv,xx,(lo,hi)in zip(variables,xs,intervals)}
  expr=sy.expand(poly.subs(sub,simultaneous=True))
  if variables:
   PP=sy.Poly(expr,*xs,extension=rt+rv);degrees=list(PP.degree_list());terms=[(ix,q4(co))for ix,co in PP.terms()]
  else:degrees=[];terms=[((),q4(expr))]
  coefficients=[]
  for index in product(*(range(nn+1)for nn in degrees)):
   value=(F(0),)*4
   for ix,co in terms:
    if all(ii<=kk for ii,kk in zip(ix,index)):
     factor=F(1)
     for ii,kk,nn in zip(ix,index,degrees):factor*=F(comb(kk,ii),comb(nn,ii))
     value=add(value,scale(co,factor))
   assert sign(value)>0,(variables,index,value)
   coefficients.append(dict(index=index,value=value))
  return dict(degrees=degrees,power_terms=[dict(index=ix,value=co)for ix,co in terms],positive_coefficients=coefficients)
 syms={str(xx):xx for xx in[e,a,b,d,c,delta]};rows=[];inherited=[];all_values=[];count=0
 for qsecond in[4,5]:
  for block in base['blocks']:
   sa,sb=block['canonical_own_supports'];active=[xx for xx,root,ss in[(e,12,sa),(a,13,sa),(b,12,sb),(d,13,sb)]if root in ss]
   zero={xx:0 for xx in[e,a,b,d]if xx not in active}
   if not block['id'].startswith('opposite'):zero[delta]=0
   for case in block['cases']:
    ia,ib=case['foreign_contacts'];ma,mb=case['third_contact_labels']
    key=dict(block_id=block['id'],foreign_contacts=[ia,ib],third_contact_labels=[ma,mb])
    key['second_middle_label']=qsecond
    if ia not in[4,5]and ib not in[4,5]:inherited.append(key);continue
    if qsecond==5 and ia not in[7,8]and ib not in[7,8]:inherited.append(key);continue
    def cost(i):return 0 if i is None else sy.Rational(1,2)if i in[7,8]and qsecond==4 else(5+rv)/10
    ua=cost(ia);vb=cost(ib)
    pp=0 if ia is None or ib is None or ((ia in[4,5])!=(ib in[4,5]))else ua if ia==ib else rt/4 if ia in[7,8]and qsecond==4 else(5+3*rv)/20
    aa=-(0 if ma==2 else sy.Rational(1,2))-((3+rv)/4 if ia in[4,5]else 0)
    bb=-(0 if mb==2 else sy.Rational(1,2))-((3+rv)/4 if ib in[4,5]else 0)
    sub=zero|{u:ua,v:vb,p:pp,ac:aa,bc:bb,Q:sy.Rational(25,19)+(1+rv)/2}
    for index,vars0 in[(ia,[e,a]),(ib,[b,d])]:
     if index is not None:
      for xx in vars0:
       if xx in active:sub[xx]=sy.Rational(1,2)
    nn=sy.expand(N.subs(sub,simultaneous=True));an=sy.expand(A.subs(sub,simultaneous=True))
    variables=[xx for xx in[e,a,b,d,c,delta]if xx in nn.free_symbols or xx in an.free_symbols]
    def expr(z):return sy.Rational(z[0])+sy.Rational(z[1])*rt
    intervals=[tuple(expr(z)for z in block['parameter_intervals'][str(xx)])for xx in variables]
    pn=bernstein(nn,variables,intervals);pa=bernstein(-an,variables,intervals)
    count+=len(pn['positive_coefficients'])+len(pa['positive_coefficients'])
    all_values.extend(rec['value']for poly in[pn,pa]for rec in poly['positive_coefficients'])
    rows.append(key|dict(free_variables=list(map(str,variables)),variable_intervals=[[q4(lo),q4(hi)]for lo,hi in intervals],
     foreign_costs=[q4(ua),q4(vb),q4(pp)],third_pairings=[q4(aa),q4(bb)],
     positive_residual_numerator=pn,negative_anchor_numerator=pa))
   print(qsecond,block['id'],'new cases so far',len(rows),flush=True)
 assert len(rows)==384 and len(inherited)==432
 minimum=all_values[0]
 for value in all_values:
  if sign(add(value,scale(minimum,-1)))<0:minimum=value
 paths=[BASE,STATE,SPARSE,REVIEW,MIDDLE,f'scripts/{TAG}.py','scripts/fresh_h9_star_cube_metric_20261014.py']
 return G.enc(dict(status='PASS_T8_ALL_THREE_INTERNAL_LABELS_GE_4_SPARSE_FAMILY_EXCLUDED',sources={pp:sha(pp)for pp in paths},
  source_mandatory_labels=state['solution']['pair_labels'],fixed_labels=fixed,finite_label_lower_bounds=bounds,
  field='Q(sqrt2,sqrt5), basis1,sqrt2,sqrt5,sqrt10',
  inherited_domain_certificate=SPARSE,independent_domain_review=REVIEW,
  domain_extension='Increasing4-5 from4 preserves every literal sparse vertex guard and lower Schur cost; wall14 norm/projections are recomputed.',
  changed_projection=dict(Q_at_first_t_lower=q4(sy.Rational(25,19)+(1+rv)/2),foreign_self_cost=q4((5+rv)/10),foreign_opposite_overlap=q4((5+3*rv)/20),third_overlap=q4((3+rv)/4)),
  label_split=dict(first4='Frozen sparse theorem',first5_second4='256 new cases with a foreign4/5 contact; without that contact Q only increases',
   first5_second5='128 new mixed-foreign cases; other cases inherit the second4 bank or old middle bank with largerQ',
   first5_second_ge6='Foreign7/8 forbidden; remaining cases already covered',
   first_ge6='Foreign4/5 forbidden; other cases inherit the old middle bank with Q>=old lower value; strict positivity persists because A<0'),
  exact_general_determinant_identities=True,case_count=len(rows),inherited_case_count=len(inherited),positive_coefficient_count=count,
  minimum_positive_coefficient=minimum,inherited_cases=inherited,cases=rows,
  scope='Compact H9 realizations of the literal T8 incidence with first1-2 finite>=7, all three other internal labels4-5,7-8,12-13 finite>=4, and14-1=2,14-2=14-4=14-5=3. All other labels free within compact incidence requirements. All nine anchor own patterns included; no finite upper label cutoff and no dotted distance fixed. Whole T incidence remains open.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if args.verify:assert data==json.loads(out.read_text())
 else:out.write_text(json.dumps(data,indent=2)+'\n')
 print(data['status'],data['case_count'],data['positive_coefficient_count'],data['minimum_positive_coefficient'],flush=True)
