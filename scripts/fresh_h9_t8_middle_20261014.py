#!/usr/bin/env python3
"""Free the second middle block I2(p), p>=4, in the complete T8 own-pattern gate.

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
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_middle_20261014'
BASE='results/fresh_h9_t8_asymmetric_20261014.json';STATE='results/fresh_h9_t8_asymmetric_20261014_search_case_08.json'
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
 base,state=[json.loads((ROOT/p).read_text())for p in[BASE,STATE]]
 for obj in[base,state]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 assert base['status']=='PASS_T8_ALL_OWN_CONTACT_PATTERNS_N_GE_4_EXCLUDED'
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 old={(a,b):m for a,b,m in base['source_mandatory_labels']};new={(a,b):m for a,b,m in state['solution']['pair_labels']}
 free={tuple(p)for p in base['free_mandatory_pairs']}|{(7,8)}
 assert new[7,8]==5 and all(new[p]==m for p,m in old.items()if p not in free)
 assert len(base['cross_right_guards'])==24 and len(base['own_support_patterns'])==9
 # All old guards use a middle block of label4; increasing that one entry
 # preserves each nonpositive relaxed Gram witness and raises Schur costs.
 # A foreign pi/3 contact to I2(p), p>=6, costs at least1, before owncost.
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
 for block in base['blocks']:
  sa,sb=block['canonical_own_supports'];active=[xx for xx,root,ss in[(e,12,sa),(a,13,sa),(b,12,sb),(d,13,sb)]if root in ss]
  zero={xx:0 for xx in[e,a,b,d]if xx not in active}
  if not block['id'].startswith('opposite'):zero[delta]=0
  for case in block['cases']:
   ia,ib=case['foreign_contacts'];ma,mb=case['third_contact_labels']
   key=dict(block_id=block['id'],foreign_contacts=[ia,ib],third_contact_labels=[ma,mb])
   if ia not in[7,8]and ib not in[7,8]:inherited.append(key);continue
   def cost(i):return 0 if i is None else sy.Rational(1,2)if i in[4,5]else(5+rv)/10
   ua=cost(ia);vb=cost(ib)
   pp=0 if ia is None or ib is None or ((ia in[4,5])!=(ib in[4,5]))else ua if ia==ib else rt/4 if ia in[4,5]else(5+3*rv)/20
   aa=-(0 if ma==2 else sy.Rational(1,2))-(sy.Rational(1,2)+rt/4 if ia in[4,5]else 0)
   bb=-(0 if mb==2 else sy.Rational(1,2))-(sy.Rational(1,2)+rt/4 if ib in[4,5]else 0)
   sub=zero|{u:ua,v:vb,p:pp,ac:aa,bc:bb,Q:sy.Rational(25,19)+rt/2}
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
  print(block['id'],'new cases so far',len(rows),flush=True)
 assert len(rows)==256 and len(inherited)==152
 minimum=all_values[0]
 for value in all_values:
  if sign(add(value,scale(minimum,-1)))<0:minimum=value
 paths=[BASE,STATE,f'scripts/{TAG}.py','scripts/fresh_h9_star_cube_metric_20261014.py']
 return G.enc(dict(status='PASS_T8_TWO_MIDDLE_AND_LAST_BLOCK_LABELS_GE_4_ALL_OWN_PATTERNS_EXCLUDED',sources={pp:sha(pp)for pp in paths},
  source_mandatory_labels=state['solution']['pair_labels'],free_mandatory_pairs=sorted(free),field='Q(sqrt2,sqrt5), basis1,sqrt2,sqrt5,sqrt10',
  inherited_domain_certificate=BASE,domain_extension='Each literal non-PD witness remains non-PD when the internal label7-8 increases from4; nonnegative inverse Schur costs increase. All old own-pattern/third-contact domains persist.',
  label_split=dict(p4='Entire frozen408-case bank',p5='256 new exact cases plus152 literally unchanged cases',p_ge6='Every nonright anchor foreign contact at7/8 costs>=1 plus positive owncost, so all such contacts are right; the152 old cases suffice.'),
  exact_general_determinant_identities=True,case_count=len(rows),inherited_case_count=len(inherited),positive_coefficient_count=count,
  minimum_positive_coefficient=minimum,inherited_cases=inherited,cases=rows,
  scope='Compact H9 realizations of this source-bound mandatory family: first label1-2 any finite>=7, internal labels7-8 and12-13 any finite>=4, four own contacts of11/17 with12/13 arbitrary finite>=2, other mandatory contacts fixed. Nine nonempty own patterns included; no finite upper label cutoff and no dotted distance fixed. Whole T incidence remains open.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if args.verify:assert data==json.loads(out.read_text())
 else:out.write_text(json.dumps(data,indent=2)+'\n')
 print(data['status'],data['case_count'],data['positive_coefficient_count'],data['minimum_positive_coefficient'],flush=True)
