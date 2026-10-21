#!/usr/bin/env python3
"""Exact full own-contact-pattern rank gate for the saved T8 mandatory family.

All four own contacts are free finite Coxeter labels, with m>=7 and n>=4
on the two variable internal edges. Six continuous polynomial boxes cover
all nine nonempty support patterns; dotted distances are never substituted.
"""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product,combinations
from math import comb
import sympy as sy
import fresh_h9_star_cube_metric_20261014 as G
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_asymmetric_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
STATE='results/fresh_h9_t8_continuation_20261014_next_case_08.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 geo,state=[json.loads((ROOT/p).read_text())for p in[GEOM,STATE]]
 for obj in[geo,state]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 geometry=geo['representatives'][8]['geometry'];V=[v['incident']for v in geometry['vertices']]
 lab={(a,b):m for a,b,m in state['solution']['pair_labels']}
 expected={(0,1):2,(0,2):3,(1,2):7,(1,14):2,(1,15):3,(1,16):3,(2,14):3,(2,15):2,(2,16):2,
  (3,4):3,(3,5):3,(4,5):4,(4,14):3,(5,14):3,(6,7):3,(6,8):3,(7,8):4,(7,15):3,(8,15):3,
  (11,12):2,(11,13):3,(12,13):7,(12,17):3,(13,17):3}
 assert lab==expected
 S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)]
 def actual(s):return next(i for i,v in enumerate(V)if set(s)<=set(v))
 for triple in[[11,12,13],[12,13,17]]:assert triple in geometry['minimal_nonfaces']
 cross=[]
 for i,j in combinations(range(4),2):
  for a0,b0 in product(blocks[i],blocks[j]):
   if i==0:
    support=[next(v for v in blocks[0]if v!=a0),a0,b0]
    cross.append(dict(pair=[a0,b0],support=support,actual_vertex=actual(support),
     witness='3-path determinant 3/4-t<0 for t=cos(pi/m)^2>3/4'))
   else:
    support=[next(v for v in blocks[i]if v!=a0),a0,b0,next(v for v in blocks[j]if v!=b0)]
    M=[[G.O if x==y else G.Z for y in range(4)]for x in range(4)]
    for x,c0 in enumerate([G.A,G.H,G.A]):M[x][x+1]=M[x+1][x]=G.neg(c0)
    w=[G.O,(F(0),F(1)),(F(0),F(1)),G.O]
    assert all(G.total(G.mul(M[x][y],w[y])for y in range(4))==G.Z for x in range(4))
    cross.append(dict(pair=[a0,b0],support=support,actual_vertex=actual(support),
     relaxed_gram=M,positive_null_vector=w,relaxation='Both internal labels lowered to4'))
 assert len(cross)==24
 anchor_budgets=[]
 for leaf,kept in product([11,17],[12,13]):
  removed=25-kept;support=sorted((set(S)-{removed})|{leaf})
  anchor_budgets.append(dict(leaf=leaf,retained_own_root=kept,support=support,actual_vertex=actual(support),
   at_most_one_foreign='owncost>=1/4; each B2 nonright contact costs>=1/2; I2(m) contact costs>1',
   foreign_label_bound='label>=4 costs>=1, plus positive owncost',
   own_label_with_foreign='ownlabel>=4 costs>=1/2 plus foreigncost>=1/2; every own nonright label must therefore be3'))
 third_budgets=[]
 for kept in[4,5]:
  support=sorted((set(S)-{1}-({4,5}-{kept}))|{14})
  third_budgets.append(dict(leaf=14,retained_own_roots=[2,kept],support=support,actual_vertex=actual(support)))
 assert F(1,4)+2*F(1,2)>1 and F(1,2)+F(1,2)==1
 assert F(1,4)+F(1,4)+F(1,2)==1
 # Every support pattern has its own literal vertex witnesses. Only the
 # resulting block algebra is renamed; no incidence automorphism is assumed.
 own_patterns=[];co2={2:G.Z,3:G.H,4:G.A}
 for sa,sb in product([(12,),(13,),(12,13)],repeat=2):
  labels=dict(lab)
  for leaf,ss in[(11,sa),(17,sb)]:
   for root in[12,13]:labels[tuple(sorted((leaf,root)))]=3 if root in ss else 2
  upper=[]
  for anchor in[11,17]:
   pair=tuple(sorted((anchor,14)))
   for vi,vertex in enumerate(V):
    if not set(pair)<=set(vertex):continue
    assert not {1,2}<=set(vertex) and not {12,13}<=set(vertex)
    test=labels|{pair:4}
    M=[[G.O if a0==b0 else G.neg(co2[test.get(tuple(sorted((a0,b0))),2)])for b0 in vertex]for a0 in vertex]
    pp=G.pivots(M)
    if G.sign(pp[-1])<=0:
     upper.append(dict(pair=pair,test_label=4,actual_vertex=vi,vertex=vertex,gram=M,
      leading_order=len(pp),last_pivot=pp[-1]));break
   else:raise AssertionError(('no fresh upper witness',sa,sb,pair))
  shared=sorted(set(sa)&set(sb));mutual=dict(pair=[11,17],actual_vertex=actual([11,17]),domain='any finite label; delta in[0,1)')
  if shared:
   tri=[11,shared[0],17];M=[[G.O if i==j else G.neg(G.H)for j in range(3)]for i in range(3)]
   assert all(G.total(row)==G.Z for row in M)
   mutual.update(domain='right',support=tri,actual_vertex=actual(tri),test_label=3,gram=M,positive_null_vector=[1,1,1])
  family='double'if len(sa)==len(sb)==2 else'asymmetric'if len(sa)+len(sb)==3 else'same'if sa==sb else'opposite'
  own_patterns.append(dict(own_supports=[sa,sb],polynomial_family=family,
   relaxed_own_labels=[[leaf,root,labels[tuple(sorted((leaf,root)))]]for leaf,root in product([11,17],[12,13])],
   anchor_pair_guard=mutual,third_contact_upper_guards=upper))
 assert len(own_patterns)==9
 # cos^2(pi/7)>81/100. U6(x)=f(x^2); f is strictly increasing above3/4.
 x=sy.symbols('x');cheb=sy.Poly(sy.chebyshevu(6,x),x)
 f=64*x**3-80*x**2+24*x-1
 assert cheb.as_expr()==f.subs(x,x*x)
 assert f.subs(x,sy.Rational(81,100))==-sy.Rational(559,15625)
 assert sy.diff(f,x).subs(x,sy.Rational(3,4))==12
 assert sy.diff(f,x,2).subs(x,sy.Rational(3,4))==128
 # Universal determinant identity. P=(e,a), C=(b,d) are own contacts on
 # roots12,13; all symbols are magnitudes. Foreign overlaps contribute p.
 e,a,b,d,c,u,v,p,delta,ac,bc,Q=sy.symbols('e a b d c u v p delta ac bc Q',real=True)
 B=1-c*c;Pcost=e*e+2*c*e*a+a*a;Ccost=b*b+2*c*b*d+d*d
 Rcross=e*b+a*d+c*(e*d+a*b);h=delta+p
 A=((1-u)*(1-v)-h*h)*B-(1-u)*Ccost-(1-v)*Pcost-2*h*Rcross+(e*d-a*b)**2
 N=-Q*A-B*((1-u)*bc*bc+(1-v)*ac*ac+2*h*ac*bc)+Ccost*ac*ac+Pcost*bc*bc-2*Rcross*ac*bc
 X=(1-u)*B-Pcost;Y=(1-v)*B-Ccost;Z=-h*B-Rcross
 # Direct polynomial determinants avoid any symbolic rational simplifier.
 scaled_det=-Q*B*(X*Y-Z*Z)-X*(B*bc)**2-Y*(B*ac)**2+2*Z*(B*ac)*(B*bc)
 assert sy.expand(scaled_det-B*B*N)==0
 assert sy.expand(X*Y-Z*Z-B*A)==0
 rt=sy.sqrt(2);K=sy.QQ.algebraic_field(rt);assert K.mod.to_list()==[1,0,-2]
 def q2(xx):
  z=[F(str(vv))for vv in K.from_sympy(sy.sympify(xx)).to_list()];z=[F(0)]*(2-len(z))+z
  return z[1],z[0]
 def bernstein(poly,variables,intervals):
  xs=sy.symbols('x0:'+str(len(variables)))
  sub={vv:lo+(hi-lo)*xx for vv,xx,(lo,hi)in zip(variables,xs,intervals)}
  expr=sy.expand(poly.subs(sub,simultaneous=True))
  if variables:
   PP=sy.Poly(expr,*xs,extension=rt);degrees=list(PP.degree_list());terms=[(ix,q2(co))for ix,co in PP.terms()]
  else:degrees=[];terms=[((),q2(expr))]
  coeff=[]
  for index in product(*(range(nn+1)for nn in degrees)):
   value=G.Z
   for ix,co in terms:
    if all(ii<=kk for ii,kk in zip(ix,index)):
     factor=F(1)
     for ii,kk,nn in zip(ix,index,degrees):factor*=F(comb(kk,ii),comb(nn,ii))
     value=G.add(value,G.scale(co,factor))
   assert G.sign(value)>0,(variables,index,value)
   coeff.append(dict(index=index,value=value))
  return dict(degrees=degrees,power_terms=[dict(index=ix,value=co)for ix,co in terms],positive_coefficients=coeff)
 foreign=[None,4,5,7,8];polynomial_blocks=[];count=0;all_values=[];total_cases=0
 definitions=[('asymmetric',[13],[12,13]),('double',[12,13],[12,13]),
  ('same_high',[13],[13]),('opposite_high',[13],[12]),('same_low',[13],[13]),('opposite_low',[13],[12])]
 for name,sa,sb in definitions:
  low=name.endswith('_low');high=name.endswith('_high');opposite=name.startswith('opposite')
  active=[vv for vv,root,ss in[(e,12,sa),(a,13,sa),(b,12,sb),(d,13,sb)]if root in ss]
  base={vv:0 for vv in[e,a,b,d]if vv not in active}
  if not opposite:base[delta]=0
  lo_own=rt/2 if low else sy.Rational(1,2);lo_c=sy.Rational(9,10)if high else rt/2
  ranges={vv:(lo_own,sy.Integer(1))for vv in active};ranges[c]=(lo_c,sy.Integer(1))
  if opposite:ranges[delta]=(sy.Integer(0),sy.Integer(1))
  rows=[]
  for ia,ib,ma,mb in product([None]if low else foreign,[None]if low else foreign,[2,3],[2,3]):
   ua=0 if ia is None else sy.Rational(1,2);vb=0 if ib is None else sy.Rational(1,2)
   pp=0 if ia is None or ib is None or ((ia in[4,5])!=(ib in[4,5]))else sy.Rational(1,2)if ia==ib else rt/4
   aa=-(0 if ma==2 else sy.Rational(1,2))-(sy.Rational(1,2)+rt/4 if ia in[4,5]else 0)
   bb=-(0 if mb==2 else sy.Rational(1,2))-(sy.Rational(1,2)+rt/4 if ib in[4,5]else 0)
   sub=base|{u:ua,v:vb,p:pp,ac:aa,bc:bb,Q:sy.Rational(25,19)+rt/2}
   if ia is not None:
    for vv in[e,a]:
     if vv in active:sub[vv]=sy.Rational(1,2)
   if ib is not None:
    for vv in[b,d]:
     if vv in active:sub[vv]=sy.Rational(1,2)
   nn=sy.expand(N.subs(sub,simultaneous=True));an=sy.expand(A.subs(sub,simultaneous=True))
   variables=[vv for vv in[e,a,b,d,c,delta]if vv in nn.free_symbols or vv in an.free_symbols]
   intervals=[ranges[vv]for vv in variables]
   pn=bernstein(nn,variables,intervals);pa=bernstein(-an,variables,intervals)
   count+=len(pn['positive_coefficients'])+len(pa['positive_coefficients'])
   all_values.extend(rec['value']for po in[pn,pa]for rec in po['positive_coefficients'])
   rows.append(dict(foreign_contacts=[ia,ib],third_contact_labels=[ma,mb],free_variables=list(map(str,variables)),
    variable_intervals=[[q2(lo),q2(hi)]for lo,hi in intervals],foreign_costs=[q2(ua),q2(vb),q2(pp)],third_pairings=[q2(aa),q2(bb)],
    positive_residual_numerator=pn,negative_anchor_numerator=pa))
  total_cases+=len(rows)
  polynomial_blocks.append(dict(id=name,canonical_own_supports=[sa,sb],parameter_intervals={str(vv):[q2(lo),q2(hi)]for vv,(lo,hi)in ranges.items()},cases=rows))
  print('block',name,'cases',len(rows),'positive coefficients so far',count,flush=True)
 assert total_cases==408
 minimum=all_values[0]
 for value in all_values:
  if G.sign(G.sub(value,minimum))<0:minimum=value
 assert G.sign(minimum)>0
 paths=[GEOM,STATE,f'scripts/{TAG}.py','scripts/fresh_h9_star_cube_metric_20261014.py']
 return G.enc(dict(status='PASS_T8_ALL_OWN_CONTACT_PATTERNS_N_GE_4_EXCLUDED',sources={pp:sha(pp)for pp in paths},
  source_mandatory_labels=state['solution']['pair_labels'],free_mandatory_pairs=[[1,2],[12,13],[11,12],[11,13],[12,17],[13,17]],
  positive_support=S,positive_support_actual_vertex=actual(S),positive_support_determinant='(1-t)*(1-c^2)/4>0',
  cross_right_guards=cross,anchor_budgets=anchor_budgets,third_budgets=third_budgets,own_support_patterns=own_patterns,
  own_support_completeness=dict(missing_triples=[[11,12,13],[12,13,17]],reason='Each missing triple is Lanner in a compact realization; both own contacts right would make it positive definite.',
   low_single_reason='For n=4,5,6 a single own label3 gives determinant3/4-c^2>=0, so each single own label is>=4. A foreign contact then violates the actual-vertex budget.'),
  cosine_threshold=dict(U6_coefficients=list(map(int,reversed(cheb.all_coeffs()))),f_at_81_over_100=F(-559,15625),fprime_at_3_over_4=12,fsecond_at_3_over_4=128,
   conclusion='cos^2(pi/7)>81/100, hence cos(pi/7)>9/10; higher finite labels increase these values'),
  exact_general_determinant_identities=True,field='Q(sqrt2); coefficients [rational,sqrt2 coefficient]',
  universal_formulas=dict(P='e^2+2*c*e*a+a^2',C='b^2+2*c*b*d+d^2',R='e*b+a*d+c*(e*d+a*b)',B='1-c^2',h='delta+p',
   A='((1-u)*(1-v)-h^2)*B-(1-u)*C-(1-v)*P-2*h*R+(e*d-a*b)^2',
   N='-Q*A-B*((1-u)*bc^2+(1-v)*ac^2+2*h*ac*bc)+C*ac^2+P*bc^2-2*R*ac*bc'),
  Q_at_lower_t=[F(25,19),F(1,2)],case_count=total_cases,positive_coefficient_count=count,minimum_positive_coefficient=minimum,blocks=polynomial_blocks,
  first_parameter_extension='Residual determinant=N/B; derivative in t is -A/[4B(1-t)^2]>0 because A<0 throughout each remaining box.',
  scope='Compact H9 Coxeter realizations of the saved T8 family with m>=7 on1-2, n>=4 on12-13, all four anchor-own finite labels arbitrary>=2, and other mandatory labels fixed. All nine nonempty own-support patterns are included. No finite upper label cutoff, no dotted distance fixed, no incidence automorphism assumed. The positive rank8 block leaves residual rank<=2, contradicted by the selected3x3 determinant. Whole T incidence remains open.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if args.verify:assert data==json.loads(out.read_text())
 else:out.write_text(json.dumps(data,indent=2)+'\n')
 print(data['status'],data['case_count'],data['positive_coefficient_count'],data['minimum_positive_coefficient'],flush=True)
