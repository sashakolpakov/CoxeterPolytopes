#!/usr/bin/env python3
"""Independent Fraction review of the opposite-root T8 family.

No producer arithmetic or SymPy imports. Rechecks literal domains, 400
quadratic Bernstein certificates, and the all-own-label identities.
"""
import argparse,json
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
from fresh_h9_star_cube_metric_20261014_cut_profiles_review import (
 F,Z,O,CO,add,sub,neg,mul,scale,total,inverse,enc,anp)
from fresh_h9_star_cube_metric_20261014_all_cut_profiles_review import sign,det_gauss
from fresh_h9_t8_metric_20261014_review import q2,padd,pneg,pmul,pdet,trim
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_t8_continuation_20261014_opposite'
SOURCE=f'results/{TAG}.json'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
STATE='results/fresh_h9_t8_continuation_20261014_case_08.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def formal_identities():
 # Nine formal indeterminates; coefficients are plain rational numbers.
 zero_key=(0,)*9
 def con(q):return {zero_key:F(q)} if q else {}
 def plus(a,b):
  d=a.copy()
  for k,v in b.items():d[k]=d.get(k,F(0))+v
  return {k:v for k,v in d.items()if v}
 def times(a,b):
  d={}
  for i,x in a.items():
   for j,y in b.items():
    k=tuple(u+v for u,v in zip(i,j));d[k]=d.get(k,F(0))+x*y
  return {k:v for k,v in d.items()if v}
 def sc(a,q):return {k:v*q for k,v in a.items()if v*q}
 def minus(a,b):return plus(a,sc(b,-1))
 def summ(xs):
  out={}
  for x in xs:out=plus(out,x)
  return out
 def prod(*xs):
  out=con(1)
  for x in xs:out=times(out,x)
  return out
 def square(x):return times(x,x)
 def deriv(x,index):
  out={}
  for k,v in x.items():
   if k[index]:
    l=list(k);l[index]-=1;out[tuple(l)]=v*k[index]
  return out
 def at_zero(x,indices):return {k:v for k,v in x.items()if all(k[i]==0 for i in indices)}
 variables=[]
 for i in range(9):
  k=list(zero_key);k[i]=1;variables.append({tuple(k):F(1)})
 a,b,c,u,v,e,A,B,Q=variables
 one=con(1);den=minus(one,square(c));U=minus(one,u);V=minus(one,v)
 L=summ([prod(minus(prod(U,V),square(e)),den),sc(prod(square(a),V),-1),
         sc(prod(square(b),U),-1),sc(prod(e,c,a,b),-2),prod(square(a),square(b))])
 N=summ([sc(prod(Q,L),-1),sc(prod(den,summ([prod(U,square(B)),prod(V,square(A)),sc(prod(e,A,B),2)])),-1),
         prod(square(a),square(B)),prod(square(b),square(A)),sc(prod(c,a,b,A,B),-2)])
 qa=minus(prod(U,den),square(a));qb=minus(prod(V,den),square(b))
 ab=sc(plus(prod(e,den),prod(c,a,b)),-1)
 ac,bc,qq=prod(A,den),prod(B,den),sc(prod(Q,den),-1)
 M=[[qa,ab,ac],[ab,qb,bc],[ac,bc,qq]]
 det={}
 for p in permutations(range(3)):
  term=prod(*(M[i][j]for i,j in enumerate(p)))
  if sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))%2:term=sc(term,-1)
  det=plus(det,term)
 assert det==prod(square(den),N)
 assert minus(prod(qa,qb),square(ab))==prod(den,L)
 derivative=summ([prod(a,plus(prod(Q,minus(V,square(b))),square(B))),prod(c,b,minus(prod(Q,e),prod(A,B)))])
 assert deriv(N,0)==sc(derivative,2)
 assert deriv(L,0)==summ([sc(prod(a,minus(V,square(b))),-2),sc(prod(e,c,b),-2)])
 right=summ([prod(Q,minus(square(c),prod(minus(one,square(a)),minus(one,square(b))))),
             square(minus(prod(a,B),prod(b,A))),sc(prod(a,b,minus(one,c),A,B),2),
             sc(prod(den,plus(square(A),square(B))),-1)])
 assert at_zero(N,[3,4,5])==right
 # Complete-square identity used when exactly one anchor has a foreign contact.
 # z²/2-z/4 = (z-1/4)²/2 -1/32.
 assert F(1,2)*F(1,4)**2==F(1,32)
 return dict(indeterminates=['alpha','beta','c','u','v','delta_plus_foreign_projection','ac','bc','Q'],
             determinant_identity=True,anchor_identity=True,alpha_derivative_identity=True,
             anchor_derivative_identity=True,right_mutual_special_identity=True)

def run():
 data=json.loads((ROOT/SOURCE).read_text());bound=set()
 assert data['status']=='PASS_T8_OPPOSITE_ROOT_ALL_OWN_LABELS_FAMILY_EXCLUDED'
 def bind(p):
  if p in bound:return
  bound.add(p)
  for s,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():
   assert sha(s)==h,s
   if s.endswith('.json'):bind(s)
 bind(SOURCE)
 geo=json.loads((ROOT/GEOM).read_text())['representatives'][8]['geometry']
 state=json.loads((ROOT/STATE).read_text());V=[v['incident']for v in geo['vertices']]
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
 assert labels=={(0,1):2,(0,2):3,(1,2):7,(1,14):2,(1,15):3,(1,16):3,(2,14):3,(2,15):2,(2,16):2,
  (3,4):3,(3,5):3,(4,5):4,(4,14):3,(5,14):3,(6,7):3,(6,8):3,(7,8):4,(7,15):3,(8,15):3,
  (11,12):2,(11,13):3,(12,13):7,(12,17):3,(13,17):2}
 assert data['source_mandatory_labels']==state['solution']['pair_labels']
 assert data['free_mandatory_pairs']==[[1,2],[12,13],[11,13],[12,17]]
 S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)]
 assert data['positive_support']==S and set(S)<=set(V[data['positive_support_actual_vertex']])
 all_cross={(a,b)for i,j in combinations(range(4),2)for a,b in product(blocks[i],blocks[j])}
 seen=set();counts={'high_cross':0,'affine_cross':0}
 for row in data['cross_right_guards']:
  a,b=row['pair'];assert (a,b)in all_cross-seen;seen.add((a,b))
  nodes=row['support'];assert set(nodes)<=set(V[row['actual_vertex']])
  i,j=S.index(a)//2,S.index(b)//2
  if i==0:
   assert nodes==[next(x for x in blocks[0]if x!=a),a,b]
   H=[[(O,),(Z,neg(O)),(Z,)],[(Z,neg(O)),(O,),(neg(CO[3]),)],[(Z,),(neg(CO[3]),),(O,)]]
   assert pdet(H)==(scale(O,F(3,4)),Z,neg(O));counts['high_cross']+=1
  else:
   assert nodes==[next(x for x in blocks[i]if x!=a),a,b,next(x for x in blocks[j]if x!=b)]
   assert labels[tuple(blocks[i])]>=4 and labels[tuple(blocks[j])]>=4
   H=[[O if x==y else Z for y in range(4)]for x in range(4)]
   for k,c in enumerate([CO[4],CO[3],CO[4]]):H[k][k+1]=H[k+1][k]=neg(c)
   null=[O,scale(CO[4],2),scale(CO[4],2),O]
   assert all(total(mul(a,b)for a,b in zip(r,null))==Z for r in H)
   assert H==[[q2(x)for x in r]for r in row['relaxed_gram']]
   assert null==list(map(q2,row['positive_null_vector']));counts['affine_cross']+=1
 assert seen==all_cross and counts=={'high_cross':12,'affine_cross':12}
 B2=[[O if i==j else neg(CO[4])if i//2==j//2 else Z for j in range(4)]for i in range(4)]
 J=inverse(B2)
 assert all(J[i][i]==scale(O,2)for i in range(4))and all(sign(x)>=0 for r in J for x in r)
 seen=set()
 for row in data['anchor_budgets']:
  leaf,kept=row['leaf'],row['retained_own_root'];seen.add((leaf,kept))
  assert (leaf,kept)in[(11,13),(17,12)]
  assert set(row['support'])==(set(S)-({12,13}-{kept}))|{leaf}<=set(V[row['actual_vertex']])
  assert labels[tuple(sorted((leaf,kept)))]==3
 assert seen=={(11,13),(17,12)}
 # Every absent I2 contact costs >1; two B2 contacts cost >=1. One
 # foreign contact costs1/2 and forces the own label to remain3.
 assert F(1,4)/(1-F(81,100))>1
 assert F(1,4)+2*F(1,2)>1 and F(1,2)+F(1,2)==1
 third=set()
 for row in data['third_budgets']:
  a,kept=row['retained_own_roots'];third.add(kept)
  assert row['leaf']==14 and a==2 and kept in[4,5]
  assert set(row['support'])=={2,kept,7,8,12,13,14}<=set(V[row['actual_vertex']])
  assert labels[2,14]==labels[kept,14]==3
 assert third=={4,5}and F(1,4)+F(1,4)+F(1,2)==1
 guard=data['anchor_mutual_upper_guard']
 assert guard['support']==[11,12,17]and set(guard['support'])<=set(V[guard['actual_vertex']])
 assert guard['test_label']==6 and guard['allowed_labels']==[2,3,4,5]
 assert labels[11,12]==2 and labels[12,17]==3
 assert 1-F(3,4)-F(1,4)==0
 upper=[]
 for row in data['third_contact_upper_guards']:
  pair=tuple(row['pair']);assert pair in[(11,14),(14,17)]and row['test_label']==4
  nodes=V[row['actual_vertex']];assert nodes==row['vertex']and set(pair)<=set(nodes)
  assert not {1,2}<=set(nodes)and not {12,13}<=set(nodes)
  lab=labels|{pair:4}
  H=[[O if a==b else neg(CO[lab.get(tuple(sorted((a,b))),2)])for b in nodes]for a in nodes]
  assert H==[[q2(x)for x in r]for r in row['gram']]
  k=row['leading_order'];d=det_gauss([r[:k]for r in H[:k]])
  before=det_gauss([r[:k-1]for r in H[:k-1]])if k>1 else O
  assert sign(before)>0 and sign(d)<=0 and mul(q2(row['last_pivot']),before)==d
  upper.append(dict(pair=list(pair),vertex=row['actual_vertex'],minor_order=k,determinant=enc(d)))
 assert {tuple(x['pair'])for x in upper}=={(11,14),(14,17)}
 # Roots of unity give f(x)=8x³+4x²-4x-1 at x=cos(2pi/7)>1/2.
 # f is increasing there and f(31/50)<0, hence cos²(pi/7)>81/100.
 f=lambda x:8*x**3+4*x*x-4*x-1
 assert f(F(31,50))==F(-559,15625)
 assert 24*F(1,2)**2+8*F(1,2)-4==6 and 48*F(1,2)+8==32
 assert (1+F(31,50))/2==F(81,100)and F(9,10)**2==F(81,100)
 assert F(data['threshold_cubic_value'])==F(-559,15625)
 assert data['chebyshev_U6']==[-1,0,24,0,-80,0,64]
 Q=add(CO[4],scale(O,F(25,19)));assert sign(sub(Q,scale(O,2)))>0
 assert anp(data['Q_at_t_lower'])==Q and anp(data['Q_minus_two'])==sub(Q,scale(O,2))
 patterns={None:[Z]*4}
 for i,wall in enumerate([4,5,7,8]):
  row=[Z]*4;row[i]=neg(CO[3]);patterns[wall]=row
 def dot(x,y):return total(mul(x[i],mul(J[i][j],y[j]))for i in range(4)for j in range(4))
 third_vector=[neg(CO[3]),neg(CO[3]),Z,Z]
 saved={(tuple(r['foreign_contacts']),r['anchor_mutual_label'],tuple(r['third_contact_labels'])):r for r in data['cases']}
 assert len(saved)==len(data['cases'])==400
 den=(O,Z,neg(O));records=[]
 def poly_scale(x,s):return tuple(mul(v,s)for v in x)
 def evaluate(x,t):
  out=Z
  for a in reversed(x):out=add(scale(out,t),a)
  return out
 def bernstein(poly):
  assert len(poly)<=3
  a=list(poly)+[Z]*(3-len(poly));lo=evaluate(a,F(9,10))
  mid=add(lo,scale(add(a[1],scale(a[2],F(9,5))),F(1,20)))
  return [lo,mid,total(a)]
 for ia,ib,mutual,ma,mb in product(patterns,patterns,[2,3,4,5],[2,3],[2,3]):
  x,y=patterns[ia],patterns[ib];u,v,p=dot(x,x),dot(y,y),dot(x,y)
  e=add(CO[mutual],p);A=sub(neg(CO[ma]),dot(x,third_vector));B=sub(neg(CO[mb]),dot(y,third_vector))
  qa=padd(poly_scale(den,sub(O,u)),(scale(O,F(-1,4)),))
  qb=padd(poly_scale(den,sub(O,v)),(scale(O,F(-1,4)),))
  ab=padd(poly_scale(den,neg(e)),(Z,scale(O,F(-1,4))))
  ac,bc,qq=poly_scale(den,A),poly_scale(den,B),poly_scale(den,neg(Q))
  M=[[qa,ab,ac],[ab,qb,bc],[ac,bc,qq]]
  row=saved[(ia,ib),mutual,(ma,mb)]
  # Compare the independent determinant expansion, not a copied formula.
  residual=tuple(map(anp,row['residual_numerator']));anchor=tuple(map(anp,row['anchor_numerator']))
  assert pdet(M)==pmul(residual,pmul(den,den))
  assert pdet([[qa,ab],[ab,qb]])==pmul(anchor,den)
  rb,abern=bernstein(residual),bernstein(anchor)
  assert all(sign(a)>0 for a in rb)and all(sign(a)<0 for a in abern)
  assert rb==list(map(anp,row['residual_bernstein']))and abern==list(map(anp,row['anchor_bernstein']))
  records.append(dict(foreign_contacts=[ia,ib],mutual_label=mutual,third_labels=[ma,mb],
   positive_residual_bernstein=list(map(enc,rb)),negative_anchor_bernstein=list(map(enc,abern))))
 formal=formal_identities()
 # Bounds for the continuous own-label extension, with all domains already
 # proved above. Foreign-free anchors alone can vary alpha,beta above1/2.
 assert F(2,8)-F(1,32)==F(7,32)>0
 assert 2*F(1,2)-F(1,2)**2==F(3,4)>0
 assert 2*(F(81,100)-F(9,16))-F(19,100)*F(1,2)==F(2,5)>0
 ext=data['all_own_extension']
 assert F(ext['one_free_one_foreign_lower'])==F(7,32)
 assert F(ext['both_free_positive_mutual_bracket_lower'])==F(3,4)
 assert F(ext['both_free_right_mutual_lower'])==F(2,5)
 paths=[SOURCE,GEOM,STATE,f'scripts/{TAG}_review.py',
  'scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py',
  'scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.py',
  'scripts/fresh_h9_t8_metric_20261014_review.py']
 return dict(status='PASS_INDEPENDENT_T8_OPPOSITE_ROOT_ALL_OWN_LABELS',sources={p:sha(p)for p in paths},
  recursively_checked_manifests=sorted(bound),cross_guards=counts,opposite_anchor_budgets=[[11,13],[17,12]],
  third_contact_failure_minors=upper,cosine_square_lower_bound='81/100',cosine_lower_bound='9/10',
  residual_cases=400,exact_cases=records,formal_all_variable_identities=formal,
  extension_bounds=dict(one_free_one_foreign='7/32',both_free_nonright_mutual_bracket='3/4',both_free_right_mutual='2/5'),
  analytic_argument='A foreign contact forces its own contact label to3. For one free variable, N_alpha/2>7/32 by completing a square; for two free variables and nonright mutual contact both derivatives are positive. For two foreign-free anchors with right mutual contact N>2/5 directly. The anchor minor decreases under every allowed own-contact increase, so the determinant also increases with t. The400 positive Bernstein certificates therefore extend to all finite own labels>=3 and both internal labels>=7.',
  scope='Only the specified opposite-root T8 mandatory family; other mandatory labels retain their saved values. All four varying labels are unbounded finite integers in their stated domains. No whole-incidence exclusion or compact H9 example is claimed.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
 out=ROOT/f'results/{TAG}_review.json'
 if args.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['residual_cases'])
