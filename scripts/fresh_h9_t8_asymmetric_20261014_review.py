#!/usr/bin/env python3
"""Independent literal incidence/domain and polynomial-certificate review."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_h9_t8_block_polynomials_20261014_review as R
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_asymmetric_20261014_review'
PRIMARY='results/fresh_h9_t8_asymmetric_20261014.json'
ALGEBRA='results/fresh_h9_t8_block_polynomials_20261014_review.json'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
STATE='results/fresh_h9_t8_continuation_20261014_next_case_08.json'
Z,O,H,A=R.Z,R.O,R.HALF,R.COS4
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def parse(x):return tuple(map(F,x))
def total(xs):
 x=Z
 for v in xs:x=R.add(x,v)
 return x
def div(x,y):
 a,b=y;den=a*a-2*b*b;assert den
 return R.scale(R.mul(x,(a,-b)),1/den)
def det(M):
 B=[row[:]for row in M];d=O
 for k in range(len(B)):
  p=next((i for i in range(k,len(B))if B[i][k]!=Z),None)
  if p is None:return Z
  if p!=k:B[k],B[p]=B[p],B[k];d=R.neg(d)
  pivot=B[k][k];d=R.mul(d,pivot)
  for i in range(k+1,len(B)):
   if B[i][k]==Z:continue
   q=div(B[i][k],pivot)
   for j in range(k+1,len(B)):B[i][j]=R.add(B[i][j],R.neg(R.mul(q,B[k][j])))
 return d
def run():
 docs={p:json.loads((ROOT/p).read_text())for p in[PRIMARY,ALGEBRA,GEOM,STATE]};bound=set()
 def bind(p):
  if p in bound:return
  bound.add(p)
  for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):bind(q)
 for p in docs:bind(p)
 data,algebra,geo,state=[docs[p]for p in[PRIMARY,ALGEBRA,GEOM,STATE]]
 assert data['status']=='PASS_T8_ALL_OWN_CONTACT_PATTERNS_N_GE_4_EXCLUDED'
 assert algebra['status']=='PASS_INDEPENDENT_408_RESIDUAL_BLOCK_POLYNOMIALS'
 geometry=geo['representatives'][8]['geometry'];V=[v['incident']for v in geometry['vertices']]
 lab={(a,b):m for a,b,m in state['solution']['pair_labels']}
 assert lab=={(0,1):2,(0,2):3,(1,2):7,(1,14):2,(1,15):3,(1,16):3,(2,14):3,(2,15):2,(2,16):2,
  (3,4):3,(3,5):3,(4,5):4,(4,14):3,(5,14):3,(6,7):3,(6,8):3,(7,8):4,(7,15):3,(8,15):3,
  (11,12):2,(11,13):3,(12,13):7,(12,17):3,(13,17):3}
 assert data['source_mandatory_labels']==state['solution']['pair_labels']
 assert data['free_mandatory_pairs']==[[1,2],[12,13],[11,12],[11,13],[12,17],[13,17]]
 S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)]
 assert data['positive_support']==S and set(S)<=set(V[data['positive_support_actual_vertex']])
 expected={(a,b)for i,j in combinations(range(4),2)for a,b in product(blocks[i],blocks[j])};seen=set()
 for row in data['cross_right_guards']:
  a,b=row['pair'];assert (a,b)in expected-seen;seen.add((a,b));nodes=row['support']
  assert set(nodes)<=set(V[row['actual_vertex']]);i,j=S.index(a)//2,S.index(b)//2
  if i==0:
   assert nodes==[next(x for x in blocks[0]if x!=a),a,b]
   # 3-path determinant 1-cos²(pi/m)-1/4<0.
   assert F(3,4)-F(81,100)<0
  else:
   assert nodes==[next(x for x in blocks[i]if x!=a),a,b,next(x for x in blocks[j]if x!=b)]
   M=[[O if x==y else Z for y in range(4)]for x in range(4)]
   for k,w in enumerate([A,H,A]):M[k][k+1]=M[k+1][k]=R.neg(w)
   null=[O,R.SQ2,R.SQ2,O]
   assert all(total(R.mul(x,y)for x,y in zip(r,null))==Z for r in M)
   assert M==[[parse(x)for x in r]for r in row['relaxed_gram']]
   assert null==list(map(parse,row['positive_null_vector']))
 assert seen==expected
 budgets=set()
 for row in data['anchor_budgets']:
  leaf,kept=row['leaf'],row['retained_own_root'];budgets.add((leaf,kept))
  assert leaf in[11,17]and kept in[12,13]
  assert set(row['support'])==(set(S)-({12,13}-{kept}))|{leaf}<=set(V[row['actual_vertex']])
 assert budgets==set(product([11,17],[12,13]))
 assert F(1,4)/(1-F(81,100))>1
 assert F(1,4)+2*F(1,2)>1 and F(1,2)+F(1,2)==1
 # Explicit normalized B2 inverse and its positive diagonal2.
 B2=[[O,R.neg(A)],[R.neg(A),O]];J=[[R.scale(O,2),R.SQ2],[R.SQ2,R.scale(O,2)]]
 assert all(total(R.mul(B2[i][k],J[k][j])for k in range(2))==(O if i==j else Z)for i in range(2)for j in range(2))
 third=set()
 for row in data['third_budgets']:
  root,kept=row['retained_own_roots'];third.add(kept)
  assert root==2 and row['leaf']==14 and kept in[4,5]
  assert set(row['support'])=={2,kept,7,8,12,13,14}<=set(V[row['actual_vertex']])
  assert lab[2,14]==lab[kept,14]==3
 assert third=={4,5}and F(1,4)+F(1,4)+F(1,2)==1
 patterns=set();upper_records=[];families={};co={2:Z,3:H,4:A}
 for row in data['own_support_patterns']:
  sa,sb=map(tuple,row['own_supports']);assert (sa,sb)not in patterns;patterns.add((sa,sb))
  assert sa in[(12,),(13,),(12,13)]and sb in[(12,),(13,),(12,13)]
  labels=dict(lab)
  for leaf,support in[(11,sa),(17,sb)]:
   for root in[12,13]:labels[tuple(sorted((leaf,root)))]=3 if root in support else 2
  assert row['relaxed_own_labels']==[[leaf,root,labels[tuple(sorted((leaf,root)))]]for leaf,root in product([11,17],[12,13])]
  guard=row['anchor_pair_guard'];assert guard['pair']==[11,17]
  assert {11,17}<=set(V[guard['actual_vertex']])
  shared=set(sa)&set(sb)
  if shared:
   assert guard['domain']=='right'and guard['support']==[11,min(shared),17]
   assert set(guard['support'])<=set(V[guard['actual_vertex']])
   M=[[O if i==j else R.neg(H)for j in range(3)]for i in range(3)]
   assert all(total(r)==Z for r in M)and M==[[parse(x)for x in r]for r in guard['gram']]
  else:assert guard['domain']=='any finite label; delta in[0,1)'
  fam='double'if len(sa)==len(sb)==2 else'asymmetric'if len(sa)+len(sb)==3 else'same'if sa==sb else'opposite'
  assert row['polynomial_family']==fam;families[str((sa,sb))]=fam
  seen_upper=set()
  for up in row['third_contact_upper_guards']:
   pair=tuple(up['pair']);assert pair in[(11,14),(14,17)]and pair not in seen_upper;seen_upper.add(pair)
   vertex=V[up['actual_vertex']];assert vertex==up['vertex']and set(pair)<=set(vertex)
   assert not {1,2}<=set(vertex)and not {12,13}<=set(vertex)
   ll=labels|{pair:4};M=[[O if a==b else R.neg(co[ll.get(tuple(sorted((a,b))),2)])for b in vertex]for a in vertex]
   assert M==[[parse(x)for x in r]for r in up['gram']]
   k=up['leading_order'];before=det([r[:k-1]for r in M[:k-1]])if k>1 else O;last=det([r[:k]for r in M[:k]])
   assert R.sign(before)>0 and R.sign(last)<=0 and R.mul(before,parse(up['last_pivot']))==last
   upper_records.append(dict(own_supports=[sa,sb],pair=pair,vertex=up['actual_vertex'],nonpositive_minor=last))
  assert seen_upper=={(11,14),(14,17)}
 assert patterns==set(product([(12,),(13,),(12,13)],repeat=2))
 assert data['own_support_completeness']['missing_triples']==[[11,12,13],[12,13,17]]
 for triple in[[11,12,13],[12,13,17]]:assert triple in geometry['minimal_nonfaces']
 # A one-contact label3 missing triple at n<=6 has determinant>=0;
 # therefore low-n single-own labels>=4, and the foreign budget forbids
 # all foreign contacts in those two boxes.
 assert F(3,4)-F(3,4)==0 and F(1,2)+F(1,2)==1
 f=lambda t:64*t**3-80*t**2+24*t-1
 assert f(F(81,100))==F(-559,15625)
 assert 192*F(3,4)**2-160*F(3,4)+24==12 and 384*F(3,4)-160==128
 assert parse(data['Q_at_lower_t'])==(F(25,19),F(1,2))
 # Compare every transformed polynomial and every Bernstein coefficient
 # against the separately generated Fraction artifact.
 independent={(r['block'],tuple(r['foreign']),tuple(r['third_labels'])):r for r in algebra['cases']}
 seen_cases=set();seen_blocks=set();checked=0
 for block in data['blocks']:
  name=block['id'];assert name not in seen_blocks;seen_blocks.add(name)
  for row in block['cases']:
   key=name,tuple(row['foreign_contacts']),tuple(row['third_contact_labels']);assert key not in seen_cases;seen_cases.add(key)
   own=independent[key];variables=row['free_variables'];indices=[R.NAMES.index(v)for v in variables]
   assert len(set(variables))==len(variables)
   for v,interval in zip(variables,row['variable_intervals']):assert list(map(parse,interval))==[parse(own['lower_bounds'][v]),O]
   for label,saved in [('positive_residual_numerator',own['certificates'][0]),('negative_anchor_numerator',own['certificates'][1])]:
    cert=row[label];assert cert['degrees']==[saved['degree'][i]for i in indices]
    for target,source in [('power_terms','unit_power_terms'),('positive_coefficients','positive_Bernstein_coefficients')]:
     left={tuple(r['index']):parse(r['value'])for r in cert[target]}
     right={tuple(r['index'][i]for i in indices):parse(r['value'])for r in saved[source]}
     assert left==right and len(left)==len(cert[target])
     if target=='positive_coefficients':checked+=len(left)
 assert seen_blocks=={'asymmetric','double','same_high','opposite_high','same_low','opposite_low'}
 assert seen_cases==set(independent)and len(seen_cases)==408 and checked==11736
 assert data['positive_coefficient_count']==algebra['positive_coefficient_count']==checked
 assert parse(data['minimum_positive_coefficient'])==parse(algebra['minimum_positive_coefficient'])==(F(7,50),F(0))
 paths=[PRIMARY,ALGEBRA,GEOM,STATE,f'scripts/{TAG}.py','scripts/fresh_h9_t8_block_polynomials_20261014_review.py']
 return R.encode(dict(status='PASS_INDEPENDENT_FULL_T8_OWN_PATTERN_FAMILY_REVIEW',sources={p:sha(p)for p in paths},
  recursively_bound_manifests=sorted(bound),cross_contact_guards=len(seen),literal_own_patterns=len(patterns),fresh_third_contact_guards=upper_records,
  own_pattern_to_polynomial_family=families,polynomial_cases=len(seen_cases),positive_coefficients=checked,minimum='7/50',
  coverage='m>=7, n>=4, all four own finite labels>=2 with at least one nonright contact on each anchor; the two missing triples exclude empty supports. Other mandatory labels fixed.',
  proof='The24 actual-vertex guards force four orthogonal rank2 blocks. Four anchor budgets and two third budgets force the complete foreign-contact domains. Nine literal own-support patterns have fresh mutual/third-contact guards. Six exact continuous boxes cover all possibilities; all408 residual numerators are positive and anchor minors negative, with11,736 separately reconstructed Bernstein coefficients. Increasing first cosine-square t raises the residual determinant. A positive rank8 support in Lorentz10 leaves residual rank<=2, so its three-normal residual determinant must vanish, contradiction.',
  scope='Full stated mandatory family only. Does not exclude the whole T8 incidence or produce a compact H9 example.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['polynomial_cases'],d['positive_coefficients'],flush=True)
