#!/usr/bin/env python3
"""Independent algebraic-field matrix replay of the first-star rank gate.

No producer imports. Reconstruct normalized Grams and their Schur complements
using SymPy's QQ(sqrt2) domain, not the producer's scalar-pair implementation.
"""
import argparse,json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
from sympy import QQ,Matrix,Rational,sqrt
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014'
K=QQ.algebraic_field(sqrt(2));ONE=K.one;ZERO=K.zero
COS={2:ZERO,3:K(Rational(1,2)),4:K.from_sympy(sqrt(2)/2)}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dm(rows):return DomainMatrix(rows,(len(rows),len(rows[0])),K)
def field(x):return K.from_sympy(Rational(x[0])+Rational(x[1])*sqrt(2))
def enc(x):
 c=x.to_list();c=[K.dom.zero]*(2-len(c))+c
 return [str(c[1]),str(c[0])]
@lru_cache(None)
def sign(x):
 e=K.to_sympy(x)
 if e==0:return 0
 if e.is_positive is True:return 1
 assert e.is_negative is True,e
 return -1
def run():
 src=f'results/{TAG}.json';data=json.loads((ROOT/src).read_text())
 paths=[src,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
 for p,h in data['sources'].items():assert sha(p)==h,p
 gp='results/fresh_h9_tree_cube_wedges_20261014.json';sp='results/fresh_h9_tree_cube_wedges_20261014_case_00.json'
 geo=json.loads((ROOT/gp).read_text());state=json.loads((ROOT/sp).read_text())
 for obj in[geo,state]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 g=geo['representatives'][0]['geometry'];V=[set(v['incident'])for v in g['vertices']]
 labels={(i,j):m for i,j,m in state['solution']['pair_labels']}
 S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
 assert data['fixed_positive_support']==S and data['positive_leaves']==P and data['cut_leaves']==C
 def gram(nodes,extra={}):
  lab=labels|extra
  return dm([[ONE if i==j else -COS[lab.get(tuple(sorted((i,j))),2)]for j in nodes]for i in nodes])
 guards=[]
 for row in data['forced_right_witnesses']+data['third_contact_upper_witnesses']:
  pair=tuple(row['pair']);nodes=sorted(V[row['actual_vertex']]);assert set(pair)<=set(nodes)
  G=gram(nodes,{pair:row['test_label']});previous=ONE
  for k in range(1,10):
   minor=G.extract(range(k),range(k)).det()
   if sign(minor)<=0:
    assert k==row['leading_order'] and minor/previous==field(row['last_pivot']);break
   previous=minor
  else:raise AssertionError('no nonpositive vertex minor')
  guards.append(dict(pair=list(pair),label=row['test_label'],prefix=k,determinant=enc(minor)))
 assert Counter(r['label']for r in guards)=={3:28,4:2}
 expected_cross={tuple(sorted((i,j)))for a in range(4)for b in range(a+1,4)for i in S[2*a:2*a+2]for j in S[2*b:2*b+2]}
 assert {tuple(r['pair'])for r in data['forced_right_witnesses']}==expected_cross|set(zip(P,C))
 assert {tuple(r['pair'])for r in data['third_contact_upper_witnesses']}=={(11,14),(14,17)}
 B=gram(S);J=B.inv();assert B*J==DomainMatrix.eye(8,K).to_dense()
 assert all(sign(B.extract(range(k),range(k)).det())>0 for k in range(1,9))
 assert [[enc(x)for x in row]for row in B.to_list()]==data['positive_gram']
 assert [[enc(x)for x in row]for row in J.to_list()]==data['positive_inverse']
 budget=[]
 for row in data['foreign_contact_budget_witnesses']:
  leaf=row['leaf'];i=next(i for i in range(4)if leaf in[P[i],C[i]])
  support=set(S)-{S[2*i]};assert support|{leaf}==set(row['support'])<=V[row['actual_vertex']]
  ordered=sorted(support);JJ=gram(ordered).inv().to_list()
  assert all(sign(x)>=0 for rr in JJ for x in rr)
  own=ordered.index(S[2*i+1]);assert JJ[own][own]==ONE
  foreign=[k for k,v in enumerate(ordered)if v not in S[2*i:2*i+2]]
  assert len(foreign)==6 and all(JJ[k][k]==K(2)for k in foreign)
  assert Fraction(1,4)+2*Fraction(1,4)*2>1 and Fraction(1,4)+2*Fraction(1,2)>1
  budget.append(dict(leaf=leaf,support=sorted(support),foreign_diagonal=2,own_cost='1/4'))
 assert len(budget)==8
 patterns=[]
 for i in range(4):
  group=[]
  for extra in[None]+[j for j in range(8)if j//2!=i]:
   p=[ZERO]*8;p[2*i]=p[2*i+1]=-COS[3]
   if extra is not None:p[extra]=-COS[3]
   group.append(p)
  patterns.append(group)
 assert [[[enc(x)for x in p]for p in arr]for arr in patterns]==data['projection_patterns']
 # Build complete 11-normal Grams in block form. The only imposed contacts
 # outside S are P3-C3=0 and C0-anchor labels2/3. Its Schur block must have
 # rank<=2 in any10-dimensional ambient form.
 def schur(pa,pb,pc,ma,mb):
  Q=dm([patterns[3][pa],patterns[3][pb],patterns[0][pc]])
  D=dm([[ONE,ZERO,-COS[ma]],[ZERO,ONE,-COS[mb]],[-COS[ma],-COS[mb],ONE]])
  return D-Q*J*Q.transpose()
 counts=Counter();rows=[]
 saved={(tuple(r['anchor_patterns']),r['third_pattern'],tuple(r['contact_labels'])):r for r in data['rank_obstructions']}
 anchor_saved={tuple(r['patterns']):r for r in data['anchors']}
 for ia,ib in product(range(7),repeat=2):
  H=schur(ia,ib,0,2,2).extract([0,1],[0,1]);det=H.det();assert sign(det)<0
  assert det==field(anchor_saved[ia,ib]['determinant'])
  for ic in range(7):
   for ma,mb in product([2,3],repeat=2):
    R=schur(ia,ib,ic,ma,mb);det=R.det();assert sign(det)!=0
    assert det==field(saved[(ia,ib),ic,(ma,mb)]['residual_determinant'])
    # Independent whole-block determinant identity in one representative
    # per anchor checks the rank11 translation without rebuilding1372.
    if ic==0 and ma==mb==2:
     Q=dm([patterns[3][ia],patterns[3][ib],patterns[0][ic]])
     D=dm([[ONE,ZERO,ZERO],[ZERO,ONE,ZERO],[ZERO,ZERO,ONE]])
     full=DomainMatrix.vstack(DomainMatrix.hstack(B,Q.transpose()),DomainMatrix.hstack(Q,D))
     assert full.det()==B.det()*det and sign(full.det())!=0
    counts[sign(det)]+=1;rows.append(enc(det))
 assert len(rows)==1372 and dict(counts)=={-1:150,1:1222}
 assert data['determinant_sign_counts']=={str(k):v for k,v in counts.items()}
 return dict(status='PASS_INDEPENDENT_FIRST_STAR_DOMAIN_AND_RANK_REPLAY',sources={p:sha(p)for p in paths},
  vertex_minor_guards=guards,schur_budget_guards=budget,anchor_count=49,
  nonzero_rank11_cases=len(rows),determinant_sign_counts={str(k):v for k,v in counts.items()},
  residual_determinants=rows,full11x11_determinants_checked=49,
  scope='No producer imports. Independent QQ(sqrt2) domain matrices check all30 actual-vertex failure minors, all8 nonnegative inverse/contact-budget premises, the complete seven projection patterns per leaf, all49 nondegenerate anchor Grams, and all1372 nonzero Schur determinants. Also verifies49 full11x11 determinant identities. Only the first saved mandatory assignment of star case0 is excluded, for all omitted labels and dotted distances; other assignments and the star incidence remain open.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['nonzero_rank11_cases'],flush=True)
