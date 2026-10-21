#!/usr/bin/env python3
"""All own contacts3/3; arbitrary finite internal labels reduce to4/5."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
from sympy import QQ,sqrt
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_zero_20261014_own_three'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
GATES='results/fresh_h9_star_zero_20261014_gates.json'
FIRST='results/fresh_h9_star_cube_metric_20261014.json'
REVIEW='results/fresh_h9_star_cube_metric_20261014_review.json'
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
K=QQ.algebraic_field(sqrt(2),sqrt(5));z=K.zero;o=K.one;h=K.convert(QQ(1,2));r2=K.from_sympy(sqrt(2))
assert K.mod.to_list()==[1,0,-14,0,9]
co={2:z,3:h,4:r2*h,5:K.from_sympy((1+sqrt(5))/4)}
def enc(x):return[str(t)for t in x.to_list()]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def sign2(a,b):
 sg=lambda x:(x>0)-(x<0)
 if not b:return sg(a)
 if not a:return sg(b)
 if sg(a)==sg(b):return sg(a)
 return sg(a)*sg(a*a-2*b*b)
def sign(x):
 p=[F(str(t))for t in x.to_list()];p=[F(0)]*(4-len(p))+p
 a3,a2,a1,a0=p;a=a0+7*a2;b=a1+17*a3;c=a1+11*a3;d=2*a2
 su=sign2(a,b);sv=sign2(c,d)
 if not sv:return su
 if not su:return sv
 if su==sv:return su
 return su*sign2(a*a+2*b*b-5*c*c-10*d*d,2*a*b-10*c*d)
def det3(a,b,c,d,e,f):return a*(d*f-e*e)-b*(b*f-c*e)+c*(b*e-c*d)
def run():
 paths=[GEOM,GATES,FIRST,REVIEW];seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,hh in d.get('sources',{}).items():
    assert sha(q)==hh,q
    if q.endswith('.json'):read(q)
  return d
 geo,gates,first,review=[read(p)for p in paths]
 V=[set(v['incident'])for v in geo['representatives'][0]['geometry']['vertices']]
 missing={frozenset(m)for m in geo['representatives'][0]['geometry']['minimal_nonfaces']}
 assert first['fixed_positive_support']==S and first['positive_leaves']==P and first['cut_leaves']==C
 assert review['status']=='PASS_INDEPENDENT_FIRST_STAR_DOMAIN_AND_RANK_REPLAY'
 own={tuple(sorted((leaf,a))):3 for i in range(4)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
 assert all(frozenset(S[2*i:2*i+2]+[leaf])in missing for i in range(4)for leaf in(P[i],C[i]))
 high=[]
 for i in range(4):
  lab=own|{tuple(S[2*j:2*j+2]):(6 if i==j else 4)for j in range(4)}
  hits=[f['id']for f in gates['families']if f['origin']=='two_dense'and all(lab.get((a,b),0)>=m for a,b,m in f['bounds'])]
  assert hits;high.append(dict(high_block=S[2*i:2*i+2],reviewed_gate=hits[0]))
 # Frozen first-state witnesses stay nonpositive under increased internal
 # cosines: every normalized off-diagonal entry only decreases.
 for row in first['forced_right_witnesses']+first['third_contact_upper_witnesses']:
  assert set(row['pair'])<=V[row['actual_vertex']]
 budget=[]
 for i in range(4):
  for leaf in(P[i],C[i]):
   support=sorted((set(S)-{S[2*i]})|{leaf});vi=next(k for k,v in enumerate(V)if set(support)<=v)
   budget.append(dict(leaf=leaf,support=support,actual_vertex=vi))
 assert F(1,4)+2*F(1,4)*2>1 and F(1,4)+F(1,2)*2>1
 patterns={}
 for group in (0,3):
  pp=[]
  for extra in [None]+[j for j in range(8)if j//2!=group]:
   v=[z]*8;v[2*group]=v[2*group+1]=-h
   if extra is not None:v[extra]=-h
   pp.append((extra,v))
  patterns[group]=pp
 def vertex_failure(labels):
  for vi,vs in enumerate(V):
   nodes=sorted(vs);B=[[o if a==b else -co[labels.get(tuple(sorted((a,b))),2)]for b in nodes]for a in nodes]
   for k in range(9):
    pivot=B[k][k]
    if sign(pivot)<=0:return dict(actual_vertex=vi,prefix_order=k+1,last_pivot=enc(pivot))
    for a in range(k+1,9):
     for b in range(k+1,9):B[a][b]-=B[a][k]*B[k][b]/pivot
  return None
 records=[];zero_cases=[];unresolved=[];n=0
 for blocks in product((4,5),repeat=4):
  inv=[o/(o-co[m]*co[m])for m in blocks]
  assert all(sign(x-K.convert(2))>=0 for x in inv)
  def pdot(x,y):
   return sum((inv[i]*(x[2*i]*y[2*i]+x[2*i+1]*y[2*i+1]+co[blocks[i]]*(x[2*i]*y[2*i+1]+x[2*i+1]*y[2*i]))for i in range(4)),z)
  labels=own|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(4)}
  ar=patterns[3];tr=patterns[0]
  aa=[o-pdot(v,v)for extra,v in ar];tt=[o-pdot(v,v)for extra,v in tr]
  ab=[[-pdot(v,w)for ex,v in ar]for ey,w in ar]
  at=[[-pdot(v,w)for ey,w in tr]for ex,v in ar]
  rows=[];anchor=[]
  for ia,ib in product(range(7),repeat=2):
   ad=aa[ia]*aa[ib]-ab[ia][ib]*ab[ia][ib];anchor.append(enc(ad));assert sign(ad)<0
   for ic in range(7):
    for ma,mb in product((2,3),repeat=2):
     ac=at[ia][ic]-co[ma];bc=at[ib][ic]-co[mb]
     d=det3(aa[ia],ab[ia][ib],ac,aa[ib],bc,tt[ic]);n+=1
     row=dict(patterns=[ia,ib,ic],contact_labels=[ma,mb],determinant=enc(d));rows.append(row)
     if d==z:
      lab=labels.copy()
      for leaf,(extra,v)in zip((P[3],C[3],C[0]),(ar[ia],ar[ib],tr[ic])):
       if extra is not None:lab[tuple(sorted((leaf,S[extra])))]=3
      lab[11,14]=ma;lab[14,17]=mb
      fail=vertex_failure(lab)
      hit=dict(block_labels=blocks,**row,actual_vertex_failure=fail);zero_cases.append(hit)
      if fail is None:unresolved.append(hit)
  records.append(dict(block_labels=blocks,anchor_determinants=anchor,cases=rows))
  print('blocks',blocks,'cases',len(rows),'zero',len(zero_cases),'unresolved',len(unresolved),flush=True)
 assert n==16*1372
 paths += [f'scripts/{TAG}.py']
 return dict(status='PASS_ALL_OWN_THREE_ARBITRARY_BLOCK_LABELS_EXCLUDED'if not unresolved else'ALL_OWN_THREE_EXACT_RANK_FRONTIER',
  sources={p:sha(p)for p in paths},case=0,all_own_contacts=sorted([list(p)+[m]for p,m in own.items()]),
  arbitrary_high_block_reduction=high,monotone_first_witness_count=30,foreign_budgets=budget,
  records=records,residual_cases=n,zero_rank_cases=zero_cases,unresolved=unresolved,
  field='Q(sqrt2,sqrt5); ANP descending rational coefficients at primitive generator sqrt2+sqrt5',
  scope='Case0,all eight leaves have both own contacts3; all four internal finite labels are arbitrary. Missing triples force internal>=4; reviewed two-dense gates exclude any internal>=6; the complete16 choices4/5 and seven conservative projections per leaf are covered. Omitted finite contacts and dotted distances remain free. This is not a whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(path.read_text())
 else:path.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['residual_cases'],len(d['zero_rank_cases']),len(d['unresolved']),flush=True)
