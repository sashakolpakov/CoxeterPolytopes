#!/usr/bin/env python3
"""Core-only domain replacement for the frozen16 finite-block rank cases."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_single_right_20261014_finite_core'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
OLD='results/fresh_h9_star_zero_20261014_own_three.json'
REVIEW='results/fresh_h9_star_zero_20261014_own_three_review.json'
S=[3,4,6,7,9,10,12,13];LEAVES=[11,14,17]
Z=(F(0),F(0));O=(F(1),F(0));H=(F(1,2),F(0));R=(F(0),F(1));RH=(F(0),F(1,2))
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def mul(x,y):return(x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def enc(x):return[str(t)for t in x]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen={}
 def read(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for q,h in d.get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):read(q)
  return d
 g=read(GEOM)['representatives'][0]['geometry'];V=[set(v['incident'])for v in g['vertices']]
 old=read(OLD);review=read(REVIEW)
 assert old['status']=='PASS_ALL_OWN_THREE_ARBITRARY_BLOCK_LABELS_EXCLUDED'
 assert review['status'].startswith('PASS_')and review['sources'][OLD]==sha(OLD)
 core=set(S+LEAVES);assert len(core)==11
 assert all(any(set(p)<=v for v in V)for p in combinations(sorted(core),2))
 assert any(set(S)<=v for v in V)
 co={3:H,4:RH};rows=[]
 def guard(path,labels,kernel,kind,test_pair,test_label,complete=False):
  n=len(path);assert len(set(path))==n and set(path)<=core
  vi=next(i for i,v in enumerate(V)if set(path)<=v)
  B=[[O if i==j else Z for j in range(n)]for i in range(n)]
  for i,m in enumerate(labels):B[i][i+1]=B[i+1][i]=neg(co[m])
  if complete:B[0][-1]=B[-1][0]=neg(H)
  assert len(kernel)==n and all(x[0]>=0 and x[1]>=0 and x!=Z for x in kernel)
  assert all(sumq==Z for sumq in [sum_add([mul(B[i][j],kernel[j])for j in range(n)])for i in range(n)])
  rows.append(dict(kind=kind,ordered_support=path,actual_vertex=vi,labels=labels,
    complete_triangle=complete,test_pair=sorted(test_pair),test_label=test_label,
    positive_kernel=[enc(x)for x in kernel],gram=[[enc(x)for x in r]for r in B]))
 for i,j in combinations(range(4),2):
  for a,b in product(S[2*i:2*i+2],S[2*j:2*j+2]):
   aa=next(x for x in S[2*i:2*i+2]if x!=a);bb=next(x for x in S[2*j:2*j+2]if x!=b)
   guard([aa,a,b,bb],[4,3,4],[O,R,R,O],'cross_block_right',[a,b],3)
 guard([11,12,17],[3,3],[O,O,O],'mutual_anchor_right',[11,17],3,True)
 # Both affineF4 witnesses lie entirely in the rank core; this replaces
 # the old upper witness that used an unused positive leaf.
 kernel=[(F(2),F(0)),(F(4),F(0)),(F(0),F(3)),(F(0),F(2)),R]
 guard([3,14,11,12,17],[3,4,3,3],kernel,'third_contact_upper',[11,14],4)
 guard([3,14,17,12,11],[3,4,3,3],kernel,'third_contact_upper',[14,17],4)
 assert len(rows)==27
 budgets=[]
 for leaf,own in [(14,[3,4]),(11,[12,13]),(17,[12,13])]:
  support=sorted((set(S)-{own[0]})|{leaf});vi=next(i for i,v in enumerate(V)if set(support)<=v)
  assert len(support)==8 and set(support)<=core
  budgets.append(dict(leaf=leaf,removed_own_root=own[0],retained_own_root=own[1],
    support=support,actual_vertex=vi,own_cost='1/4',foreign_inverse_diagonal_lower_bound='2'))
 assert F(1,4)+2*F(1,4)*2>1 and F(1,4)+F(1,2)*2>1
 assert old['zero_rank_cases']==[]and old['unresolved']==[]
 record_keys=[];count=0
 expected={(a,b,c,m,n)for a,b,c in product(range(7),repeat=3)for m,n in product((2,3),repeat=2)}
 for record in old['records']:
  block=tuple(record['block_labels']);assert block in set(product((4,5),repeat=4));record_keys.append(block)
  keys=set()
  for row in record['cases']:
   key=tuple(row['patterns'])+tuple(row['contact_labels']);assert key not in keys;keys.add(key)
   assert any(F(x)!=0 for x in row['determinant']);count+=1
  assert keys==expected
 assert sorted(record_keys)==list(product((4,5),repeat=4))and count==21952
 fixed=sorted([sorted((leaf,r))+[3]for leaf,roots in [(14,[3,4]),(11,[12,13]),(17,[12,13])]for r in roots])
 families=[dict(id='finite_core_'+''.join(map(str,b)),origin='core_only_finite_block_rank',
   fixed=sorted(fixed+[S[2*i:2*i+2]+[b[i]]for i in range(4)]),bounds=[])for b in product((4,5),repeat=4)]
 sources=[GEOM,OLD,REVIEW,'scripts/fresh_h9_star_zero_20261014_own_three.py',f'scripts/{TAG}.py']
 return dict(status='PASS_CORE_ONLY_FINITE_INTERNAL_LABEL_FAMILY',case=0,
  sources={p:sha(p)for p in sources},rank_core=sorted(core),core_own_equalities=fixed,
  internal_label_choices=[4,5],positive_kernel_guards=rows,foreign_budgets=budgets,
  inherited_rank_cases=count,inherited_zero_cases=0,families=families,
  proof=['Every cross-block contact is right: a trial label3 contains the actual affine4-3-4 path. Larger labels or other acute contacts preserve its positive-kernel obstruction.',
   'The anchor mutual contact is right by an actual affineA2 triangle. Both anchor-third labels are at most3 by the new literal affineF4 five-node paths, using no noncore leaf.',
   'For each core leaf, the retained-own-root support has projection cost1/4. Every foreign block inverse is entrywise nonnegative with diagonal>=2. Thus two foreign nonright contacts, or one of label>=4, would make an actual vertex nonpositive. Only zero or one foreign label3 remains.',
   'These are exactly the seven conservative projection patterns and two third-contact labels in the frozen21952-case rank bank. Every determinant is nonzero, so the eleven-core Gram has rank11, impossible in a ten-dimensional Lorentz space.',
   'No noncore own-contact value or high-internal reduction is used. The explicit internal4/5 restriction is essential to this inherited finite bank.'],
  scope='Star0,all four internal labels4or5,core leaves14/11/17 own3/3. Every noncore own contact,other finite meeting contact and dotted distance free. No whole-incidence exclusion; no claim for internal>=6 from this certificate.')
def sum_add(xs):
 ans=Z
 for x in xs:ans=add(ans,x)
 return ans
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['positive_kernel_guards']),d['inherited_rank_cases'])
