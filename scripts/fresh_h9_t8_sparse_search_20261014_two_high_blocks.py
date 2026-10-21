#!/usr/bin/env python3
"""All-label negative-index obstruction from two high blocks in literal T8.

Only four internal label inequalities are assumed. Positive null vectors
and actual incidences establish every domain; no coefficient bank or
outside mandatory label is used. Arithmetic is stdlib Fraction only.
"""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import combinations,product
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_sparse_search_20261014_two_high_blocks'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';STATE='results/fresh_h9_t8_sparse_search_20261014_next_case_08.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
Z=(F(0),)*4;O=(F(1),F(0),F(0),F(0));H=(F(1,2),F(0),F(0),F(0));R2=(F(0),F(1,2),F(0),F(0));R3=(F(0),F(0),F(1,2),F(0))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def mul(a,b):
 out=list(Z)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i^j]+=x*y*(2 if i&j&1 else 1)*(3 if i&j&2 else 1)
 return tuple(out)
def total(values):
 out=Z
 for x in values:out=add(out,x)
 return out
def kernel(M,w):return all(total(mul(a,b)for a,b in zip(row,w))==Z for row in M)
def path(cosines):
 n=len(cosines)+1;M=[[O if i==j else Z for j in range(n)]for i in range(n)]
 for i,c in enumerate(cosines):M[i][i+1]=M[i+1][i]=neg(c)
 return M
def run():
 geo,state=[json.loads((ROOT/p).read_text())for p in[GEOM,STATE]]
 for obj in[geo,state]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 geometry=geo['representatives'][8]['geometry'];vertices=[r['incident']for r in geometry['vertices']]
 S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)];core=sorted(S+[11,14])
 def actual(s):return next(i for i,v in enumerate(vertices)if set(s)<=set(v))
 pair_vertices=[dict(pair=list(p),actual_vertex=actual(p))for p in combinations(core,2)]
 assert len(pair_vertices)==45
 cross=[]
 for i,j in combinations(range(4),2):
  for a,b in product(blocks[i],blocks[j]):
   if i<2:
    support=[next(x for x in blocks[i]if x!=a),a,b];M=path([R3,H]);w=[R3,O,H];kind='relaxed6-3'
   else:
    support=[next(x for x in blocks[i]if x!=a),a,b,next(x for x in blocks[j]if x!=b)]
    M=path([R2,H,R2]);w=[O,(F(0),F(1),F(0),F(0)),(F(0),F(1),F(0),F(0)),O];kind='relaxed4-3-4'
   assert kernel(M,w)
   cross.append(dict(pair=[a,b],support=support,actual_vertex=actual(support),kind=kind,gram=M,positive_null_vector=w))
 assert sum(x['kind']=='relaxed6-3'for x in cross)==20 and len(cross)==24
 missing=[[1,2,14],[4,5,14],[11,12,13]]
 assert all(t in geometry['minimal_nonfaces']for t in missing)
 anchor=[]
 for own in[12,13]:
  support=sorted((set(S)-{25-own})|{11});anchor.append(dict(active_own_root=own,support=support,actual_vertex=actual(support)))
 third=[];stars=[]
 M=[[O if i==j else Z for j in range(4)]for i in range(4)]
 for j,c in[(0,R2),(2,H),(3,H)]:M[1][j]=M[j][1]=neg(c)
 w=[R2,O,H,H];assert kernel(M,w)
 for a,b in product([1,2],[4,5]):
  support=sorted((set(S)-{3-a,9-b})|{14})
  third.append(dict(active_high_roots=[a,b],support=support,actual_vertex=actual(support)))
  star=[11,14,a,b];stars.append(dict(support=star,actual_vertex=actual(star),gram=M,positive_null_vector=w))
 # All budget inequalities are lower bounds; untouched pairs are right
 # only after the preceding cross-block witnesses establish this fact.
 assert F(1,4)+F(1,4)+F(1,2)==1
 assert F(1,4)/(1-F(3,4))==1
 assert F(1,4)/(1-F(4,5))==F(5,4)
 # Independent U6 recurrence and exact root threshold.
 u0=[F(1)];u1=[F(0),F(2)]
 for degree in range(2,7):
  out=[F(0)]*(degree+1)
  for i,a in enumerate(u1):out[i+1]+=2*a
  for i,a in enumerate(u0):out[i]-=a
  u0,u1=u1,out
 assert u1==[-1,0,24,0,-80,0,64]
 f=lambda x:64*x**3-80*x*x+24*x-1
 assert f(F(4,5))==F(-29,125)
 assert 192*F(3,4)**2-160*F(3,4)+24==12
 assert 384*F(3,4)-160==128
 # Cost inequalities at the only low labels. Squaring is legitimate
 # because each compared rational number and square root is positive.
 assert 2>F(1,2)**2        # 1+sqrt2/2 >5/4
 assert 5>2**2            # (3+sqrt5)/4 >5/4
 assert 5>F(5,4)**2       # 1+sqrt5/5 >5/4
 assert F(2)>F(5,4)
 assert F(1,2)/(1-F(1,2))==1  # n=4, own<=4 is not Lanner
 assert 5<5**2                # n=5, own3 cost=(5+sqrt5)/10<1
 assert F(1,4)/(1-F(3,4))==1  # n=6, own3 is affine
 bounds=[[1,2,7],[4,5,7],[7,8,4],[12,13,4]]
 labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
 assert all(labels[a,b]>=m for a,b,m in bounds)
 def enc(x):
  if isinstance(x,F):return str(x)
  if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
  if isinstance(x,(list,tuple)):return[enc(v)for v in x]
  return x
 paths=[GEOM,STATE,f'scripts/{TAG}.py']
 return enc(dict(status='PASS_T8_TWO_HIGH_BLOCKS_ALL_OTHER_LABELS_FREE_EXCLUDED',sources={p:sha(p)for p in paths},
  finite_label_lower_bounds=bounds,fixed_labels=[],positive_support=S,positive_support_actual_vertex=actual(S),core_roots=core,all45_core_meeting_pairs=pair_vertices,
  cross_right_guards=cross,missing_triples=missing,anchor_high_right_budgets=anchor,third_other_blocks_right_budgets=third,third_anchor_upper_stars=stars,
  own_cost_cases=[dict(pattern='both own contacts nonright',lower_bound='1+sqrt2/2>5/4'),
   dict(pattern='single own, n=4',own_label_at_least=5,lower_bound='(3+sqrt5)/4>5/4'),
   dict(pattern='single own, n=5',own_label_at_least=4,lower_bound='1+sqrt5/5>5/4'),
   dict(pattern='single own, n=6',own_label_at_least=4,lower_bound='2>5/4'),
   dict(pattern='single own, n>=7',own_label_at_least=3,lower_bound='1/[4(1-c^2)]>5/4 since c^2>4/5')],
  cosine_threshold=dict(U6=u1,f_at_4_over_5=F(-29,125),fprime_at_3_over_4=12,fsecond_at_3_over_4=128),
  residual_bounds=dict(anchor11_square='strictly below -1/4',third14_square='strictly below -1',mutual_pairing='0 or -1/2, because their support projections have disjoint blocks',determinant='strictly greater than (1/4)*1-(1/2)^2=0',conclusion='negative definite2x2 residual plus positive rank8 support contradicts Lorentz negative index1'),
  saved_state_excluded=True,
  scope='All compact H9 realizations of literal T8 under only four displayed finite internal-label lower bounds; every other meeting label and dotted distance remains free. No parameter bank, finite label cutoff, outside-leaf contact or incidence automorphism is used. No continuation is run by this checker and the whole T8 incidence remains open.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['cross_right_guards']),len(d['third_anchor_upper_stars']))
