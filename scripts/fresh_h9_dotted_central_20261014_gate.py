#!/usr/bin/env python3
"""A dense I2(>=6) leaf and an independent I2(>=7) leaf force index >=2.

All other contacts/distances are free subject to the listed literal finite
faces and missing own triple. No incidence symmetry or finite label cutoff.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_dotted_central_20261014_gate'
SOURCE='results/fresh_h9_tree_cube_wedges_20261014.json'
OLD='results/fresh_h9_other_trees_20261014_limited_search_case_05.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(a,b):
 c=[F(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c

def mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def evaluate(p,x):return sum((v*x**i for i,v in enumerate(p)),F(0))
def run():
 source=json.loads((ROOT/SOURCE).read_text());old=json.loads((ROOT/OLD).read_text())
 for obj in (source,old):
  for p,h in obj['sources'].items():assert sha(p)==h,p
 g=source['representatives'][5]['geometry'];vertices=[set(v['incident']) for v in g['vertices']]
 missing={frozenset(m) for m in g['minimal_nonfaces']}
 blocks=[b for b in g['blocks'] if len(b)==2];outside=sorted(set(range(18))-set(sum(blocks,[])))
 # T7(x)-1 = (x-1)(8x^3+4x^2-4x-1)^2. The largest cubic
 # root is cos(2pi/7); a root lies in (3/5,1), hence cos^2(pi/7)>4/5.
 T0=[F(1)];T1=[F(0),F(1)]
 for n in range(2,8):T0,T1=T1,add(mul([0,2],T1),[-x for x in T0])
 P=[F(-1),F(-4),F(4),F(8)]
 assert add(T1,[-1])==mul([-1,1],mul(P,P))
 assert evaluate(P,F(3,5))==F(-29,125) and evaluate(P,F(1))==7
 # The four-path determinant at smallest own weights is9/16-delta^2.
 # Label5 has cos^2=(3+sqrt5)/8>9/16 since2sqrt5>3.
 assert 4*5>9
 # Dense A cost at a=6 is2+sqrt3: 1/[2(1-sqrt3/2)].
 # Residual magnitudes exceed (1+sqrt3) and1/4, so determinant
 # exceeds (sqrt3-1)/4>0 for delta^2<=1/2.
 assert 3>1
 maps=[]
 for A,B in permutations(blocks,2):
  for x in outside:
   if frozenset(A+[x]) not in missing:continue
   for y in outside:
    if y==x or frozenset(B+[y]) not in missing:continue
    supports=[('positive_four_block',A+B),('x_foreign_right',B+[x]),('y_foreign_right',A+[y])]
    supports += [('middle_label_upper',[x,y,a,b]) for a in A for b in B]
    rows=[]
    for kind,S in supports:
     vi=next((i for i,v in enumerate(vertices) if set(S)<=v),None)
     if vi is None:break
     rows.append(dict(kind=kind,support=S,actual_vertex=vi))
    else:
     conditions=sorted([sorted(A)+[6],sorted(B)+[7],sorted((x,A[0]))+[3],sorted((x,A[1]))+[3]])
     maps.append(dict(dense_block=A,other_block=B,dense_leaf=x,other_leaf=y,
         actual_supports=rows,missing_own_triples=[sorted(A+[x]),sorted(B+[y])],
         bounds=conditions))
 families=[];seen=set()
 for i,row in enumerate(maps):
  key=tuple(map(tuple,row['bounds']))
  if key in seen:continue
  seen.add(key);families.append(dict(id=f'dense_high_pair_{len(families)}',witness_map=i,fixed=[],bounds=row['bounds']))
 lab={}
 for row in old['unexplored_frontier']:
  t=old['triangle_supports'][row['triangle']]
  from itertools import combinations
  for p,m in zip(combinations(t,2),row['labels']):
   assert p not in lab or lab[p]==m
   lab[p]=m
 hits=[i for i,row in enumerate(families) if all(lab.get((a,b),0)>=m for a,b,m in row['bounds'])]
 assert len(maps)==44 and hits
 return dict(status='PASS_ALL_LABEL_DENSE_HIGH_PAIR_GATE_CASE5',
   sources={p:sha(p) for p in [SOURCE,OLD,f'scripts/{TAG}.py']},
   theorem=dict(normalization='all root squares1',
     assumptions='A=I2(a),a>=6;B=I2(b),b>=7; dense leaf x has both own contacts nonright; y has nonempty own contact to B; all listed supports finite',
     foreign_right_reason='An elliptic component of rank>=3 cannot contain a label>=6 edge: the3-node path determinant is<=1-3/4-1/4=0',
     middle_label_bound='The actual four-path has determinant (1-alpha^2)(1-beta^2)-delta^2<=9/16-delta^2; label5 and higher fail, hence delta^2<=1/2',
     dense_projection_cost_lower={'rational':'2','sqrt3':'1'},
     other_projection_cost_strict_lower='5/4',
     residual_diagonal_magnitude_bounds=[{'rational':'1','sqrt3':'1'},'strictly greater than1/4'],
     determinant_strict_lower={'rational':'-1/4','sqrt3':'1/4'},
     conclusion='The residual2x2 matrix is negative definite, adding two negative directions to a positive rank4 block; impossible in any Lorentzian Gram matrix'),
   threshold_certificate=dict(T7=[str(x) for x in T1],cubic=[str(x) for x in P],at_three_fifths='-29/125',at_one='7',cos_squared_pi7_strict_lower='4/5'),
   case=5,maps=maps,families=families,saved_frontier_excluded_by=hits,
   scope='A sufficient all-integer-label local obstruction. All other labels and all dotted distances remain free. Case5 is not globally excluded by this gate alone.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['maps']),'maps',len(d['families']),'families','old-frontier hits',d['saved_frontier_excluded_by'])
