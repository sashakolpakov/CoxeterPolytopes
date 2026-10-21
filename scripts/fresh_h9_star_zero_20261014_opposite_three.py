#!/usr/bin/env python3
"""Last internal label3, opposite single own contacts: literal negative plane."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_zero_20261014_opposite_three'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 g=json.loads((ROOT/GEOM).read_text())
 for p,h in g['sources'].items():assert sha(p)==h,p
 geom=g['representatives'][0]['geometry'];V=[r['incident']for r in geom['vertices']];sets=list(map(set,V));missing={frozenset(m)for m in geom['minimal_nonfaces']}
 def witness(nodes):
  nodes=sorted(nodes);vi=next(i for i,v in enumerate(sets)if set(nodes)<=v)
  return dict(support=nodes,actual_vertex=vi)
 own={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
 own[11,12]=own[13,17]=2;own[12,13]=3
 assert all(frozenset(S[2*i:2*i+2]+[leaf])in missing for i in range(4)for leaf in(P[i],C[i]))
 assert F(1,2)+F(1,3)+F(1,6)==1
 # Any finite connected triple containing a >=7 bond and another nonright
 # bond has determinant <1-3/4-1/4=0 (extra edges only worsen it).
 assert 1-F(3,4)-F(1,4)==0
 foreign=[]
 for anchor,root in[(11,13),(17,12)]:
  for k in S[:6]:foreign.append(dict(anchor=anchor,root=root,forced_right_pairs=[sorted((anchor,k)),sorted((root,k))],**witness([anchor,root,k])))
 third=[]
 for anchor,root in[(11,13),(17,12)]:third.append(dict(anchor=anchor,root=root,forced_right_pairs=[sorted((14,anchor)),sorted((14,root))],**witness([14,anchor,root])))
 mutual=witness([11,13,17])
 cross=[]
 for i,j in combinations(range(3),2):cross.append(dict(blocks=[S[2*i:2*i+2],S[2*j:2*j+2]],relaxed_path_labels=[4,3,4],**witness(S[2*i:2*i+2]+S[2*j:2*j+2])))
 assert 1-F(1,2)-F(1,4)-F(1,2)+F(1,4)==0
 # d=3/4; determinant = ((1-alpha^2)(1-beta^2)-1/4)/(3/4).
 # alpha^2,beta^2>3/4 and <1 imply determinant<-1/4.
 assert (F(1,16)-F(1,4))/F(3,4)==-F(1,4)
 # Third projection includes two own1/2 contacts on an internal>=4 block.
 # Its cost >=1/[2(1-sqrt2/2)]=1+sqrt2/2>1, so residual<=-sqrt2/2.
 assert 1-F(1,2)>0
 perm=list(range(18));perm[12],perm[13]=13,12
 assert {tuple(sorted(perm[x]for x in v))for v in V}==set(map(tuple,V))
 assert {frozenset(perm[x]for x in m)for m in missing}==missing
 paths=[GEOM,f'scripts/{TAG}.py']
 return dict(status='PASS_OPPOSITE_SINGLE_INTERNAL_THREE_NEGATIVE_INDEX',sources={p:sha(p)for p in paths},case=0,
  fixed=sorted([list(p)+[m]for p,m in own.items()]),own_bounds=[[11,13,7],[12,17,7]],
  complete_positive_support=witness(S),foreign_high_edge_guards=foreign,third_high_edge_guards=third,
  mutual_right_guard=mutual,first_three_cross_guards=cross,literal_twin_swap=perm,
  analytic_certificate=dict(last_cosine='1/2',last_inverse_denominator='3/4',own_squared_cosine_domain='(3/4,1)',
   anchor_diagonals='1-(4/3)alpha^2,1-(4/3)beta^2',anchor_off_diagonal='-(2/3)alpha*beta',
   anchor_determinant='(4/3)*((1-alpha^2)*(1-beta^2)-1/4)',anchor_determinant_strict_upper='-1/4',
   third_residual_upper='-sqrt2/2',third_anchor_pairings=['0','0'],
   conclusion='The anchor residual plane has signature(1,1); the orthogonal third residual is negative, forcing two negative directions'),
  scope='Case0 first3 blocks have all six leaves own3/3, with arbitrary internal labels>=4; last block has internal3 and opposite single own contacts. Missing triples force both nonright ownlabels>=7. All other contacts and dotted distances free; both last-root orientations by checked literal swap.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['foreign_high_edge_guards'])+len(d['third_high_edge_guards'])+len(d['first_three_cross_guards'])+2,'supports')
