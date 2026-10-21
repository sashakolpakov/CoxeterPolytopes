#!/usr/bin/env python3
"""Independent Fraction biquadratic replay of the monotone cut family.

No SymPy or producer imports. Multiplication uses the square-root basis
1,sqrt2,sqrt5,sqrt10; ANP output is evaluated at sqrt2+sqrt5 only for comparison.
"""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product,permutations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014_cut_profiles'
Z=(F(0),)*4;O=(F(1),F(0),F(0),F(0));R2=(F(0),F(1),F(0),F(0));TH=(F(0),F(1),F(1),F(0))
CO={2:Z,3:(F(1,2),F(0),F(0),F(0)),4:(F(0),F(1,2),F(0),F(0)),5:(F(1,4),F(0),F(1,4),F(0))}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def sub(x,y):return add(x,neg(y))
def scale(x,a):return tuple(v*a for v in x)
def mul(x,y):
 q=[F(0)]*4
 for i,a in enumerate(x):
  if not a:continue
  for j,b in enumerate(y):
   if b:q[i^j]+=a*b*(2 if i&j&1 else 1)*(5 if i&j&2 else 1)
 return tuple(q)
def inv(x):
 y=O
 for mask in[1,2,3]:y=mul(y,tuple((-v if (i&mask).bit_count()%2 else v)for i,v in enumerate(x)))
 d=mul(x,y);assert d[0]and d[1:]==(0,0,0);return scale(y,1/d[0])
def total(xs):
 out=Z
 for x in xs:out=add(out,x)
 return out
def inverse(A):
 n=len(A);a=[list(row)+[O if i==j else Z for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  k=next(k for k in range(j,n)if a[k][j]!=Z);a[j],a[k]=a[k],a[j]
  t=inv(a[j][j]);a[j]=[mul(x,t)for x in a[j]]
  for k in range(n):
   if k!=j:
    t=a[k][j]
    if t!=Z:a[k]=[sub(x,mul(t,y))for x,y in zip(a[k],a[j])]
 return [row[n:]for row in a]
def determinant(A):
 n=len(A);out=Z
 for p in permutations(range(n)):
  q=O
  for i,j in enumerate(p):q=mul(q,A[i][j])
  if sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))%2:q=neg(q)
  out=add(out,q)
 return out
def anp(code):
 out=Z
 for x in code:out=add(mul(out,TH),(F(x),F(0),F(0),F(0)))
 return out
def enc(x):return list(map(str,x))
def run():
 src=f'results/{TAG}.json';data=json.loads((ROOT/src).read_text())
 assert data['status']=='PASS_MONOTONE_CUT_FAMILY_EXCLUDED'
 for p,h in data['sources'].items():assert sha(p)==h,p
 firstpath='results/fresh_h9_star_cube_metric_20261014.json';reviewpath='results/fresh_h9_star_cube_metric_20261014_review.json'
 first=json.loads((ROOT/firstpath).read_text());oldreview=json.loads((ROOT/reviewpath).read_text())
 assert oldreview['status']=='PASS_INDEPENDENT_FIRST_STAR_DOMAIN_AND_RANK_REPLAY'
 for obj in[first,oldreview]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 gp='results/fresh_h9_tree_cube_wedges_20261014.json';g=json.loads((ROOT/gp).read_text())['representatives'][0]['geometry']
 V=[set(v['incident'])for v in g['vertices']];S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
 assert first['fixed_positive_support']==S and first['positive_leaves']==P and first['cut_leaves']==C
 B=[[O if i==j else neg(CO[4])if i//2==j//2 else Z for j in range(8)]for i in range(8)];J=inverse(B)
 assert all(total(mul(B[i][k],J[k][j])for k in range(8))==(O if i==j else Z)for i in range(8)for j in range(8))
 assert J==[[(scale(O,2)if i==j else R2)if i//2==j//2 else Z for j in range(8)]for i in range(8)]
 guards=[]
 for i in range(4):
  for j in range(2):
   root=S[2*i+j];support=(set(S)-{S[2*i+1-j]})|{C[i]};tri={P[i],root,C[i]}
   vi=next(k for k,v in enumerate(V)if support<=v);ti=next(k for k,v in enumerate(V)if tri<=v)
   assert len(support)==8 and F(1,2)+2*F(1,4)==1
   assert 1-F(1,4)-F(3,4)==0
   guards.append(dict(cut=C[i],increased_own_root=root,eight_support=sorted(support),actual_vertex=vi,
    label6_path=sorted(tri),path_actual_vertex=ti,label6_determinant='1-1/4-3/4=0'))
 assert guards==data['own_contact_guards']
 @lru_cache(None)
 def pdot(x,y):return total(mul(x[i],mul(J[i][j],y[j]))for i in range(8)for j in range(8)if x[i]!=Z and J[i][j]!=Z and y[j]!=Z)
 def patterns(group,profile):
  ans=[]
  for extra in([None]+[j for j in range(8)if j//2!=group]if profile==(3,3)else[None]):
   row=[Z]*8;row[2*group]=neg(CO[profile[0]]);row[2*group+1]=neg(CO[profile[1]])
   if extra is not None:row[extra]=neg(CO[3])
   ans.append((extra,tuple(row)))
  return ans
 saved={}
 for row in data['records']:
  left,right=tuple(row['C0_profile']),tuple(row['C3_profile'])
  for case in row['cases']:saved[left,right,tuple(case['foreign_contacts']),tuple(case['anchor_labels'])]=anp(case['residual_determinant'])
 assert len(saved)==6300
 nonzero=0;zeros=[];distinct=set();anchors=set();profiles=list(product([3,4,5],repeat=2))
 for left,right in product(profiles,repeat=2):
  for ia,x in patterns(3,(3,3)):
   for ib,y in patterns(3,right):
    qx=sub(O,pdot(x,x));qy=sub(O,pdot(y,y));xy=neg(pdot(x,y));ad=determinant([[qx,xy],[xy,qy]])
    assert ad!=Z;anchors.add((right,ia,ib))
    for ic,z in patterns(0,left):
     qz=sub(O,pdot(z,z))
     for ma,mb in product([2,3],repeat=2):
      xz=sub(neg(CO[ma]),pdot(x,z));yz=sub(neg(CO[mb]),pdot(y,z))
      d=determinant([[qx,xy,xz],[xy,qy,yz],[xz,yz,qz]])
      assert d==saved[left,right,(ia,ib,ic),(ma,mb)]
      distinct.add(d)
      if d!=Z:nonzero+=1
      else:
       assert left==(3,3)and right==(4,4)and ia in[0,1]and ib is None and ic in[6,7]and(ma,mb)==(2,3)
       # Actual affine A3 four-cycle C0--S[ia]--P3--S[ic]--C0.
       nodes=[C[0],S[ia],P[3],S[ic]];vertex=next(k for k,v in enumerate(V)if set(nodes)<=v)
       H=[[O if i==j else neg(CO[3])if (i-j)%4 in[1,3]else Z for j in range(4)]for i in range(4)]
       assert all(total(row)==Z for row in H)and determinant(H)==Z
       zeros.append(dict(C0_profile=list(left),C3_profile=list(right),foreign_contacts=[ia,ib,ic],anchor_labels=[ma,mb],
        cycle=nodes,actual_vertex=vertex,gram=[[enc(v)for v in row]for row in H],positive_null_vector=[1,1,1,1]))
 assert nonzero==6296 and len(zeros)==4 and len(distinct)==765 and len(anchors)==105
 assert data['residual_cases']==6300 and data['profile_pairs']==81 and data['anchor_zero_cases']==0
 assert data['rank_excluded_cases']==nonzero and len(data['affine_vertex_excluded_cases'])==4
 byforeign={tuple(r['foreign_contacts']):r for r in zeros}
 for row in data['affine_vertex_excluded_cases']:
  z=byforeign[tuple(row['rank_survivor']['foreign_contacts'])]
  assert set(row['support'])==set(z['cycle'])<=V[row['actual_vertex']]
  H=[[anp(x)for x in r]for r in row['gram']]
  assert all(total(r)==Z for r in H)and row['null_vector']==[1,1,1,1]
 # All four rank-zero possibilities have now been rejected by an actual
 # spherical-vertex obstruction; rank alone would leave them unresolved.
 paths=[src,f'scripts/{TAG}.py',firstpath,reviewpath,gp,f'scripts/{TAG}_review.py']
 return dict(status='PASS_INDEPENDENT_MONOTONE_CUT_FAMILY_EXCLUSION',sources={p:sha(p)for p in paths},
  inherited_independently_checked_vertex_gates=30,new_own_contact_guards=guards,
  profile_pairs=81,exact_residual_cases=6300,nonzero_rank_cases=nonzero,
  distinct_residual_determinants=len(distinct),distinct_nondegenerate_anchor_cases=len(anchors),
  rank_zero_cases_closed_by_actual_affine_cycles=zeros,
  scope='Separate Fraction arithmetic in the four-element square-root basis, no SymPy or producer imports. Independently reconstructs all6300 residual determinants and compares every result with the primary primitive-element encoding, verifies both own-root support orientations and all8 label6 paths, and rejects all4 zero-determinant cases by literal actual affineA3 vertex minors. Inherits the30 first-state contact guards through their separate frozen review and monotonicity. Only the fixed-B2/positive-3/own-cut-nonright family is excluded; ownright cuts, changed positive contacts and changed B2 labels remain outside.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['nonzero_rank_cases'],len(d['rank_zero_cases_closed_by_actual_affine_cycles']),flush=True)
