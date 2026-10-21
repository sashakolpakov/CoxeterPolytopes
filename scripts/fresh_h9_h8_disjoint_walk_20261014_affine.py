#!/usr/bin/env python3
"""Complete irreducible affine rank-eight one-wall extensions of the62supports."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_disjoint_walk_20261014_affine'
def run():
 path='results/fresh_h9_h8_disjoint_walk_20261014.json';d=json.loads((ROOT/path).read_text())
 for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
 out=[];counts=Counter()
 for row in d['queue']:
  z=row['profile_data'];G=[[A.Q(*x)for x in r]for r in z['normalized_Gram']];J=[[A.Q(*x)for x in r]for r in z['normalized_inverse']]
  S=row['support'];cs=z['support_components'];where={i:k for k,c in enumerate(cs)for i in c}
  if any(x.b for r in G for x in r):counts['supports_with_old_label5']+=1;continue
  # Connected Euclidean Coxeter diagrams of rank8 have only3/4edges,
  # new-node degree<=4, and are trees or the affine A8 cycle.
  found=[]
  for n in range(1,5):
   for ids in combinations(range(8),n):
    if len({where[S[i]]for i in ids})!=len(cs):continue
    if n not in[len(cs),len(cs)+1]:continue
    for labels in product([3,4],repeat=n):
     counts['exact_cost_trials']+=1;c0=A.Q();c1=A.Q()
     for k,i in enumerate(ids):c0+=J[i][i]*(1 if labels[k]==3 else 2)
     for k,l in combinations(range(n),2):
      t=2*J[ids[k]][ids[l]]
      if labels[k]==labels[l]:c0+=t*(1 if labels[k]==3 else 2)
      else:c1+=t
     if c1 or c0!=2:continue
     # Kernel is (J*c,1), where c has entries1 or sqrt2.
     kernel=[]
     for i in range(8):
      aa=sum((J[i][j]for j,m in zip(ids,labels)if m==3),A.Q())
      bb=sum((J[i][j]for j,m in zip(ids,labels)if m==4),A.Q())
      assert aa.sign()>=0 and bb.sign()>=0 and (aa or bb)
      kernel.append([aa,bb])
     found.append(dict(contacts=[[S[i],m]for i,m in zip(ids,labels)],projection_square=2,
        positive_kernel_first8_as_a_plus_b_sqrt2=kernel,last_kernel_coefficient=1))
  out.append(dict(germ=row['id'],support=S,affine_profiles=found))
  counts['profiles']+=len(found);counts['eligible_supports']+=1
 paths=[path,f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
 return A.enc(dict(status='PASS_COMPLETE_LOCAL_AFFINE_ENDPOINT_PROFILES',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  counts=dict(counts),supports=out,scope='One new wall completing the8PD edge support to a bounded Euclidean cusp link. The full9-node Gram must be irreducible affine: the8old roots are independent, so a disconnected leftover finite component would leave the link unbounded. Affine rank8 classification permits onlylabels3/4 and degree<=4; exact Schur cost2 and strict positive kernel check every accepted profile. No global old-wall feasibility claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],flush=True)
