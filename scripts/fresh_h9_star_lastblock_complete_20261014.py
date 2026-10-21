#!/usr/bin/env python3
"""Logical all-label composition: first3all-own3, last internal>=4."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_lastblock_complete_20261014'
OLD='results/fresh_h9_star_shared_high_20261014_gates.json';BASE='results/fresh_h9_star_zero_20261014_gates.json';GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
NAMES=['fresh_h9_star_lastblock_20261014','fresh_h9_star_shared_high_20261014','fresh_h9_star_opposite_high_20261014']
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,h in d.get('sources',{}).items():
    assert sha(q)==h,q
    if q.endswith('.json'):read(q)
  return d
 paths=[OLD,BASE,GEOM]+[f'results/{n}{s}.json'for n in NAMES for s in('','_review')]
 old,base,geo=[read(p)for p in paths[:3]];proofs=[]
 for n in NAMES:
  p=f'results/{n}.json';r=f'results/{n}_review.json';d=read(p);review=read(r)
  assert d['status'].startswith('PASS_')and review['status'].startswith('PASS_')
  assert review['sources'][p]==sha(p)
  proofs.append(d)
 finite,shared,opposite=proofs
 assert finite['status']=='PASS_COMPLETE_LAST_BLOCK_4_5_OWN_PROFILE_EXCLUSION'
 assert shared['status']=='PASS_SHARED_SINGLE_OWN_ROOT_ALL_HIGH_LABELS_EXCLUDED'
 assert opposite['status']=='PASS_OPPOSITE_SINGLE_PI3_ALL_HIGH_LABELS_EXCLUDED'
 common={tuple(sorted((lf,r))):3 for i in range(3)for lf in(P[i],C[i])for r in S[2*i:2*i+2]}
 commonrows=sorted([list(p)+[m]for p,m in common.items()])
 assert finite['fixed_first_three_own']==commonrows and finite['last_internal_values']==[4,5]
 firstbounds=[S[2*i:2*i+2]+[4]for i in range(3)]
 assert finite['internal_lower_bounds']==firstbounds
 for d in(shared,opposite):
  assert {tuple((a,b)):m for a,b,m in d['fixed_own_contacts']if a not in(11,12,13,17)and b not in(11,12,13,17)}==common
  assert d['internal_lower_bounds']==firstbounds+[[12,13,7]]
 missing={frozenset(m)for m in geo['representatives'][0]['geometry']['minimal_nonfaces']}
 assert all(frozenset(S[2*i:2*i+2]+[lf])in missing for i in range(4)for lf in(P[i],C[i]))
 assert F(1,3)+F(1,3)+F(1,3)==1
 lower=common|{tuple(S[2*i:2*i+2]):(4 if i<3 else 6)for i in range(4)}
 strong=[];dense=[]
 def hit(lab,origin):
  hits=[f for f in base['families']if f['origin']==origin and not f['fixed']and all(lab.get((a,b),0)>=m for a,b,m in f['bounds'])]
  assert hits
  f=hits[0];return dict(family=f['id'],bounds=f['bounds'])
 for lf in(11,17):
  for root in(12,13):strong.append(dict(leaf=lf,own_root=root,threshold=4,**hit(lower|{tuple(sorted((lf,root))):4},'strong_contact')))
 for lf in(11,17):dense.append(dict(leaf=lf,**hit(lower|{tuple(sorted((lf,r))):3 for r in(12,13)},'two_dense')))
 # With all last own labels<=3, these are complete Lanner possibilities.
 low={n:[list(pr)for pr in product((2,3),repeat=2)if F(1,n)+sum(F(1,m)for m in pr)<1]for n in(6,7)}
 assert low=={6:[[3,3]],7:[[2,3],[3,2],[3,3]]}
 # For every m>=7 the same test remains true: (2,2) always has sum>1;
 # (2,3) and(3,2) have sum5/6+1/m<1; (3,3) has sum2/3+1/m<1.
 assert F(5,6)+F(1,7)<1 and F(2,3)+F(1,7)<1
 single=[]
 for a,b in product(((2,3),(3,2)),repeat=2):
  d=shared if a==b else opposite
  wanted=common|{tuple(sorted((lf,r))):m for lf,pr in zip((11,17),(a,b))for r,m in zip((12,13),pr)}
  matches=[]
  for swap in(False,True):
   pm=d['literal_twin_swap']if swap else list(range(18))
   got={tuple(sorted((pm[x],pm[y]))):m for x,y,m in d['fixed_own_contacts']}
   if got==wanted:matches.append(swap)
  assert matches
  single.append(dict(anchor_profiles=[a,b],proof='shared'if a==b else'opposite',swap=matches[0]))
 fs=old['families'][:]
 fs.append(dict(id='first_three_own_three_last_internal_ge4',origin='complete_all_label_last_block_composition',fixed=commonrows,bounds=firstbounds+[[12,13,4]]))
 paths+=[f'scripts/{TAG}.py']
 return dict(status='PASS_STAR_ZERO_COMPLETE_LAST_INTERNAL_GE4_COMPOSITION',case=0,sources={p:sha(p)for p in paths},families=fs,fixed_first_three_own=commonrows,last_internal_lower_bound=4,strong_own_threshold_cases=strong,dense_own_cases=dense,low_own_profiles=low,remaining_single_cases=single,scope='First3blocks all-own3; last internal arbitrary finite>=4 and all last own contacts arbitrary. Other finite labels and dotted distances remain free. Last internal2/3 and changed first3own contacts remain open; this is not whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
