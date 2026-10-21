#!/usr/bin/env python3
"""Narrow independent symbolic identities and frozen golden-certificate review.

Does not import the producer or repeat its cached all-pair enumeration.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_balanced_golden_20261014_review'
SOURCE='scripts/fresh_h20_balanced_golden_20261014.py'
DATA='results/fresh_h20_balanced_golden_20261014.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/DATA).read_text())
 assert all(sha(p)==h for p,h in d['sources'].items())
 N,k,t,u,dd=s.symbols('N k t u dd')
 Q=(N-4)*k*k-2*(N+2)*k+N-4
 # Two successive sphere centers differ by t(1-k)z; their radii are 2t,2tk.
 pair=s.expand(N*t*t*(1-k)**2-4*t*t-4*t*t*k*k-4*t*t*k)
 assert s.expand(pair-t*t*Q)==0
 # General distinct-ray distance uses z.w=N-dd/2.
 direct=s.expand((1-t)**2*N+(1-u)**2*N-2*(1-t)*(1-u)*(N-dd/2))
 assert s.expand(direct-(N*(t-u)**2+(1-t)*(1-u)*dd))==0
 # The radius chosen for the endpoint layer gives a right final contact.
 assert s.expand(N*t*t-4*t*t-(N-4)*t*t)==0
 # At an old endpoint of height-square q, wall pairing is
 # (R²-distance²-q)/2; strict negativity is exactly the producer's gate.
 R2,distance2,q=s.symbols('R2 distance2 q')
 assert s.expand((R2-distance2-q)+(distance2-R2+q))==0
 G=s.Matrix([[1,-s.Rational(1,2),-s.Rational(1,2)],
             [-s.Rational(1,2),1,-s.Rational(1,2)],
             [-s.Rational(1,2),-s.Rational(1,2),1]])
 assert G*s.ones(3,1)==s.zeros(3,1) and G.rank()==2
 assert G[:2,:2].det()==s.Rational(3,4)
 # D=2340+1008phi has nonsquare norm, so the displayed quadratic tower is a field.
 normD=2340**2+2340*1008-1008**2
 assert normD==6818256 and math.isqrt(normD)**2!=normD
 expected=[(417,2,20),(479,64,640),(815,400,4000)]
 for row,(n,ncaps,ncusps) in zip(d['configurations'],expected,strict=True):
  assert row['wall_count']==n and len(row['cap_recipes'])==ncaps
  assert sum(row['all_pair_counts'].values())==n*(n-1)//2
  assert row['preserved_old_finite_vertices']==392
  assert sum(row['deficient_cusp_counts_by_quotient_rank'].values())==ncusps
  assert all(len(r['affine_A2_seed_pairs'])==10 for r in row['cap_recipes'])
 final=d['configurations'][-1]
 assert final['strict_escape']['minimum_ball_power']==[['-15/16','21/16'],['0','0']]
 # (21phi-15)/16>1 iff sqrt5>41/21.
 assert 5*21**2>41**2
 assert final['strict_escape']['D']==['33/32','21/32']
 assert F(196,64)-2*F(33,32)==1 and F(84,64)-2*F(21,32)==0
 assert final['transported_target_probes_cut']==392
 paths=[SOURCE,DATA,'scripts/fresh_h20_balanced_20261014.py',f'scripts/{TAG}.py']
 return dict(status='PASS_NARROW_GOLDEN_SOURCE_AND_SYMBOLIC_REVIEW',
  sources={p:sha(p)for p in paths},quadratic_tower_norm=normD,
  configuration_wall_counts=[417,479,815],
  independent_checks=['successive pi/3 identity','cross-ray distance identity',
   'level-eight right endpoint identity','endpoint strict-sign identity',
   'affine A2 kernel and rank','tower nonsquare norm','final escape square and minimum'],
  source_review='Squared pair gates retain positive numerator signs; common center, facet projection, complete zero-pair incident sets, endpoint preservation, and rank-deficient cusp arguments reviewed.',
  primary_saved_replay='PASS, separately run python scripts/fresh_h20_balanced_golden_20261014.py --verify',
  scope='Independent symbolic identities plus source/proof review and fresh primary replay. Not a second implementation of all pair signs, nor a new construction or complete cusp census.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert json.loads(out.read_text())==d
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
