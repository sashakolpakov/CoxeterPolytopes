#!/usr/bin/env python3
"""Independent integer-partition and rational-inequality review; no producer import."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='dimension_bound_vertex_components_20261015_review'
PRIMARY='results/dimension_bound_vertex_components_20261015.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/PRIMARY).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 # Unbounded integer knapsack by allowed part sizes, distinct from the producer's composition recurrence.
 best=[0]+[-10**9]*100;maxima={};checks=0
 for r in range(1,30):
  for n in range(r,101):best[n]=max(best[n],best[n-r]+r*(r-1)//2)
  for n in range(4,101):
   q,s=divmod(n,r);assert best[n]==q*r*(r-1)//2+s*(s-1)//2
   maxima[n,r]=best[n];checks+=1
 assert checks==d['exact_partition_checks']
 partitions=[1]+[0]*29
 for r in range(1,30):
  for n in range(r,30):partitions[n]+=partitions[n-r]
 thresholds={}
 for n in range(4,101):
  nik=F(4)+F(4,n-2 if n%2==0 else n-1)
  threshold=F(n*(n-1),2)*(F(5)/nik-1);thresholds[n]=threshold
  if n<=29:
   row=next(r for r in d['moments']if r['dimension']==n)
   assert F(row['nikulin_strict_upper'])==nik
   assert F(row['strict_lower_average_connected_pairs'])==threshold
   assert F(row['three_free_strict_upper_average_components'])==n-threshold
   assert row['rank_partitions']==partitions[n]
  for r in range(1,30):
   if n>=4*r+1:assert maxima[n,r]<=threshold
 for row in d['dimension_rank_thresholds']:
  n,r=row['dimension'],row['maximum_component_rank'];mn=(n+r-1)//r
  assert row['maximum_connected_pairs']==maxima[n,r]
  assert row['general_compact_excluded']==(maxima[n,r]<=thresholds[n])
  assert row['minimum_vertex_components']==mn and row['maximum_diagram_edges']==n-mn
  assert row['three_free_excluded']==(n-mn<=thresholds[n])
 for row in d['class_bounds']:
  r=row['maximum_component_rank']
  allowed=[n for n in range(4,30)if maxima[n,r]>thresholds[n]and n-(n+r-1)//r>thresholds[n]]
  assert row['three_free_bound']==max(allowed)
  assert row['general_bound']==4*r and row['bound_combined_with_Vinberg']==min(4*r,29)
 for row in d['rank_distribution']:
  n,r=row['dimension'],row['rank_cutoff'];a=maxima[n,r];b=n*(n-1)//2
  assert row['low_vertex_maximum']==a
  assert F(row['strict_lower_fraction_with_larger_component'])==(thresholds[n]-a)/(b-a)
 for row in d['odd_boundary_checks']:
  r=row['r'];assert row['first_excluded_dimension']==4*r+1
  assert row['exact_partition_maximum']==4*comb(r,2)==2*r*(r-1)
 # Universal tail: symbolic coefficient identities evaluated in a separate polynomial basis.
 # d=4r+2+t (even) gives d*t/8 >=0; d=4r+3+t (odd) gives ((4r+3+t)(1+t)+5)/8>0.
 for r in range(1,101):
  for t in range(0,9,2):
   n=4*r+2+t;assert F(n*(n-6),8)-F(n*(r-1),2)==F(n*t,8)
   n=4*r+3+t;assert F((n-1)*(n-5),8)-F(n*(r-1),2)==F((4*r+3+t)*(1+t)+5,8)>0
 assert 'every missing face to have size2' in d['theorem_scope']
 paths=[PRIMARY,'scripts/dimension_bound_vertex_components_20261015.py',f'scripts/{TAG}.py']
 return dict(status='PASS_INDEPENDENT_COMPACT_DIMENSION_PARTITION_INEQUALITIES',sources={p:sha(p)for p in paths},
  integer_knapsack_cases=checks,moment_rows=len(d['moments']),partition_counts_verified_through=29,
  class_bounds=d['class_bounds'],distribution_rows=len(d['rank_distribution']),
  scope='Independent arithmetic and partition extremum review only. The geometric Nikulin and flag quadrilateral-charge premises require the separately sourced proof. The unrestricted29/995bounds are not improved; flag means every missing face has size2, not merely no triangular2faces. No producer arithmetic imports.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['integer_knapsack_cases'])
