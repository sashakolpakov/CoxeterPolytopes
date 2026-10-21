#!/usr/bin/env python3
"""Exact Gale realization and literal old-catalogue comparison of bounded proposals."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
import numpy as np
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import clutter
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_gale_frontier_20261014'
FILES=[f'results/{TAG}_{s}.json'for s in ['walk','low4']]
OLD=['results/fresh_h9_bugaenko_surgery_20261014_incidence.json','results/fresh_h9_bugaenko_remaining_20261014_newfacet.json']
B=list(combinations(range(14),4));BI={b:i for i,b in enumerate(B)};C=list(combinations(range(14),5));IX=np.array(B)
CI=np.array([[BI[c[:i]+c[i+1:]]for i in range(5)]for c in C]);CS=np.array([1,-1,1,-1,1]);MASK=[sum(1<<i for i in c)for c in C]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def minors(V):
 V=np.array(V,dtype=np.int64);assert V.shape==(14,4)and np.max(np.abs(V))<=8192
 assert 24*8192**4<2**63
 a=V[IX];out=np.zeros(len(B),dtype=np.int64)
 for p in permutations(range(4)):
  sign=(-1)**sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))
  out+=sign*a[:,0,p[0]]*a[:,1,p[1]]*a[:,2,p[2]]*a[:,3,p[3]]
 return list(map(int,out))
def build():
 oldkeys={};old_sources={}
 for f in OLD:
  d=json.loads((ROOT/f).read_text());old_sources[f]=len(d['cases'])
  for row in d['cases']:oldkeys.setdefault(tuple(row['canonical_minimal_nonfaces']),[]).append([f,row['case']])
 records=[];seen=set();stats=Counter();inputs=[]
 for f in FILES:
  data=json.loads((ROOT/f).read_text());inputs.append(dict(path=f,counts=data['counts'],distinct_rows=len(data['rows'])))
  for idx,row in enumerate(data['rows']):
   V=row['V'];det=minors(V);assert all(det)
   circuits=[];weights=[0]*14
   for c in C:
    co=[(-1)**i*det[BI[c[:i]+c[i+1:]]]for i in range(5)]
    if all(v>0 for v in co)or all(v<0 for v in co):
     circuits.append(sum(1<<i for i in c))
     for i,v in zip(c,co):weights[i]+=abs(v)
   assert all(w>0 for w in weights)
   assert all(sum(weights[i]*V[i][j]for i in range(14))==0 for j in range(4))
   facets=sorted(16383^c for c in circuits)
   assert facets==sorted(row['vertices'])and all(v.bit_count()==9 for v in facets)
   assert all(any(v>>i&1 for v in facets)for i in range(14))
   mm=clutter(facets,14);assert list(mm)==row['missing']
   assert all(m.bit_count()in[2,4,5]for m in mm)
   assert sum(m.bit_count()==2 for m in mm)>=2
   can,order=canonical(mm,14);assert list(can)==row['canonical']and list(order)==row['permutation']
   assert can not in seen;seen.add(can)
   inv={v:i for i,v in enumerate(order)}
   canonV=sorted(sum(1<<inv[j]for j in range(14)if v>>j&1)for v in facets)
   ridges=Counter(v^(1<<i)for v in facets for i in range(14)if v>>i&1);assert set(ridges.values())=={2}
   prior=oldkeys.get(can,[]);stats['old_catalogue_duplicate'if prior else'outside_both_closed_surgery_catalogues']+=1
   records.append(dict(case=len(records),source=f,source_row=idx,gale=V,positive_kernel=weights,
    nonzero_four_minors=len(det),minimum_absolute_minor=min(map(abs,det)),determinant_first_four=det[0],
    facets=facets,minimal_nonfaces=mm,canonical_minimal_nonfaces=can,canonical_vertices=canonV,
    canonical_permutation=order,missing_size_counts=dict(Counter(m.bit_count()for m in mm)),
    old_catalogue_occurrences=prior,ridge_count=len(ridges)))
 paths=FILES+OLD+[f'scripts/{TAG}_geometry.py',f'scripts/{TAG}_walk.py',f'scripts/{TAG}_low4.py',
  'scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py','scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py']
 paths += [f'src/coxeter_polytopes/combinatorics/h9_14_record{k}.py'for k in[0,1,3,4,5,7,9,31]]
 return dict(status='PASS_EXACT_REALIZABLE_GALE_NERVES_BOUNDED_SAMPLE',sources={p:sha(p)for p in paths},inputs=inputs,
  counts=dict(stats),cases=records,old_catalogue_sizes=old_sources,
  proof='All four-row determinants are nonzero. The sum of all oriented positive five-circuit dependencies gives a strictly positive integral kernel vector w. Rescaling row i by w_i gives row sum zero; its transpose kernel has dimension10 and contains1, hence defines a rational affine14-point configuration in dimension9. Positive-circuit complements are precisely facets by Gale duality, all of order9. Every label occurs in a facet, so all14 points are vertices. Exact ridge-degree2 is also checked. This proves actual polytopality, not merely chirotope feasibility.',
  scope='Only the saved matrices from two bounded proposal walks. No exhaustiveness among Gale configurations, open missing-edge shadows, or compact H9 polytopes. Numerical SVD supplies initial proposals only; every archived realization and incidence is integer/rational exact.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(build()));p=ROOT/f'results/{TAG}_geometry.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],len(out['cases']),out['counts'])
