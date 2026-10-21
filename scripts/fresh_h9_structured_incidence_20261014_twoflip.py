#!/usr/bin/env python3
"""Controlled two-flip diagnostic through one large temporary missing face.

Source-vertex carriers restrict this to the fixed product-section mechanism.
Abstract sphere outputs need a separate realization if any gate survives.
"""
import argparse,json
from itertools import product
from collections import Counter
from pathlib import Path
from fresh_h9_structured_incidence_20261014 import ROOT,sha,bits,Gates,components
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
TAG='fresh_h9_structured_incidence_20261014_twoflip';SOURCE='results/fresh_h9_structured_incidence_20261014.json'

def subsets(m):
 o=m
 while m:yield m;m=(m-1)&o
 yield 0

def altered_clutter(M,counts,removed,added):
 delta=Counter()
 for f in removed:
  for s in subsets(f):delta[s]-=1
 for f in added:
  for s in subsets(f):delta[s]+=1
 def face(s):return counts.get(s,0)+delta.get(s,0)>0
 candidates=set(M)|{s for s,d in delta.items()if d<0 and not face(s)}
 return sorted(s for s in candidates if not face(s)and all(face(s^(1<<i))for i in bits(s)))

def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h
 gates=Gates();stats=Counter();rows=[];first_rows=[];unique={};original={}
 for r in d['cases']:original[r['facet_count'],tuple(r['canonical_minimal_nonfaces'])]=r['case']
 for r in d['cases']:
  co=r['coefficients'];n=sum(map(len,co));assert n==r['facet_count']
  off=[];p=0
  for c in co:off.append(p);p+=len(c)
  carriers={}
  for v in product(*[range(len(c))for c in co]):
   A=sum(1<<(o+j)for c,o,i in zip(co,off,v)for j in range(len(c))if j!=i)
   assert A.bit_count()==10;carriers[A]=list(v)
  V={x['incident']for x in r['vertex_records']};M=r['minimal_nonfaces'];facecounts=Counter(s for v in V for s in subsets(v))
  for A in sorted(carriers):
   P=sum(1<<i for i in bits(A)if A^(1<<i)in V);N=A^P
   if not 2<=P.bit_count()<=4 or P not in M:continue
   old={v for v in V if v&N==N};required={A^(1<<i)for i in bits(P)}
   if old!=required:continue
   new={A^(1<<i)for i in bits(N)};V1=(V-old)|new;stats['first_exact_large_face_flips']+=1
   first=dict(parent=r['case'],carrier=A,old_missing=P,temporary_missing=N,failed_second_carriers=[])
   for B in sorted(carriers):
    if B==A or B&N!=N:continue
    Q=B^N
    old2={v for v in V1 if v&Q==Q};required2={B^(1<<i)for i in bits(N)}
    stats['candidate_second_carriers']+=1
    if old2!=required2:
     witness=min(old2^required2)
     first['failed_second_carriers'].append(dict(carrier=B,witness=witness,unexpected=witness in old2))
     continue
    # N is the new minimal nonface after the first exact bistellar flip.
    # The second flip removes it. Its new Q has size2..4.
    new2={B^(1<<i)for i in bits(Q)};V2=(V1-old2)|new2
    mm=altered_clutter(M,facecounts,V-V2,V2-V);stats['two_flip_recipes']+=1
    row=dict(parent=r['case'],first_carrier=A,old_missing=P,temporary_missing=N,second_carrier=B,new_missing=Q,
             source_vertices=[carriers[A],carriers[B]],removed_facets=sorted(V-V2),added_facets=sorted(V2-V))
    if max(map(int.bit_count,mm))>5:row.update(verdict='LARGE_FINAL_MINIMAL_NONFACE',witness=next(m for m in mm if m.bit_count()>5))
    else:
     stats['small_final_clutter']+=1;can,order=canonical(mm,n);key=n,tuple(can)
     if key in original:row.update(verdict='ORIGINAL_OBSTRUCTED_TYPE',original_case=original[key],permutation=list(order))
     elif key in unique:row.update(verdict='DUPLICATE_NEW_SPHERE',representative=unique[key],permutation=list(order))
     else:
      unique[key]=len(rows);row.update(minimal_nonfaces=mm,canonical_minimal_nonfaces=list(can),canonical_order=list(order),
        actual_sphere_facets=sorted(V2),product_factors=components(mm,n))
      gate=gates.check(mm)
      if gate:row.update(verdict='EXACT_INHERITED_OBSTRUCTION',gate=gate)
      else:row.update(verdict='UNRESOLVED_REALIZATION_AND_LABEL_FRONTIER')
    stats[row['verdict']]+=1;rows.append(row)
   first_rows.append(first)
  print('parent',r['case'],dict(stats),flush=True)
 paths=[SOURCE,f'scripts/{TAG}.py','scripts/fresh_h9_structured_incidence_20261014.py',
        'scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py']+gates.paths
 return dict(status='TERMINAL_EXACT_TWO_FLIP_DIAGNOSTIC',sources={p:sha(p)for p in paths},counts=dict(stats),first_flips=first_rows,records=rows,
  scope='For the16 saved realized section nerves, exhaust source-vertex-carrier first flips2/8,3/7,4/6 followed by a source-carrier inverse-size flip that removes the temporary large missing face, with a distinct carrier. Exact star equalities certify sphere moves. Newly formed abstract spheres are not claimed Euclidean-polytopal. Inherited obstructions are exact when present; all ungated outputs remain unresolved, not excluded.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();p=ROOT/f'results/{TAG}.json';d=json.loads(json.dumps(d))
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'])
