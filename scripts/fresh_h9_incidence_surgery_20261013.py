#!/usr/bin/env python3
"""Affine-incidence constraints for three previously unmutated H9 sources."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse,gzip,hashlib,json
import numpy as np
import networkx as nx
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/alexeev_parked_rank10_h9_examples.json'
OUT=ROOT/'results/fresh_h9_incidence_surgery_20261013.json'
IDS=['A14_8_1_C02','A14_10_1_C02','A16_14_1_C01']
def kernel(A):
 M=[list(r)for r in A];n=len(M);piv=[];k=0
 for j in range(n):
  found=next((i for i in range(k,n)if M[i][j]),None)
  if found is None:continue
  M[k],M[found]=M[found],M[k];v=M[k][j];M[k]=[x/v for x in M[k]]
  for i in range(n):
   if i!=k:
    v=M[i][j];M[i]=[x-v*y for x,y in zip(M[i],M[k])]
  piv.append(j);k+=1
 if k!=n-1:return None
 free=next(j for j in range(n)if j not in piv);v=[F(0)]*n;v[free]=F(1)
 for i,j in enumerate(piv):v[j]=-M[i][free]
 if all(x<0 for x in v):v=[-x for x in v]
 assert all(x>0 for x in v) and all(sum(x*y for x,y in zip(row,v))==0 for row in A)
 return v

def cover_solutions(components,edges,parallel,limit=300):
 remaining=[c for c in components if not any(i in c and j in c for i,j in parallel)]
 masks=[sum(1<<k for k,c in enumerate(remaining)if i in c and j in c)for i,j in edges]
 full=(1<<len(remaining))-1
 # Dominated columns can be removed when finding the minimum cardinality:
 # replacing one by a superset never worsens any cover.
 unique=sorted(set(masks)-{0})
 dominant=[a for a in unique if not any(a!=b and a|b==b for b in unique)]
 options=[[m for m in dominant if m>>i&1]for i in range(len(remaining))]
 calls=0
 @lru_cache(None)
 def feasible(rem,k):
  nonlocal calls;calls+=1
  if not rem:return ()
  if not k:return None
  useful=[m&rem for m in dominant]
  maximum=max(x.bit_count()for x in useful)
  if (rem.bit_count()+maximum-1)//maximum>k:return None
  i=min((i for i in range(len(remaining))if rem>>i&1),key=lambda i:len(options[i]))
  for m in sorted(options[i],key=lambda m:-(m&rem).bit_count()):
   tail=feasible(rem&~m,k-1)
   if tail is not None:return (m,)+tail
  return None
 minimum=0
 while feasible(full,minimum)is None:minimum+=1
 proof=feasible(full,minimum)
 # Exhaust finite-edge supports directly. For all-pair supports retain only
 # a concrete dominant-column minimum witness; no exhaustive support claim.
 sols=[];total=0;complete=len(edges)<30
 if complete:
  for ids in combinations(range(len(edges)),minimum):
   m=0
   for i in ids:m|=masks[i]
   if m==full:
    total+=1
    if len(sols)<limit:sols.append([list(edges[i])for i in ids])
 else:
  sols=[sorted([list(edges[masks.index(m)])for m in proof])]
 return {'remaining_circuits':len(remaining),'minimum_changes':minimum,
  'minimum_supports':sols,'total_minimum_supports':total if complete else None,
  'enumeration_complete':complete and total<=limit,'enumeration_limit':limit,
  'minimum_proof':'exact memoized set-cover with dominated-column replacement',
  'dynamic_states':calls,'dominant_column_count':len(dominant)}

def build():
 d=json.loads(SOURCE.read_text());rows=[]
 for source in d['cases']:
  if source['id']not in IDS:continue
  G=[[F(x)for x in r]for r in source['section']['restricted_mirrors']['positive_wall_gram']]
  n=len(G);assert n in[15,16]
  path=ROOT/source['artifacts']['ray_archive_path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==source['artifacts']['ray_archive_sha256']
  with gzip.open(path,'rt')as f:
   header=json.loads(next(f));rays=[json.loads(x)for x in f]
  affine={};cusps=[]
  for ray,active,q,kind in rays:
   if kind!='ideal':continue
   active=[i-1 for i in active];g=nx.Graph();g.add_nodes_from(active)
   g.add_edges_from((i,j)for i,j in combinations(active,2)if G[i][j])
   comps=[]
   for c in nx.connected_components(g):
    ids=tuple(sorted(c));A=[[G[i][j]for j in ids]for i in ids];v=kernel(A);assert v is not None
    affine[ids]=v;comps.append(list(ids))
   cusps.append({'active':active,'components':comps})
  finite=[];parallel=[];zero=[];hist={}
  for i,j in combinations(range(n),2):
   y=G[i][j]*G[i][j]/(G[i][i]*G[j][j]);assert G[i][j]<=0
   label={F(0):2,F(1,4):3,F(1,2):4,F(1):0}[y]
   hist[str(label)]=hist.get(str(label),0)+1
   (parallel if label==0 else zero if label==2 else finite).append((i,j))
  components=list(affine)
  row={'id':source['id'],'facet_count':n,'source_positive_gram':[[str(x)for x in r]for r in G],
   'ray_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'ideal_vertices':cusps,
   'affine_circuits':[{'nodes':list(c),'positive_kernel':[str(x)for x in affine[c]]}for c in sorted(affine)],
   'old_pair_histogram':hist,'parallel_pairs':[list(p)for p in parallel],
   'finite_edge_cover':cover_solutions(components,finite,parallel),
   'all_pair_cover':cover_solutions(components,finite+zero,parallel)}
  print(row['id'],'circuits',len(affine),'parallel',len(parallel),'finite gate',row['finite_edge_cover']['minimum_changes'],len(row['finite_edge_cover']['minimum_supports']),row['finite_edge_cover']['enumeration_complete'],'all gate',row['all_pair_cover']['minimum_changes'],flush=True)
  rows.append(row)
 return {'status':'EXACT_AFFINE_CIRCUIT_OBLIGATIONS_WITH_FINITE_SET_COVER_OPTIMIZATION','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'cases':rows,
 'scope':'Every preserved affine circuit obstructs compactness with retained actual facets. Old parallel pairs are allowed to vary freely; remaining circuits must have an internally changed pair. Set-cover minima use exact bitset recursion; geometry/circuit kernels are exact Fraction.'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=build()
 if a.verify:assert d==json.loads(OUT.read_text());print('PASS replay')
 else:OUT.write_text(json.dumps(d,indent=2)+'\n');print('WROTE',OUT)
if __name__=='__main__':main()
