#!/usr/bin/env python3
"""Separate integer/rational Gale-to-primal replay and literal gate audit.

No Gale producer, determinant expansion, clutter or nauty code is imported.
Global label searches reuse traversal but replace finite and modular arithmetic.
"""
import argparse,json,copy
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from flint import fmpz_mat,fmpq_mat,fmpq
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_gale_frontier_20261014'
GEOM=f'results/{TAG}_geometry.json';LABELS=f'results/{TAG}_labels.json'
PAT='results/fresh_h9_newfacet_frontier_20261014_patterns.json';REL='results/fresh_h9_newfacet_frontier_20261014_relations.json'
LOCAL='results/fresh_h9_newfacet_frontier_20261014_local_independent.json'
BLOCK='results/fresh_h9_bugaenko_remaining_20261014_block_gate_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d

def bits(m):return [i for i in range(m.bit_length())if m>>i&1]
def face_table(facets):
 flag=bytearray(1<<14)
 for f in facets:
  s=f
  while s:flag[s]=1;s=(s-1)&f
 flag[0]=1
 missing=[]
 for s in range(1,1<<14):
  if flag[s]:continue
  t=s;ok=True
  while t:
   b=t&-t;t-=b
   if not flag[s^b]:ok=False;break
  if ok:missing.append(s)
 return flag,missing

def geometry(row):
 G=row['gale'];minors={c:int(fmpz_mat([G[i]for i in c]).det())for c in combinations(range(14),4)}
 assert all(minors.values());weights=[0]*14;circuits=[]
 for c in combinations(range(14),5):
  co=[(-1)**i*minors[c[:i]+c[i+1:]]for i in range(5)]
  if all(v>0 for v in co)or all(v<0 for v in co):
   pos=[abs(v)for v in co];circuits.append((c,pos))
   for i,v in zip(c,pos):weights[i]+=v
 assert weights==row['positive_kernel']and all(v>0 for v in weights)
 assert all(sum(weights[i]*G[i][j]for i in range(14))==0 for j in range(4))
 facets=sorted(16383^sum(1<<i for i in c)for c,v in circuits)
 assert facets==row['facets'];assert all(any(f>>i&1 for f in facets)for i in range(14))
 # Explicit rational homogeneous primal coordinates. H is a transpose-kernel
 # basis whose ten free columns are the identity. Replace its last row by w,
 # then divide column i by w_i to make the first coordinate equal to one.
 B=fmpq_mat([[G[i][j]for i in range(4)]for j in range(4)])
 X=-B.inv()*fmpq_mat([[G[i][j]for i in range(4,14)]for j in range(4)])
 H=fmpq_mat(10,14)
 for j in range(10):
  for i in range(4):H[j,i]=X[i,j]
  H[j,j+4]=1
 assert H*fmpq_mat(G)==fmpq_mat(10,4)
 A=fmpq_mat([[1]*14]+[[H[j,i]/weights[i]for i in range(14)]for j in range(9)])
 assert A.rank()==10
 assert A*fmpq_mat([[weights[i]*v for v in G[i]]for i in range(14)])==fmpq_mat(10,4)
 for c,co in circuits:
  values=[0]*14
  for i,v in zip(c,co):values[i]=v
  first=fmpq(values[13],weights[13]);ell=[first]+[values[j+4]-first*weights[j+4]for j in range(9)]
  actual=fmpq_mat([ell])*A
  assert all(actual[0,i]==fmpq(values[i],weights[i])for i in range(14))
  f=[i for i in range(14)if i not in c]
  assert fmpq_mat([[A[j,i]for i in f]for j in range(10)]).rank()==9
  assert all(actual[0,i]>0 for i in c)and all(actual[0,i]==0 for i in f)
 flag,missing=face_table(facets);assert missing==row['minimal_nonfaces']
 perm=row['canonical_permutation'];assert sorted(perm)==list(range(14));inv={v:i for i,v in enumerate(perm)}
 def trans(m):return sum(1<<inv[i]for i in bits(m))
 assert sorted(map(trans,missing))==row['canonical_minimal_nonfaces']
 assert sorted(map(trans,facets))==row['canonical_vertices']
 assert min(map(abs,minors.values()))==row['minimum_absolute_minor']
 ridge=Counter(f^(1<<i)for f in facets for i in bits(f));assert set(ridge.values())=={2}
 return flag,len(facets),len(ridge)

def run():
 d=load(GEOM);labels=load(LABELS);pat=load(PAT);rel=load(REL);local=load(LOCAL);load(BLOCK)
 assert local['obstruction_count']==44 and local['six_node_relation_rows']==96
 P={r['pattern']:r for r in pat['patterns']if r['verdict']=='EXACT_LOCAL_OBSTRUCTION'};P['seven']=rel['seven_node_obstruction']
 counts=Counter();rows=[];allkeys=set();searches=[]
 for row,gate in zip(d['cases'],labels['cases']):
  assert row['case']==gate['case']==len(rows)
  flag,nfac,nridge=geometry(row);key=tuple(row['canonical_minimal_nonfaces']);assert key not in allkeys;allkeys.add(key)
  if gate['gate']=='local':
   p=P[gate['pattern']];perm=gate['permutation'];assert len(perm)==len(set(perm))==p['order'];union=sum(1<<i for i in perm)
   lift=lambda m:sum(1<<perm[i]for i in bits(m))
   assert {m for m in row['minimal_nonfaces']if m&union==m}=={lift(m)for m in p['minimal_nonfaces']}
   for f in range(1<<p['order']):
    if all(f&m!=m for m in p['minimal_nonfaces']):assert flag[lift(f)];counts['local_actual_face_checks']+=1
  elif gate['gate']=='two_block':
   a,b=gate['blocks'];assert not a&b and a.bit_count()in[4,5]and b.bit_count()in[4,5]
   assert a in row['minimal_nonfaces']and b in row['minimal_nonfaces']
   for i in bits(a):
    for j in bits(b):assert flag[(a|b)^(1<<i)^(1<<j)];counts['block_actual_face_checks']+=1
  else:assert gate['gate']=='complete_label_search';searches.append((row,gate['search']))
  counts['four_by_four_determinants']+=1001;counts['rational_primal_facets']+=nfac
  counts[gate['gate']]+=1;rows.append(dict(case=row['case'],facets=nfac,ridges=nridge,gate=gate['gate']))
  if len(rows)%100==0:print('independent rational primal',len(rows),flush=True)
 # Separate finite-tree and finite-field arithmetic; search control flow shared.
 import fresh_h9_bugaenko_closure_20261014 as C
 import fresh_h9_newfacet_frontier_20261014_engine as E
 from fresh_h9_bugaenko_closure_20261014_large_first_audit import finite,determinant
 C.finite=E.finite=finite;C.residue=determinant
 for row,saved in searches:
  got=json.loads(json.dumps(E.search(row,labels['node_limit'])))
  assert got==saved and got['verdict']=='EXACT_COMPLETE_OBSTRUCTION'
  counts['replayed_search_states']+=got['nodes']
 for kind in ['gale','positive_kernel','facet']:
  forged=copy.deepcopy(d['cases'][0])
  if kind=='gale':forged['gale'][0][0]+=1
  elif kind=='positive_kernel':forged['positive_kernel'][0]+=1
  else:forged['facets'][0]^=1
  try:geometry(forged)
  except AssertionError:counts['forged_'+kind+'_rejected']+=1
  else:raise AssertionError('forged geometry accepted: '+kind)
 paths=[GEOM,LABELS,PAT,REL,LOCAL,BLOCK,f'scripts/{TAG}_independent.py',
  'scripts/fresh_h9_newfacet_frontier_20261014_engine.py','scripts/fresh_h9_bugaenko_closure_20261014.py',
  'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py','scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
 return dict(status='PASS_SEPARATE_GALE_PRIMAL_AND_GATE_REPLAY',sources={p:sha(p)for p in paths},
  nerve_count=len(rows),counts=dict(counts),cases=rows,
  scope='FLINT integer determinants and rational linear algebra independently reconstruct all saved Gale circuits, positive kernels and explicit affine primal coordinates, every supporting facet and its affine rank8, all minimal nonfaces and literal local-gate embeddings. No primary Gale/clutter/canonicalization code is imported. Canonical relabeling maps are checked, but nauty uniqueness is inherited. All23 complete searches are replayed with separate finite-tree/modular arithmetic; DFS and relation propagation are reused. Previous44 local obstruction and two-block proof certificates are inherited source-bound theorems, not independently reproved here. The180000 proposal moves are not re-enumerated by this verifier, and no whole14-facet classification is claimed.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}_independent.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['nerve_count'],out['counts'])
