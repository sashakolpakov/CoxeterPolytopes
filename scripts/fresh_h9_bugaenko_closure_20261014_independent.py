#!/usr/bin/env python3
"""Separate case-zero proof: explicit automorphisms, finite gates, signatures."""
import argparse,json
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
import sympy as S
from fresh_h9_bugaenko_surgery_20261014_rank_independent import spherical
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_closure_20261014'
OLD='fresh_h9_bugaenko_surgery_20261014';PAIRS=list(combinations(range(14),2));INDEX={p:i for i,p in enumerate(PAIRS)}
def tup(m):return tuple(i for i in range(14)if m>>i&1)
def run():
 paths=[f'results/{TAG}.json',f'results/{OLD}_incidence.json',f'results/{OLD}_labels_0_25.json']
 d,inc,joins=[json.loads((ROOT/p).read_text())for p in paths]
 for source in[d,inc,joins]:
  for p,h in source['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 assert d['status']=='PASS_COMPLETE_CASE0_CLOSURE'and d['metric_candidate_count']==0
 case=inc['cases'][0];M=list(map(tup,case['canonical_minimal_nonfaces']));V=list(map(tup,case['canonical_vertices']))
 row=joins['cases'][0];assert not row['state_limit_reached']and row['unprocessed_missing_faces']==0
 states=row['terminal_states'];assert len(states)==864
 blocks=d['twin_blocks'];assert sorted(i for b in blocks for i in b)==list(range(14))
 for block in blocks:
  assert len({tuple(i in m for m in M)for i in block})==1
 group=[]
 for choices in product(*(tuple(permutations(b))for b in blocks)):
  p=list(range(14))
  for b,c in zip(blocks,choices):
   for i,j in zip(b,c):p[i]=j
  assert {tuple(sorted(p[i]for i in m))for m in M}==set(M)
  assert {tuple(sorted(p[i]for i in v))for v in V}==set(V)
  group.append(p)
 assert group==d['verified_literal_group']and len(group)==288
 domain=[tuple(p)for p,m in states[0]];loc={p:i for i,p in enumerate(domain)}
 assert all([tuple(p)for p,m in st]==domain for st in states)
 vectors=[tuple(m for p,m in st)for st in states];lookup={v:i for i,v in enumerate(vectors)};assert len(lookup)==864
 covered=set();records=[];minors={}
 for orbit in d['orbits']:
  ii=orbit['representative'];values=vectors[ii];members=set()
  for p in group:
   moved=tuple(values[loc[tuple(sorted((p[a],p[b])))]]for a,b in domain);assert moved in lookup;members.add(lookup[moved])
  assert sorted(members)==orbit['members']and len(members)==144 and not covered&members;covered|=members
  labels=[2]*91;assigned={}
  for p,m in states[ii]:k=INDEX[tuple(p)];labels[k]=m;assigned[k]=m
  assert all(spherical(tuple(labels[INDEX[p]]for p in combinations(v,2)))for v in V)
  assert len(orbit['tree'])==1 and orbit['tree'][0]['kind']=='ORTHOGONAL_LANNER'
  A,B=[M[i]for i in orbit['tree'][0]['missing_faces']];assert not set(A)&set(B)
  negative=[]
  for block in[A,B]:
   if len(block)==2:
    assert block in M;negative.append(dict(support=block,type='dotted_pair',determinant='1-x^2<0 for every x>1'));continue
   code=tuple(labels[INDEX[p]]for p in combinations(block,2));assert all(INDEX[p]in assigned for p in combinations(block,2))
   if code not in minors:
    G=S.eye(len(block))
    for(i,j),m in zip(combinations(range(len(block)),2),code):G[i,j]=G[j,i]=-S.cos(S.pi/m)
    ds=[S.simplify(G[:k,:k].det())for k in range(1,len(block)+1)]
    assert all(x.is_positive is True for x in ds[:-1])and ds[-1].is_negative is True
    minors[code]=list(map(str,ds))
   negative.append(dict(support=block,type='exact_Lorentzian_Gram',labels=code,leading_minors=minors[code]))
  cross=[]
  for a in A:
   for b in B:
    p=tuple(sorted((a,b)));k=INDEX[p]
    if k in assigned:assert labels[k]==2;cross.append(dict(pair=p,reason='assigned_right'));continue
    labels[k]=3;bad=next((j for j,v in enumerate(V)if not spherical(tuple(labels[INDEX[p]]for p in combinations(v,2)))),None);labels[k]=2
    assert bad is not None
    cross.append(dict(pair=p,reason='label3_relaxation_nonspherical',witness_vertex=bad))
  records.append(dict(representative=ii,members=sorted(members),negative_blocks=negative,forced_cross_pairs=cross))
 assert covered==set(range(864))
 files=paths+[f'scripts/{TAG}.py',str(Path(__file__).resolve().relative_to(ROOT)),f'scripts/{OLD}_rank_independent.py']
 return dict(status='PASS_INDEPENDENT_COMPLETE_CASE0_NEGATIVE_INDEX_OBSTRUCTION',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},partial_labelings_covered=864,literal_group_order=288,orbits=records,signature_negative_index_at_least=2,scope='Independent orbit reconstruction, exact SymPy leading-minor signatures, and separate finite-diagram classifier for every forced orthogonality. The complete864 initial Lanner join remains an explicitly source-bound upstream completeness dependency. No metric optimization or arbitrary distance bound is used.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}_independent.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
