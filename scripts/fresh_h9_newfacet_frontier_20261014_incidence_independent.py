#!/usr/bin/env python3
"""Independent literal recipe replay for all5536 changed-order nerves."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import clutter
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_newfacet_frontier_20261014_incidence_independent'
INC='results/fresh_h9_bugaenko_remaining_20261014_newfacet.json';RAW='data/known_examples/bugaenko_compact_h8.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def cut(vertices,face,new):
 f=frozenset(face);assert 2<=len(f)<=5 and any(f<=v for v in vertices)
 return {v for v in vertices if not f<=v}|{frozenset((v-{a})|{new})for v in vertices if f<=v for a in f}
def run():
 d=json.loads((ROOT/INC).read_text());raw=json.loads((ROOT/RAW).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 base={frozenset(v)for v in raw['geometry']['vertices']};records=[];canon=set()
 for row in d['cases']:
  first=cut(base,row['first_face'],11);wedged=set()
  for v in first:
   if 11 in v:wedged.add(v|{12})
   else:wedged.update([v|{11},v|{12}])
  vertices=cut(wedged,row['second_face'],13)
  assert all(len(v)==9 for v in vertices)
  masks=[sum(1<<i for i in v)for v in vertices]
  mm=clutter(masks,14)
  assert max(map(int.bit_count,mm))<=5
  order=row['canonical_order'];assert sorted(order)==list(range(14));mapping={a:i for i,a in enumerate(order)}
  def move(m):return sum(1<<mapping[i]for i in range(14)if m>>i&1)
  assert sorted(map(move,masks))==row['canonical_vertices']
  assert sorted(map(move,mm))==row['canonical_minimal_nonfaces']
  assert tuple(row['canonical_minimal_nonfaces'])not in canon;canon.add(tuple(row['canonical_minimal_nonfaces']))
  degrees=Counter(frozenset(v-{i})for v in vertices for i in v);assert set(degrees.values())=={2}
  records.append(dict(case=row['case'],vertices=len(vertices),minimal_nonfaces=len(mm),ridges=len(degrees)))
 assert len(records)==5536
 files=[INC,RAW,f'scripts/{TAG}.py','scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py']
 return dict(status='PASS_SEPARATE_ALL5536_CHANGED_ORDER_INCIDENCES',sources={p:sha(p)for p in files},cases=len(records),actual_vertices=sum(r['vertices']for r in records),records=records,scope='Every saved changed recipe is rebuilt by literal primal vertex surgeries. Minimal nonfaces are recovered independently from the full downward complex; canonical maps and ridge degrees are checked. No producer surgery formula or nauty canonicalizer is imported. This replays every saved type, but the safe operation-orbit completeness and isomorphism deduplication remain producer/source-review premises. Polytopality follows from wedge and sufficiently small face truncations.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['actual_vertices'],flush=True)
