#!/usr/bin/env python3
"""Independent full-face incidence replay; imports no surgery producer."""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_surgery_20261014'
def bits(m):
 while m:
  b=m&-m;yield b;m-=b
@lru_cache(None)
def downward(v):
 if not v:return 1
 b=v&-v;a=downward(v^b)
 return a|(a<<b)
@lru_cache(None)
def coordinate_masks(n):
 return tuple(sum(1<<j for j in range(1<<n)if j&(1<<i))for i in range(n))
def clutter(vertices,n):
 faces=0
 for v in vertices:faces|=downward(v)
 allbits=(1<<(1<<n))-1;missing=allbits^faces
 for i,has in enumerate(coordinate_masks(n)):
  b=1<<i
  missing&=(faces<<b)|(allbits^has)
 return tuple(b.bit_length()-1 for b in bits(missing))
def cut(vertices,face,new):
 assert 2<=len(face)<=5
 f=frozenset(face);assert any(f<=v for v in vertices)
 return {v for v in vertices if not f<=v}|{(v-{a})|{new}for v in vertices if f<=v for a in f}
def masks(vertices):return sorted(sum(1<<i for i in v)for v in vertices)
def run():
 source=f'data/known_examples/bugaenko_compact_h8.json';path=f'results/{TAG}_incidence.json'
 d=json.loads((ROOT/path).read_text());base=json.loads((ROOT/source).read_text())
 for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
 old={frozenset(v)for v in base['geometry']['vertices']};assert len(old)==41 and all(len(v)==8 for v in old)
 assert list(clutter(masks(old),11))==d['base_minimal_nonfaces']
 rows=[];total_vertices=0;canonical=set()
 for r in d['cases']:
  c=r['wedge_carrier'];wedged=set()
  for v in old:
   if c in v:wedged.add(v|{11})
   else:wedged.update([v|{c},v|{11}])
  v=cut(cut(wedged,r['first_face'],12),r['second_face'],13)
  mm=masks(v);assert mm==r['actual_vertices']
  C=clutter(mm,14);assert list(C)==r['minimal_nonfaces'] and max(map(int.bit_count,C))<=5
  assert all(len(w)==9 for w in v)
  ridges=Counter(w-{i}for w in v for i in w);assert set(ridges.values())=={2}
  order=r['canonical_order'];assert sorted(order)==list(range(14));mp={a:i for i,a in enumerate(order)}
  cv=masks({frozenset(mp[i]for i in w)for w in v})
  cm=sorted(sum(1<<mp[b.bit_length()-1]for b in bits(m))for m in C)
  assert cv==r['canonical_vertices'] and cm==r['canonical_minimal_nonfaces']
  assert tuple(cm)not in canonical;canonical.add(tuple(cm));total_vertices+=len(v)
  rows.append(dict(case=r['case'],vertices=len(v),minimal_nonfaces=len(C),eight_faces=len(ridges)))
 assert len(rows)==2953 and sum(r['recipe_orbit_multiplicity']for r in d['cases'])==5986
 files=[source,path,f'scripts/{TAG}_incidence.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return dict(status='PASS_SEPARATE_FULL_FACE_INCIDENCE_REPLAY',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},cases=len(rows),actual_vertices_checked=total_vertices,records=rows,
  scope='Every saved recipe is rebuilt by primal vertex wedge and dual stellar truncation. Minimal nonfaces are independently recovered from the entire downward face complex, not the producer clutter formula. All canonical relabelings and ridge degrees are checked. Orbit-coverage and final isomorphism deduplication are source-reviewed and producer-replayed, not separately re-enumerated here. Polytopality follows from wedge and sufficiently small face truncations of the certified compact source.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');o=a.parse_args();d=run();p=ROOT/f'results/{TAG}_incidence_independent.json'
 if o.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['cases'],d['actual_vertices_checked'])
