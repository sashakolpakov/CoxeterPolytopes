#!/usr/bin/env python3
"""Transfer proved predicates through a literal automorphism of star0."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_single_right_20261014_symmetry'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
PRIMARY='results/fresh_h9_star_single_right_20261014.json'
REVIEW='results/fresh_h9_star_single_right_20261014_review.json'
PERM=[1,0,17,12,13,16,9,10,15,6,7,14,3,4,11,8,5,2]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen={}
 def read(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for q,h in d.get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):read(q)
  return d
 g=read(GEOM)['representatives'][0]['geometry'];primary=read(PRIMARY);review=read(REVIEW)
 assert primary['case']==0 and primary['status']=='PASS_TWO_NONCORE_CUT_LEAVES_UNRESTRICTED'
 assert review['status'].startswith('PASS_')and review['sources'][PRIMARY]==sha(PRIMARY)
 p=PERM;assert sorted(p)==list(range(18))and all(p[p[i]]==i for i in range(18))
 V=[tuple(sorted(v['incident']))for v in g['vertices']]
 M=[tuple(sorted(m))for m in g['minimal_nonfaces']]
 vm={v:i for i,v in enumerate(V)};mm={m:i for i,m in enumerate(M)}
 assert len(vm)==len(V)==378 and len(mm)==len(M)
 vp=[vm[tuple(sorted(p[i]for i in v))]for v in V]
 mp=[mm[tuple(sorted(p[i]for i in m))]for m in M]
 assert sorted(vp)==list(range(len(V)))and sorted(mp)==list(range(len(M)))
 def key(f):return tuple(tuple(sorted(tuple(row)for row in f[k]))for k in ('fixed','bounds'))
 families=primary['families'][:];lookup={key(f):f['id']for f in families};transfers=[]
 assert len(lookup)==len(families)==161
 for f in primary['families']:
  image={k:sorted([sorted((p[a],p[b]))+[m]for a,b,m in f[k]])for k in ('fixed','bounds')}
  assert all(0<=a<b<18 and 2<=m<=7 for k in image for a,b,m in image[k])
  target=lookup.get(key(image));added=target is None
  if added:
   target='reversal_'+f['id'];lookup[key(image)]=target
   families.append(dict(id=target,origin='literal_star_zero_reversal',**image))
  transfers.append(dict(source=f['id'],target=target,added=added,**image))
 assert len(families)==173 and sum(r['added']for r in transfers)==12
 sources=[GEOM,PRIMARY,REVIEW,f'scripts/{TAG}.py']
 return dict(status='PASS_LITERAL_STAR_REVERSAL_GATE_TRANSFER',case=0,
  sources={p:sha(p)for p in sources},permutation=p,vertex_permutation=vp,
  minimal_nonface_permutation=mp,transfers=transfers,families=families,
  original_predicates=161,new_predicates=12,
  proof='The involution permutes all378 literal actual vertices and the complete minimal-nonface clutter. Relabeling any compact Coxeter Gram realization by it preserves every incidence, label, dotted parameter, rank and signature. Thus every mapped predicate is excluded by its original hash-bound theorem. Each equality and lower-bound pair is explicitly mapped; no finite-label cutoff is introduced.',
  scope='Adds twelve literal automorphism images to161 reviewed predicates on star0. Does not itself classify the remaining label assignments or claim any automorphism of other catalogue cases.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
