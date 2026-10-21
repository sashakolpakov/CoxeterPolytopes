#!/usr/bin/env python3
"""Complete small affine-component maps for the new distance-eight overlap domain.

Retains the frozen ideal5 pair sets and checks every newly selected pattern.
No optimization or dimension assertion.
"""
import argparse,json,hashlib
from collections import Counter
from pathlib import Path
import fresh_ideal5_small_factor_payments_20261017 as S
P=S.P;ROOT=P.ROOT
TAG='fresh_distance8_ideal5_small_factor_20261017'
OVERLAP='fresh_distance8_high_overlap_20261017'

def key(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['match'])))
def run():
 data=S.read(OVERLAP);old=S.read(S.TAG);prior=S.read(P.OLD)
 diagrams={r['index']:r for r in S.read(P.CENSUS)['diagrams']}
 patterns=[]
 assert len(data['isolation']['cases'])==2614
 for row in data['isolation']['cases']:
  if row.get('component_order_bound',0)<14:continue
  di,_=row['left'];dj,_=row['right'];gram,_=P.R.unite(P.R.bond_graph(diagrams[di]),P.R.bond_graph(diagrams[dj]),tuple(map(tuple,row['match'])))
  ids=row['marked_component'];labels=[[int(2-gram[i][j])if i!=j and gram[i][j]else 0 for j in ids]for i in ids]
  patterns.append(dict(index=row['index'],left=row['left'],right=row['right'],match=row['match'],nodes=ids,labels=labels,
    marks=[ids.index(x)for x in row['preferred_corner']],order_bound=row['component_order_bound'],method=row['method']))
 assert len(patterns)==1153
 assert Counter(p['order_bound']for p in patterns)=={14:313,15:354,16:346,17:140}
 oldkey={key(r):i for i,r in enumerate(prior['isolation']['cases'])};cache={};records=[]
 for rec in old['records']:
  m=rec['order'];kind=rec['kind'];target=rec['target'];paid=set(map(tuple,rec['paid_pairs']));drop=set(map(tuple,rec['dropped_pairs']))
  union=set();rows=[];affected=[]
  for p in patterns:
   if p['order_bound']<m or len(p['nodes'])>m:continue
   cachekey=(m,kind,str(p['labels']))
   if cachekey not in cache:cache[cachekey]=P.embeddings(p['labels'],target)
   maps=cache[cachekey];pairs={tuple(sorted(mp[i]for i in p['marks']))for mp in maps}
   union|=pairs
   if maps:rows.append(dict(case=p['index'],embeddings=maps,marked_images=[list(x)for x in sorted(pairs)]))
   if pairs&drop:affected.append(p)
  assert union<=paid|drop
  correspondence=[]
  if kind=='D'and m in(16,17):
   for p in affected:
    assert key(p)in oldkey
    correspondence.append(dict(new_case=p['index'],old_case=oldkey[key(p)],left=p['left'],right=p['right'],match=p['match']))
   assert sorted(r['old_case']for r in correspondence)==rec['dropped_affected_cases']
   assert len(correspondence)==9
  records.append(dict(order=m,kind=kind,target=target,preferred_pairs=[list(x)for x in sorted(union)],
    retained_paid_pairs=rec['paid_pairs'],retained_cost=rec['cost'],dropped_pairs=rec['dropped_pairs'],case_embeddings=rows,
    dropped_affected_cases=[p['index']for p in affected],old_sealed_case_correspondence=correspondence))
  print('NEW_IDEAL_MAPS',m,kind,len(union),len(paid),len(affected),flush=True)
 stems=[OVERLAP,S.TAG,P.OLD,P.CENSUS,'fresh_affine_d13_long_pair_transfer_20261017',
  'fresh_affine_d14_long_pair_transfer_20261017','fresh_order16_finite_cross_payments_20261017',P.TAG,
  'fresh_cycle14_ideal5_transfer_20261017','fresh_finite8_ideal5_local_replacement_20261017',
  'fresh_finite8_cycle_replacement_20261017']
 paths=[Path(__file__),ROOT/'scripts'/f'{S.TAG}.py',ROOT/'scripts'/f'{P.TAG}.py']
 paths += [ROOT/f'{folder}/{s}.{ext}'for s in stems for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
 return dict(status='PASS_COMPLETE_DISTANCE8_SMALL_FACTOR_PAIR_INCLUSION',patterns=patterns,
  order_bound_counts=dict(Counter(p['order_bound']for p in patterns)),records=records,
  sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
  scope='All1153 selected patterns with possible affine order14..17 use the frozen ideal5 pair sets or the validated universal/literal replacements. No new ideal cost and no dimension claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],flush=True)
