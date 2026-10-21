#!/usr/bin/env python3
"""Independent literal17-state closure and geometry/filter review.
Reuses only the new independent graph/tuple-Fraction engines read-only.
The all-rank single-chain seed-reduction theorem is an explicit dependency.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_extended_forest_weights_20261015_review as R
import fresh_sqrt5_sealed_subcore_closure_20261015_review as S
import fresh_sqrt5_forest_sensitive_weights_20261015_review as W
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_singular_attachment_closure_20261015_review'

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=read('results/fresh_sqrt5_singular_attachment_closure_20261015.json');first=read('results/fresh_sqrt5_singular_extension_frontier_20261015.json');chain=read('results/fresh_sqrt5_single_chain_completion_20261015.json');old=read('results/fresh_sqrt5_quad_contacts_20261015.json');base=read('results/fresh_sqrt5_quad_exhaustion_20261015.json')['labelled_union_supports'];extended=read('results/fresh_sqrt5_extended_forest_weights_20261015.json');original=read('results/fresh_sqrt5_forest_sensitive_weights_20261015.json')
 assert d['complete']and d['processed']==d['marked_states']==17 and len(d['records'])==17
 weights=R.old_weights(original)
 for v in extended['upgraded_orbits']:weights[W.decode(v['key'])]=tuple(map(F,v['new']))
 roots=[('old',r['class_index'],base[r['union_index']])for r in old['records']if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and r['singular']and not r['accepted']]+[('new',r['index'],r['row'])for r in first['records']if not r['accepted']]
 assert len(roots)==158
 pool=defaultdict(list);targets={};sizes=set()
 for fam,i,row in roots:
  M=S.matrix(row);iv=S.invariant(M);assert all(S.isomorphism(M,N)is None for _,_,N in pool[iv]);pool[iv].append((fam,i,M));targets[fam,i]=M;sizes.add(len(row['white']))
 seeds=[]
 for rec in first['records']:
  for c in rec['accepted']:seeds.append((dict(kind='seven_continuations',source_index=rec['index'],contacts=c),R.extend(rec['row'],c)))
 for k,r in enumerate(chain['records']):
  if r.get('length',0)>1:
   assert r['status']=='UNIQUE_LENGTH_UNPAID_PRINCIPAL_STAR';seeds.append((dict(kind='longer_A_chain',source_record=k,length=r['length']),r['row']))
 assert len(seeds)==23 and Counter(r[0]['kind']for r in seeds)=={'seven_continuations':7,'longer_A_chain':16}
 # Unique full literal origin matching; source order is not trusted.
 statepool=defaultdict(list);origins=[];originset=set();seedmaps=[]
 for rec in d['records']:
  M=S.matrix(rec['row']);iv=S.invariant(M);assert all(S.isomorphism(M,N)is None for N in statepool[iv]);statepool[iv].append(M)
  for o in rec['origins']:
   ident={k:v for k,v in o.items()if k!='to_representative'}
   j=next(j for j,(meta,_)in enumerate(seeds)if meta==ident);assert j not in originset;originset.add(j)
   SM=S.matrix(seeds[j][1]);assert S.isomorphism(SM,M)is not None;S.check_map(SM,M,o['to_representative']);seedmaps.append((j,rec['index'],o['to_representative']))
 assert originset==set(range(23))
 for j,saved in enumerate(d['seed_maps']):
  actual=next(r for r in seedmaps if r[0]==j);assert saved['target']==actual[1]and saved['map']==actual[2]and saved['kind']==seeds[j][0]['kind']
 stops=[]
 for rec in first['records']:
  if rec['accepted']:
   corners,sums=R.payments(rec['row'],weights);assert min(sums)>=1;saved=next(r for r in d['initial_five_zero_further_contact_cases']if r['index']==rec['index']);assert sums==list(map(F,saved['sums']));stops.append(dict(index=rec['index'],sums=list(map(str,sums))))
 assert len(stops)==5
 count=Counter();trials=Counter();subsets=0;records=[]
 for rec in d['records']:
  row=rec['row'];M=S.matrix(row);matches=[]
  for size in sorted(sizes):
   for white in combinations(row['white'],size):
    ids=list(range(4))+list(white);SM=[[M[i][j]for j in ids]for i in ids];subsets+=1
    for fam,i,T in pool.get(S.invariant(SM),[]):
     mp=S.isomorphism(SM,T)
     if mp is None:continue
     S.check_map(SM,T,mp);outside=set(row['white'])-set(white)
     assert all(M[a][b]==2 for a in outside for b in ids);matches.append((fam,i,tuple(ids)))
  expected=[]
  for m in rec['known_subcores']:
   ids=m['nodes'];SM=[[M[i][j]for j in ids]for i in ids];S.check_map(SM,targets[m['family'],m['index']],m['to_target']);assert m['outside_white']==sorted(set(row['white'])-set(ids));expected.append((m['family'],m['index'],tuple(ids)))
  assert sorted(matches)==sorted(expected)
  if matches:
   assert rec['status']=='REDUCES_TO_KNOWN_158_SEALED_CATALOGUE';out=dict(index=rec['index'],status=rec['status'],subcores=[dict(family=fam,index=i,nodes=list(ids))for fam,i,ids in sorted(matches)])
  else:
   mapping={'PAID_H4':'BECOMES_PAID_ENTIRE_H4','REJECT_FIELD':'REJECT_FINITE_CYCLE_FIELD','REJECT_RELATION':'REJECT_PHYSICAL_RELATION_AND_CONJUGATE_KERNEL','REJECT_AFFINE':'REJECT_AFFINE','REJECT_TWO_NEGATIVE':'REJECT_TWO_NEGATIVE','REJECT_CONJUGATE':'REJECT_CONJUGATE_NOT_PSD','ACCEPT_UNPAID':'PASS_UNPAID_EXTRA_WHITE'}
   transformed=dict(index=rec['index'],row=row,profiles=len(rec['trials']),accepted=rec['accepted'],trials=[dict(t,status=mapping[t['status']])for t in rec['trials']])
   check,unused=R.audit_domain(transformed);trials.update(t['status']for t in rec['trials']);assert not check['accepted']and not rec['accepted']
   corners,sums=R.payments(row,weights);assert min(sums)>=1 and sums==list(map(F,rec['payment']['sums_without_with_external_a1']))
   assert rec['status']=='NEW_SEALED_PAID';out=dict(index=rec['index'],status=rec['status'],domain=check,payment_sums=list(map(str,sums)))
  count[rec['status']]+=1;records.append(out)
 assert subsets==3054 and dict(count)==d['counts']and dict(trials)==d['trial_counts']
 assert count=={'REDUCES_TO_KNOWN_158_SEALED_CATALOGUE':6,'NEW_SEALED_PAID':11}and sum(trials.values())==138
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_singular_attachment_closure_20261015.py','results/fresh_sqrt5_singular_attachment_closure_20261015.json','scripts/fresh_sqrt5_extended_forest_weights_20261015_review.py','results/fresh_sqrt5_extended_forest_weights_20261015_review.json','scripts/fresh_sqrt5_sealed_subcore_closure_20261015_review.py','scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py','results/fresh_sqrt5_single_chain_completion_20261015.json','scripts/fresh_sqrt5_single_chain_completion_20261015_review.py','results/fresh_sqrt5_single_chain_completion_20261015_review.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_INDEPENDENT_COMPLETE_17_STATE_ATTACHMENT_QUEUE',sources={p:sha(p)for p in paths},literal_seeds=23,distinct_states=17,known_sealed_catalogue=158,known_subcore_subset_checks=subsets,counts=dict(sorted(count.items())),trial_counts=dict(sorted(trials.items())),initial_five_zero_contact_cases=stops,records=records,queue_complete=True,cutoff_reached=False,method='Independent colored-graph isomorphism and all eligible induced-white subsets; separate tuple-Fraction physical/conjugate inertia and exactPD-pruned contact recursion; every literal seed and origin rebuilt.',scope='Complete closure of the23 explicit seeds, conditional on the separately stated all-rank single-chain seed-reduction theorem. No arbitrary nonsingular/dotted-frontier closure or new dimension bound.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],out['trial_counts'])
