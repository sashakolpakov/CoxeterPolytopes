#!/usr/bin/env python3
"""Raw23 aggregate-vertex relative roles, before an incidence/charge claim.

All ten finite supports and all235 source classes are inventoried. A rank13
host containing both active forests realizes every component-overlap pattern.
Identical recovered native stars are one aggregate, not coexistence events.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement,combinations
from pathlib import Path
import fresh_raw357_relative_marker_pilot_20261018 as V
G=V.G;C=V.C;A=V.A;R=V.R;L=V.L;T=V.T;H=V.H;ROOT=V.ROOT
TAG='fresh_raw23_same_family_relative_pilot_20261018'
WITNESS='fresh_raw23_surviving_coface_witness_20261018'
INC='fresh_sqrt5_prism_simplex_base_incidence_20261017'
RAW='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017'
QUAD='fresh_sqrt5_rank13_packing_quad_separation_20261017'
def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())

def native_maps(g,ids,om,target,kind):
 J,_,_,tens=T.projection_data(target)
 maps=R.active_embeddings(C.active_graph(g,ids),target);assert maps
 heights=[]
 anchor=maps[0]
 cov=[tuple(sorted((anchor[j],tuple(g[i][j]))for j in ids if g[i][j]!=A.RIGHT))for i in om]
 for cv in cov:
  field=L.field_for((cv,),tens);z=L.projection(cv,cv,J,tens,field)
  assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;heights.append(h)
 field=L.field_for(cov,tens);z=L.projection(cov[0],cov[1],J,tens,field)
 product=field.negative_contact(g[om[0]][om[1]])
 assert field.signs_difference(tuple(x-y for x,y in zip(z,product)),heights[0]*heights[1])==0
 out=[];seen=set()
 for mp in maps:
  assert set(mp)==set(ids)and len(set(mp.values()))==5
  assert all(set(target[mp[i]])<=set(mp.values())for i in ids)
  assert all(target[mp[i]].get(mp[j],2)==g[i][j][0]for i,j in combinations(ids,2))
  cv=[tuple(sorted((mp[j],tuple(g[i][j]))for j in ids if g[i][j]!=A.RIGHT))for i in om]
  order=sorted(range(2),key=lambda i:cv[i])
  native=dict(contacts=[[[j,list(e)]for j,e in cv[i]]for i in order],
   squared_timelike_heights=[heights[i].out()for i in order],
   native_omitted_product=list(g[om[0]][om[1]]),occupied_pairs=[])
  key=H.frozen(native['contacts'])
  assert key not in seen;seen.add(key)
  out.append(dict(index=len(out),kind=kind,active_embedding=sorted(mp.items()),native=native))
 return out

def run(pairs=False):
 started=time.monotonic();raw=read(RAW);quad=read(QUAD);inc=read(INC);w=read(WITNESS)
 for p,h in w['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 c=quad['canonical_cores'][19];assert c['source_core']==23 and c['required']==10
 g=A.V.raw_graph(raw['records'][23]);assert len(g)==7
 scaled,sc=C.scaled(g);assert C.P.congruence(scaled)['inertia']==[5,1,1]
 supports=[]
 for om in combinations(range(7),2):
  ids=sorted(set(range(7))-set(om));inertia=C.P.congruence(C.P.sub(scaled,ids))['inertia']
  if inertia==[5,0,0]:supports.append(dict(kind=len(supports),omitted=list(om),support=ids,forest=C.forest_type(g,ids)))
 assert [r['omitted']for r in supports]==c['pd_deletions']and len(supports)==10
 occurrence=Counter(tuple(sorted(r['vertex_omitted']))for r in c['aggregate_corner_occurrences'])
 assert set(occurrence)=={tuple(r['omitted'])for r in supports}and set(occurrence.values())=={4}
 autos=C.V.automorphisms(g)
 # Unique complete source contexts, retaining all zero-mass types too.
 actual=read('fresh_sqrt5_actual_provenance_dual_20261017')
 prior=read('fresh_all_marker_families_dual_continuation_20261018')
 explicit={actual['provenance_masses'][pi]['source_context']:f for pi,f in prior['prepared']['provenance_to_face']}
 whole=defaultdict(list)
 for ctx,f,ci,k,bs in inc['contexts']:whole[explicit.get(ctx,f)].append(ctx)
 alpha={f:F(a)for f,a in w['candidate']['face_multipliers']};family=[]
 for ctx,f,ci,k,bs in inc['contexts']:
  if ci!=19:continue
  f=explicit.get(ctx,f);assert whole[f]==[ctx]and k==5
  family.append(dict(source_context=ctx,face=f,bases=bs,required_quads=10,current_mass=str(alpha.get(f,F()))))
 assert len(family)==235 and len({r['face']for r in family})==235
 assert sum(F(r['current_mass'])>0 for r in family)==50
 models=dict(R.C.models());cases={};counts=Counter();hosts=[]
 for a,b in combinations_with_replacement(supports,2):
  forest=sorted(a['forest']+b['forest']+['A1']*3)
  assert sum(len(models[s])for s in forest)==13
  target=R.P.forest_graph(forest,models)
  left=native_maps(g,a['support'],a['omitted'],target,a['kind'])
  right=left if a['kind']==b['kind']else native_maps(g,b['support'],b['omitted'],target,b['kind'])
  anchor=left[0];local=Counter()
  for other in right:
   if H.frozen(anchor['native']['contacts'])==H.frozen(other['native']['contacts']):
    counts['identical_aggregate_presentations']+=1;continue
   key,pattern,payload=V.signature(target,a['kind'],anchor,b['kind'],other)
   counts['anchor_pairs']+=1;local['anchor_pairs']+=1
   if key in cases:assert cases[key]['payload']==V.normal(payload)
   else:cases[key]=dict(kinds=[a['kind'],b['kind']],omissions=[a['omitted'],b['omitted']],pattern=pattern,
    payload=V.normal(payload),forest=forest,finite_frame=target,anchor=anchor,other=other)
  hosts.append(dict(kinds=[a['kind'],b['kind']],forest=forest,left_orientations=len(left),right_orientations=len(right),**local))
  print('HOST',a['kind'],b['kind'],'cases',len(cases),'seconds',round(time.monotonic()-started,2),flush=True)
 assert len(hosts)==55
 records=[]
 for key,z in sorted(cases.items()):
  if pairs:
   assert len(cases)<=1200,'Bounded pilot threshold exceeded; inspect preflight first.'
   J,_,_,tens=T.projection_data(z['finite_frame'])
   z['result']=L.pair_test(J,tens,z['anchor']['native'],z['other']['native'])
  records.append(z)
 paths={f'scripts/{TAG}.py'}|{f'results/{s}.json'for s in(RAW,QUAD,INC,WITNESS,
  'fresh_sqrt5_actual_provenance_dual_20261017','fresh_all_marker_families_dual_continuation_20261018')}
 paths|={str(Path(m.__file__).relative_to(ROOT))for m in(V,G,C,A,R,L,T,H)}
 return V.normal(dict(status='PASS_RAW23_SAME_FAMILY_RELATIVE_NATIVE_PILOT'if pairs else'PASS_RAW23_SAME_FAMILY_RELATIVE_PREFLIGHT',
  sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
  geometry=dict(native_gram=g,scaled_gram=C.P.out(scaled),squared_scales=[x.out()for x in sc],
   native_automorphisms=autos,finite_vertex_types=supports,quad_partitions=c['quad_partitions'],
   quadrilateral_corner_occurrences=len(c['aggregate_corner_occurrences']),actual_vertex_occurrences=10),
  complete_source_family=family,hosts=hosts,relative_cases=records,
  counts=dict(counts,source_classes=235,positive_source_classes=50,finite_supports=10,
   host_pairs=55,relative_cases=len(records),native_pair_tests=len(records)if pairs else 0,
   retained_relative_cases=sum(r['result']['allowed']for r in records)if pairs else None),
  old_graph_outcomes_changed=False,exported_columns=[],optimizer_invocations=0,
  scope='Complete relative-role domain for two raw23 aggregate occurrences at one finite vertex, with distinct ordered whole native roles and conservative native gates if requested. One aggregate contributes one event at each of10 vertices and4 quadrilateral corners there, not an orientation multiplier. No extra-vertex incidence, graph capacity, charge, matrix change or optimized bound is claimed.'))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--pairs',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 path=ROOT/f'results/{TAG}{""if args.pairs else"_preflight"}.json'
 assert args.verify or not path.exists();out=run(args.pairs)
 if args.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],flush=True)
