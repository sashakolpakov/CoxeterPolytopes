#!/usr/bin/env python3
"""Independent graph/Fraction replay of unique-A1/E7 cross weights.
Only the frozen separate224-review graph engine is imported, no producers.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_enlarged_a1_weights_20261015_review as R
ROOT=Path(__file__).resolve().parents[1];PRIMARY='fresh_sqrt5_parameter_weight_stress_20261015_cross';TAG=PRIMARY+'_review'
def load(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def inputs():
 rows=[(fam,i,row,None)for fam,i,row in R.literal_inputs()]
 families=load('results/fresh_sqrt5_subdivision_families_20261015.json')
 assert len(families['rank_at_most29_patterns'])==191
 for r in families['rank_at_most29_patterns']:rows.append(('subdivision',str(r['source_index'])+':'+str(r['subdivisions']),r['row'],None))
 one=load('results/fresh_sqrt5_rooted_component_completion_20261015.json');assert len(one['rank29_instances'])==1715
 for r in one['rank29_instances']:rows.append(('one_touch',r['index'],r['row'],None))
 dotted=load('results/fresh_sqrt5_dotted_rank29_exhaustion_20261015.json')
 assert dotted['processed']==dotted['marked_states']==438 and dotted['status']=='COMPLETE_DOTTED_ARITHMETIC_ATTACHMENT_QUEUE_THROUGH_DIMENSION29'
 selected=[r for r in dotted['records']if r['result']['singular']]
 assert[r['index']for r in selected]==dotted['weight_obligations']and len(selected)==264
 rows.extend(('dotted',r['index'],r['row'],r['ratio'])for r in selected)
 assert len(rows)==2397;return rows

def corners(row):
 M=R.label_matrix(row);es=[(i,j,M[i][j])for i,j in combinations(range(len(M)),2)];out=[]
 for a,b in product((0,2),(1,3)):
  parts=R.components(R.graph(row['white']+[a,b],es));left=next(g for g in parts if a in g);right=next(g for g in parts if b in g)
  count=sum(len(g)==1 for g in parts);base=F(0);key=None;cross=None
  if b in left:base=R.baseline(left,a,b);key=R.tree_key(left,(a,b))
  elif count==1:
   if len(left)==1 and R.type_of(right)=='E7':cross=R.tree_key(right,(b,))
   if len(right)==1 and R.type_of(left)=='E7':cross=R.tree_key(left,(a,))
  out.append(dict(marked=[a,b],types=sorted(R.type_of(g)for g in parts),A1=count,baseline=base,key=key,cross=cross))
 return out

def run():
 P=load('results/'+PRIMARY+'.json');x=list(map(F,P['solution']['exact_primal']));zi=len(x)-1;assert zi==922 and min(x)>=0
 connected={};cross={};label_to_key={}
 for r in P['connected_orbits']:
  g=R.graph(r['nodes'],[(a,b,m)for(a,b),m in zip(combinations(r['nodes'],2),r['labels'])]);key=R.tree_key(g,r['marked'])
  assert key not in connected and R.baseline(g,*r['marked'])==0 and R.type_of(g)==R.alias(r['type'])
  connected[key]=r['key'];label_to_key[r['key']]=('connected',key)
 for r in P['cross_orbits']:
  g=R.graph(r['nodes'],[(a,b,m)for(a,b),m in zip(combinations(r['nodes'],2),r['labels'])]);key=R.tree_key(g,(r['marked_root'],))
  assert R.type_of(g)=='E7'and key not in cross;cross[key]=r['key'];label_to_key[r['key']]=('cross',key)
 assert len(connected)==306 and len(cross)==4
 index={}
 for j,r in enumerate(P['variable_order']):
  kind,key=label_to_key[r['orbit']];context=r['A1_count_capped2'];k=kind,key,context
  assert k not in index and context in range(3)and(kind!='cross'or context==1);index[k]=j
 assert set(index)=={('connected',k,c)for k in connected for c in range(3)}|{('cross',k,1)for k in cross}
 raw=inputs();assert len(P['literal_maps'])==2397 and len(P['faces'])==1465
 matrices=[R.label_matrix(r['row'])for r in P['faces']];pool=defaultdict(list);comparisons=0
 for fi,(row,M)in enumerate(zip(P['faces'],matrices)):
  parameter=tuple(row['dotted_ratio'])if row['dotted_ratio']is not None else None;inv=R.face_invariant(M),parameter
  for N in pool[inv]:comparisons+=1;assert not R.face_iso(M,N)
  pool[inv].append(M)
 mapped=Counter()
 for(fam,i,row,parameter),mp in zip(raw,P['literal_maps']):
  assert(mp['family'],mp['index'])==(fam,i);fi=mp['representative'];target=P['faces'][fi]
  assert target['dotted_ratio']==parameter;R.check_bijection(R.label_matrix(row),matrices[fi],mp['to_representative']);mapped[fi]+=1
  assert dict(family=fam,index=i,to_representative=mp['to_representative'])in target['origins']
 assert len(mapped)==1465 and sum(len(f['origins'])for f in P['faces'])==2397
 unique={};payments=[];case=0;cross_instances=0;connected_seen=set();cross_seen=set()
 for fi,face in enumerate(P['faces']):
  cs=corners(face['row']);rank=face['row']['nodes']-2;assert rank<=29
  for c,s,sc in zip(cs,face['corners'],face['cross_corners']):
   assert c['marked']==s['marked']==sc['marked']and c['types']==sorted(R.alias(t)for t in s['core_forest'])and c['baseline']==F(s['baseline'])and c['A1']==sc['core_A1_count']
   if s['orbit']is not None:assert connected[c['key']]==s['orbit'];connected_seen.add(c['key'])
   if sc['cross_orbit']is not None:assert cross[c['cross']]==sc['cross_orbit'];cross_seen.add(c['cross'])
  sums=[]
  # Re-derive the feasible external-count cases using actual rank28/29 forests.
  feasible=[]
  for capped in range(3):
   okay=False
   for ambient in(28,29):
    remaining=ambient-rank
    for actual in range(max(0,remaining)+1):
     nonA1=remaining-actual
     if min(2,actual)==capped and(nonA1==0 or nonA1>=2):okay=True
   if okay:feasible.append(capped)
  assert feasible==list(range(min(2,29-rank)+1))
  for external in feasible:
   co=Counter();base=sum(c['baseline']for c in cs)
   for c in cs:
    context=min(2,c['A1']+external)
    if c['key']in connected and c['baseline']==0:co[index['connected',c['key'],context]]+=1
    if c['cross']in cross and c['A1']+external==1:co[index['cross',c['cross'],1]]+=1;cross_instances+=1
   val=base+sum(v*x[j]for j,v in co.items());assert val>=1;sums.append(dict(external_A1_count_capped2=external,payment=str(val)))
   expected=dict(face=fi,external_A1_count_capped2=external,baseline=str(base),coefficients=[[j,v]for j,v in sorted(co.items())]);assert expected==P['face_constraints'][case];case+=1
   unique.setdefault((tuple(sorted(co.items())),base-1),[]).append(dict(face=fi,origins=face['origins'],external_A1_count_capped2=external))
  payments.append(dict(face=fi,payments=sums))
 assert case==4212 and len(unique)==2925 and connected_seen==set(connected)and cross_seen==set(cross)
 rows=[];rhs=[];meta=[]
 for(co,b),w in sorted(unique.items()):rows.append({j:-v for j,v in co});rhs.append(b);meta.append(dict(kind='face',witnesses=w))
 totals={};baselines={};coefficients={};root_multiplicity=Counter()
 for name,g in R.MODELS.items():
  base=F(0);coef=[Counter()for _ in range(3)]
  for a,b in combinations(g,2):
   base+=R.baseline(g,a,b);key=R.tree_key(g,(a,b))
   if key in connected:
    assert R.baseline(g,a,b)==0
    for c in range(3):coef[c][index['connected',key,c]]+=1
  if name=='E7':
   for root in g:
    key=R.tree_key(g,(root,))
    if key in cross:coef[1][index['cross',key,1]]+=1;root_multiplicity[key]+=1
  totals[name]=[base+sum(v*x[j]for j,v in row.items())for row in coef];baselines[name]=base;coefficients[name]=coef
 for r in P['cross_orbits']:
  kind,key=label_to_key[r['key']];assert kind=='cross'and root_multiplicity[key]==r['node_multiplicity']and r['variable']==index['cross',key,1]
 assert set(root_multiplicity)==set(cross)
 names=[r['type']for r in P['exact_certificate']['component_totals']];assert len(names)==len(set(names))==92 and set(names)==set(R.MODELS)
 for r in P['exact_certificate']['component_totals']:
  name=r['type'];assert r['rank']==len(R.MODELS[name])and F(r['baseline'])==baselines[name]and list(map(F,r['totals']))==totals[name]
 for cut in P['solution']['forest_constraints']:
  c=cut['context'];forest=[names[i]for i in cut['components']];assert sum(len(R.MODELS[t])for t in forest)==28 and min(2,forest.count('A1'))==c
  co=Counter({zi:-1});b=F(0)
  for name in forest:co.update(coefficients[name][c]);b-=baselines[name]
  rows.append(dict(co));rhs.append(b);meta.append(cut)
 slacks=[b-sum(v*x[j]for j,v in row.items())for row,b in zip(rows,rhs)];assert min(slacks)>=0
 y=[F(0)]*len(rows)
 for d in P['solution']['nonzero_dual']:
  i=d['row'];assert y[i]==0 and d['constraint']==meta[i];y[i]=F(d['multiplier']);assert y[i]<0
 cols=[F(0)]*(zi+1)
 for mult,row in zip(y,rows):
  for j,a in row.items():cols[j]+=mult*a
 assert all(v<=int(j==zi)for j,v in enumerate(cols))
 dual=sum(v*b for v,b in zip(y,rhs));assert dual==x[-1]==F(P['solution']['optimum'])==F(230,3)
 dp=R.forest_dp(totals);forests=[]
 for n,fronts in[(28,P['exact_certificate']['fronts']),(29,P['rank29_same_rule']['fronts'])]:
  for r in fronts:
   c=r['context'];v=dp[c][n,c];assert v==F(r['maximum']);forest=r['maximizing_forest']
   assert sum(len(R.MODELS[t])for t in forest)==n and min(2,forest.count('A1'))==c and sum(totals[t][c]for t in forest)==v
   if n==28:assert r['rank_maxima']==[str(dp[c][k,c])if(k,c)in dp[c]else None for k in range(29)]
   forests.append(dict(rank=n,context=c,maximum=str(v),witness=forest))
 assert max(dp[c][28,c]for c in range(3))==dual and max(dp[c][29,c]for c in range(3))==84
 paths=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json','scripts/fresh_sqrt5_enlarged_a1_weights_20261015_review.py','results/fresh_sqrt5_enlarged_a1_weights_20261015_review.json','results/fresh_sqrt5_parameter_weight_stress_20261015.json','results/fresh_sqrt5_subdivision_families_20261015.json','results/fresh_sqrt5_rooted_component_completion_20261015.json','results/fresh_sqrt5_dotted_rank29_exhaustion_20261015.json','results/fresh_sqrt5_dotted_parameters_20261015_review.json','results/fresh_sqrt5_dotted_rank29_exhaustion_20261015_review.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in load(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_SEPARATE_UNIQUE_A1_E7_GRAPH_DUAL_DP_REPLAY',sources={p:sha(p)for p in paths},literal_inputs=2397,distinct_parameter_marked_graphs=1465,marked_graph_comparisons=comparisons,face_cases=4212,distinct_constraints=2925,connected_orbits=306,cross_orbits=4,cross_corner_occurrences=cross_instances,variables=922,finite_component_types=92,primal=str(x[-1]),dual=str(dual),minimum_primal_slack=str(min(slacks)),face_payments=payments,forest_maxima=forests,cross_root_multiplicities=[dict(key=cross[k],multiplicity=v)for k,v in sorted(root_multiplicity.items())],
  proof='A unique entireA1 gives exactlyone disconnected pair with each E7 node, so rooted-node multiplicities add to the E7 total only in context1. Both other contexts have zero cross contribution. Reverse-type exact forest DP therefore covers every forest. Dimension-feasible spectator cases independently reconstructed from possible residual ranks at28/29.',
  scope='Separate implementation using the earlier independent graph engine, no producer/optimizer imports. Sealing, arithmetic queue completeness and family classifications remain explicit bound dependencies. The family instances are conservative no-further-contact constraints, not a complete finite-extension theorem. No dimension bound. Same construction lane authorship is disclosed.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';out=json.loads(json.dumps(run()))
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['primal'],out['forest_maxima'])
