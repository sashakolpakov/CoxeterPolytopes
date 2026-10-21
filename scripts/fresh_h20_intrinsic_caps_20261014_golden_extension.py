#!/usr/bin/env python3
"""Exact half-vertex extension of the frozen intrinsic golden H20 seed."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from fresh_h20_intrinsic_caps_20261014_golden_seed import ROOT,ZERO,ONE,q,add,sub,mul,scale,sign,enc,sha
TAG='fresh_h20_intrinsic_caps_20261014_golden_extension'
SRC='results/fresh_h20_intrinsic_caps_20261014_golden_seed.json'
HELPER='scripts/fresh_h20_intrinsic_caps_20261014_golden_seed.py'
def dec(x):return tuple(map(F,x))
def run():
 src=json.loads((ROOT/SRC).read_text())
 assert src['status']=='PASS_EXACT_275_FACET_INTRINSIC_H20_GOLDEN_SEED_WITH_STRICT_ESCAPE'
 for p,h in src['sources'].items():assert sha(p)==h,p
 B=[[F(x)for x in r]for r in src['spatial_Gram']]
 def sd(a,b):
  out=ZERO
  for i,row in enumerate(B):
   for j,x in enumerate(row):
    if x:out=add(out,scale(mul(a[i],b[j]),x))
  return out
 def ld(a,b):return sub(sub(sd(a,b),mul(a[19],b[20])),mul(a[20],b[19]))
 old=[[dec(x)for x in r]for r in src['walls']]
 centers=[r[:19]for r in old[23:]]
 norms=[sd(z,z)for z in centers]
 dots=[[sd(a,b)for b in centers]for a in centers]
 assert all(dots[i][i]==norms[i]for i in range(252))
 new=[[scale(x,F(1,2))for x in z]+[scale(sub(scale(n,F(1,4)),ONE),F(1,2)),ONE]for z,n in zip(centers,norms)]
 walls=old+new;assert len(walls)==527 and all(ld(r,r)==ONE for r in new)
 counts=Counter();mins={}
 def lower(key,val):
  if key not in mins or sign(sub(val,mins[key]))<0:mins[key]=val
 def distance(i,j,a=F(1),b=F(1)):
  return sub(add(scale(norms[i],a*a),scale(norms[j],b*b)),scale(dots[i][j],2*a*b))
 ids=src['source_wall_order'];tri={23,50,200}
 for i,r in enumerate(new):
  omitted=set(src['finite_endpoints'][i]['omitted_source_walls'])
  for j in range(22):
   p=ld(r,old[j]);b=old[j][19]
   if ids[j]not in omitted:
    assert p==scale(b,F(-1,2));counts['label5'if ids[j]in tri else'label3']+=1
   else:assert sign(add(p,ONE))<0;counts['dotted']+=1
  pc=scale(sub(q(5),scale(norms[i],F(1,4))),F(1,2))
  assert ld(r,old[22])==pc and sign(add(pc,q(2)))<0;counts['dotted']+=1
  for j in range(252):
   d=distance(i,j,F(1,2),F(1));lower('midpoint_to_vertex_square',d)
   assert sign(sub(d,q(4)))>0;counts['dotted']+=1
  for j in range(i):
   d=distance(i,j,F(1,2),F(1,2));lower('midpoint_separation_square',d)
   assert sign(sub(d,q(4)))>0;counts['dotted']+=1
 assert dict(counts)=={'label3':4284,'label5':504,'dotted':96138}
 assert sum(counts.values())==527*526//2-275*274//2
 center=[dec(x)for x in src['common_strict_center']['point']]
 cps=[ld(center,r)for r in new];assert all(sign(p)<0 for p in cps)
 # Existing actual finite endpoints and central facet remain feasible.
 for i,row in enumerate(src['finite_endpoints']):
  for j in range(252):assert sign(distance(i,j,F(1),F(1,2)))>0
 central=[dec(x)for x in src['central_facet_point']['point']]
 assert all(sign(ld(central,r))<0 for r in new)
 # A strict point on every new facet proves that none is redundant.
 facet_points=[]
 for i,r in enumerate(new):
  x=r[:19]+[add(r[19],ONE),ONE]
  assert ld(x,x)==q(-1)and sign(ld(x,center))<0 and ld(x,r)==ZERO
  assert all(sign(ld(x,w))<0 for w in old[:23])
  assert all(sign(distance(i,j,F(1,2),F(1)))>0 for j in range(252))
  assert all(sign(distance(i,j,F(1,2),F(1,2)))>0 for j in range(252)if i!=j)
  facet_points.append({'point':x,'norm':q(-1),'only_incident_wall':275+i})
 target=[dec(x)for x in src['strict_remaining_escape']['point']]
 targetpairs=[ld(target,r)for r in new]
 assert targetpairs[0]==ONE
 # A quarter-vertex point remains outside every hemisphere.
 z=[scale(x,F(1,4))for x in centers[0]];n=scale(norms[0],F(1,16))
 escape=z+[scale(sub(n,ONE),F(1,2)),ONE]
 assert ld(escape,escape)==ONE and sign(ld(escape,center))<0
 ep=[ld(escape,r)for r in walls];assert all(sign(p)<0 for p in ep)
 paths=[SRC,HELPER,str(Path(__file__).resolve().relative_to(ROOT))]
 full=Counter(src['pair_counts']);full.update(counts)
 return enc({'status':'PASS_EXACT_527_FACET_GOLDEN_EXTENSION_WITH_STRICT_ESCAPE','sources':{p:sha(p)for p in paths},'field':src['field'],'coordinates':src['coordinates'],'Lorentz_signature':[20,1],'spatial_Gram':B,'retained_wall_count':275,'new_walls':new,'new_norms':[ONE]*252,'new_pair_counts':dict(counts),'all527_pair_counts':dict(full),'minimum_distances':mins,'common_center':{'point':center,'norm':q(-9),'new_pairings':cps},'old252_actual_finite_endpoints_all_retained':True,'new_actual_facet_points':facet_points,'previous_escape_cut':{'point':target,'new_pairings':targetpairs,'wall275_pairing':ONE},'strict_remaining_escape':{'point':escape,'norm':ONE,'center_pairing':ld(escape,center),'all527_pairings':ep},'scope':'All275 old facets retained and252 new actual facets. Exact full Coxeter pair gate and252 retained finite endpoints. New caps cut the previous half-vertex escape, but a strict quarter-vertex escape proves infinite volume. No complete corner census or unrestricted obstruction claimed.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['new_pair_counts'])
