#!/usr/bin/env python3
"""Partial weight optimization allowing unique-A1 / rooted-E7 cross pairs.
The nonnegative baseline and all frozen2397 literal inputs remain in force.
"""
import argparse,json
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_parameter_weight_stress_20261015 as P

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_parameter_weight_stress_20261015_cross'
C=P.C;G=P.G
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def cross_corners(row,orbits):
 es=dict(zip(combinations(range(row['nodes']),2),row['labels']));out=[]
 for a,b in product((0,2),(1,3)):
  parts=G.pieces(G.graph(row['white']+[a,b],es));left=next(g for g in parts if a in g);right=next(g for g in parts if b in g);key=None
  count=sum(len(g)==1 for g in parts)
  if b not in left and count==1:
   component,root=(right,b)if len(left)==1 else(left,a)if len(right)==1 else(None,None)
   if component is not None and G.family(component)=='E7':
    key='A1xE7:'+json.dumps(G.code(component,(root,)),separators=(',',':'));ids=sorted(component)
    record=dict(key=key,type='E7',nodes=ids,labels=[component[u].get(v,2)for u,v in combinations(ids,2)],marked_root=root)
    if key not in orbits:orbits[key]=record
  out.append(dict(marked=[a,b],core_A1_count=count,cross_orbit=key))
 return out

def build():
 data=P.build();cross={};newfaces=[]
 for face in data['faces']:
  r=dict(face);r['cross_corners']=cross_corners(face['row'],cross);newfaces.append(r)
 data['faces']=newfaces;oldvars=len(data['variables']);keys=sorted(cross)
 data['variables']=data['variables']+[(k,1)for k in keys];index={k:oldvars+i for i,k in enumerate(keys)};data['zi']=len(data['variables'])
 unique={};constraints=[]
 for fc in data['faceconstraints']:
  fi=fc['face'];face=data['faces'][fi];external=fc['external_A1_count_capped2'];co=Counter(dict(fc['coefficients']))
  for c in face['cross_corners']:
   if c['cross_orbit']is not None and c['core_A1_count']+external==1:co[index[c['cross_orbit']]]+=1
  baseline=F(fc['baseline']);key=tuple(sorted(co.items())),baseline-1
  unique.setdefault(key,[]).append(dict(face=fi,origins=face['origins'],external_A1_count_capped2=external))
  constraints.append(dict(face=fi,external_A1_count_capped2=external,baseline=str(baseline),coefficients=[[j,v]for j,v in sorted(co.items())]))
 rows=[];rhs=[];meta=[]
 for(co,b),witnesses in sorted(unique.items()):rows.append({j:-v for j,v in co});rhs.append(b);meta.append(dict(kind='face',witnesses=witnesses))
 data.update(faceconstraints=constraints,rows=rows,rhs=rhs,meta=meta);components=deepcopy(data['components'])
 e7=next(g for name,g in G.models()if name=='E7');multiplicities=Counter()
 for root in e7:
  key='A1xE7:'+json.dumps(G.code(e7,(root,)),separators=(',',':'))
  if key in index:multiplicities[key]+=1
 assert set(multiplicities)==set(keys)
 for component in components:
  if component['type']=='E7':
   for key,m in multiplicities.items():component['coeffs'][1][index[key]]+=m
 data.update(components=components,cross_orbits=[dict(cross[k],node_multiplicity=multiplicities[k],variable=index[k])for k in keys],connected_variables=oldvars)
 return data

def run(saved=None):
 data=build();print('faces',len(data['faces']),'connectedvars',data['connected_variables'],'crossvars',len(data['cross_orbits']),'constraints',len(data['rows']),flush=True)
 solution=C.solve(data)if saved is None else saved['solution'];proof=C.verify(data,solution);rank29=P.rank29_check(data,solution)
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_parameter_weight_stress_20261015.py','results/fresh_sqrt5_parameter_weight_stress_20261015.json','scripts/fresh_sqrt5_a1_count_weights_20261015.py','results/fresh_sqrt5_enlarged_a1_weights_20261015_review.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 x=list(map(F,solution['exact_primal']))
 return dict(status='PASS_EXACT_PARTIAL_UNIQUE_A1_E7_CROSS_WEIGHT_OPTIMUM',sources={p:sha(p)for p in paths},literal_inputs=data['literal_count'],distinct_marked_face_graphs=len(data['faces']),literal_maps=data['literal_maps'],face_cases=len(data['faceconstraints']),distinct_constraints=len(data['rows']),connected_orbits=data['orbits'],cross_orbits=data['cross_orbits'],variables=data['zi'],variable_order=[dict(orbit=k,A1_count_capped2=c,kind='connected'if j<data['connected_variables']else'unique_A1_E7')for j,(k,c)in enumerate(data['variables'])],faces=data['faces'],face_constraints=data['faceconstraints'],solution=solution,exact_certificate=proof,target=77,meets_target=F(solution['optimum'])<=77,rank29_same_rule=rank29,cross_weights=[dict(orbit=r['key'],variable=r['variable'],node_multiplicity=r['node_multiplicity'],weight=str(x[r['variable']]))for r in data['cross_orbits']],
  rule='Existing connected-pair baseline/bonuses plus nonnegative weights on disconnected marked A1/E7 pairs only when the entire vertex has exactlyoneA1. Each such variable depends on the rooted E7 orbit. All other disconnected pairs have weight0.',
  separation_proof='With exactlyoneA1, every root of each E7 component forms exactlyone cross pair with that unique A1. Thus its cross contribution is additive over E7 components and equals the rooted-node multiplicity sum. Context0 and>=2 have no cross term. The original exact A1-count forest DP remains complete.',
  scope='Conservative partial weight system on the same2397 literal inputs as the preceding model. Family-instance constraints are not arbitrary-extension sealing proofs; finite-label attachment coverage remains incomplete. No dimension theorem.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';saved=json.loads(p.read_text())if a.verify else None;out=json.loads(json.dumps(run(saved)))
 if a.verify:assert out==saved
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['solution']['optimum'],'rank29',out['rank29_same_rule']['maximum'])
