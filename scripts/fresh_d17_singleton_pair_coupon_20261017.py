#!/usr/bin/env python3
"""Exact D17 payment-support diagnostic and one sealed A4/D13 transfer.

D17 denotes the affine diagram's order17; its standard affine type is D16.
This local package does not assert a global finite budget or dimension bound.
"""
import argparse,json,hashlib
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import fresh_ideal5_small_factor_payments_20261017 as S
import fresh_h17_half_coupon_transfer_20261017_review as H
P=S.P;ROOT=P.ROOT
TAG='fresh_d17_singleton_pair_coupon_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
OVERLAP='fresh_distance8_high_overlap_20261017'
OLDH='fresh_nonincident_carrier_incidence_20261017'

def key(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['match'])))
def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def run():
 x=read(SMALL);old=read(S.TAG);over=read(OVERLAP)
 rec=next(r for r in x['records']if r['order']==17 and r['kind']=='D')
 prior=next(r for r in old['records']if r['order']==17 and r['kind']=='D')
 target=rec['target'];autos=P.embeddings(target,target);assert len(autos)==8
 adj,ds,features=P.distance_data(target)
 paid=set(map(tuple,rec['retained_paid_pairs']));basic=set(map(tuple,prior['basic_pairs']))
 assert len(paid)==127 and len(basic)==110 and basic<=paid
 extra=paid-basic;orbits=[]
 while extra:
  p=min(extra);orbit={tuple(sorted(a[i]for i in p))for a in autos};assert orbit<=extra;extra-=orbit
  indices=[r['case']for r in rec['case_embeddings']if set(map(tuple,r['marked_images']))&orbit]
  cases=[]
  for i in indices:
   row=over['isolation']['cases'][i]
   cases.append({k:row[k]for k in('index','left','right','match','preferred_corner','component_order_bound','method')})
  orbits.append(dict(representative=list(p),distance=ds[p[0]][p[1]],pairs=[list(q)for q in sorted(orbit)],contexts=cases))
 assert [len(r['pairs'])for r in orbits]==[2,2,1,4,4,4]
 assert [len(r['contexts'])for r in orbits]==[2,2,1,4,4,16]
 single=next(r for r in orbits if len(r['pairs'])==1)
 assert single['pairs']==[[1,11]] and single['distance']==10 and len(single['contexts'])==1
 case=single['contexts'][0];assert case['index']==2518
 assert key(case)==(25,(3,4),25,(3,4),((9,9),))
 core=over['isolation']['cores'][over['isolation']['cases'][2518]['core_key']]
 assert core['gram']==[[str(v)for v in row]for row in H.G]
 assert H.inertia(H.G)==[17,1,1]
 assert all(sum(H.G[i][j]*H.Z[j]for j in H.N)==0 for i in H.N)
 affine=sorted(H.N-{3,13});sig=H.inertia(H.sub(affine));assert sig==[16,0,1]
 K={i:(1 if len(H.ADJ[i]&set(affine))==1 else 2)for i in affine}
 assert all(sum(H.G[i][j]*K[j]for j in affine)==0 for i in affine)
 assert H.Z[3]==-2 and H.Z[13]==2
 black={3,4,13,14};white=sorted(H.N-black);assert H.inertia(H.sub(white))==[15,0,0]
 # Every actual edge principal is strict, including the omitted black roles.
 for b in black:assert H.inertia(H.sub(sorted(set(white)|{b})))==[16,0,0]
 rows=[]
 for a,b in product((3,4),(13,14)):
  ids=sorted(set(white)|{a,b});signature=H.inertia(H.sub(ids));ca=H.component(ids,a);cb=H.component(ids,b)
  finite=signature==[17,0,0];parts=[]
  if not finite:
   assert[a,b]==[4,14]and ids==affine and signature==[16,0,1]
   assert ca[b]==10
   row=dict(marks=[a,b],signature=signature,kind='affineD16_order17',distance=10,removed_pair=[1,11],finite_unit=0,half_coupon='0')
  else:
   assert b not in ca
   for comp,mark in((ca,a),(cb,b)):
    name,role=H.name_and_role(comp,mark);nodes=sorted(comp);degree={i:len(H.ADJ[i]&set(nodes))for i in nodes}
    if name.startswith('A'):
     role='endpoint'if degree[mark]==1 else'other'
    else:
     fork=next(i for i in nodes if degree[i]==3)
     tail_neighbors=[i for i in H.ADJ[fork]&set(nodes)if degree[i]>1]
     assert len(tail_neighbors)==1
     if mark==tail_neighbors[0]:role='long_arm_neighbor_of_fork'
    parts.append(dict(nodes=nodes,type=name,mark=mark,role=role))
   coupon=F(0)
   for i,j in((0,1),(1,0)):
    if(parts[i]['type'],parts[i]['role'],parts[j]['type'],parts[j]['role'])==('A4','endpoint','D13','long_arm_neighbor_of_fork'):coupon=F(1,2)
   row=dict(marks=[a,b],signature=signature,kind='finite',components=parts,finite_unit=0,half_coupon=str(coupon))
  rows.append(row)
 assert [r['half_coupon']for r in rows]==['0','1/2','1/2','0']
 assert sum(F(r['half_coupon'])for r in rows)==1
 # Match the genuine quadrilateral in the independently frozen H17 chamber ledger.
 face=next(f for f in read(OLDH)['face_ledger']['two_faces']if f['polygon']==4 and f['side']==[3,4,13,14])
 assert face['white']==white and len(face['vertices'])==4 and not face['compact']
 assert sorted(c['marked']for c in face['corners'])==sorted(r['marks']for r in rows)
 # The singleton is fixed by every automorphism and occurs in no basic/triangle pair.
 reduced=paid-{(1,11)};assert len(reduced)==126 and basic<=reduced
 assert all({tuple(sorted(a[i]for i in p))for p in reduced}==reduced for a in autos)
 triangle=set()
 for p in combinations(range(17),2):
  if ds[p[0]][p[1]]<=5 or(ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features)):triangle.add(p)
 assert triangle<=reduced
 # Every newly unpaid emitted image has exactly the local core just certified.
 affected=[]
 for r in rec['case_embeddings']:
  if (1,11)in set(map(tuple,r['marked_images'])):affected.append(r['case'])
 assert affected==[2518]
 stems=[SMALL,OVERLAP,S.TAG,OLDH,'fresh_h17_half_coupon_transfer_20261017_review',
        'fresh_finite8_ideal5_combined_payments_20261017']
 paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in stems for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
 return dict(status='PASS_D17_SINGLETON_PAIR_LOCAL_TRANSFER',sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
   affine_type='D16',diagram_order=17,automorphisms=autos,basic_pair_count=110,old_cost=127,new_cost=126,
   extra_pair_orbits=orbits,removed_pair=[1,11],affected_overlap_cases=affected,
   literal_core=dict(gram=[[str(v)for v in row]for row in H.G],signature=[17,1,1],radical=H.Z,white=white,
      affine_nodes=affine,affine_positive_kernel=[K[i]for i in affine],omitted_radical_coefficients={'3':-2,'13':2},corners=rows),
   actual_h17_quadrilateral=face,new_paid_pairs=[list(p)for p in sorted(reduced)],triangle_pairs_retained=[list(p)for p in sorted(triangle)],
   new_finite_coupon=dict(weight='1/2',small_type='A4',small_mark='endpoint',large_type='D13',large_mark='unique long-arm neighbor of the fork center'),
   scope='Complete saved support/orbit diagnostic and one actual-star local transfer. Finite vertex budget and new global dimension comparison remain separate obligations.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run();path=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'six orbits, one core, two half-coupons; D17 cost127->126',flush=True)
