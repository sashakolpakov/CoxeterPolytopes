#!/usr/bin/env python3
"""Complete finite8/ideal5 local payments and the unchanged-envelope obstruction.

No improved dimension bound is asserted.
"""
import argparse,json
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import fresh_ideal5_combined_payments_20261017 as C
import fresh_ideal5_piecewise_dimension77_20261017 as B
ROOT=C.ROOT;S=C.S
TAG='fresh_finite8_ideal5_combined_payments_20261017'
OVERLAP='fresh_distance8_high_overlap_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
LOW='fresh_distance8_low_core_20261017'
HIGH='fresh_distance8_high_disjoint_20261017'

def fixed_frontier():
 rows=[]
 for n in range(70,78):
  d=n-1;shapes=set()
  for a in range(d//12+1):
   for b in range((d-12*a)//14+1):
    for c in range((d-12*a-14*b)//16+1):
     for e in range((d-12*a-14*b-16*c)//17+1):
      for p,q,low in product((0,13),(0,15),range(12)):
       high=d-12*a-14*b-16*c-17*e-p-q-low
       if high==0 or high>=18:shapes.add(tuple(sorted([12]*a+[14]*b+[16]*c+[17]*e+([p]if p else[])+([q]if q else[])+([low]if low else[])+([high]if high else[]))))
  profiles=[B.profile(rs)for rs in sorted(shapes)];smallest=min(profiles,key=lambda r:F(r['margin']))
  beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
  assert F(smallest['margin'])<0
  rows.append(dict(ambient_dimension=n,reduced_profiles=len(profiles),minimum_margin=smallest['margin'],minimizing_ranks=smallest['ranks'],finite8_slack=str(beta-8*n)))
 boundary=B.profile((12,)+(16,)*4);deficit=-F(boundary['margin'])
 assert deficit==F(16664,1085773)and boundary['link_vertices']==1085773
 return dict(diagnostic_minima=rows,boundary_at_n77=boundary,
   unnormalized_cusp_deficit=16664,finite8_slack_at_n77='68',strict_global_shape_constant_at_n77=380,
   required_normalized_single_rank12_cost_decrease=str(deficit),
   required_each_rank16_normalized_decrease_if_shared_equally=str(deficit/4),
   required_order13_pair_cost_decrease=str(13*deficit),required_order17_pair_cost_decrease_each=str(17*deficit/4),
   scope='Formal cusp-profile obstruction to the unchanged envelope only; no geometric realization or new nonexistence bound asserted.')

def run():
 stems=(OVERLAP,SMALL,LOW,HIGH,C.TAG,S.TAG,C.TRI,C.MIXED,C.FINITE,S.D13,S.D14,S.OLD16,S.P.TAG,
       'fresh_cycle14_ideal5_transfer_20261017','fresh_finite8_cycle_replacement_20261017',
       'fresh_distance8_adaptive_leaf_overlap_20261017')
 inputs={s:S.read(s)for s in stems}
 for s,d in inputs.items():
  assert d['status'].startswith('PASS')
  for p,h in d['sources'].items():assert C.digest(ROOT/p)==h,(s,p)
 ov=inputs[OVERLAP];new=inputs[SMALL];old=inputs[S.TAG]
 assert ov['census']['consistent_matchings']==82197 and len(ov['isolation']['cases'])==2614
 assert len(new['patterns'])==1153 and len(new['records'])==16
 assert inputs[LOW]['extension_closure']['component_order_bound']==12
 assert len(inputs[LOW]['extension_closure']['cases'])==66
 assert inputs[C.FINITE]['dimension_maximum']==80
 assert all(F(r['excess'])<=0 for r in inputs[C.FINITE]['dimension_endpoints'])
 triangle=[];short_reduction=[]
 for newrow,oldrow in zip(new['records'],old['records']):
  assert(newrow['order'],newrow['kind'])==(oldrow['order'],oldrow['kind'])
  assert newrow['retained_paid_pairs']==oldrow['paid_pairs']
  assert newrow['retained_cost']==oldrow['cost']
  target=oldrow['target'];adj,ds,features=S.P.distance_data(target);paid=set(map(tuple,oldrow['paid_pairs']))
  near={i for i in range(len(target))if any(ds[i][f]<=5 for f in features)}
  required={p for p in combinations(range(len(target)),2)if ds[p[0]][p[1]]<=5 or(ds[p[0]][p[1]]<=8 and set(p)&near)}
  assert required<=paid
  missing=[p for p in combinations(range(len(target)),2)if ds[p[0]][p[1]]<=8 and p not in paid]
  assert not missing or oldrow['kind']=='A'
  triangle.append(dict(order=oldrow['order'],kind=oldrow['kind'],required_pairs=[list(p)for p in sorted(required)],retained=True))
  short_reduction.append(dict(order=oldrow['order'],kind=oldrow['kind'],unpaid_short8_pairs=[list(p)for p in missing],
      transfer='universal cycle14'if oldrow['order']==14 else'universal cycle>=15'if missing else'none'))
 # Whole-small alternatives and literal terminal coupons remain finite8-safe.
 assert {r['order']for r in inputs[S.D13]['formulas']['whole_affine12_payments']}=={12}
 assert {r['order']for r in inputs[S.D14]['profile_formulas']['small_affine_rows']}=={13}
 oldcases=inputs[S.OLD16]['exceptional_D15']['cases']+inputs[S.P.TAG]['order17']['replacement_cases']
 assert len(oldcases)==18
 for r in oldcases:
  assert sum(F(c['weight'])for c in r['corners']if c['finite'])>=1
  assert all(F(c['weight'])==0 for c in r['corners']if not c['finite'])
  for c in r['corners']:
   if c['finite'] and F(c.get('base','0'))>0:
    assert c['distance'] is not None and (c['distance']<=8 or max(c['component_orders'])<=17)
 equality_count=0
 for stem,ek,uk in ((S.D13,'extreme_equality','unit'),(S.D14,'equality','finite_small17_unit')):
  for row in inputs[stem][ek]['rows']:
   equality_count+=1;payment=F(0)
   for c in row['corners']:
    finite=c['signature'][1:]==[0,0]
    unit=any(len(b['marks'])==2 and len(b['nodes'])<=17 for b in c['components'])if finite else False
    assert int(unit)==c[uk]
    payment+=int(unit)+F(c['half_coupon'])
   assert payment>=1
 assert equality_count==144
 costs={}
 for r in new['records']:costs.setdefault(str(r['order']),{})[r['kind']]=r['retained_cost']
 assert costs==inputs[C.TAG]['small_order_costs']
 paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in stems for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
 paths += [ROOT/'scripts'/f'{B.TAG}.py',ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
 return dict(status='PASS_COMPLETE_FINITE8_IDEAL5_LOCAL_PAYMENT_AND_BUDGET',sources={str(p.relative_to(ROOT)):C.digest(p)for p in paths},
   finite_vertex_budget='8n',finite_budget_dimension_maximum=80,triangle_minimum='2',quadrilateral_minimum='1',
   high_overlap_contexts=2614,low_overlap_contexts=66,low_component_order_bound=12,selected_small_patterns=1153,
   small_order_costs=costs,ideal_rank_envelope=C.ENVELOPE,triangle_pair_retention=triangle,short8_reduction=short_reduction,
   full_new_overlap_domain=True,old_order16_17_affected_keys_matched=True,fixed_envelope_frontier=fixed_frontier(),
   scope='Complete finite8/ideal5 face-payment theorem with the frozen finite coupons for dimensions<=80. The unchanged ideal envelope still obstructs the same scalar method at n77. No new global dimension bound.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run();path=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],flush=True)
