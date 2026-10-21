#!/usr/bin/env python3
"""Bounded exact node3 D17 alternate-corner and current-payment diagnostic.

No weight theorem, forest optimization, or dimension conclusion.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_finite_d17_a4_d13_half_transfer_20261018 as L
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as C
ROOT=L.ROOT;T=L.T;A=L.A;Q=L.Q;D=L.D
TAG='fresh_finite_d17_node3_alternate_support_scout_20261018'
BALANCE='fresh_center_reversion_rank16_balance_ambient60_frontier_20261018'
CHANGE=L.HALF
SPECS=[dict(name='A4end_A13fourth',small='A4',sm=[0,3],large='A13',lm=[3,9],weight='1/2',cost='2'),
       dict(name='A2_A14fourth',small='A2',sm=[0,1],large='A14',lm=[3,10],weight='1/2',cost='2'),
       dict(name='A4end_A13end',small='A4',sm=[0,3],large='A13',lm=[0,12],weight='1/2',cost='2'),
       dict(name='A1_A16fourth',small='A1',sm=[0],large='A16',lm=[3,12],weight='1/2',cost='1')]

def match(g,blocks,marks):
 a,b=[next(b for b in blocks if m in b['nodes'])for m in marks];out=[]
 if a is b:return out
 for aa,bb,m,n in((a,b,*marks),(b,a,marks[1],marks[0])):
  for s in SPECS:
   if aa['type']!=s['small']or bb['type']!=s['large']:continue
   x=T.mark_pairs(g,aa['nodes'],[m],len(aa['nodes']),'A');y=T.mark_pairs(g,bb['nodes'],[n],len(bb['nodes']),'A')
   if x==[(i,)for i in s['sm']]and y==[(i,)for i in s['lm']]:out.append(dict(spec=s['name'],small=aa,small_images=x,large=bb,large_images=y,increment=s['weight']))
 assert len(out)<=1
 return out

def run():
 stems=tuple(dict.fromkeys((*C.STEMS,C.TAG,L.TAG,L.A.OVER,BALANCE)))
 docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
 old=docs[L.TAG];over=docs[A.OVER]
 assert old['half_pairs']==[list(p)for p in sorted(CHANGE)]
 proofs=[]
 for i in L.CASES:
  c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
  proof=L.seal(c,core);assert L.enc(proof)==next(p for p in old['selected_proofs']if p['index']==i)
  g=A.matrix(proof['gram']);corners=[]
  for r in proof['corners']:
   matches=match(g,r['components'],r['marks']);debit=F(1,2)if r['selected']else F()
   payment=F(r['payment'])-debit+sum((F(x['increment'])for x in matches),F())
   corners.append(dict(marks=r['marks'],components=r['components'],old_payment=r['payment'],D17_debit=debit,new_coupons=matches,new_payment=payment))
  assert sum(r['new_payment']for r in corners)==1
  assert sum(bool(r['new_coupons'])for r in corners)==1
  proofs.append(dict(index=i,corners=corners,total='1',whole_seal_source=L.TAG))
 baseline={C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
 literal=[];hits=Counter();bad=[]
 for stem,index,g,white,cs in T.literal_inputs(docs):
  before=baseline[C.R.key(stem,index)];assert len(before['corners'])==len(cs);rows=[]
  for j,(marks,*_)in enumerate(cs):
   assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));delta=F();images=[];matches=[]
   if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
    blocks=D.finite_blocks(g,ids);shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
    if shared and shared['type']=='D17':
     images=T.mark_pairs(g,shared['nodes'],marks,17)
     if set(images)&CHANGE:assert set(images)<=CHANGE;delta-=F(1,2);hits['D17_debits']+=1
    matches=match(g,blocks,marks);delta+=sum((F(x['increment'])for x in matches),F());hits.update(x['spec']for x in matches)
   value=F(before['corners'][j]['retained_payment_lower_bound'])+delta
   rows.append(dict(marks=marks,images=images,coupons=matches,previous=before['corners'][j]['retained_payment_lower_bound'],delta=delta,payment=value))
  total=sum(r['payment']for r in rows);rec=dict(source=stem,case=index,corners=rows,total=total);literal.append(rec)
  if total<1 or any(r['payment']<0 for r in rows):bad.append(rec)
 assert len(literal)==216
 # Narrow leaf-only scenario: two already-retained whole A17 units replace
 # the unnecessary A4/A13 endpoint coupon, so only A1/A16 is new.
 leaf={(3,15),(3,16)}
 wa={tuple(r['pair']):F(r['weight'])for r in docs[BALANCE]['local_proofs']['finite_weighted_pairs']['A17']}
 assert wa[0,13]==wa[3,16]==1
 leafrows=[];leafbad=[];leafhits=Counter();selected=[]
 for p in proofs:
  if p['index']not in(1944,2032,2041,2337):continue
  oldp=next(z for z in old['selected_proofs']if z['index']==p['index']);g=A.matrix(oldp['gram']);corners=[]
  for z,c in zip(oldp['corners'],p['corners']):
   val=F(z['payment'])-(F(1,2)if z['selected']else F());im=[]
   block=next((b for b in z['components']if b['type']=='A17'and set(z['marks'])<=set(b['nodes'])),None)
   if block:
    im=T.mark_pairs(g,block['nodes'],z['marks'],17,'A');assert set(im)=={(0,13),(3,16)};val=max(val,min(wa[q]for q in im))
   cc=[x for x in c['new_coupons']if x['spec']=='A1_A16fourth'];val+=sum((F(x['increment'])for x in cc),F())
   corners.append(dict(marks=z['marks'],A17_images=im,new_coupons=cc,payment=val))
  assert sum(c['payment']for c in corners)>=1
  selected.append(dict(index=p['index'],corners=corners,total=sum(c['payment']for c in corners)))
 for stem,index,g,white,cs in T.literal_inputs(docs):
  before=baseline[C.R.key(stem,index)];rows=[]
  for j,(marks,*_)in enumerate(cs):
   ids=sorted(set(white)|set(marks));val=F(before['corners'][j]['retained_payment_lower_bound']);di=[];ai=[];credits=[]
   if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
    blocks=D.finite_blocks(g,ids);shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
    if shared and shared['type']=='D17':
     di=T.mark_pairs(g,shared['nodes'],marks,17)
     if set(di)&leaf:assert set(di)<=leaf;val-=F(1,2);leafhits['D17_debits']+=1
    if shared and shared['type']=='A17':
     ai=T.mark_pairs(g,shared['nodes'],marks,17,'A')
     if set(ai)=={(0,13),(3,16)}:
      if val<1:leafhits['retained_A17_unit_strengthening']+=1
      val=max(val,F(1))
    credits=[x for x in match(g,blocks,marks)if x['spec']=='A1_A16fourth']
    val+=sum((F(x['increment'])for x in credits),F());leafhits.update(x['spec']for x in credits)
   rows.append(dict(marks=marks,D17_images=di,A17_images=ai,new_coupons=credits,payment=val))
  rr=dict(source=stem,case=index,corners=rows,total=sum(c['payment']for c in rows));leafrows.append(rr)
  if rr['total']<1 or any(c['payment']<0 for c in rows):leafbad.append(rr)
 assert len(leafrows)==216
 rep=C.P.Z.B.B.representative;types=[t for t in C.P.Z.B.B.types()if t[1]<=59]
 def added(a,b):
  ans=F()
  for s in SPECS:
   if {a,b}=={(s['small'][0],int(s['small'][1:])),(s['large'][0],int(s['large'][1:]))}:ans+=F(s['cost'])
  return ans
 violations=[];costs=[]
 for a in types:
  for b in types:
   actual=F(C.actual_cross(a,b),16)+added(a,b);upper=F(C.upper_cross(rep(a),rep(b)),16)
   if actual>upper:violations.append(dict(types=[a,b],new_actual=actual,current_upper=upper,excess=actual-upper))
   leafactual=F(C.actual_cross(a,b),16)+(1 if {a,b}=={('A',1),('A',16)}else 0)
   assert leafactual<=upper
 for a,b in((('A',4),('A',13)),(('A',2),('A',14)),(('A',1),('A',16))):
  costs.append(dict(types=[a,b],old_actual=F(C.actual_cross(a,b),16),added=added(a,b),new_actual=F(C.actual_cross(a,b),16)+added(a,b),current_upper=F(C.upper_cross(rep(a),rep(b)),16)))
 paths={Path(__file__),Path(L.__file__),Path(C.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__)}
 paths.update(ROOT/f'{f}/{s}.{e}'for s in stems for f,e in(('scripts','py'),('results','json'),('results','md')))
 return L.enc(dict(status='PASS_D17_NODE3_ALTERNATE_SUPPORT_SCOUT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
  proposed_zero_pairs=sorted(CHANGE),proposed_coupons=SPECS,selected_proofs=proofs,literal_configurations=216,literal_hits=dict(hits),literal_rows=literal,literal_failures=bad,
  simultaneous_type_pairs=len(types)**2,cost_rows=costs,unchanged_upper_domination_failures=violations,
  current_D17_cost=str(F(C.actual_base(('D',17)),16)),proposed_D17_cost=str(F(C.actual_base(('D',17)),16)-2),
  leaf_only=dict(selected_proofs=selected,retained_whole_A17_unit_pairs=[[0,13],[3,16]],new_coupon=SPECS[-1],
   literal_rows=leafrows,literal_hits=dict(leafhits),literal_failures=leafbad,all_current_type_pairs_dominated=True,
   proposed_D17_cost=str(F(C.actual_base(('D',17)),16)-1)),
  geometric_theorem_proved=False,finite_forest_count_performed=False,new_dimension_bound=None,
  scope='Eight inherited whole-core seals and exact current216 literal/coupon costs only. Universal/later payment-use inventory and a changed full finite count remain required.'))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'literal',out['literal_hits'],'failures',len(out['literal_failures']),'costs',out['cost_rows'],'typefailures',len(out['unchanged_upper_domination_failures']),flush=True)
