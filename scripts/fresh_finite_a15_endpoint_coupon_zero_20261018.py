#!/usr/bin/env python3
"""A15 cost103, paying for two zeros by upgrading a whole A2/A13 coupon.

This local theorem adds two units per A2/A13 component pair. No forest
budget or ambient bound is inferred from the local identities.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_rank15_existing_half_transfer_20261018 as R
import fresh_finite_d17_a4_d13_half_transfer_20261018 as E

ROOT=R.ROOT; T=R.T; A=R.A; Q=R.Q; D=R.D
TAG='fresh_finite_a15_endpoint_coupon_zero_20261018'
DEPTH='fresh_depth8_finite8_zero_transfer_20261018'
ZERO={(0,13),(1,14)}
STEMS=tuple(dict.fromkeys((*R.STEMS,R.TAG,E.TAG,DEPTH)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    prior=docs[R.TAG]; latest=docs[E.TAG]
    assert latest['new_cost']=='118' and latest['all_other_pair_weights_unchanged']
    assert docs[DEPTH]['feature_depth']==8 and docs[DEPTH]['quadrilateral_minimum']==1
    rec=next(r for r in prior['rules'] if r['kind']=='A')
    old={tuple(r['pair']):F(r['weight']) for r in rec['weighted_pairs']}
    assert len(old)==105 and sum(old.values())==104 and all(old[p]==F(1,2) for p in ZERO)
    weights={p:(F() if p in ZERO else w) for p,w in old.items()}
    assert sum(weights.values())==103
    assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    for mp in rec['automorphisms']:
        assert {tuple(sorted(mp[i] for i in p)) for p in ZERO}==ZERO
    scout=next(r for r in docs[R.R.TAG]['records'] if r['kind']=='A')
    over=docs[A.OVER]
    support=R.validate_support(scout,over)
    assert support==dict(kind='A',patterns=536,maps=482,affected_cases=[1156,1616])
    coupon=next(c for c in docs[R.S.TAG]['coupons'] if c['name']=='A2_A13end')
    assert coupon['small_nodes']==[0,1] and coupon['large_nodes']==[0,12]
    assert coupon['weight']=='1/2' and coupon['pair_cost']=='2'
    proofs=[]
    for i in support['affected_cases']:
        c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        proof=R.core_proof(c,core,'A',{})
        baseline=next(p for p in prior['selected_proofs'] if p['index']==i)
        assert enc(proof)==baseline
        selected=alternative=0
        for corner in proof['corners']:
            if corner['reason']=='selected_rank15_half':
                assert corner['payment']==F(1,2)
                corner['payment']=F();selected+=1
            elif corner['reason']=='existing_whole_A2_A13_endpoint_half':
                assert corner['payment']==F(1,2)
                corner['payment']=F(1)
                corner['reason']='upgraded_whole_A2_A13_endpoint_unit';alternative+=1
        assert selected==alternative==1 and sum(c['payment'] for c in proof['corners'])==1
        proofs.append(proof)
    previous={(r['source'],E.N.normalize(r['case'])):r for r in latest['literal_rows']}
    assert len(previous)==216
    literal=[];hits=0
    for stem,index,g,white,corners in T.literal_inputs(docs):
        base=previous[stem,E.N.normalize(index)];rows=[]
        assert len(corners)==len(base['corners'])
        for j,(marks,*_) in enumerate(corners):
            assert marks==base['corners'][j]['marks']
            ids=sorted(set(white)|set(marks));images=[];debit=F()
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                b=next((b for b in D.finite_blocks(g,ids) if set(marks)<=set(b['nodes'])),None)
                if b and b['type']=='A15':
                    images=T.mark_pairs(g,b['nodes'],marks,15,'A')
                    debits={old[p]-weights[p] for p in images};assert len(debits)==1
                    debit=debits.pop()
            value=F(base['corners'][j]['retained_payment'])-debit
            assert value>=0;hits+=bool(debit)
            rows.append(dict(marks=marks,whole_A15_images=images,debit=debit,
                             retained_payment_lower_bound=value))
        total=sum(r['retained_payment_lower_bound'] for r in rows);assert total>=1
        literal.append(dict(source=stem,case=index,corners=rows,total_lower_bound=total))
    assert len(literal)==216 and hits==0
    # Positive replacement contributions in the other recent finite theorems
    # cannot spend one of the new A15 zeros. The coupon increase only helps.
    other=[p for p in prior['selected_proofs'] if p['kind']!='A']
    other+=docs[R.E.TAG]['selected_proofs']+latest['selected_proofs']
    for proof in other:
        for c in proof['corners']:
            if F(c['payment']):
                assert not any(b['type']=='A15' and set(c['marks'])<=set(b['nodes'])
                               for b in c['components'])
    paths=[Path(__file__),Path(R.__file__),Path(E.__file__),Path(T.__file__),
           Path(A.__file__),Path(Q.__file__),Path(D.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return enc(dict(status='PASS_FINITE_A15_ENDPOINT_COUPON_ZERO_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        target=rec['target'],automorphisms=rec['automorphisms'],zero_pairs=sorted(ZERO),
        old_cost='104',new_cost='103',weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        retained_basic_pairs=rec['retained_basic_pairs'],complete_support_validation=support,
        selected_proofs=proofs,literal_rows=literal,literal_configurations=216,literal_debited_corners=hits,
        upgraded_coupon=dict(name='A2_A13end',small='whole A2',small_marks=[0,1],
            large='whole A13',large_marks=[0,12],old_weight='1/2',new_weight='1',
            marked_pairs_per_component_pair=4,old_component_pair_cost='2',
            new_component_pair_cost='4',component_pair_cost_increase='2'),
        inherited_D17_cost='118',compatible_feature_depth=8,
        all_other_finite_and_ideal_weights_unchanged=True,triangle_minimum=2,quadrilateral_minimum=1,
        finite_budget_computed=False,new_dimension_bound=None,
        scope='Complete saved536-pattern/482-map validation, two inherited all-label radical seals reconstructed, latest216 literal interfaces after D17=118, and separate recent-transfer noninterference. The real extra cost is2 per whole A2/A13 component pair and must be included in any future forest bound. No new ambient theorem or finite-volume example.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'A15',result['new_cost'],'coupon_increase',
          result['upgraded_coupon']['component_pair_cost_increase'],flush=True)
