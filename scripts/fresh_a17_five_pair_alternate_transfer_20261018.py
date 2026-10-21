#!/usr/bin/env python3
"""Promote five whole-A17 alternate-corner reductions with full interfaces.

The bound pilot supplies exact selected geometry/literal/type arithmetic;
the report proves arbitrary-support and complete current composition.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_a17_alternate_corner_pilot_20261018 as S
import fresh_finite_d17_leaf_a16_half_transfer_20261018 as LEAF
import fresh_d17_existing_three_corner_transfer_20261018 as NODE2
import fresh_order15_bd_basic_fractional_transfer_20261018 as IDEAL

ROOT=S.ROOT;X=S.X;T=S.T;A=S.A;D=S.D
TAG='fresh_a17_five_pair_alternate_transfer_20261018'
STEMS=tuple(dict.fromkeys((*S.STEMS,S.TAG,LEAF.TAG,NODE2.TAG,IDEAL.TAG)))


def actual_base(t):return X.actual_base(t)-(64 if t==('A',17) else 0)
def upper_base(t):return X.upper_base(t)-(64 if t==('A',17) else 0)


def actual_cross(a,b):
    pair={a,b}
    return X.actual_cross(a,b)+(64 if pair=={('D',5),('A',12)} else 32 if pair=={('A',1),('A',14)} else 0)


def positive_interfaces(docs):
    examined=[];hits=[]
    for stem,d in docs.items():
        for name in ('selected_proofs','radical_transfers','whole_lorentz_cases','local_proofs'):
            ps=d.get(name,[])
            if isinstance(ps,dict):ps=ps.get('proofs',[])
            if not isinstance(ps,list):continue
            for p in ps:
                if 'gram' not in p:continue
                g=A.matrix(p['gram'])
                for c in p.get('corners',[]):
                    pay=c.get('payment',c.get('used_new_payment',c.get('retained_payment_lower_bound','0')))
                    try:positive=F(pay)>0
                    except (TypeError,ValueError):continue
                    if not positive or 'components' not in c:continue
                    b=next((b for b in c['components'] if b['type']=='A17' and set(c['marks'])<=set(b['nodes'])),None)
                    if not b:continue
                    images=T.mark_pairs(g,b['nodes'],c['marks'],17,'A')
                    r=dict(source=stem,array=name,index=p.get('index'),marks=c['marks'],images=images,payment=str(pay));examined.append(r)
                    if set(images)&S.CHANGE.keys():hits.append(r)
    assert not hits,hits
    return dict(examined_positive_whole_A17_corners=examined,changed_later_positive_corners=hits,
        scope='Stored direct positive-corner arrays with physical Grams in the bound local packets. Supplements the complete written current-output inventory; not a new universal geometric census.')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    fresh=S.run();assert fresh==docs[S.TAG]
    assert not fresh['literal_changed'] and not fresh['literal_failures'] and not fresh['domination_failures']
    assert fresh['selected_support']['complete_eligible_patterns']==140 and fresh['selected_support']['complete_maps']==48
    assert len(fresh['selected_core_diagnostics'])==5
    over=docs[A.OVER]
    for p in fresh['selected_core_diagnostics']:
        c=over['isolation']['cases'][p['index']];core=over['isolation']['cores'][c['core_key']]
        assert core['gauge']==[1]*19 and core['gram']==p['gram']
        assert all(F(p['gram'][i][i])==2 for i in range(19))
    old={tuple(r['pair']):F(r['weight']) for r in docs[S.B.TAG]['local_proofs']['finite_weighted_pairs']['A17']}
    new={p:S.CHANGE.get(p,w) for p,w in old.items()};assert sum(new.values())==F(451,4)
    assert all(new[tuple(p)]==1 for p in docs[LEAF.TAG]['retained_whole_A17_unit_pairs'])
    assert set(S.CHANGE).isdisjoint(S.B.PAIRS['A17']) and new[1,15]==1
    assert docs[NODE2.TAG]['new_cost']=='927/8' and docs[IDEAL.TAG]['new_coupon']['canonical_large_node']==9
    assert actual_cross(('D',5),('A',12))==108 and actual_cross(('A',1),('A',14))==32
    types=[t for t in X.P.Z.B.B.types() if t[1]<=59];rep=X.P.Z.B.B.representative
    assert all(upper_base(rep(t))>=actual_base(t)>=0 for t in types)
    for a in types:
        for b in types:assert actual_cross(a,b)==actual_cross(b,a)<=X.upper_cross(rep(a),rep(b))
    paths={Path(__file__),Path(S.__file__),Path(X.__file__),Path(LEAF.__file__),Path(NODE2.__file__),Path(IDEAL.__file__)}
    paths.update(ROOT/f'{f}/{s}.{e}' for s in STEMS for f,e in (('scripts','py'),('results','json')))
    paths.update(ROOT/f'results/{s}.md' for s in STEMS if s!=S.TAG)
    return S.enc(dict(status='PASS_A17_FIVE_PAIR_ALTERNATE_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        changed_pairs=fresh['proposed_pairs'],old_cost='467/4',new_cost='451/4',component_saving='4',
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(new.items())],
        selected_support=fresh['selected_support'],selected_proofs=fresh['selected_core_diagnostics'],
        current_literal_interface=dict(configurations=216,changed=[],rows=fresh['literal_records']),
        positive_interface_audit=positive_interfaces(docs),
        retained_A17_unit_pairs=docs[LEAF.TAG]['retained_whole_A17_unit_pairs'],
        new_coupons=[dict(small='whole_D5_short_leaf',large='whole_A12_third_node',small_nodes=[3,4],large_nodes=[2,9],weight='1',actual_component_increment='4'),
            dict(small='whole_A1',large='whole_A14_third_or_fourth_node',large_nodes=[2,3,10,11],weight='1/2',actual_component_increment='2')],
        actual_component_costs=dict(D5_A12='27/4',A1_A14='2'),
        simultaneous_actual_type_pairs=len(types)**2,all_cross_upper_costs_unchanged=True,
        all_ideal_weights_unchanged=True,triangle_minimum=2,quadrilateral_minimum=1,
        finite_forest_count_performed=False,new_dimension_bound=None,
        scope='Fresh replay of complete140/48 selected pilot, five exact arbitrary-support radical seals and216 unchanged literal interfaces; current positive-array audit and explicit later/universal composition. Every new coupon charged with simultaneous typed domination. No new finite forest optimization or ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_cost'],'later_A17',len(out['positive_interface_audit']['examined_positive_whole_A17_corners']),flush=True)
