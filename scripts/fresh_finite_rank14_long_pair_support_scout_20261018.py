#!/usr/bin/env python3
"""Complete selected/literal rank14 orbit inventory, not a payment theorem."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d17_five_unit_reduction_20261018 as N
import fresh_small_d16_coupon_and_ideal17_rebalance_20261018 as CURRENT

ROOT=N.ROOT;T=N.T;A=N.A;P=N.P;Q=N.Q
TAG='fresh_finite_rank14_long_pair_support_scout_20261018'
IDEAL='fresh_order14_ambient62_ideal_margin_closure_20261018'
STEMS=tuple(dict.fromkeys((*N.STEMS,N.TAG,IDEAL,CURRENT.TAG)))


def key(stem,index):
    return stem,tuple(index) if isinstance(index,list) else index


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[N.TAG]['new_cost']=='120' and docs[N.TAG]['new_finite_coupons']==[]
    old={key(r['source'],r['case']):r for r in docs[CURRENT.TAG]['combined_literal_interface']['rows']}
    assert len(old)==216
    over=docs[A.OVER]
    patterns=[c for c in over['isolation']['cases'] if (c.get('component_order_bound') or 0)>=14 and len(c['marked_component'])<=14]
    assert patterns and len({c['index'] for c in patterns})==len(patterns)
    graphs={kind:T.target(14,kind) for kind in 'ABD'}
    literal=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=old[key(stem,index)];corners=[]
        assert len(cs)==len(prev['corners'])
        for (marks,unit,_),baseline in zip(cs,prev['corners']):
            assert marks==baseline['marks']
            ids=sorted(set(white)|set(marks));shared=None;images={}
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                shared=next((p for p in A.components(g,ids) if set(marks)<=set(p)),None)
            if shared and len(shared)==14:
                for kind,graph in graphs.items():
                    maps=T.T.mapped(g,shared,marks,graph)
                    if maps:images[kind]=sorted({tuple(m['marked_images']) for m in maps})
                assert len(images)==1
                assert baseline['retained_payment_lower_bound']=='1'
            corners.append(dict(marks=marks,baseline_payment=baseline['retained_payment_lower_bound'],whole_rank14_images=images))
        assert sum(F(c['baseline_payment']) for c in corners)==F(prev['total_lower_bound'])
        literal.append(dict(source=stem,case=index,corners=corners,baseline_total=prev['total_lower_bound']))
    assert len(literal)==216
    records=[]
    for kind,graph in graphs.items():
        adj,dist,_=P.distance_data(graph);autos=P.embeddings(graph,graph)
        pairs=set(combinations(range(14),2));assert len(pairs)==91
        basic={p for p in pairs if dist[p[0]][p[1]]<=8 or dist[p[0]][p[1]]==9 and any(len(adj[i])==1 for i in p)}
        long={p for p in pairs if dist[p[0]][p[1]]>=10};todo=set(long);orbits=[]
        while todo:
            p=min(todo);orb={tuple(sorted(mp[i] for i in p)) for mp in autos}
            assert orb<=long and not orb&basic;orbits.append(sorted(orb));todo-=orb
        selected=[];union=set();maps_count=0
        for c in patterns:
            core=over['isolation']['cores'][c['core_key']];g=A.matrix(core['gram'])
            maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],graph)
            images=sorted({tuple(m['marked_images']) for m in maps});union.update(images);maps_count+=len(maps)
            selected.append(dict(index=c['index'],method=c['method'],core_key=c['core_key'],
                component_order_bound=c['component_order_bound'],displayed_order=len(c['marked_component']),
                maps=maps,marked_images=images))
        orbit_rows=[]
        for orb in orbits:
            os=set(map(tuple,orb));hits=[]
            for r in selected:
                affected=[m for m in r['maps'] if tuple(m['marked_images']) in os]
                if affected:hits.append(dict(index=r['index'],method=r['method'],core_key=r['core_key'],
                    displayed_order=r['displayed_order'],affected_maps=affected))
            trials=[]
            for weight in (F(1,2),F()):
                rows=[]
                for r in literal:
                    debit=F();which=[]
                    for j,c in enumerate(r['corners']):
                        images=set(map(tuple,c['whole_rank14_images'].get(kind,[])))
                        if images&os:
                            assert images<=os;debit+=1-weight;which.append(j)
                    value=F(r['baseline_total'])-debit
                    if which:rows.append(dict(source=r['source'],case=r['case'],corners=which,
                        baseline_total=r['baseline_total'],debit=str(debit),hypothetical_total=str(value)))
                trials.append(dict(proposed_weight=str(weight),factor_saving=str(len(orb)*(1-weight)),
                    affected_literal_records=rows,minimum_total=str(min([F(r['baseline_total']) for r in literal if not any(key(r['source'],r['case'])==key(x['source'],x['case']) for x in rows)]+[F(r['hypothetical_total']) for r in rows])),
                    literal_violations=[r for r in rows if F(r['hypothetical_total'])<1],
                    geometric_coverage_proved=False))
            orbit_rows.append(dict(pairs=orb,distance=dist[orb[0][0]][orb[0][1]],
                selected_contexts=len(hits),selected_maps=sum(len(r['affected_maps']) for r in hits),
                method_counts=dict(Counter(r['method'] for r in hits)),
                displayed_order_counts=dict(Counter(str(r['displayed_order']) for r in hits)),
                complete_selected_obligations=hits,literal_trials=trials))
        assert maps_count>0
        records.append(dict(kind=kind,rank=14,current_cost='91',target=graph,automorphisms=autos,
            retained_basic_pairs=sorted(basic),complete_eligible_patterns=len(patterns),complete_maps=maps_count,
            selected_pair_union=sorted(union),unsupported_long_pairs=sorted(long-union),
            selected_patterns=selected,long_pair_orbits=orbit_rows))
    paths=[Path(CURRENT.__file__),Path(__file__),Path(N.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(P.__file__),Path(Q.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE_RANK14_SELECTED_LITERAL_SUPPORT_SCOUT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        records=records,current_literal_baseline=CURRENT.TAG,complete_literal_configurations=216,
        literal_interfaces=literal,new_weight_rule_proved=False,finite_DP_run=False,new_dimension_bound=None,
        unresolved_universal_scope='Every proposed deletion still needs arbitrary-support actual-component geometry and complete nonselected/later coverage. Selected absence and passing literal trials alone do not prove a transfer.',
        scope='Complete finite A/B/D14 selected labelled embeddings and all216 current sparse-D16/ideal17 literal interfaces; failed selected pattern searches retained. Orbit trials debit one proposed orbit at a time; no simultaneous theorem or finite optimization is inferred.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],[(r['kind'],r['complete_maps'],len(r['long_pair_orbits'])) for r in out['records']],flush=True)
