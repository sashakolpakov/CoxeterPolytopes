#!/usr/bin/env python3
"""Exact fixed-witness barrier for a stated inherited-charge extension class.

This does not build new columns or constrain arbitrary stronger inequalities.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw357_inherited_charge_barrier_20261018'
MODEL='fresh_all_marker_families_dual_continuation_20261018'
INCIDENCE='fresh_raw357_extra_vertex_incidence_pilot_20261018'


def run():
    paths={f'results/{MODEL}.json',f'scripts/{MODEL}.py',
        'results/fresh_all_marker_families_full_model_20261018.md',
        *(f'{folder}/{INCIDENCE}.{ext}'for folder,ext in
          [('scripts','py'),('results','json'),('results','md')]),
        str(Path(__file__).relative_to(ROOT))}
    model=json.loads((ROOT/f'results/{MODEL}.json').read_text())
    inc=json.loads((ROOT/f'results/{INCIDENCE}.json').read_text())
    for data in (model,inc):
        for p,h in data['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']}
    mu={v:F(x)for v,x in dual['vertex_multipliers']}
    assert all(x>=0 for x in alpha.values())and all(x>=0 for x in mu.values())
    num=F(dual['objective_numerator']);den=sum(mu.values(),F())
    assert den==F(dual['normalizer'])and num/den==F(dual['bound'])
    family=inc['complete_family_inventory'];domains=inc['extra_marker_inventory']
    assert len(family)==47 and len({r['face']for r in family})==47
    assert len(domains)==141
    source={r['face']:r for r in family}
    seen=set();counts=Counter();demand=Counter()
    for row in domains:
        f=row['face'];r=source[f];key=(f,tuple(row['omitted']))
        assert key not in seen;seen.add(key);counts[f]+=1
        assert row['actual_extra_occurrences']==r['required_quads']==1
        assert row['source_context']==r['source_context']
        assert r['raw_family']==357 and r['extra_vertex_occurrences']==3
        assert F(row['current_mass'])==F(r['current_mass'])==alpha.get(f,F())
        assert tuple(row['omitted'])in((2,4),(2,5),(2,6))
        demand[row['target_vertex']]+=alpha.get(f,F())
    assert set(counts)==set(source)and set(counts.values())=={3}
    assert len(demand)==138
    mass=sum((alpha.get(f,F())for f in source),F())
    positive=[[f,str(alpha[f])]for f in sorted(source)if alpha.get(f,F())]
    assert len(positive)==9 and mass==F(83284683,10**10)
    extra=sum(demand.values(),F());assert extra==3*mass
    extended=dict(mu)
    for v,e in demand.items():extended[v]=mu.get(v,F())+e
    newden=sum(extended.values(),F());assert newden==den+extra
    multi=[]
    for c in model['prepared']['inherited_multi_vertex_columns']:
        costs=c['integer_vertex_capacities'];assert costs and all(n>0 for v,n in costs)
        load=sum((alpha.get(f,F())*n for f,n in c['integer_face_coefficients']),F())
        old=sum((mu.get(v,F())*n for v,n in costs),F())
        new=sum((extended.get(v,F())*n for v,n in costs),F())
        assert load<=old<=new
        multi.append(dict(variable=c['variable'],load=str(load),old_cost=str(old),new_cost=str(new)))
    assert len(multi)==14
    ratio=num/newden;assert ratio>12
    return dict(status='PASS_RAW357_INHERITED_CHARGE_TEMPLATE_BARRIER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        hypotheses=[
          'Every new row has one unit of vertex budget at its target.',
          'After removing raw357 marker terms, its nonnegative face coefficient vector is dominated coordinatewise by a convex combination of existing normalized unit-budget columns at that same target.',
          'Every raw357 marker coefficient is nonnegative and at most one.',
          'Raw357 occurrence coefficients are exactly one at each of its three extra vertices; no orientation multiplier is used.',
          'All existing columns, including signed and multi-owner columns, remain unchanged.'],
        source_classes=47,pointed_domains=141,target_forests=138,
        positive_source_faces=positive,total_source_mass=str(mass),total_added_demand=str(extra),
        target_demand=[[v,str(e)]for v,e in sorted(demand.items())],
        numerator=str(num),old_normalizer=str(den),old_bound=str(num/den),
        dominating_normalizer=str(newden),feasible_dual_lower_bound=str(ratio),
        margin_above12=str(ratio-12),remaining_normalizer_gap_to12=str(num/12-newden),
        preserved_multi_owner_costs=multi,exported_columns=[],optimizer_invocations=0,
        scope='A conditional-template barrier: any extension satisfying all stated hypotheses retains this explicit >12 feasible dual. No new joint charge, actual least repair, arbitrary stronger-cut barrier, geometric realization or dimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists(),path
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['feasible_dual_lower_bound'],flush=True)
