#!/usr/bin/env python3
"""Complete literal screening for newly discounted quadrilateral branches.

This is a finite graph/provenance certificate. Arithmetic one-white decisions
are inherited from the separately reviewed all-8298 contact census, not
recomputed or promoted to arbitrary attachment closure here.
"""
import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import fresh_sqrt5_forest_sensitive_weights_20261015_review as G

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_adjacent_discount_frontier_20261016'
BASE = 'fresh_sqrt5_all_quad_unions_20261015'
CONTACT = 'fresh_sqrt5_all_quad_contacts_20261015'
ACCEPT = {'PASS_EXTRA_WHITE_NECESSARY_TESTS', 'DOTTED_PARAMETER_UNFILTERED'}


def read(tag):
    return json.loads((ROOT / f'results/{tag}.json').read_text())


def analysis(row):
    n = row['nodes']
    assert row['white'] == list(range(4, n))
    edges = dict(zip(combinations(range(n), 2), row['labels']))
    result = []
    for a, b in product((0, 2), (1, 3)):
        components = G.pieces(G.graph(row['white'] + [a, b], edges))
        c = next(c for c in components if a in c)
        family = G.family(c)
        connected = b in c
        base = int(G.base(c, a, b)) if connected else 0
        adjacent = c[a].get(b) == 3
        eligible = adjacent and len(c) >= 6
        # Complete finite-type extension fact: H3,H4,F4 are maximal under
        # connected finite extension. Small A/B/D components extend to rank6.
        potential = adjacent and (len(c) >= 6 or family.startswith(('A', 'B', 'D'))
                                  or family == 'I2(3)')
        assert not adjacent or base == 1
        assert not eligible or potential
        result.append(dict(marked=[a, b], baseline=base, adjacent3=adjacent,
                           eligible_now=eligible, eligible_after_finite_extension=potential,
                           marked_component=sorted(c), component_type=family,
                           component_rank=len(c),
                           corner_forest=sorted(G.family(t) for t in components),
                           distance=G.distances(c, a).get(b)))
    B = sum(c['baseline'] for c in result)
    current = sum(c['eligible_now'] for c in result)
    potential = sum(c['eligible_after_finite_extension'] for c in result)
    affected = (B == 1 and potential == 1) or (B == 2 and potential == 2)
    actual = (B == 1 and current == 1) or (B == 2 and current == 2)
    return dict(baseline=B, adjacent3_count=sum(c['adjacent3'] for c in result),
                eligible_now=current, eligible_after_finite_extension=potential,
                potentially_affected=affected, affected_if_components_entire=actual,
                corners=result)


def extend(row, contacts):
    n = row['nodes']
    assert len(contacts) == n
    edges = dict(zip(combinations(range(n), 2), row['labels']))
    edges.update({(i, n): m for i, m in enumerate(contacts)})
    return dict(nodes=n+1, white=list(range(4, n+1)),
                labels=[edges[p] for p in combinations(range(n+1), 2)])


def transition_check(parent, child):
    assert child['baseline'] >= parent['baseline']
    for p, c in zip(parent['corners'], child['corners']):
        assert p['marked'] == c['marked'] and p['adjacent3'] == c['adjacent3']
        # An existing path in a finite forest never changes under extension.
        if p['distance'] is not None:
            assert p['distance'] == c['distance']
        if p['component_type'] in ('H4', 'F4'):
            assert c['component_type'] == p['component_type']
        if p['component_type'] == 'H3':
            assert c['component_type'] in ('H3', 'H4')
        assert not c['eligible_after_finite_extension'] or p['eligible_after_finite_extension']
    # Universal safe pruning persists under every finite extension.
    if parent['baseline'] > 0 and not parent['potentially_affected']:
        assert not child['potentially_affected']


def generate():
    base, classes, domains = read(BASE), read(BASE+'_classes'), read(CONTACT)
    assert base['status'] == 'COMPLETE_ALL_QUAD_UNIONS'
    assert classes['marked_classes'] == 8298
    assert domains['status'] == 'COMPLETE_ALL_QUAD_ONE_WHITE_DOMAINS'
    assert len(domains['records']) == 8298
    rows = base['labelled_union_supports']
    counters, affected_status, transitions = Counter(), Counter(), Counter()
    records, immediate, outgoing = [], [], []
    all_transition_checks = 0
    for index, cl in enumerate(classes['marked_isomorphism_classes']):
        ui = cl['representative']
        row, old = rows[ui], domains['records'][index]
        assert old['class_index'] == index and old['union_index'] == ui
        a = analysis(row)
        counters[f"B{a['baseline']}_potential{a['eligible_after_finite_extension']}"] += 1
        trials = old.get('trials', [])
        accepted = old.get('accepted', [])
        assert sorted(accepted) == sorted(t['contacts'] for t in trials if t['status'] in ACCEPT)
        children = []
        for contact in accepted:
            child = extend(row, contact)
            ca = analysis(child)
            transition_check(a, ca)
            all_transition_checks += 1
            if a['potentially_affected']:
                reason = ('STILL_AFFECTED' if ca['potentially_affected'] else
                          'STABLY_SAFE_BASELINE')
                transitions[(old['dotted'], reason)] += 1
                children.append(dict(contacts=contact, reason=reason, analysis=ca))
                if ca['potentially_affected']:
                    outgoing.append(dict(parent_class=index, contacts=contact, row=child,
                                         parameter_type='dotted' if old['dotted'] else 'finite',
                                         analysis=ca))
        if a['potentially_affected']:
            affected_status[old['status']] += 1
            rec = dict(class_index=index, union_index=ui, row=row, analysis=a,
                       inherited_status=old['status'], singular=old.get('singular'),
                       domain_profiles=old.get('structural_profile_count', 0),
                       accepted_profiles=len(accepted), children=children)
            records.append(rec)
            if old['status'] == 'SEALED_SINGULAR_FINITE_CORE' and a['affected_if_components_entire']:
                immediate.append(dict(class_index=index, row=row, analysis=a))
    assert len(records) == 1216 and len(immediate) == 35
    paths = [f'scripts/{TAG}.py', f'scripts/{G.__name__}.py',
             f'results/{BASE}.json', f'results/{BASE}_classes.json',
             f'results/{BASE}_review_complete.json', f'results/{CONTACT}.json',
             f'results/{CONTACT}_review.json',
             'scripts/fresh_sqrt5_marked_discount_census_20261016.py',
             'scripts/fresh_sqrt5_marked_compensation_proposal_20261016.py']
    sources = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    return dict(status='PASS_COMPLETE_MINIMAL_UNION_ADJACENT_DISCOUNT_SCREEN',
                sources=sources, minimal_classes=8298, dimension_target=21,
                discount_hypotheses=dict(distance2='0', adjacent3='0 <= delta < 2/3',
                                         eligible='entire finite component rank >=6 and marked edge label3'),
                baseline_distribution=dict(counters), potentially_affected_classes=len(records),
                inherited_status_counts=dict(affected_status),
                all_accepted_extension_transition_checks=all_transition_checks,
                affected_transitions={f"{'dotted' if d else 'finite'}_{s}":v
                                      for (d,s),v in transitions.items()},
                immediate_sealed_weight_obligations=immediate,
                outgoing_literal_roots=len(outgoing), outgoing_roots=outgoing, records=records,
                scope='Complete screening of all8298 reviewed minimal unions and inherited first-white domains. The35 sealed singular rows require every retained facet mark and common spectator context. Outgoing finite and symbolic-dotted roots still need full all-attachment closure and exact payment. No new dimension bound, no claim that previous B0 queues cover these new branches, and no arithmetic contact replay here.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(generate()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == data
        print('PASS saved literal frontier replay')
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], data['potentially_affected_classes'],
          data['affected_transitions'], data['outgoing_literal_roots'])
