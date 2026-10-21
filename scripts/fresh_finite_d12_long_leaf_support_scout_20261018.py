#!/usr/bin/env python3
"""Complete selected/literal D12 endpoint-to-fork-leaf support diagnostic.

No local theorem: arbitrary supporting-white completions and nonselected
payment use must still be proved, even if a literal trial is successful.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_rank14_long_pair_support_scout_20261018 as R
import fresh_a15_d15_ambient59_forest_closure_20261018 as C

ROOT = R.ROOT
TAG = 'fresh_finite_d12_long_leaf_support_scout_20261018'
T, A, P, Q = R.T, R.A, R.P, R.Q
STEMS = tuple(dict.fromkeys((*R.STEMS, C.TAG)))
PAIRS = {(0, 10), (0, 11)}


def run():
    docs = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem, doc in docs.items():
        assert doc['status'].startswith('PASS'), stem
        for path, digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (stem, path)
    over = docs[A.OVER]
    patterns = [c for c in over['isolation']['cases'] if
                (c.get('component_order_bound') or 0) >= 12 and len(c['marked_component']) <= 12]
    assert len(patterns) == len({c['index'] for c in patterns})
    target = T.target(12, 'D')
    adjacency, distances, _ = P.distance_data(target)
    autos = P.embeddings(target, target)
    assert {p for p in combinations(range(12), 2) if distances[p[0]][p[1]] >= 10} == PAIRS
    assert all({tuple(sorted(mp[i] for i in p)) for p in PAIRS} == PAIRS for mp in autos)
    selected, hits = [], []
    for c in patterns:
        core = over['isolation']['cores'][c['core_key']]
        g = A.matrix(core['gram'])
        maps = T.T.mapped(g, c['marked_component'], c['preferred_corner'], target)
        affected = [m for m in maps if tuple(m['marked_images']) in PAIRS]
        row = dict(index=c['index'], method=c['method'], core_key=c['core_key'],
                   component_order_bound=c['component_order_bound'], displayed_order=len(c['marked_component']),
                   maps=maps, affected_maps=affected)
        selected.append(row)
        if affected:
            hits.append(row)
    print('SELECTED', len(patterns), 'maps', sum(len(r['maps']) for r in selected),
          'contexts', len(hits), flush=True)
    current = {R.key(r['source'], r['case']): r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    assert len(current) == 216
    literal = []
    for stem, index, g, white, corners in T.literal_inputs(docs):
        prior = current[R.key(stem, index)]
        assert len(corners) == len(prior['corners'])
        entries = []
        for (marks, unit, unused), previous in zip(corners, prior['corners']):
            assert marks == previous['marks']
            ids = sorted(set(white) | set(marks))
            images = []
            if Q.inertia(A.sub(g, ids)) == [len(ids), 0, 0]:
                shared = next((part for part in A.components(g, ids) if set(marks) <= set(part)), None)
                if shared and len(shared) == 12:
                    mapped = T.T.mapped(g, shared, marks, target)
                    images = sorted({tuple(m['marked_images']) for m in mapped})
            changed = bool(set(images) & PAIRS)
            if changed:
                assert set(images) <= PAIRS and F(previous['retained_payment_lower_bound']) == 1
            entries.append(dict(marks=marks, current_lower=previous['retained_payment_lower_bound'],
                                whole_D12_images=images, candidate_debit=changed))
        total = sum(F(r['current_lower']) for r in entries)
        assert total == F(prior['total_lower_bound'])
        literal.append(dict(source=stem, case=index, corners=entries, current_total=str(total)))
    assert len(literal) == 216
    trials = []
    for weight in (F(1, 2), F()):
        rows = []
        for old in literal:
            affected = [i for i, r in enumerate(old['corners']) if r['candidate_debit']]
            if not affected:
                continue
            debit = len(affected)*(1-weight)
            total = F(old['current_total'])-debit
            rows.append(dict(source=old['source'], case=old['case'], corners=affected,
                             previous_total=old['current_total'], debit=str(debit), proposed_total=str(total)))
        trials.append(dict(proposed_pair_weight=str(weight), factor_saving=str(2*(1-weight)),
                           affected=rows, violations=[r for r in rows if F(r['proposed_total']) < 1],
                           geometric_coverage_proved=False))
    assert C.actual_base(('D', 12)) == 16*66
    witness = [('A', 1)]*7+[('D', 5)]*3+[('D', 12)]*3
    assert F(C.value(witness), 16) == F(3063, 8)
    paths = {Path(__file__), Path(R.__file__), Path(C.__file__), Path(T.__file__), Path(T.T.__file__),
             Path(A.__file__), Path(P.__file__), Path(Q.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS
                 for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return dict(status='PASS_D12_LONG_LEAF_SELECTED_LITERAL_SCOUT',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
                target=target, automorphisms=autos, proposed_pairs=sorted(PAIRS), current_D12_cost='66',
                selected_patterns=selected, complete_eligible_patterns=len(patterns),
                complete_maps=sum(len(r['maps']) for r in selected),
                affected_contexts=hits, affected_maps=sum(len(r['affected_maps']) for r in hits),
                method_counts=dict(Counter(r['method'] for r in hits)),
                displayed_order_counts=dict(Counter(str(r['displayed_order']) for r in hits)),
                current_literal_baseline=C.TAG, literal_interfaces=literal, literal_trials=trials,
                formal_witness=witness, formal_witness_current_cost='3063/8',
                hypothetical_witness_cost_with_two_zero_pairs='3015/8',
                new_weight_rule_proved=False, finite_DP_run=False, new_dimension_bound=None,
                scope='Complete selected labelled D12 embeddings and216 current literal interfaces. Failed searches retained. Literal trial failures are obligations, not counterexamples to all possible replacement rules; no universal completion or output-use proof is asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], 'literal violations', [len(r['violations']) for r in out['literal_trials']], flush=True)
