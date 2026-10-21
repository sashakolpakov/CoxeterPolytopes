#!/usr/bin/env python3
"""Complete provenance diagnostic for deleting only the finite A17 endpoint pair.

This is a local payment support audit. No new ambient bound is asserted.
Frozen certificates supply geometric completeness and actual-component sealing;
this checker enumerates the complete interfaces affected by the proposed change.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_selective_distance10_dimension116_20261017_review as Q

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite_a17_endpoint_payment_support_20261017'
OVER = 'fresh_distance8_high_overlap_20261017'
D13 = 'fresh_affine_d13_long_pair_transfer_20261017'
D14 = 'fresh_affine_d14_long_pair_transfer_20261017'
OLD16 = 'fresh_order16_finite_cross_payments_20261017'
OLD17 = 'fresh_ideal7_finite_cross_payments_20261017'
MIXED = 'fresh_finite8_ideal5_local_replacement_20261017'
ADAPT = 'fresh_distance8_adaptive_leaf_overlap_20261017'
SUCCESSORS = (
    'fresh_d17_singleton_pair_coupon_20261017',
    'fresh_d17_partial_weight_transfer_20261017',
    'fresh_d16_partial_weight_transfer_20261017',
    'fresh_order13_affine_b_singleton_transfer_20261017',
    'fresh_affine_d12_leaf_pair_transfer_20261017',
)
STEMS = (OVER, D13, D14, OLD16, OLD17, MIXED, ADAPT, *SUCCESSORS,
    'fresh_finite8_ideal5_combined_payments_20261017',
    'fresh_distance8_ideal5_small_factor_20261017',
    'fresh_distance8_low_core_20261017',
    'fresh_distance8_high_disjoint_20261017',
    'fresh_ideal5_finite8_triangle_payments_20261017',
    'fresh_finite8_cycle_replacement_20261017',
    'fresh_cycle14_ideal5_transfer_20261017',
    'fresh_cycle13_ideal5_transfer_20261017',
    'fresh_order15_affine_b_pair_reduction_20261017',
    'fresh_b14_b15_singleton_completion_transfer_20261017',
    'fresh_ideal5_weighted_dimension71_20261017',
    'fresh_ideal5_weighted_dimension69_20261017',
    'fresh_finite8_d12_b10_coupon_budget_20261017',
)


def matrix(x):
    return [[F(v) for v in row] for row in x]


def sub(g, ids):
    return [[g[i][j] for j in ids] for i in ids]


def components(g, ids):
    adj = {i: {j for j in ids if i != j and g[i][j]} for i in ids}
    unseen = set(ids)
    parts = []
    while unseen:
        seen = {min(unseen)}
        todo = list(seen)
        while todo:
            i = todo.pop()
            for j in adj[i]-seen:
                seen.add(j)
                todo.append(j)
        unseen -= seen
        parts.append(sorted(seen))
    return parts


def simple_path(g, ids):
    """Return a deterministic native-simple path, or None, in any positive gauge."""
    adj = {i: sorted(j for j in ids if i != j and g[i][j]) for i in ids}
    if any(g[i][j] >= 0 or 4*g[i][j]**2 != g[i][i]*g[j][j]
           for i in ids for j in adj[i]):
        return None
    if len(ids) == 1:
        return ids[:]
    ends = [i for i in ids if len(adj[i]) == 1]
    if len(ends) != 2 or any(len(adj[i]) > 2 for i in ids):
        return None
    path = [min(ends)]
    while True:
        nxt = [i for i in adj[path[-1]] if i not in path]
        if not nxt:
            break
        assert len(nxt) == 1
        path.append(nxt[0])
    return path if len(path) == len(ids) else None


def endpoint_debit(g, parts, marks):
    shared = next((ns for ns in parts if set(marks) <= set(ns)), None)
    if shared is None or len(shared) != 17:
        return 0, None
    path = simple_path(g, shared)
    if path is None:
        return 0, None
    positions = sorted(path.index(m) for m in marks)
    return int(positions == [0, 16]), dict(path=path, marked_positions=positions)


def all_dicts(value, path='root'):
    if isinstance(value, dict):
        yield path, value
        for key, val in value.items():
            yield from all_dicts(val, path+'.'+key)
    elif isinstance(value, list):
        for index, val in enumerate(value):
            yield from all_dicts(val, path+f'[{index}]')


def selected_support(over):
    rows = []
    routing = Counter()
    pair_union = set()
    cases = over['isolation']['cases']
    assert len(cases) == 2614
    assert over['census']['consistent_matchings'] == 82197
    for c in cases:
        method = c['method']
        bound = c.get('component_order_bound')
        if bound is None:
            assert method in ('finite_leaf9_or_ideal_feature5', 'nine_context_adaptive_finite_payment')
            if method == 'finite_leaf9_or_ideal_feature5':
                assert c['distance'] == 9 and c['protected_marked_leaves']
            routing[method] += 1
            continue
        assert bound <= 17
        if bound < 17:
            routing['whole_component_at_most16'] += 1
            continue
        core = over['isolation']['cores'][c['core_key']]
        g = matrix(core['gram'])
        ns, marks = c['marked_component'], c['preferred_corner']
        assert set(marks) <= set(ns)
        assert len(components(g, ns)) == 1
        path = simple_path(g, ns)
        maps = []
        if path is not None:
            assert Q.inertia(sub(g, ns)) == [len(ns), 0, 0]
            for start in range(18-len(path)):
                for ordered in (path, path[::-1]):
                    mp = {node: start+i for i, node in enumerate(ordered)}
                    # The complete embeddings of a connected path into a path
                    # are its consecutive intervals in the two orientations.
                    for i, j in combinations(ns, 2):
                        assert bool(g[i][j]) == (abs(mp[i]-mp[j]) == 1)
                    marked = sorted(mp[m] for m in marks)
                    assert marked != [0, 16]
                    pair_union.add(tuple(marked))
                    maps.append(dict(images=[mp[i] for i in ns], marked_images=marked))
        outcome = 'all_A17_maps_avoid_endpoint_pair' if maps else 'no_induced_A17_embedding'
        routing[outcome] += 1
        rows.append(dict(index=c['index'], method=method, component_order_bound=bound,
            core_key=c['core_key'], nodes=ns, marks=marks, path=path,
            outcome=outcome, maps=maps))
    assert len(rows) == 140 and sum(bool(r['maps']) for r in rows) == 15
    assert sum(len(r['maps']) for r in rows) == 48 and len(pair_union) == 10
    assert routing == Counter(whole_component_at_most16=2391,
        finite_leaf9_or_ideal_feature5=74, nine_context_adaptive_finite_payment=9,
        no_induced_A17_embedding=125, all_A17_maps_avoid_endpoint_pair=15)
    return dict(contexts=len(cases), routing=dict(routing), eligible_bound17_contexts=len(rows),
        path_patterns=15, complete_induced_maps=48, marked_pair_union=sorted(pair_union), cases=rows)


def legacy_replacements(inputs):
    over = inputs[OVER]
    g17 = matrix(over['isolation']['cores']['(25, 25, ((9, 9),))']['gram'])
    g16 = matrix(inputs[OLD16]['exceptional_D15']['gram'])
    rows = []
    for stem, g, saved in (
        (OLD16, g16, inputs[OLD16]['exceptional_D15']['cases']),
        (OLD17, g17, inputs[OLD17]['order17']['replacement_cases']),
    ):
        assert len(saved) == 9
        for r in saved:
            white = sorted(set(range(len(g)))-set(r['black']))
            corners = []
            for c in r['corners']:
                ids = sorted(set(white)|set(c['marks']))
                parts = components(g, ids)
                assert sorted(parts) == sorted(b['nodes'] for b in c['blocks'])
                sig = Q.inertia(sub(g, ids))
                assert c['finite'] == (sig == [len(ids), 0, 0])
                debit, a17 = endpoint_debit(g, parts, c['marks']) if c['finite'] else (0, None)
                assert not debit or F(c['base']) == 1
                retained = F(c['weight'])-debit
                assert retained >= 0
                shared = next((p for p in parts if set(c['marks']) <= set(p)), None)
                shared_type = None
                if c['finite'] and shared:
                    if simple_path(g, shared) is not None:
                        shared_type = 'A'+str(len(shared))
                    else:
                        adj = {i: [j for j in shared if i != j and g[i][j]] for i in shared}
                        forks = [i for i in shared if len(adj[i]) == 3]
                        assert len(forks) == 1
                        assert sum(len(adj[i]) == 1 for i in adj[forks[0]]) == 2
                        shared_type = 'D'+str(len(shared))
                corners.append(dict(marks=c['marks'], signature=sig,
                    old_weight=c['weight'], endpoint_debit=debit,
                    retained_weight=str(retained), shared_finite_component_type=shared_type,
                    A17_component=a17))
            old = sum(F(c['old_weight']) for c in corners)
            new = sum(F(c['retained_weight']) for c in corners)
            assert old >= 1 and new >= 1
            rows.append(dict(source=stem, case=r['case'], corners=corners,
                old_total=str(old), retained_total=str(new)))
    affected = [r for r in rows if r['old_total'] != r['retained_total']]
    assert len(affected) == 1 and affected[0]['case'] == 662
    assert affected[0]['old_total'] == '3' and affected[0]['retained_total'] == '2'
    assert [c['shared_finite_component_type'] for c in affected[0]['corners']
        if F(c['retained_weight'])] == ['D17', 'D17']
    return dict(cases=rows, affected_case=662, remaining_total_at_affected_case='2')


def d13_gram(row):
    g = [[F(2 if i == j else 0) for j in range(19)] for i in range(19)]
    L, R = row['leaf_attachments']
    edges = [(i, i+1) for i in range(9)]+[(0, 10), (0, 11), (9, 12), (9, 13)]
    edges += [(14, L), (15, R), (16, row['A3_attached_to']), (16, 17), (17, 18)]
    for i, j in edges:
        g[i][j] = g[j][i] = F(-1)
    return g


def universal_D_equalities(inputs):
    out = []
    a17 = []
    for stem, key, unit_key in ((D13, 'extreme_equality', 'unit'), (D14, 'equality', 'finite_small17_unit')):
        saved = inputs[stem][key]['rows']
        assert len(saved) == 72
        for ri, r in enumerate(saved):
            g = d13_gram(r) if stem == D13 else matrix(r['gram'])
            assert Q.inertia(g) == r['signature']
            corners = []
            for ci, c in enumerate(r['corners']):
                ids = sorted(set(r['white'])|set(c['marks']))
                parts = components(g, ids)
                assert sorted(parts) == sorted(p['nodes'] for p in c['components'])
                assert Q.inertia(sub(g, ids)) == c['signature']
                debit, path = endpoint_debit(g, parts, c['marks']) if c[unit_key] else (0, None)
                assert not debit
                if path:
                    a17.append(dict(source=stem, row=ri, corner=ci, marks=c['marks'], **path))
                corners.append(dict(marks=c['marks'], unit=c[unit_key], half=c['half_coupon'],
                    endpoint_debit=debit))
            total = sum(F(c['unit'])+F(c['half']) for c in corners)
            assert total >= 1
            out.append(dict(source=stem, row=ri, corners=corners, retained_total=str(total)))
    assert len(out) == 144 and len(a17) == 8
    assert {tuple(r['marked_positions']) for r in a17} <= {(0, 13), (3, 16)}
    return dict(rows=out, A17_corners=a17, endpoint_debits=0)


def successor_interfaces(inputs):
    rows = []
    for stem in SUCCESSORS:
        corners = [(path, c) for path, c in all_dicts(inputs[stem])
            if 'components' in c and 'marks' in c and ('finite_unit' in c or 'unit' in c)]
        assert corners
        units = []
        for path, c in corners:
            unit = c.get('finite_unit', c.get('unit', 0))
            if not unit:
                continue
            marked = [p for p in c['components'] if set(c['marks']) <= set(p['nodes'])]
            if marked:
                assert len(marked) == 1 and marked[0]['type'] != 'A17'
                typ = marked[0]['type']
            else:
                assert c['full_A1_D12_upgrade']
                assert sorted(p['type'] for p in c['components']) == ['A1', 'D12']
                typ = 'separate A1/D12 full coupon'
            units.append(dict(path=path, marks=c['marks'], component_type=typ))
        rows.append(dict(source=stem, explicitly_typed_corners=len(corners), finite_units=units,
            A17_endpoint_units=0))
    # The B15 transfer selects an already finite D14/A13 component, or excludes
    # all proper affine completions. Its checked outcome metadata is recorded.
    b15 = inputs['fresh_order15_affine_b_pair_reduction_20261017']
    assert len(b15['cases']) == 7
    b15_rows = []
    for r in b15['cases']:
        o = r['outcome']
        assert o['maximum_order'] <= 14 and r['finite_payment'] == 1
        assert o.get('whole_finite_type') == 'D14' or o.get('base_finite_type') == 'A12'
        b15_rows.append(dict(index=r['index'], maximum_whole_order=o['maximum_order'],
            outcome=o['method']))
    bzero = inputs['fresh_b14_b15_singleton_completion_transfer_20261017']
    assert all(r['map_outcomes'] == {'projection_exclusion': r['affected_map_count']}
        for r in bzero['records'])
    return dict(typed_unit_interfaces=rows,
        affine_B15_selected_cases=len(b15['cases']),
        affine_B15_whole_component_alternatives=b15_rows,
        affine_B14_B15_zero_transfer_maps=sum(r['affected_map_count'] for r in bzero['records']),
        scope='Explicit finite corner interfaces plus inherited universal branch proofs; no search over arbitrary spectator rank.')


def inherited_domain(inputs):
    low = inputs['fresh_distance8_low_core_20261017']['extension_closure']
    assert low['component_order_bound'] == 12 and len(low['cases']) == 66
    adaptive = []
    for r in inputs[ADAPT]['cases']:
        distances = [x['distance'] for _, x in all_dicts(r['payment']) if 'distance' in x]
        bounds = [x['whole_component_bound'] for _, x in all_dicts(r['payment']) if 'whole_component_bound' in x]
        assert all(d <= 9 for d in distances) and all(b <= 16 for b in bounds)
        adaptive.append(dict(index=r['index'], method=r['method'], retained_leaf_distances=distances,
            whole_component_bounds=bounds))
    eq = inputs[MIXED]['cycle18_equality']['corners']
    coupons = []
    for c in eq:
        if F(c['coupon']):
            assert c['kind'] == 'finite' and c['unit'] == 0 and c['marked_distance'] is None
            types = sorted(p['type'] for p in c['marked_components'])
            assert types == ['A1', 'A17'] and F(c['coupon']) == F(1, 2)
            a17 = next(p for p in c['marked_components'] if p['type'] == 'A17')
            assert a17['endpoint_distances'] == [5, 11]
            coupons.append(dict(marks=c['marks'], types=types, A17_positions=[5, 11], retained_coupon=c['coupon']))
    assert len(coupons) == 2
    return dict(low_overlap_contexts=66, low_component_bound=12, adaptive_contexts=adaptive,
        cycle18_cross_component_coupons=coupons,
        triangle_route='Finite distance at most8; ideal whole factors at most9; unchanged.',
        short_path_route='Finite8 directly. Universal cycles and trees use finite8, actual leaf9, or unchanged cross-component coupons.',
        high_disjoint_route='Complete inherited all-native-label exclusion, independent of whole A17 payment.',
        D13_D14_immediate_route='Whole affine orders12/13; later D12/B12 refinements checked by successor interfaces.',
        spectral_route='Cycle13/14/16/18 finite branches use distance at most8, actual leaf9, or separate-component coupons; none requires an endpoint-pair unit.')


def run():
    inputs = {}
    for stem in STEMS:
        d = json.loads((ROOT/f'results/{stem}.json').read_text())
        assert d['status'].startswith('PASS'), stem
        for path, expected in d['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == expected, (stem, path)
        inputs[stem] = d
    pairs = set(combinations(range(17), 2))
    basic = {p for p in pairs if p[1]-p[0] <= 8 or (p[1]-p[0] == 9 and (0 in p or 16 in p))}
    assert len(pairs) == 136 and len(basic) == 102
    drop = (0, 16)
    assert tuple(sorted(16-i for i in drop)) == drop and drop not in basic
    out = dict(status='PASS_FINITE_A17_ENDPOINT_COMPLETE_SUPPORT_DIAGNOSTIC',
        candidate_type='whole finite A17', canonical_path=list(range(17)),
        candidate_removed_pair=drop, old_same_component_cost=136, candidate_same_component_cost=135,
        short8_and_actual_leaf9_pairs=sorted(basic), basic_pair_count=102,
        automorphism_invariant=True, cross_component_coupons_unchanged=True,
        selected_support=selected_support(inputs[OVER]),
        legacy_replacements=legacy_replacements(inputs),
        universal_D_equalities=universal_D_equalities(inputs),
        inherited_domain=inherited_domain(inputs), successor_interfaces=successor_interfaces(inputs),
        new_ambient_dimension_bound=None, finite_forest_maximum_recomputed=False,
        scope='Complete affected local support/provenance audit of one A17 endpoint deletion. Inherited geometric classifications, all-label domain exhaustiveness and sealing are bound and reviewed, not rerun. No new global dimension bound or finite-forest optimum asserted.')
    paths = [Path(__file__), Path(Q.__file__)]
    paths += [ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    out['sources'] = {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}
    return json.loads(json.dumps(out))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], '15 patterns; 48 maps; legacy662 retains2; 144 D equalities retained', flush=True)
