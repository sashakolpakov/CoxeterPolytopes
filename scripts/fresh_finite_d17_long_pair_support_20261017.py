#!/usr/bin/env python3
"""Complete finite A17/B17/D17 long-pair payment support and local deletion.

Frozen inputs supply the all-label geometric census and arbitrary-support
sealing proofs. This checker exhausts the interfaces affected by the deletion.
It does not calculate any finite-forest maximum or ambient dimension bound.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_a17_endpoint_payment_support_20261017 as A
import fresh_ideal7_finite_cross_payments_20261017 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite_d17_long_pair_support_20261017'
PARTIAL = 'fresh_d17_partial_weight_transfer_20261017'
NEW_CYCLE = 'fresh_cycle12_ideal5_transfer_20261017'
NEW_D = 'fresh_affine_d11_leaf_pair_transfer_20261017'
NEW_B = 'fresh_order12_affine_b_half_transfer_20261017'
NEW_C = 'fresh_order12_affine_c_endpoint_reduction_20261017'


def labels(g, ids):
    out = [[0]*len(ids) for _ in ids]
    for i, j in combinations(range(len(ids)), 2):
        x, y = ids[i], ids[j]
        if not g[x][y]:
            continue
        assert g[x][y] < 0
        ratio = 4*g[x][y]**2/(g[x][x]*g[y][y])
        if ratio not in (1, 2):
            return None
        out[i][j] = out[j][i] = 3 if ratio == 1 else 4
    return out


def target(kind):
    edges = [(i, i+1, 3) for i in range(16)]
    if kind == 'B':
        edges[-1] = (15, 16, 4)
    if kind == 'D':
        edges = [(i, i+1, 3) for i in range(14)]+[(14,15,3),(14,16,3)]
    return P.labelled(17, edges)


def mapped(g, ids, marks, graph):
    p = labels(g, ids)
    if p is None:
        return []
    maps = P.embeddings(p, graph)
    mi = [ids.index(m) for m in marks]
    return [dict(images=m, marked_images=sorted(m[i] for i in mi)) for m in maps]


def selected(doc, graph):
    out = []
    for c in doc['isolation']['cases']:
        if c.get('component_order_bound') != 17:
            continue
        g = A.matrix(doc['isolation']['cores'][c['core_key']]['gram'])
        ids, marks = c['marked_component'], c['preferred_corner']
        assert len(A.components(g, ids)) == 1
        maps = mapped(g, ids, marks, graph)
        out.append(dict(index=c['index'], method=c['method'], core_key=c['core_key'],
            nodes=ids, marks=marks, normalized_native_labels=labels(g,ids), maps=maps))
    assert len(out) == 140
    return out


def corner_record(g, white, c, graph, unit, weight):
    ids = sorted(set(white)|set(c['marks']))
    parts = A.components(g, ids)
    shared = next((p for p in parts if set(c['marks']) <= set(p)), None)
    maps = mapped(g, shared, c['marks'], graph) if shared and len(shared) == 17 else []
    if maps:
        assert A.Q.inertia(A.sub(g, ids)) == [len(ids), 0, 0]
    return dict(marks=c['marks'], unit=unit, old_weight=str(weight),
                shared_nodes=shared, maps=maps)


def literals(data, graph):
    rows = []
    g17 = A.matrix(data[A.OVER]['isolation']['cores']['(25, 25, ((9, 9),))']['gram'])
    g16 = A.matrix(data[A.OLD16]['exceptional_D15']['gram'])
    for stem, g, cases in (
        (A.OLD16, g16, data[A.OLD16]['exceptional_D15']['cases']),
        (A.OLD17, g17, data[A.OLD17]['order17']['replacement_cases']),
    ):
        for r in cases:
            white = sorted(set(range(len(g)))-set(r['black']))
            corners = [corner_record(g, white, c, graph, int(c['finite'] and F(c['base']) == 1), F(c['weight']))
                       for c in r['corners']]
            rows.append(dict(source=stem, case=r['case'], corners=corners))
    for stem, key, unit_key in ((A.D13, 'extreme_equality', 'unit'),
                               (A.D14, 'equality', 'finite_small17_unit')):
        for ri, r in enumerate(data[stem][key]['rows']):
            g = A.d13_gram(r) if stem == A.D13 else A.matrix(r['gram'])
            corners = [corner_record(g, r['white'], c, graph, c[unit_key],
                                    F(c[unit_key])+F(c['half_coupon'])) for c in r['corners']]
            rows.append(dict(source=stem, case=ri, corners=corners))
    for r in data[PARTIAL]['cases']:
        g = A.matrix(r['gram'])
        corners = [corner_record(g, r['white'], c, graph, c['finite_unit'],
                    F(c['finite_unit'])+F(c['finite_half'])+F(c['ideal_weight'])) for c in r['corners']]
        rows.append(dict(source=PARTIAL, case=r['index'], corners=corners))
    assert len(rows) == 174
    return rows


def later_interfaces(data):
    rows = []
    for stem in (*A.SUCCESSORS, NEW_CYCLE, NEW_D, NEW_B, NEW_C):
        corners = []
        for path, c in A.all_dicts(data[stem]):
            if not ('components' in c and 'marks' in c and
                    ('unit' in c or 'finite_unit' in c)):
                continue
            unit = c.get('finite_unit', c.get('unit',0))
            if not unit:
                continue
            shared = next((p for p in c['components'] if set(c['marks']) <= set(p['nodes'])), None)
            if shared:
                typ = shared['type']
                if len(shared['nodes']) == 17:
                    assert stem == PARTIAL and typ == 'D17'
                else:
                    assert len(shared['nodes']) < 17
            else:
                if c.get('full_A1_D12_upgrade'):
                    typ = 'separate A1/D12 full coupon'
                else:
                    assert c.get('full_A1_D11_upgrade')
                    typ = 'separate A1/D11 full coupon'
            corners.append(dict(path=path, marks=c['marks'], type=typ))
        rows.append(dict(source=stem, typed_unit_corners=corners,
                         type_counts=dict(Counter(c['type'] for c in corners))))
    assert next(r for r in rows if r['source']==PARTIAL)['type_counts'] == {'D17':4}
    # The remaining B11 preferred finite components are declared by their
    # whole-component outcome, before any hypothetical affine completion.
    b11 = [c for c in data[NEW_B]['cases'] if c['outcome']=='whole_finite_B11_before_any_affine_completion']
    assert len(b11) == 2
    # Cycle12 arbitrary single-profile outputs are finite paths of length<=6;
    # its equality outputs use separate A11 coupons or ideal/finite paths<=4.
    single = data[NEW_CYCLE]['single_profiles']
    return dict(typed_interfaces=rows, inherited=A.successor_interfaces(data),
        new_B11_preferred_components=len(b11), cycle12_single_profiles=len(single),
        universal_scope='Cycle12 finite paths<=6, whole A1/A11 and A2/A11 coupons; D11 one-sided whole D12 or A1/D11; D11 two-sided finite corners of order16. No whole rank17 same-component debit.',
        C12_scope='No new finite alternative; complete provenance deletion only.')


def simultaneous_literal_totals(records):
    out = []
    for i, row in enumerate(records[0]['literal_contexts']):
        corners = []
        for ci, old in enumerate(row['corners']):
            debit = 0
            matched = []
            for rec in records:
                r = rec['literal_contexts'][i]
                assert (r['source'],r['case']) == (row['source'],row['case'])
                c = r['corners'][ci]
                assert (c['marks'],c['old_weight'],c['unit']) == (old['marks'],old['old_weight'],old['unit'])
                images = {tuple(m['marked_images']) for m in c['maps']}
                if images and c['unit']:
                    removed = set(map(tuple,rec['removed_pairs']))
                    assert images <= removed or not (images&removed)
                    matched.append(dict(type=rec['kind']+'17',marked_images=sorted(images)))
                    if images <= removed:
                        debit += 1
            assert debit <= 1
            retained = F(old['old_weight'])-debit
            assert retained >= 0
            corners.append(dict(marks=old['marks'], old_weight=old['old_weight'],
                debit=debit, retained_weight=str(retained), whole_rank17_matches=matched))
        total = sum(F(c['retained_weight']) for c in corners)
        assert total >= 1
        out.append(dict(source=row['source'],case=row['case'],corners=corners,retained_total=str(total)))
    affected = [r for r in out if any(c['debit'] for c in r['corners'])]
    assert len(affected) == 1
    assert affected[0]['source'] == A.OLD17 and affected[0]['case'] == 662
    assert affected[0]['retained_total'] == '2'
    return out


def run():
    stems = tuple(dict.fromkeys((A.TAG,*A.STEMS,NEW_CYCLE,NEW_D,NEW_B,NEW_C,
             'fresh_ideal5_weighted_dimension68_20261017')))
    data = {s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem, doc in data.items():
        assert doc['status'].startswith('PASS')
        for path, h in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == h, (stem,path)
    records = []
    for kind in 'ABD':
        graph = target(kind)
        adj, ds, _ = P.distance_data(graph)
        pairs = set(combinations(range(17), 2))
        basic = {p for p in pairs if ds[p[0]][p[1]] <= 8 or
            (ds[p[0]][p[1]] == 9 and any(len(adj[x]) == 1 for x in p))}
        high = selected(data[A.OVER], graph)
        high_pairs = {tuple(m['marked_images']) for r in high for m in r['maps']}
        literal = literals(data, graph)
        literal_pairs = {tuple(m['marked_images']) for r in literal for c in r['corners']
            if c['unit'] for m in c['maps']}
        required = basic|high_pairs|literal_pairs
        unused = pairs-required
        autos = P.embeddings(graph, graph)
        todo = set(unused); orbits = []
        while todo:
            p = min(todo)
            orbit = {tuple(sorted(m[i] for i in p)) for m in autos}
            assert orbit <= todo
            todo -= orbit
            orbits.append(dict(pairs=sorted(orbit), distance=ds[p[0]][p[1]]))
        removed = unused|({(0,16)} if kind=='A' else set())
        retained = pairs-removed
        assert basic <= retained and high_pairs <= retained
        assert all({tuple(sorted(m[x] for x in p)) for m in autos} <= removed for p in removed)
        rec = dict(kind=kind, graph=graph, basic_pairs=sorted(basic),
            automorphisms=autos,
            selected_high_pairs=sorted(high_pairs), selected_high_contexts=high,
            literal_unit_pairs=sorted(literal_pairs), literal_contexts=literal,
            conservatively_required_pairs=sorted(required),
            unsupported_pairs=sorted(unused), unsupported_orbits=orbits,
            prior_deleted_pairs=[(0,16)] if kind=='A' else [],
            removed_pairs=sorted(removed), retained_pairs=sorted(retained),
            old_cost=135 if kind=='A' else 136, new_cost=len(retained))
        records.append(rec)
        print(kind, 'basic',len(basic),'high contexts',sum(bool(r['maps']) for r in high),
            'maps',sum(len(r['maps']) for r in high),'high pairs',len(high_pairs),
            'literal pairs',len(literal_pairs),'unused',len(unused),'cost',len(required),flush=True)
        print('UNSUPPORTED',kind,sorted(unused),flush=True)
    assert [r['new_cost'] for r in records] == [114,105,127]
    totals = simultaneous_literal_totals(records)
    domain = A.inherited_domain(data)
    routing = data[A.TAG]['selected_support']['routing']
    assert sum(routing.values()) == 2614
    assert routing['whole_component_at_most16'] == 2391
    assert routing['finite_leaf9_or_ideal_feature5'] == 74
    assert routing['nine_context_adaptive_finite_payment'] == 9
    paths = [Path(__file__),Path(A.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in stems
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE_RANK17_COMPLETE_LONG_PAIR_LOCAL_DELETION',
        records=records, simultaneous_literal_totals=totals,
        inherited_domain=domain, later_interfaces=later_interfaces(data),
        selected_domain_routing=dict(whole_component_at_most16=2391,
            protected_leaf9=74,adaptive_leaf9_or_order16=9,eligible_bound17=140),
        new_payment_theorem=True, new_dimension_bound=None,
        cross_component_coupons_unchanged=True, ideal_weights_unchanged=True,
        finite_forest_maximum_recomputed=False,
        scope='Simultaneous whole finite A17/B17/D17 same-component long-pair deletion, preserving all short8/actual-leaf9 pairs and every coupon. Complete affected selected and literal interfaces reconstructed; geometric census, all-label classification and arbitrary-support sealing inherited. No ambient bound or finite-forest maximum.',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(p.read_text())
    else:
        p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
