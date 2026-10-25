#!/usr/bin/env python3
"""SOURCE PREPARATION: lower-order/cycle terminal bindings, actual ambient58.

No legacy mathematical worker is imported. Fixed physical receipts, canonical
typed outputs and geometric/priority premises are reported separately. Root
must review and launch each named batch through the existing serial guard.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
import fcntl
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
ENGINE = 'scripts/finite57_terminal_literal_engine_small_affine_20261025.py'
ENGINE_SHA = '39dd55d2b2abef3db045e646edcc95bc1c8363f483d9f853d8c1d4e1b65d0e1f'
LIB = 'scripts/finite57_literal_weights_audit_20261024.py'
LIB_SHA = 'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0'
MANIFEST = 'results/finite57_g4_lower_cycle_source_inputs_20261025.json'
MANIFEST_SHA = '9f26c3f91c8c47126e8506e23d7a77497735974e012209a582da1344ddeb232b'
MAX_FILE = 4 * 1024**2
MAX_ROWS = 10000
STEMS = {
    'order10-cycle': 'fresh_order10_cycle_half_transfer_20261018',
    'order10-b': 'fresh_order10_b_terminal_half_transfer_20261018',
    'order10-d': 'fresh_order10_d_opposite_fork_fractional_transfer_20261018',
    'cycle11': 'fresh_cycle11_half_pair_transfer_20261018',
    'order11-bcd': 'fresh_order11_bcd_half_pair_transfer_20261018',
    'order11-bc': 'fresh_order11_bc_zero_transfer_20261018',
    'order12-center': 'fresh_order12_b_center_terminal_transfer_20261018',
    'order12-near': 'fresh_order12_bc_near_terminal_half_transfer_20261018',
    'order12-selected': 'fresh_order12_bc_unit_coupon_transfer_20261017',
    'order13-sparse': 'fresh_order13_bc_sparse_terminal_transfer_20261018',
    'cycle12': 'fresh_cycle12_ideal5_transfer_20261017',
    'cycle13': 'fresh_cycle13_ideal5_transfer_20261017',
    'cycle14': 'fresh_cycle14_ideal5_transfer_20261017',
    'cycle16': 'fresh_cycle16_half_coupon_transfer_20261017',
    'cycle18-tree': 'fresh_finite8_ideal5_local_replacement_20261017',
    'large-d': 'fresh_large_finite_d_leaf_half_transfer_20261017',
    **{f'native-corridors-{i}': 'fresh_large_finite_d_leaf_half_transfer_20261017' for i in range(4)},
}
PINS = {}; SIZES = {}; READ = {}; W = E = None; FUEL = [0]; ROWS = []; PRIORITIES = []


def need(condition, message):
    if not condition:
        raise ValueError(message)


def raw(relative, cap=MAX_FILE):
    path = (ROOT / relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded repository input path')
    with path.open('rb') as stream:
        data = stream.read(cap + 1)
    need(0 < len(data) <= cap, 'bounded input size: ' + relative)
    return data


def unique(pairs):
    result = {}
    for key, value in pairs:
        need(key not in result, 'duplicate JSON key')
        result[key] = value
    return result


def decode(data):
    return json.loads(data, object_pairs_hook=unique,
        parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON number')))


def captured(relative):
    data = raw(relative)
    need(relative in PINS and len(data) == SIZES[relative]
         and sha256(data).hexdigest() == PINS[relative], 'source identity: ' + relative)
    READ[relative] = PINS[relative]
    return data


def packet(batch):
    stem = STEMS[batch]
    for folder, ext in (('scripts', 'py'), ('results', 'md')):
        captured(f'{folder}/{stem}.{ext}')
    value = decode(captured(f'results/{stem}.json'))
    need(type(value) is dict and str(value.get('status', '')).startswith('PASS'), 'inherited packet schema')
    return value


def load(name, relative, data):
    module = types.ModuleType(name)
    module.__file__ = str(ROOT / relative)
    exec(compile(data, module.__file__, 'exec'), module.__dict__)
    return module


def exact(value):
    need(type(value) in (str, int, F), 'exact rational schema')
    text = str(value)
    need(len(text) <= 128, 'rational length')
    return F(text)


def domain(rows, key, expected):
    need(type(rows) is list and len(rows) == len(expected), 'complete finite row count')
    actual = [key(r) for r in rows]
    need(len(actual) == len(set(actual)) and set(actual) == set(expected), 'exact finite model key domain')
    return rows


def gram(n, edges=()):
    need(type(n) is int and 1 <= n <= 64, 'model order')
    g = [['2' if i == j else '0' for j in range(n)] for i in range(n)]
    for i, j, magnitude in edges:
        need(0 <= i < n and 0 <= j < n and i != j and g[i][j] == '0', 'distinct model edge')
        g[i][j] = g[j][i] = '-(' + str(magnitude) + ')'
    return g


def enlarged(base, n, edges):
    need(len(base) <= n, 'extension order')
    out = gram(n)
    for i, row in enumerate(base):
        for j, value in enumerate(row):
            out[i][j] = str(value)
    for i, j, magnitude in edges:
        need(0 <= i < n and 0 <= j < n and i != j and out[i][j] == '0', 'new physical edge')
        out[i][j] = out[j][i] = '-(' + str(magnitude) + ')'
    return out


def path_edges(nodes, magnitude='1'):
    return [(i, j, magnitude) for i, j in zip(nodes, nodes[1:])]


def affine(kind, order):
    target = E.graph(kind, order, W, affine=True)
    magnitudes = {3: '1', 4: 'sqrt(2)', 5: '(1+sqrt(5))/2', 6: 'sqrt(3)', 99: '2'}
    return gram(order, [(i, j, magnitudes[label]) for i, row in enumerate(target)
                       for j, label in row.items() if i < j])


def one(key, g, white, marks, minimum, premise='whole marked factors from the named source proof'):
    E.charge(FUEL)
    need(len(ROWS) < MAX_ROWS, 'output row cap')
    result = E.corner(g, list(white), list(marks), W, FUEL)
    need(exact(result['weight']) >= exact(minimum), 'current literal below chosen output: ' + key)
    ROWS.append(dict(key=key, chosen_lower=str(exact(minimum)), current=result,
                     whole_factor_premise=premise))
    return exact(result['weight'])


def receipt(key, g, white, rows, minimum=1):
    seen = set(); total = F()
    for i, (marks, lower) in enumerate(rows):
        pair = tuple(sorted(marks))
        need(len(pair) == 2 and pair not in seen, 'distinct actual receipt corners')
        seen.add(pair)
        total += one(key + '/corner/' + str(i), g, white, marks, lower)
    need(total >= exact(minimum), 'current distinct-corner total below required receipt: ' + key)
    return dict(key=key, current_total=str(total), required_total=str(exact(minimum)),
                distinct_corners=len(seen))


def tail(base, root, rank, attachment):
    n = len(base)
    return enlarged(base, n + 1 + rank,
        [(n, attachment, '1')] + path_edges([n] + list(range(n + 1, n + 1 + rank))))


def immediate_profiles(key, base_kind, order, selected, rows):
    """Rebuild the exact zero-own immediate opposite from its native profile.

    Zero current intrinsic values are named priority LINKS, never units or a
    successful numerical receipt. Their independent whole-face theorem remains
    a separately listed geometric premise.
    """
    base = affine(base_kind, order)
    graph = E.graph(base_kind, order, W, affine=True)
    selected_value = W.affine_weight(base_kind, order, *selected)
    expected = {(mark, name) for mark in selected for name in ('end4', 'next3')}
    need(len(rows) == 4 and {(r['mark'], r['profile']) for r in rows} == expected,
         'complete native immediate profile domain')
    for row in rows:
        mark = row['mark']; other = next(i for i in selected if i != mark)
        neighbor = next(iter(graph[mark]))
        if len(graph[neighbor]) == 3:
            first = next(i for i in graph[neighbor] if i != mark and len(graph[i]) == 1)
            second = neighbor
        else:
            first = neighbor
            second = next(i for i in graph[neighbor] if i != mark)
        contact = first if row['profile'] == 'end4' else second
        ids = [i for i in range(order) if i != mark]
        old = [[base[i][j] for j in ids] for i in ids]
        g = enlarged(old, order, [(order - 1, ids.index(contact),
                                  'sqrt(2)' if row['profile'] == 'end4' else '1')])
        marks = [order - 1, ids.index(other)]
        value = one(key + '/' + str(mark) + '/' + row['profile'], g,
                    [i for i in range(order) if i not in marks], marks, 0)
        current = ROWS[-1]['current']
        need(current['kind'] == 'ideal' and current['forest'] == [[row['target'], order]],
             'actual full immediate factor versus stored target type')
        if value == 0 or selected_value == 0:
            priority = ('fresh_order10_c_endpoint_zero_transfer_20261018' if order == 10
                        else 'fresh_order11_bc_zero_transfer_20261018')
            for folder, ext in (('scripts', 'py'), ('results', 'md'), ('results', 'json')):
                captured(f'{folder}/{priority}.{ext}')
            PRIORITIES.append(dict(key=key, profile=row['profile'], mark=mark,
                intrinsic_selected=str(selected_value), intrinsic_distinct=str(value),
                required_independent_source=priority,
                link='Apply independent current B/C zero whole-face priority first; no intrinsic unit is asserted.'))
        else:
            need(selected_value + value >= 1, 'immediate positive intrinsic receipt below one')


def order10(batch, doc):
    receipts = []
    if batch == 'order10-cycle':
        base = affine('A', 10)
        cycle_profiles(batch, doc)
        rows = domain(doc['whole_coupon_corners'], lambda r: r['spectator_rank'], range(4))
        for r in rows:
            rank = r['spectator_rank']; g = tail(base, 10, rank, 0)
            one('cycle10/coupon/' + str(rank), g, [i for i in range(len(g)) if i not in (0, 5, 10)],
                [10, 5], r['half_weight'])
        mixed = domain(doc['mixed_whole_corners'], lambda r: (r['spectator'], r['orientation']),
                       [('A7', 0), ('D6', 0), ('D6', 1)])
        for i, r in enumerate(mixed):
            g = r['gram']; white = [j for j in range(len(g)) if j not in (0, 5, 10, 11)]
            one('cycle10/mixed/' + str(i), g, white, [10, 11], r['payment'])
    elif batch == 'order10-b':
        base = affine('B', 10)
        immediate_profiles('B10/immediate', 'B', 10, (8, 7), doc['physical']['immediate'])
        rows = domain(doc['finite_corners'], lambda r: (r['extreme_attachment'], r['spectator_path_rank']),
                      list(product((8, 9), range(4))))
        for i, r in enumerate(rows):
            g = tail(base, 10, r['spectator_path_rank'], r['extreme_attachment'])
            one('B10/coupon/' + str(i), g, [j for j in range(len(g)) if j not in (8, 7, 10)],
                [10, 7], r['payment'])
    else:
        base = affine('D', 10)
        immediate_profiles('D10/immediate', 'D', 10, (6, 8), doc['physical']['immediate'])
        rows = domain(doc['whole_finite_coupons'], lambda r: (r['old_attachment'], r['spectator_rank']),
                      list(product((6, 7), range(8))))
        for i, r in enumerate(rows):
            g = tail(base, 10, r['spectator_rank'], r['old_attachment'])
            one('D10/coupon/' + str(i), g, [j for j in range(len(g)) if j not in (6, 8, 10)],
                r['marks'], r['payment'])
        keys = [(kind, ori, a, b) for kind in ('A15', 'D7')
                for ori in range(1 if kind == 'A15' else 2) for a, b in product((6, 7), (8, 9))]
        mixed = domain(doc['mixed_whole_corners'],
            lambda r: (r['spectator'], r['orientation'], *r['old_attachments']), keys)
        for i, r in enumerate(mixed):
            if not r['viable']:
                need(r['payment'] is None, 'excluded model cannot earn a payment')
                continue
            payment = r['payment']; c = r['corners'][payment['corner']]
            one('D10/mixed/' + str(i), r['gram'],
                [j for j in range(len(r['gram'])) if j not in (6, 8, 10, 11)],
                c['marks'], payment['weight'])
        need(sum(bool(r['viable']) for r in mixed) == 6, 'complete viable mixed model count')
    return receipts


def cycle_profiles(batch, doc):
    order = {'order10-cycle': 10, 'cycle11': 11, 'cycle12': 12, 'cycle13': 13, 'cycle14': 14, 'cycle16': 16}[batch]
    rows = doc['profiles']['single_profiles'] if order in (10, 11) else doc['single_profile_checks'] if order == 16 else doc['single_profiles']
    arc_pairs = [(5, 5)] if order == 10 else [(5, 6)] if order == 11 else [(6, 6)] if order == 12 else [(6, 7)] if order == 13 else [(8, 8)] if order == 16 else [(6, 8), (7, 7)]
    need(len(rows) == 6 * len(arc_pairs), 'complete cycle single-profile count')
    expected = {(a, b, name, side) for a, b in arc_pairs for name in ('end3', 'end4', 'next3') for side in (0, 1)}
    seen = set()
    for i, r in enumerate(rows):
        name = r['profile']; side = int(name.rsplit('_', 1)[1]); family = name.rsplit('_', 1)[0]
        need(family in ('end3', 'end4', 'next3') and side in (0, 1), 'cycle profile name')
        a, b = (r['P'], r['Q']) if order == 14 else arc_pairs[0]
        key = (a, b, family, side); need(key in expected and key not in seen, 'cycle profile key')
        seen.add(key)
        position = 2 if family == 'next3' else 1
        contact = position if side == 0 else order - position
        g = gram(order, path_edges(list(range(order - 1))) +
                 [(order - 1, contact - 1, 'sqrt(2)' if family == 'end4' else '1')])
        marks = [order - 1, a - 1]
        one(batch + '/single/' + str(i), g, [j for j in range(order) if j not in marks], marks, 1)
    need(seen == expected, 'complete cycle profile set')
    if order in (10, 11):
        if order == 11:
            g = doc['profiles']['both_ends_Gram']
        else:
            g = gram(10, path_edges(list(range(9))) + [(9, 0, '1'), (9, 8, '1')])
        marks = [order - 1, 4]
        one(batch + '/both-ends', g, [j for j in range(order) if j not in marks], marks, '1/2')


def cycles(batch, doc):
    receipts = []
    if batch == 'cycle18-tree':
        for name, order in (('cycle18_equality', 18), ('cycle16_end3_equality', 16)):
            r = doc[name]
            need(r['order'] == order and len(r['corners']) == 4, 'exact large-cycle equality context')
            receipts.append(receipt(name, r['gram'], r['white'],
                    [(c['marks'], exact(c['unit']) + exact(c['coupon'])) for c in r['corners']]))
        tree = doc['tree']
        profiles = ('isolated', 'outer_end3', 'outer_end4', 'outer_next3',
                    'inner_end3', 'inner_next3', 'both_ends3', 'inner_end4')
        rows = domain(tree['rows'], lambda r: (r['P'], r['c'], r['d']),
                      list(product((6, 7, 8, 9), profiles, profiles)))
        counts = {name: sum(r['kind'] == name for r in rows)
                  for name in ('immediate_payment', 'positive_gap', 'sealed_finite_A7')}
        need(counts == dict(immediate_payment=45, positive_gap=210, sealed_finite_A7=1), 'complete tree-profile partition')
        one('tree/sealed-A7', gram(7, path_edges(list(range(7)))), list(range(1, 6)), [0, 6], 1)
        outer = tree['added_outer_exception']
        g = gram(9, [(i, j, '1') for i, j in outer['edges']])
        typ, images, ideal = E.classify(E.matrix(g), list(range(9)), [0, 8], W, FUEL)
        need(typ == ('E', 9) and ideal, 'outer D8 prolongation is exact affineE8')
        ROWS.append(dict(key='tree/outer-D8-exclusion', graph_type=list(typ),
                         scope='Native affine type binding ONLY: source proves this block lies in a PD actual edge; no payment credited.'))
        return receipts
    cycle_profiles(batch, doc)
    if batch == 'cycle12':
        names = ['isolated_1_sqrt2', 'isolated_sqrt2_1', 'mixed_isolated_both', 'mixed_both_isolated',
                 'two_A1', 'one_sqrt2_A1', 'one_two_A1', 'one_A3_middle', 'one_B4', 'one_D4',
                 'twins_A3', 'twins_D4', 'twins_D5', 'twins_D6', 'twins_D7']
        rows = domain(doc['literal_cases'], lambda r: r['name'], names)
        for r in rows:
            if r['name'] in ('twins_D6', 'twins_D7'):
                need(r['outcome'] == 'excluded by a negative direction at the proposed actual corner', 'excluded cycle12 schema')
                continue
            if r['outcome'] == 'retained short path':
                marks = [12, 13]
                one('cycle12/' + r['name'], r['gram'], r['white'], marks, 1)
            else:
                selected = [(c['black_pair'], sum(exact(q['weight']) for q in c['coupons'])) for c in r['corners']]
                receipts.append(receipt('cycle12/' + r['name'], r['gram'], r['white'], selected))
    elif batch == 'cycle14':
        g = gram(17, [(i, (i + 1) % 14, '1') for i in range(14)] +
                 [(0, 14, '1'), (7, 15, '1'), (14, 16, '1')])
        q = doc['equality']
        need(q['white'] == [i for i in range(14) if i not in (0, 7)] + [16], 'cycle14 physical white support')
        need([r['marks'] for r in q['corners']] == [[0, 7], [7, 14], [14, 15], [15, 0]], 'cycle14 distinct corner domain')
        # Old halves were deliberately redistributed to41/64+23/64.
        receipts.append(receipt('cycle14/equality-current', g, q['white'], [(r['marks'], 0) for r in q['corners']]))
    elif batch == 'cycle16':
        r = doc['equality_core']
        need(len(r['gram']) == 18 and len(r['corners']) == 4, 'exact cycle16 equality domain')
        receipts.append(receipt('cycle16/equality', r['gram'], r['white'], [(c['marks'], c['coupon']) for c in r['corners']]))
    return receipts


def order11(batch, doc):
    if batch == 'order11-bc':
        rows = domain(doc['B_equality_corners'], lambda r: (r['heavy_profile'], r['fork_attachment']),
                      list(product(('isolated', 'next3'), (9, 10))))
        return [receipt('BC11/equality/' + str(i), r['gram'], r['white'],
                [(r['first_corner'], r['first_payment']), (r['second_corner'], r['second_payment'])])
                for i, r in enumerate(rows)]
    immediate = doc['formulas']['whole_affine_half_alternatives']
    need(len(immediate) == 12, 'complete BCD11 immediate set')
    for kind, selected in (('B', (8, 9)), ('C', (0, 10)), ('D', (7, 9))):
        immediate_profiles('BCD11/' + kind + '/immediate', kind, 11, selected,
                           [r for r in immediate if r['kind'] == kind])
    small = domain(doc['selected_whole_small_corners'],
        lambda r: (r['kind'], r['attached_spectator_A_rank'], r['old_attachment']),
        [('B', rank, attach) for rank, attach in product(range(2), (9, 10))] +
        [('D', rank, attach) for rank, attach in product(range(4), (7, 8))])
    for i, r in enumerate(small):
        kind = r['kind']; old, far = (9, 8) if kind == 'B' else (7, 9)
        g = tail(affine(kind, 11), 11, r['attached_spectator_A_rank'], r['old_attachment'])
        need(r['selected_corner_marks'] == [11, far], 'BCD11 marked-output source')
        one('BCD11/small/' + str(i), g, [j for j in range(len(g)) if j not in (old, far, 11)],
            r['selected_corner_marks'], r['weight'])
    rows = domain(doc['sealed_far_spectator_corners'],
        lambda r: (r['spectator'], r['orientation'], *r['old_attachments']),
        [(kind, ori, a, b) for kind in ('A7', 'D6') for ori in range(1 if kind == 'A7' else 2)
         for a, b in product((7, 8), (9, 10))])
    for i, r in enumerate(rows):
        if r['outcome'] == 'excluded by actual edge/corner':
            continue
        need(r['outcome'] == 'retained finite unit', 'BCD11 mixed outcome')
        selected = [c for c in r['corners'] if c['retained_finite_unit']]
        need(bool(selected), 'nonempty current-output predecessor')
        white = [j for j in range(len(r['gram'])) if j not in (7, 9, 11, 12)]
        for j, c in enumerate(selected):
            one('BCD11/mixed/' + str(i) + '/' + str(j), r['gram'], white, c['marks'], 1)
    return []


def order12(batch, doc):
    receipts = []
    if batch == 'order12-center':
        rows = domain(doc['literal_cases'], lambda r: (r['center_leaf'], r['heavy_profile']),
                      list(product((10, 11), ('isolated', 'end3', 'next3'))))
        for i, r in enumerate(rows):
            receipts.append(receipt('B12/center/' + str(i), r['gram'], r['supporting_white'],
                [(c['marks'], c['payment']) for c in r['corners']]))
    elif batch == 'order12-near':
        rows = domain(doc['whole_finite_corners'], lambda r: r['profile'], ('isolated', 'mate3'))
        for r in rows:
            one('BC12/finite/' + r['profile'], r['gram'], r['white'],
                r['paying_corner'], r['distinct_payment'])
        immediate = domain(doc['immediate_affine_roles'], lambda r: r['profile'], ('end4', 'next3'))
        for r in immediate:
            nodes = r['nodes']; marks = [nodes.index(i) for i in (10, 12)]
            one('BC12/immediate/' + r['profile'], r['gram'],
                [j for j in range(len(nodes)) if j not in marks], marks, r['current_weight'])
            if exact(r['current_weight']) == 0:
                priority = 'fresh_order11_bc_zero_transfer_20261018'
                for folder, ext in (('scripts', 'py'), ('results', 'md'), ('results', 'json')):
                    captured(f'{folder}/{priority}.{ext}')
                need(r['payment'] == 'independent_BC11_zero_theorem_first', 'exact BC12 zero dispatch')
                PRIORITIES.append(dict(key='BC12/immediate/' + r['profile'],
                    intrinsic_distinct='0', required_independent_source=priority,
                    status='NAMED_PRIORITY_SOURCE_NOT_REPLAYED',
                    link='Independent BC11 zero whole-face theorem first; no intrinsic unit is asserted.'))
    else:
        expected = [98, 100, 119, 121, 164, 167, 171, 1452, 1728, 1734, 1745, 2121, 2127, 2138]
        rows = domain(doc['B']['cases'], lambda r: r['index'], expected)
        count = 0
        for r in rows:
            core = r['inherited_geometry']; g = core['gram']; white = core['white']
            if not r['payment_variants']:
                need(core['outcome'] == 'whole_finite_B11_before_any_affine_completion', 'proper-only B12 output')
                one('BC12/selected/' + str(r['index']), g, white, core['preferred_corner'], 1)
            for variant in r['payment_variants']:
                count += 1
                extended = variant['variant'] == 'one_extra_root_payment'
                gg = core['one_extra_root_gram'] if extended else g
                ww = white + [len(g)] if extended else white
                receipts.append(receipt('BC12/selected/' + str(r['index']) + '/' + variant['variant'], gg, ww,
                    [(c['marks'], c['new_weight']) for c in variant['corners']]))
        need(count == 16, 'complete selected B12 physical variant count')
        crows = domain(doc['C']['cases'], lambda r: r['index'], (164, 1452))
        b_by_index = {r['index']: r for r in rows}
        need(sum(len(r['maps']) for r in crows) == 4, 'complete selected C12 endpoint map count')
        for r in crows:
            need(r['end_projection_certificate'] == b_by_index[r['index']]['inherited_geometry']
                 and not b_by_index[r['index']]['payment_variants']
                 and exact(r['retained_finite_unit']) == 1,
                 'C12 proper-completion exclusion uses the SAME bound B11 whole-corner receipt')
    return receipts


def order13_sparse(doc):
    records = domain(doc['records'], lambda r: r['kind'], ('B', 'C'))
    receipts = []; seals = set()
    for record in records:
        expected = [157, 174, 952, 2474] if record['kind'] == 'B' else [625, 1556, 1642]
        cases = domain(record['cases'], lambda r: r['index'], expected)
        for r in cases:
            proof = r['proof']
            if 'inherited_full_core_seal' not in proof:
                need('maps' in proof and bool(proof['maps']), 'nonempty proper-completion exclusion input')
                continue
            need(record['kind'] == 'B' and r['index'] in (174, 2474), 'sparse B13 radical source domain')
            seals.add(r['index']); core = proof['inherited_full_core_seal']
            receipts.append(receipt('BC13/sparse/' + str(r['index']), core['gram'], core['white'],
                [(c['marks'], c['payment']) for c in proof['corners']]))
    need(seals == {174, 2474}, 'complete sparse radical seal set')
    return receipts


def typed(batch, doc):
    need(batch == 'large-d', 'typed family')
    # All types of finite rank<=5; I2(m),m>=5 is ONE symbolic weight class.
    types_small = [('A', n) for n in range(1, 6)] + [('B', n) for n in range(2, 6)] + \
                  [('D', 4), ('D', 5), ('F', 4), ('H', 3), ('H', 4), ('I', 2)]
    for typ in types_small:
        for mark in range(typ[1]):
            for leaf in (8, 9):
                E.charge(FUEL)
                value = W.finite_cross_weight(typ, mark, ('D', 10), leaf)
                need(value >= F(1, 2), 'universal small/D10 short-leaf ordinary half')
                ROWS.append(dict(key='largeD/small/' + str((typ, mark, leaf)),
                    typed_input=list(typ), mark=mark, D10_leaf=leaf, current=str(value),
                    domain='Every finite type of rank<=5; I denotes all I2(m),m>=5, same printed weight.'))
    for rank in range(20, 59):
        for leaf in (rank - 2, rank - 1):
            value = W.finite_same_weight(('D', rank), rank - 11, leaf)
            need(value == F(1, 2), 'largeD intrinsic short-leaf9 half')
            ROWS.append(dict(key='largeD/rank/' + str((rank, leaf)), type=['D', rank],
                marks=[rank - 11, leaf], current=str(value), actual_ambient58_rank_bound=True))
    need(len(doc['coupon_cases']) == 16, 'complete sixteen source coupon geometries')
    for r in doc['coupon_cases']:
        need(r['small_component_order_bound'] in (1, 5), 'source whole small-order bound')
        need(len(r['exposed_D10']) == 2, 'both actual D marked orientations')
        for i, q in enumerate(r['exposed_D10']):
            E.charge(FUEL)
            g = E.matrix(r['core_gram']); nodes = q['whole_D10_nodes']
            typ, images, is_affine = E.classify(g, nodes, [q['D10_mark']], W, FUEL)
            need(typ == ('D', 10) and not is_affine and set(images) == {(8,), (9,)}, 'physical exposed D10 short-leaf role')
            ROWS.append(dict(key='largeD/exposed/' + str((r['index'], i)), type=list(typ),
                all_mark_images=[list(p) for p in images], small_rank_bound=r['small_component_order_bound'],
                whole_seal_premise='Source all-rank tail/indirect-contact proof, not the displayed D20 representative alone.'))
    cycle = domain(doc['cycle_D_leaf9_orders'], lambda r: r['shorter_arc'], (6, 7, 8, 9))
    for r in cycle:
        need(r['longer_arc'] == 10 and r['finite_D_order'] == r['shorter_arc'] + 10 <= 19,
             'symbolic cycle longer-arc forced D-order bound')
    return []


def native_corridors(batch):
    """Complete ambient58 marked-role overlist, not an all-rank geometry proof.

    An extension of a displayed finite path can change its whole type/rank.
    The fixed no-spectator models do NOT cover these possibilities. This
    overlist checks all actual classical ranks and current zeros explicitly;
    a zero is dispatched to the separate G6 deep-path theorem, never counted
    as an intrinsic unit. Leaf9 short-fork halves retain their other half.
    """
    shard = int(batch.rsplit('-', 1)[1]); need(0 <= shard < 4, 'corridor shard')
    providers = {6: ('finite57.current.deep_finite6.v1', 'fresh_depth6_finite6_zero_transfer_20261018'),
                 7: ('finite57.current.deep_finite7.v1', 'fresh_depth6_finite7_adaptive_transfer_20261018'),
                 8: ('finite57.current.deep_finite8.v1', 'fresh_depth6_finite8_zero_transfer_20261018')}
    for provider, stem in providers.values():
        for folder, ext in (('scripts', 'py'), ('results', 'md'), ('results', 'json')):
            captured(f'{folder}/{stem}.{ext}')
    for stem in ('fresh_cycle11_half_pair_transfer_20261018',
                 'fresh_cycle12_ideal5_transfer_20261017',
                 'fresh_cycle13_ideal5_transfer_20261017',
                 'fresh_cycle14_ideal5_transfer_20261017'):
        captured('scripts/' + stem + '.py'); captured('results/' + stem + '.md')
    for rank in range(9, 59):
        if rank % 4 != shard:
            continue
        for kind in ('A', 'B', 'D'):
            graph = W.classical_graph(kind, rank)
            distances = {}
            for start in range(rank):
                ds = {start: 0}; queue = [start]
                for node in queue:
                    for other in graph[node]:
                        if other not in ds:
                            ds[other] = ds[node] + 1; queue.append(other)
                need(len(ds) == rank, 'connected typed target')
                distances[start] = ds
            features = [i for i in range(rank) if len(graph[i]) != 2
                        or any(label != 3 for label in graph[i].values())]
            units = zeros = halves = 0
            for a in range(rank):
                for b in range(a + 1, rank):
                    distance = distances[a][b]
                    leaf9 = distance == 9 and (len(graph[a]) == 1 or len(graph[b]) == 1)
                    if distance > 8 and not leaf9:
                        continue
                    E.charge(FUEL)
                    value = W.finite_same_weight((kind, rank), a, b)
                    if value == 0:
                        depths = [min(distances[mark][f] for f in features) for mark in (a, b)]
                        need(distance in (6, 7, 8) and
                             min(depths) >= 6,
                             'zero corridor needs its EXACT deep-feature priority')
                        zeros += 1
                        need(len(PRIORITIES) < MAX_ROWS, 'bounded outgoing G6 provider-key count')
                        provider, source = providers[distance]
                        PRIORITIES.append(dict(provider_id=provider, ambient_dimension=58,
                            factor_kind=kind, factor_rank=rank, canonical_marks=[a, b],
                            marked_distance=distance, feature_depths=depths,
                            current_intrinsic_weight='0', source_family=source,
                            status='UNBOUND_G6_CURRENT_OUTPUT_BRIDGE'))
                    elif leaf9 and kind == 'D' and rank >= 20 and {a, b} & {rank - 2, rank - 1}:
                        need(value == F(1, 2), 'largeD leaf9 current half')
                        halves += 1
                    else:
                        need(value == 1, 'retained current short corridor unit')
                        units += 1
            ROWS.append(dict(key='native-corridor/' + str((kind, rank)), type=[kind, rank],
                unit_marked_roles=units, zero_deep_marked_roles=zeros, largeD_leaf9_halves=halves,
                full_marked_role_overlist='All pairs at distance<=8 or distance9 with a marked leaf',
                mathematical_premise='Finite classification at rank>=9 is A/B/D; actual ambient58 bounds rank.'))
            if halves:
                PRIORITIES.append(dict(key='native-corridor/' + str((kind, rank)), count=halves,
                    link='Large finiteD short-leaf9 theorem needs its distinct whole-small/D10 short-leaf half; this batch checks no geometric seal.'))
    return []


def main():
    global PINS, W, E
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=('index', *STEMS), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    lock_path = ROOT / 'formal/.proof-batches/exclusive.lock'
    need(lock_path.is_file(), 'shared proof lock is missing')
    with lock_path.open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
            raise ValueError('shared parent proof lock must already be held')
    manifest_bytes = raw(MANIFEST)
    need(sha256(manifest_bytes).hexdigest() == MANIFEST_SHA, 'frozen G4 source manifest identity')
    manifest = decode(manifest_bytes)
    entries = manifest['inputs']
    need(type(entries) is list and len(entries) <= 100, 'bounded source manifest')
    for entry in entries:
        need(entry['path'] not in PINS, 'duplicate input source path')
        PINS[entry['path']] = entry['sha256']
        SIZES[entry['path']] = entry['bytes']
    need(PINS[ENGINE] == ENGINE_SHA and PINS[LIB] == LIB_SHA, 'frozen engine/literal versions')
    source_name = str(Path(__file__).resolve().relative_to(ROOT)); own_source = raw(source_name)
    captured('manuscripts/finite_volume57.tex')
    E = load('finite57_g4_exact_engine', ENGINE, captured(ENGINE))
    W = load('finite57_g4_literal_weights', LIB, captured(LIB))
    receipts = []
    if args.batch == 'index':
        for relative in PINS:
            captured(relative)
        status = 'SOURCE_INDEX_ONLY_NO_TERMINAL_CHECK'
    else:
        doc = packet(args.batch)
        if args.batch.startswith('order10-'):
            receipts = order10(args.batch, doc)
        elif args.batch.startswith('cycle'):
            receipts = cycles(args.batch, doc)
        elif args.batch.startswith('order11-'):
            receipts = order11(args.batch, doc)
        elif args.batch.startswith('order12-'):
            receipts = order12(args.batch, doc)
        elif args.batch == 'order13-sparse':
            receipts = order13_sparse(doc)
        elif args.batch.startswith('native-corridors-'):
            receipts = native_corridors(args.batch)
        else:
            receipts = typed(args.batch, doc)
        need(bool(ROWS), 'empty terminal binding cannot earn PASS')
        status = 'PASS_G4_CURRENT_TERMINALS_DECLARED_SCOPE'
    for relative, pin in READ.items():
        need(sha256(raw(relative)).hexdigest() == pin, 'input changed during batch')
    need(raw(MANIFEST) == manifest_bytes and raw(source_name) == own_source, 'controller sources changed during batch')
    print(json.dumps(dict(status=status, batch=args.batch, source_sha256=sha256(own_source).hexdigest(),
        manifest_sha256=MANIFEST_SHA, read_inputs=READ, rows=ROWS, physical_receipts=receipts,
        fuel_used=FUEL[0], unmarked_witnesses_credited=False,
        ambient_premise='Actual dimension58: finite whole rank<=58; cusp whole affine order<=58. No earlier59–65 stage is rebound.',
        remaining_premises=['Complete original geometric/native parameter reductions and selected-model domains',
            'Actual whole marked-factor seals and physical-corner incidence',
            'Inherited inertia/proper-completion exclusions, not recomputed here',
            'Zero intrinsic B/C outputs use the named independent earlier-priority/short-path theorems; no intrinsic unit is inferred'],
        earlier_priority_links=PRIORITIES,
        scope='Selected current terminal values and explicitly recomputed distinct-corner receipts only. No universal all-rank geometry, complete G4/face routing or finite57 PASS.'), sort_keys=True))


if __name__ == '__main__':
    main()
