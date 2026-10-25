#!/usr/bin/env python3
"""SOURCE ONLY: bind the64 quiet mate4 keys to repaired affineB13 priority.

This is the exact previously unbound consumer, with current direct receipts
for all four old long-orbit promotion cases. Its source-proved universal
applicability and geometric seals remain explicit premises.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
LOADER = 'scripts/index_finite57_funded_selected_whole_providers_20261025.py'
LOADER_SHA = '0a245deec21b283f0d47af0323cf9557fd6ba4465f6e25a2b7c47897f2a46e86'
EXTRA = 'results/finite57_quiet_b13_priority_extra_inputs_20261025.json'
EXTRA_SHA = 'd2d142fba95fb0db7e91137ed044529f2b3bbc10a44660075b56e7a7669f6fdc'
OLD = 'fresh_order13_affine_b_singleton_transfer_20261017'
PROMOTE = 'fresh_ideal5_weighted_dimension71_20261017'
U = 'fresh_affine_d12_leaf_pair_transfer_20261017'
PAIRS = [[10, 11], [10, 12]]
MAX_OUTPUT = 500 * 1024


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    path = (ROOT / relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded quiet bridge source path')
    with path.open('rb') as stream:
        data = stream.read(2 * 1024**2 + 1)
    need(0 < len(data) <= 2 * 1024**2, 'bounded quiet bridge source read')
    return data


def prepare():
    encoded = raw(LOADER)
    need(sha256(encoded).hexdigest() == LOADER_SHA, 'whole-provider loader identity')
    helper = types.ModuleType('finite57_quiet_current_loader')
    helper.__file__ = str(ROOT / LOADER)
    exec(compile(encoded, str(ROOT / LOADER), 'exec'), helper.__dict__)
    r, c = helper.prepare()
    data = raw(EXTRA)
    need(sha256(data).hexdigest() == EXTRA_SHA, 'quiet extra manifest identity')
    rows = json.loads(data, object_pairs_hook=c.unique)['inputs']
    expected = {f'{folder}/{stem}.{suffix}' for stem in (OLD, PROMOTE)
                for folder, suffix in (('scripts', 'py'), ('results', 'json'))}
    need(len(rows) == 4 and {p['path'] for p in rows} == expected, 'exact four extra priority inputs')
    for p in rows:
        if p['path'] in c.PINS:
            need(c.PINS[p['path']] == p['sha256'], 'extra/base input collision identity')
        c.PINS[p['path']] = p['sha256']
        need(c.digest(p['path']) == p['sha256'], 'quiet priority extra source/input pin')
        c.READ.add(p['path'])
    return helper, r, c


def promotion_receipts(c):
    old = c.field(OLD, ('conditional_long_orbit',))
    need(old['pairs'] == PAIRS and old['removed'] is False, 'exact diagnostic old long orbit')
    cases = {p['index']: p for p in old['cases']}
    need(len(cases) == len(old['cases']) == 4 and set(cases) == {166, 170, 1744, 2137},
         'all four actual long-orbit promotion cases')
    promoted = c.field(PROMOTE, ('order13_interfaces',))
    rows = promoted['promoted_B_distance11_cases']
    need({p['case'] for p in rows} == set(cases) and len(rows) == 4,
         'exact original promotion provider keys')
    need(len(promoted['records']) == 4 and {p['kind'] for p in promoted['records']} == {'A', 'B', 'C', 'D'},
         'exact four order13 promotion family records')
    b = next(p for p in promoted['records'] if p['kind'] == 'B')
    need(b['order'] == 13 and all(pair in b['removed_pairs'] for pair in PAIRS),
         'promoted theorem actually removed both current long-orbit pairs')
    repairs = c.b13_repairs()['receipts']
    fixed = {p['index']: p for p in repairs}
    need(set(fixed) == {166, 1744}, 'exact current finite B13 endpoint repairs')
    result = []
    for index, p in cases.items():
        current = c.receipt(OLD, index, p['gram'], p['white'], p['corners'],
                            'accepted positive-kernel/full-core-radical long-orbit seal')
        need(current['direct_current_receiver_condition'] and not current['legacy_pointwise_bounds_used'],
             'promotion uses only direct CURRENT receiver condition')
        if index in fixed:
            need(current['context_sha256'] == fixed[index]['context_sha256'],
                 'repaired endpoint consumer is exactly this old physical model')
            need({tuple(q['marks']) for q in current['corners']} ==
                 {tuple(q['marks']) for q in fixed[index]['corners']},
                 'current repair/old source actual corner-domain association')
        result.append(dict(case=index, current_receipt=current,
                           current_endpoint_repair=fixed.get(index),
                           saved_old_promotion_payment_used=False))
    return result


def quiet_consumers(c):
    need(all(c.W.affine_weight('B', 13, *pair) == 0 for pair in PAIRS),
         'current affineB13 long-orbit intrinsic values are zero')
    source = c.field(U, ('nonextreme_zero_cases',))
    all_repairs = c.field(c.FINAL, ('quiet_repairs',))
    need(len(all_repairs) == 24 and len({p['index'] for p in all_repairs}) == 24,
         'complete exact24 quiet source keys')
    chosen = [p for p in all_repairs if p['method'] == 'independent_current_affine_B13_priority']
    need({p['index'] for p in chosen} == set(range(4, 47, 6)) and len(chosen) == 8,
         'exact eight original quiet mate4 consumer keys')
    result = []
    for p in chosen:
        original = source[p['index']]
        need(original['extreme_profile'] == 'mate3' and original['nonextreme_profile'] == p['profile'] == 'mate4',
             'quiet original/profile association')
        value = c.value(original['gram'], original['white'], p['paying_corner'])
        need(value['kind'] == 'ideal' and value['forest'] == [['B', 13]], 'whole actual affineB13 priority factor')
        need(value['marked_images'] == p['images'] == PAIRS and F(value['weight']) == 0,
             'complete canonical long-orbit images; zero is not an intrinsic unit')
        for types in product(('A', 'D'), repeat=3):
            result.append(dict(consumer_key=[p['index'], *types],
                provider_key=dict(theorem='current affineB13 independent whole-face replacement',
                    kind='B', order=13, marked_images=PAIRS, quiet_source_index=p['index']),
                physical_context_sha256=sha256(c.canonical([original['gram'], original['white']]).encode()).hexdigest(),
                original_paying_marks=p['paying_corner'], intrinsic_current_weight='0',
                conclusion='Whole receiving face is paid by the repaired independent long-orbit theorem; no intrinsic unit is added.',
                applicability_premise='Accepted universal affineB13 full-face priority with four current long-orbit promotion receipts; its source DAG precedes the final D13 guard and permits arbitrary orthogonal spectator forests.'))
    need(len(result) == 64 and len({c.canonical(p['consumer_key']) for p in result}) == 64,
         'exact64 quiet rank15-choice consumer joins')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=('quiet-b13',), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative)
    helper, r, c = prepare()
    receipts = promotion_receipts(c)
    consumers = quiet_consumers(c)
    for path in c.READ:
        need(c.digest(path) == c.PINS[path], 'quiet priority input/source changed')
    need(raw(relative) == before and sha256(raw(LOADER)).hexdigest() == LOADER_SHA,
         'quiet bridge/loader source changed')
    need(sha256(raw(EXTRA)).hexdigest() == EXTRA_SHA, 'quiet extra manifest changed')
    need(sha256(helper.raw(helper.ROUTES)).hexdigest() == helper.ROUTES_SHA,
         'whole-provider route definitions changed')
    need(sha256(helper.raw(helper.COLLECTOR)).hexdigest() == helper.COLLECTOR_SHA,
         'direct collector source changed')
    need(sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA, 'base110 input manifest changed')
    payload = json.dumps(dict(status='PASS_NAMED_QUIET_B13_CURRENT_PRIORITY_BRIDGE', batch=args.batch,
        source_sha256=sha256(before).hexdigest(), current_terminal_bridge_complete=False,
        consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)}, fuel_used=c.FUEL[0],
        current_promotion_receipts=receipts, exact_quiet_consumer_joins=consumers,
        scope='Only64 quiet mate4 consumer keys and their four repaired direct-current promotion receipts. Universal source theorem/finite-data premises remain explicit; no complete G6/global finite57 PASS.'), sort_keys=True)
    need(len(payload.encode()) <= MAX_OUTPUT, 'quiet priority output cap; refusal not truncated keys')
    print(payload)


if __name__ == '__main__':
    main()
