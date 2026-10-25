#!/usr/bin/env python3
"""One original quiet case, eight ordered witness contexts, at explicit premises.

The unchanged parent supplied192 applications. This worker selects ONE of
its24 physical cases before matrix/corner work. It never runs the parent
final_quiet function, and no child asserts completion of the24-case union.
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
BASE = 'scripts/check_finite57_current_universal_receivers_20261025.py'
BASE_SHA = '29ddbcf65cae2a7c7ff7703a1b4a81bad3ca1fdb0ade4bdbffda26fb9c8a1631'
CASE_INDICES = tuple(3+6*i+j for i in range(8) for j in range(3))
BATCHES = tuple(f'quiet-case-{i:02d}' for i in CASE_INDICES)
MAX_OUTPUT = 500*1024


def need(condition, message):
    if not condition:
        raise ValueError(message)


def raw(name):
    path = (ROOT/name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'existing repository source')
    with path.open('rb') as inp:
        data = inp.read(4*1024**2+1)
    need(0 < len(data) <= 4*1024**2, 'bounded source read')
    return data


def quiet_case(u, c, index):
    need(type(index) is int and index in CASE_INDICES, 'one exact original quiet case')
    original = c.field(u.U, ('nonextreme_zero_cases',))
    repairs = c.field(c.FINAL, ('quiet_repairs',))
    need(len(original) == 48, 'complete original order13 nonextreme48 domain')
    by = u.domain(c, repairs, lambda p: p['index'], list(CASE_INDICES),
                  'exact finalguard24 original quiet indices')
    # Select BEFORE quiet_physical, E.matrix, current, role or classify.
    p = by[index]
    q = original[index]
    old = u.ORDER13_MARKS[index//12]
    side = ('left', 'right')[(index%12)//6]
    profile = ('isolated4', 'mate4', 'center3')[index%3]
    need(tuple(q['marks']) == old and q['nonextreme_side'] == side and
         q['extreme_profile'] == 'mate3' and
         q['nonextreme_profile'] == p['profile'] == profile,
         'exact original ordered quiet-profile source association')
    g = u.quiet_physical(c, q)
    white = sorted(set(range(13))-set(old))
    need(c.E.matrix(g) == c.E.matrix(q['gram']) and q['white'] == white,
         'complete actual quiet Gram and original white identity')
    quiet = [14, old[0]] if side == 'left' else [old[1], 13]
    paying = (([old[1], 13] if side == 'left' else [14, old[0]])
              if profile != 'center3' else [13, 14])
    need(p['quiet_corner'] == quiet and p['paying_corner'] == paying and
         set(quiet) != set(paying), 'two DISTINCT actual quiet/pay corners from SAME physical star')
    out = []
    for kinds in product(('A', 'D'), repeat=3):
        extra = tuple((kind, 15) for kind in kinds)
        first = u.current(c, g, white, quiet, '1/2', extra, True)
        need(first['current']['complete'] and not first['current']['ignored_unsealed_components'] and
             sum(t[1] for t in first['current']['forest']) == 58,
             'complete guardTRUE finite rank58 quiet corner')
        premises = [u.assumed('accepted_full_star_zero_energy_seal', source=u.BRIDGE,
            assertion='R12 nonextreme8+4 gives E0; fullZ inverse SPD and disjoint columns force BOTH entire columns zero. Same three labelled whole A15/D15 factors are orthogonal to the entire star and remain whole at BOTH receivers.'),
            u.assumed('complete_original_native_profiles_and_actual_factor_fiber_geometry', source=u.U)]
        if profile == 'mate4':
            second = u.current(c, g, white, paying, '0')
            need(second['current']['kind'] == 'ideal' and second['current']['forest'] == [['B', 13]] and
                 second['current']['marked_images'] == p['images'] == [[10, 11], [10, 12]],
                 'exact WHOLE current B13zero priority image')
            premises.append(u.assumed('independent_quiet_B13_current_priority',
                theorem='finite57.current.quiet_B13_long_orbit.whole_face.v1',
                worker='scripts/check_finite57_quiet_b13_current_priority_bridge_20261025.py',
                batch='quiet-b13', consumer_key=[index, *kinds], marked_images=p['images']))
            total = None
        else:
            second = u.current(c, g, white, paying, '1/2' if profile == 'center3' else '1', extra, True)
            need(second['current']['complete'] and not second['current']['ignored_unsealed_components'] and
                 sum(t[1] for t in second['current']['forest']) == 58,
                 'same complete actual rank58 paying forest')
            total = str(F(first['current']['weight'])+F(second['current']['weight']))
        out.append(u.application([index, *kinds], 'finite57.current.final_guard_quiet_receiver.whole_face.v1',
            [c.FINAL, 'quiet_repairs', repairs.index(p)], premises,
            original_source_selector=[u.U, 'nonextreme_zero_cases', index],
            original_source_context_sha256=u.context(c, q['gram'], white),
            physical_context_sha256=u.context(c, g, white),
            ordered_whole_witnesses=[dict(copy=i, type=[k, 15]) for i, k in enumerate(kinds)],
            same_three_witness_copies_at_both_corners=True, quiet_receiver=first, paying_receiver=second,
            distinct_current_total=total, intrinsic_zero_credited_as_unit=False))
    expected = {(index, *kinds) for kinds in product(('A', 'D'), repeat=3)}
    keys = [tuple(p['application_key']) for p in out]
    need(len(keys) == len(set(keys)) == 8 and set(keys) == expected,
         'exact eight ordered witness contexts for THIS physical case')
    return dict(case_index=index, original24_case_domain=list(CASE_INDICES),
        theorem_applications=out, case_complete_at_declared_premises=True,
        quiet24_complete_at_declared_premises=False,
        independent_quiet_B13_runtime_verified=False, finalguard_one_sided_eight_not_recomputed=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=BATCHES, required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    own = str(Path(__file__).resolve().relative_to(ROOT))
    before, base = raw(own), raw(BASE)
    need(sha256(base).hexdigest() == BASE_SHA, 'exact unchanged universal parent')
    u = types.ModuleType('finite57_quiet_case_parent')
    u.__file__ = str(ROOT/BASE)
    exec(compile(base, u.__file__, 'exec'), u.__dict__)
    a, r, c, extra, ledger = u.prepare()  # Enforces the already-held shared lock.
    need(BASE not in c.PINS or c.PINS[BASE] == BASE_SHA, 'parent pin collision')
    c.PINS[BASE] = BASE_SHA
    c.READ.add(BASE)
    result = quiet_case(u, c, int(args.batch.rsplit('-', 1)[1]))
    for name in c.READ:
        need(c.digest(name) == c.PINS[name], 'consumed quiet-case source/input changed')
    need(raw(own) == before and raw(BASE) == base and raw(u.LEDGER) == ledger and
         sha256(raw(u.LOADER)).hexdigest() == u.LOADER_SHA,
         'case worker/parent/shared loader/ledger changed')
    need(raw(a.EXTRA) == extra and sha256(raw(a.ROUTE)).hexdigest() == a.ROUTE_SHA and
         sha256(raw(r.COLLECTOR)).hexdigest() == r.COLLECTOR_SHA and
         sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA,
         'shared frozen loader chain changed')
    payload = dict(status='PASS_CURRENT_FINAL_QUIET_CASE_AT_DECLARED_PREMISES',
        batch=args.batch, source_sha256=sha256(before).hexdigest(), base_universal_source_sha256=BASE_SHA,
        result=result, consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)}, fuel_used=c.FUEL[0],
        scope='Only one exact original physical quiet case and its eight ordered whole witness contexts. Complete native profiles, saved matrix/spectrum inputs and whole-star geometry remain explicit premises; remaining cases and quiet-B13 priority require separate actual receipts.',
        independent_S8_verified=False, current_protected74_verified=False,
        current_terminal_bridge_complete=False, global_finite57_verified=False)
    encoded = json.dumps(payload, sort_keys=True)
    need(len(encoded.encode())+1 <= MAX_OUTPUT, 'complete case output cap; never truncate')
    print(encoded)


if __name__ == '__main__':
    main()
