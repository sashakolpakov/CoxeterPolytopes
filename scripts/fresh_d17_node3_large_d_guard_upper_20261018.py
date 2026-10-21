#!/usr/bin/env python3
"""Complete rank57 increasing-family upper for a proposed spectator guard.

This is an arithmetic diagnostic, not a local geometric theorem. It generates
new small additive tables; it does not read the giant ambient57 certificate.
Run with --verify for separate saved-table replay.
"""
import argparse
import json
import resource
import signal
import sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import fresh_ambient57_payment_frontier_20261018 as C

ROOT = C.ROOT
TAG = 'fresh_d17_node3_large_d_guard_upper_20261018'
N = 57


def memory_bytes():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak if sys.platform == 'darwin' else 1024*peak


def memory_guard():
    assert memory_bytes() < 1024**3, '1GiB memory limit reached'


def digest(path):
    h = sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def table(weights, saved=None):
    if saved is None:
        v = [0] + [None] * N
        for r, w in weights.items():
            for n in range(r, N + 1):
                if v[n-r] is not None:
                    z = v[n-r] + w
                    if v[n] is None or z > v[n]:
                        v[n] = z
        last = [None] + [next((r for r, w in weights.items()
            if r <= n and v[n-r] is not None and v[n] == v[n-r]+w), None)
            for n in range(1, N+1)]
    else:
        assert saved['weights'] == [[r, w] for r, w in weights.items()]
        v, last = saved['values'], saved['last']
    assert len(v) == len(last) == N+1 and v[0] == 0 and last[0] is None
    inequalities = predecessors = 0
    for n in range(1, N+1):
        possible = []
        for r, w in weights.items():
            if r <= n and v[n-r] is not None:
                assert v[n] is not None and v[n] >= v[n-r]+w
                inequalities += 1
                possible.append(r)
        if possible:
            r = last[n]
            assert r in possible and v[n] == v[n-r]+weights[r]
            predecessors += 1
        else:
            assert v[n] is None and last[n] is None
    return dict(weights=list(weights.items()), values=v, last=last), inequalities, predecessors


def run(saved=None):
    priors = {(tuple(r['counts']), r['k5'], r['a12']): r
        for r in saved['classes']} if saved else {}
    classes = []
    checks = predecessors = options = quarter_options = reserves = 0
    for cs, used, key in C.fixed_domain():
        memory_guard()
        # Only these fixed ranks can receive either proposed new coupon.
        if not (cs[4] or cs[5]):
            continue
        for k5 in range((N-used)//5+1):
            rem = N-used-5*k5
            if rem < 18:
                continue
            fixed = C.fixed_components(cs, key) + [('A', 5)]*k5
            fc = C.upper_value(fixed)
            old = priors.get((cs, k5, key[9]))
            if saved:
                assert old is not None
            weights = {r: C.item(r, key, k5) + (128*cs[4] if r == 4 else 0)
                + (128*cs[5] if r == 2 else 0) for r in C.ITEMS}
            tabs = []
            for no2 in (False, True):
                ws = {r:w for r,w in weights.items() if not no2 or r != 2}
                tab, a, b = table(ws, old['tables'][int(no2)] if old else None)
                checks += a
                predecessors += b
                tab['exclude2'] = no2
                tabs.append(tab)
            rows = []
            for large_rank in range(18, rem+1):
                left_rank = rem-large_rank
                # A whole D>=18 factor has no old cross/context interaction.
                # Reserve one such actual factor; all others stay in the table.
                ld = C.actual_base(('D', large_rank))
                ordinary = [(tabs[0]['values'][left_rank], 0, 0, left_rank)]
                for reserve in (1, 2):
                    oldguard = reserve == 2 and bool(cs[1])
                    newguard = bool(cs[0] and k5)
                    if left_rank < reserve or not (oldguard or newguard) or (reserve == 2 and not oldguard):
                        continue
                    for a2 in range((left_rank-reserve)//2+1):
                        rest = left_rank-reserve-2*a2
                        bonus = (36*a2*(a2-1) if oldguard else 0) + (64*a2*k5 if newguard else 0)
                        v = tabs[1]['values'][rest]
                        if v is not None:
                            ordinary.append((reserve*weights[1]+a2*weights[2]+v+bonus,
                                reserve, a2, rest))
                ordinary = [x for x in ordinary if x[0] is not None]
                assert ordinary
                options += len(ordinary)
                best, a1, a2, rest = max(ordinary)
                # Retain the current order17 alternate-half transfer, with
                # D5/A12 and A11/A1 witness types relaxed to their fixed fibers.
                delta17 = 16*k5*key[9] + (8*k5*left_rank if cs[2] else 0)
                upper = fc + ld + best + delta17
                quarter = None
                if cs[1] and left_rank >= 7:
                    for qa1 in range(1, left_rank-5):
                        for qa4 in range(1, (left_rank-qa1-2)//4+1):
                            for qa2 in range(1, (left_rank-qa1-4*qa4)//2+1):
                                qr = left_rank-qa1-4*qa4-2*qa2
                                v = tabs[1]['values'][qr]
                                if v is None:
                                    continue
                                bonus = (36*qa2*(qa2-1) if qa1 >= 2 else 0) + (64*qa2*k5 if cs[0] else 0) + 64*qa2*qa4
                                z = fc+ld+qa1*weights[1]+qa2*weights[2]+qa4*weights[4]+v+bonus+delta17
                                quarter_options += 1
                                if quarter is None or z > quarter[0]:
                                    quarter = (z, qa1, qa2, qa4, qr)
                    if quarter is not None:
                        upper = max(upper, quarter[0])
                rows.append(dict(reserved_D_rank=large_rank, residual_rank=left_rank,
                    ordinary_upper=F(fc+ld+best+delta17,64),
                    ordinary_explicit=[a1,a2,rest], quarter=quarter,
                    order17_delta_upper=F(delta17,64), upper=F(upper,64)))
                reserves += 1
            classes.append(dict(counts=cs, k5=k5, a12=key[9], key=key,
                remaining_rank=rem, tables=tabs, reserved_rows=rows,
                upper=max(r['upper'] for r in rows)))
    if saved:
        assert len(classes) == len(priors)
    maximum = max(r['upper'] for r in classes)
    memory_guard()
    return C.enc(dict(status='PASS_D17_NODE3_LARGE_D_GUARD_UPPER',
        sources={str(p.relative_to(ROOT)):digest(p) for p in
            (Path(__file__), Path(C.__file__), Path(C.C.__file__), Path(C.P.__file__))},
        scale=64, ambient_rank=N, complete_classes=len(classes),
        reserved_large_D_options=reserves, Bellman_inequalities=checks,
        attaining_predecessors=predecessors, ordinary_options=options,
        quarter_options=quarter_options, classes=classes,
        increasing_family_upper=maximum,
        all_increasing_forests_below_curvature=maximum < 364,
        all_increasing_forests_below_inherited_maximum=maximum < 372,
        geometric_theorem_proved=False, full_finite_maximum_recounted=False,
        scope='Complete upper for all forests with whole D>=18 and an A4/A13 or A2/A14 coupon; presence/A1 guards relaxed; current order17 increases included. D17 decreases ignored.'))


if __name__ == '__main__':
    # This macOS sandbox rejects lowering RLIMIT_AS/DATA. Check resident
    # storage after every small fixed class and read no giant JSON.
    signal.alarm(59)
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'
    assert args.verify or not path.exists()
    saved = json.loads(path.read_text()) if args.verify else None
    out = run(saved)
    if args.verify:
        assert saved == out
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['complete_classes'], out['increasing_family_upper'])
