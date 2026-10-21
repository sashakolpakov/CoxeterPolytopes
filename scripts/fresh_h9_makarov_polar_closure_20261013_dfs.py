#!/usr/bin/env python3
"""Bounded numerical clique/coverage search in the polar bank."""
import argparse
import json
import time
import numpy as np
from fresh_h9_makarov_polar_closure_20261013 import setup, ROOT, TAG


def bits(row):
    return sum(1 << int(i) for i in np.flatnonzero(row))


def indexes(mask):
    while mask:
        a = mask & -mask
        yield a.bit_length()-1
        mask ^= a


def run(seconds=90):
    start = time.time()
    source, Q, allR, c, refs, points = setup()
    R = allR[27:]
    cov = R@Q@points.T > 1e-7
    pp = R@Q@R.T
    ok = pp < -1+2e-6
    for m in range(2, 41):
        ok |= abs(pp+np.cos(np.pi/m)) < 2e-6
    np.fill_diagonal(ok, False)
    covers = list(map(bits, cov))
    compatible = list(map(bits, ok))
    cutters = list(map(bits, cov.T))
    full = (1 << len(R))-1
    uncovered = (1 << len(points))-1
    failures = set()
    nodes = 0
    best = len(points)
    answer = None
    timed_out = False
    trace = []

    def dfs(available, missing, selected):
        nonlocal nodes, best, answer, timed_out
        nodes += 1
        if nodes % 1024 == 0 and time.time()-start > seconds:
            timed_out = True
            return False
        if not missing:
            answer = selected
            return True
        if (available, missing) in failures:
            return False
        count = missing.bit_count()
        if count < best:
            best = count
        options = None
        witness = None
        smallest = len(R)+1
        for j in indexes(missing):
            opts = available & cutters[j]
            n = opts.bit_count()
            if not n:
                trace.append({'selected': selected, 'witness_index': j, 'options': []})
                failures.add((available, missing))
                return False
            if n < smallest:
                smallest, options, witness = n, opts, j
                if n == 1:
                    break
        order = sorted(indexes(options), key=lambda i: (covers[i] & missing).bit_count(), reverse=True)
        trace.append({'selected': selected, 'witness_index': witness, 'options': order})
        for i in order:
            if dfs(available & compatible[i], missing & ~covers[i], selected+[i]):
                return True
            if timed_out:
                return False
            available &= ~(1 << i)
        failures.add((available, missing))
        return False

    dfs(full, uncovered, [])
    out = {'status': 'NUMERICAL_COMPATIBLE_POLAR_COVER' if answer is not None else
                     'TIME_LIMIT' if timed_out else 'NUMERICAL_EXHAUSTED_POLAR_COVER',
           'nodes': nodes, 'cache_size': len(failures), 'best_uncovered_count': best,
           'seconds': time.time()-start,
           'branch_trace': trace,
           'selected': None if answer is None else [refs[27+i] for i in answer],
           'selected_normals': None if answer is None else R[answer].tolist(),
           'scope': 'Numerical proposal search only. Exhaustion here is not an exact Coxeter nonexistence certificate.'}
    (ROOT/'results'/f'{TAG}_dfs.json').write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], nodes, 'nodes;', 'best', best, 'seconds', out['seconds'], flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--seconds', type=int, default=90)
    run(ap.parse_args().seconds)
