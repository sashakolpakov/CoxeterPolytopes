#!/usr/bin/env python3
"""Exact higher-dimensional cyclic-base triage before Coxeter label search."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_highbase_construction_20261014_cyclic'
OUT = ROOT/f'results/{TAG}.json'


def compositions(n, total, prefix=()):
    if n == 1:
        yield prefix+(total,)
    else:
        for x in range(1, total-n+2):
            yield from compositions(n-1, total-x, prefix+(x,))


def base(n, d):
    vertices = [tuple(t**k for k in range(1,d+1)) for t in range(n)]
    facets, supports = [], []
    for f in combinations(range(n), d):
        polynomial = [1]
        for root in f:
            nxt = [0]*(len(polynomial)+1)
            for k,a in enumerate(polynomial):
                nxt[k] -= root*a
                nxt[k+1] += a
            polynomial = nxt
        values = [sum(a*t**k for k,a in enumerate(polynomial)) for t in range(n)]
        off = [values[t] for t in range(n) if t not in f]
        if min(off) > 0 or max(off) < 0:
            sign = 1 if min(off)>0 else -1
            assert [t for t,v in enumerate(values) if v == 0] == list(f)
            assert sign*sum(values) > 0  # vertex average is strictly inside each facet
            facets.append(f)
            supports.append([sign*a for a in polynomial])
    faces = {()}
    for f in facets:
        for k in range(1,d+1):
            faces.update(combinations(f,k))
    missing = []
    for k in range(2,d+2):
        for a in combinations(range(n),k):
            if a not in faces and all(a[:i]+a[i+1:] in faces for i in range(k)):
                missing.append(a)
    assert set(Counter(e for f in facets for e in combinations(f,d-1)).values()) == {2}
    return dict(dimension=d, vertex_count=n, moment_curve_vertices=vertices,
                dual_facets=facets, exact_support_polynomials=supports,
                minimal_nonfaces=missing)


def run():
    records = []
    for n,d in [(10,4),(11,5),(12,6)]:
        b = base(n,d)
        total = n+9-d
        candidates = list(compositions(n,total))
        assert len(candidates) == comb(total-1,n-1)
        histogram = Counter()
        valid = []
        for w in candidates:
            order = max(sum(w[i] for i in m) for m in b['minimal_nonfaces'])
            histogram[order] += 1
            if order <= 5:
                valid.append(w)
        records.append(dict(base=b, dimension_after_wedges=9, facets_after_wedges=total,
                            positive_weight_compositions=len(candidates),
                            maximum_missing_order_histogram=dict(histogram),
                            all_admissible_weights=valid))
        print(n,d,len(candidates),'admissible',len(valid),flush=True)
    return json.loads(json.dumps(dict(
        status='EXACT_CYCLIC_BASE_MISSING_FACE_FILTER',
        sources={f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()},
        cases=records,
        scope='All positive multiwedge weights raising the three displayed cyclic bases to dimension nine. '
              'Exact convex moment-curve facets and complete minimal nonfaces. '
              'A missing face of order above five excludes compact Coxeter realization; other bases are not classified.')))


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run()
    if a.verify:
        assert out==json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
