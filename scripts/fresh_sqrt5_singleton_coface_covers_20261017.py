#!/usr/bin/env python3
"""Exact singleton mixtures on every exported individual coface slice."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_singleton_coface_covers_20261017'
SOURCE = 'fresh_sqrt5_labelled_coface_base_packing_20261017'


def run():
    src = json.loads((ROOT / f'results/{SOURCE}.json').read_text())
    assert src['status'] == 'PASS_EXACT_LABELLED_COFACE_SUBSET_PACKING_DIAGNOSTIC'
    alpha = {i: F(x) for i, x in src['face_multipliers']}
    faces = {}
    for b in src['actual_bases']:
        p, f = b['provenance'], b['face']
        assert p not in faces or faces[p] == f
        faces[p] = f
    rows, hist = [], Counter()
    for s in src['coface_slice_demands']:
        masses = []
        for mask, terms in s['block_coefficients']:
            assert mask.bit_count() >= 2
            assert mask < 1 << len(s['vertex_roots'])
            d = sum((F(c) * alpha.get(faces[p], F()) for p, c in terms), F())
            assert d >= 0
            if d:
                masses.append([mask, str(d)])
        assert masses == s['block_marginals']
        total = sum((F(d) for _, d in masses), F())
        capacity = s['coface_multiplicity'] * alpha.get(s['simplex_face'], F())
        assert total == F(s['total_block_mass'])
        assert capacity == F(s['capacity']) and total <= capacity
        if total:
            hist[len(masses)] += 1
        rows.append(dict(slice=s['index'], simplex_face=s['simplex_face'],
                         singleton_mixture=masses, empty_mass=str(capacity-total)))
    assert len(rows) == 34520 and sum(hist.values()) == 870
    paths = [f'scripts/{TAG}.py', f'results/{SOURCE}.json']
    return dict(status='PASS_EXACT_SINGLETON_COVERS_ON_ALL_INDIVIDUAL_COFACE_SLICES',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                counts=dict(slices=len(rows), positive_slices=sum(hist.values()),
                            positive_block_count_histogram=sorted(hist.items())),
                covers=rows, input_bound=src['input_bound'], optimizer_invocations=0,
                scope='Every exported averaged individual-coface marginal is a mixture of empty and singleton-base configurations. Any compatibility relaxation retaining these individual options admits these marginals. This does not certify common ambient realizability or simultaneous cross-coface compatibility; the new cross-coface cap constraints violate the same old witness.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT/f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'])
