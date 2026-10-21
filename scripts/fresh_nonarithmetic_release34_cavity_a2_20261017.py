#!/usr/bin/env python3
"""Bounded seven-profile cavity-cluster discovery, all native gates retained."""
import argparse
import json
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import numpy as np
import fresh_nonarithmetic_release34_cavity_20261017 as B

ROOT = B.ROOT
TAG = 'fresh_nonarithmetic_release34_cavity_a2_20261017'
INPUT = ROOT / 'results/fresh_nonarithmetic_release34_cavity_20261017_white_input.json'
ALPHABET = list(range(2, 17)) + [24, 32, 64]


def build_input():
    import fresh_nonarithmetic_s6_release34_exact_20261017 as X
    z = X.geometry()
    n, pair, h, s = z['normals'], z['pair'], z['common_center'], X.S
    old = X.D.E.labels()
    def literal(i, j):
        if i == j:
            return X.I(2)
        if max(i, j) < 55 and old[i][j]:
            return -X.D.theta(old[i][j])
        if i in X.NEW or j in X.NEW:
            return X.I(0)  # All four new walls are right to the six whites.
        return X.U.rnd(pair(n[i], n[j]))
    return dict(white_support=s,
                white_gram=[[literal(i, j).enc() for j in s] for i in s],
                old_white_contacts={str(i): [literal(i, j).enc() for j in s] for i in sorted(n)},
                center_white_products=[X.U.rnd(pair(h, n[i])).enc() for i in s],
                source_sha256=sha256(Path(X.__file__).read_bytes()).hexdigest())


def compute():
    z, g, transform, labels, normals, h, sides = B.setup()
    inp = json.loads(INPUT.read_text())
    wg = np.array([[B.mid(x) for x in row] for row in inp['white_gram']])
    wi = np.linalg.inv(wg)
    oldw = np.array([[B.mid(x) for x in inp['old_white_contacts'][str(i)]] for i in labels])
    hw = np.array([B.mid(x) for x in inp['center_white_products']])
    hw /= np.sqrt(-B.mid(z['projected_center_square']))
    profiles = [np.zeros(6)]
    descriptions = [dict(kind='right_S6')]
    for j in [4, 5]:
        for m in [3, 4, 5]:
            c = np.zeros(6)
            c[j] = -B.theta(m)
            profiles.append(c)
            descriptions.append(dict(kind='one_A2_bond', white=inp['white_support'][j], label=m))
    index = {x: i for i, x in enumerate(labels)}
    records, counts = [], Counter()
    for profile, white in enumerate(profiles):
        coeff = wi @ white
        norm = 2 - white @ coeff
        offsets = oldw @ coeff
        for group, frame in enumerate(B.FRAMES):
            inds = [index[i] for i in frame]
            f = normals[inds]
            gram = f @ B.METRIC @ f.T
            inv = np.linalg.inv(gram)
            kernel = np.linalg.svd(f @ B.METRIC)[2][-1]
            kernel /= np.sqrt(kernel @ B.METRIC @ kernel)
            for ms in product(ALPHABET, repeat=3):
                native = -np.array([B.theta(m) for m in ms])
                c = native - offsets[inds]
                if any(np.linalg.eigvalsh(np.array([
                    [gram[0, 0], gram[0, j], c[0]],
                    [gram[0, j], gram[j, j], c[j]], [c[0], c[j], norm]]))[0] <= 1e-9
                       for j in [1, 2]):
                    continue
                counts['finite_profile'] += 1
                p = c @ inv @ f
                rem = norm - p @ B.METRIC @ p
                if rem < 1e-10:
                    counts['no_real_or_degenerate'] += 1
                    continue
                for sign in [-1, 1]:
                    n = p + sign * np.sqrt(rem) * kernel
                    products = normals @ B.METRIC @ n + offsets
                    assert max(abs(products[inds] - native)) < 1e-7
                    products[inds] = native
                    center = n @ B.METRIC @ h + coeff @ hw
                    if center >= -1e-8:
                        counts['wrong_center'] += 1
                        continue
                    if max(products) > 2e-7:
                        counts['obtuse'] += 1
                        continue
                    relations = [B.relation(x) for x in products]
                    if any(r is None or r == 'parallel_or_unresolved' for r in relations):
                        counts['noncoxeter_or_nearparallel'] += 1
                        continue
                    counts['accepted_numeric'] += 1
                    records.append(dict(index=len(records), profile=profile, group=group,
                        frame=frame, labels=list(ms), branch=sign, normal=n.tolist(),
                        white=white.tolist(), center_product=float(center),
                        finite_contacts=[[i, r] for i, r in zip(labels, relations) if isinstance(r, int)]))
        print('profile', profile, 'cumulative candidates', len(records), flush=True)
    if not records:
        return dict(status='EMPTY_BOUNDED_NUMERICAL_PROFILE_BANK', counts=dict(counts))
    n = np.array([r['normal'] for r in records])
    w = np.array([r['white'] for r in records])
    groups = np.array([r['group'] for r in records])
    finite, neighbors = [], defaultdict(set)
    for start in range(0, len(records), 256):
        dots = n[start:start+256] @ B.METRIC @ n.T + w[start:start+256] @ wi @ w.T
        mask = (dots > -2+1e-7) & (dots < 2e-7)
        for ii, j in np.argwhere(mask):
            i = int(start + ii)
            j = int(j)
            if i >= j or groups[i] == groups[j]:
                continue
            x = float(dots[ii, j])
            rel = B.relation(x)
            if not isinstance(rel, int):
                continue
            finite.append(dict(pair=[i, j], relation=rel, product=x))
            neighbors[i].add(j)
            neighbors[j].add(i)
    print('finite compatible pairs', len(finite), flush=True)
    def dot(i, j):
        return float(n[i] @ B.METRIC @ n[j] + w[i] @ wi @ w[j])
    selections = {tuple(r['pair']) for r in finite}
    for a, near in neighbors.items():
        for b, c in combinations(sorted(near), 2):
            t = tuple(sorted([a, b, c]))
            if len(set(groups[list(t)])) != 3:
                continue
            rel = B.relation(dot(b, c))
            if rel is not None and rel != 'parallel_or_unresolved':
                selections.add(t)
    batches = []
    for ids in sorted(selections):
        ns = np.vstack([normals[sides], n[list(ids)]])
        ls = [labels[i] for i in sides] + [1000+i for i in ids]
        batches.append(dict(candidates=list(ids), boundary=B.ledger(ns, ls)))
    batches.sort(key=lambda r: (r['boundary']['open_ends'], r['boundary']['exterior']))
    baseline = B.ledger(normals[sides], [labels[i] for i in sides])
    assert baseline['open_ends'] == 17 and baseline['finite'] == 23
    return dict(status='BOUNDED_NUMERICAL_SEVEN_PROFILE_CAVITY_CLUSTER_DIAGNOSTIC',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(B.__file__), B.PARENT, INPUT]},
        profiles=descriptions, alphabet=ALPHABET, frames=B.FRAMES, counts=dict(counts),
        candidates=records, finite_pairs=finite, batches=batches, baseline=baseline,
        scope='Seven possible actual-facet white profiles; bounded framing labels and only three selected two-endpoint groups. Native offsets retained in every old gate and new-new gate. Numerical only; not an exhaustive completion or exact geometry claim.')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    p.add_argument('--rebuild-input', action='store_true')
    args = p.parse_args()
    if args.rebuild_input or not INPUT.exists():
        out = build_input()
        if args.verify:
            assert json.loads(INPUT.read_text()) == out
        else:
            INPUT.write_text(json.dumps(out, indent=2)+'\n')
    out = compute()
    path = ROOT/'results'/(TAG+'.json')
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['counts'], flush=True)
    if out.get('batches'):
        print('best', out['batches'][0], flush=True)


if __name__ == '__main__':
    main()
