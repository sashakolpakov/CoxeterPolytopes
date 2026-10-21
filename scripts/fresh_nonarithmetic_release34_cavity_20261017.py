#!/usr/bin/env python3
"""Bounded numerical endpoint-cluster scout for the deleted34 cavity.

Consumes frozen exact residual coordinates; no numerical result is an exact
Coxeter or completeness certificate. All old native contacts are tested.
"""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_release34_cavity_20261017'
PARENT = ROOT / 'results/fresh_nonarithmetic_release34_s6_boundary_20261017.json'
ALPHABET = list(range(2, 33)) + [48, 64, 128]
FRAMES = [[7, 12, 47], [8, 12, 52], [30, 47, 52]]
METRIC = np.diag([-1., 1., 1., 1.])


def mid(x):
    return float((Fraction(x[0]) + Fraction(x[1])) / 2)


def theta(m):
    return 0. if m == 2 else 2 * np.cos(np.pi / m)


def relation(x, tol=2e-7):
    if x < -2 - tol:
        return 'dotted'
    if abs(x + 2) <= tol:
        return 'parallel_or_unresolved'
    if abs(x) <= tol:
        return 2
    if x > 0 or x < -2:
        return None
    m = int(round(np.pi / np.arccos(-x / 2)))
    return m if m >= 3 and abs(x + theta(m)) <= tol else None


def setup():
    z = json.loads(PARENT.read_text())['after65']
    g = np.array([[mid(x) for x in row] for row in z['residual_gram']])
    g = (g + g.T) / 2
    vals, vecs = np.linalg.eigh(g)
    assert vals[0] < 0 < vals[1]
    transform = vecs @ np.diag(1 / np.sqrt(abs(vals)))
    assert np.max(abs(transform.T @ g @ transform - METRIC)) < 1e-12
    labels = [r['wall'] for r in z['projections']]
    cov = np.array([[mid(x) for x in r['covector']] for r in z['projections']])
    normals = cov @ transform @ METRIC
    h = np.linalg.solve(transform, np.array([mid(x) for x in z['projected_common_center']]))
    h /= np.sqrt(-h @ METRIC @ h)
    sides = [labels.index(i) for i in z['actual_sides']]
    return z, g, transform, labels, normals, h, sides


def endpoint_candidates(labels, normals, h):
    all_records, groups, counts = [], [], Counter()
    index = {x: i for i, x in enumerate(labels)}
    for group, frame in enumerate(FRAMES):
        f = normals[[index[i] for i in frame]]
        gram = f @ METRIC @ f.T
        gram[np.abs(gram) < 1e-12] = 0
        # Projected old7 has square4/3, not its native square2. This matters
        # for actual S6 endpoint existence, although Coxeter contacts below
        # continue to use the native normalization.
        for j in range(3):
            rational = float(Fraction(float(gram[j, j])).limit_denominator(100))
            assert abs(rational - gram[j, j]) < 1e-9
            gram[j, j] = rational
        inv = np.linalg.inv(gram)
        kernel = np.linalg.svd(f @ METRIC)[2][-1]
        square = kernel @ METRIC @ kernel
        assert square > 0
        kernel /= np.sqrt(square)
        accepted = []
        for ms in product(ALPHABET, repeat=3):
            c = -np.array([theta(m) for m in ms])
            # Actual finite endpoints at the first/second and first/third
            # old frame pairs are required; no other corner is assumed.
            if any(np.linalg.eigvalsh(np.array([
                [gram[0, 0], gram[0, j], c[0]],
                [gram[0, j], gram[j, j], c[j]], [c[0], c[j], 2.]]))[0] <= 1e-9
                   for j in [1, 2]):
                continue
            counts['finite_profile'] += 1
            p = c @ inv @ f
            rem = 2 - p @ METRIC @ p
            if rem < -1e-9:
                counts['no_real_norm_solution'] += 1
                continue
            if rem <= 1e-9:
                counts['degenerate_norm_solution'] += 1
                continue
            for sign in [-1, 1]:
                n = p + sign * np.sqrt(rem) * kernel
                assert abs(n @ METRIC @ n - 2) < 1e-7
                products = normals @ METRIC @ n
                assert max(abs(products[[index[i] for i in frame]] - c)) < 1e-7
                products[[index[i] for i in frame]] = c
                center = n @ METRIC @ h
                if center >= -1e-8:
                    counts['wrong_center'] += 1
                    continue
                if max(products) > 2e-7:
                    counts['obtuse'] += 1
                    continue
                relations = [relation(x) for x in products]
                if any(r is None or r == 'parallel_or_unresolved' for r in relations):
                    counts['noncoxeter_or_nearparallel'] += 1
                    continue
                counts['accepted_numeric'] += 1
                row = dict(index=len(all_records), group=group, frame=frame,
                           labels=list(ms), branch=sign, normal=n.tolist(),
                           center_product=float(center),
                           finite_contacts=[[i, r] for i, r in zip(labels, relations) if isinstance(r, int)],
                           minimum_dotted_margin=float(min([-2 - x for x, r in zip(products, relations) if r == 'dotted'], default=0)))
                all_records.append(row)
                accepted.append(row['index'])
        groups.append(accepted)
    return all_records, groups, counts


def ledger(normals, labels):
    """Numerical complete triple-intersection diagnostic, not certification."""
    cov = normals @ METRIC
    cov = cov / np.linalg.norm(cov, axis=1)[:, None]
    chart = cov.sum(axis=0)
    supports = list(combinations(range(len(labels)), 3))
    matrices = np.array([np.vstack([cov[list(s)], chart]) for s in supports])
    dets = np.linalg.det(matrices)
    nonsingular = np.abs(dets) > 1e-10
    vec = np.zeros((len(supports), 4))
    vec[nonsingular] = np.linalg.solve(matrices[nonsingular], np.tile([0., 0., 0., -1.], (sum(nonsingular), 1))[:, :, None])[:, :, 0]
    vec[nonsingular] /= np.linalg.norm(vec[nonsingular], axis=1)[:, None]
    products = vec @ cov.T
    good = np.where(nonsingular & (products.max(axis=1) < 2e-10))[0]
    rays, bysupport = [], {}
    for k in good:
        active = tuple(np.where(abs(products[k]) < 2e-10)[0])
        if active in bysupport:
            continue
        bysupport[active] = len(rays)
        v = vec[k]
        rays.append((active, v, float(v @ METRIC @ v)))
    pair_index = defaultdict(list)
    for i, (support, v, sq) in enumerate(rays):
        for pair in combinations(support, 2):
            pair_index[pair].append(i)
    counts = Counter()
    edges = []
    for pair, ends in pair_index.items():
        if len(ends) != 2:
            counts['degenerate_face'] += 1
            continue
        a, b = [rays[i] for i in ends]
        finite = (a[2] < -1e-12) + (b[2] < -1e-12)
        ab = float(a[1] @ METRIC @ b[1])
        det = a[2] * b[2] - ab * ab
        if finite == 2:
            kind = 'compact'
        elif finite == 1:
            kind = 'one_open'
        elif ab < -1e-12 and det < -1e-15:
            kind = 'two_open'
        else:
            kind = 'outside'
        counts[kind] += 1
        edges.append(dict(sides=[labels[j] for j in pair], kind=kind))
    return dict(rays=len(rays), finite=sum(r[2] < -1e-12 for r in rays),
                exterior=sum(r[2] > 1e-12 for r in rays), edges=dict(counts),
                open_ends=counts['one_open'] + 2 * counts['two_open'],
                noncompact_edges=[e for e in edges if e['kind'] != 'compact'],
                singular_triples=int(sum(~nonsingular)))


def compute():
    z, g, transform, labels, normals, h, sides = setup()
    records, groups, counts = endpoint_candidates(labels, normals, h)
    print('endpoint candidates', [len(x) for x in groups], dict(counts), flush=True)
    compatibility, adjacent, finite_adjacent = [], set(), set()
    for a, b in combinations(records, 2):
        if a['group'] == b['group']:
            continue
        x = float(np.array(a['normal']) @ METRIC @ np.array(b['normal']))
        rel = relation(x)
        if rel is None or rel == 'parallel_or_unresolved':
            continue
        adjacent.add((a['index'], b['index']))
        if isinstance(rel, int):
            finite_adjacent.add((a['index'], b['index']))
        compatibility.append(dict(pair=[a['index'], b['index']], product=x, relation=rel))
    print('compatible pairs', len(compatibility), 'finite pairs', len(finite_adjacent), flush=True)
    batches = []
    selections = sorted(finite_adjacent)
    neighbors = defaultdict(set)
    for a, b in finite_adjacent:
        neighbors[a].add(b)
        neighbors[b].add(a)
    triples = set()
    for a, near in neighbors.items():
        for b, c in combinations(sorted(near), 2):
            t = tuple(sorted([a, b, c]))
            if len({records[i]['group'] for i in t}) == 3 and tuple(sorted([b, c])) in adjacent:
                triples.add(t)
    selections += sorted(triples)
    print('connected finite-contact batch selections', len(selections), flush=True)
    for batch in selections:
        ns = np.vstack([normals[sides], [records[i]['normal'] for i in batch]])
        labs = [labels[i] for i in sides] + [1000+i for i in batch]
        score = ledger(ns, labs)
        batches.append(dict(candidates=list(batch), boundary=score))
    batches.sort(key=lambda r: (r['boundary']['open_ends'], r['boundary']['exterior'], -r['boundary']['finite']))
    return dict(status='BOUNDED_NUMERICAL_CAVITY_CLUSTER_DIAGNOSTIC',
                sources={str(PARENT.relative_to(ROOT)): sha256(PARENT.read_bytes()).hexdigest()},
                alphabet=ALPHABET, frames=FRAMES, old_labels=labels,
                residual_basis=z['residual_basis'], transform=transform.tolist(),
                counts=dict(counts), candidates=records, groups=groups,
                compatible_pair_count=len(compatibility),
                compatible_finite_pairs=[r for r in compatibility if isinstance(r['relation'], int)], batches=batches,
                baseline=ledger(normals[sides], [labels[i] for i in sides]),
                scope='Finite bounded three-contact endpoint families; score only batches connected by finite new-new contacts. Numeric tolerance gates are discovery, not exact acceptance; other labels and contact patterns remain open.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT / 'results' / (TAG + '.json')
    out = compute()
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ['status', 'counts', 'baseline']}), flush=True)
    if out['batches']:
        print('best batch', json.dumps(out['batches'][0]), flush=True)


if __name__ == '__main__':
    main()
