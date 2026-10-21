#!/usr/bin/env python3
"""Coherent multi-cap configurations and exact spacelike volume certificates."""
import argparse
from fractions import Fraction as Q
from itertools import combinations
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.spatial import HalfspaceIntersection

from scripts.fresh_h9_joint_completion_20261013 import setup, read, ROOT
from scripts.fresh_h8_threefree_verify_20261013 import gate as old_gate, pair, qsign
from scripts.fresh_h8_alternate_caps_20261013 import real
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, Z, O, parse, enc, scale, dot, new_pair_gate)

SOURCE = 'data/enumerations/fresh_h9_quartic_completion_20261013_targeted.json'
OUTPUT = ROOT / 'data/enumerations/fresh_h9_quartic_completion_20261013_geometry.json'


def numeric(field, x):
    if field is K:
        return (mp.mpf(x[0].numerator)/x[0].denominator +
                mp.mpf(x[1].numerator)/x[1].denominator*(1+mp.sqrt(5))/2)
    return numeric(field.base, x[0])+mp.sqrt(numeric(field.base, field.d))*numeric(field.base, x[1])


def run(limit=8):
    G, I, center, old, oldpolys, oldcn, In, hn, _ = setup()
    source = read(SOURCE)
    E = Quadratic(K, parse(source['anchor_radicand']))
    pool = [(r['label'], tuple(map(Q, r['coefficients_descending'])),
             *map(Q, r['largest_root_interval'])) for r in read(source['angle_source'])]
    fields, coefficients, pairings, keep, cut = [], [], [], [], []
    cap_rows = list(source['caps'])
    references = [{'source': SOURCE, 'index': i} for i in range(len(cap_rows))]
    if source.get('required_cap_source'):
        other = read(source['required_cap_source'])
        ids = source['required_cap_indices']
        assert other['anchor_radicand'] == source['anchor_radicand']
        cap_rows = [other['caps'][i] for i in ids]+cap_rows
        references = [{'source': source['required_cap_source'], 'index': i} for i in ids]+references
    for c in cap_rows:
        F = Quadratic(E, E.parse(c['discriminant_over_E']))
        fields.append(F)
        coefficients.append(tuple(map(F.parse, c['seed_coefficients'])))
        pairings.append(tuple(map(F.parse, c['seed_pairings'])))
        keep.append(set(c['retained_nodes']))
        cut.append(set(c['cut_target_indices']))
    mp.mp.dps = 90
    cn = np.array([[float(numeric(F, x)) for x in c] for F, c in zip(fields, coefficients)])
    checks, adjacency = {}, [set() for _ in fields]
    for i, j in combinations(range(len(fields)), 2):
        row = new_pair_gate(fields[i], pairings[i], fields[j].d, coefficients[j], pool)
        checks[f'{i},{j}'] = row
        if row['compatible']:
            adjacency[i].add(j)
            adjacency[j].add(i)
    seeds = []
    for i in range(len(fields)):
        if 5 not in cut[i]:
            continue
        for j in adjacency[i]:
            if 22 not in cut[j]:
                continue
            common = keep[i] & keep[j]
            score = (len(common), len(cut[i] | cut[j]), -i, -j)
            seeds.append((score, tuple(sorted((i, j)))))
    seeds.sort(reverse=True)
    cases, seen = [], set()
    for _, seed in seeds:
        selected = list(seed)
        retained = set.intersection(*(keep[i] for i in selected))
        eligible = set.intersection(*(adjacency[i] for i in selected))
        # Preserve the coherent shared core while filling with further caps.
        while eligible:
            allowed = [i for i in eligible if retained <= keep[i]]
            if not allowed:
                break
            covered = set.union(*(cut[i] for i in selected))
            j = max(allowed, key=lambda j: (len(cut[j]-covered),
                                           len(adjacency[j] & eligible), -j))
            selected.append(j)
            eligible &= adjacency[j]
        selected.sort()
        key = tuple(selected), tuple(sorted(retained))
        if key in seen:
            continue
        seen.add(key)
        oldids = sorted(retained)
        oldchecks = {}
        for i, j in combinations(oldids, 2):
            row = old_gate(old[i], old[j], oldpolys)
            assert row['compatible']
            oldchecks[f'{i},{j}'] = row
        assert source['anchor'] in oldids
        A = np.vstack((np.eye(10), oldcn[oldids], cn[selected]))
        half = np.c_[A[:, 9, None]-A[:, :9], -A[:, 9]]
        half /= np.linalg.norm(half[:, :9], axis=1)[:, None]
        vv = HalfspaceIntersection(half, hn[:9], qhull_options='Qx').intersections
        ys = np.c_[vv, 1-vv.sum(axis=1)]
        qq = np.einsum('bi,ij,bj->b', ys, In, ys)
        order = sorted(np.flatnonzero(qq > 1e-8), key=lambda j: -qq[j])
        escapes, masks = [], set()
        for vi in order:
            certificate = None
            for epsilon in (1e-3, 1e-5, 1e-7):
                y0 = ys[vi]+epsilon*hn
                for denominator in (10**9, 10**12, 10**15):
                    y = tuple(scale(O, Q(round(float(x)*denominator), denominator)) for x in y0)
                    if not all(K.sign(x) > 0 for x in y):
                        continue
                    norm, cp = pair(y, I, y), pair(y, I, center)
                    if K.sign(norm) <= 0 or K.sign(cp) >= 0:
                        continue
                    if any(qsign((dot(K, old[i]['a'], y), dot(K, old[i]['b'], y)), old[i]['d']) <= 0 for i in oldids):
                        continue
                    values = [dot(F, c, tuple(F.lift(E.lift(t)) for t in y))
                              for F, c in zip(fields, coefficients)]
                    signs = [F.sign(value) for F, value in zip(fields, values)]
                    if any(signs[i] <= 0 for i in selected):
                        continue
                    cutters = [i for i, sign in enumerate(signs) if sign < 0]
                    oldcutters = [i for i, n in enumerate(old)
                                  if qsign((dot(K, n['a'], y), dot(K, n['b'], y)), n['d']) < 0]
                    certificate = {'y': list(map(enc, y)), 'norm': enc(norm),
                                   'center_pairing': enc(cp),
                                   'new_dictionary_cutters': cutters,
                                   'old_dictionary_cutters': oldcutters}
                    break
                if certificate:
                    break
            if certificate and tuple(certificate['new_dictionary_cutters']) not in masks:
                masks.add(tuple(certificate['new_dictionary_cutters']))
                escapes.append(certificate)
                # A whole-dictionary escape already certifies this chamber;
                # do not inspect thousands of other vertices for duplicate masks.
                if not certificate['new_dictionary_cutters']:
                    break
            if len(escapes) >= 6:
                break
        case = {'selected_new_caps': selected, 'retained_old_caps': oldids,
                'facet_count': 10+len(oldids)+len(selected),
                'covered_saved_targets': sorted(set.union(*(cut[i] for i in selected),
                                                         {8})),
                'old_pair_checks': oldchecks, 'numeric_vertices': len(ys),
                'numeric_spacelike_vertices': int((qq > 1e-8).sum()),
                'escapes': escapes,
                'status': 'EXACT_ESCAPES' if escapes else 'UNRESOLVED_VOLUME_CANDIDATE'}
        cases.append(case)
        result = {'source': SOURCE, 'cap_references': references,
                  'all_new_pair_checks': checks,
                  'compatible_cross_gap_pairs': len(seeds),
                  'cases': cases,
                  'scope': 'Exact compatible configurations with retained coherent cores. '
                           'Numerical vertices propose exact strict future spacelike witnesses; '
                           'absence of an escape would not certify finite volume.'}
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
        print({'case': len(cases)-1, 'new': selected, 'old': len(oldids),
               'facets': case['facet_count'], 'vertices': len(ys),
               'escapes': len(escapes)}, flush=True)
        if not escapes or len(cases) >= limit:
            break


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=8)
    p.add_argument('--source', default=SOURCE)
    p.add_argument('--output')
    a = p.parse_args()
    SOURCE = a.source
    if a.output:
        OUTPUT = ROOT/a.output
    run(a.limit)
