#!/usr/bin/env python3
"""Actual right-gap coupling, permitting at most two declared cap replacements."""
import argparse
from itertools import combinations, product
import json
import numpy as np
import mpmath as mp
from scripts.fresh_h9_gap_completion_20261013_base import *


def run(limit=40, maximum_fixed=2, focus=0, drop_budget=0, stage='45', anchor_choice='right'):
    G, I, h, ps, target_source, target_case = configuration(stage)
    anchor_ref = ANCHOR if anchor_choice == 'right' else {
        'source': 'data/enumerations/fresh_h9_replaced_core_20261013_polar_caps.json', 'index': 0}
    anchor_name = 'cap2' if anchor_choice == 'right' else 'required0'
    E, a = canonical_anchor(anchor_ref)
    mp.mp.dps = 90
    an = np.array([float(value(E, t)) for t in a])
    bn = np.array([float(value(K, t)) for t in ps['old348']['c']])
    In = np.array([[float(value(K, t)) for t in row] for row in I])
    hn = np.array([float(value(K, t)) for t in h])
    ys = [tuple(map(parse, t['certificate']['y']))
          for t in read(target_source)['cases'][target_case]['geometry']['escapes']]
    Yn = np.array([[float(value(K, t)) for t in y] for y in ys])
    names = list(ps)
    removable = [n for n in names if (n.startswith(('required', 'gap')) or n in ('cap2', 'cap5', 'swap_new6'))
                 and n != anchor_name]
    protected = [i for i, n in enumerate(names) if n not in removable]
    cn = np.array([[float(value(ps[n]['field'], t)) for t in ps[n]['c']] for n in names])
    qn = np.array([float(value(K, ps[n]['q'])) for n in names])
    allowed = np.array([np.cos(np.pi/m)**2 for m in (2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30)])
    proposals = []
    counts = {'linear_families': 0, 'seed_center_roots': 0, 'focus_cutters': 0,
              'coherent_by_drops': {}}
    for free in combinations(range(10), 3):
        best = max((abs(an[j]*bn[k]-an[k]*bn[j]), i, j, k)
                   for i in free for j, k in combinations([z for z in free if z != i], 2))
        det, u, v, w = best
        if det < 1e-9:
            continue
        det = an[v]*bn[w]-an[w]*bn[v]
        remaining = [i for i in range(10) if i not in free]
        specs = [((), ())]
        for count in range(1, maximum_fixed+1):
            specs += [(s, labels) for s in combinations(remaining, count)
                      for labels in product((3, 5), repeat=count)]
        fixed = np.zeros((len(specs), 10))
        for j, (positions, labels) in enumerate(specs):
            for i, label in zip(positions, labels):
                fixed[j, i] = 1 if label == 3 else float(value(K, P))
        line1 = np.zeros(10)
        line1[u] = 1
        line1[v] = (an[w]*bn[u]-an[u]*bn[w])/det
        line1[w] = (an[u]*bn[v]-an[v]*bn[u])/det
        A = line1 @ In @ line1
        if abs(A) < 1e-9:
            continue
        for am, pm in product((2, 3, 5), (2, 10)):
            M, N = {2: 0, 3: 1, 5: float(value(K, P))}[am], 0 if pm == 2 else 1
            r0, r1 = M-fixed@an, N-fixed@bn
            line0 = fixed.copy()
            line0[:, v] = (r0*bn[w]-an[w]*r1)/det
            line0[:, w] = (an[v]*r1-r0*bn[v])/det
            B = 2*line0 @ In @ line1
            C = np.einsum('bi,ij,bj->b', line0, In, line0)-2
            D = B*B-4*A*C
            counts['linear_families'] += len(specs)
            for branch in (-1, 1):
                t = (-B+branch*np.sqrt(np.maximum(D, 0)))/(2*A)
                ws = line0+t[:, None]*line1
                ratio = ws*ws/4
                near = np.min(abs(ratio[:, :, None]-allowed), axis=2) < 1e-7
                ok = (D > 1e-8) & np.all((ws >= -1e-7) & ((ws > 2+1e-7) | near), axis=1)
                ok &= ws @ In @ hn < -1e-7
                counts['seed_center_roots'] += int(ok.sum())
                cuts = ws @ In @ Yn.T
                ok &= cuts[:, focus] > 1e-7
                counts['focus_cutters'] += int(ok.sum())
                for j in np.flatnonzero(ok):
                    pairs = -cn @ ws[j]
                    rr = pairs*pairs/(2*qn)
                    good = (pairs <= 1e-7) & ((rr > 1+1e-7) |
                            (np.min(abs(rr[:, None]-allowed), axis=1) < 1e-7))
                    if not all(good[protected]):
                        continue
                    dropped = [names[i] for i in np.flatnonzero(~good)]
                    counts['coherent_by_drops'][str(len(dropped))] = counts['coherent_by_drops'].get(str(len(dropped)), 0)+1
                    if len(dropped) > drop_budget:
                        continue
                    positions, labels = specs[j]
                    score = (-len(dropped), int((cuts[j] > 1e-7).sum()), min(cuts[j]), cuts[j, focus])
                    proposals.append((score, (u, v, w, positions, labels, am, pm, branch)))
    proposals.sort(reverse=True)
    print({'counts': counts, 'proposals': len(proposals)}, flush=True)
    pool = angle_pool(ANGLE_SOURCE)
    caps, rejected, seen = [], [], set()
    for score, spec in proposals:
        if len(caps)+len(rejected) >= limit:
            break
        try:
            row, record = exact_cap(spec, G, I, h, ps, pool, anchor_ref)
        except AssertionError:
            rejected.append({'spec': spec, 'reason': 'exact_norm_seed_gate'})
            continue
        nk = tuple(round(float(value(row['field'], p)), 10) for p in row['p'])
        if nk in seen:
            continue
        seen.add(nk)
        checks = {n: relation(row, ps[n], pool) for n in names}
        dropped = [n for n in names if checks[n]['kind'] not in ('finite', 'dotted')]
        if len(dropped) > drop_budget or not set(dropped) <= set(removable):
            rejected.append({'spec': spec, 'reason': 'exact_core_gate', 'dropped': dropped})
            continue
        F = row['field']
        signs = [F.sign(dot(F, row['c'], tuple(lift_K(F, t) for t in y))) for y in ys]
        assert signs[focus] < 0
        record.update({'anchor_name': anchor_name, 'pair_checks': checks, 'dropped': dropped,
                       'target_signs': signs, 'cuts': [i for i, s in enumerate(signs) if s < 0]})
        caps.append(record)
        print({'accepted': len(caps), 'dropped': dropped, 'anchor_label': record['anchor_label'],
               'polar_label': record['polar_label'], 'cuts': record['cuts']}, flush=True)
    result = {'status': 'EXACT_NEW_ANCHOR_GAP_CAPS', 'source': target_source, 'source_case': target_case,
              'configuration': stage, 'cap_name_prefix': 'gap' if stage == '45' else 'add',
              'anchor_reference': anchor_ref, 'focus_target': focus, 'drop_budget': drop_budget,
              'removable': removable, 'angle_source': ANGLE_SOURCE,
              'proposal_counts': counts, 'caps': caps, 'rejected': rejected,
              'scope': 'Bounded full-new-anchor construction; exact accepted caps and '
                       'declared replacements. Numerical screening is not exhaustive.'}
    suffix = '' if stage == '45' and anchor_choice == 'right' else f'_{stage}_{anchor_choice}'
    (ROOT/(PREFIX+suffix+f'_caps_drop{drop_budget}.json')).write_text(json.dumps(result, indent=2)+'\n')
    print({'accepted': len(caps), 'rejected': len(rejected)}, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=40)
    p.add_argument('--max-fixed', type=int, default=2)
    p.add_argument('--focus', type=int, default=0)
    p.add_argument('--drop-budget', type=int, default=0)
    p.add_argument('--configuration', default='45', choices=('45', '47', 'swapped47'))
    p.add_argument('--anchor', default='right', choices=('left', 'right'))
    a = p.parse_args()
    run(a.limit, a.max_fixed, a.focus, a.drop_budget, a.configuration, a.anchor)
