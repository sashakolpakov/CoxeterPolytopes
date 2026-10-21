#!/usr/bin/env python3
"""Transport existing mirrors by the exact seed6/seed7 symmetry."""
import json
import numpy as np
import mpmath as mp
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate

STAGE = PREFIX+'_geometry.json'


def stage():
    G, I, h, ps = current()
    data = read(STAGE)
    bank = read(data['source'])
    case = data['cases'][0]
    for i in case['selected_caps']:
        ps[f'gap{i}'] = cap_profile(bank['caps'][i], f'gap{i}')
    assert list(ps) == case['walls']
    return G, I, h, ps


def run():
    G, I, h, ps = stage()
    perm = list(range(10))
    perm[6], perm[7] = 7, 6
    assert all(G[i][j] == G[perm[i]][perm[j]] for i in range(10) for j in range(10))
    assert all(h[i] == h[perm[i]] for i in range(10))
    _, _, _, originals = profiles()
    candidates = originals|ps
    targets = [tuple(map(parse, r['certificate']['y']))
               for r in read(STAGE)['cases'][0]['geometry']['escapes']]
    pool = angle_pool(ANGLE_SOURCE)
    mp.mp.dps = 90
    names = list(ps)
    cp = np.array([[float(value(r['field'], t)) for t in r['p']] for r in ps.values()])
    q = np.array([float(value(K, r['q'])) for r in ps.values()])
    allowed = np.array([np.cos(np.pi/m)**2 for m in (2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30)])
    proposals = []
    counts = {'transported': 0, 'any_target_cut': 0, 'all_compatible': 0}
    for source, original in candidates.items():
        f = original['field']
        if all(f.sign(f.sub(original['c'][i], original['c'][perm[i]])) == 0 for i in range(10)):
            continue
        row = dict(original, name='swap_'+source,
                   c=tuple(original['c'][i] for i in perm),
                   p=tuple(original['p'][i] for i in perm))
        counts['transported'] += 1
        signs = [f.sign(dot(f, row['c'], tuple(lift_K(f, t) for t in y))) for y in targets]
        if not any(s < 0 for s in signs):
            continue
        counts['any_target_cut'] += 1
        c = np.array([float(value(f, t)) for t in row['c']])
        pair = cp@c
        ratio = pair*pair/(float(value(K, row['q']))*q)
        good = (pair <= 1e-7) & ((ratio > 1+1e-7) |
                               (np.min(abs(ratio[:, None]-allowed), axis=1) < 1e-7))
        drops = [names[i] for i in np.flatnonzero(~good)]
        if len(drops) > 2 or any(n.startswith('seed') for n in drops):
            continue
        checks = {n: relation(row, r, pool) for n, r in ps.items()}
        drops = [n for n, r in checks.items() if r['kind'] not in ('finite', 'dotted')]
        if len(drops) > 2 or any(n.startswith('seed') for n in drops):
            continue
        counts['all_compatible'] += not drops
        proposals.append((source, row, drops, checks, signs))
        print({'source': source, 'dropped': drops, 'cuts': [i for i, s in enumerate(signs) if s < 0]}, flush=True)
    cases = []
    for source, row, drops, checks, signs in proposals:
        kept = [r for n, r in ps.items() if n not in drops]
        selected = kept+[row]
        ns = [r['name'] for r in selected]
        pairs = dict(read(STAGE)['cases'][0]['pair_checks'])
        pairs.update({'|'.join(sorted((row['name'], n))): r for n, r in checks.items()})
        pairs = {k: r for k, r in pairs.items() if all(n in ns for n in k.split('|'))}
        census = affine_census(ns, pairs)
        case = {'transported_source': source, 'permutation': perm, 'dropped': drops,
                'walls': ns, 'pair_checks': pairs, 'target_signs': signs, 'affine_census': census}
        if census['affine_circuits']:
            case['status'] = 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
        else:
            case['geometry'] = volume_gate(I, h, selected)
            case['status'] = case['geometry']['status']
        cases.append(case)
        print({'case': len(cases)-1, 'source': source, 'status': case['status'],
               'maximum_numeric_vertex_norm': case.get('geometry', {}).get('maximum_numeric_vertex_norm')}, flush=True)
    data = {'source': STAGE, 'source_case': 0, 'permutation': perm, 'counts': counts,
            'cases': cases, 'scope': 'Existing source normals transported by a seed isometry; '
            'at most two declared cap replacements. Numeric screening is not a complete '
            'finite-angle exclusion. Actual accepted geometry has exact certificates.'}
    (ROOT/(PREFIX+'_swap.json')).write_text(json.dumps(data, indent=2)+'\n')
    print(counts, flush=True)


if __name__ == '__main__':
    run()
