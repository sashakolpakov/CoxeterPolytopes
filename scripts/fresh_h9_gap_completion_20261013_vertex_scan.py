#!/usr/bin/env python3
"""Numerical proposals forcing a new finite vertex on an elliptic rank8 face."""
import json
from itertools import combinations
import numpy as np
import mpmath as mp
from scipy.linalg import null_space
from scripts.fresh_h9_gap_completion_20261013_base import *


def run():
    G, I, h, ps, source, case = configuration('48')
    names = list(ps)
    pairs = read(source)['cases'][case]['pair_checks']
    mp.mp.dps = 90
    C = np.array([[float(value(r['field'], t)) for t in r['c']] for r in ps.values()])
    Pn = np.array([[float(value(r['field'], t)) for t in r['p']] for r in ps.values()])
    Qn = np.array([float(value(K, r['q'])) for r in ps.values()])
    In = np.array([[float(value(K, t)) for t in row] for row in I])
    hn = np.array([float(value(K, t)) for t in h])
    gram = (C@Pn.T)/np.sqrt(Qn[:, None]*Qn[None, :])
    neighbors = [set() for _ in names]
    for i, j in combinations(range(len(names)), 2):
        r = pairs['|'.join(sorted((names[i], names[j])))]
        if r['kind'] == 'finite':
            neighbors[i].add(j)
            neighbors[j].add(i)
    blocks, visited = [], 0
    def walk(ids, eligible):
        nonlocal visited
        if len(ids) == 8:
            blocks.append(ids)
            return
        if len(ids)+len(eligible) < 8:
            return
        for i in sorted(eligible):
            if len(ids)+sum(j >= i for j in eligible) < 8:
                break
            trial = ids+(i,)
            visited += 1
            if np.linalg.eigvalsh(gram[np.ix_(trial, trial)])[0] <= 1e-8:
                continue
            walk(trial, {j for j in eligible & neighbors[i] if j > i})
    walk((), set(range(len(names))))
    print({'elliptic_rank8_blocks': len(blocks), 'tested_extensions': visited}, flush=True)
    grid = np.arange(-14, 14.0001, .04)
    proposals = []
    counts = {'blocks': len(blocks), 'dotted_only_solutions': 0, 'finite_attachment_solutions': 0}
    labels = (2, 3, 4, 5, 6, 10)
    allowed = np.array([np.cos(np.pi/m)**2 for m in labels])
    for block in blocks:
        N = null_space(C[list(block)])
        assert N.shape == (10, 2)
        eig, vec = np.linalg.eigh(N.T@In@N)
        assert eig[0] < 0 < eig[1]
        timelike, spatial = N@vec[:, 0]/np.sqrt(-eig[0]), N@vec[:, 1]/np.sqrt(eig[1])
        # Every non-block wall must be ultraparallel in this first route.
        pp = np.sqrt(2)*(np.cosh(grid)[:, None]*spatial + np.sinh(grid)[:, None]*timelike)
        pp = np.r_[pp, -pp]
        normalized = (pp@C.T)/np.sqrt(2*Qn)
        outside = [i for i in range(len(names)) if i not in block]
        ok = (pp@In@hn > 1e-8) & np.all(normalized[:, outside] < -1-1e-7, axis=1)
        if np.any(ok):
            j = min(np.flatnonzero(ok), key=lambda j: np.linalg.norm(pp[j]))
            proposals.append({'kind': 'eight_right_others_dotted', 'block': list(block),
                              'block_names': [names[i] for i in block],
                              'seed_pairings_numeric': pp[j].tolist()})
            counts['dotted_only_solutions'] += 1
        # A finite attachment adds one exact scalar equation to the Lorentz plane.
        local_seen = set()
        for anchor in outside:
            aa = C[anchor]@N
            if np.linalg.norm(aa) < 1e-8:
                continue
            direction = np.array([-aa[1], aa[0]])
            normform = N.T@In@N
            A = direction@normform@direction
            if abs(A) < 1e-10:
                continue
            for label in labels:
                target = -np.sqrt(2*Qn[anchor])*np.cos(np.pi/label)
                origin = aa*target/(aa@aa)
                B = 2*origin@normform@direction
                D = B*B-4*A*(origin@normform@origin-2)
                if D < 1e-10:
                    continue
                for branch in (-1, 1):
                    x = (-B+branch*np.sqrt(D))/(2*A)
                    p = N@(origin+x*direction)
                    if p@In@hn <= 1e-7:
                        continue
                    pair = C@p
                    ratio = pair*pair/(2*Qn)
                    good = (pair <= 1e-7) & ((ratio > 1+1e-7) |
                        (np.min(abs(ratio[:, None]-allowed), axis=1) < 1e-7))
                    if not np.all(good):
                        continue
                    key = tuple(np.round(p, 8))
                    if key in local_seen:
                        continue
                    local_seen.add(key)
                    proposals.append({'kind': 'eight_right_finite_anchor', 'block': list(block),
                        'block_names': [names[i] for i in block], 'anchor_name': names[anchor],
                        'anchor_label': label, 'seed_pairings_numeric': p.tolist()})
                    counts['finite_attachment_solutions'] += 1
    data = {'source': source, 'source_case': case, 'wall_names': names, 'counts': counts,
            'numeric_elliptic_blocks': [[names[i] for i in block] for block in blocks],
            'proposals': proposals,
            'scope': 'Numerical rank-eight positive-definite blocks and compatible finite-vertex '
                     'cap proposals. No exact positive conclusion or exhaustive exclusion.'}
    (ROOT/(PREFIX+'_vertex_proposals.json')).write_text(json.dumps(data, indent=2)+'\n')
    print({'counts': counts, 'proposals': len(proposals)}, flush=True)


if __name__ == '__main__':
    run()
