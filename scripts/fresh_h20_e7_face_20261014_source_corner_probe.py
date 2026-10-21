#!/usr/bin/env python3
"""Exact diagram probe of the published nineteen-wall f_2^13 chamber.

The tabulated source chamber has no A11+A2 corner for direct face transfer.
This neither excludes other chambers/lattices nor proves nonreflectivity.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_face_20261014_source_corner_probe'


def run():
    N = []
    for i in range(1, 13):
        r = [0]*14; r[i] = -1; r[i+1] = 1; N.append(r)
    N.append([0]*13+[-1])
    # McLeod2013 TableE.2, n=13, plus the B13 chamber at the basepoint.
    N.extend([[1, 2]+[0]*12, [1]+[1]*4+[0]*9, [2]+[1]*10+[0]*3,
              [3, 3]+[1]*11+[0], [3, 2, 2]+[1]*11, [5]+[2]*13])
    assert len(N) == 19 and all(len(r) == 14 for r in N)
    G = [[-2*r[0]*s[0]+sum(a*b for a, b in zip(r[1:], s[1:])) for s in N] for r in N]
    norms = [G[i][i] for i in range(19)]
    assert norms == [2]*12+[1]+[2]*4+[1, 2]
    counts = Counter()
    for i in range(19):
        for j in range(i):
            g = G[i][j]; assert g <= 0
            ratio = F(g*g, norms[i]*norms[j])
            if ratio >= 1: label = 'parallel' if ratio == 1 else 'dotted'
            else:
                assert ratio in [0, F(1, 4), F(1, 2), F(3, 4)]
                label = {F(0): 'right', F(1, 4): 'pi3', F(1, 2): 'pi4', F(3, 4): 'pi6'}[ratio]
            counts[label] += 1
    # Every node in A11+A2 has a pi3 neighbor. Integer pairings rule out
    # norm1 nodes, even next to a norm2 node (required squared pairing1/2).
    indices = [i for i, q in enumerate(norms) if q == 2]
    tested = 0; hits = []
    for J in combinations(indices, 13):
        tested += 1
        if any(G[i][j] not in [0, -1] for i, j in combinations(J, 2)): continue
        adj = {i: {j for j in J if G[i][j] == -1} for i in J}
        unseen = set(J); sizes = []; valid = True
        while unseen:
            todo = [min(unseen)]; C = set()
            while todo:
                i = todo.pop()
                if i in C: continue
                C.add(i); todo.extend(adj[i]-C)
            unseen -= C; degrees = [len(adj[i]) for i in C]
            if max(degrees) > 2 or sum(degrees) != 2*(len(C)-1): valid = False; break
            sizes.append(len(C))
        if valid and sorted(sizes) == [2, 11]: hits.append(list(J))
    assert tested == comb(17, 13) == 2380 and not hits
    script = f'scripts/{TAG}.py'
    return dict(status='PASS_TABULATED_F2_H13_CHAMBER_HAS_NO_A11_A2_CORNER', sources={script: sha256((ROOT/script).read_bytes()).hexdigest()},
        literature=dict(author='John Angus McLeod', title='Arithmetic Hyperbolic Reflection Groups', year=2013,
            url='https://etheses.durham.ac.uk/id/eprint/7743/1/thesis31072013.pdf', locations='Section3.1.4 printedp54; TableE.2 printedp149/PDFpage156 zero-based',
            transcription_scope='Six nonzero-height roots at n13 are transcribed from TableE.2. Reflectivity/completeness are literature inputs, not independently reproved by this diagram probe.'),
        diagonal_form=[-2]+[1]*13, normals=N, gram=G, norms=norms, pair_counts=dict(counts), thirteen_root_subsets=tested, matching_A11_A2_subsets=hits,
        scope='No direct A11+A2-corner matching to this particular published19wall fundamental chamber. Does not exclude another lattice, another reflection subgroup, a different chamber construction, rational isometry of forms, or completion of the current H13 face.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); a = ap.parse_args()
    d = run(); p = ROOT/f'results/{TAG}.json'
    if a.verify: assert d == json.loads(p.read_text())
    else: p.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['thirteen_root_subsets'])
