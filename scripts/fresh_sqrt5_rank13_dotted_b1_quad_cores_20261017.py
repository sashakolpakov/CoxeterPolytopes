#!/usr/bin/env python3
"""Exact local arithmetic compact quadrilaterals from completed B1 dot stops.

Checks all singular cores with <=15 roots. Attachment-closure completeness
remains the separate producer replay; this script proves each local geometry.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_sqrt5_global_b1_dotted_fast_20261017 as P
import fresh_sqrt5_rank13_five_core_spectators_20261017 as S

A, C, K = S.A, S.C, P.K
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017'
INPUT = f'results/{P.TAG}.json'
LEMMA = 'results/fresh_singular_quad_local_realization_20261017.md'


def dump(x):
    return json.dumps(x, separators=(',', ':'))


def sha(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def run():
    prior = json.loads((ROOT/INPUT).read_text())
    assert prior['complete_queue'] and prior['marked_states'] == 25535
    assert prior['sources'] == P.manifest()
    selected = [r for r in prior['records'] if r['result']['singular'] and r['row']['nodes'] <= 15]
    assert len(selected) == 1399
    models = [(name, g) for name, g in C.models() if len(g) <= 13]
    finite = {C.code(g): name for name, g in models}
    component_proofs = [S.finite_component(name, g) for name, g in models]
    records, counts = [], Counter()
    for ri, source in enumerate(selected):
        row, parameters = source['row'], source['parameters']
        n = row['nodes']
        assert row['white'] == list(range(4, n))
        M = P.field_matrix(row, parameters)
        assert M is not None
        assert all(M[i][i].sign() > 0 and M[i][i].conj().sign() > 0 for i in range(n))
        assert all(M[i][j].sign() <= 0 for i, j in combinations(range(n), 2))
        cl = source['result']['classification']
        assert cl['physical'] == [n-2, 1, 1] and cl['conjugate'] == [n-1, 0, 1]
        kernel = [K(*x) for x in cl['kernel']]
        assert len(kernel) == n
        assert all(not sum((x*y for x, y in zip(r, kernel)), K()) for r in M)
        signs = [x.sign() for x in kernel]
        assert signs[0] == signs[2] and signs[1] == signs[3] and signs[0]*signs[1] == -1
        # The actual pair products are rebuilt from their exact K-valued gauges.
        edges = P.Q.edges(row)
        squares = {2: K(), 3: K(1), 4: K(2), 5: P.Q.PHI**2, 6: K(3), 10: 2+P.Q.PHI}
        dots = []
        for (i, j), label in edges.items():
            expected = K(*parameters[P.edge_key(i, j)]) if label == 0 else squares[label]
            assert 4*M[i][j]**2/(M[i][i]*M[j][j]) == expected
            if label == 0:
                assert (expected-4).sign() > 0
                dots.append([i, j])
        assert dots
        corners = []
        for pair in [(0, 1), (1, 2), (2, 3), (0, 3)]:
            graph = C.graph(row['white']+list(pair), edges)
            pieces = C.pieces(graph)
            # Recognition against exactly checked finite component models,
            # together with the exact normalized products above, proves PD.
            forest = tuple(sorted(finite[C.code(g)] for g in pieces))
            code = sorted(C.code(g, pair) for g in pieces)
            corners.append(dict(black=list(pair), forest=list(forest), marked_forest_code=code))
        # PD corner rank n-2, a dotted indefinite pair, and a nonzero radical
        # force physical inertia(n-2,1,1), independently of the saved label.
        # Deleting black3 gives a full basis since its kernel coefficient!=0.
        basis = [i for i in range(n) if i != 3]
        conjugate = [[A.K(M[i][j].a, -M[i][j].b) for j in basis] for i in basis]
        assert A.inertia(conjugate) == [n-1, 0, 0]
        cycles = P.D.cycles(M)
        assert all(r['integral'] for r in cycles)
        gram = [[x.enc() for x in r] for r in M]
        records.append(dict(index=ri, source_index=source['index'], nodes=n,
                            row=row, parameters=parameters, gram_sha256=sha256(dump(gram).encode()).hexdigest(),
                            physical_inertia=[n-2, 1, 1], conjugate_inertia=[n-1, 0, 1],
                            quotient_basis=basis, kernel=[x.enc() for x in kernel],
                            black_kernel_signs=signs[:4], dotted_pairs=dots,
                            integral_cycles=len(cycles), corners=corners, spectator_rank=15-n))
        counts['cores'] += 1
        counts['integral_cycles'] += len(cycles)
        counts['nodes_'+str(n)] += 1
        if ri % 50 == 0:
            print('LOCAL_DOTTED_QUAD', ri, source['index'], n, flush=True)
    paths = [f'scripts/{TAG}.py', INPUT, f'scripts/{P.TAG}.py',
             f'scripts/{P.D.TAG}.py', f'scripts/{S.TAG}.py', f'scripts/{A.TAG}.py',
             f'scripts/{C.TAG}.py', LEMMA]
    return dict(status='PASS_EXACT_RANK13_ELIGIBLE_DOTTED_B1_LOCAL_QUADRILATERALS',
                sources={p: sha(p) for p in paths}, counts=dict(counts),
                finite_component_certificates=component_proofs, records=records,
                proof='All exact field Grams are rebuilt. The four corner blocks match separately checked finite Coxeter components by their normalized products. A PD corner, a strictly dotted pair and the displayed nonzero radical force the physical inertia. A separately computed conjugate principal basis is positive definite. Every simple-cycle product is integral. The common-center, compact-four-ray face and stable integral reflection module follow from the explicit local-realization lemma.',
                scope='Every singular stopping state of the completed B1 dotted queue with at most15nodes is individually certified as a local admissible Qsqrt5 compact quadrilateral star. Additional orthogonal finite white spectators of rank15-n place it in ambientH13. No globally compact Coxeter chamber or dimension bound is claimed. This application does not independently reprove attachment closure or the producer hereditary payment pruning.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
        print('PASS_SAVED_DOTTED_B1_LOCAL_QUAD_REPLAY')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'], flush=True)
