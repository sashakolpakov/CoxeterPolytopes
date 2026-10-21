#!/usr/bin/env python3
"""Exact source and Gram-table obligations for the all-real H4 obstruction.

The geometric/arithmetic proof is in the matching report. This program
verifies its finite arithmetic premises, not the general theorems it cites.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import hashlib
import json
from fresh_h20_h4_edge_sections_20261013 import K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_h4_face_closure_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'


def sha(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def inertia(A):
    A = [[F(x) for x in row] for row in A]
    signs = []
    while A:
        if not any(A[i][i] for i in range(len(A))):
            pair = next(((i, j) for i in range(len(A)) for j in range(i)
                         if A[i][j]), None)
            if pair is None:
                signs.extend([0]*len(A))
                break
            i, j = pair
            # Congruence replacing basis vector i by i+j.
            for k in range(len(A)):
                A[k][i] += A[k][j]
            for k in range(len(A)):
                A[i][k] += A[j][k]
        i = next(i for i in range(len(A)) if A[i][i])
        A[0], A[i] = A[i], A[0]
        for row in A:
            row[0], row[i] = row[i], row[0]
        pivot = A[0][0]
        signs.append(1 if pivot > 0 else -1)
        A = [[A[i][j] - A[i][0]*A[0][j]/pivot
              for j in range(1, len(A))] for i in range(1, len(A))]
    return dict(Counter(map(str, signs)))


def components(E, G):
    todo = set(E)
    out = []
    while todo:
        cc = [min(todo)]
        todo.remove(cc[0])
        for i in cc:
            more = [j for j in sorted(todo) if G[i][j]]
            cc.extend(more)
            todo.difference_update(more)
        out.append(cc)
    return out


def run():
    raw = json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q = [[-int(x) for x in row] for row in raw['basis_gram']]
    W = [[int(x) for x in row] for row in raw['wall_coordinates_as_columns_transposed']]
    C = [[sum(Q[i][j]*w[j] for j in range(22)) for i in range(22)] for w in W]
    G = [[sum(x*y for x, y in zip(w, c)) for c in C] for w in W]
    assert len(W) == 210 and all(len(w) == 22 for w in W)
    signature = inertia(Q)
    assert signature == {'1': 21, '-1': 1}
    assert all(F(2*x, G[i][i]).denominator == 1
               for i, c in enumerate(C) for x in c)
    values = defaultdict(set)
    for i in range(210):
        for j in range(i):
            values[tuple(sorted([G[i][i], G[j][j]]))].add(G[i][j])
    assert dict(values) == {(2, 2): {-1, 0}, (2, 4): {-2, 0},
                            (4, 4): {-8, -6, -4, -2}}
    qe = K(F(1, 4), F(1, 4))
    assert (qe-F(4, 5)).sign() > 0
    short_modified = K(2) - qe
    h4_cosine = 1/(2*(K(2)-1/qe))
    assert h4_cosine == K(F(3, 8), F(1, 8))
    # Complete squared-Coxeter-cosine list lying in Q(sqrt5).
    allowed = {2: K(), 3: K(F(1, 4)), 4: K(F(1, 2)),
               5: K(F(3, 8), F(1, 8)), 6: K(F(3, 4)),
               10: K(F(5, 8), F(1, 8))}
    types = [('rational_short', 2, 0), ('rational_long', 4, 0),
             ('modified_short', 2, -1), ('modified_long', 4, -2)]
    table = []
    for (name1, d1, k1), (name2, d2, k2) in combinations_with_replacement(types, 2):
        if k1 == k2 == 0:
            continue
        n1, n2 = K(d1)-qe*k1*k1, K(d2)-qe*k2*k2
        assert n1.sign() > 0 and n2.sign() > 0
        for old in sorted(values[tuple(sorted([d1, d2]))]):
            pair = K(old) - qe*k1*k2
            determinant = n1*n2 - pair*pair
            row = {'types': [name1, name2], 'old_pairing': old,
                   'projected_norms': [n1.enc(), n2.enc()],
                   'projected_pairing': pair.enc(), 'rank2_determinant': determinant.enc()}
            if not pair:
                assert (k1 == 0) != (k2 == 0) and old == 0
                row['verdict'] = 'ORTHOGONAL_TO_UNMODIFIED_RATIONAL_ROOT'
            elif determinant.sign() <= 0:
                row['verdict'] = 'NO_FINITE_INTERSECTION'
            else:
                assert pair.sign() < 0
                cosine = pair*pair/(n1*n2)
                assert all(cosine != v for v in allowed.values())
                row['verdict'] = 'FORBIDDEN_COXETER_ANGLE'
                row['squared_cosine'] = cosine.enc()
            table.append(row)
    assert Counter(row['verdict'] for row in table)['FORBIDDEN_COXETER_ANGLE'] == 3
    walk = json.loads((ROOT/WALK).read_text())
    assert walk['status'] == 'PASS_COMPLETE_CORNER_ORBIT_CENSUS'
    assert not walk['pending_orbit_ids'] and walk['sources'][RAW] == sha(RAW)
    edges = set()
    count = 0
    a3_count = 0
    for corner in walk['corner_orbits']:
        for edge in corner['edge_orbits']:
            E = edge['edge_walls']
            if tuple(E) in edges:
                continue
            edges.add(tuple(E))
            for cc in components(E, G):
                if (len(cc) == 3 and all(G[i][i] == 2 for i in cc)
                    and sorted(sum(G[i][j] != 0 for j in cc)-1 for i in cc) == [1, 1, 2]
                    and all(G[i][j] in [0, -1] for i in cc for j in cc if i != j)):
                    a3_count += 1
                if (len(cc) == 4 and all(G[i][i] == 2 for i in cc)
                    and sorted(sum(G[i][j] != 0 for j in cc)-1 for i in cc) == [1, 1, 2, 2]
                    and all(G[i][j] in [0, -1] for i in cc for j in cc if i != j)):
                    count += 1
    assert len(edges) == 6492 and count == 1660
    # A connected spherical extension of the H3 chain is H4, attaching
    # at the endpoint opposite its label5 edge. Orthogonality to the
    # anchor forces (e,r_l)=q*k. The required new label3 then gives
    # d-q*k² = 2*G_lc². Irrationality forces k=0; the source pools
    # allow only a norm2 root paired -1 to that far endpoint.
    h3_extensions = []
    for d, anchor_pool, far_pool in [(2, [0, -1], [0, -1]),
                                     (4, [0, -2], [0, -2])]:
        for k in anchor_pool:
            for g in far_pool:
                if g < 0 and K(d)-qe*k*k == K(2*g*g):
                    h3_extensions.append([d, k, g])
    assert h3_extensions == [[2, 0, -1]]
    return {'status': 'PASS_EXACT_PREMISES_OF_ALL_REAL_A4_TO_H4_OBSTRUCTION',
            'sources': {p: sha(p) for p in [RAW, WALK, 'scripts/'+Path(__file__).name,
                        'scripts/fresh_h20_h4_edge_sections_20261013.py']},
            'source_signature': signature, 'integral_root_reflections': 210,
            'source_pairing_pools': {str(k): sorted(v) for k, v in values.items()},
            'normal_squared_norm': qe.enc(), 'H4_squared_cosine': h4_cosine.enc(),
            'complete_Qsqrt5_squared_Coxeter_cosines': {str(m): v.enc() for m, v in allowed.items()},
            'modified_facet_adjacency_table': table,
            'table_counts': dict(Counter(row['verdict'] for row in table)),
            'literal_edge_supports': len(edges), 'A4_components': count,
            'continuous_one_contact_families': 2*count,
            'A3_components': a3_count,
            'continuous_A3_to_H3_families': 2*a3_count,
            'H3_connected_extension_solutions_norm_anchor_far': h3_extensions,
            'scope': 'Exact finite arithmetic premises. The matching report proves that an actual one-contact A4-to-H4 section face forces a compact rational arithmetic H16 or H15 chamber. The A3-to-H3 extension forces that forbidden H4 face or a compact rational arithmetic H17 or H16 chamber. All real parameters in the3320 A4 and4230 A3 saved source-edge families are excluded; no floating-point census is used in that proof.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == result
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['table_counts'])
