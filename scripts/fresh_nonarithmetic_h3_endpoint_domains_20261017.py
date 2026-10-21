#!/usr/bin/env python3
"""Complete continuous endpoint families for all open edges of the S6 face.

This is a domain certificate, not a claim that the remaining transcendental
Coxeter contact equations or coordinated compact closure have been solved.
"""
import argparse
import itertools
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import sympy as sy
import fresh_nonarithmetic_a2_escape_partition_20261017 as P

U, B, I, C53 = P.U, P.B, P.I, P.C53
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_h3_endpoint_domains_20261017'
S = P.S
THETA = {2: sy.Integer(0), 3: sy.Integer(1), 4: sy.sqrt(2),
         5: (1 + sy.sqrt(5)) / 2, 6: sy.sqrt(3),
         10: sy.sqrt((5 + sy.sqrt(5)) / 2)}


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def labels():
    old = B.read('results/fresh_h9_persistent_cusp_screen_20261016.json')
    matrix = [r + [0] * 6 for r in old['pair_matrix']] + [[0] * 55 for _ in range(6)]
    rows = B.read(P.C53.CL)['caps'] + [dict(index=52, relations=B.read(P.CAP52)['relations'])]
    rows += B.read(P.CAPS55)['caps']
    for cap in rows:
        i = cap['index']
        matrix[i][i] = 2
        for relation in cap['relations']:
            j = relation['wall']
            matrix[i][j] = matrix[j][i] = relation.get('label', 0)
    return matrix


def classify_support(ids, matrix):
    if any(matrix[i][j] == 0 for i, j in itertools.combinations(ids, 2)):
        return None, ['contains_dotted_pair']
    gram = sy.Matrix([[2 if i == j else -THETA[matrix[i][j]] for j in ids] for i in ids])
    minors = [sy.simplify(gram[:k, :k].det()) for k in range(1, 5)]
    if not all(x.is_positive is True for x in minors):
        assert any(x.is_nonpositive is True for x in minors)
        return None, [str(x) for x in minors]
    return gram, [str(x) for x in minors]


def profiles(gram):
    inverse = sy.simplify(gram.inv())
    assert all(x.is_nonnegative is True for x in inverse)
    isolated = [i for i in range(4) if all(i == j or gram[i, j] == 0 for j in range(4))]
    for i in range(4):
        if i in isolated:
            assert inverse[i, i] == sy.Rational(1, 2)
        else:
            assert sy.simplify(inverse[i, i] - sy.Rational(2, 3)).is_nonnegative is True
    free = [i for i in range(4) if sy.simplify(inverse[i, i] - 2).is_negative is True]
    assert all(i in free or sy.simplify(inverse[i, i] - 2).is_nonnegative is True
               for i in range(4))
    accepted, rejected, null = [], 0, []
    for assignment in itertools.product(range(2, 6), repeat=len(free)):
        ls = [2] * 4
        for i, m in zip(free, assignment):
            ls[i] = m
        c = sy.Matrix([THETA[m] for m in ls])
        cost = sy.expand((c.T * inverse * c)[0])
        residual = sy.simplify(2 - cost)
        if residual == 0:
            null.append(ls)
        elif residual.is_positive is True:
            accepted.append(dict(labels=ls, projection_cost=str(cost), residual_square=str(residual)))
        else:
            assert residual.is_negative is True
            rejected += 1
    return dict(inverse=[[str(sy.simplify(inverse[i, j])) for j in range(4)] for i in range(4)],
                inverse_diagonal=[str(sy.simplify(inverse[i, i])) for i in range(4)],
                isolated_A1_slots=isolated, free_slots=free,
                finite_label_profiles=accepted, zero_residual_profiles=null,
                positive_residual_profiles=len(accepted), negative_residual_count=rejected,
                tested_vectors=4 ** len(free),
                tails=[dict(slot=i, labels='m>=6 at this slot; all others2',
                            magnitude='c=2cos(pi/m)', continuous_relaxation='sqrt(3)<=c<2',
                            projection_cost='c^2/2', residual_square='2-c^2/2') for i in isolated])


def run():
    _, qi, normals, h, _, _, _, cone = P.physical_data()
    pair = lambda x, y: U.idot(x, U.imv(qi, y))
    matrix = labels()
    rays = [(kind, tuple(ids)) for kind, groups in cone['complete_census'].items()
            for ids in groups if set(S) <= set(ids)]
    assert Counter(kind for kind, _ in rays) == dict(finite_vertices=13, hyperideal_rays=7)
    assert all(len(ids) == 9 for _, ids in rays)
    # Every raw one-face is recognized by its complete common-zero ray set.
    raw = sorted({tuple(sorted(set(a) & set(b))) for (_, a), (_, b)
                  in itertools.combinations(rays, 2) if len(set(a) & set(b)) == 8})
    assert len(raw) == 30
    vertex_ids = S + [7, 8, 11]
    assert vertex_ids in cone['complete_census']['finite_vertices']
    vg = [[pair(normals[i], normals[j]) for j in vertex_ids] for i in vertex_ids]
    P.pd(vg)
    vi, _ = C53.inverse(vg)
    v = P.sub(h, P.combine(U.imv(vi, [pair(h, normals[i]) for i in vertex_ids]),
                           [normals[i] for i in vertex_ids]))
    assert pair(v, v).h < 0 and pair(v, h).h < 0
    records, inventory = [], []
    for active in raw:
        endpoints = [(kind, ids) for kind, ids in rays if set(active) <= set(ids)]
        assert len(endpoints) == 2
        sides = [i for i in active if i not in S]
        ids = [5, 6] + sides
        assert all(matrix[i][j] == 2 for i in S[:4] for j in ids)
        gram4, minors = classify_support(ids, matrix)
        item = dict(active=list(active), sides=sides,
                    endpoint_kinds=[kind for kind, _ in endpoints],
                    endpoint_supports=[list(ss) for _, ss in endpoints], residual_leading_minors=minors,
                    geodesic=gram4 is not None)
        if gram4 is None:
            inventory.append(item)
            continue
        open_ends = sum(kind == 'hyperideal_rays' for kind, _ in endpoints)
        item['open_ends'] = open_ends
        inventory.append(item)
        if not open_ends:
            continue
        prof = profiles(gram4)
        full = S[:4] + ids
        g = [[pair(normals[i], normals[j]) for j in full] for i in full]
        pivots = P.pd(g)
        gi, inverse_pivots = C53.inverse(g)

        def project(x):
            return P.sub(x, P.combine(U.imv(gi, [pair(x, normals[i]) for i in full]),
                                     [normals[i] for i in full]))

        hp = project(h)
        hp2 = pair(hp, hp)
        assert hp2.h < 0
        hscale = U.sqrt(-hp2)
        k0 = [U.rnd(x / hscale) for x in hp]
        assert pair(k0, h).h < 0
        # An exact coordinate vector, projected and made perpendicular to k0,
        # supplies the positive direction of the two-dimensional kernel.
        for coordinate in range(10):
            e = [I(int(i == coordinate)) for i in range(10)]
            ep = project(e)
            alpha = pair(ep, k0)
            w = [U.rnd(x + alpha * y) for x, y in zip(ep, k0)]
            w2 = pair(w, w)
            if w2.l > 0:
                break
        else:
            raise AssertionError('No positive kernel direction')
        ws = U.sqrt(w2)
        k1 = [U.rnd(x / ws) for x in w]
        # Duals d_i satisfy <d_i,n_j>=delta_ij in the last four slots.
        # With negative contact vector -c, p=-sum c_i d_i exactly.
        duals = [P.combine([gi[j][4 + i] for j in range(8)],
                           [normals[j] for j in full]) for i in range(4)]
        response = [[pair(d, n) for d in duals] for n in normals]
        b = [pair(k1, n) for n in normals]
        d = [pair(k0, n) for n in normals]
        records.append(dict(sides=sides, full_support=full, residual_nodes=ids,
                            open_ends=open_ends, profiles=prof,
                            positive_gram_pivots=pivots, inverse_pivots=inverse_pivots,
                            full_gram=[P.enc(row) for row in g],
                            kernel_future_unit=P.enc(k0), kernel_positive_unit=P.enc(k1),
                            kernel_coordinate_seed=coordinate,
                            unnormalized_future_square=hp2.enc(),
                            unnormalized_positive_square=w2.enc(),
                            native_response_matrix=[P.enc(row) for row in response],
                            native_kernel_positive_pairings=P.enc(b),
                            native_kernel_future_pairings=P.enc(d),
                            common_interior_vertex_response=P.enc([pair(x, v) for x in duals]),
                            common_interior_vertex_positive=pair(k1, v).enc(),
                            common_interior_vertex_future=pair(k0, v).enc()))
    assert sum(row['geodesic'] for row in inventory) == 28
    assert len(records) == 13 and sum(row['open_ends'] for row in records) == 17
    finite_count = sum(row['profiles']['positive_residual_profiles'] for row in records)
    tail_count = sum(len(row['profiles']['tails']) for row in records)
    assert finite_count == 296 and tail_count == 14
    paths = [f'scripts/{TAG}.py', f'scripts/{P.TAG}.py', f'results/{P.TAG}.json',
             P.U.CAP, P.C53.CL, P.CAP52, P.CAPS55, P.CONE55,
             'results/fresh_h9_persistent_cusp_screen_20261016.json',
             'results/fresh_nonarithmetic_nonright_s6_20261016_cone_review.json']
    for path in paths:
        if path.endswith('.json'):
            B.read(path)
    return dict(status='PASS_COMPLETE_ALL_LABEL_CONTINUOUS_OPEN_EDGE_ENDPOINT_DOMAINS',
                sources={path: sha(path) for path in paths}, retained_walls=55,
                face_support=S, face_ray_counts=dict(Counter(kind for kind, _ in rays)),
                raw_one_faces=30, geodesic_edges=28, open_edges=13, open_edge_ends=17,
                edge_inventory=inventory, finite_profiles=finite_count,
                continuous_high_label_tails=tail_count, endpoint_domains=records,
                common_interior_vertex_support=vertex_ids, common_interior_vertex=P.enc(v),
                exact_normal_formula='For c=(2cos(pi/m_i)) with c_i=0 for label2, let p=-sum_i c_i*d_i, rho=2-c^T J c>0. Every root with these eight prescribed contacts and norm2 is r=p+s*sqrt(rho)*(cosh(u)*k1+sinh(u)*k0), s in {-1,1}, u real. For high-label tails only one isolated-A1 coordinate c is nonzero and sqrt3<=c<2; actual labels are integers m>=6.',
                all55_native_gate_formula='For every j=0..54, <r,n_j>=-response[j]*c+s*sqrt(rho)*(b[j]*cosh(u)+d[j]*sinh(u)). Require this to belong to {0,-2cos(pi/m):integer m>=3} union(-infinity,-2]. No unspecified finite label is silently rejected. The same formula against the saved finite vertex must be strictly negative.',
                completeness='The H4 block forces all four contacts right. In each residual rank4 support every nonisolated inverse diagonal is at least2/3, so a label>=6 cannot occur there. A label>=6 at isolatedA1 costs at least3/2; any second nonzero contact costs at least1/2, so all others must be zero. Thus label2..5 vectors plus isolatedA1 tails are complete. Positive Schur residual is necessary and sufficient for a finite endpoint Gram.',
                interior_equivalence='Under all55 acuteness inequalities, pairing strictly negatively with the saved old finite vertex is equivalent to admitting some common strict future interior: necessity follows finite-vertex persistence, sufficiency by perturbing that vertex toward the retained strict center. The retained chosen center itself is not imposed on new walls.',
                scope='Complete unavoidable endpoint normal families for the13 open edges, including all high-label tails. This certificate provides exact native offsets and parameter formulas; it does not solve their remaining all55 Coxeter gates, mutual new-wall compatibility, or coordinated compact closure. No finite-volume or no-completion claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
        print('PASS saved all-label endpoint-domain replay')
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['finite_profiles'], data['continuous_high_label_tails'])
