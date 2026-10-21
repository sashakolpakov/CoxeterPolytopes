#!/usr/bin/env python3
"""Two affine-E8 cusp completions of the frozen77-wall H20 chamber.

Exact Q(sqrt2)(sqrt79) arithmetic; a common integer shift is selected from
1,2,4,... by strict algebraic gates. No additional cap or frontier search.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_contact_closure_20261014_replacement_review as K
import fresh_h12_a12_source_20261014 as Qlin

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e8_cusp_completion_20261014'
BASE = 'results/fresh_h20_contact_closure_20261014_h12_replacement.json'
Z = K.ZERO


def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def vec(row): return [tuple(map(F, x)) for x in row]


def run():
    d = json.loads((ROOT/BASE).read_text())
    for p, h in d['sources'].items(): assert sha(p) == h, p
    Q = [[F(x) for x in row] for row in d['form']]
    old = list(map(vec, d['normals'])); assert len(old) == 77
    H = vec(d['strict_center']); escape = vec(d['strict_escape']['point'])
    cache = {}

    def dot(v, w):
        key = tuple(w)
        if key not in cache:
            cache[key] = [tuple(sum(Q[i][j]*w[j][k] for j in range(21)) for k in range(4)) for i in range(21)]
        out = Z
        for a, b in zip(v, cache[key]): out = K.add(out, K.mul(a, b))
        return out

    def linear(vectors, coeff):
        return [tuple(sum(c*v[j][k] for c, v in zip(coeff, vectors)) for k in range(4)) for j in range(21)]

    norms = list(map(F, d['norms']))
    assert all(dot(r, r) == K.scalar(q) and q > 0 for r, q in zip(old, norms))
    E = d['E8_support']; I = [[F(x) for x in row] for row in d['E8_inverse']]
    assert len(E) == 8 and E[0] == d['E8_minimum_weight_node'] == 12
    B = [[dot(old[i], old[j]) for j in E] for i in E]
    assert all(x[1:] == (0, 0, 0) for row in B for x in row)
    B = [[x[0] for x in row] for row in B]; Qlin.positive_definite(B)
    assert all(sum(I[i][k]*B[k][j] for k in range(8)) == int(i == j) for i in range(8) for j in range(8))
    assert [i for i in range(8) if I[i][i] == 2] == [0]
    assert all(x > 0 for row in I for x in row)
    minus_weight = linear([old[i] for i in E], [-x for x in I[0]])
    assert dot(minus_weight, minus_weight) == K.scalar(2)
    assert all(dot(minus_weight, old[j]) == K.scalar(-int(i == 0)) for i, j in enumerate(E))
    rays = [vec(row['point']) for row in d['H12_face_ideal_vertices']]
    assert len(rays) == 2 and K.sign(dot(*rays)) < 0
    old_incident = []
    for nu, row in zip(rays, d['H12_face_ideal_vertices']):
        assert dot(nu, nu) == Z and K.sign(dot(nu, H)) < 0
        assert dot(minus_weight, nu) == Z and minus_weight == vec(row['E8_link_recession'])
        vals = [dot(nu, r) for r in old]
        inc = [i for i, x in enumerate(vals) if x == Z]
        assert inc == row['incident'] and all(K.sign(x) <= 0 for x in vals)
        assert all(dot(minus_weight, old[i]) == K.scalar(-int(i == E[0])) for i in inc)
        old_incident.append(inc)

    def roots(t):
        return [[K.add(x, K.scale(y, t)) for x, y in zip(minus_weight, nu)] for nu in rays]

    trials = []; chosen = None
    for exponent in range(9):
        t = 2**exponent; caps = roots(t); failures = []
        for capindex, (r, inc) in enumerate(zip(caps, old_incident)):
            assert dot(r, r) == K.scalar(2)
            if K.sign(dot(r, H)) >= 0: failures.append(['center', capindex])
            for i, w in enumerate(old):
                g = dot(r, w)
                if i in inc:
                    assert g == K.scalar(-int(i == E[0]))
                elif not (K.sign(g) < 0 and K.sign(K.add(K.mul(g, g), K.scalar(-2*norms[i]))) > 0):
                    failures.append(['not_strictly_dotted', capindex, i])
        g = dot(*caps)
        if not (K.sign(g) < 0 and K.sign(K.add(K.mul(g, g), K.scalar(-4))) > 0):
            failures.append(['caps_not_strictly_dotted'])
        trials.append(dict(t=t, failures=failures))
        if not failures: chosen = t; break
    assert chosen == 8
    new = roots(chosen); N = old+new; norms += [F(2), F(2)]
    counts = Counter(); pair_rows = []
    for i in range(77, 79):
        for j in range(i):
            g = dot(N[i], N[j]); assert K.sign(g) <= 0
            if g == Z: kind = 'right'
            elif K.mul(g, g) == K.scalar(norms[i]*norms[j]/4): kind = 'pi3'
            else:
                assert K.sign(K.add(K.mul(g, g), K.scalar(-norms[i]*norms[j]))) > 0
                kind = 'dotted'
            counts[kind] += 1; pair_rows.append(dict(walls=[i, j], pairing=g, relation=kind))
    assert dict(counts) == {'right': 39, 'pi3': 2, 'dotted': 114}
    assert len(pair_rows) == 155
    assert K.sign(dot(H, H)) < 0 and all(K.sign(dot(H, r)) < 0 for r in N)
    # Verify actual facets directly, not merely retention of a halfspace list.
    facets = []
    for i, r in enumerate(N):
        coefficient = K.scale(dot(H, r), -1/norms[i])
        point = [K.add(a, K.mul(coefficient, b)) for a, b in zip(H, r)]
        assert K.sign(dot(point, point)) < 0 and K.sign(dot(point, H)) < 0
        assert all(dot(point, w) == Z if j == i else K.sign(dot(point, w)) < 0 for j, w in enumerate(N))
        facets.append(point)

    links = []
    for number, (nu, oldJ) in enumerate(zip(rays, old_incident)):
        vals = [dot(nu, r) for r in N]
        J = [i for i, x in enumerate(vals) if x == Z]
        assert J == oldJ+[77+number] and all(K.sign(x) <= 0 for x in vals)
        Gram = [[dot(N[i], N[j]) for j in J] for i in J]
        assert all(x[1:] == (0, 0, 0) for row in Gram for x in row)
        G = [[x[0] for x in row] for row in Gram]
        adj = {i: {j for j in J if i != j and dot(N[i], N[j]) != Z} for i in J}
        unseen = set(J); blocks = []; quotient_rank = 0
        while unseen:
            todo = [min(unseen)]; component = set()
            while todo:
                i = todo.pop()
                if i in component: continue
                component.add(i); todo.extend(adj[i]-component)
            unseen -= component; C = sorted(component)
            A = [[dot(N[i], N[j])[0] for j in C] for i in C]
            RR, piv = Qlin.rref(A); assert len(piv) == len(C)-1
            free = next(i for i in range(len(C)) if i not in piv)
            weights = [F(0)]*len(C); weights[free] = 1
            for row, j in zip(RR, piv): weights[j] = -row[free]
            if all(x < 0 for x in weights): weights = [-x for x in weights]
            assert all(x > 0 for x in weights)
            minor = [[A[i][j] for j in range(len(C)) if j != free] for i in range(len(C)) if i != free]
            pivots = Qlin.positive_definite(minor)
            total = linear([N[i] for i in C], weights)
            a, b = next((total[j][k], nu[j][k]) for j in range(21) for k in range(4) if nu[j][k])
            mu = a/b
            assert mu > 0 and total == [K.scale(x, mu) for x in nu]
            quotient_rank += len(C)-1
            blocks.append(dict(incident=C, gram=A, positive_kernel=weights,
                               positive_principal_pivots=pivots, null_multiple=mu))
        assert quotient_rank == 19
        e8block = next(row for row in blocks if 77+number in row['incident'])
        assert e8block['incident'] == sorted(E+[77+number]) and len(e8block['incident']) == 9
        links.append(dict(point=nu, incident=J, gram=G, affine_blocks=blocks,
                          euclidean_quotient_rank=quotient_rank, nonincident_strict=True))

    # The caps project to future null normals on the E8-perpendicular H12
    # face; its35 finite vertices are therefore unchanged (also checked).
    preserved = []
    for row in d['H12_face_finite_vertices']:
        p = vec(row['point']); values = [dot(p, r) for r in new]
        assert K.sign(dot(p, p)) < 0 and all(K.sign(x) < 0 for x in values)
        assert all(dot(p, old[i]) == Z for i in row['incident'])
        preserved.append(dict(incident=row['incident'], new_pairings=values))
    assert len(preserved) == 35
    escape_norm = dot(escape, escape)
    assert escape_norm == tuple(map(F, d['strict_escape']['norm']))
    assert K.sign(escape_norm) > 0 and K.sign(dot(escape, H)) < 0
    escape_pairs = [dot(escape, r) for r in N]
    assert all(K.sign(x) < 0 for x in escape_pairs)
    paths = [BASE, 'scripts/fresh_h20_contact_closure_20261014_h12_replacement.py',
             'scripts/fresh_h20_contact_closure_20261014_replacement_review.py',
             'scripts/fresh_h12_a12_source_20261014.py', f'scripts/{TAG}.py']
    return Qlin.encode(dict(status='PASS_79_WALL_TWO_E8_CUSP_COMPLETIONS_SAME_STRICT_ESCAPE',
        sources={p: sha(p) for p in paths}, field_basis=d['field_basis'], form=Q,
        retained_walls=77, normals=N, norms=norms, E8_support=E,
        E8_weight_node=E[0], negative_fundamental_weight=minus_weight,
        cusp_rays=rays, common_shift=chosen, doubling_trials=trials,
        new_pair_counts=dict(counts), new_pairs=pair_rows, strict_center=H,
        strict_facet_points=facets, completed_cusp_links=links,
        preserved_H12_finite_vertices=preserved,
        strict_escape=dict(point=escape, norm=escape_norm, pairings=escape_pairs,
                           new_pairings=escape_pairs[-2:], center_pairing=dot(escape, H)),
        conditional_nonquasi=d['conditional_nonquasi'],
        scope='Retains all77 actual facets and adds exactly2 norm2 caps -omega+8nu. All155 new pairs checked. The two selected ambient cusp links now have only positive affine blocks spanning quotientrank19;35 specified finite vertices and the full H12 face are retained. Same exact strict spacelike escape survives. No whole chamber cusp/vertex census or finite-volume claim; no further search.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    d = run(); path = ROOT/f'results/{TAG}.json'
    if args.verify: assert d == json.loads(path.read_text())
    else: path.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['new_pair_counts'], d['strict_escape']['norm'])
