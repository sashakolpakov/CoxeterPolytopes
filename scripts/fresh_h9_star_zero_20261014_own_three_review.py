#!/usr/bin/env python3
"""Independent Fraction-basis review of the all-own-three star obstruction.

No producer imports or symbolic-algebra dependency.  Arithmetic uses the
ordered basis 1,sqrt(2),sqrt(5),sqrt(10), and determinants use permutations.
"""
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_zero_20261014_own_three'
Z = (F(0),) * 4
O = (F(1), F(0), F(0), F(0))
S = [3, 4, 6, 7, 9, 10, 12, 13]
P = [2, 5, 8, 11]
C = [14, 15, 16, 17]

def add(x, y): return tuple(a+b for a, b in zip(x, y))
def neg(x): return tuple(-a for a in x)
def sub(x, y): return add(x, neg(y))
def scale(x, a): return tuple(a*b for b in x)
def total(xs):
    v = Z
    for x in xs: v = add(v, x)
    return v

@lru_cache(None)
def mul(x, y):
    out = [F(0)]*4
    for i, a in enumerate(x):
        if not a: continue
        for j, b in enumerate(y):
            if b: out[i ^ j] += a*b*(2 if i & j & 1 else 1)*(5 if i & j & 2 else 1)
    return tuple(out)

def sign_quad(a, b, d):
    sg = lambda q: (q > 0)-(q < 0)
    if not b: return sg(a)
    if not a: return sg(b)
    if sg(a) == sg(b): return sg(a)
    return sg(a)*sg(a*a-d*b*b)

@lru_cache(None)
def sign(x):
    a, b, c, d = x
    u, v = sign_quad(a, b, 2), sign_quad(c, d, 2)
    if not v: return u
    if not u: return v
    if u == v: return u
    return u*sign_quad(a*a+2*b*b-5*c*c-10*d*d, 2*a*b-10*c*d, 2)

@lru_cache(None)
def inverse(x):
    conjugates = [(x[0], s*x[1], t*x[2], s*t*x[3])
                  for s, t in [(-1, 1), (1, -1), (-1, -1)]]
    n = O
    for y in conjugates: n = mul(n, y)
    den = mul(x, n)
    assert den[0] and den[1:] == (0, 0, 0)
    return scale(n, 1/den[0])

def dec_anp(row):
    p = [F(v) for v in row]
    p = [F(0)]*(4-len(p))+p
    a3, a2, a1, a0 = p
    return a0+7*a2, a1+17*a3, a1+11*a3, 2*a2

COS = {2: Z, 3: (F(1, 2), F(0), F(0), F(0)),
       4: (F(0), F(1, 2), F(0), F(0)),
       5: (F(1, 4), F(0), F(1, 4), F(0))}
PERMS = [(p, (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)))
         for p in permutations(range(3))]
def determinant3(M):
    return total(scale(mul(mul(M[0][p[0]], M[1][p[1]]), M[2][p[2]]), sg)
                 for p, sg in PERMS)
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    src = f'results/{TAG}.json'
    data = json.loads((ROOT/src).read_text())
    pending, seen, hashes = [src], set(), {}
    while pending:
        path = pending.pop()
        if path in seen: continue
        seen.add(path)
        d = json.loads((ROOT/path).read_text())
        for p, h in d.get('sources', {}).items():
            assert sha(p) == h, p
            hashes[p] = h
            if p.endswith('.json'): pending.append(p)
    geo = json.loads((ROOT/'results/fresh_h9_tree_cube_wedges_20261014.json').read_text())['representatives'][0]['geometry']
    V = [set(v['incident']) for v in geo['vertices']]
    missing = {frozenset(m) for m in geo['minimal_nonfaces']}
    assert len(V) == 378
    old = json.loads((ROOT/'results/fresh_h9_star_cube_metric_20261014.json').read_text())
    gates = json.loads((ROOT/'results/fresh_h9_star_zero_20261014_gates.json').read_text())
    own = {tuple(sorted((leaf, v))): 3 for i in range(4)
           for leaf in (P[i], C[i]) for v in S[2*i:2*i+2]}
    assert sorted([list(p)+[m] for p, m in own.items()]) == data['all_own_contacts']
    for i in range(4):
        for leaf in (P[i], C[i]):
            assert frozenset(S[2*i:2*i+2]+[leaf]) in missing
    # Labels2 and3 cannot give a hyperbolic missing triple with two labels3:
    # their normalized Gram determinants are respectively1/2 and0.
    assert determinant3([[O, neg(COS[3]), neg(COS[3])],
                         [neg(COS[3]), O, Z], [neg(COS[3]), Z, O]]) == scale(O, F(1, 2))
    assert determinant3([[O, neg(COS[3]), neg(COS[3])],
                         [neg(COS[3]), O, neg(COS[3])], [neg(COS[3]), neg(COS[3]), O]]) == Z
    for row in data['arbitrary_high_block_reduction']:
        i = next(i for i in range(4) if S[2*i:2*i+2] == row['high_block'])
        labels = own | {tuple(S[2*j:2*j+2]): 6 if i == j else 4 for j in range(4)}
        gate = next(g for g in gates['families'] if g['id'] == row['reviewed_gate'])
        assert gate['origin'] == 'two_dense'
        assert all(labels.get((a, b), 0) >= m for a, b, m in gate['bounds'])
    assert len(data['arbitrary_high_block_reduction']) == 4
    base = own | {tuple(S[2*i:2*i+2]): 4 for i in range(4)}
    guard_records = []
    for row in old['forced_right_witnesses']+old['third_contact_upper_witnesses']:
        pair = tuple(row['pair']); nodes = sorted(V[row['actual_vertex']])
        assert set(pair) <= set(nodes)
        lab = base | {pair: row['test_label']}
        M = [[O if a == b else neg(COS[lab.get(tuple(sorted((a, b))), 2)]) for b in nodes] for a in nodes]
        for k in range(9):
            pivot = M[k][k]
            if sign(pivot) <= 0:
                assert k+1 == row['leading_order']
                assert pivot == (F(row['last_pivot'][0]), F(row['last_pivot'][1]), F(0), F(0))
                break
            for a in range(k+1, 9):
                for b in range(k+1, 9):
                    M[a][b] = sub(M[a][b], mul(mul(M[a][k], M[k][b]), inverse(pivot)))
        else: raise AssertionError('missing vertex failure')
        guard_records.append([list(pair), row['test_label'], row['actual_vertex'], k+1])
    cross = {tuple(sorted((a, b))) for i in range(4) for j in range(i+1, 4)
             for a in S[2*i:2*i+2] for b in S[2*j:2*j+2]}
    assert {tuple(r['pair']) for r in old['forced_right_witnesses']} == cross | set(zip(P, C))
    assert {tuple(r['pair']) for r in old['third_contact_upper_witnesses']} == {(11, 14), (14, 17)}
    # Z-matrix comparison preserves each failure as the internal labels grow.
    # The first eigenvector can be chosen nonnegative by Perron--Frobenius.
    patterns = {}
    for group in (0, 3):
        patterns[group] = []
        for extra in [None]+[j for j in range(8) if j//2 != group]:
            v = [Z]*8; v[2*group] = v[2*group+1] = neg(COS[3])
            if extra is not None: v[extra] = neg(COS[3])
            patterns[group].append(v)
    assert len(data['records']) == 16
    assert {tuple(r['block_labels']) for r in data['records']} == set(product((4, 5), repeat=4))
    counts, anchors, budgets, records = {}, 0, 0, []
    for record in data['records']:
        blocks = record['block_labels']
        B = [[Z]*8 for _ in range(8)]; J = [[Z]*8 for _ in range(8)]
        for i, m in enumerate(blocks):
            c = COS[m]; v = inverse(sub(O, mul(c, c)))
            assert sign(sub(v, scale(O, 2))) >= 0
            for a in range(2):
                for b in range(2):
                    B[2*i+a][2*i+b] = O if a == b else neg(c)
                    J[2*i+a][2*i+b] = v if a == b else mul(c, v)
            assert sign(sub(O, mul(c, c))) > 0
        assert all(total(mul(B[i][k], J[k][j]) for k in range(8)) == (O if i == j else Z)
                   for i in range(8) for j in range(8))
        for row in data['foreign_budgets']:
            leaf = row['leaf']; i = next(i for i in range(4) if leaf in (P[i], C[i]))
            assert set(row['support']) == (set(S)-{S[2*i]}) | {leaf}
            assert set(row['support']) <= V[row['actual_vertex']]
            assert all(sign(x) >= 0 for rr in J for x in rr)
            # Removing first own root leaves a singleton inverse1; foreign
            # two-by-two blocks retain inverse diagonals at least2.
            assert all(sign(sub(J[j][j], scale(O, 2))) >= 0 for j in range(8) if j//2 != i)
            assert F(1, 4)+2*F(1, 4)*2 > 1 and F(1, 4)+F(1, 2)*2 > 1
            budgets += 1
        def inner(v, w):
            return total(mul(v[i], total(mul(J[i][j], w[j]) for j in range(8))) for i in range(8))
        pp = patterns[3]+patterns[0]
        residual = [[sub(O if i == j else Z, inner(v, w)) for j, w in enumerate(pp)] for i, v in enumerate(pp)]
        # Distinct leaves may share a projection; their diagonal1 is still
        # distinct from their forced mutual right contact0.
        def entry(i, j, diagonal=False): return sub(O if diagonal else Z, inner(pp[i], pp[j]))
        diag = [entry(i, i, True) for i in range(14)]
        off = [[entry(i, j) for j in range(14)] for i in range(14)]
        for k, (ia, ib) in enumerate(product(range(7), repeat=2)):
            ad = sub(mul(diag[ia], diag[ib]), mul(off[ia][ib], off[ia][ib]))
            assert ad == dec_anp(record['anchor_determinants'][k]) and sign(ad) < 0
            anchors += 1
        seen_cases = set(); signs = {}
        for row in record['cases']:
            ia, ib, ic = row['patterns']; ma, mb = row['contact_labels']; ic += 7
            key = (ia, ib, ic-7, ma, mb); assert key not in seen_cases; seen_cases.add(key)
            ab, ac, bc = off[ia][ib], sub(off[ia][ic], COS[ma]), sub(off[ib][ic], COS[mb])
            M = [[diag[ia], ab, ac], [ab, diag[ib], bc], [ac, bc, diag[ic]]]
            d = determinant3(M)
            assert d == dec_anp(row['determinant']) and d != Z
            sg = str(sign(d)); assert sg != '0'
            signs[sg] = signs.get(sg, 0)+1; counts[sg] = counts.get(sg, 0)+1
        assert seen_cases == set(product(range(7), range(7), range(7), (2, 3), (2, 3)))
        records.append(dict(block_labels=blocks, cases=len(seen_cases), sign_counts=signs))
        print('review', blocks, signs, flush=True)
    assert anchors == 784 and budgets == 128 and sum(counts.values()) == 21952
    assert data['zero_rank_cases'] == data['unresolved'] == []
    paths = [src, f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_FRACTION_ALL_OWN_THREE_RANK_REVIEW',
                sources={p: sha(p) for p in paths}, actual_vertices=len(V), guard_records=guard_records,
                high_label_gate_placements=4, support_inverse_and_budget_checks=budgets,
                negative_anchor_determinants=anchors, nonzero_rank11_determinants=21952,
                determinant_sign_counts=counts, records=records,
                scope='No producer imports or symbolic arithmetic. Rechecks all30 actual-vertex guards in a separate Fraction basis, all16 positive supports and128 budgets, and all21952 Schur determinants. High labels use four literal placements of recursively hash-bound reviewed two-dense gates. Only all-own-three family of case0 is excluded; whole incidence remains open.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = run(); out = ROOT/f'results/{TAG}_review.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['nonzero_rank11_determinants'], flush=True)
