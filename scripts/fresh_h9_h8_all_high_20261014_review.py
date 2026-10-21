#!/usr/bin/env python3
"""Separate algebraic-field and 192-bit proof of the all-high obstruction."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
from fresh_h9_h8_facet_cover_20261014_line_review import I, BITS

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_all_high_20261014_review'
PROOF = 'results/fresh_h9_h8_all_high_20261014.json'
BASE = 'results/fresh_h9_mixed_d4_coupling_20261014.json'
SOURCE = 'data/known_examples/bugaenko_compact_h8.json'
SQ5 = s.sqrt(5); RR = s.sqrt((16+4*SQ5)/11)
K = s.QQ.algebraic_field(SQ5, RR); Z, O = K.zero, K.one

def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def decode(x):
    return K.from_sympy(s.Rational(x[0][0])+s.Rational(x[0][1])*SQ5+
                       (s.Rational(x[1][0])+s.Rational(x[1][1])*SQ5)*RR)
def source_decode(x):
    return K.from_sympy(s.Rational(x['a'])+s.Rational(x['b'])*(1+SQ5)/2)
def ev(expr):
    if not isinstance(expr, s.Basic): expr = K.to_sympy(expr)
    if expr.is_Rational: return I(F(int(expr.p), int(expr.q)))
    if expr.is_Add: return sum((ev(x) for x in expr.args), I(0))
    if expr.is_Mul:
        out = I(1)
        for x in expr.args: out *= ev(x)
        return out
    if expr.is_Pow:
        base, power = expr.args
        if power == s.Rational(1, 2): return ev(base).sqrt()
        if power == -1: return 1/ev(base)
        if power == s.Rational(-1, 2): return 1/ev(base).sqrt()
        if power.is_Integer and power > 0:
            out = I(1)
            for _ in range(int(power)): out *= ev(base)
            return out
    raise AssertionError(expr)
def dot(a, b): return sum((x*y for x, y in zip(a, b)), Z)
def mv(M, v): return [dot(row, v) for row in M]
def pair(v, w, Q): return dot(v, mv(Q, w))
def encode(M):
    if isinstance(M, (list, tuple)): return [encode(x) for x in M]
    return str(K.to_sympy(M))
def eliminate(M):
    M = [row[:] for row in M]; pivots = []
    for j in range(len(M[0])):
        k = len(pivots)
        p = next((i for i in range(k, len(M)) if M[i][j]), None)
        if p is None: continue
        M[k], M[p] = M[p], M[k]
        v = M[k][j]; M[k] = [x/v for x in M[k]]
        for i in range(len(M)):
            if i != k:
                v = M[i][j]; M[i] = [x-v*y for x, y in zip(M[i], M[k])]
        pivots.append(j)
        if len(pivots) == len(M): break
    return M, pivots
def inverse(M):
    n = len(M)
    R, pivots = eliminate([row+[O if i == j else Z for j in range(n)] for i, row in enumerate(M)])
    assert pivots == list(range(n))
    return [row[n:] for row in R]
def null_one(M):
    R, pivots = eliminate(M); n = len(M[0]); assert len(pivots) == n-1
    free = next(i for i in range(n) if i not in pivots)
    v = [Z]*n; v[free] = O
    for i, p in enumerate(pivots): v[p] = -R[i][free]
    assert all(not x for x in mv(M, v))
    return v
def ldl(M):
    M = [row[:] for row in M]; pivots = []
    for k in range(len(M)):
        p = M[k][k]; assert p; pivots.append(p)
        for i in range(k+1, len(M)):
            for j in range(k+1, len(M)): M[i][j] -= M[i][k]*M[k][j]/p
    return pivots

def run():
    pending, seen = [PROOF], set()
    while pending:
        p = pending.pop()
        if p in seen: continue
        seen.add(p); d = json.loads((ROOT/p).read_text())
        for path, h in d.get('sources', {}).items():
            assert sha(path) == h, path
            if path.endswith('.json'): pending.append(path)
    proof = json.loads((ROOT/PROOF).read_text())
    base = json.loads((ROOT/BASE).read_text()); src = json.loads((ROOT/SOURCE).read_text())
    Q = [[decode(x) for x in row] for row in base['form']]
    N = [[decode(x) for x in row] for row in base['final_normals']]
    ambient_pivots = [ev(x) for x in ldl(Q)]
    assert sum(x.l > 0 for x in ambient_pivots) == 9
    assert sum(x.h < 0 for x in ambient_pivots) == 1
    norms = [pair(n, n, Q) for n in N]; assert all(ev(x).l > 0 for x in norms)
    assert all(x == 2*O for x in norms[:11])
    SQ = [[source_decode(x) for x in row] for row in src['lattice_gram']]
    SN = [[source_decode(x) for x in row] for row in src['simple_roots']]
    SH = [source_decode(x) for x in src['geometry']['interior_point']]
    vertices = []
    for vi in (32, 33):
        v = null_one([mv(SQ, SN[i]) for i in src['geometry']['vertices'][vi]])
        if ev(pair(v, SH, SQ)).l > 0: v = [-x for x in v]
        assert ev(pair(v, SH, SQ)).h < 0
        assert ev(pair(v, v, SQ)).h < 0
        assert all(not pair(v, n, SQ) or ev(pair(v, n, SQ)).h < 0 for n in SN)
        vertices.append(v)
    B = K.from_sympy(SQ5-1); DB = K.from_sympy(17+8*SQ5)
    g = pair(*vertices, SQ); n0, n1 = [pair(v, v, SQ) for v in vertices]
    assert ev(g).h < 0 and g*g == DB*DB*n0*n1
    # This identity provides the exact diameter threshold for C=2.
    assert 2*(O+4*O/B)**2-O == DB
    bi, di = ev(B), ev(DB); L2 = 2*(B+2*O)
    assert (di-(1+4/bi)).l > 0
    dotted, C = [], {}
    for i, j in combinations(range(12), 2):
        g = pair(N[i], N[j], Q)
        assert not g or ev(g).h < 0
        c = -2*ev(g)/(ev(norms[i])*ev(norms[j])).sqrt()
        if (c-2).l > 0: dotted.append([i, j])
        else: assert (2-c).l > 0
        C[i, j] = c
    assert dotted == [[4, 11], [5, 6], [8, 11]] == proof['required_dotted_pairs']
    def union_cosh(x, y, c):
        d = (((x+2)*(y+2)).sqrt()+c)/(x*y).sqrt()
        assert d.l > 1
        radius_sum_cosh = ((x+2)*(y+2)/(x*y)).sqrt()+2/(x*y).sqrt()
        radius_sum_sinh = (2*(y+2)/(x*y)).sqrt()+(2*(x+2)/(x*y)).sqrt()
        return d*radius_sum_cosh+(d*d-1).sqrt()*radius_sum_sinh
    records = []
    for pair_ids, support, cutoff in [((4, 11), [0, 2, 9], 4),
                                      ((5, 6), [0, 1, 2, 3], 3),
                                      ((8, 11), [0, 2, 5], 3)]:
        i, j = pair_ids; c = C[i, j]
        assert all(not pair(N[k], N[l], Q) for k in support for l in pair_ids)
        G = [[pair(N[k], N[l], Q) for l in support] for k in support]
        pivots = ldl(G); assert all(ev(p).l > 0 for p in pivots)
        J = inverse(G)
        assert all(not x or ev(x).l > 0 for row in J for x in row)
        inverse_sum = sum((x for row in J for x in row), Z)
        expected_sum = K.from_sympy(s.Rational(3, 2) if len(support) == 3 else 166+72*SQ5)
        assert inverse_sum == expected_sum
        x = I(F(cutoff*cutoff, 2)-2); assert (x-bi).l > 0
        diameter_margin = di-union_cosh(x, bi, c); assert diameter_margin.l > 0
        lower = ev(L2*inverse_sum)-2*cutoff*cutoff/(c-2)
        assert lower.l > 2
        source_row = next(row for row in proof['three_pair_obstructions'] if row['pair'] == [i, j])
        assert source_row['positive_support'] == support and source_row['anchor_contact_upper'] == cutoff
        records.append(dict(pair=[i, j],support=support,positive_pivots=encode(pivots),
                            inverse_entry_sum=encode(inverse_sum),diameter_margin=diameter_margin.enc(),
                            cap_norm_lower=lower.enc()))
    # Monotonicity is uniform over the entire positive domain, not sampled.
    x, y, c, a, b = s.symbols('x y c a b', positive=True)
    d = (s.sqrt((x+2)*(y+2))+c)/s.sqrt(x*y)
    derivative = -(2*s.sqrt(y+2)/s.sqrt(x+2)+c)/(2*x*s.sqrt(x*y))
    assert s.simplify(s.diff(d, x)-derivative) == 0
    G = s.Matrix([[2, -c], [-c, 2]]); v = s.Matrix([a, b])
    penalty = 2*(a*a+b*b+c*a*b)/(c*c-4)
    assert s.simplify((v.T*G.inv()*v)[0]+penalty) == 0
    assert s.simplify(penalty.subs(b, a)-2*a*a/(c-2)) == 0
    paths = [PROOF, BASE, SOURCE, f'scripts/{TAG}.py',
             'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_ALL_HIGH_THREE_PAIR_PROJECTION_REVIEW',
                sources={p: sha(p) for p in paths}, interval_bits=BITS,
                reconstructed_source_vertices=encode(vertices),pair_comparisons=66,
                required_dotted_pairs=dotted,three_pair_obstructions=records,
                exact_threshold_identity='2*(1+4/(sqrt5-1))^2-1=17+8sqrt5',
                scope='No producer or Fraction-field/128-bit producer arithmetic imports. Separate SymPy algebraic field, physical nullspaces and rational192-bit intervals reconstruct both source vertices, the66old pairs, all three positive supports and strict norm contradictions. Uniform monotonicity and Lorentz inverse identities checked symbolically. Only coverage by the fixed twelve old-wall balls is excluded; extra nonperpendicular covering walls remain outside scope.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    result = run(); out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], len(result['three_pair_obstructions']))
