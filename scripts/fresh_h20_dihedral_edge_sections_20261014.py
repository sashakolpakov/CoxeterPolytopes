#!/usr/bin/env python3
"""Exact algebra for the all-label m>=7 one-contact A2 section obstruction.

No numerical producer is imported. The separate Markdown gives the compact
face and invariant-lattice argument; this file checks its literal source,
uniform pair table, and the only exceptional cyclotomic equality.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from collections import Counter
import argparse, hashlib, json
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_dihedral_edge_sections_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'
OUT = ROOT / f'results/{TAG}.json'
q = S.Symbol('q')
lo, hi = S.Rational(2, 3), S.Rational(3, 4)


def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def positive(expr, a=lo, b=hi):
    """Prove a rational function positive on the entire open interval."""
    num, den = S.fraction(S.cancel(expr))
    assert (num*den).subs(q, (a+b)/2) > 0
    for poly in [num, den]:
        p = S.Poly(poly, q)
        for endpoint in [a, b]:
            while p.degree() and p.eval(endpoint) == 0:
                p = p.exquo(S.Poly(q-endpoint, q))
        assert p.count_roots(a, b) == 0
    return str(S.factor(expr))


def run():
    raw = json.loads((ROOT / RAW).read_text())['coxeter_input']
    Q = -S.Matrix(raw['basis_gram']).applyfunc(S.Rational)
    W = S.Matrix(raw['wall_coordinates_as_columns_transposed']).applyfunc(S.Rational)
    G, C = W*Q*W.T, W*Q
    assert W.shape == (210, 22) and all(x.q == 1 for x in W)
    norms = [G[i, i] for i in range(210)]
    assert Counter(norms) == {2: 42, 4: 168}
    assert all(G[i, j] <= 0 for i in range(210) for j in range(i))
    assert all((2*C[i, j]/norms[i]).q == 1 for i in range(210) for j in range(22))
    pools = {}
    for d1, d2 in [(2, 2), (2, 4), (4, 4)]:
        pools[(d1, d2)] = sorted({G[i, j] for i in range(210) for j in range(i)
                                 if sorted([norms[i], norms[j]]) == [d1, d2]})
    assert pools == {(2, 2): [-1, 0], (2, 4): [-2, 0], (4, 4): [-8, -6, -4, -2]}
    factors = S.factor_list(Q.charpoly().as_poly())[1]
    signature = [sum(m*p.count_roots(0, S.oo) for p, m in factors),
                 sum(m*p.count_roots(-S.oo, 0) for p, m in factors)]
    assert signature == [21, 1]
    A = S.Matrix([[2, -1], [-1, 2]])
    assert A.inv()[0, 0] == lo
    c = S.Symbol('c')
    qm = 2*c/(4*c-1)
    assert S.simplify(1/(2*(2-1/qm))-c) == 0
    kinds = {'rational_short': (2, 0), 'rational_long': (4, 0),
             'modified_short': (2, -1), 'modified_long': (4, -2)}
    for d, k in kinds.values():
        positive(d-q*k*k)
    positive(4*q-2)
    positive(16*q-4)
    rows = []
    for ka, kb in combinations_with_replacement(kinds, 2):
        d, a = kinds[ka]; e, b = kinds[kb]
        if a == b == 0:
            continue
        for g in pools[tuple(sorted([d, e]))]:
            n1, n2, p = d-q*a*a, e-q*b*b, g-q*a*b
            ratio = S.cancel(p*p/(n1*n2))
            row = {'types': [ka, kb], 'old_pairing': int(g),
                   'projected_norms': [str(n1), str(n2)],
                   'projected_pairing': str(p), 'squared_cosine': str(ratio)}
            if p == 0:
                assert (a == 0) != (b == 0)
                row['status'] = 'RIGHT_TO_ORIGINAL_RATIONAL_ROOT'
            elif (ka, kb, g) == ('rational_long', 'modified_long', -2):
                assert S.cancel(ratio-1/(4*(1-q))) == 0
                positive(ratio-S.Rational(3, 4)); positive(1-ratio)
                row['status'] = 'ONLY_POSSIBLE_NONRIGHT_CASE'
            else:
                try:
                    positive(ratio-1)
                    row['status'] = 'STRICTLY_DISJOINT'
                except AssertionError:
                    if S.cancel(ratio-1/(2-q)) == 0:
                        lower, upper = S.Rational(3, 4), S.Rational(4, 5)
                    else:
                        lower, upper = S.Rational(1, 4), S.Rational(1, 2)
                    positive(ratio-lower); positive(upper-ratio)
                    row.update(status='STRICT_GAP_BETWEEN_COXETER_ANGLES',
                               strict_bounds=[str(lower), str(upper)])
            rows.append(row)
    # The exceptional finite contact is equivalent to
    # cos(2pi/m)*cos(2pi/n)=1/2. For m,n>=7 its sole solution is (8,8).
    x, y = S.symbols('x y')
    assert S.cancel((1/(4*(1-qm))-(1+y)/2).subs(c, (1+x)/2)) == (1-2*x*y)/(4*x)
    f7 = 8*x**3+4*x**2-4*x-1
    f9 = 8*x**3-6*x+1
    assert S.diff(f7, x).subs(x, S.Rational(1, 2)) > 0
    assert S.diff(f9, x).subs(x, S.Rational(1, 2)) == 0
    assert S.diff(f7, x, 2).subs(x, S.Rational(1, 2)) > 0
    assert S.diff(f9, x, 2).subs(x, S.Rational(1, 2)) > 0
    assert f7.subs(x, S.Rational(3, 5)) < 0 < f7.subs(x, S.Rational(5, 8))
    assert f9.subs(x, S.Rational(1, 2)) < 0 < f9.subs(x, S.Rational(4, 5))
    r = (S.sqrt(5)-1)/2
    assert S.simplify(f7.subs(x, r)) == 4*S.sqrt(5)-9
    assert 80 < 81 and S.simplify(r*(1+S.sqrt(5))/4) == S.Rational(1, 2)
    assert S.Rational(5, 8)*S.Rational(4, 5) == S.Rational(1, 2)
    q8 = S.simplify(qm.subs(c, (2+S.sqrt(2))/4))
    assert S.simplify(1/(4*(1-q8))-(2+S.sqrt(2))/4) == 0
    # Every right neighbor of a modified long root is an unchanged short root.
    ml_right = [r for r in rows if 'modified_long' in r['types'] and r['status'].startswith('RIGHT')]
    assert len(ml_right) == 1 and ml_right[0]['types'] == ['rational_short', 'modified_long']
    walk = json.loads((ROOT/WALK).read_text())
    assert walk['status'] == 'PASS_COMPLETE_CORNER_ORBIT_CENSUS'
    seen, components = set(), 0
    for corner in walk['corner_orbits']:
        for edge in corner['edge_orbits']:
            E = tuple(edge['edge_walls'])
            if E in seen:
                continue
            seen.add(E); left = set(E)
            while left:
                cc = [min(left)]; left.remove(cc[0])
                for i in cc:
                    nn = [j for j in sorted(left) if G[i, j]]
                    cc += nn; left.difference_update(nn)
                if len(cc) == 2 and all(G[i, i] == 2 for i in cc) and G[cc[0], cc[1]] == -1:
                    components += 1
    assert len(seen) == 6492 and components == 2524
    sources = [RAW, WALK, f'scripts/{TAG}.py']
    return {'status': 'PASS_ALL_LABEL_M_GE_7_ONE_CONTACT_A2_FACE_GATE',
            'sources': {p: sha(p) for p in sources}, 'source_signature': list(map(int, signature)),
            'all210_source_reflections_preserve_Z22': True,
            'old_pairing_pools': {str(k): list(map(int, v)) for k, v in pools.items()},
            'section_norm_formula': str(qm), 'open_section_norm_interval': ['2/3', '3/4'],
            'modified_facet_pair_table': rows, 'pair_status_counts': dict(Counter(r['status'] for r in rows)),
            'only_nonright_label_exception': [8, 8],
            'cyclotomic_ordering_certificate': {'f7': str(f7), 'f9': str(f9),
                'cos_2pi_over7_bounds': ['(sqrt(5)-1)/2', '5/8'],
                'cos_2pi_over9_upper_bound': '4/5',
                'cos_squared_pi_over7_strict_lower_bound': '4/5'},
            'literal_source_edges': len(seen), 'norm2_A2_components': components,
            'anchored_continuous_families': 2*components,
            'excluded_labels': 'Every integer m>=7', 'not_excluded_here': [4, 5, 6],
            'scope': 'Exact uniform algebra plus the compact-face and rational invariant-lattice proof in the report. No enumeration or exclusion of m=4,5,6 is claimed by this certificate.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    out = run()
    if args.verify:
        assert out == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status']); print(out['pair_status_counts'])
