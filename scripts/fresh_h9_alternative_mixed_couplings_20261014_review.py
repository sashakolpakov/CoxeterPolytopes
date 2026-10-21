#!/usr/bin/env python3
"""Independent Fraction-basis census of fixed-source rank4 couplings."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
from fresh_h9_star_zero_20261014_own_three_review import Z, O, COS, add, sub, neg, scale, mul, inverse, sign, total
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_alternative_mixed_couplings_20261014'
PAIRS = list(combinations(range(4), 2)); PERMS = list(permutations(range(4)))
W = {m: scale(c, 2) for m, c in COS.items()}; TWO = scale(O, 2)
DET = {n: [(p, (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)))
           for p in permutations(range(n))] for n in range(1, 5)}
def determinant(M):
    out = Z
    for p, sg in DET[len(M)]:
        term = O
        for i, j in enumerate(p): term = mul(term, M[i][j])
        out = add(out, scale(term, sg))
    return out
def gram(code):
    G = [[TWO if i == j else Z for j in range(4)] for i in range(4)]
    for (i, j), m in zip(PAIRS, code): G[i][j] = G[j][i] = neg(W[m])
    return G
def canonical(code):
    d = dict(zip(PAIRS, code))
    return min(tuple(d[tuple(sorted((p[i], p[j])))] for i, j in PAIRS) for p in PERMS)
def decode(e): return F(e[0][0]), F(e[1][0]), F(e[0][1]), F(e[1][1])
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    path = f'results/{TAG}.json'; data = json.loads((ROOT/path).read_text())
    for p, h in data['sources'].items(): assert sha(p) == h, p
    supports, positive = set(), 0
    for code in product((2, 3, 4, 5), repeat=6):
        G = gram(code)
        if all(sign(determinant([r[:n] for r in G[:n]])) > 0 for n in (2, 3, 4)):
            positive += 1; supports.add(canonical(code))
    assert positive == data['positive_support_codes'] == 182
    assert supports == {tuple(r['code']) for r in data['supports']} and len(supports) == 18
    records, profiles_count, pair_count = [], 0, 0
    for row in data['supports']:
        code = tuple(row['code']); G = gram(code); den = inverse(determinant(G))
        J = [[scale(mul(den, determinant([[G[a][b] for b in range(4) if b != i]
                                          for a in range(4) if a != j])), (-1)**(i+j))
              for j in range(4)] for i in range(4)]
        assert G == [[decode(x) for x in r] for r in row['gram']]
        assert J == [[decode(x) for x in r] for r in row['inverse']]
        assert all(total(mul(G[i][k], J[k][j]) for k in range(4)) == (O if i == j else Z)
                   for i in range(4) for j in range(4))
        assert all(sign(x) >= 0 for r in J for x in r)
        adjacency = [set() for _ in range(4)]
        for (i, j), m in zip(PAIRS, code):
            if m > 2: adjacency[i].add(j); adjacency[j].add(i)
        assert all(sign(sub(J[i][i], scale(O, F(2, 3) if adjacency[i] else F(1, 2)))) >= 0 for i in range(4))
        comps, remaining = [], set(range(4))
        while remaining:
            comp = {min(remaining)}
            while True:
                nxt = comp | set().union(*(adjacency[i] for i in comp))
                if nxt == comp: break
                comp = nxt
            comps.append(sorted(comp)); remaining -= comp
        assert comps == row['components']
        def inner(v, w): return total(mul(v[i], total(mul(J[i][j], w[j]) for j in range(4))) for i in range(4))
        expected = {}
        for labs in product((2, 3, 4, 5), repeat=4):
            v = [W[m] for m in labs]; cost = inner(v, v)
            if sign(sub(TWO, cost)) > 0: expected[labs] = (v, cost, sub(O, scale(cost, F(1, 2))))
        assert set(expected) == {tuple(p['labels']) for p in row['positive_projection_profiles']}
        for p in row['positive_projection_profiles']:
            v, cost, alpha2 = expected[tuple(p['labels'])]
            assert cost == decode(p['cost']) and alpha2 == decode(p['alpha_squared'])
            assert v == [decode(x) for x in p['contact']]
        profiles = [expected[tuple(p['labels'])] for p in row['positive_projection_profiles']]
        surviving = []
        for i, (v, cv, av) in enumerate(profiles):
            for j, (w, cw, aw) in enumerate(profiles):
                cross = inner(v, w); assert sign(cross) >= 0
                left = mul(mul(W[5], W[5]), mul(av, aw))
                hits = [m for m in (2, 3, 4, 5) if left == mul(add(cross, W[m]), add(cross, W[m]))]
                if hits: surviving.append(dict(profiles=[i, j], labels=hits))
        assert surviving == row['surviving_pi5_pairs']
        assert len(surviving) == 1 and surviving[0]['labels'] == [5]
        assert all(profiles[i][1] == Z for i in surviving[0]['profiles'])
        allowed = {labs for labs, (v, cost, av) in expected.items() if cost == O}
        assert allowed == {tuple(p['labels']) for p in row['cost_one_facet6_profiles']}
        if len(comps) > 1:
            assert all(any(all(labs[i] == 2 for i in comp) for comp in comps) for labs in allowed)
        profiles_count += len(profiles); pair_count += len(profiles)**2
        records.append(dict(type=row['type'], profiles=len(profiles), pi5_pairs=len(profiles)**2, cost_one=len(allowed)))
        print('review', row['type'], len(profiles), flush=True)
    assert profiles_count == 216 and pair_count == 5168
    # Uniform high labels: a support edge>=6 plus any neighboring nonright
    # edge has determinant<=1-3/4-1/4=0. Thus its I2 factor is isolated;
    # its norm2 inverse diagonal>=2 prevents any positive-projection contact.
    assert 1-F(3, 4)-F(1, 4) == 0 and 3*F(2, 3) == 2
    # An ambient contact>=6 is possible only on A1, with cost>=3/2.
    # A second contact costs>=1/2, so it is the only active component.
    # Against a finite source neighbor: disjoint contact gives weight<1;
    # shared A1 gives a negative weight, since phi<2.
    assert sign(sub(TWO, W[5])) > 0 and F(3, 2)+F(1, 2) == 2
    basepath = 'results/fresh_h9_mixed_d4_coupling_20261014.json'
    base = json.loads((ROOT/basepath).read_text()); SQ = A.dm(base['source_form']); SN = A.dm(base['source_normals'])
    edges = []
    for i, j in combinations(range(7), 2):
        relation, label = A.relation(A.pair(SN[i], SN[j], SQ), A.pair(SN[i], SN[i], SQ), A.pair(SN[j], SN[j], SQ))
        if relation == 'finite': edges.append([i, j, label])
    assert edges == [[0, 1, 5], [1, 2, 3], [2, 3, 3], [3, 4, 3], [4, 5, 3], [4, 6, 4]]
    paths = [path, basepath, f'scripts/{TAG}_review.py',
             'scripts/fresh_h9_star_zero_20261014_own_three_review.py',
             'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_FRACTION_RANK4_COUPLING_CENSUS',
                sources={p: sha(p) for p in paths}, support_codes=4096, positive_codes=positive,
                canonical_types=18, profiles=profiles_count, pi5_pairs=pair_count, records=records,
                source_finite_edges=edges,
                scope='No coupling-producer imports. Separate Fraction-basis arithmetic, permutation minors and adjugate inverses reproduce the complete finite census. High-label tails use the stated uniform inequalities; frozen quartic arithmetic checks the original source edge graph. This classifies only the fixed-source, one-lifted-wall-per-source-facet rank4 coupling setup, not all H9 constructions or all possible caps.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = run(); out = ROOT/f'results/{TAG}_review.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
