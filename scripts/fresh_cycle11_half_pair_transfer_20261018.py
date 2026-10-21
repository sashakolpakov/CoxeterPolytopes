#!/usr/bin/env python3
"""Universal affine cycle11 distance5 half transfer and triangle exclusion.

Exact scalar/matrix corroboration for the accompanying all-rank proof.
No project arithmetic helpers or bounded spectator census are imported.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_cycle11_half_pair_transfer_20261018'
FAR = 'fresh_spectator_five_quarters_transfer_20261017'
LOW = 'fresh_spectator_energy_below_one_20261017'
CYCLE = 'fresh_cycle12_ideal5_transfer_20261017'
CURRENT = 'fresh_deep_finite8_dimension64_20261018'
STEMS = (FAR, LOW, CYCLE, CURRENT)
TARGET = F(14, 11)


def positive(expr):
    expr = s.simplify(expr)
    assert expr.is_positive is not None, expr
    return bool(expr.is_positive)


def gram(n, edges):
    g = 2*s.eye(n)
    for i, j, magnitude in edges:
        g[i, j] = g[j, i] = -magnitude
    return g


def distances(g, start):
    ds = {start: 0}
    queue = [start]
    for i in queue:
        for j in range(g.rows):
            if i != j and g[i, j] != 0 and j not in ds:
                ds[j] = ds[i]+1
                queue.append(j)
    return ds


def pd_pivots(g):
    g = g.copy()
    pivots = []
    while g.rows:
        p = s.simplify(g[0, 0])
        assert positive(p), p
        pivots.append(str(p))
        g = (g[1:, 1:]-g[1:, :1]*g[:1, 1:]/p).applyfunc(s.simplify)
    return pivots


def one_sided_spectrum():
    # Native magnitude >=phi is too costly even on the smallest A1.
    phi = (1+s.sqrt(5))/2
    assert positive(phi**2/2-s.Rational(14, 11))
    # If magnitude=sqrt2, a diagonal<=7/11<1 forces a whole A_r leaf.
    # r/(r+1)<=7/11 implies r<=7/4, so only A1 and energy1 remain.
    assert F(7, 11)/(1-F(7, 11)) == F(7, 4)
    # A-path interior: if both endpoint distances+1>=3, diagonal>=3/2.
    # Otherwise position2 has diagonal2-4/N; the upper bound forces N<=11/2.
    assert F(3, 2) > TARGET
    assert F(4)/(2-TARGET) == F(11, 2)
    interiors = [dict(N=N, energy=str(F(2*(N-2), N))) for N in (4, 5)]
    assert [r['energy'] for r in interiors] == ['1', '6/5']
    # B inverse diagonals are integers or r/2; only1 fits the bound.
    assert F(3, 2) > TARGET
    # D inverse diagonals are integers or r/4; onlyr4,5 fit.
    assert 4*TARGET == F(56, 11)
    d_leaves = [dict(rank=r, energy=str(F(r, 4))) for r in (4, 5)]
    exceptional = []
    for n in (6, 7, 8):
        edges = []
        node = 1
        for length in (n-4, 2, 1):
            previous = 0
            for _ in range(length):
                edges.append((previous, node, 1))
                previous, node = node, node+1
        assert node == n
        g = gram(n, edges)
        inv = g.inv()
        assert all(positive(v-s.Rational(14, 11)) for v in inv.diagonal())
        exceptional.append(dict(kind=f'E{n}', gram=[[str(x) for x in row] for row in g.tolist()],
                                inverse_diagonal=list(map(str, inv.diagonal()))))
    for kind, labels in (('F4', [1, s.sqrt(2), 1]),
                         ('H3', [phi, 1]), ('H4', [phi, 1, 1])):
        g = gram(len(labels)+1, [(i, i+1, q) for i, q in enumerate(labels)])
        inv = g.inv().applyfunc(s.simplify)
        assert all(positive(v-s.Rational(14, 11)) for v in inv.diagonal())
        exceptional.append(dict(kind=kind, gram=[[str(x) for x in row] for row in g.tolist()],
                                inverse_diagonal=list(map(str, inv.diagonal()))))
    # Every I2(m),m>=5 has at least the label5 inverse diagonal.
    i2_min = s.simplify(2/(4-phi**2))
    assert positive(i2_min-s.Rational(14, 11))
    # With two active components, each is<=14/11-1/2=17/22<1.
    cap = TARGET-F(1, 2)
    assert cap == F(17, 22) and cap/(1-cap) == F(17, 5)
    endpoint_values = [F(r, r+1) for r in range(1, 4)]
    sums = [dict(ranks=[r, t], energy=str(F(r, r+1)+F(t, t+1)))
            for r, t in combinations_with_replacement(range(1, 4), 2)]
    assert all(F(row['energy']) != TARGET for row in sums)
    assert TARGET not in (F(1), F(6, 5), F(5, 4)) and TARGET > 1
    assert F(4, 3) > TARGET and F(3, 2) > TARGET
    return dict(bound=str(TARGET), one_component_spectrum=['r/(r+1), integer r>=1', '1', '6/5', '5/4'],
                nonunit_label5_minimum=str(s.simplify(phi**2/2)),
                nonunit_label4_forced_rank_maximum='7/4',
                A_interior_cases=interiors, D_short_leaf_cases=d_leaves,
                exceptional_exclusions=exceptional, I2_label5_inverse_minimum=str(i2_min),
                two_component_cap=str(cap), two_component_endpoint_energies=list(map(str, endpoint_values)),
                two_component_sums=sums,
                theorem='At difference energy14/11 under separate strict edge extensions, disjoint contacts force a mixed-component path of length<=5 between the new roots. Otherwise one or two one-sided component energies cannot sum to14/11.',
                arbitrary_rank='Classical all-rank inverse formulas and the inherited mixed-component4/3 lower bound; exceptional finite types checked exactly; all I2 labels>=5 excluded monotonically.')


def cycle_profiles():
    # Old cycle0..10, a=0,b=5; newc=11. The displayed actual corner
    # has the10 old nodes other than a plus c.
    records = []
    for side in (0, 1):
        for name, pos, q in (('end3', 1, 1), ('end4', 1, s.sqrt(2)), ('next3', 2, 1)):
            contact = pos if side == 0 else 11-pos
            ids = list(range(1, 12))
            g = gram(11, [(i-1, i, 1) for i in range(1, 10)]+[(10, contact-1, q)])
            ds = distances(g, 10)
            assert ds[4] <= 6
            leaf = 9 if side == 0 else 0
            neighbours = [j for j in range(11) if j != leaf and g[leaf, j]]
            assert len(neighbours) == 1
            near = neighbours[0]
            assert sum(j != near and g[near, j] != 0 for j in range(11)) == 2
            records.append(dict(profile=f'{name}_{side}', native_contact=str(q),
                                marked_distance=ds[4], fixed_leaf=ids[leaf],
                                fixed_degree_two_neighbor=ids[near], positive_pivots=pd_pivots(g)))
    both = gram(11, [(i, i+1, 1) for i in range(10)]+[(10, 0, 1)])
    assert both*s.ones(11, 1) == s.zeros(11, 1)
    assert both.rank() == 10 and distances(both, 10)[4] == 5
    # In an actual PSD corner the positive kernel forces every external
    # contact to this cycle to vanish. Its whole factor is therefore11.
    return dict(single_profiles=records, both_ends_Gram=[[str(x) for x in row] for row in both.tolist()],
                both_ends_positive_kernel=[1]*11, both_ends_marked_distance=5,
                both_ends_payments=['original cycle11 half', 'different actual cycle11 half'])


def height_cases():
    alpha = s.Rational(11, 30)
    h, k, p, q, r = s.symbols('h k p q r', nonzero=True)
    core = s.Matrix([[alpha, -alpha, -h, 0], [-alpha, alpha, 0, -k],
                     [-h, 0, 2-p, -q], [0, -k, -q, 2-r]])
    formula = 1-2*alpha/h**2-2*alpha/k**2+alpha*(p/h**2-2*q/(h*k)+r/k**2)
    assert s.factor(core.det()/(h*h*k*k)-formula) == 0
    assert 1-2*alpha == s.Rational(4, 15) > 0  # both heights>=sqrt2
    assert 1-2*alpha-s.Rational(2, 3)*alpha == s.Rational(1, 45) > 0
    phi = (1+s.sqrt(5))/2
    Dphi = s.simplify(2+2/phi**2-1/alpha)
    assert Dphi == s.Rational(25, 11)-s.sqrt(5)
    assert positive(Dphi) and positive(1/(2*phi**2)-Dphi)
    Dsqrt = s.simplify(2+1-1/alpha)
    assert Dsqrt == s.Rational(3, 11) < s.Rational(1, 2)
    forced_energy = 2*Dsqrt
    forced_rank = forced_energy/(1-forced_energy)
    assert forced_energy == s.Rational(6, 11) and forced_rank == s.Rational(6, 5)
    assert 4-1/alpha == s.Rational(14, 11)
    return dict(alpha='11/30', determinant_identity=str(formula),
                both_heights_at_least_sqrt2_gap='4/15', one_height_at_least_sqrt3_gap='1/45',
                cases=[dict(heights=['1', '1'], required_energy='14/11',
                            resolution='mixed spectator component supplies opposite-black path<=5'),
                       dict(heights=['1', 'sqrt(2)'], required_energy='3/11',
                            resolution='height1 support empty by coordinate Cauchy; remaining one-sided unscaled energy6/11 would force nonintegralA rank6/5'),
                       dict(heights=['1', 'phi'], required_energy=str(Dphi),
                            resolution='required energy below every nonempty coordinate Cauchy bound, but positive')],
                phi_coordinate_gap=str(s.simplify(1/(2*phi**2)-Dphi)))


def triangle():
    a, b = 0, 5
    allowed_a = {(a+d) % 11 for d in (-2, -1, 1, 2)}
    allowed_b = {(b+d) % 11 for d in (-2, -1, 1, 2)}
    assert allowed_a == {1, 2, 9, 10} and allowed_b == {3, 4, 6, 7}
    assert not allowed_a & allowed_b and b not in allowed_a and a not in allowed_b
    return dict(old_marks=[a, b], allowed_after_delete_a=sorted(allowed_a),
                allowed_after_delete_b=sorted(allowed_b), common_allowed_contacts=[],
                proof='Complete PSD one-root extensions of each A10 allow only its endpoint/second-node contacts or both endpoints. Their allowed sets are disjoint; the other marked node is interior. The third side is orthogonal to the entire affine factor, contradicting its positive cusp height.',
                conclusion='No actual triangle has varying affine cycle11 and marked distance5.')


def run():
    data = {stem: json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem, packet in data.items():
        assert packet['status'].startswith('PASS'), stem
        for path, digest in packet['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (stem, path)
    assert data[FAR]['same_component_far_lower_bound'] == '4/3'
    spectrum = one_sided_spectrum()
    profiles = cycle_profiles()
    heights = height_cases()
    triangle_record = triangle()
    weights = [dict(pair=[i, j], distance=min(j-i, 11-j+i),
                    weight='1/2' if min(j-i, 11-j+i) == 5 else '1')
               for i, j in combinations(range(11), 2)]
    assert sum(F(row['weight']) for row in weights) == F(99, 2)
    assert sum(row['weight'] == '1/2' for row in weights) == 11
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md',
              ROOT/'results/fresh_distance8_ideal_corner_replacement_20261017.md']
    return dict(status='PASS_CYCLE11_DISTANCE5_HALF_TRANSFER_AND_TRIANGLE_EXCLUSION',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                spectrum=spectrum, profiles=profiles, isolated_heights=heights,
                triangle=triangle_record, weighted_pairs=weights,
                old_total_cost='55', new_total_cost='99/2', added_finite_cost='0',
                other_ideal_types_unchanged=True, order11_all_family_envelope='55',
                theorem='Halve all11 distance5 pairs in the actual affine cycle11. Every affected quadrilateral has another retained unit or a distinct half; no triangle contains an affected corner. All finite weights/coupons and other ideal weights stay unchanged.',
                scope='Universal local transfer with arbitrary supporting spectators and native labels. Classical/exceptional finite and affine classification plus inherited strict-edge energy bounds are explicit inputs. No new all-family order11 envelope, finite budget, ambient dimension bound or chamber is asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], 'new cost', out['new_total_cost'], flush=True)
