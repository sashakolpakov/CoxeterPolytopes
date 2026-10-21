#!/usr/bin/env python3
"""ONE bounded exact short-root pencil; run through the shared-lock wrapper.

Use run_bounded_research_check_20261018.py --parent-lock -- THIS --write
or --verify. The wrapper owns resource/RSS/parent-death gates; this process
adds a45second alarm. No scientific producer, numerical package or optimizer.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import signal
from check_arithmetic_dn_opposed_fork_models_20261019 import field, inverse

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_omega3_short_pencil_20261019.json'
Q = field(5)
PINS = {
    'results/h19_twenty_six_wall_seed_20261019.json':
    'e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b',
    'scripts/check_h19_twenty_six_wall_seed_20261019.py':
    'f0ac1d28225a37be17e9b61f2e436ea2fc3b1511efb32b02212b33539fecb363',
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
}


def native(a=0, b=0):
    return Q(F(a)+F(b)/2, F(b)/2)


def decode(pair):
    return Q(F(pair[0]), F(pair[1]))


def dot(left, right):
    return sum((x*y for x, y in zip(left, right)), Q())


def mv(matrix, vector):
    return [dot(row, vector) for row in matrix]


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            value.update(block)
    return value.hexdigest()


def rational_sqrt(value):
    if value < 0:
        return None
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    return F(a, b) if a*a == value.numerator and b*b == value.denominator else None


def field_sqrt(value):
    """Return the positive K-square root, or None; no floating approximation."""
    if value.sign() <= 0:
        return Q() if value == 0 else None
    norm_sqrt = rational_sqrt(value.a*value.a-5*value.b*value.b)
    if norm_sqrt is None:
        return None
    for sign in (1, -1):
        rational_part = rational_sqrt((value.a+sign*norm_sqrt)/2)
        if rational_part is None:
            continue
        if rational_part:
            candidate = Q(rational_part, value.b/(2*rational_part))
        elif value.b == 0:
            irrational_part = rational_sqrt(value.a/5)
            if irrational_part is None:
                continue
            candidate = Q(0, irrational_part)
        else:
            continue
        if candidate*candidate == value:
            return candidate if candidate.sign() > 0 else -candidate
    return None


def inputs():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    with (ROOT/'results/h19_twenty_six_wall_seed_20261019.json').open() as stream:
        saved = json.load(stream)
    assert saved['walls'] == 26 and saved['wall_pairs'] == 325
    assert saved['checker_sha256'] == PINS['scripts/check_h19_twenty_six_wall_seed_20261019.py']
    for name, expected in saved['sources'].items():
        assert digest(ROOT/name) == expected, name
    phi = native(0, 1)
    G = [[Q() for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, phi, phi, Q(2)
    for i in range(8):
        G[i+2][i+2] = Q(2)
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
        G[i+2][j+2] = G[j+2][i+2] = Q(-1)
    gi = inverse(G, Q)
    with (ROOT/'data/candidates/h19_vinberg_prefix.json').open() as stream:
        old = json.load(stream)['output']['roots_in_vinberg_order']
    roots = {f'old{i}': [native(v['a'], v['b']) for v in row]
             for i, row in enumerate(old) if i not in (10, 20, 21)}
    for name, encoded in saved['new_covectors'].items():
        assert len(encoded) == 10
        roots[name] = mv(gi, [decode(v) for v in encoded])
    assert len(roots) == 26 and set(roots) == set(saved['norms'])
    assert all(len(v) == 10 for v in roots.values())
    cov = {name: mv(G, v) for name, v in roots.items()}
    norms = {name: dot(v, cov[name]) for name, v in roots.items()}
    for name in roots:
        assert norms[name] == decode(saved['norms'][name])
        assert norms[name].sign() > 0
        assert all(v.integral() for v in roots[name])
        assert all((2*v/norms[name]).integral() for v in cov[name])
    # omega3 is reconstructed from G^-1, not a duplicated producer root list.
    omega = [gi[i+2][5] for i in range(8)]
    assert omega == list(map(Q, [8, 16, 24, 20, 15, 10, 5, 12]))
    assert dot(omega, mv([row[2:] for row in G[2:]], omega)) == 20
    p = [Q(2), -phi]+[Q()]*8
    center = list(p)
    center[1] += F(1, 10)
    for i in range(8):
        center[i+2] += sum(gi[i+2][j+2] for j in range(8))/10
    assert dot(p, mv(G, p)).sign() < 0 and dot(p, mv(G, center)).sign() < 0
    assert all(dot(p, g).sign() <= 0 for g in cov.values())
    assert dot(center, mv(G, center)).sign() < 0
    assert all(dot(center, g).sign() < 0 for g in cov.values())
    edge = next(e for e in saved['initial_p_edge_first_hits'] if e['freed_wall'] == 'old4')
    assert edge['kind'] == 'spacelike'
    parameter = decode(edge['parameter'])
    ray = p[:2]+[parameter*v for v in omega]
    assert dot(ray, mv(G, ray)) == decode(edge['square'])
    assert dot(ray, mv(G, ray)).sign() > 0
    assert all(dot(ray, g).sign() <= 0 for g in cov.values())
    return phi, G, roots, norms, p, center, ray, omega, parameter


def run():
    phi, G, roots, norms, p, center, ray, omega, parameter = inputs()
    D, stages, solutions, survivors = 1+5*phi, Counter(), [], []
    # If0<u<=300, |u'|<=2, then b=(u-u')/sqrt5 lies in(-1,151).
    # Also a=u'-b*phi' lies in[-2,b+2]. This finite box is conservative.
    for b in range(151):
        for a in range(-2, b+3):
            stages['coefficient_pairs'] += 1
            u = native(a, b)
            up = u.conjugate()
            if not (u.sign() > 0 and (u-300).sign() <= 0
                    and (up+2).sign() >= 0 and (up-2).sign() <= 0):
                continue
            stages['u_in_requested_domain'] += 1
            t = field_sqrt((1+D*u*u)/10)
            if t is None:
                continue
            stages['field_square_solutions'] += 1
            record = {'u_phi_coefficients': [a, b], 'u': u.encoded(), 't': t.encoded()}
            if not t.integral():
                record['rejected'] = ['nonintegral_t']; solutions.append(record)
                continue
            stages['integral_square_solutions'] += 1
            assert t.sign() > 0 and 10*t*t-D*u*u == 1
            assert (a+3*b) % 5 in (2, 3)
            root = [2*u, -phi*u]+[t*v for v in omega]
            covector = mv(G, root)
            assert covector == [-D*u, Q()]+[t if i == 3 else Q() for i in range(8)]
            assert dot(root, covector) == 2
            assert all(v.integral() for v in root+covector)
            products, categories, invalid = {}, Counter(), []
            for name, old in roots.items():
                value = dot(old, covector)
                square = value*value/(2*norms[name])
                if value == 0:
                    kind = 'right'
                elif value.sign() > 0:
                    kind = 'reject_positive'
                elif square.sign() > 0 and (square-1).sign() < 0:
                    kind = next((f'angle{m}' for m, v in [
                        (3, Q(F(1, 4))), (4, Q(F(1, 2))),
                        (5, phi*phi/4), (6, Q(F(3, 4)))] if square == v), 'reject_finite_angle')
                else:
                    kind = 'ultraparallel' if (square-1).sign() > 0 else 'reject_parallel'
                categories[kind] += 1
                products[name] = [value.encoded(), kind]
                if kind.startswith('reject'):
                    invalid.append(name)
            assert sum(categories.values()) == 26
            cp, cut = dot(center, covector), dot(ray, covector)
            assert cut == -2*D*u+20*t*parameter
            pp, epsilon = dot(p, covector), F(1)
            assert pp == -2*D*u and pp.sign() < 0
            while (pp+epsilon*cp).sign() >= 0:
                epsilon /= 2
            new_center = [v+epsilon*w for v, w in zip(p, center)]
            center_square = dot(new_center, mv(G, new_center))
            assert center_square.sign() < 0 and dot(p, mv(G, new_center)).sign() < 0
            assert all(dot(new_center, mv(G, v)).sign() < 0 for v in roots.values())
            new_product = dot(new_center, covector)
            assert new_product == pp+epsilon*cp and new_product.sign() < 0
            reasons = ([] if cut.sign() > 0 else ['does_not_cut_old4_ray'])
            reasons += ([] if not invalid else ['incompatible_native_pairs'])
            record.update(old_center_product=cp.encoded(), old_center_strict=cp.sign() < 0,
                          new_center_epsilon=str(epsilon), new_center_square=center_square.encoded(),
                          new_center_root_product=new_product.encoded(), ray_product=cut.encoded(),
                          invalid_walls=invalid, pair_categories=dict(categories), rejected=reasons)
            solutions.append(record)
            if not reasons:
                stages['survivors'] += 1
                survivors.append(record | {'coordinates': [v.encoded() for v in root],
                    'covector': [v.encoded() for v in covector], 'wall_products': products,
                    'strict_center': [v.encoded() for v in new_center]})
    assert stages['coefficient_pairs'] == 12080
    return {'status': 'PASS_BOUNDED_OMEGA3_SHORT_PENCIL_SEARCH',
        'source_sha256': digest(Path(__file__)), 'sources': PINS,
        'domain': {'u_positive_at_most': 300, 'absolute_conjugate_u_at_most': 2,
                   'b_box_inclusive': [0, 150], 'a_box': '-2<=a<=b+2',
                   'equation': '10*t^2-D*u^2=1; u,t integral; t positive'},
        'stage_counts': dict(stages), 'exact_field_square_solutions': solutions,
        'survivors': survivors,
        'scope': 'Only the stated bounded omega3 norm-two pencil and exact requested native '
                 'labels2/3/4/5/6 or strictly ultraparallel; parallel pairs rejected. '
                 'All26 old wall products are tested for every integral solution. '
                 'The old center is never a rejection filter: p+epsilon*c supplies a strict center. '
                 'No other pencil/norm class, root completeness, vertex preservation, '
                 'finite-volume, compactness, nonreflectivity or dimension conclusion.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.alarm(45)
    result = run()
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    assert len(encoded.encode()) < 1024**2
    if args.verify:
        assert json.loads(CERT.read_text()) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream:
            stream.write(encoded)
    signal.alarm(0)
    print(json.dumps({'status': result['status'], 'counts': result['stage_counts'],
                      'survivors': len(result['survivors'])}))


if __name__ == '__main__':
    main()
