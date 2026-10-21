#!/usr/bin/env python3
"""One bounded exact short-root/profile search; coordinator execution only.

Use run_bounded_research_check_20261018.py --parent-lock -- THIS --write
or --verify. No inner lock, numerical packages, producer imports or solver.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path
import signal
from check_arithmetic_dn_opposed_fork_models_20261019 import field, inverse

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_omega2_finite_profiles_20261019.json'
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
    'results/compact_h9_discriminant19_edge_pencil_obstructions_20261019.md':
    '41e78966ecb8137714347e79ddc5e434044409a85d5d98bcecd437350e8ace0a',
    'results/compact_h9_discriminant19_edge_profile_congruences_20261019.md':
    '3fac696b113f386e60ae15debeb2ff3f643b8681034086cd306023c6b4f6cdff',
}


def native(a=0, b=0):
    return Q(F(a)+F(b)/2, F(b)/2)


def decode(pair):
    return Q(F(pair[0]), F(pair[1]))


def dot(left, right):
    return sum((a*b for a, b in zip(left, right)), Q())


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
    if value.sign() <= 0:
        return Q() if value == 0 else None
    norm_sqrt = rational_sqrt(value.a*value.a-5*value.b*value.b)
    if norm_sqrt is None:
        return None
    for sign in (1, -1):
        c = rational_sqrt((value.a+sign*norm_sqrt)/2)
        if c is None:
            continue
        if c:
            candidate = Q(c, value.b/(2*c))
        elif value.b == 0:
            e = rational_sqrt(value.a/5)
            if e is None:
                continue
            candidate = Q(0, e)
        else:
            continue
        if candidate*candidate == value:
            return candidate if candidate.sign() > 0 else -candidate
    return None


def integer_coefficients(value):
    a, b = value.a-value.b, 2*value.b
    assert a.denominator == b.denominator == 1
    return int(a), int(b)


def inputs():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    saved = json.loads((ROOT/'results/h19_twenty_six_wall_seed_20261019.json').read_text())
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
    ci = [row[2:] for row in gi[2:]]
    old = json.loads((ROOT/'data/candidates/h19_vinberg_prefix.json').read_text())
    roots = {f'old{i}': [native(v['a'], v['b']) for v in row]
             for i, row in enumerate(old['output']['roots_in_vinberg_order'])
             if i not in (10, 20, 21)}
    for name, encoded in saved['new_covectors'].items():
        roots[name] = mv(gi, [decode(v) for v in encoded])
    assert len(roots) == 26 and set(roots) == set(saved['norms'])
    covectors = {name: mv(G, v) for name, v in roots.items()}
    norms = {name: dot(v, covectors[name]) for name, v in roots.items()}
    for name, root in roots.items():
        assert norms[name] == decode(saved['norms'][name]) and norms[name].sign() > 0
        assert all(v.integral() for v in root)
        assert all((2*v/norms[name]).integral() for v in covectors[name])
    # Fixed roots are negative coordinate vectors; eta below is a magnitude.
    assert roots['old0'] == [Q(), Q(-1)]+[Q()]*8
    for j in range(8):
        assert roots[f'old{j+1}'] == [Q(), Q()]+[Q(-int(i == j)) for i in range(8)]
    omega = [row[2] for row in ci]
    assert omega == list(map(Q, [10, 20, 30, 24, 18, 12, 6, 15]))
    p = [Q(2), -phi]+[Q()]*8
    c = p[:2]+[sum(row)/10 for row in ci]
    c[1] += F(1, 10)
    assert dot(p, mv(G, p)).sign() < 0 and dot(c, mv(G, c)).sign() < 0
    assert dot(p, mv(G, c)).sign() < 0
    assert all(dot(p, v).sign() <= 0 and dot(c, v).sign() < 0 for v in covectors.values())
    edge = next(e for e in saved['initial_p_edge_first_hits'] if e['freed_wall'] == 'old3')
    assert edge['kind'] == 'spacelike'
    parameter = decode(edge['parameter'])
    ray = p[:2]+[parameter*v for v in omega]
    assert parameter.sign() > 0 and dot(ray, mv(G, ray)) == decode(edge['square'])
    assert dot(ray, mv(G, ray)).sign() > 0
    assert all(dot(ray, v).sign() <= 0 for v in covectors.values())
    return phi, G, ci, roots, norms, p, c, ray, parameter


def profiles(phi, ci):
    # W=A1(binary)+A2(nodes0,1)+A4(nodes3..6)+A1(node7).
    choices = lambda nodes: [(None, Q())]+[(i, v) for i in nodes for v in (Q(1), phi)]
    accepted, summaries, counts = [], [], Counter()
    for selection in product(choices([-1]), choices([0, 1]), choices([3, 4, 5, 6]), choices([7])):
        counts['component_choices'] += 1
        beta, eta, contacts = Q(), [Q()]*8, {}
        for node, value in selection:
            if node is None:
                continue
            contacts['old0' if node == -1 else f'old{node+1}'] = '1' if value == 1 else 'phi'
            if node == -1:
                beta = value
            else:
                eta[node] = value
        q = beta*beta/2
        q += 2*(eta[0]*eta[0]+eta[1]*eta[1])/3
        q += 4*(eta[3]*eta[3]+eta[6]*eta[6])/5
        q += 6*(eta[4]*eta[4]+eta[5]*eta[5])/5+eta[7]*eta[7]/2
        if not contacts or (q-2).sign() >= 0:
            continue
        ident = f'profile_{len(summaries)+1:02d}'
        Q0, R = dot(eta, mv(ci, eta)), dot(ci[2], eta)
        K = 60-15*beta*beta-30*Q0+R*R
        assert K == 30*(2-q) and K.sign() > 0 and K.integral()
        a, b = integer_coefficients(K)
        residue3, residue19 = (a % 3, b % 3), (a+15*b) % 19
        ra, rb = integer_coefficients(R)
        r3 = (ra % 3, rb % 3)
        # K=R^2 mod3 is vacuous. Use the ORIGINAL norm equation instead:
        # R nonzero lets T vary freely; R=0 forces D*x^2=beta^2+2Q0-4.
        ha, hb = integer_coefficients(beta*beta+2*Q0-4)
        x2mod3 = ((2*hb) % 3, (2*(ha+hb)) % 3)  # D inverse =2phi in F9
        pass3 = r3 != (0, 0) or x2mod3 in {(0, 0), (1, 0), (2, 0), (1, 1), (2, 2)}
        pass19 = residue19 in {i*i % 19 for i in range(19)}
        counts['schur_profiles'] += 1
        counts['passes_mod3'] += int(pass3)
        counts['passes_both_mod3_mod19'] += int(pass3 and pass19)
        summaries.append({'id': ident, 'contact_magnitudes': contacts, 'q': q.encoded(),
                          'Q0': Q0.encoded(), 'R': R.encoded(), 'K': K.encoded(),
                          'K_mod3': list(residue3), 'K_mod19': residue19,
                          'R_mod3': list(r3), 'required_x_square_mod3_if_R_zero':
                          list(x2mod3) if r3 == (0, 0) else None,
                          'passes_mod3': pass3, 'passes_mod19': pass19})
        if pass3 and pass19:
            accepted.append((ident, beta, eta, R, K))
    assert counts == {'component_choices': 405, 'schur_profiles': 51,
                      'passes_mod3': 40, 'passes_both_mod3_mod19': 34}
    return accepted, summaries, counts


def pair_kind(value, norm, phi):
    square = value*value/(2*norm)
    if value == 0:
        return 'right'
    if value.sign() > 0:
        return 'reject_positive'
    if (square-1).sign() > 0:
        return 'ultraparallel'
    if square == 1:
        return 'reject_parallel'
    return next((f'angle{m}' for m, v in [(3, Q(F(1, 4))), (4, Q(F(1, 2))),
                 (5, phi*phi/4), (6, Q(F(3, 4)))] if square == v), 'reject_finite_angle')


def run():
    phi, G, ci, roots, norms, p, c, ray, parameter = inputs()
    selected, summaries, counts = profiles(phi, ci)
    D, xs, solutions, survivors = 1+5*phi, [], [], []
    for b in range(301):
        # m=floor(b*(sqrt5-1)/2); six integers cover |x'|<=2 exactly.
        m = (isqrt(5*b*b)-b)//2
        for a in range(m-2, m+4):
            counts['x_coefficient_pairs'] += 1
            x = native(a, b)
            if (x.sign() > 0 and (x-600).sign() <= 0
                    and (x.conjugate()+2).sign() >= 0 and (x.conjugate()-2).sign() <= 0):
                xs.append((a, b, x))
    assert counts['x_coefficient_pairs'] == 1806
    counts['x_in_requested_domain'] = len(xs)
    for ident, beta, eta, R, K in selected:
        for a, b, x in xs:
            counts['profile_x_square_tests'] += 1
            positive_A = field_sqrt(K+15*D*x*x)
            if positive_A is None:
                continue
            counts['field_square_solutions'] += 1
            for A in (positive_A, -positive_A):
                counts['signed_A_branches'] += 1
                T, y = (A-R)/30, (beta-phi*x)/2
                if not (T.integral() and y.integral()):
                    counts['reject_T_or_y_nonintegral'] += 1
                    continue
                root = [x, y]+mv(ci, [eta[j]+(T if j == 2 else Q()) for j in range(8)])
                covector = mv(G, root)
                assert all(v.integral() for v in root+covector)
                assert dot(root, covector) == 2 and A*A-15*D*x*x == K
                assert dot(roots['old0'], covector) == -beta
                assert all(dot(roots[f'old{j+1}'], covector) == -eta[j] for j in range(8) if j != 2)
                counts['integral_norm_two_roots'] += 1
                products, invalid, categories = {}, [], Counter()
                for name, old in roots.items():
                    value = dot(old, covector)
                    kind = pair_kind(value, norms[name], phi)
                    products[name] = [value.encoded(), kind]
                    categories[kind] += 1
                    if kind.startswith('reject'):
                        invalid.append(name)
                assert sum(categories.values()) == 26
                cut, cp, pp = dot(ray, covector), dot(c, covector), dot(p, covector)
                assert cut == -D*x+parameter*A and pp == -D*x and pp.sign() < 0
                epsilon = F(1)
                while (pp+epsilon*cp).sign() >= 0:
                    epsilon /= 2
                center = [v+epsilon*w for v, w in zip(p, c)]
                center_square = dot(center, mv(G, center))
                assert center_square.sign() < 0 and dot(p, mv(G, center)).sign() < 0
                assert dot(center, covector).sign() < 0
                assert all(dot(center, mv(G, v)).sign() < 0 for v in roots.values())
                reasons = ([] if cut.sign() > 0 else ['does_not_cut_old3_ray'])
                reasons += ([] if not invalid else ['incompatible_native_pairs'])
                record = {'profile': ident, 'x_phi_coefficients': [a, b], 'A': A.encoded(),
                          'T': T.encoded(), 'y': y.encoded(), 'ray_product': cut.encoded(),
                          'old_center_product': cp.encoded(), 'old_center_strict': cp.sign() < 0,
                          'new_center_epsilon': str(epsilon), 'new_center_square': center_square.encoded(),
                          'invalid_walls': invalid, 'pair_categories': dict(categories), 'rejected': reasons}
                solutions.append(record)
                if not reasons:
                    tau = D*x/A
                    endpoint = p[:2]+[tau*row[2] for row in ci]
                    endpoint_square = dot(endpoint, mv(G, endpoint))
                    assert tau.sign() > 0 and (parameter-tau).sign() > 0
                    assert endpoint_square == -2*D*K/(A*A) and endpoint_square.sign() < 0
                    assert dot(p, mv(G, endpoint)).sign() < 0
                    endpoint_products = {n: dot(endpoint, mv(G, v)) for n, v in roots.items()}
                    assert all(v.sign() <= 0 for v in endpoint_products.values())
                    white = {'old0'} | {f'old{j+1}' for j in range(8) if j != 2}
                    assert {n for n, v in endpoint_products.items() if v == 0} == white
                    assert dot(endpoint, covector) == 0
                    survivors.append(record | {'coordinates': [v.encoded() for v in root],
                        'covector': [v.encoded() for v in covector], 'wall_products': products,
                        'strict_center': [v.encoded() for v in center],
                        'endpoint_parameter': tau.encoded(), 'endpoint_square': endpoint_square.encoded(),
                        'endpoint': [v.encoded() for v in endpoint]})
    assert counts['profile_x_square_tests'] == 34*len(xs)
    return {'status': 'PASS_BOUNDED_OMEGA2_FINITE_PROFILE_SEARCH',
        'source_sha256': digest(Path(__file__)), 'sources': PINS,
        'domain': {'x_positive_at_most': 600, 'absolute_conjugate_x_at_most': 2,
                   'b_box_inclusive': [0, 300],
                   'a_box': 'm-2<=a<=m+3; m=(isqrt(5*b*b)-b)//2',
                   'equation': 'A^2-15D*x^2=K; T=(A-R)/30; y=(beta-phi*x)/2 integral'},
        'stage_counts': dict(counts), 'profiles': summaries,
        'integral_norm_two_solutions': solutions, 'survivors': survivors,
        'scope': 'Only nonzero finite short contact profiles for the fixed old3 white star, '
                 'at most one 1/phi contact per W component and Schur q<2, in the stated x box. '
                 'Both signs of A and all26 wall products are checked for every integral root. '
                 'Old-center positivity never rejects a candidate. No other star/pencil, '
                 'root completeness, vertex preservation, finite volume, nonreflectivity or dimension claim.'}


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
