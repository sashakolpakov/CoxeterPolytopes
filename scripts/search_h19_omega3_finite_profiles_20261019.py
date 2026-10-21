#!/usr/bin/env python3
"""One bounded exact omega3 short-profile search; coordinator execution only.

Run through run_bounded_research_check_20261018.py --parent-lock -- THIS
--write or --verify. No inner lock. The pinned omega2 module supplies
only pure arithmetic helpers: its inputs/profiles/run are NEVER invoked.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json
from math import isqrt
from pathlib import Path
import signal

from search_h19_omega2_finite_profiles_20261019 import (
    Q, native, decode, dot, mv, digest, field_sqrt, integer_coefficients,
    pair_kind, inverse)

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_omega3_finite_profiles_20261019.json'
PINS = {
    'results/h19_twenty_seven_wall_seed_20261019.json':
    'b5b529b485a01b9315be6f31bbb122c170652273d3932b1eecece85680b5ce7d',
    'scripts/check_h19_twenty_seven_wall_seed_20261019.py':
    '4fed900c65b747dbf99226d16fd0e99b7e1f3106fd3269e0ebce9882455f9cf1',
    'results/h19_twenty_six_wall_seed_20261019.json':
    'e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b',
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'scripts/search_h19_omega2_finite_profiles_20261019.py':
    '02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'results/compact_h9_discriminant19_omega3_profile_congruences_20261019.md':
    '0c9a29fc83a7c1abac1330b72b771612423c0f0ea9009a63daf6f71d47a038fd',
}


def validate_record_sources(record):
    for name, expected in record['sources'].items():
        assert digest(ROOT/name) == expected, name


def fixed_inputs():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    prior = json.loads((ROOT/'results/h19_twenty_six_wall_seed_20261019.json').read_text())
    current = json.loads((ROOT/'results/h19_twenty_seven_wall_seed_20261019.json').read_text())
    found = json.loads((ROOT/'results/h19_omega2_finite_profiles_20261019.json').read_text())
    assert prior['status'] == 'PASS_EXACT_26_WALL_LOCAL_ALGEBRA'
    assert current['status'] == 'PASS_EXACT_27_WALL_LOCAL_ALGEBRA'
    assert found['status'] == 'PASS_BOUNDED_OMEGA2_FINITE_PROFILE_SEARCH'
    assert current['walls'] == 27 and current['wall_pairs'] == 351
    assert current['checker_sha256'] == PINS['scripts/check_h19_twenty_seven_wall_seed_20261019.py']
    for record in (prior, current, found):
        validate_record_sources(record)
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
    for name, encoded in prior['new_covectors'].items():
        roots[name] = mv(gi, [decode(v) for v in encoded])
    assert len(roots) == 26 and set(roots) == set(prior['norms'])
    roots['q'] = mv(gi, [decode(v) for v in current['new_covector']])
    assert roots['q'] == [decode(v) for v in current['new_root']]
    assert roots['q'] == [native(1, 3), native(-1, -2)]+[
        native(a, b) for a, b in [(7, 10), (14, 20), (20, 30), (16, 24),
                                 (12, 18), (8, 12), (4, 6), (10, 15)]]
    covectors = {name: mv(G, root) for name, root in roots.items()}
    norms = {name: dot(root, covectors[name]) for name, root in roots.items()}
    assert len(roots) == 27
    for name, root in roots.items():
        expected = Q(2) if name == 'q' else decode(prior['norms'][name])
        assert norms[name] == expected and expected.sign() > 0
        assert all(v.integral() for v in root)
        assert all((2*v/expected).integral() for v in covectors[name])
    assert roots['old0'] == [Q(), Q(-1)]+[Q()]*8
    for j in range(8):
        assert roots[f'old{j+1}'] == [Q(), Q()]+[Q(-int(i == j)) for i in range(8)]
    omega = [row[3] for row in ci]
    assert omega == list(map(Q, [8, 16, 24, 20, 15, 10, 5, 12]))
    p = [Q(2), -phi]+[Q()]*8
    c = p[:2]+[sum(row, Q())/10 for row in ci]
    c[1] += F(1, 10)
    assert dot(p, mv(G, p)).sign() < 0 and dot(c, mv(G, c)).sign() < 0
    assert dot(p, mv(G, c)).sign() < 0
    assert all(dot(p, v).sign() <= 0 and dot(c, v).sign() < 0 for v in covectors.values())
    edge = next(e for e in current['initial_p_edge_first_hits'] if e['freed_wall'] == 'old4')
    assert edge['kind'] == 'spacelike'
    parameter = decode(edge['parameter'])
    ray = p[:2]+[parameter*v for v in omega]
    assert parameter.sign() > 0 and dot(ray, mv(G, ray)) == decode(edge['square'])
    assert dot(ray, mv(G, ray)).sign() > 0
    assert all(dot(ray, v).sign() <= 0 for v in covectors.values())
    return phi, G, ci, roots, norms, p, c, ray, parameter


def finite_profiles(phi, ci):
    # W=A1(binary)+A4(nodes0,1,2,7)+A3(nodes4,5,6).
    def choices(nodes):
        return [(None, Q())]+[(i, v) for i in nodes for v in (Q(1), phi)]
    accepted, summaries, counts = [], [], Counter()
    square19 = {i*i % 19 for i in range(19)}
    square5 = {0, 1, 4}
    for selection in product(choices([-1]), choices([0, 1, 2, 7]), choices([4, 5, 6])):
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
        qW = beta*beta/2
        qW += 4*(eta[0]*eta[0]+eta[7]*eta[7])/5
        qW += 6*(eta[1]*eta[1]+eta[2]*eta[2])/5
        qW += 3*(eta[4]*eta[4]+eta[6]*eta[6])/4+eta[5]*eta[5]
        if (qW-2).sign() >= 0:
            continue
        ident = f'omega3_profile_{len(summaries)+1:02d}'
        Q0, R = dot(eta, mv(ci, eta)), dot(ci[3], eta)
        K = 40-10*beta*beta-20*Q0+R*R
        assert K == 20*(2-qW) and K.sign() > 0 and K.integral()
        ka, kb = integer_coefficients(K)
        residue19 = (ka+15*kb) % 19
        pass19 = residue19 in square19
        # If R/5 is integral then A=20T+R lies in5O. Divide the
        # norm equation by5 BEFORE reducing modulo sqrt5 (phi=3).
        divided_test = (R/5).integral()
        required_x2, pass5 = None, True
        if divided_test:
            assert (K/5).integral()
            a, b = integer_coefficients(K/5)
            required_x2 = (2*(a+3*b)) % 5  # (-2)^(-1)=2 in F5.
            pass5 = required_x2 in square5
        counts['schur_profiles_including_zero'] += 1
        counts['nonzero_schur_profiles'] += int(bool(contacts))
        counts['passes_prime_D'] += int(pass19)
        counts['passes_prime_D_and_divided_5'] += int(pass19 and pass5)
        summaries.append({'id': ident, 'contact_magnitudes': contacts,
                          'white_schur_energy': qW.encoded(), 'Q0': Q0.encoded(),
                          'R': R.encoded(), 'K': K.encoded(), 'K_mod19': residue19,
                          'passes_prime_D': pass19, 'R_in_5O': divided_test,
                          'required_x_square_mod_sqrt5': required_x2,
                          'passes_divided_5_test': pass5})
        if pass19 and pass5:
            assert contacts and all(v == '1' for v in contacts.values())
            accepted.append((ident, beta, eta, R, K))
    assert counts == {'component_choices': 189, 'schur_profiles_including_zero': 29,
                      'nonzero_schur_profiles': 28, 'passes_prime_D': 22,
                      'passes_prime_D_and_divided_5': 18}
    return accepted, summaries, counts


def run():
    phi, G, ci, roots, norms, p, c, ray, parameter = fixed_inputs()
    selected, summaries, counts = finite_profiles(phi, ci)
    D, xs, solutions, survivors = 1+5*phi, [], [], []
    for b in range(301):
        m = (isqrt(5*b*b)-b)//2
        for a in range(m-2, m+4):
            counts['x_coefficient_pairs'] += 1
            x = native(a, b)
            if (x.sign() > 0 and (x-600).sign() <= 0
                    and (x.conjugate()+2).sign() >= 0 and (x.conjugate()-2).sign() <= 0):
                xs.append((a, b, x))
    assert counts['x_coefficient_pairs'] == 1806
    counts['x_in_requested_domain'] = len(xs)
    retained_bytes = 0
    for ident, beta, eta, R, K in selected:
        for a, b, x in xs:
            counts['profile_x_square_tests'] += 1
            positive_A = field_sqrt(K+10*D*x*x)
            if positive_A is None:
                continue
            counts['field_square_solutions'] += 1
            for A in (positive_A, -positive_A):
                counts['signed_A_branches'] += 1
                T, y = (A-R)/20, (beta-phi*x)/2
                if not (T.integral() and y.integral()):
                    counts['reject_T_or_y_nonintegral'] += 1
                    continue
                root = [x, y]+mv(ci, [eta[j]+(T if j == 3 else Q()) for j in range(8)])
                covector = mv(G, root)
                assert all(v.integral() for v in root+covector)
                assert dot(root, covector) == 2 and A*A-10*D*x*x == K
                assert dot(roots['old0'], covector) == -beta
                assert all(dot(roots[f'old{j+1}'], covector) == -eta[j] for j in range(8) if j != 3)
                counts['integral_norm_two_roots'] += 1
                products, invalid, categories = {}, [], Counter()
                for name, old in roots.items():
                    value = dot(old, covector)
                    kind = pair_kind(value, norms[name], phi)
                    products[name] = [value.encoded(), kind]
                    categories[kind] += 1
                    if kind.startswith('reject'):
                        invalid.append(name)
                assert sum(categories.values()) == 27
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
                reasons = ([] if cut.sign() > 0 else ['does_not_cut_old4_ray'])
                reasons += ([] if not invalid else ['incompatible_native_pairs'])
                record = {'profile': ident, 'x_phi_coefficients': [a, b], 'A': A.encoded(),
                          'T': T.encoded(), 'y': y.encoded(), 'ray_product': cut.encoded(),
                          'old_center_product': cp.encoded(), 'old_center_strict': cp.sign() < 0,
                          'new_center_epsilon': str(epsilon), 'new_center_square': center_square.encoded(),
                          'invalid_walls': invalid, 'pair_categories': dict(categories), 'rejected': reasons}
                retained_bytes += len(json.dumps(record, sort_keys=True).encode())
                assert retained_bytes < 512*1024, 'Bounded retained-record limit reached'
                solutions.append(record)
                if not reasons:
                    tau = D*x/A
                    endpoint = p[:2]+[tau*row[3] for row in ci]
                    endpoint_square = dot(endpoint, mv(G, endpoint))
                    assert tau.sign() > 0 and (parameter-tau).sign() > 0
                    assert endpoint_square == -2*D*K/(A*A) and endpoint_square.sign() < 0
                    assert dot(p, mv(G, endpoint)).sign() < 0
                    endpoint_products = {n: dot(endpoint, mv(G, v)) for n, v in roots.items()}
                    assert all(v.sign() <= 0 for v in endpoint_products.values())
                    white = {'old0'} | {f'old{j+1}' for j in range(8) if j != 3}
                    assert {n for n, v in endpoint_products.items() if v == 0} == white
                    assert dot(endpoint, covector) == 0
                    survivor = record | {'coordinates': [v.encoded() for v in root],
                        'covector': [v.encoded() for v in covector], 'wall_products': products,
                        'strict_center': [v.encoded() for v in center],
                        'endpoint_parameter': tau.encoded(), 'endpoint_square': endpoint_square.encoded(),
                        'endpoint': [v.encoded() for v in endpoint],
                        'cap_pair_kinds': {n: products[n][1] for n in ('r', 'w')}}
                    retained_bytes += len(json.dumps(survivor, sort_keys=True).encode())
                    assert retained_bytes < 512*1024, 'Bounded retained-record limit reached'
                    survivors.append(survivor)
    assert counts['profile_x_square_tests'] == 18*len(xs)
    return {'status': 'PASS_BOUNDED_OMEGA3_FINITE_PROFILE_SEARCH',
        'source_sha256': digest(Path(__file__)), 'sources': PINS,
        'domain': {'x_positive_at_most': 600, 'absolute_conjugate_x_at_most': 2,
                   'b_box_inclusive': [0, 300],
                   'a_box': 'm-2<=a<=m+3; m=(isqrt(5*b*b)-b)//2',
                   'equation': 'A^2-10D*x^2=K; T=(A-R)/20; y=(beta-phi*x)/2 integral'},
        'stage_counts': dict(counts), 'profiles': summaries,
        'integral_norm_two_solutions': solutions, 'survivors': survivors,
        'scope': 'Only finite short contact profiles for the fixed old4 white star, '
                 'at most one 1/phi contact per W component and Schur q<2, in the stated x box. '
                 'Both A signs and all27 native wall products are checked for every integral root. '
                 'The old center is never a rejection filter. No other star, root completeness, '
                 'vertex preservation, finite volume, nonreflectivity or dimension claim.'}


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
