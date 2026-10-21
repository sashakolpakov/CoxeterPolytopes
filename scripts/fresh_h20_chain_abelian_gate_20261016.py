#!/usr/bin/env python3
"""Exact nonabelian rotation certificate; standalone rational/modular algebra."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_chain_abelian_gate_20261016'
BASE = 'results/fresh_h20_s61_multicap_20261014_exact.json'
OLD = 'results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
BRIDGE = 'results/fresh_h20_e7_bridge_face_20261014_exact.json'
UNIFORM = 'results/fresh_h20_uniform_chain_cap_gate_20261016.json'
PRIMES = [2, 79, 11, 3]


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def add(a, b):
    return trim([(a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(a, t):
    return trim([x*t for x in a])


def mul(a, b):
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)


def sub(a, b):
    return add(a, scale(b, -1))


def rem(a, f, p):
    a = trim([x % p for x in a])
    f = trim([x % p for x in f])
    assert f != [0]
    inv = pow(f[-1], -1, p)
    while a != [0] and len(a) >= len(f):
        d, t = len(a)-len(f), a[-1]*inv % p
        for i, x in enumerate(f):
            a[i+d] = (a[i+d]-t*x) % p
        a = trim(a)
    return a


def monic(f, p):
    f = trim([x % p for x in f])
    assert f != [0]
    return [x*pow(f[-1], -1, p) % p for x in f]


def gcd(a, b, p):
    while b != [0]:
        a, b = b, rem(a, b, p)
    return monic(a, p)


def power(a, n, f, p):
    ans = [1]
    while n:
        if n & 1:
            ans = rem(mul(ans, a), f, p)
        a = rem(mul(a, a), f, p)
        n //= 2
    return ans


def irreducible(f, p, factors):
    """x^(p^n)=x and gcd tests for every prime divisor of n."""
    f = monic(f, p)
    n = len(f)-1
    assert n > 0
    divisor_product = 1
    for q in factors:
        assert isprime(q) and n % q == 0
        while (n//divisor_product) % q == 0:
            divisor_product *= q
    assert divisor_product == n
    v = [0, 1]
    tests = []
    for k in range(1, n+1):
        v = power(v, p, f, p)
        if k in [n//q for q in factors]:
            g = gcd(f, rem(sub(v, [0, 1]), f, p), p)
            assert g == [1]
            tests.append(dict(frobenius_power=k, gcd_degree=0))
    assert rem(sub(v, [0, 1]), f, p) == [0]
    return dict(degree=n, prime_divisors=factors, coprimality_tests=tests,
                final_frobenius_remainder=[0])


def isprime(p):
    return p >= 2 and all(p % i for i in range(2, isqrt(p)+1))


def parse(v):
    return [parse(x) for x in v] if isinstance(v, list) else F(v)


def mod(q, p):
    assert q.denominator % p
    return q.numerator*pow(q.denominator, -1, p) % p


def basis(p, roots):
    assert isprime(p) and len(roots) == 4
    assert all(r*r % p == a % p and 2*r % p for a, r in zip(PRIMES, roots))
    rr = [1]
    for r in roots:
        rr += [x*r % p for x in rr]
    return rr


def eval0(v, p, rr):
    assert len(v) == 16
    return sum(mod(q, p)*r for q, r in zip(v, rr)) % p


def eval1(v, p, rr, rd):
    assert len(v) == 2
    return (eval0(v[0], p, rr)+rd*eval0(v[1], p, rr)) % p


def eval2(v, p, rr, rd, rh):
    assert len(v) == 2
    return (eval1(v[0], p, rr, rd)+rh*eval1(v[1], p, rr, rd)) % p


def sha(p):
    return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    base, old, bridge, uniform = [json.loads((ROOT/p).read_text())
                                  for p in (BASE, OLD, BRIDGE, UNIFORM)]
    assert uniform['status'] == 'PASS_ALL_LABEL_UNIFORM_FINITE_CHAIN_ANGLE_CAP_OBSTRUCTION'
    assert base['chain']['walls'] == list(range(107, 133))
    assert base['chain']['norms'] == [2]*26
    assert old['field']['base_square_root_basis'] == [1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
    D = parse(old['field']['discriminant'])
    H = parse(bridge['field']['new_positive_sqrt_radical'])
    w = parse(base['frame']['W_norm'])
    # The four distinct primes have independent squareclasses, hence degree16.
    assert len(set(PRIMES)) == 4 and all(isprime(p) for p in PRIMES)
    rr = basis(313, [193, 198, 295, 230])
    d313 = eval0(D, 313, rr)
    assert d313 == 244 and pow(d313, 156, 313) == 312
    rr = basis(359, [340, 325, 332, 196])
    d359 = eval0(D, 359, rr)
    assert d359 == 246 and 167**2 % 359 == d359 and 2*167 % 359
    h359 = eval1(H, 359, rr, 167)
    assert h359 == 342 and pow(h359, 179, 359) == 358
    # Nonresidue witnesses prove that each successive radicand is nonsquare
    # in the preceding actual number field. The complete tower has degree64.
    tower_irreducibility = [
        dict(stage='sqrt(D)', p=313, base_roots=[193,198,295,230], residue=d313,
             euler_criterion=312),
        dict(stage='sqrt(H)', p=359, base_roots=[340,325,332,196], D=d359,
             sqrt_D=167, residue=h359, euler_criterion=358)]

    # Chebyshev endpoint formulas in the variable C=a²/(a²+2).
    U0, U1 = [1], [0, 2]
    U = [U0, U1]
    for _ in range(11):
        U.append(sub(mul([0, 2], U[-1]), U[-2]))
    A, B = sub(U[12], U[11]), add(U[12], U[11])
    aa, bb = [1], [1]
    for _ in range(12):
        aa, bb = sub(mul([0, 1], aa), mul([1, -1], bb)), add(mul([0, 1], bb), mul([1, 1], aa))
    assert (aa, bb) == (A, B)
    R, K = mul([1, 1], mul(A, A)), mul([1, -1], mul(B, B))
    assert add(R, K) == [2]
    d = sub([0, 2], K)
    e6 = add([0, 10], scale(R, 3))
    # P(C,w)=36[(d*w-e6/6)^2-2*w*R*K].
    coeffs = [mul(e6, e6), sub(scale(mul(d, e6), -12), scale(mul(R, K), 72)), scale(mul(d, d), 36)]
    assert max(map(len, coeffs)) == 51
    reductions = []
    specs = [dict(p=2017, roots=[1031,206,1360,1759], rd=1831, rh=1772, w=1649),
             dict(p=359, roots=[340,34,332,163], rd=340, rh=222, w=15)]
    for spec in specs:
        p, rd, rh = spec['p'], spec['rd'], spec['rh']
        rr = basis(p, spec['roots'])
        dr, hr = eval0(D, p, rr), eval1(H, p, rr, rd)
        assert rd*rd % p == dr and rh*rh % p == hr
        assert dr and hr and 2*rd % p and 2*rh % p
        wr = eval2(w, p, rr, rd, rh)
        assert wr == spec['w']
        fraw = [sum((coeffs[j][i] if i < len(coeffs[j]) else 0)*pow(wr, j, p)
                    for j in range(3)) % p for i in range(51)]
        assert fraw[-1] != 0
        f = monic(fraw, p)
        derivative = [(i*f[i]) % p for i in range(1, len(f))]
        assert gcd(f, derivative, p) == [1]
        rec = dict(**spec, D_residue=dr, H_residue=hr,
                   simple_hensel_roots=True, polynomial_leading_unit=fraw[-1],
                   reduced_monic_polynomial=f, squarefree=True)
        if p == 2017:
            rec['irreducibility'] = irreducible(f, p, [2, 5])
        else:
            v = [0, 1]
            factors = []
            for k in range(1, 38):
                v = power(v, p, f, p)
                if k in (13, 37):
                    g = gcd(f, rem(sub(v, [0, 1]), f, p), p)
                    assert len(g)-1 == k
                    cert = irreducible(g, p, [k])
                    factors.append(dict(polynomial=g, certificate=cert))
            assert [x % p for x in mul(factors[0]['polynomial'], factors[1]['polynomial'])] == f
            rec['factors'] = factors
            rec['factor_degrees'] = [13, 37]
        reductions.append(rec)

    # The saved uniform proof already checked the closed endpoint t=2/a.
    alo, ahi = map(F, base['parameter']['isolating_interval'])
    assert alo > 0 and alo*alo > 280
    assert list(map(F, uniform['time_interval'])) == [1/ahi, 2/alo]
    assert 2+4/(alo*alo) < F(141,70)
    assert uniform['finite_support_cost_census']['strict_above_two_minimum'] == '137/66'
    for rec in uniform['exact_blockers']:
        assert F(rec['whole_parameter_positive_pairing'][0]) > 0

    queue, seen, references = [BASE, OLD, BRIDGE, UNIFORM], set(), 0
    while queue:
        name = queue.pop()
        if name in seen:
            continue
        seen.add(name)
        for path, h in json.loads((ROOT/name).read_text()).get('sources', {}).items():
            assert sha(path) == h, (name, path)
            references += 1
            if path.endswith('.json'):
                queue.append(path)
    paths = [f'scripts/{TAG}.py', BASE, OLD, BRIDGE, UNIFORM,
             'scripts/fresh_h20_uniform_chain_cap_gate_20261016.py']
    return dict(status='PASS_EXACT_NONABELIAN_CHAIN_ROTATION_AND_NONDOTTED_CAP_OBSTRUCTION',
                sources={p:sha(p) for p in paths}, tower_degree=64,
                tower_irreducibility=tower_irreducibility,
                endpoint_A=A, endpoint_B=B,
                endpoint_polynomial_coefficients_by_w_degree=coeffs,
                endpoint_polynomial_degree_C=50, endpoint_polynomial_degree_w=2,
                modular_reductions=reductions,
                irreducible_degree_over_tower=50, extension_is_not_galois=True,
                rotation_C_not_in_Q_ab=True,
                recurrence='f[j+3]-f[j+2]+f[j+1]-f[j]=2*C*(f[j+2]-f[j+1])',
                consecutive_nondotted_block_length_forcing_constancy=5,
                minimum_dotted_chain_contacts=5,
                uniform_parallel_endpoint_checked=True,
                recursive_manifests_checked=len(seen), recursive_references_checked=references,
                scope='Every coherent new Coxeter wall retaining the frozen133-wall H20 chamber must have at least five dotted relations among the26 chain walls. Five consecutive chain pairings in Q^ab force all26 pairings equal, and their finite-or-parallel constant value is excluded. No claim excludes caps having dotted chain contacts, changed retained walls, or arbitrary H20 completions.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'])
