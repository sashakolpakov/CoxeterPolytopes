#!/usr/bin/env python3
"""Review thin alphabets and calibrate the complete two-block census at H=30.

This is a candidate-enumeration calibration, not an independent geometry replay.
The coefficient-box oracle does not use thin_alphabet's interval pruning.  The
full Host oracle independently joins all four blocks and the scalar coordinate.
Both reuse the repository's exact cubic arithmetic.
"""
from fractions import Fraction as F
from itertools import product, permutations
import json
from pathlib import Path
import time

from fresh_cubic7_dihedral_h9_20261013 import (
    Host, K, T, C, D, BETA, BLOCK, Z, O, alphabet, enc, pair, mv, ray_key,
)
from fresh_cubic7_target_alphabet_20261013 import (
    thin_alphabet, DUAL, BOUNDS, INTERVALS,
)

ROOT = Path(__file__).resolve().parents[1]


def ceil(x):
    return -((-x.numerator)//x.denominator)


def coefficient_oracle(caps):
    """Exhaust every coefficient in a proved box, without interval pruning."""
    if any(c.sign(e) < 0 for e, c in enumerate(caps)):
        return ()
    ceilings = []
    for e, c in enumerate(caps):
        m = 0
        while (K(m*m)-c).sign(e) < 0:
            m += 1
        ceilings.append(m)
    bounds = tuple(ceil(sum(b*m for b, m in zip(row, ceilings))) for row in BOUNDS)
    answer = []
    for a, b, c in product(*(range(-n, n+1) for n in bounds)):
        x = K(a, b, c)
        if all((cap-x*x).sign(e) >= 0 for e, cap in enumerate(caps)):
            answer.append(x)
    return tuple(sorted(answer))


def target_roots(height):
    """Direct two-block reconstruction, without calling its search producer."""
    roots = set()
    records = []
    alpha = 2*T
    for q in (K(2), 2*BETA):
        times = tuple(x for x in thin_alphabet((height*q/alpha, -q/alpha, -q/alpha))
                      if x.sign() > 0 and (2*alpha*x/q).is_integral)
        weights = tuple(w for w in thin_alphabet((2*q*(height+1), 2*q, 2*q))
                        if w.sign() >= 0 and (2*w/q).is_integral)
        blocks = []
        for u in weights:
            for w in weights:
                if (w-u).sign() < 0:
                    continue
                y = ((2*u+C*w)/BETA, (C*u+2*w)/BETA)
                if not all(a.is_integral for a in y):
                    continue
                energy = pair(y, BLOCK, y)
                if (q*(height+1)-energy).sign() < 0:
                    continue
                if any((q-energy).sign(e) < 0 for e in (1, 2)):
                    continue
                blocks.append((energy, y))
        norm_roots = set()
        # A direct Cartesian product, independent of the producer's hash join.
        for x, (e1, y1), (e2, y2) in product(times, blocks, blocks):
            if e1+e2 != q+alpha*x*x:
                continue
            for i, j in permutations(range(4), 2):
                r = [x]+[Z]*9
                r[1+2*i:3+2*i] = y1
                r[1+2*j:3+2*j] = y2
                norm_roots.add(tuple(r))
        roots.update(norm_roots)
        records.append({'norm': q.to_json(), 'controls': len(times),
                        'weights': len(weights), 'blocks': len(blocks),
                        'distinct_roots': len(norm_roots)})
    return roots, records


def main():
    started = time.monotonic()
    basis = (O, T, T*T)
    trace_duality = [[(d*b).trace() for b in basis] for d in DUAL]
    assert trace_duality == [[F(i == j) for j in range(3)] for i in range(3)]
    assert all((T-lo).sign(e) > 0 and (hi-T).sign(e) > 0
               for e, (lo, hi) in enumerate(INTERVALS))
    assert all((K(bound)-d).sign(e) > 0 and (K(bound)+d).sign(e) > 0
               for d, row in zip(DUAL, BOUNDS) for e, bound in enumerate(row))
    assert BETA == 2-T and BETA.norm() == 7 and C.norm() == -1
    assert all(BETA.sign(e) > 0 for e in range(3))
    # In F_2[theta]/(theta^3+theta^2+1), x^2=0 only for x=0.
    residue_squares = []
    for coeffs in product(range(2), repeat=3):
        square = K(*coeffs)**2
        mod = tuple(int(c) % 2 for c in (square.a, square.b, square.c))
        assert bool(any(mod)) == bool(any(coeffs))
        residue_squares.append({'element': list(coeffs), 'square': list(mod)})
    h = 30
    caps_cases = [(Z, Z, Z), (-O, O, O)]
    for q in (K(2), 2*BETA):
        caps_cases.extend(((h*q/(2*T), -q/(2*T), -q/(2*T)),
                           (2*q*(h+1), 2*q, 2*q)))
    comparisons = []
    for caps in caps_cases:
        thin = thin_alphabet(caps)
        brute = coefficient_oracle(caps)
        assert thin == brute
        comparisons.append({'caps': enc(caps), 'count': len(thin)})
    host = Host(2*T, 4, True)
    for k in range(4):
        first, second = host.initial()[2*k:2*k+2]
        # The inequalities are -u<=0 and D(u-w)<=0, with D>0.
        assert D.sign() > 0
        assert first[1+2*k] == -O
        assert second[1+2*k] == D and second[2+2*k] == -D
    targeted, census = target_roots(h)
    full = set()
    full_controls = []
    for q in host.norms:
        controls = tuple(x for x in alphabet((h*q/(2*T), -q/(2*T), -q/(2*T)))
                         if x.sign() > 0 and (4*T*x/q).is_integral)
        count = 0
        for x in controls:
            all_roots, _ = host.candidates(x, q)
            for r in all_roots:
                supports = sum(any(r[1+2*i:3+2*i]) for i in range(4))
                if r[-1] == Z and supports <= 2:
                    assert q in (K(2), 2*BETA)
                    full.add(r)
                    count += 1
        full_controls.append({'norm': q.to_json(), 'controls': len(controls),
                              'scalar_zero_at_most_two_blocks': count})
    assert targeted == full
    support_counts = {str(k): sum(sum(any(r[1+2*i:3+2*i]) for i in range(4)) == k
                                 for r in full) for k in range(3)}
    for r in full:
        q = pair(r, host.gram, r)
        assert all((2*a/q).is_integral for a in mv(host.gram, r))
        assert r[-1] == Z
    output = {
        'status': 'PASS', 'height': h,
        'trace_duality': [[str(x) for x in row] for row in trace_duality],
        'mod2_square_injectivity': residue_squares,
        'thin_vs_unpruned_coefficient_box': comparisons,
        'two_block_census': census, 'full_host_controls': full_controls,
        'matched_distinct_roots': len(full), 'support_counts': support_counts,
        'matched_roots': [enc(r) for r in sorted(full, key=ray_key)],
        'scope': 'Static completeness proof plus independent-enumeration H30 calibration; '
                 'project cubic arithmetic is shared. No H1000 dictionary rerun or geometry replay.',
        'seconds': round(time.monotonic()-started, 3),
    }
    path = ROOT/'results/fresh_cubic7_two_block_completeness_audit_20261013.json'
    path.write_text(json.dumps(output, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: output[k] for k in ('status', 'matched_distinct_roots', 'support_counts', 'seconds')}))


if __name__ == '__main__':
    main()
