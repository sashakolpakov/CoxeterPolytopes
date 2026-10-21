#!/usr/bin/env python3
"""Replay the conditional obstruction to quasi-arithmetic completion.

This verifies finite algebraic data. The implication uses Vinberg's criterion,
as explained in the accompanying report; it does not prove finite volume.
"""
import argparse
from fractions import Fraction as F
import fresh_h22_incidence_prism_20261013_verify as V

OUT = 'results/' + V.PRE + '_arithmetic_obstruction.json'


def run():
    saved = V.read(OUT)
    assert V.sha(V.OUT) == saved['source_sha256']
    M, roots, _ = V.setup()
    short = [i for i, r in enumerate(roots[:210]) if V.pair(M, r, r) == 2]
    assert len(short) == 42 and short == saved['old_short_indices']
    gram = [[V.pair(M, roots[i], roots[j]) for j in short] for i in short]
    sums = [sum(row) for row in gram]
    assert sums == [-3] * 42 == saved['unnormalized_Gram_row_sums']
    assert sum(sums) == -126 == saved['ones_quadratic_value']

    # An actual retained pair has equal norm q and pairing -q*T/2.
    recipe = V.read(V.REC)['recipes'][0]
    anchor = recipe['anchor']
    q = F(recipe['norm'])
    u, v = list(map(F, recipe['u'])), list(map(F, recipe['v']))
    assert V.pair(M, roots[anchor], roots[anchor]) == q
    assert V.pair(M, u, roots[anchor]) == 1
    assert V.pair(M, v, roots[anchor]) == 0
    assert V.pair(M, u, v) == 0
    qu, qv = V.pair(M, u, u), V.pair(M, v, v)
    assert qv > 0
    rad = (q/qv, F(0), -q*q*qu/(4*qv))
    assert V.ka(V.ks(V.TWO, q*q*qu/4), V.ks(rad, qv)) == V.ks(V.ONE, q)
    assert V.ki(rad)[0] > 0
    assert V.read(V.OUT)['new_caps'][0]['finite_contacts'][str(anchor)] == 7

    # The length-two cycle of 2*normalized-Gram equals T^2.
    x = V.TWO
    poly = V.ka(V.ka(V.km(V.km(x, x), x), V.ks(V.km(x, x), -5)),
                V.ka(V.ks(x, 6), V.ks(V.ONE, -1)))
    assert poly == V.ZERO
    assert saved['cycle_polynomial_ascending'] == [-1, 6, -5, 1]
    f = lambda t: t**3 - 5*t*t + 6*t - 1
    assert f(1) == 1 and f(-1) == -13  # All possible rational roots.
    print(saved['status'])


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true', help='Replay the frozen artifact (default).')
    p.parse_args()
    run()
