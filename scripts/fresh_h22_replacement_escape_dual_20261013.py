#!/usr/bin/env python3
"""Exact positive-functional obstruction and necessary nonright contacts."""
import argparse
import json
from fractions import Fraction as F
import fresh_h22_core_replacement_20261013 as S

TAG = 'fresh_h22_replacement_escape_dual_20261013'
OUT = 'results/' + TAG + '.json'


def run():
    d = S.read(S.OUT)
    for path, digest in d['sources'].items():
        assert S.sha(path) == digest
    Q, W, C, h, M, R = S.setup()
    p = list(map(F, d['projection_p']))
    delta = list(map(F, d['cusp_ray']))
    I = d['twenty_retained_incident']
    E = list(map(F, d['strict_escape']['coordinates']))
    assert S.pair(M, E, E) == 1586
    assert S.pair(M, E, h + [0]) == -814
    assert all(S.pair(M, E, R[i]) < 0 for i in d['retained_base_indices'])
    basis = [p + [0], delta + [0], [F(0)]*22 + [F(1)]]
    witnesses = [5, 188, 229]
    matrix = [[S.pair(M, b, R[i]) for b in basis] for i in witnesses]
    assert matrix == [[-4, -1, 0], [-192, -2, 0], [-423, -8, 2]]
    det = 2 * (8 - 192)
    assert det == -368
    coefficients = [F(20608, 23), F(12995, 23), F(40)]
    functional = [S.pair(M, b, E) for b in basis]
    assert functional == [-128984, -2346, 80]
    assert [sum(a*row[j] for a, row in zip(coefficients, matrix))
            for j in range(3)] == functional
    assert all(a > 0 for a in coefficients)
    remainder = [E[j] - sum(a*R[i][j] for a, i in zip(coefficients, witnesses))
                 for j in range(23)]
    gram = [[S.pair(M, R[i], R[j]) for j in I] for i in I]
    pivots = S.positive_pivots(gram)
    old_coefficients = S.solve(gram, [S.pair(M, remainder, R[i]) for i in I])
    assert remainder == [sum(a*R[i][j] for a, i in zip(old_coefficients, I))
                         for j in range(23)]
    negative = [i for i, a in zip(I, old_coefficients) if a < 0]
    positive = [i for i, a in zip(I, old_coefficients) if a > 0]
    assert negative == [7, 14, 16, 17, 21, 23, 25, 28, 29, 30, 31, 48, 50, 63]
    assert positive == [2, 22, 33, 40, 45, 168]
    assert len(negative) + len(positive) == 20
    assert all(S.pair(M, b, R[i]) == 0 for b in basis for i in I)
    # I is PD20, and the three witness restrictions have rank three on I^perp.
    # Hence I together with these witnesses spans the full ambient space.
    # If a nonzero acute normal is right to the fourteen negative terms,
    # its pairing with E is strictly negative: equality would kill all
    # nine remaining positive terms as well and force the normal to zero.
    u = list(map(F, d['old_projection_numerator_u']))
    U = F(d['u_square'])
    a = S.pair(Q, u, E[:22])
    assert a < 0 and a*a > 2*E[-1]*E[-1]*U
    # Two sharper identities, proposed by numerical LP but checked here
    # as literal rational vector identities without any LP dependency.
    single_negative = {
        7: {4:'77/2',9:'39/2',11:'123/4',14:'77/2',17:'125/2',
            18:'9/2',21:'203/4',22:'203/4',23:'203/4',24:'123/4',
            25:'123/4',26:'39/2',28:'45/2',34:'123/4',39:'125/2',
            40:'125/2',48:'38',49:'43/4',244:'40',7:'-161/2'},
        168: {6:'69/4',14:'57/8',17:'589/8',21:'69/4',23:'337/4',
              24:'64',25:'75/4',26:'45/8',28:'71/2',29:'41/8',
              39:'65/8',40:'1',48:'45/8',49:'207/8',50:'743/8',
              143:'443/16',188:'57/16',200:'131/4',327:'3',
              428:'219/8',490:'77/8',168:'-50'}
    }
    for target, recipe in single_negative.items():
        cs = {i:F(c) for i,c in recipe.items()}
        assert set(cs) <= set(d['retained_base_indices'])
        assert [i for i,c in cs.items() if c < 0] == [target]
        assert all(c for c in cs.values())
        assert [sum(c*R[i][j] for i,c in cs.items()) for j in range(23)] == E
    return {
        'status': 'PASS_EXACT_ALL_REAL_PLANE_OBSTRUCTION_AND_14_CONTACT_GATE',
        'sources': {path: S.sha(path) for path in
                    (S.OUT, 'scripts/' + S.PRE + '.py', 'scripts/' + TAG + '.py')},
        'escape': list(map(str, E)),
        'escape_square': '1586',
        'basis': [list(map(str, b)) for b in basis],
        'witness_walls': witnesses,
        'witness_restrictions': [list(map(str, row)) for row in matrix],
        'restriction_determinant': str(det),
        'positive_witness_coefficients': list(map(str, coefficients)),
        'escape_restriction': list(map(str, functional)),
        'old_incident': I,
        'old_incident_coefficients': list(map(str, old_coefficients)),
        'positive_gram_pivots': list(map(str, pivots)),
        'necessary_nonright_contact_set': negative,
        'other_incident_walls': positive,
        'single_negative_identities': {
            str(target): {str(i):c for i,c in recipe.items()}
            for target,recipe in single_negative.items()},
        'mandatory_strictly_negative_pairs_for_any_escape_cutter': [7,168],
        'scope': 'Any nonzero real normal acute to the retained 811 parent '
                 'normals and perpendicular to the fourteen listed old walls '
                 'has strictly negative pairing with the saved escape. In '
                 'particular the entire three-dimensional complement of all '
                 'twenty incident old walls is excluded as a gap-cutting '
                 'family, with arbitrary heights, norms and Coxeter labels. '
                 'Finite additions from these restricted families cannot '
                 'give finite volume. General nonright additions and other '
                 'wall deletions remain outside the obstruction.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.verify:
        assert result == S.read(OUT)
    else:
        (S.P / OUT).write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])
