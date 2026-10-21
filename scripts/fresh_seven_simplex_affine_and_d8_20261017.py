#!/usr/bin/env python3
"""Exact four hard seven-simplex cases; arbitrary outside whites are in the proof.

Uses the frozen standalone Fraction Q(sqrt2) checker, not producer SymPy
arithmetic. The classical attachment census supplies the complete profiles.
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

from fresh_finite_edge_cusp_tradeoff12_20261017 import (
    Q, gram, inertia, inverse, kernel1, matrix_data, mm, sub, transpose,
)

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_seven_simplex_affine_and_d8_20261017'
CENSUS = 'fresh_high_rank_quasilanner_census_20261017'


def cost_matrices(B, profiles):
    out = []
    for b in range(len(B)):
        ids = [i for i in range(len(B)) if i != b]
        H = sub(B, ids)
        assert inertia(H) == [len(B)-1, 0, 0]
        J = inverse(H)
        assert all(x.sign() >= 0 for row in J for x in row)
        P = [[profile[i] for profile in profiles] for i in ids]
        out.append(mm(mm(transpose(P), J), P))
    return out


def residual(K, chosen, native=Q()):
    H = sub(K, chosen)
    return [[Q(2 if i == j else 0) - H[i][j] -
             (native if i != j else Q()) for j in range(len(H))]
            for i in range(len(H))]


def det2(H):
    return H[0][0]*H[1][1]-H[0][1]*H[1][0]


def finite_black_rank9(diagrams):
    rows = []
    for L in diagrams:
        if L['order'] != 9:
            continue
        G = gram(9, L['edges'])
        for w in range(9):
            if w in L['ideal_deletions']:
                continue
            black = [i for i in range(9) if i != w]
            assert inertia(sub(G, black)) == [8, 0, 0]
            ideals = L['ideal_deletions']
            rec = dict(diagram=L['index'], internal_white=w)
            if len(ideals) >= 2:
                rec.update(method='outside_whites_orthogonal', maximum_ideals=len(ideals))
            else:
                a = ideals[0]
                thresholds = []
                for b in black:
                    if b == a:
                        continue
                    ids = [i for i in range(9) if i != b]
                    J = inverse(sub(G, ids))
                    t = (1/J[ids.index(a)][ids.index(a)]).rational()
                    thresholds.append([b, str(t)])
                low = min(F(t) for _, t in thresholds)
                ties = [b for b, t in thresholds if F(t) == low]
                rec.update(method='one_black_diagonal_decrease', minimum=str(low),
                           thresholds=thresholds, ties=ties, maximum_ideals=1+len(ties))
            assert rec['maximum_ideals'] <= 3
            rows.append(rec)
    assert len(rows) == 31
    return rows


def affine_cases(source):
    out = []
    for index, source_name, support in (
        (20, 'B7', [(0,), (1,)]),
        (21, 'D7', [(0,), (1,), (6,), (7,)]),
        (22, 'E7', [(4,), (7,), (4, 7)]),
    ):
        L = source['diagrams'][index]
        G = gram(9, L['edges'])
        B = sub(G, list(range(8)))
        assert inertia(B) == [7, 0, 1]
        kernel = kernel1(B)
        profiles = [[Q(i in supp) for i in range(8)] for supp in support]
        provenances = [p for p in source['attachment_provenances']
                       if p['affine_source'] == source_name]
        assert sorted(p['contacts'] for p in provenances) == sorted(
            [3 if i in supp else 2 for i in range(8)] for supp in support)
        # Every provenance map must recover this same literal affine source.
        for p in provenances:
            d = source['diagrams'][p['diagram']]
            D = gram(d['order'], d['edges'])
            mapping = dict(p['isomorphism'])
            assert sub(D, [mapping[i] for i in range(8)]) == B
        K = cost_matrices(B, profiles)
        singles = []
        duplicates = []
        for i in range(len(profiles)):
            costs = [k[i][i].rational() for k in K]
            high = max(costs)
            assert high > 1
            b = costs.index(high)
            repeated = [[Q(2*(r == c))-K[b][i][i] for c in range(2)] for r in range(2)]
            assert det2(repeated).sign() < 0
            duplicates.append(dict(profile=i, black=b, determinant=det2(repeated).data()))
            low = 2-high
            assert low >= 0
            ties = [b for b, c in enumerate(costs) if c == high]
            singles.append(dict(profile=i, costs=list(map(str, costs)),
                                outside_diagonal_maximum=str(low),
                                possible_ideal_blacks=ties, maximum_ideals=len(ties)))
        pairs = []
        for i, j in itertools.combinations(range(len(profiles)), 2):
            H0 = [residual(k, [i, j]) for k in K]
            bad = [b for b, h in enumerate(H0) if inertia(h)[1]]
            if bad:
                b = bad[0]
                assert det2(H0[b]).sign() < 0
                pairs.append(dict(profiles=[i, j], allowed_native_magnitudes=[],
                                  exclusion_black=b, determinant_at_right=det2(H0[b]).data()))
                continue
            allowed = [0]
            H1 = [residual(k, [i, j], Q(1)) for k in K]
            if all(inertia(h)[1] == 0 for h in H1):
                allowed.append(1)
                witnesses = [b for b, h in enumerate(H1) if not det2(h)]
                assert witnesses
                witness = witnesses[0]
                cap = 1
            else:
                witnesses = [b for b, h in enumerate(H0) if not det2(h)]
                assert witnesses
                witness = witnesses[0]
                cap = 0
            assert K[witness][i][j].sign() > 0
            # Increasing the nonnegative native contact magnitude decreases
            # this determinant strictly, so these caps cover every Coxeter label.
            rec = dict(profiles=[i, j], allowed_native_magnitudes=allowed,
                       all_label_cap_black=witness, all_label_cap=str(cap), statuses=[])
            for native in allowed:
                H = H0 if native == 0 else H1
                ideals = [b for b, h in enumerate(H) if not det2(h)]
                if ideals:
                    b = ideals[0]
                    ker = kernel1(H[b])
                    method = 'positive_null_vector_forces_entrywise_nonnegative_D_zero'
                    status = dict(native=native, original_ideal_blacks=ideals,
                                  positive_null_at=b, positive_null=[x.data() for x in ker],
                                  method=method, maximum_ideals=len(ideals))
                else:
                    assert index == 22 and [i, j] == [0, 1] and native == 0
                    status = dict(native=0, original_ideal_blacks=[],
                                  method='E7_right_pair_indirect_white_inequality',
                                  maximum_ideals=2, possible_ideal_blacks=[4, 7])
                assert status['maximum_ideals'] <= 3
                rec['statuses'].append(status)
            pairs.append(rec)
        # Every three distinct profiles contain an excluded pair.
        permitted = {tuple(p['profiles']) for p in pairs if p['allowed_native_magnitudes']}
        for triple in itertools.combinations(range(len(profiles)), 3):
            assert any(pair not in permitted for pair in itertools.combinations(triple, 2))
        out.append(dict(diagram=index, affine_source=source_name, black_gram=matrix_data(B),
                        positive_kernel=[x.data() for x in kernel], profile_supports=support,
                        cost_matrices=[matrix_data(k) for k in K], single_profiles=singles,
                        repeated_profile_exclusions=duplicates, pair_profiles=pairs,
                        maximum_ideals=max([x['maximum_ideals'] for x in singles] +
                            [r['maximum_ideals'] for p in pairs for r in p.get('statuses', [])])))
    assert [r['maximum_ideals'] for r in out] == [1, 2, 3]
    return out


def e7_right_pair(diagrams):
    B = sub(gram(9, diagrams[22]['edges']), list(range(8)))
    profiles = [[Q(i == node) for i in range(8)] for node in (4, 7)]
    K = cost_matrices(B, profiles)
    H = [residual(k, [0, 1]) for k in K]
    assert H[4] == [[Q(2), Q()], [Q(), Q(F(1, 2))]]
    assert H[7] == [[Q(F(1, 2)), Q()], [Q(), Q(2)]]
    # H4-D and H7-D PSD imply x,y<=1/2. D PSD implies z<=sqrt(xy)<=1/2.
    strict = []
    for b in (0, 2, 3, 5, 6):
        assert not H[b][0][1]
        lower = (H[b][0][0]-Q(F(1, 2)))*(H[b][1][1]-Q(F(1, 2)))-Q(F(1, 4))
        assert lower.sign() > 0
        strict.append(dict(black=b, determinant_lower_bound=lower.data()))
    assert H[1] == [[Q(F(9, 8)), Q(F(-1, 8))], [Q(F(-1, 8)), Q(F(9, 8))]]
    # Its diagonal entries after D are >=5/8 and |offdiag|<=5/8.
    # A zero determinant can only have x=y=z=1/2, which violates H4-D PSD.
    half = [[Q(F(1, 2))]*2 for _ in range(2)]
    endpoint = [[H[4][i][j]-half[i][j] for j in range(2)] for i in range(2)]
    assert det2(endpoint) == Q(F(-1, 4))
    return dict(residuals=[matrix_data(h) for h in H],
                outside_D_bounds=dict(x='[0,1/2]', y='[0,1/2]', z='[0,sqrt(xy)]'),
                strict_determinant_bounds=strict, black1_equality_requires='x=y=z=1/2',
                equality_black4_determinant=det2(endpoint).data(), possible_ideal_blacks=[4, 7])


def d8_case(diagrams):
    L = gram(10, diagrams[23]['edges'])
    black = list(range(8))
    assert inertia(sub(L, black)) == [8, 0, 0]
    standard_to_L = [2, 0, 1, 9, 3, 4, 5, 6, 7]
    E = sub(L, standard_to_L)
    standard = sub(gram(10, diagrams[25]['edges']), list(range(9)))
    assert E == standard and inertia(E) == [8, 0, 1]
    assert [L[8][i] for i in standard_to_L] == [Q()]*8+[-Q(0, 1)]
    profiles = [[(Q(1) if kind == 'A' else Q(0, 1)) if i == 8 else Q()
                 for i in range(9)] if kind in ('A', 'B') else
                [Q(i == 7) for i in range(9)] for kind in ('A', 'B', 'C')]
    K = cost_matrices(E, profiles)
    exclusions = []
    for p in range(3):
        paircost = [[K[1][i][j] for j in (1, p)] for i in (1, p)]
        H = [[Q(2*(i == j))-paircost[i][j] for j in range(2)] for i in range(2)]
        determinant = det2(H)
        assert determinant.sign() < 0
        exclusions.append(dict(additional_profile='ABC'[p], determinant_at_right=determinant.data()))
    assert [r['determinant_at_right'] for r in exclusions] == [
        Q(F(-4, 3)).data(), Q(F(-28, 9)).data(), Q(F(-10, 9)).data()]
    costs = [[j, str(K[j][1][1].rational())] for j in range(9) if j != 3]
    high = max(F(c) for _, c in costs)
    assert high == F(16, 9) and [j for j, c in costs if F(c) == high] == [1]
    return dict(black_type='finite_D8', original_internal_whites=[8, 9],
                E8_standard_to_L=standard_to_L, E8_internal_white_standard=3,
                existing_direct_profile='B', delete1_cost_gram=matrix_data(K[1]),
                excludes_every_second_direct_white=exclusions, existing_profile_costs=costs,
                unique_maximum_black=1, maximum_ideals_if_all_B_plus_u_PSD=1,
                maximum_ideals_including_rank9_alternative=3)


def run():
    source = json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    easy = finite_black_rank9(source['diagrams'])
    affine = affine_cases(source)
    right_pair = e7_right_pair(source['diagrams'])
    d8 = d8_case(source['diagrams'])
    paths = [f'scripts/{TAG}.py', f'scripts/{CENSUS}.py', f'results/{CENSUS}.json',
             'scripts/fresh_finite_edge_cusp_tradeoff12_20261017.py',
             'results/fresh_finite_edge_cusp_tradeoff12_20261017.json',
             'results/fresh_finite_edge_cusp_tradeoff12_20261017.md',
             'results/fresh_multi_ideal_triangle_scope_20261017.md']
    return dict(status='PASS_FOUR_HARD_SEVEN_SIMPLEX_CASES_AT_MOST_THREE_IDEAL',
                sources={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                finite_black_rank9_alternatives=easy, affine_black_cases=affine,
                E7_right_pair_outside_white_certificate=right_pair, finite_D8_case=d8,
                scope='Four hard minimal-star cases only; no global seven-simplex theorem or ambient dimension bound by this package alone.',
                arbitrary_white_reduction='For affine B every directly touching white has a classified profile. Repeated profiles and all-label pair inequalities leave only the listed U. Z is orthogonal to B, and D=R^tG_Z^-1R is PSD and entrywise nonnegative. Single and positive-null residuals, plus the explicit E7 pair inequality, control every black vertex. The D8 case first reduces any indefinite B+u to a separately checked finite-B rank9 case; otherwise it uses affineE8 with internal white3 and the existing B-profile.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], 'saved_replay=', args.verify)
