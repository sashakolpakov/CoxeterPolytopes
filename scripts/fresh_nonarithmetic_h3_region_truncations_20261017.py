#!/usr/bin/env python3
"""All-label single-wall tests for five patches of the retained S6 face.

This is a necessary-condition obstruction for three patches, not an obstruction
to coordinated multi-facet compactification. All sign tests use outward
rational intervals for the frozen exact physical normal recipes.
"""
import argparse
import itertools
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_h3_endpoint_domains_20261017 as E

P, U, I = E.P, E.U, E.I
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_h3_region_truncations_20261017'
PATCHES = [
    dict(sides=[7, 30, 49], edges=[[7, 30], [7, 49], [30, 49]],
         vertices=[[7, 30, 49]]),
    dict(sides=[7, 33, 49], edges=[[7, 33], [7, 49], [33, 49]],
         vertices=[[7, 33, 49]]),
    dict(sides=[30, 33, 49], edges=[[30, 33], [30, 49], [33, 49]],
         vertices=[[30, 33, 49]]),
    dict(sides=[8, 30, 44, 52], edges=[[8, 44], [8, 52], [30, 44], [30, 52]],
         vertices=[[8, 44, 52], [30, 44, 52]]),
    dict(sides=[30, 34, 47, 52], edges=[[30, 47], [30, 52], [34, 47], [34, 52]],
         vertices=[[30, 47, 52], [34, 47, 52]]),
]


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def profile_intersection(patch, domains):
    accepted = []
    for a2 in [(2, 2)] + [(m, 2) for m in (3, 4, 5)] + [(2, m) for m in (3, 4, 5)]:
        # Code6 denotes the full continuous tail, not the single label6.
        for side_labels in itertools.product(range(2, 7), repeat=len(patch['sides'])):
            labels = {5: a2[0], 6: a2[1], **dict(zip(patch['sides'], side_labels))}
            for edge in patch['edges']:
                p = domains[tuple(edge)]['profiles']
                ls = [labels[i] for i in [5, 6] + edge]
                if 6 in ls:
                    good = (ls.count(6) == 1 and ls.index(6) in p['isolated_A1_slots']
                            and all(x == 2 for x in ls if x != 6))
                else:
                    good = ls in [r['labels'] for r in p['finite_label_profiles']]
                if not good:
                    break
            else:
                assert sum(m == 6 for m in labels.values()) <= 1
                accepted.append(labels)
    return accepted


def run():
    endpoint = E.B.read(f'results/{E.TAG}.json')
    domains = {tuple(r['sides']): r for r in endpoint['endpoint_domains']}
    _, qi, normals, h, _, _, _, cone = P.physical_data()
    pair = lambda x, y: U.idot(x, U.imv(qi, y))
    old_hyper = sorted(tuple(i for i in ids if i not in P.S)
                       for ids in cone['complete_census']['hyperideal_rays']
                       if set(P.S) <= set(ids))
    assert old_hyper == sorted(tuple(v) for p in PATCHES for v in p['vertices'])
    edge_counts = Counter(tuple(e) for p in PATCHES for e in p['edges'])
    assert edge_counts == Counter({tuple(r['sides']): r['open_ends']
                                   for r in endpoint['endpoint_domains']})
    # The two quadrilateral patches join pairs only across raw cone edges that
    # have no geodesic portion. Their four listed boundary edges are geodesic.
    raw = {tuple(r['sides']): r for r in endpoint['edge_inventory']}
    for sides in [(44, 52), (47, 52)]:
        assert not raw[sides]['geodesic']
        assert raw[sides]['endpoint_kinds'] == ['hyperideal_rays'] * 2

    records = []
    expected = [(2, 0), (5, 0), (8, 0), (84, 4), (84, 4)]
    for number, patch in enumerate(PATCHES):
        profiles = profile_intersection(patch, domains)
        finite = sum(6 not in row.values() for row in profiles)
        tails = len(profiles) - finite
        assert (finite, tails) == expected[number]
        ids = P.S + patch['sides']
        frame = [normals[i] for i in ids]
        gram = [[pair(x, y) for y in frame] for x in frame]
        inverse, pivots = P.C53.inverse(gram)
        item = dict(patch_index=number, **patch, frame=ids,
                    full_gram=[P.enc(row) for row in gram], inverse_pivots=pivots,
                    finite_profiles=finite, high_label_tails=tails,
                    profiles=[{str(i): m for i, m in row.items()} for row in profiles])
        if number in (1, 2):
            item['status'] = 'DOMAIN_ONLY_NOT_EXCLUDED_OR_ACCEPTED'
            records.append(item)
            continue
        if len(ids) == 9:
            kernel = P.sub(h, P.combine(U.imv(inverse, [pair(h, n) for n in frame]), frame))
            kernel_square = pair(kernel, kernel)
            assert kernel_square.l > 0
            item.update(kernel=P.enc(kernel), kernel_square=kernel_square.enc())
        decisions = []
        for index, labels in enumerate(profiles):
            tail = [i for i in ids if labels.get(i) == 6]
            contact = [-U.theta(labels.get(i, 2)) if labels.get(i) != 6 else I(0) for i in ids]
            p = P.combine(U.imv(inverse, contact), frame)
            psquare = pair(p, p)
            if tail:
                assert len(ids) == 10 and len(tail) == 1
                direction = P.combine(U.imv(inverse, [I(-int(i == tail[0])) for i in ids]), frame)
                # q(p+c*direction)-2 = a*c^2+b*c+d; c>=sqrt3.
                a, b, d = pair(direction, direction), U.rnd(2 * pair(p, direction)), U.rnd(psquare - 2)
                assert a.l > 0
                if all(x.l == x.h == 0 for x in contact):
                    assert all(x.l == x.h == 0 for x in p)
                    assert b.l == b.h == 0
                else:
                    assert b.l > 0
                at_sqrt3 = U.rnd(3 * a + U.sqrt(I(3)) * b + d)
                assert at_sqrt3.l > 0
                decisions.append(dict(profile=index, status='TAIL_NORM_ALWAYS_ABOVE2',
                                      tail_wall=tail[0], continuous_domain='sqrt3<=c<2',
                                      polynomial=P.enc([a, b, d]), value_at_sqrt3=at_sqrt3.enc(),
                                      reason='a>0,b>=0 and polynomial(sqrt3)>0'))
            elif len(ids) == 10:
                difference = U.rnd(psquare - 2)
                assert difference.h < 0 or difference.l > 0
                decisions.append(dict(profile=index, status='UNIQUE_NORMAL_HAS_WRONG_NORM',
                                      norm_difference=difference.enc()))
            else:
                delta = U.rnd((2 - psquare) / kernel_square)
                assert delta.l > 0
                scale = U.sqrt(delta)
                for sign in (-1, 1):
                    r = [U.rnd(x + sign * scale * y) for x, y in zip(p, kernel)]
                    if sign == 1:
                        product = pair(r, normals[4])
                        assert product.l > 0
                        result = dict(status='OBTUSE_TO_OLD_WALL', wall=4, product=product.enc())
                    else:
                        product = pair(r, normals[50])
                        if index == 0:
                            assert product.h < 0 and product.l > -1
                            bracket = [2, 3]
                        else:
                            assert index == 1 and product.h < -1
                            assert (product + U.theta(4)).l > 0
                            bracket = [3, 4]
                        result = dict(status='ALL_LABEL_COXETER_GAP', wall=50,
                                      product=product.enc(), consecutive_labels=bracket)
                    decisions.append(dict(profile=index, branch=sign, delta=delta.enc(), **result))
        item['status'] = 'ALL_REAL_SINGLE_WALL_TRUNCATION_EXCLUDED'
        item['decisions'] = decisions
        records.append(item)
    assert Counter(r['status'] for r in records) == dict(
        ALL_REAL_SINGLE_WALL_TRUNCATION_EXCLUDED=3, DOMAIN_ONLY_NOT_EXCLUDED_OR_ACCEPTED=2)
    counts = Counter(d['status'] for r in records for d in r.get('decisions', []))
    assert counts == dict(OBTUSE_TO_OLD_WALL=2, ALL_LABEL_COXETER_GAP=2,
                          UNIQUE_NORMAL_HAS_WRONG_NORM=168, TAIL_NORM_ALWAYS_ABOVE2=8)
    paths = [f'scripts/{TAG}.py', f'scripts/{E.TAG}.py', f'results/{E.TAG}.json',
             f'scripts/{P.TAG}.py', f'results/{P.TAG}.json', P.U.CAP, P.C53.CL,
             P.CAP52, P.CAPS55, P.CONE55,
             'results/fresh_nonarithmetic_nonright_s6_20261016_cone_review.json']
    return dict(status='PASS_THREE_ALL_LABEL_SINGLE_WALL_PATCH_OBSTRUCTIONS',
                sources={p: sha(p) for p in paths}, retained_walls=55, face_support=P.S,
                hyperideal_vertices=7, patches=5, boundary_open_ends=17,
                label_code6='The full tail c=2cos(pi/m),integer m>=6; proof relaxes to sqrt3<=c<2.',
                patch_records=records, decision_counts=dict(counts),
                conclusion='No single actual S6-facet Coxeter wall can terminate all three selected open edge ends of patch[7,30,49], or all four selected ends of either patch[8,30,44,52] or[30,34,47,52]. All fields and all finite labels are covered.',
                unresolved='The other two patches have only a complete necessary-profile list here. Coordinated multiple-wall compactification of any patch is not excluded. No new cap or finite-volume claim.',
                proof='A compact completion must terminate every listed old geodesic edge at a finite vertex. Each selected endpoint has a positive-definite rank9 ambient Gram, so the frozen complete endpoint profiles apply. Their intersections give all listed profiles and tails. Nine prescribed contacts leave an affine line whose norm equation has exactly the two recorded roots; ten independent contacts leave a unique normal or the recorded one-parameter affine family. The exact outward sign bounds exhaust these branches without assuming a finite angle alphabet on unspecified contacts.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
        print('PASS saved all-label single-wall patch replay')
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['decision_counts'])
