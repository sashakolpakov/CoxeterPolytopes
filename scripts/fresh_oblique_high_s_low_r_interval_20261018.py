#!/usr/bin/env python3
"""Continuous high-S exclusion for the 172 frozen, numerically selected R recipes.

Physical interval reconstruction specifies each recipe exactly. This does not
recertify completeness of the earlier numerical native-wall bank.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import fresh_oblique_high_s_polar_interval_20261018 as P

ROOT, TAG = P.ROOT, 'fresh_oblique_high_s_low_r_interval_20261018'
V, add, sub, mul, div, neg, sqrt = P.V, P.add, P.sub, P.mul, P.div, P.neg, P.sqrt
dot = P.I.dot


def restore(g, wc, profiles, row):
    """Restore a hash-bound, previously checked physical interval extension."""
    wall = row['wall']
    assert wall not in g
    products = {int(i): V(v) for i, v in row['new_products'].items()}
    assert set(products) == set(g)
    for i, value in products.items():
        g[i][wall] = value
    g[wall] = {**products, wall: V(row['projected_square'])}
    wc[wall] = profiles[row['profile']]


def continuous(g, inherited):
    frame = [52, 68, 71]
    inv, det = P.inverse([[g[a][b] for b in frame] for a in frame])
    rhs, cv = [V(-1), V(0), V(0)], [g[i][8] for i in frame]
    alpha = neg(dot(rhs, inv, rhs))
    q8 = sub(g[8][8], dot(cv, inv, cv))
    assert det[1] < 0 and alpha[0] > 0 and q8[0] > 0
    assert P.I.enc(alpha) == inherited['alpha']
    assert P.I.enc(q8) == inherited['residual8_square']
    contacts = [g[i][70] for i in frame]
    a = dot(rhs, inv, contacts)
    b = div(sub(g[8][70], dot(cv, inv, contacts)), sqrt(q8))
    z = (div(P.theta(6), sqrt(add(V(2), mul(V(3), alpha))))[0],
         div(V(2), sqrt(add(V(2), mul(V(4), alpha))))[1])
    assert P.I.enc(z) == inherited['transformed_z_domain']
    # S has white profile zero, so this projected product is the full product.
    derivative = sub(a, mul(mul(alpha, z), b))
    assert derivative[1] < 0
    def evaluate(theta):
        return sub(mul(theta, a), mul(sqrt(add(V(2), mul(alpha, mul(theta, theta)))), b))
    at6, limit = evaluate(P.theta(6)), evaluate(V(2))
    full_range = limit[0], at6[1]
    assert full_range[1] < -2
    return dict(S_frame=frame, S_frame_inverse=inv, S_frame_determinant=det,
        alpha=alpha, residual8_square=q8, A_dot_R=a, B_dot_R=b,
        theta_domain=[P.theta(6)[0], F(2)], transformed_z_domain=z,
        derivative_range=derivative, at_label6=at6, at_limit=limit,
        full_product_range=full_range, monotonicity='strictly_decreasing')


def run():
    old = json.loads((ROOT/f'results/{P.H.TAG}.json').read_text())
    prior = json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    for document in (old, prior):
        for path, digest in document['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    assert prior['status'] == 'PASS_FIXED_HIGH_N_HIGH_S_HIGH_RM_PHYSICAL_INTERVAL_OBSTRUCTION'
    assert len(old['records']) == len(prior['parents']) == 2
    assert len(prior['records']) == 20
    profiles = [[V(q) for q in row] for row in prior['profiles']]
    wi = [[V(q) for q in row] for row in prior['primary_projection']['white_inverse']]
    g0 = {int(a): {int(b): V(q) for b, q in row.items()}
          for a, row in prior['primary_projection']['projected_Gram'].items()}
    wc0 = {int(a): [V(q) for q in row]
           for a, row in prior['primary_projection']['old_white_covectors'].items()}
    numerical_old = {int(k): np.array(v) for k, v in old['numerical_snapshot']['old_normals'].items()}
    records, domain = [], []
    for source, parent in zip(old['records'], prior['parents']):
        assert source['index'] == parent['index']
        assert source['instruction_key'] == parent['instruction_key']
        g, wc = {a: row.copy() for a, row in g0.items()}, wc0.copy()
        numerical = {**numerical_old, **{int(i): np.array(v) for i, v in source['prefix_normals'].items()}}
        for step in parent['steps']:
            restore(g, wc, profiles, step)
        assert len(source['caps']) == 10
        for cap in source['caps']:
            matches = [row for row in prior['records']
                       if (row['parent'], row['N_index']) == (source['index'], cap['index'])]
            assert len(matches) == 1
            previous = matches[0]
            assert V(previous['continuous']['positive_S_sheet_old8'])[0] > 0
            assert [r for r in cap['S_bank']['symbolic_high_families'] if r[1][2] == 2] == [[0, [6, 2, 2], [0]]]
            gn, wn = {a: row.copy() for a, row in g.items()}, wc.copy()
            restore(gn, wn, profiles, previous['N_extension'])
            # The ancestral R constructed M; it is not the R being tested.
            assert 70 in gn
            del gn[70], wn[70]
            for row in gn.values():
                del row[70]
            nums = {**numerical, 73: np.array(cap['normal'])}
            assert [r['index'] for r in cap['R_bank']['caps']] == list(range(len(cap['R_bank']['caps'])))
            for r in cap['R_bank']['caps']:
                assert r['labels'][0] == r['labels'][2] == 2 and 2 <= r['labels'][1] <= 5
                assert 0 <= r['profile'] < 7 and r['branch'] in (-1, 1)
                key = [source['index'], cap['index'], r['index']]
                domain.append(key)
                full, white = {a: row.copy() for a, row in gn.items()}, wn.copy()
                proof = P.extension(full, white, wi, profiles,
                    {i: nums[i] for i in gn}, 70, [69, 71, 73], r['profile'], r['labels'], np.array(r['normal']))
                proof['saved_numerical_branch'] = r['branch']
                record = dict(index=len(records), key=key, N_labels=cap['labels'],
                    R_profile=r['profile'], R_labels=r['labels'], R_branch=r['branch'],
                    R_extension=proof, continuous=continuous(full, previous['continuous']))
                records.append(record)
            print('CONTINUOUS_LOW_R', source['index'], cap['labels'], len(cap['R_bank']['caps']), flush=True)
    assert len(domain) == len({tuple(k) for k in domain}) == len(records) == 172
    assert len({tuple(k[:2]) for k in domain}) == 20
    debit = old['exact_deletion_debit']
    assert debit['lost_finite_count'] == 19 and debit['compact_faces_by_dimension'] == {'2': 12, '3': 1}
    paths = {Path(__file__), Path(P.__file__), Path(P.I.__file__), Path(P.H.__file__)}
    paths.update(ROOT/path for document in (old, prior) for path in document['sources'])
    paths.update(ROOT/f'{folder}/{stem}.{ext}' for stem in (P.TAG, P.H.TAG)
                 for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return P.I.enc(dict(status='PASS_SELECTED_172_LOW_R_CONTINUOUS_HIGH_S_PHYSICAL_INTERVAL_EXCLUSION',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        interval_scale=P.I.PREC, selected_domain=domain, records=records,
        inherited_physical_projection=prior['primary_projection'], inherited_profiles=prior['profiles'],
        exact_deletion_debit=debit, deleted_walls=old['deleted_walls'], topology=old['topology'],
        necessary_conditions=dict(S_profile=0, S52_label_at_least=6, SL=2, SM=2, RX=2, RN=2,
            RM_labels=[2, 3, 4, 5], R_profile_domain=list(range(7))),
        summary=dict(parent_classes=2, parent_origins=4, N_geometries=20, selected_R_recipes=172,
            physical_R_extensions=172, excluded_positive_S_sheets=20, strictly_decreasing_products=172,
            continuously_dotted_joins=172,
            maximum_product_upper=max(r['continuous']['full_product_range'][1] for r in records),
            maximum_derivative_upper=max(r['continuous']['derivative_range'][1] for r in records),
            unresolved=0, accepted=0, new_cones=0),
        scope='Physical rational interval exclusion of every high-S parameter for exactly172 frozen numerically selected R recipes on20 fixed N geometries. The prior numerical native-bank rejection domain is not recertified; this is not an all-low-R, all-N, all-incidence or chamber-volume theorem.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'
    assert args.verify or not path.exists()
    result = run()
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['summary'], flush=True)
