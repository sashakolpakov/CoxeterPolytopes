#!/usr/bin/env python3
"""Remaining eight symbolic high-S families, with independent high parameters."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_oblique_high_s_complete_low_r_interval_20261018 as C

P, Q, ROOT = C.P, C.Q, C.ROOT
TAG = 'fresh_oblique_nonright_high_s_interval_20261018'
V, add, sub, mul, div, neg, sqrt, dot = C.V, C.add, C.sub, C.mul, C.div, C.neg, C.sqrt, C.dot
FAMILIES = [[0, [k, 2, 6], [2]] for k in range(2, 6)] + [[0, [6, 2, k], [0]] for k in range(3, 6)] + [[0, [6, 2, 6], [0, 2]]]


def rectangle(basis, rcoeffs, family):
    profile, labels, high = family
    assert profile == 0 and labels[1] == 2 and labels[2] >= 3
    theta = (P.theta(6)[0], F(2)) if 0 in high else P.theta(labels[0])
    tau = (P.theta(6)[0], F(2)) if 2 in high else P.theta(labels[2])
    aa, ac, cc, a8, c8, sigma = (basis[k] for k in ('AA', 'AC', 'CC', 'A8', 'C8', 'residual8_square'))
    rho = sub(V(2), add(add(mul(aa, mul(theta, theta)),
        mul(V(2), mul(ac, mul(theta, tau)))), mul(cc, mul(tau, tau))))
    assert rho[0] > 0
    radical = sqrt(rho)
    positive8 = add(add(mul(theta, a8), mul(tau, c8)), mul(radical, sqrt(sigma)))
    assert positive8[0] > 0
    joins = []
    for r in rcoeffs:
        value = sub(add(mul(theta, r['AR']), mul(tau, r['CR'])), mul(radical, r['BR']))
        assert value[1] < -2
        joins.append(dict(R_index=r['index'], full_SR_range=value))
    return dict(family=family, theta_domain=theta, tau_domain=tau, norm_remainder_range=rho,
        positive_S_sheet_old8=positive8, negative_S_R_joins=joins)


def run():
    old = json.loads((ROOT/f'results/{P.H.TAG}.json').read_text())
    prior = json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    cube = json.loads((ROOT/f'results/{C.TAG}.json').read_text())
    for doc in (old, prior, cube):
        for path, digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    assert cube['summary']['exact_R_sheets'] == 1120 and cube['summary']['positive_old8_R_sheets'] == 560
    assert len(old['records']) == len(prior['parents']) == 2
    assert len(cube['records']) == len(prior['records']) == 20
    profiles = [[V(q) for q in row] for row in prior['profiles']]
    g0 = {int(a): {int(b): V(q) for b, q in row.items()}
          for a, row in prior['primary_projection']['projected_Gram'].items()}
    wc0 = {int(a): [V(q) for q in row]
           for a, row in prior['primary_projection']['old_white_covectors'].items()}
    rows = []
    for parent, source in zip(prior['parents'], old['records']):
        assert parent['index'] == source['index'] and parent['instruction_key'] == source['instruction_key']
        g, wc = {a: row.copy() for a, row in g0.items()}, wc0.copy()
        for step in parent['steps']:
            Q.restore(g, wc, profiles, step)
        assert len(source['caps']) == 10
        for cap in source['caps']:
            key = source['index'], cap['index']
            pm = [r for r in prior['records'] if (r['parent'], r['N_index']) == key]
            cm = [r for r in cube['records'] if (r['parent'], r['N_index']) == key]
            assert len(pm) == len(cm) == 1
            assert [f for f in cap['S_bank']['symbolic_high_families'] if f[1] != [6, 2, 2]] == FAMILIES
            gn, wn = {a: row.copy() for a, row in g.items()}, wc.copy()
            Q.restore(gn, wn, profiles, pm[0]['N_extension'])
            frame = [52, 68, 71]
            gram = [[gn[a][b] for b in frame] for a in frame]
            inv, det = P.inverse(gram)
            a, c, cv = [V(-1), V(0), V(0)], [V(0), V(0), V(-1)], [gn[i][8] for i in frame]
            sigma = sub(gn[8][8], dot(cv, inv, cv))
            assert det[1] < 0 and gram[0][0][0] > 0 and sigma[0] > 0
            basis = dict(frame=frame, Gram=gram, inverse=inv, determinant=det, old8_covector=cv,
                old8_square=gn[8][8], residual8_square=sigma,
                AA=dot(a, inv, a), AC=dot(a, inv, c), CC=dot(c, inv, c),
                A8=dot(a, inv, cv), C8=dot(c, inv, cv))
            coeffs = []
            assert len(cm[0]['recipes']) == 28
            for r in cm[0]['recipes']:
                contacts = [V(r['negative_R_sheet_products'][str(i)]) for i in frame]
                ar, cr = dot(a, inv, contacts), dot(c, inv, contacts)
                br = div(sub(V(r['negative_R_sheet_products']['8']), dot(cv, inv, contacts)), sqrt(sigma))
                assert P.I.enc(ar) == r['continuous']['A_dot_R'] and P.I.enc(br) == r['continuous']['B_dot_R']
                assert V(r['positive_R_sheet_old8'])[0] > 0
                coeffs.append(dict(index=r['index'], profile=r['profile'], RM_label=r['RM_label'], AR=ar, CR=cr, BR=br))
            assert [(r['profile'], r['RM_label']) for r in coeffs] == list(product(range(7), range(2, 6)))
            families = [rectangle(basis, coeffs, family) for family in FAMILIES]
            rows.append(dict(index=len(rows), parent=key[0], N_index=key[1], N_labels=cap['labels'],
                S_basis=basis, R_coefficients=coeffs, families=families))
            print('NONRIGHT_HIGH_S', key, '8families224rectangles_PASS', flush=True)
    assert len(rows) == len({(r['parent'], r['N_index']) for r in rows}) == 20
    families = [f for r in rows for f in r['families']]
    joins = [j for f in families for j in f['negative_S_R_joins']]
    assert len(families) == 160 and len(joins) == 4480
    debit = old['exact_deletion_debit']
    assert debit['lost_finite_count'] == 19 and debit['compact_faces_by_dimension'] == {'2': 12, '3': 1}
    paths = {Path(__file__), Path(C.__file__), Path(Q.__file__), Path(P.__file__), Path(P.I.__file__)}
    paths.update(ROOT/path for doc in (old, prior, cube) for path in doc['sources'])
    paths.update(ROOT/f'{folder}/{stem}.{ext}' for stem in (P.H.TAG, P.TAG, C.TAG)
                 for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return P.I.enc(dict(status='PASS_EIGHT_NONRIGHT_HIGH_S_FAMILIES_PHYSICAL_INTERVAL_EXCLUSION',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        interval_scale=P.I.PREC, symbolic_families=FAMILIES, records=rows,
        exact_deletion_debit=debit, deleted_walls=old['deleted_walls'], topology=old['topology'],
        high_RM_nonright_SM_obstruction=dict(RM_square_at_least=3, SM_square_at_least=1,
            maximum_three_root_determinant=0, strict_positive_definiteness_impossible=True),
        summary=dict(parent_classes=2, parent_origins=4, N_geometries=20, symbolic_family_types=8,
            continuous_S_families=160, exact_R_recipes=560, inherited_R_sheets=1120,
            independent_two_high_families=20, positive_S_old8_sheets=160,
            continuously_dotted_rectangles=4480, combined_R_S_sheet_pairs=17920,
            minimum_S_norm_remainder=min(f['norm_remainder_range'][0] for f in families),
            minimum_positive_S_old8=min(f['positive_S_sheet_old8'][0] for f in families),
            maximum_product_upper=max(j['full_SR_range'][1] for j in joins),
            unresolved=0, accepted=0, new_cones=0),
        scope='Physical interval exclusion for eight remaining symbolic high-S families, all7x4 low-R recipes and both S/R sheets on20 fixed N geometries. Both high parameters remain independent. Uses frozen complete physical R cube, not a finite-corner filter or numerical bank rejection. Other N/parent geometries and incidences remain open; no volume claim.'))


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
