#!/usr/bin/env python3
"""Both exact R sheets for all7 profiles x4 low RM labels on20 fixed N roots.

No numerical R-bank rejection or sampled high-S label is used in this proof.
The parent and N recipes retain their frozen residual-witness definitions.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_oblique_high_s_low_r_interval_20261018 as Q

P, ROOT = Q.P, Q.ROOT
TAG = 'fresh_oblique_high_s_complete_low_r_interval_20261018'
V, add, sub, mul, div, neg, sqrt, dot = Q.V, Q.add, Q.sub, Q.mul, Q.div, Q.neg, Q.sqrt, Q.dot


def complete_cube(g, wc, wi, profiles, inherited):
    frame, witness = [69, 71, 73], 10
    gram = [[g[a][b] for b in frame] for a in frame]
    inv, det = P.inverse(gram)
    assert det[1] < 0 and gram[0][0][0] > 0
    cv = [g[i][witness] for i in frame]
    sigma = sub(g[witness][witness], dot(cv, inv, cv))
    assert sigma[0] > 0
    rows = []
    for profile, label in product(range(7), range(2, 6)):
        labels = [2, label, 2]
        rhs = [sub(neg(P.theta(m)), dot(profiles[profile], wi, wc[i]))
               for i, m in zip(frame, labels)]
        square = sub(V(2), dot(profiles[profile], wi, profiles[profile]))
        rho = sub(square, dot(rhs, inv, rhs))
        assert rho[0] > 0
        factor = sqrt(div(rho, sigma))
        sheets = {}
        for sign in (-1, 1):
            contacts = {}
            for j in g:
                cj = [g[i][j] for i in frame]
                contacts[j] = rhs[frame.index(j)] if j in frame else add(dot(rhs, inv, cj),
                    mul(mul(V(sign), factor), sub(g[witness][j], dot(cv, inv, cj))))
            sheets[sign] = contacts
        positive8 = add(sheets[1][8], dot(profiles[profile], wi, wc[8]))
        assert positive8[0] > 0
        full = {a: row.copy() for a, row in g.items()}
        for i, value in sheets[-1].items():
            full[i][70] = value
        full[70] = {**sheets[-1], 70: square}
        continuous = Q.continuous(full, inherited)
        rows.append(dict(index=len(rows), profile=profile, RM_label=label, R_labels=labels,
            rhs=rhs, projected_square=square, norm_remainder=rho,
            positive_R_sheet_old8=positive8,
            negative_R_sheet_products={i: sheets[-1][i] for i in (8, 52, 68, 71)},
            continuous=continuous))
    assert [(r['profile'], r['RM_label']) for r in rows] == list(product(range(7), range(2, 6)))
    return dict(R_frame=frame, frame_Gram=gram, inverse=inv, determinant=det,
        witness=witness, witness_covector=cv, witness_residual=sigma, recipes=rows)


def run():
    old = json.loads((ROOT/f'results/{P.H.TAG}.json').read_text())
    prior = json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    selected = json.loads((ROOT/f'results/{Q.TAG}.json').read_text())
    for doc in (old, prior, selected):
        for path, digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    assert selected['summary']['selected_R_recipes'] == 172
    assert len(old['records']) == len(prior['parents']) == 2 and len(prior['records']) == 20
    profiles = [[V(q) for q in row] for row in prior['profiles']]
    wi = [[V(q) for q in row] for row in prior['primary_projection']['white_inverse']]
    g0 = {int(a): {int(b): V(q) for b, q in row.items()}
          for a, row in prior['primary_projection']['projected_Gram'].items()}
    wc0 = {int(a): [V(q) for q in row]
           for a, row in prior['primary_projection']['old_white_covectors'].items()}
    records = []
    for source, parent in zip(old['records'], prior['parents']):
        assert source['index'] == parent['index'] and source['instruction_key'] == parent['instruction_key']
        g, wc = {a: row.copy() for a, row in g0.items()}, wc0.copy()
        for step in parent['steps']:
            Q.restore(g, wc, profiles, step)
        assert len(source['caps']) == 10
        for cap in source['caps']:
            matches = [r for r in prior['records']
                       if (r['parent'], r['N_index']) == (source['index'], cap['index'])]
            assert len(matches) == 1
            pr = matches[0]
            assert [r for r in cap['S_bank']['symbolic_high_families'] if r[1][2] == 2] == [[0, [6, 2, 2], [0]]]
            assert V(pr['continuous']['positive_S_sheet_old8'])[0] > 0
            gn, wn = {a: row.copy() for a, row in g.items()}, wc.copy()
            Q.restore(gn, wn, profiles, pr['N_extension'])
            del gn[70], wn[70]
            for row in gn.values():
                del row[70]
            cube = complete_cube(gn, wn, wi, profiles, pr['continuous'])
            records.append(dict(index=len(records), parent=source['index'], N_index=cap['index'],
                N_labels=cap['labels'], positive_S_sheet_old8=pr['continuous']['positive_S_sheet_old8'], **cube))
            print('COMPLETE_CUBE', source['index'], cap['labels'], '28recipes56sheets_PASS', flush=True)
    assert len(records) == len({(r['parent'], r['N_index']) for r in records}) == 20
    rows = [r for parent in records for r in parent['recipes']]
    assert len(rows) == 560
    debit = old['exact_deletion_debit']
    assert debit['lost_finite_count'] == 19 and debit['compact_faces_by_dimension'] == {'2': 12, '3': 1}
    paths = {Path(__file__), Path(Q.__file__), Path(P.__file__), Path(P.I.__file__)}
    paths.update(ROOT/path for doc in (old, prior, selected) for path in doc['sources'])
    paths.update(ROOT/f'{folder}/{stem}.{ext}' for stem in (P.H.TAG, P.TAG, Q.TAG)
                 for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return P.I.enc(dict(status='PASS_COMPLETE_FIXED_LOW_R_CUBE_CONTINUOUS_HIGH_S_INTERVAL_EXCLUSION',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        interval_scale=P.I.PREC, records=records,
        inherited_physical_projection=prior['primary_projection'], inherited_profiles=prior['profiles'],
        exact_deletion_debit=debit, deleted_walls=old['deleted_walls'], topology=old['topology'],
        necessary_conditions=dict(N_profile=0, NX_label_at_least=6, S_profile=0,
            S52_label_at_least=6, SL=2, SM=2, RX=2, RN=2,
            R_profile_domain=list(range(7)), RM_labels=[2, 3, 4, 5]),
        summary=dict(parent_classes=2, parent_origins=4, N_geometries=20, complete_R_recipes=560,
            exact_R_sheets=1120, positive_norm_remainders=560, excluded_positive_S_sheets=20,
            positive_old8_R_sheets=560, continuously_dotted_negative_R_sheets=560,
            minimum_norm_remainder=min(r['norm_remainder'][0] for r in rows),
            minimum_positive_R_old8=min(r['positive_R_sheet_old8'][0] for r in rows),
            maximum_product_upper=max(r['continuous']['full_product_range'][1] for r in rows),
            maximum_derivative_upper=max(r['continuous']['derivative_range'][1] for r in rows),
            unresolved=0, accepted=0, new_cones=0),
        scope='Complete physical interval exclusion for all7 white profiles x4 low RM labels xboth norm sheets on20 fixed N recipes, with S profile0 and S52>=6, SL=SM=right. No numerical R-bank rejection is used. Other N parameters, high-S families with SM nonright and other incidences remain outside the theorem; no chamber-volume conclusion.'))


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
