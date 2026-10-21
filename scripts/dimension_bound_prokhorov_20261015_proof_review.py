#!/usr/bin/env python3
"""Small exact checks and bindings for the separate geometric/source audit.

The geometric conclusions are proved in the matching report. This is not a
second verification of Prokhorov's published topological lemmas.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'dimension_bound_prokhorov_20261015_proof_review'

def run():
    inputs = [
        'scripts/dimension_bound_prokhorov_20261015.py',
        'results/dimension_bound_prokhorov_20261015.json',
        'scripts/dimension_bound_prokhorov_20261015_review.py',
        'results/dimension_bound_prokhorov_20261015_review.json',
        'scripts/dimension_bound_vertex_components_20261015.py',
        'results/dimension_bound_vertex_components_20261015.json',
    ]
    primary = json.loads((ROOT / inputs[1]).read_text())
    arithmetic = json.loads((ROOT / inputs[3]).read_text())
    assert primary['derived_upper_dimension'] == 990
    assert arithmetic['connected_graph_checks'] == 160
    for data in (primary, arithmetic):
        assert data['status'].startswith('PASS_')
        for path, digest in data['sources'].items():
            assert sha256((ROOT / path).read_bytes()).hexdigest() == digest

    # The two closed graph-distance formulas supply the only large-rank cases.
    ksum = 2 * sum(range(1, 8)) + sum(range(1, 16))
    assert ksum == 176
    assert ksum - (2 + 1) == 173
    small = []
    for r in range(1, 17):
        deficit = 29 * r - F(3 * r * (r - 1), 2)
        assert deficit >= F(13 * r, 2) > 0
        small.append([r, str(deficit)])
    assert F(13 * 27, 2) > 173

    # Polynomial identity, not a finite-dimensional test as a tail argument.
    assert [991**2 - 997*991 + 6532, 2*991 - 997, 1] == [586, 985, 1]
    rows = []
    for n in (990, 991, 992):
        c = F(29, 3) - F(173, 3*(n-1))
        p = n*n - 997*n + 6532
        assert F(n) - (96*c + 68) == F(p, n-1)
        rows.append(dict(n=n, C=str(c), necessary_strict_upper=str(96*c+68),
                         polynomial=p))
    paths = [f'scripts/{TAG}.py'] + inputs
    return dict(
        status='PASS_PRIMARY_GEOMETRY_AND_SCOPE_AUDIT',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        exact_checks=dict(path_deficit=176, D_deficit=173,
                          small_component_crude_bounds=small,
                          excluded_tail_shift_coefficients=[586,985,1],
                          boundary_values=rows),
        published_inputs=[
            'Prokhorov 1987 p408: every edge-angle scheme is elliptic of rank n-1.',
            'Prokhorov 1987 p410 Proposition5: unordered marked-node pairs bound edge incidences.',
            'Prokhorov 1987 p406 Lemma2: positive C with unchanged good-face weights gives n<96C+68.',
            'Finite irreducible real Coxeter classification: only paths and D have rank>=17.',
        ],
        conclusion='The bound n<=990 follows for finite-sided finite-volume real hyperbolic Coxeter polytopes from the explicitly cited published inputs.',
        not_claimed='No independent reproof of all published cusp/topology cases; no publication priority, accepted-best-bound, or verified Nikulin300 proof claim.',
        literature_status=dict(
            nikulin_primary='ICM1986 vol1 printed654,657-659; CHRG means all finite-covolume reflection groups, not arithmetic only.',
            nikulin_announcement='Printed659 asserts dimension<300 with an outline; the extracted report gives an affine-A cusp example, not a complete numerical treatment of all cusp types.',
            source_url='https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf',
            source_pdf_sha256='0375d371245a73ed40ebd95ffc66846b2990d54fd25a447f75ecd358f4d1fcff',
            source_pdf_pages=974,
            nikulin_pdf_pages_one_based=[756,773],
            dolgachev_url='https://sites.lsa.umich.edu/idolga/wp-content/uploads/sites/1334/2024/08/reflections08.pdf',
            resolution='990 is an audited explicit refinement of Prokhorov; the older stronger announcement prevents a new-current-best claim.'),
        conditional_compact_audit='PASS source/proof audit: rank-r component bound d<=4r; quadrilateral charge additionally requires full flag/3-free missing-face condition, not only absence of triangular two-faces.',
    )

if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true')
    args=parser.parse_args(); data=run(); out=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'])
