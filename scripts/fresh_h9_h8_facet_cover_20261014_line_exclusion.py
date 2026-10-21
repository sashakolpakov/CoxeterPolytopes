#!/usr/bin/env python3
"""Exclude the fixed H8 attachment on the entire all-dotted H--E line."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_facet_cover_20261014_line_exclusion'
SPLIT = F('.9967505')

def run():
    paths = ['results/fresh_h9_h8_facet_cover_20261014_parameter_gate.json',
             'results/fresh_h9_h8_facet_cover_20261014_gap_gate.json',
             'results/fresh_h9_mixed_d4_coupling_20261014.json',
             'results/fresh_h9_mixed_d4_global_caps_20261014.json']
    docs = [json.loads((ROOT/p).read_text()) for p in paths]
    pending, seen = list(zip(paths, docs)), set()
    while pending:
        path, doc = pending.pop()
        if path in seen: continue
        seen.add(path)
        for p, h in doc.get('sources', {}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
            if p.endswith('.json'): pending.append((p, json.loads((ROOT/p).read_text())))
    param, gap, base, cap = docs
    assert param['status'] == 'PASS_EXACT_MONOTONE_PARAMETER_DIAMETER_GATE'
    assert gap['status'] == 'PASS_FIXED_LITERAL_CAP_H8_SOURCE_ATTACHMENT_EXCLUDED'
    small = A.E(A.Q(*gap['exceptional_minimum_norm']))
    bound = A.E(A.Q(*gap['ordinary_minimum_norm']))
    assert bound == 4*small
    assert sorted(r['source_vertex'] for r in gap['exceptional_rays']) == [0, 32, 33]
    assert all(any(A.Q(*x) == 1 for x in r['contacts']) for r in gap['exceptional_rays'])
    lower, upper = map(F, param['necessary_open_t_interval'])
    assert lower < SPLIT < upper
    Q, N, H = A.dm(base['form']), A.dm(base['final_normals']), A.dv(base['center'])
    E = [A.E(F(x)) for x in cap['literal_escape']]
    z = [h*(1-SPLIT)+e*SPLIT for h, e in zip(H, E)]
    q = A.pair(z, z, Q); assert q.sign() > 0
    g = [A.pair(z, n, Q) for n in N]
    norms = [A.pair(n, n, Q) for n in N]
    kappa = [g[i]*g[i]/q-norms[i] for i in range(12)]
    assert all(v.sign() > 0 for v in kappa)
    assert all(norms[i] == 2 for i in (0, 1, 2, 3, 7))
    low7 = A.dec(param['lower_projected_norm_squares'][7])
    assert (2*small-low7).sign() > 0 and (kappa[7]-small).sign() > 0
    assert all((bound-kappa[i]).sign() > 0 for i in range(4))
    # Strictly negative off-diagonal products suffice for distinct projected
    # center directions throughout every all-dotted part of the line.
    # Even original right pairs give cosh distance>1, since gamma_i>1.
    pairs = []
    for i in range(4):
        for j in range(i+1, 4):
            old_pair = A.pair(N[i], N[j], Q); assert old_pair.sign() <= 0
            pairs.append(dict(pair=[i, j], old_product=A.enc(old_pair)))
    paths += [f'scripts/{TAG}.py', 'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return dict(status='PASS_ENTIRE_DOTTED_CAP_LINE_FIXED_H8_ATTACHMENT_EXCLUDED',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                necessary_t_interval=param['necessary_open_t_interval'], split_parameter=str(SPLIT),
                small_norm=A.enc(small), ordinary_bound=A.enc(bound),
                lower_wall7_margin=A.enc(2*small-low7),
                split_wall7_margin=A.enc(kappa[7]-small),
                split_four_center_margins=A.enc([bound-kappa[i] for i in range(4)]),
                split_projected_norm_squares=A.enc(kappa), distinct_center_pairs=pairs,
                early_interval='lower<t<=split: a<kappa7^2<2a forces an exceptional source-vertex ray with a nonzero contact strictly between1 andsqrt2, impossible for Coxeter angles.',
                late_interval='split<=t<upper: four distinct centers0,1,2,3 all have kappa_i^2<4a and must occupy three exceptional source-vertex rays, impossible.',
                scope='Entire all-dotted part of the specified H--E segment, all isometric placements of the fixed11-facet compact BugaenkoH8 source perpendicular to the new cap, retaining all12 old walls and imposing acute Coxeter pairs and finite volume. Different cap directions, source polytopes, nonperpendicular sidewalls, extra source facets and changed retained walls remain outside the theorem.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = run(); out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['split_parameter'])
