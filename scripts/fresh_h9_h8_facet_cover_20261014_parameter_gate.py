#!/usr/bin/env python3
"""Exact necessary diameter bounds for one H8-source attachment line."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate import I, ei, SCALE, BITS
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_facet_cover_20261014_parameter_gate'
LEFT = F('.99673989')
NULL_UPPER = F('.996739892')
LOWER = F('.99674949')
UPPER = F('.996752404')

def run():
    source = 'data/known_examples/bugaenko_compact_h8.json'
    base = 'results/fresh_h9_mixed_d4_coupling_20261014.json'
    cap = 'results/fresh_h9_mixed_d4_global_caps_20261014.json'
    src, d, cc = [json.loads((ROOT/p).read_text()) for p in (source, base, cap)]
    def sf(x): return A.E(A.Q(F(x['a'])+F(x['b'])/2, F(x['b'])/2))
    SQ = [[sf(x) for x in r] for r in src['lattice_gram']]
    SN = [[sf(x) for x in r] for r in src['simple_roots']]
    SH = [sf(x) for x in src['geometry']['interior_point']]
    assert A.pair(SH, SH, SQ).sign() < 0
    V = []
    for ids in src['geometry']['vertices']:
        v = A.kernel_one([A.mv(SQ, SN[i]) for i in ids]); assert v is not None
        if A.pair(v, SH, SQ).sign() > 0: v = [-x for x in v]
        norm = A.pair(v, v, SQ); assert norm.sign() < 0
        products = [A.pair(v, n, SQ) for n in SN]
        assert all(p.sign() <= 0 for p in products)
        assert [i for i, p in enumerate(products) if not p] == list(ids)
        V.append((v, norm))
    assert len(V) == 41
    distances = []
    for i in range(41):
        for j in range(i+1, 41):
            v, qv = V[i]; w, qw = V[j]
            cosh = -ei(A.pair(v, w, SQ))/(ei(qv*qw).sqrt())
            distances.append((i, j, cosh))
    source_lower = next(x for i, j, x in distances if (i, j) == (32, 33))
    source_upper = I.raw(max(x.lo for i, j, x in distances), max(x.hi for i, j, x in distances))
    Q, N, H = A.dm(d['form']), A.dm(d['final_normals']), A.dv(d['center'])
    E = [A.E(F(x)) for x in cc['literal_escape']]
    hh, he, ee = [A.pair(x, y, Q) for x, y in ((H, H), (H, E), (E, E))]
    a, b, c = hh+ee-2*he, 2*(he-hh), hh
    q = lambda t: a*t*t+b*t+c
    assert a.sign() < 0 and (2*a+b).sign() > 0
    assert q(LEFT).sign() < 0 and q(NULL_UPPER).sign() > 0
    G = [[A.pair(n, m, Q) for m in N] for n in N]
    gh = [A.pair(H, n, Q) for n in N]; ge = [A.pair(E, n, Q) for n in N]
    dg = [e-h for h, e in zip(gh, ge)]
    monotone = []
    for i in range(12):
        assert gh[i].sign() < 0 and ge[i].sign() < 0 and G[i][i].sign() > 0
        assert all(G[i][j].sign() <= 0 for j in range(12) if j != i)
        # gamma' = -(2g'q-gq')/(2 sqrt(nn) q^(3/2)).
        lin, const = dg[i]*b-2*a*gh[i], 2*dg[i]*c-gh[i]*b
        values = [lin*LEFT+const, lin+const]
        assert all(x.sign() > 0 for x in values)
        monotone.append(dict(wall=i, linear=A.enc(lin), constant=A.enc(const), endpoint_values=A.enc(values)))
    def geometry(t):
        qq = q(t); assert qq.sign() > 0
        g = [h+v*t for h, v in zip(gh, dg)]
        assert all((g[i]*g[i]-qq*G[i][i]).sign() > 0 for i in range(12))
        PP = [[G[i][j]-g[i]*g[j]/qq for j in range(12)] for i in range(12)]
        k = [(-ei(PP[i][i])).sqrt() for i in range(12)]
        sr = [ei(G[i][i]).sqrt()/k[i] for i in range(12)]
        cr = [-ei(g[i])/(ei(qq).sqrt()*k[i]) for i in range(12)]
        rows = []
        for i in range(12):
            for j in range(i, 12):
                dc = I(1) if i == j else -ei(PP[i][j])/(k[i]*k[j])
                rad = dc.square()-1
                ds = I.raw(max(0, rad.lo), rad.hi).sqrt()
                bound = dc*(cr[i]*cr[j]+sr[i]*sr[j])+ds*(sr[i]*cr[j]+cr[i]*sr[j])
                rows.append(dict(pair=[i, j], center_cosh=dc.enc(), ball_union_cosh=bound.enc()))
        return rows, [-PP[i][i] for i in range(12)]
    lower_rows, lower_kappa = geometry(LOWER)
    upper_rows, upper_kappa = geometry(UPPER)
    assert all(int(r['ball_union_cosh'][1]) < source_lower.lo for r in lower_rows)
    far = next(r for r in upper_rows if r['pair'] == [3, 7])
    assert int(far['center_cosh'][0]) > source_upper.hi
    paths = [source, base, cap, f'scripts/{TAG}.py',
             'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
             'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py']
    return dict(status='PASS_EXACT_MONOTONE_PARAMETER_DIAMETER_GATE',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                interval_bits=BITS, interval_scale=str(SCALE), null_bracket=[str(LEFT), str(NULL_UPPER)],
                necessary_open_t_interval=[str(LOWER), str(UPPER)],
                source_vertices=41, source_vertex_pairs=len(distances),
                source_diameter_cosh_lower_pair=[32, 33], source_diameter_cosh_lower=source_lower.enc(),
                source_diameter_cosh_upper=source_upper.enc(), gamma_derivative_sign_witnesses=monotone,
                lower_endpoint_pairs=lower_rows, upper_endpoint_far_pair=far,
                lower_projected_norm_squares=A.enc(lower_kappa), upper_projected_norm_squares=A.enc(upper_kappa),
                scope='Exact necessary interval for compact Bugaenko H8 perpendicular attachment to all-dotted caps on the fixed H--E segment, with acuteness and finite volume. All center distances and radii increase because every gamma decreases and all old mutual products are nonpositive. The lower bound uses the diameter of the ball union; the upper uses the diameter of the forced center hull. No placement, angle compatibility, coverage or finite-volume example is asserted. Arithmetic reuses the frozen quartic and dyadic libraries.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = run(); out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['necessary_open_t_interval'])
