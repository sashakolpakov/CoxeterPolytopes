#!/usr/bin/env python3
"""Numerical necessary-diameter diagnostics along one dotted-cap segment.

The sampled trends and bisection crossings are NOT interval certificates.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import fresh_h9_mixed_d4_coupling_20261014_review as A
import fresh_h9_h8_facet_cover_20261014_geometry as B
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_facet_cover_20261014_parameter_scout'

def run():
    mp.mp.dps = 80
    base = B.data()
    source_diameter = max(base['D'])
    D = json.loads((ROOT/'results/fresh_h9_mixed_d4_coupling_20261014.json').read_text())
    C = json.loads((ROOT/'results/fresh_h9_mixed_d4_global_caps_20261014.json').read_text())
    Q = A.dm(D['form']); N = A.dm(D['final_normals']); H = A.dv(D['center'])
    E = [A.E(A.F(x)) for x in C['literal_escape']]
    def f(x): return mp.mpf(x.numerator)/x.denominator
    def ev(x):
        return f(x.a.a)+f(x.a.b)*mp.sqrt(5)+(f(x.b.a)+f(x.b.b)*mp.sqrt(5))*mp.sqrt((16+4*mp.sqrt(5))/11)
    G = mp.matrix([[ev(A.pair(x, y, Q)) for y in N] for x in N])
    eh, hh, ee = [ev(A.pair(x, y, Q)) for x, y in ((E, H), (H, H), (E, E))]
    gh = mp.matrix([ev(A.pair(H, n, Q)) for n in N]); ge = mp.matrix([ev(A.pair(E, n, Q)) for n in N])
    a, b, c = hh+ee-2*eh, 2*(eh-hh), hh
    def roots(a, b, c):
        d = b*b-4*a*c
        return [(-b-mp.sqrt(d))/(2*a), (-b+mp.sqrt(d))/(2*a)] if d >= 0 else []
    t0 = min(x for x in roots(a, b, c) if 0 < x < 1)
    bounds = []
    for i in range(12):
        d, nn = ge[i]-gh[i], G[i, i]
        bounds.extend((x, i) for x in roots(d*d-nn*a, 2*gh[i]*d-nn*b, gh[i]**2-nn*c) if t0 < x < 1)
    t1, wall = min(bounds)
    def geom(u, full=False):
        t = t0+u*(t1-t0); q = a*t*t+b*t+c; g = gh+t*(ge-gh)
        PP = G-g*g.T/q; k = [mp.sqrt(-PP[i, i]) for i in range(12)]
        R = [mp.asinh(mp.sqrt(G[i, i])/k[i]) for i in range(12)]
        ds = [[mp.acosh(max(1, -PP[i, j]/(k[i]*k[j]))) for j in range(12)] for i in range(12)]
        cd, ci = max((ds[i][j], (i, j)) for i in range(12) for j in range(i, 12))
        ud, ui = max((ds[i][j]+R[i]+R[j], (i, j)) for i in range(12) for j in range(i, 12))
        row = dict(u=mp.nstr(u, 30), t=mp.nstr(t, 60), center_diameter=mp.nstr(cd, 40),
                   center_pair=ci, union_diameter=mp.nstr(ud, 40), union_pair=ui,
                   maximum_radius=mp.nstr(max(R), 40), passes_necessary_diameters=cd <= source_diameter <= ud)
        if full: row.update(radii=[mp.nstr(x, 40) for x in R])
        return row, cd, ud
    samples = [geom(mp.mpf(i)/1000)[0] for i in range(1, 1000)]
    def bisect(index):
        lo, hi = mp.mpf('1e-12'), 1-mp.mpf('1e-12')
        for _ in range(100):
            mid = (lo+hi)/2
            if geom(mid)[index] < source_diameter: lo = mid
            else: hi = mid
        return geom((lo+hi)/2, True)[0]
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_h9_h8_facet_cover_20261014_geometry.py',
             'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
             'data/known_examples/bugaenko_compact_h8.json',
             'results/fresh_h9_mixed_d4_coupling_20261014.json',
             'results/fresh_h9_mixed_d4_global_caps_20261014.json']
    return dict(status='NUMERICAL_PARAMETER_DIAMETER_DIAGNOSTIC_ONLY',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, precision_digits=80,
                source_vertex_count=41, source_diameter=mp.nstr(source_diameter, 60),
                t_null=mp.nstr(t0, 70), t_first_parallel=mp.nstr(t1, 70), parallel_wall=wall,
                literal_cap=geom((mp.mpf(C['interpolation_parameter'])-t0)/(t1-t0), True)[0],
                lower_numeric_crossing=bisect(2), upper_numeric_crossing=bisect(1),
                samples=samples, passing_samples=sum(r['passes_necessary_diameters'] for r in samples),
                scope='Only one H--E line and one compact H8 source. The two inequalities are necessary geometric conditions, but the999 samples and two bisection crossings are numerical diagnostics. No monotonicity theorem, uniform exclusion, valid placement, Coxeter compatibility or ball coverage is claimed.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = json.loads(json.dumps(run())); out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['passing_samples'], result['lower_numeric_crossing']['u'], result['upper_numeric_crossing']['u'])
