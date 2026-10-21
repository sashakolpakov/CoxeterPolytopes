#!/usr/bin/env python3
"""Exclude the remaining mixed all-D5 quad without any half-weight assumption."""
import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact21_d5_fork_quad_single_basic_unit_20261018 as B

P = B.P
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_d5_mixed_fork_quad_exclusion_20261018'


def run():
    x,y = B.Poly({(1,0):1}),B.Poly({(0,1):1})
    r = [[Q(5,4),Q(-3,4),x-Q(3,4),Q(-1,2)],
         [Q(-3,4),Q(5,4),Q(-3,4),y-Q(1,2)],
         [x-Q(3,4),Q(-3,4),Q(5,4),Q(-1,2)],
         [Q(-1,2),y-Q(1,2),Q(-1,2),1]]
    f = (2*x+1)*(1+y-y*y)+3*y-5
    assert 2*B.det_poly(r) == (2-x)*f
    s = P.K((0,0,1,0))
    phi = (1+s)/2
    cases = []
    for yy,xx in ((P.K(-1),P.K(Q(-9,2))),(-phi,-5*phi/4)):
        assert (2*xx+1)*(1+yy-yy*yy)+3*yy-5 == 0
        assert not B.integral(xx)
        cases.append(P.enc([yy,xx]))
    disc = 34+14*s
    assert (3-s)*(3-s)+4*s*(5+s) == disc
    obstructions = [B.certify_nonsquare(v) for v in (P.K(28),disc)]
    sources = [f'scripts/{TAG}.py',
               'scripts/fresh_compact21_d5_fork_quad_single_basic_unit_20261018.py',
               'scripts/fresh_compact29_small_faces_20261015.py',
               'results/fresh_compact21_d5_fork_quad_triangle_priority_20261018.json',
               'results/fresh_compact21_single_field_20261016.json']
    pins = {}
    for path in sources:
        with (ROOT/path).open('rb') as stream:
            pins[path] = file_digest(stream,'sha256').hexdigest()
    return {'status':'PASS_MIXED_ALL_D5_QUAD_EXCLUDED_NO_HALF_WEIGHTS',
            'sources':pins,'corank_polynomial':f.encode(),
            'finite_y_nonintegral_x':cases,'finite_x_nonsquare_discriminants':obstructions,
            'D5_fork_priority_rule_needs_half_units':False,
            'all_distance_exclusion_proved_in_companion_report':True,
            'new_dimension_bound':None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    signal.alarm(59)
    rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert rss() < 1024**3
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        result = run()
        assert rss() < 1024**3
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == result
        else:
            assert not target.exists()
            target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'peak_rss_bytes':rss(),
                      'new_dimension_bound':None}))
