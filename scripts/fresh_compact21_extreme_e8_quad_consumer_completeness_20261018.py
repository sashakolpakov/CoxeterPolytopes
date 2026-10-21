#!/usr/bin/env python3
"""Exact bounded coefficients for the complete extreme E8 quad proof.

Only finite endpoint diagrams and sparse determinant identities are
enumerated. All-distance completeness is proved in the companion report.
"""
import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import product
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact21_d5_fork_quad_single_basic_unit_20261018 as B

P = B.P
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_extreme_e8_quad_consumer_completeness_20261018'


def gram(ids, edges):
    return [[P.K(2 if i == j else -1 if tuple(sorted((i,j))) in edges else 0)
             for j in ids] for i in ids]


def unit_interval(phi, prime, lo, hi):
    def power(r):
        value = P.K(1)
        for _ in range(abs(r)):
            value *= phi if r > 0 else 1/phi
        return value
    answer = []
    for a in ([0] if prime == 1 else [0,1]):
        # These exact tail inequalities, with monotonicity by phi^2,
        # exclude every odd exponent outside [-3,3], not just a sample.
        assert (prime**a*power(5)-hi).sign() > 0
        assert (prime**a*power(-5)-lo).sign() < 0
        for r in (-3,-1,1,3):
            conjugate_magnitude = prime**a*power(-r)
            if (conjugate_magnitude-lo).sign() >= 0 and (hi-conjugate_magnitude).sign() >= 0:
                answer.append(prime**a*power(r))
    return answer


def run():
    e8_edges = {(i,i+1) for i in range(6)}|{(2,7)}
    d5_edges = {(0,1),(1,2),(2,3),(2,4)}
    e8 = gram(list(range(8)), e8_edges)
    inverse = P.inverse(e8)
    entries = {(0,0):4,(0,1):7,(0,7):5,(1,1):14,(7,7):8,
               (5,5):6,(5,6):3,(6,6):2,(0,6):2}
    assert all(inverse[i][j] == n for (i,j),n in entries.items())
    old = P.schur(e8,[1,2,3,4,5,7],[0,6])
    assert old == [[P.K(Q(1,2)),P.K(Q(-1,2))],[P.K(Q(-1,2)),P.K(1)]]
    endpoints = {}
    for missing, marked, expected in ((0,6,[None,1,7]),(6,0,[None,5])):
        ids = [i for i in range(8) if i != missing]+[8]
        accepted = []
        for contact in [None]+[i for i in range(8) if i != missing]:
            edges = e8_edges|({(contact,8)} if contact is not None else set())
            g = gram(ids,edges)
            if P.psd_status(g) == 1 and not B.basic_payment(g,ids.index(marked),ids.index(8)):
                accepted.append(contact)
        assert accepted == expected
        endpoints[str(missing)] = accepted
    attachments = []
    for i in range(5):
        if P.psd_status(gram(list(range(6)),d5_edges|{(i,5)})) == 1:
            attachments.append(i)
    assert attachments == [0,3,4]
    profiles = []
    for c,d in product([None,0,3,4],repeat=2):
        edges = d5_edges|({(c,5)} if c is not None else set())|({(d,6)} if d is not None else set())
        if P.psd_status(gram(list(range(7)),edges)) == 1:
            profiles.append([c,d])
    assert len(profiles) == 11
    assert [r for r in profiles if None not in r] == [[0,3],[0,4],[3,0],[4,0]]

    x,y = B.Poly({(1,0):1}),B.Poly({(0,1):1})
    polys = {}
    for delta in (Q(0),Q(1),Q(5,4)):
        a = 2-delta
        rnl = [[Q(1,2)-x*(1/a),Q(-1,2),Q(-1,2)],
               [Q(-1,2),1,y-1],[Q(-1,2),y-1,1]]
        fnl = (2*x-a)*(-y)-a
        assert 2*a*B.det_poly(rnl) == (2-y)*fnl
        for name,u,factor,target in (
                ('S',x-Q(3,2),2-x,(1-x)*(y-a)-a/2),
                ('T',x-1,Q(3,2)-x,(1-2*x)*(y-a)-a)):
            rsn = [[Q(1,2),Q(-1,2),u],[Q(-1,2),1-y*(1/a),Q(-1,2)],
                   [u,Q(-1,2),Q(1,2)]]
            assert (a if name == 'S' else 2*a)*B.det_poly(rsn) == factor*target
        polys[f'N-L_delta={delta}'] = fnl.encode()
    for name,u,fac,target in (
            ('S-L',x-Q(3,2),2-x,(1-x)*(-y)-1),
            ('T-L',x-1,Q(3,2)-x,(1-2*x)*(-y)-2)):
        r = [[Q(1,2),Q(-1,2),u,Q(-1,2)],
             [Q(-1,2),1,Q(-1,2),y-1],
             [u,Q(-1,2),Q(1,2),Q(-1,2)],
             [Q(-1,2),y-1,Q(-1,2),1]]
        assert (1 if name == 'S-L' else 2)*B.det_poly(r) == fac*(2-y)*target
        polys[name] = target.encode()
    # C=0 factorization of the whole-old-vertex height determinant.
    for dc,dd in ((0,0),(1,0),(Q(5,4),0),(0,1),(0,Q(5,4))):
        a,b = Q(2)-dc,Q(2)-dd
        height_det = (4*x-a)*(2*y-b)-4*x*y
        assert height_det == 2*((2*x-a)*(y-b)-a*b/2)
    # Both-contact branch: use y=t and substitute Y=t^2/X as a
    # Laurent polynomial, retaining exact identity rather than samples.
    invx = B.Poly({(-1,0):1})
    yy = y*y*invx
    for dc,dd,expected in (
            (Q(1),Q(5,4),8*y*y-4*y+1-6*x-4*yy),
            (Q(5,4),Q(1),8*y*y-4*y+1-8*x-3*yy)):
        a,b = 2-dc,2-dd
        height_det = (4*x-a)*(2*yy-b)-(2*y+Q(1,2))**2
        assert 2*height_det == expected

    s = P.K((0,0,1,0))
    phi = (1+s)/2
    intervals = []
    for prime,lo,hi,expected in (
            (1,Q(1,2),Q(1),[phi]),(1,Q(1),Q(2),[1/phi]),
            (2,Q(1),Q(2),[1/phi,2*phi]),
            (3,Q(3,2),Q(3),[1/phi,3*phi]),
            (3,Q(1),Q(2),[1/phi,3*phi])):
        found = unit_interval(phi,prime,lo,hi)
        assert set(tuple(v.v) for v in found) == set(tuple(v.v) for v in expected)
        intervals.append({'product':prime,'minus_conjugate_interval':[str(lo),str(hi)],
                          'complete_U_values':P.enc(found),'odd_exponent_tails_checked':[-5,5]})
    bad_x = [(1+phi)/2,(3+1/phi)/8,(3+3*phi)/8]
    assert all(not B.integral(v) for v in bad_x)
    bad_alpha = [(1-2*phi)/2,(1-3*phi)/2]
    assert all(not B.integral(v) for v in bad_alpha)
    assert ((1-1/phi)/2).sign() > 0
    disc = 20+8*s
    obstruction = B.certify_nonsquare(disc)
    native = [P.K(1),P.K(2),phi*phi,P.K(3),(5+s)/2]
    assert all((B.sigma(v)-Q(3,16)).sign() > 0 for v in native)
    only_y = [v for v in native if (v-1).sign() > 0 and (Q(1,2)-B.sigma(v)).sign() >= 0]
    assert only_y == [phi*phi]
    Y = phi*phi
    aa,bb,cc = 4-4/Y,P.K(-2),Q(1,2)-Q(3,2)*Y
    assert bb*bb-4*aa*cc == disc
    assert (phi*phi-1)*(phi*phi-2) == 1
    sources = [f'scripts/{TAG}.py',
               'scripts/fresh_compact21_d5_fork_quad_single_basic_unit_20261018.py',
               'scripts/fresh_compact29_small_faces_20261015.py',
               'results/fresh_compact21_d5_two_e8_triangle_coupling_20261018.json',
               'results/fresh_compact21_single_field_20261016.json']
    pins = {}
    for path in sources:
        with (ROOT/path).open('rb') as stream:
            pins[path] = file_digest(stream,'sha256').hexdigest()
    return {'status':'PASS_COMPLETE_EXTREME_E8_QUAD_CONSUMER_COEFFICIENTS',
            'sources':pins,'E8_inverse_entries':{str(k):v for k,v in entries.items()},
            'zero_basic_endpoint_contacts_by_deleted_E8_root':endpoints,
            'spherical_D5_single_attachments':attachments,'spherical_D5_double_profiles':profiles,
            'corank_polynomials':polys,'complete_short_unit_intervals':intervals,
            'nonintegral_X_values':P.enc(bad_x),'nonintegral_alpha_values':P.enc(bad_alpha),
            'two_contact_final_nonsquare_discriminant':obstruction,
            'remaining_opposite_squares':P.enc([phi*phi,phi*phi]),
            'all_opposite_distances_covered_in_written_proof':True,
            'specific_E8_pair':[0,6],'all_E8_long_pairs_classified':False,
            'new_forest_budget_rows':['old_total-1<=z','old_total-2*b_06<=z'],
            'old_model_optimized':False,'new_dimension_bound':None}


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
            assert json.loads(target.read_text()) == json.loads(json.dumps(result))
        else:
            assert not target.exists()
            target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'D5_profiles':11,'peak_rss_bytes':rss(),
                      'new_dimension_bound':None}))
