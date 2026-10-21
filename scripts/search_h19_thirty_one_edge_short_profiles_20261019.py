#!/usr/bin/env python3
"""One bounded62-profile short-root search at the CURRENT31-wall H19 edge.

Coordinator only: bounded parent-lock wrapper, --write or --verify.
Only pure arithmetic helpers are imported; no old input/check/run is called.
No optimizer, unrestricted queue, subprocess or inner lock is used.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import signal

from search_h19_omega2_finite_profiles_20261019 import (
    Q, native as N, decode, dot, mv, digest, field_sqrt, pair_kind)
from check_arithmetic_dn_opposed_fork_models_20261019 import inverse, inertia

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_thirty_one_edge_short_profiles_20261019.json'
PINS = {
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'results/h19_twenty_six_wall_seed_20261019.json':
    'e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b',
    'results/h19_twenty_seven_wall_seed_20261019.json':
    'b5b529b485a01b9315be6f31bbb122c170652273d3932b1eecece85680b5ce7d',
    'results/h19_twenty_nine_wall_seed_20261019.json':
    'd6e7129224a47621c30737d70c87b713b818e8b4154e0b8de2ea5f3ac49cf031',
    'results/h19_thirty_wall_seed_20261019.json':
    'd308740335fabc8cae80224670d777d24d8bcfb2c3484fec4b0ed0bb48758176',
    'results/h19_thirty_one_wall_seed_20261019.json':
    '5a3b25aa139793f304d535f56d6033ad103d88cffe07c20278bc7b42c2712d57',
    'results/h19_omega3_finite_profiles_20261019.json':
    'aa85c4c7c5373f66f124702563a55f9c88a7b7a28fe962f78db5bbe01af54c1d',
    'scripts/check_h19_thirty_one_wall_seed_20261019.py':
    'e81cd38a68c5262acc6421c16f8d8114e2480f5f453edde503a32de30710c76e',
    'scripts/search_h19_omega2_finite_profiles_20261019.py':
    '02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'results/h19_thirty_one_escape_finite_profile_reduction_20261019.md':
    '2be3d769f288ab25bfdb5831115c3002c1f8022890d30b4940d527ed4c01336c',
}


def enc(vector):
    return [v.encoded() for v in vector]


def add(left, right, scale=1):
    return [a+scale*b for a, b in zip(left, right)]


def fixed_inputs():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    read = lambda stem: json.loads((ROOT/'results'/(stem+'_20261019.json')).read_text())
    old, qsave, jsave, ksave, current, hsave = map(read, (
        'h19_twenty_six_wall_seed', 'h19_twenty_seven_wall_seed',
        'h19_twenty_nine_wall_seed', 'h19_thirty_wall_seed',
        'h19_thirty_one_wall_seed', 'h19_omega3_finite_profiles'))
    assert current['status'] == 'PASS_EXACT_31_WALL_AND_BOUNDED_EDGE_DIAGNOSTIC'
    assert current['checker_sha256'] == PINS['scripts/check_h19_thirty_one_wall_seed_20261019.py']
    assert current['walls'] == 31 and current['wall_pairs'] == 465
    for record in (old, qsave, jsave, ksave, current, hsave):
        assert record['status'].startswith('PASS')
        for name, expected in record['sources'].items():
            assert digest(ROOT/name) == expected, name
    phi, D = N(0, 1), N(1, 5)
    u, t = N(10, 17), N(12, 19)
    assert 7*t*t-D*u*u == 1
    G = [[Q() for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, phi, phi, Q(2)
    for i in range(8):
        G[i+2][i+2] = Q(2)
    for i, j in ((0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)):
        G[i+2][j+2] = G[j+2][i+2] = Q(-1)
    assert inertia(G) == [9,1,0] and inertia(G, True) == [10,0,0]
    gi = inverse(G, Q)
    ci = [row[2:] for row in gi[2:]]
    assert ci[1] == list(map(Q, (7,14,20,16,12,8,4,10)))
    assert ci[2] == list(map(Q, (10,20,30,24,18,12,6,15)))
    prefix = json.loads((ROOT/'data/candidates/h19_vinberg_prefix.json').read_text())
    roots = {f'old{i}': [N(v['a'], v['b']) for v in row]
        for i, row in enumerate(prefix['output']['roots_in_vinberg_order'])
        if i not in (10,20,21)}
    for name, cv in old['new_covectors'].items():
        roots[name] = mv(gi, [decode(v) for v in cv])
    roots['q'] = [decode(v) for v in qsave['new_root']]
    assert len(hsave['survivors']) == 1
    roots['h'] = [decode(v) for v in hsave['survivors'][0]['coordinates']]
    for name, record in (('j', jsave), ('K', ksave), ('L', current)):
        roots[name] = [decode(v) for v in record['new_root']]
    assert len(roots) == 31
    cov = {n: mv(G, v) for n, v in roots.items()}
    norms = {n: dot(v, cov[n]) for n, v in roots.items()}
    for name, root in roots.items():
        expected = Q(2) if name in ('q','h','j','K','L') else decode(old['norms'][name])
        assert norms[name] == expected and expected.sign() > 0
        assert all(v.integral() for v in root)
        assert all((2*v/expected).integral() for v in cov[name])
    pairs, categories = [], Counter()
    for a, b in combinations(roots, 2):
        value = dot(roots[a], cov[b])
        if norms[a] == 2:
            kind = pair_kind(value, norms[b], phi)
        elif norms[b] == 2:
            kind = pair_kind(value, norms[a], phi)
        else:
            assert value == 0 or (value.sign() < 0 and (value*value-norms[a]*norms[b]).sign() > 0)
            kind = 'right' if value == 0 else 'ultraparallel'
        assert not kind.startswith('reject')
        pairs.append([a,b,value.encoded(),kind])
        categories[kind] += 1
    pair_hash = hashlib.sha256(json.dumps(pairs,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert len(pairs) == 465 and pair_hash == current['wall_pairs_sha256']
    assert dict(categories) == current['wall_pair_categories']
    center = [decode(v) for v in current['strict_center']]
    square = lambda v: dot(v, mv(G, v))
    assert square(center).sign() < 0 and all(dot(center,v).sign() < 0 for v in cov.values())
    V = [Q(1),-phi/2]+[D*u*v/(14*t) for v in ci[1]]
    delta = [Q(),Q()]+[a-F(10,7)*b for a,b in zip(ci[2],ci[1])]
    horizon = Q(F(-42,220), F(27,220))
    E = add(V, delta, horizon)
    witness = current['bounded_edge_walk']['witness']
    assert witness['from'] == 3 and witness['freed'] == 'old3' and witness['first_hit'] == ['q']
    assert horizon == decode(witness['parameter']) and E == [decode(v) for v in witness['endpoint']]
    assert V == [decode(v) for v in next(r for r in current['bounded_edge_walk']['vertices'] if r['id']==3)['point']]
    assert square(V) == -D/(14*t*t) and square(delta) == F(10,7)
    assert dot(V,mv(G,delta)) == 0 and dot(delta,cov['old3']) == -1
    assert square(E) == Q(F(6227,4840),F(-691,1210)) and square(E).sign() > 0
    white = ['old0','old1']+[f'old{i}' for i in range(4,9)]+['s']
    assert {n for n in roots if dot(V,cov[n]) == 0} == set(white)|{'old3'}
    assert {n for n in roots if dot(V,cov[n]) == 0 and dot(delta,cov[n]) == 0} == set(white)
    assert all(dot(v,cv).sign() <= 0 for v in (V,E) for cv in cov.values())
    assert dot(V,mv(G,center)).sign() < 0
    assert inertia([[dot(roots[a],cov[b]) for b in white] for a in white]) == [8,0,0]
    for j in range(8):
        assert roots[f'old{j+1}'] == [Q(),Q()]+[Q(-int(i==j)) for i in range(8)]
    assert roots['old0'] == [Q(),Q(-1)]+[Q()]*8
    assert cov['s'] == [-D*u,Q(),Q(),t]+[Q()]*6
    return phi,D,u,t,G,ci,roots,cov,norms,center,V,delta,E,horizon,white,current


def finite_profiles(phi, ci, roots, cov, white):
    raw = []
    # The four isolated white roots; -1 means binary old0, -2 means s.
    isolates = (-1,0,7,-2)

    def insert(contacts, family):
        beta, epsilon, eta = Q(), Q(), [Q()]*8
        for node, value in contacts.items():
            if node == -1:
                beta = value
            elif node == -2:
                epsilon = value
            else:
                eta[node] = value
        raw.append((beta,epsilon,eta,family))

    insert({}, 'isolated')
    for i in isolates:
        for value in (Q(1),phi):
            insert({i:value}, 'one_isolate')
    for i,j in combinations(isolates,2):
        for a,b in ((Q(1),Q(1)),(Q(1),phi),(phi,Q(1))):
            insert({i:a,j:b}, 'two_isolates')
    for nodes in combinations(isolates,3):
        insert({i:Q(1) for i in nodes}, 'three_isolates')
    assert len(raw) == 31
    for node in (3,6,4,5):
        max_isolates = 2 if node in (3,6) else 1
        for k in range(max_isolates+1):
            for nodes in combinations(isolates,k):
                insert({node:Q(1)}|{i:Q(1) for i in nodes}, 'touches_A4')
    assert len(raw) == 63
    summaries, selected, seen = [], [], set()
    W = [[dot(roots[a],cov[b]) for b in white] for a in white]
    for beta,epsilon,eta,family in raw:
        code = tuple(tuple(v.encoded()) for v in [beta,epsilon]+eta)
        assert code not in seen
        seen.add(code)
        q = (beta*beta+epsilon*epsilon+eta[0]*eta[0]+eta[7]*eta[7])/2
        chain = eta[3:7]
        ai = [[Q(F(v,5)) for v in row] for row in ((4,3,2,1),(3,6,4,2),(2,4,6,3),(1,2,3,4))]
        q += dot(chain,mv(ai,chain))
        assert (q-2).sign() < 0
        contact_vector = [beta,eta[0]]+eta[3:8]+[epsilon]
        extension = [row+[-contact_vector[i]] for i,row in enumerate(W)]
        extension.append([-v for v in contact_vector]+[Q(2)])
        assert inertia(extension) == [9,0,0]
        R1,R2 = dot(ci[1],eta),dot(ci[2],eta)
        Q0 = dot(eta,mv(ci,eta))
        B = 7*R2-10*R1
        K = B*B+140-35*(beta*beta+epsilon*epsilon)-70*Q0+5*R1*R1
        assert K == 70*(2-q) and K.sign() > 0 and K.integral()
        ident = f'current31_profile_{len(summaries):02d}'
        zero = all(v == 0 for v in [beta,epsilon]+eta)
        row = dict(id=ident,family=family,beta=beta.encoded(),epsilon=epsilon.encoded(),
            eta=enc(eta),white_projection_square=q.encoded(),R1=R1.encoded(),R2=R2.encoded(),
            Q0=Q0.encoded(),B=B.encoded(),K=K.encoded(),excluded_zero_profile_mod5=zero)
        summaries.append(row)
        if not zero:
            selected.append((ident,beta,epsilon,eta,R1,B,K))
    assert len(selected) == 62
    assert 3 not in {i*i%5 for i in range(5)}
    return selected,summaries


def run():
    phi,D,u,t,G,ci,roots,cov,norms,oldcenter,V,delta,E,horizon,white,current = fixed_inputs()
    selected,summaries = finite_profiles(phi,ci,roots,cov,white)
    counts, xs = Counter(), []
    for b in range(301):
        m = (isqrt(5*b*b)-b)//2
        for a in range(m-2,m+4):
            counts['x_coefficient_pairs'] += 1
            x = N(a,b)
            if (x.sign() > 0 and (x-600).sign() <= 0
                    and (x.conjugate()+2).sign() >= 0 and (x.conjugate()-2).sign() <= 0):
                xs.append((a,b,x))
    assert counts['x_coefficient_pairs'] == 1806 and len(xs) == 1073
    counts['x_in_requested_domain'] = len(xs)
    square = lambda v: dot(v,mv(G,v))
    rejected, cuts, survivors, outcome_hash = [], [], [], hashlib.sha256()
    retained_bytes = 0

    def retain(target, record):
        nonlocal retained_bytes
        retained_bytes += len(json.dumps(record,sort_keys=True).encode())
        assert retained_bytes < 850000, 'Retained-certificate cap reached; no completion claim'
        target.append(record)

    for ident,beta,epsilon,eta,R1,B,K in selected:
        for aa,bb,x in xs:
            counts['profile_x_pairs'] += 1
            n = (x+u*epsilon)/t
            y = (beta-phi*x)/2
            if not n.integral():
                counts['reject_nonintegral_n'] += 1
                continue
            if not y.integral():
                counts['reject_nonintegral_y'] += 1
                continue
            counts['field_square_tests'] += 1
            positive_A = field_sqrt(5*D*n*n+K)
            if positive_A is None:
                continue
            counts['field_square_solutions'] += 1
            for A in (positive_A,-positive_A):
                counts['signed_A_branches'] += 1
                U = (A-B)/10
                T = (D*u*n-7*t*epsilon-20*U-R1)/14
                if not U.integral() or not T.integral():
                    counts['reject_nonintegral_T_U'] += 1
                    continue
                z = mv(ci,[eta[j]+(T if j==1 else U if j==2 else Q()) for j in range(8)])
                root = [x,y]+z
                cv = mv(G,root)
                assert all(v.integral() for v in root+cv)
                assert dot(root,cv) == 2 and A*A-5*D*n*n == K
                assert all((2*v/2).integral() for v in cv)
                assert dot(root,cov['old0']) == -beta and dot(root,cov['s']) == -epsilon
                assert all(dot(root,cov[f'old{j+1}']) == -eta[j] for j in (0,3,4,5,6,7))
                counts['integral_primitive_norm_two_roots'] += 1
                products,invalid,categories = {},[],Counter()
                for name, old in roots.items():
                    value = dot(old,cv)
                    kind = pair_kind(value,norms[name],phi)
                    products[name] = [value.encoded(),kind]
                    categories[kind] += 1
                    if kind.startswith('reject'):
                        invalid.append(name)
                    elif norms[name] != 2:
                        assert kind in ('right','ultraparallel'), 'Forbidden finite mixed bond'
                assert sum(categories.values()) == 31
                vp,dp,cut,cp = dot(V,cv),dot(delta,cv),dot(E,cv),dot(oldcenter,cv)
                assert vp == -D*(x+u*epsilon)/(14*t*t) and vp.sign() < 0
                assert dp == (10*U+8*eta[3]+6*eta[4]+4*eta[5]+2*eta[6]+5*eta[7])/7
                reasons = ([] if not invalid else ['incompatible_native_pairs'])
                if cut.sign() <= 0:
                    reasons.append('does_not_cut_current_escape')
                record = dict(profile=ident,x_phi_coefficients=[aa,bb],n=n.encoded(),A=A.encoded(),
                    U=U.encoded(),T=T.encoded(),y=y.encoded(),V_product=vp.encoded(),
                    edge_slope=dp.encoded(),escape_product=cut.encoded(),old_center_product=cp.encoded(),
                    invalid_walls=invalid,pair_categories=dict(categories),rejected=reasons)
                outcome_hash.update((json.dumps(record,sort_keys=True,separators=(',',':'))+'\n').encode())
                if reasons:
                    if len(rejected) < 128:
                        retain(rejected,record)
                    else:
                        counts['rejected_records_omitted_from_preview'] += 1
                    continue
                assert dp.sign() > 0
                tau = -vp/dp
                endpoint = add(V,delta,tau)
                ep = {name:dot(endpoint,oldcv) for name,oldcv in cov.items()}
                finite = square(endpoint).sign() < 0
                assert tau.sign() > 0 and (horizon-tau).sign() > 0
                assert all(v.sign() <= 0 for v in ep.values()) and dot(endpoint,cv) == 0
                if not finite:
                    counts['compatible_cuts_with_nonfinite_endpoint'] += 1
                    retain(cuts,record|dict(endpoint_parameter=tau.encoded(),endpoint_square=square(endpoint).encoded(),
                        outcome='compatible_escape_cut_but_not_a_finite_endpoint'))
                    continue
                assert (D/(20*t*t)-tau*tau).sign() > 0
                assert dot(V,mv(G,endpoint)).sign() < 0
                assert {name for name,value in ep.items() if value == 0} == set(white)
                active_roots = [roots[name] for name in white]+[root]
                assert inertia([[dot(a,mv(G,b)) for b in active_roots] for a in active_roots]) == [9,0,0]
                center_scale = F(1)
                while (vp+center_scale*cp).sign() >= 0:
                    center_scale /= 2
                center = add(V,oldcenter,center_scale)
                assert square(center).sign() < 0 and dot(V,mv(G,center)).sign() < 0
                assert dot(center,cv).sign() < 0 and all(dot(center,oldcv).sign() < 0 for oldcv in cov.values())
                witness = add(center,root,-dot(center,cv)/2)
                assert square(witness).sign() < 0 and dot(witness,cv) == 0
                assert all(dot(witness,oldcv).sign() < 0 for oldcv in cov.values())
                retention = {name:dot([decode(v) for v in saved['point']],cv)
                    for name,saved in current['finite_vertex_witnesses'].items()}
                assert len(retention) == 7
                caps = {}
                for name in ('r','w'):
                    cap_point = add(center,roots[name],-dot(center,cov[name])/norms[name])
                    assert square(cap_point).sign() < 0 and dot(cap_point,cov[name]) == 0
                    caps[name] = dict(pair_kind=products[name][1],product=products[name][0],
                        strict_cap_point_product=dot(cap_point,cv).encoded(),
                        full_cap_preservation_certified=products[name][1]=='ultraparallel' and dot(cap_point,cv).sign()<0)
                named_retained = all(v.sign() < 0 for v in retention.values())
                caps_certified = all(v['full_cap_preservation_certified'] for v in caps.values())
                all_retained = named_retained and caps_certified
                full = record|dict(coordinates=enc(root),covector=enc(cv),wall_products=products,
                    endpoint=enc(endpoint),endpoint_parameter=tau.encoded(),endpoint_square=square(endpoint).encoded(),
                    endpoint_old_products={name:value.encoded() for name,value in ep.items()},
                    strict_center=enc(center),center_scale=str(center_scale),new_facet_witness=enc(witness),
                    old_named_vertex_products={name:value.encoded() for name,value in retention.items()},
                    old_cap_facets=caps,all_seven_named_vertices_strictly_retained=named_retained,
                    both_full_caps_certified=caps_certified,required_retention_certified=all_retained,
                    primitive_reason='Integral norm2 in an even O-lattice: a nonunit coordinate divisor would divide1 after squaring.',
                    outcome='compatible_finite_endpoint_all_retained' if all_retained else 'compatible_finite_endpoint_retention_not_certified')
                counts['compatible_finite_endpoint_roots'] += 1
                if all_retained:
                    retain(survivors,full)
                else:
                    retain(cuts,full)
    assert counts['profile_x_pairs'] == 62*len(xs)
    assert all(digest(ROOT/name)==expected for name,expected in PINS.items())
    return dict(status='PASS_BOUNDED_CURRENT31_62_SHORT_PROFILE_SEARCH',source_sha256=digest(Path(__file__)),
        sources=PINS,old_walls=31,fixed_white_walls=white,old_escape_source=dict(vertex=3,freed='old3',first_hit='q'),
        domain=dict(x_positive_at_most=600,absolute_conjugate_x_at_most=2,b_inclusive=[0,300],
            a_rule='m-2..m+3; m=(isqrt(5*b*b)-b)//2',
            equation='A^2-5D*((x+u*epsilon)/t)^2=K; both A signs; T,U,y integral'),
        profiles=summaries,finite_short_profiles=63,retained_short_profiles=62,
        zero_short_profile_all_height_excluded_mod5=True,long_profile_searched=False,
        stage_counts=dict(counts),integral_outcomes_sha256=outcome_hash.hexdigest(),
        rejected_integral_root_preview=rejected,other_compatible_cuts=cuts,survivors=survivors,
        global_finite_volume_proved=False,
        scope='All62 necessary nonzero SHORT finite profiles at the ACTUAL31-wall V/E edge in the declared x box. '
              'Both field-square signs, integral primitive norm2 roots, ALL31 signed native pairs, actual finite endpoint, '
              'seven named vertex retention and two cap orientations are checked. Old-center positivity is never a rejection. '
              'Full-cap certification here requires an ultraparallel pair and a strict cap point; right-cap cases are '
              'retained as compatible cuts with cap retention uncertified. Rejected-record preview is capped; '
              'the hash binds every integral outcome. No long search, all-height census, '
              'simultaneous choice of survivors,32-wall check, full volume, reflectivity or dimension theorem.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    signal.alarm(45)
    result = run()
    payload = json.dumps(result,sort_keys=True,indent=2)+'\n'
    assert len(payload.encode()) < 1024**2
    if args.verify:
        assert json.loads(CERT.read_text()) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream:
            stream.write(payload)
    signal.alarm(0)
    print(json.dumps(dict(status=result['status'],counts=result['stage_counts'],survivors=len(result['survivors']))))


if __name__ == '__main__':
    main()
