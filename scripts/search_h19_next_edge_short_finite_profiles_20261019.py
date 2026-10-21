#!/usr/bin/env python3
"""Bounded11-profile norm2 search on the fixed28-wall H19 edge.

Coordinator only: run_bounded_research_check_20261018.py --parent-lock --
THIS --write or --verify. No inner lock or imported producer/checker body.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import signal
from search_h19_omega2_finite_profiles_20261019 import (
    Q, native, decode, dot, mv, digest, field_sqrt, pair_kind)
from check_arithmetic_dn_opposed_fork_models_20261019 import inverse, inertia

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_next_edge_short_finite_profiles_20261019.json'
PINS = {
    'scripts/search_h19_omega2_finite_profiles_20261019.py':
    '02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'results/h19_twenty_six_wall_seed_20261019.json':
    'e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b',
    'results/h19_twenty_seven_wall_seed_20261019.json':
    'b5b529b485a01b9315be6f31bbb122c170652273d3932b1eecece85680b5ce7d',
    'results/h19_twenty_eight_wall_seed_20261019.json':
    '3996be94b07f251685a15650e50321456a4f3ed5548e0ad927d94abba1525b83',
    'results/h19_omega3_finite_profiles_20261019.json':
    'aa85c4c7c5373f66f124702563a55f9c88a7b7a28fe962f78db5bbe01af54c1d',
    'scripts/check_h19_twenty_eight_wall_seed_20261019.py':
    'f9ccdb1dbb5fe6a810bf456be56b3f15d1b83be69fd0d080931713eae3738600',
    'results/h19_next_edge_sixteen_profiles_research_20261019.md':
    '6dc87d20c5cf207d90aabe28d3f1aa1b8c38596484f8fd6ab98c6532d7ff95cb',
    'results/review_h19_next_edge_long_profiles_mod19_finite57_20261019.md':
    '0aec792ba6c29c852ed95b1eeec56ddd3fb4ff2b606dbfd95f4ce6b08673cbda',
    'results/h19_twenty_eight_wall_seed_execution_and_escape_20261019.md':
    'f0087b974aac69c758747cec515440888735bc1ac214e10217e858e8ec3b29f7',
}


def add(left, right, scale=1):
    return [a+scale*b for a, b in zip(left, right)]


def enc(vector):
    return [v.encoded() for v in vector]


def fixed_inputs():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    read = lambda stem: json.loads((ROOT/'results'/(stem+'_20261019.json')).read_text())
    old, middle, saved, found = map(read, ('h19_twenty_six_wall_seed',
        'h19_twenty_seven_wall_seed', 'h19_twenty_eight_wall_seed', 'h19_omega3_finite_profiles'))
    assert old['walls'] == 26 and middle['walls'] == 27 and saved['walls'] == 28
    assert saved['status'] == 'PASS_EXACT_28_WALL_AND_BOUNDED_EDGE_DIAGNOSTIC'
    assert saved['wall_pairs'] == 378 and saved['checker_sha256'] == PINS[
        'scripts/check_h19_twenty_eight_wall_seed_20261019.py']
    for record in (old, middle, saved, found):
        assert str(record['status']).startswith('PASS')
        for name, expected in record['sources'].items():
            assert digest(ROOT/name) == expected, name
    phi = native(0, 1)
    G = [[Q() for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, phi, phi, Q(2)
    for i in range(8): G[i+2][i+2] = Q(2)
    for i, j in [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)]:
        G[i+2][j+2] = G[j+2][i+2] = Q(-1)
    assert inertia(G) == [9,1,0] and inertia(G, True) == [10,0,0]
    gi = inverse(G, Q); ci = [row[2:] for row in gi[2:]]
    prefix = json.loads((ROOT/'data/candidates/h19_vinberg_prefix.json').read_text())
    roots = {f'old{i}': [native(v['a'],v['b']) for v in row]
             for i, row in enumerate(prefix['output']['roots_in_vinberg_order']) if i not in (10,20,21)}
    for name, covector in old['new_covectors'].items():
        roots[name] = mv(gi, [decode(v) for v in covector])
    roots['q'] = [decode(v) for v in middle['new_root']]
    assert len(found['survivors']) == 1
    h = found['survivors'][0]; assert h['x_phi_coefficients'] == [22,36]
    roots['h'] = [decode(v) for v in h['coordinates']]
    cov = {n: mv(G,v) for n,v in roots.items()}
    norms = {n: dot(v,cov[n]) for n,v in roots.items()}
    assert len(roots) == 28 and cov['h'] == [decode(v) for v in h['covector']]
    for name, root in roots.items():
        assert norms[name] == (2 if name in ('q','h') else decode(old['norms'][name]))
        assert all(v.integral() for v in root)
        assert all((2*v/norms[name]).integral() for v in cov[name])
    pairs, categories = [], Counter()
    for a, b in combinations(roots, 2):
        v = dot(roots[a],cov[b])
        if norms[a] == 2: kind = pair_kind(v,norms[b],phi)
        elif norms[b] == 2: kind = pair_kind(v,norms[a],phi)
        else:
            assert v == 0 or (v.sign() < 0 and (v*v-norms[a]*norms[b]).sign() > 0)
            kind = 'right' if v == 0 else 'ultraparallel'
        assert not kind.startswith('reject')
        pairs.append([a,b,v.encoded(),kind]); categories[kind] += 1
    code = json.dumps(pairs,sort_keys=True,separators=(',',':')).encode()
    assert hashlib.sha256(code).hexdigest() == saved['wall_pairs_sha256']
    assert dict(categories) == saved['wall_pair_categories']
    p = [Q(2),-phi]+[Q()]*8
    center = add(p,[Q(),Q(F(1,10))]+[sum(row,Q())/10 for row in ci])
    omega0, omega1 = ci[0], ci[1]
    V = [Q(1),-phi/2]+[phi*phi*phi*v/4 for v in omega0]
    direction = [Q(),Q()]+[v-F(7,4)*u for u,v in zip(omega0,omega1)]
    E = add(V,direction,F(1,5)); f = 2*phi-3
    square = lambda v: dot(v,mv(G,v))
    assert square(center).sign() < 0 and dot(p,mv(G,center)).sign() < 0
    assert all(dot(center,v).sign() < 0 for v in cov.values())
    assert square(V) == -f/4 and square(direction) == F(7,4)
    assert dot(V,mv(G,direction)) == 0 and square(E) == (41-25*phi)/50
    assert square(E).sign() > 0 and dot(p,mv(G,E)).sign() < 0
    assert all(dot(v,c).sign() <= 0 for v in (V,E) for c in cov.values())
    white = ['old0']+[f'old{j}' for j in range(3,9)]+['z']
    assert [n for n in roots if dot(V,cov[n]) == 0] == ['old0']+[f'old{j}' for j in range(2,9)]+['z']
    assert {n for n in roots if dot(direction,cov[n]) == 0 and dot(V,cov[n]) == 0} == set(white)
    assert inertia([[dot(roots[a],cov[b]) for b in white] for a in white]) == [8,0,0]
    assert dot(direction,cov['old2']) == -1 and dot(E,cov['old18']) == 0
    witness = saved['bounded_edge_walk']['witness']
    assert saved['bounded_edge_walk']['stop'] == 'spacelike_endpoint'
    assert E == [decode(v) for v in witness['endpoint']]
    assert square(E) == decode(witness['square'])
    for j in range(8):
        assert roots[f'old{j+1}'] == [Q(),Q()]+[Q(-int(i == j)) for i in range(8)]
    return phi,G,ci,roots,cov,norms,center,V,direction,E,white


def profiles(phi, ci):
    accepted, summaries = [], []
    chain = [7,2,3,4,5,6]
    for beta, choice in product((Q(),Q(1),phi), [(None,Q())]+[(j,v) for j in chain for v in (Q(1),phi)]):
        node,value = choice; eta = [Q()]*8
        if node is not None: eta[node] = value
        energy = beta*beta/2+(Q(chain.index(node)+1)*(7-chain.index(node)-1)*value*value/7 if node is not None else 0)
        if (energy-2).sign() >= 0: continue
        R0,R1,Q0 = dot(ci[0],eta),dot(ci[1],eta),dot(eta,mv(ci,eta))
        L = 4*R1-7*R0
        H = 56-14*beta*beta-28*Q0+7*R0*R0+L*L
        assert H == 28*(2-energy) and H.sign() > 0
        excluded = node is None and beta in (Q(),Q(1))
        if excluded:
            rhs = 8-2*beta*beta
            assert ((rhs.a*rhs.a-5*rhs.b*rhs.b)/(-1)) % 7 == 6
        ident = f'next_edge_profile_{len(summaries)+1:02d}'
        summaries.append(dict(id=ident,beta=beta.encoded(),A6_node=node,
            A6_magnitude=value.encoded(),schur_energy=energy.encoded(),H=H.encoded(),
            R0=R0.encoded(),R1=R1.encoded(),Q0=Q0.encoded(),L=L.encoded(),excluded_mod7=excluded))
        if not excluded: accepted.append((ident,beta,eta,R0,L,H))
    assert len(summaries) == 13 and len(accepted) == 11
    assert sum(all(v == 0 for v in eta) for _,_,eta,_,_,_ in accepted) == 1
    # Independent reviewed long branch: phi=15, 5u=6a, 5a^2=1 mod19.
    assert (2*15**3)%19 == 5 and pow(15, -2, 19) == 6
    assert all((5*a*a)%19 != 1 for a in (0,1,15))
    return accepted,summaries


def run():
    phi,G,ci,roots,cov,norms,oldcenter,V,direction,E,white = fixed_inputs()
    selected,summaries = profiles(phi,ci)
    f,xs,solutions,survivors,counts,retained = 2*phi-3,[],[],[],Counter(),0
    for b in range(301):
        m = (isqrt(5*b*b)-b)//2
        for a in range(m-2,m+4):
            counts['x_coefficient_pairs'] += 1
            x = native(a,b)
            if x.sign() > 0 and (x-600).sign() <= 0 and (x.conjugate()+2).sign() >= 0 and (x.conjugate()-2).sign() <= 0:
                xs.append((a,b,x))
    assert counts['x_coefficient_pairs'] == 1806
    counts['x_in_requested_domain'] = len(xs)
    square = lambda v: dot(v,mv(G,v))
    for ident,beta,eta,R0,L,H in selected:
        for a,b,x in xs:
            counts['profile_x_square_tests'] += 1
            positive = field_sqrt(H+7*f*x*x)
            if positive is None: continue
            counts['field_square_solutions'] += 1
            for B in (positive,-positive):
                counts['signed_B_branches'] += 1
                gamma = (B-L)/7
                alpha,y = (phi*phi*phi*x-R0-7*gamma)/4,(beta-phi*x)/2
                if not all(v.integral() for v in (gamma,alpha,y)):
                    counts['reject_reconstruction_nonintegral'] += 1; continue
                spatial = mv(ci,[alpha,gamma]+eta[2:]); root = [x,y]+spatial
                cv = mv(G,root)
                assert all(v.integral() for v in root+cv) and dot(root,cv) == 2
                assert B*B-7*f*x*x == H and dot(root,cov['z']) == 0
                assert dot(root,cov['old0']) == -beta
                assert all(dot(root,cov[f'old{j+1}']) == -eta[j] for j in range(2,8))
                counts['integral_primitive_norm_two_roots'] += 1
                products,invalid,categories = {},[],Counter()
                for n,r in roots.items():
                    v = dot(r,cv); kind = pair_kind(v,norms[n],phi)
                    products[n] = [v.encoded(),kind]; categories[kind] += 1
                    if kind.startswith('reject'): invalid.append(n)
                vp,dp,cut,cp = dot(V,cv),dot(direction,cv),dot(E,cv),dot(oldcenter,cv)
                assert vp == -f*x/4 and vp.sign() < 0 and dp == B/4 and cut == (B-5*f*x)/20
                epsilon = F(1)
                while (vp+epsilon*cp).sign() >= 0: epsilon /= 2
                center = add(V,oldcenter,epsilon)
                assert square(center).sign() < 0 and dot(V,mv(G,center)).sign() < 0
                assert dot(center,cv).sign() < 0 and all(dot(center,c).sign() < 0 for c in cov.values())
                reasons = ([] if cut.sign() > 0 else ['does_not_cut_E'])+([] if not invalid else ['incompatible_native_pairs'])
                record = dict(profile=ident,x_phi_coefficients=[a,b],B=B.encoded(),gamma=gamma.encoded(),
                    alpha=alpha.encoded(),y=y.encoded(),V_product=vp.encoded(),E_product=cut.encoded(),
                    old_center_product=cp.encoded(),old_center_strict=cp.sign()<0,center_epsilon=str(epsilon),
                    invalid_walls=invalid,pair_categories=dict(categories),rejected=reasons)
                solutions.append(record)
                retained += len(json.dumps(record).encode()); assert retained < 900000
                if reasons: continue
                tau = f*x/B; endpoint = add(V,direction,tau); ep = {n: dot(endpoint,c) for n,c in cov.items()}
                assert tau.sign() > 0 and (Q(F(1,5))-tau).sign() > 0
                assert square(endpoint) == -f*H/(4*B*B) and square(endpoint).sign() < 0
                assert dot(V,mv(G,endpoint)).sign() < 0 and dot(endpoint,cv) == 0
                assert {n for n,v in ep.items() if v == 0} == set(white)
                assert all(v.sign() < 0 for n,v in ep.items() if n not in white)
                assert inertia([[dot(roots[a0],cov[b0]) for b0 in white] for a0 in white]) == [8,0,0]
                active_roots = [roots[n] for n in white]+[root]
                assert inertia([[dot(u,mv(G,v)) for v in active_roots] for u in active_roots]) == [9,0,0]
                projected = add(center,root,-dot(center,cv)/2)
                assert square(projected).sign() < 0 and dot(projected,cv) == 0
                assert all(dot(projected,c).sign() < 0 for c in cov.values())
                caps = {n:dict(product=products[n][0],kind=products[n][1],
                    whole_facet_preservation=products[n][1]=='ultraparallel') for n in ('r','w')}
                full = record|dict(coordinates=enc(root),covector=enc(cv),wall_products=products,
                    strict_center=enc(center),new_facet_witness=enc(projected),endpoint_parameter=tau.encoded(),
                    endpoint=enc(endpoint),endpoint_square=square(endpoint).encoded(),
                    endpoint_old_products={n:v.encoded() for n,v in ep.items()},old_cap_facets=caps)
                survivors.append(full); retained += len(json.dumps(full).encode()); assert retained < 900000
    assert counts['profile_x_square_tests'] == 11*len(xs)
    return dict(status='PASS_BOUNDED_NEXT_EDGE_11_SHORT_PROFILE_SEARCH',source_sha256=digest(Path(__file__)),
        sources=PINS,domain=dict(x_positive_at_most=600,absolute_conjugate_x_at_most=2,
        b_inclusive=[0,300],a_rule='m-2..m+3; m=(isqrt(5*b*b)-b)//2',
        equation='B^2-7F*x^2=H; both B signs; gamma,alpha,y in O'),
        old_walls=28,fixed_white_walls=white,profiles=summaries,retained_short_profiles=11,
        long_profiles_excluded_mod19=3,stage_counts=dict(counts),integral_norm_two_solutions=solutions,
        survivors=survivors,global_finite_volume_proved=False,
        scope='Bounded x box and all11 necessary short finite profiles at the fixed V,E edge only. '
        'All28 signed native pairs,integrality,actual finite endpoint and common center are checked. '
        'Old-center positivity never rejects. Primitive norm2 follows from evenness of G. '
        'Full H8 facet counts remain inherited; preservation is recorded only for ultraparallel caps. '
        'No other edge,all-height short census,completion,nonreflectivity or dimension claim.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true'); group.add_argument('--verify',action='store_true')
    args = parser.parse_args(); signal.alarm(45)
    result = run(); payload = json.dumps(result,sort_keys=True,indent=2)+'\n'
    assert len(payload.encode()) < 1024**2
    if args.verify: assert json.loads(CERT.read_text()) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream: stream.write(payload)
    signal.alarm(0)
    print(json.dumps(dict(status=result['status'],counts=result['stage_counts'],survivors=len(result['survivors']))))


if __name__ == '__main__':
    main()
