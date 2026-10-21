#!/usr/bin/env python3
"""Independent geometry and simultaneous polar batch review.

Uses frozen independent nested-radical and interval engines, never a new
producer import. Helpers also support the independent complete-cone review.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_nonarithmetic_escape35_20261016_review as A
import fresh_nonarithmetic_escape35_20261016_cone_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_joint_regions_20261016_batch_review'
BATCH = 'results/fresh_nonarithmetic_joint_regions_20261016_batch.json'
POLAR = 'results/fresh_nonarithmetic_joint_regions_20261016_polar_gate.json'
OLD = 'results/fresh_nonarithmetic_escape35_20261016.json'
OLD_REVIEW = 'results/fresh_nonarithmetic_escape35_20261016_review.json'
CONE_REVIEW = 'results/fresh_nonarithmetic_escape35_20261016_cone_review.json'
T, I, K, Z, O = A.T, R.I, A.K, A.Z, A.O
THETA = R.interval_expression(K.ext.as_expr())


def digest(p): return sha256((ROOT / p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT / p).read_text())
def db(raw): return K.new([K.dom.convert(F(v)) for v in raw])


def dec(raw):
    assert len(raw) == 8
    return A.dt([[[raw[4*b2+2*b1+b0] for b0 in range(2)]
                  for b1 in range(2)] for b2 in range(2)])


@lru_cache(None)
def ib(x):
    out = I()
    for c in x.to_list(): out = out * THETA + I(F(str(c)))
    return out


SQ = [ib(A.P.P.D).sqrt(), ib(A.P.D).sqrt(), ib(A.D).sqrt()]


def iv(x):
    x = T(x)
    def low(v): return ib(v.a) + SQ[0] * ib(v.b)
    def mid(v): return low(v.a) + SQ[1] * low(v.b)
    return mid(x.a) + SQ[2] * mid(x.b)


def sg(x):
    if x == 0: return 0
    v = iv(x)
    if v.lo > 0: return 1
    if v.hi < 0: return -1
    return x.sign()


def dot(a, b): return A.dot(a, b)
def mv(q, v): return A.mv(q, v)
def pair(a, b, q): return dot(a, mv(q, b))
def idot(a, b): return R.idot(a, b)
def imv(a, b): return R.imv(a, b)
def ids(m): return R.ids(m)


def geometry():
    batch, old, oldrev, conerev = map(read, [BATCH, OLD, OLD_REVIEW, CONE_REVIEW])
    assert oldrev['sources'][OLD] == digest(OLD)
    assert conerev['status'] == 'PASS_INDEPENDENT_COMPLETE_36_WALL_CONE_BY_RANK_CLIPPING'
    assert conerev['sources'][OLD] == digest(OLD)
    Q = [[db(v) for v in row] for row in batch['form']]
    assert Q == [[db(v) for v in row] for row in old['form']]
    H = list(map(dec, batch['center']))
    assert H == [T(db(v)) for v in old['center']]
    N = [list(map(dec, row)) for row in batch['batch_normals']]
    assert N[:36] == [A.vt(row) for row in old['normals']]
    assert len(N) == 48
    cov = [mv(Q, n) for n in N]
    QI = [[ib(v) for v in row] for row in Q]
    HI, NI, covI = [iv(v) for v in H], [[iv(x) for x in n] for n in N], [[iv(x) for x in c] for c in cov]
    hcov = mv(Q, H)
    old_rays = {tuple(r['active']): r for r in conerev['final_rays']}
    old_cross = {s['wall']: {tuple(c['new_support']): c for c in s['crossings']}
                 for s in conerev['clipping_stages']}
    basecov = [[A.base(v) for v in row] for row in cov[:32]]

    @lru_cache(None)
    def exact_old(support):
        """An independently reconstructed positive representative; no fixed scale."""
        assert support in old_rays or max(support) <= 35
        last = max(support)
        if last < 32:
            rows = [basecov[i][:] for i in support]
            pivots = []
            for col in reversed(range(10)):
                k = len(pivots)
                j = next((j for j in range(k, len(rows)) if rows[j][col]), None)
                if j is None: continue
                rows[k], rows[j] = rows[j], rows[k]
                lead = rows[k][col]
                rows[k] = [x / lead for x in rows[k]]
                for j in range(len(rows)):
                    if j != k and rows[j][col]:
                        lead = rows[j][col]
                        rows[j] = [x - lead*y for x, y in zip(rows[j], rows[k])]
                pivots.append(col)
            assert len(pivots) == 9
            free = next(i for i in range(10) if i not in pivots)
            v = [Z] * 10
            v[free] = O
            for k, col in enumerate(pivots): v[col] = -rows[k][free]
            f = -sum((sum((x*y for x, y in zip(v, c)), Z) for c in basecov), Z)
            assert f
            return tuple(T(x/f) for x in v)
        row = old_cross[last][support]
        a = exact_old(tuple(row['negative_support']))
        b = exact_old(tuple(row['positive_support']))
        ga, gb = dot(a, cov[last]), dot(b, cov[last])
        assert sg(ga) < 0 < sg(gb)
        return tuple(gb*x-ga*y for x, y in zip(a, b))

    return dict(batch=batch, old=old, inherited=conerev, Q=Q, QI=QI, H=H, HI=HI,
                hcov=hcov, hI=[iv(v) for v in hcov], N=N, NI=NI, cov=cov,
                covI=covI, old_rays=old_rays, exact_old=exact_old)


def relation(g, product_norm):
    gs = sg(g)
    if gs == 0: return {'kind': 'right', 'label': 2}
    if gs > 0: return {'kind': 'obtuse'}
    square = g*g
    if sg(square-product_norm) > 0: return {'kind': 'dotted'}
    sq5 = T(A.R.BASIS[1])
    for label, value in [(3, T(1)), (4, T(2)), (5, (T(3)+sq5)/2),
                         (6, T(3)), (10, (T(5)+sq5)/2)]:
        if sg(4*square-value*product_norm) == 0:
            return {'kind': 'finite', 'label': label}
    assert sg(4*square-product_norm) > 0 and sg(2*product_norm-4*square) > 0
    return {'kind': 'forbidden_finite_gap', 'neighbor_labels': [3, 4]}


def run():
    G = geometry()
    batch, Q, H, N, cov = [G[k] for k in ('batch', 'Q', 'H', 'N', 'cov')]
    hh = pair(H, H, Q)
    assert sg(hh) < 0
    caps, norms, records = [], [], []
    for index, row in enumerate(batch['candidates']):
        v = list(map(dec, row['normal']))
        norm, hp = pair(v, v, Q), dot(v, G['hcov'])
        assert sg(norm) > 0 and sg(hp) < 0
        values = [dot(v, c) for c in cov[:36]]
        active = [i for i, x in enumerate(values) if sg(x) == 0]
        assert active == row['right_support']
        inherited = G['old_rays'][tuple(active)]
        assert inherited['kind'] == 'hyperideal_rays'
        # The inherited active rank9 gives a unique projective ray. The
        # negative center product fixes its orientation, independently of
        # the producer's crossing reconstruction or normalization.
        relations = [relation(x, norm*pair(n, n, Q)) for x, n in zip(values, N[:36])]
        assert relations == row['old_relations']
        assert all(r['kind'] in ('right', 'finite', 'dotted') for r in relations)
        assert norm == dec(row['norm']) and hp == dec(row['center_pairing'])
        assert values == list(map(dec, row['old_pairings']))
        if index < 12: assert v == N[36+index]
        caps.append(v); norms.append(norm)
        records.append(dict(index=index, ray_index=row['ray_index'], active=active,
                            norm_bounds=iv(norm).saved(), center_pair_bounds=iv(hp).saved(),
                            relations=relations))
        print('independent polar contacts', index, flush=True)
    pairs, bad = [], []
    for i, j in combinations(range(15), 2):
        g = pair(caps[i], caps[j], Q)
        rel = relation(g, norms[i]*norms[j])
        expected = next(p for p in batch['mutual_pairs'] if p['pair'] == [i, j])
        assert rel == {k:v for k,v in expected.items() if k not in ('pair', 'pairing')}
        assert g == dec(expected['pairing'])
        if rel['kind'] not in ('right', 'finite', 'dotted'): bad.append([i, j, rel['kind']])
        pairs.append(dict(pair=[i,j], **rel))
    assert bad == batch['incompatible_pairs']
    assert all(max(row[:2]) >= 12 for row in bad)
    hp = [dot(H, c) for c in cov]
    q = [dot(n,c) for n,c in zip(N,cov)]
    assert all(sg(v) < 0 for v in hp) and all(sg(v) > 0 for v in q)
    facets = []
    for i in range(48):
        v = [q[i]*h-hp[i]*n for h,n in zip(H,N[i])]
        assert v == list(map(dec,batch['facet_witnesses'][i]['point']))
        norm = q[i]*q[i]*hh-q[i]*hp[i]*hp[i]
        assert sg(norm) < 0 and sg(q[i]*hh) < 0 and dot(v,cov[i]) == 0
        # Every off-wall value q_i*h_j-h_i*g_ij is strictly negative.
        for j in range(48):
            if i != j: assert sg(dot(N[i],cov[j])) <= 0
        facets.append(dict(wall=i,norm_bounds=iv(norm).saved()))
    print('independent simultaneous48 actual facets PASS', flush=True)
    escape = list(map(dec,batch['common_strict_escape']['point']))
    ray = list(map(dec,batch['common_strict_escape']['ray']))
    assert dot(ray,G['hcov']) == -1
    assert escape == [v+T(K.convert(F(1,8192)))*h for v,h in zip(ray,H)]
    assert sg(pair(escape,escape,Q)) > 0 and sg(dot(escape,G['hcov'])) < 0
    for n in N[:36]+caps: assert sg(pair(escape,n,Q)) < 0
    preserve = 0
    for old in G['old_rays'].values():
        if old['kind'] != 'finite_vertices': continue
        point = [I(*map(F,interval)) for interval in old['normalized_coordinate_bounds']]
        assert all(idot(point,c).strict_sign() < 0 for c in G['covI'][36:])
        preserve += 1
    assert preserve == 100
    extra = read(POLAR)['additional_exact_polars']
    assert len(extra) == 1 and extra[0]['ray_index'] == 8
    extra_n = list(map(dec, extra[0]['normal']))
    extra_q = pair(extra_n,extra_n,Q)
    extra_rel = [relation(pair(extra_n,n,Q),extra_q*pair(n,n,Q)) for n in N[:36]]
    assert extra_rel == extra[0]['relations']
    assert all(r['kind'] in ('right','finite','dotted') for r in extra_rel)
    assert [i for i,c in enumerate(cov[:36]) if sg(dot(extra_n,c)) == 0] == extra[0]['active']
    assert tuple(extra[0]['active']) in G['old_rays']
    assert sg(pair(extra_n,escape,Q)) < 0 and sg(pair(extra_n,H,Q)) < 0
    assert all(relation(pair(extra_n,n,Q),extra_q*q)['kind'] == 'dotted' for n,q in zip(caps,norms))
    paths = [BATCH,POLAR,OLD,OLD_REVIEW,CONE_REVIEW,f'scripts/{TAG}.py',
             'scripts/fresh_nonarithmetic_escape35_20261016_review.py',
             'scripts/fresh_nonarithmetic_escape35_20261016_cone_review.py']
    return dict(status='PASS_INDEPENDENT_SIXTEEN_POLARS_AND_SIMULTANEOUS_48_FACETS',
                sources={p:digest(p) for p in paths},recursive_integrity=R.integrity(paths),
                polar_candidates=records, mutual_pairs=pairs,incompatible_pairs=bad,
                facet_certificates=facets,strict_old_finite_vertices_preserved=preserve,
                common_escape_norm_bounds=iv(pair(escape,escape,Q)).saved(),
                additional_polar=dict(ray_index=8,relations=extra_rel,common_escape_pair_bounds=iv(pair(extra_n,escape,Q)).saved()),
                scope='Independent nested radical exact contacts and simultaneous48facet coherence,100oldfinite vertices retained,16acceptedpolars sharing a strict spacelike escape. The complete initial36cone is inherited from its independent review. Completeness of the all-label polar classification and the48cone census are separate reviews; no numerical candidate bank proves either.')


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify: assert result==read(str(p.relative_to(ROOT)))
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
