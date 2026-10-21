#!/usr/bin/env python3
"""Independent rank-based clipping review with certified interval ray DAGs.

No producer arithmetic is imported. Starting rays are reconstructed in a
separate quartic primitive field. New rays are exact positive convex
combinations, with independent512-bit intervals proving all strict signs.
Adjacency uses common-active rank8, not the producer's face-ray count test.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sympy as sp
import fresh_h9_joint_corner_20261015_bridge_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_escape35_20261016_cone_review'
PRIMARY = 'results/fresh_nonarithmetic_escape35_20261016_cone.json'
WALLS = 'results/fresh_nonarithmetic_escape35_20261016.json'
LOCAL_REVIEW = 'results/fresh_nonarithmetic_escape35_20261016_review.json'
START = 'results/fresh_h9_joint_pair_20261015_two_extension.json'
START_REVIEW = 'results/fresh_h9_joint_pair_20261015_two_review.json'
BITS = 512
DEN = 1 << BITS


class I:
    __slots__ = ('lo', 'hi')
    def __init__(self, lo=0, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi
    @staticmethod
    def get(v): return v if isinstance(v, I) else I(v)
    @staticmethod
    def out(lo, hi):
        return I(F((lo * DEN).__floor__(), DEN), F((hi * DEN).__ceil__(), DEN))
    def __add__(self, v):
        v = I.get(v)
        return I.out(self.lo + v.lo, self.hi + v.hi)
    __radd__ = __add__
    def __neg__(self): return I(-self.hi, -self.lo)
    def __sub__(self, v): return self + -I.get(v)
    def __rsub__(self, v): return I.get(v) + -self
    def __mul__(self, v):
        v = I.get(v)
        products = [x * y for x in (self.lo, self.hi) for y in (v.lo, v.hi)]
        return I.out(min(products), max(products))
    __rmul__ = __mul__
    def __truediv__(self, v):
        v = I.get(v)
        assert v.lo > 0 or v.hi < 0
        return self * I.out(1 / v.hi, 1 / v.lo)
    def sqrt(self):
        assert self.lo > 0
        lo = isqrt((self.lo * DEN * DEN).__floor__())
        hi = isqrt((self.hi * DEN * DEN).__floor__()) + 1
        assert F(lo * lo, DEN * DEN) <= self.lo and F(hi * hi, DEN * DEN) > self.hi
        return I(F(lo, DEN), F(hi, DEN))
    def strict_sign(self):
        assert self.lo > 0 or self.hi < 0, ('unresolved strict sign', self.lo, self.hi)
        return 1 if self.lo > 0 else -1
    def saved(self):
        den = 1 << 192
        lo, hi = F((self.lo * den).__floor__(), den), F((self.hi * den).__ceil__(), den)
        if self.lo > 0 and lo <= 0 or self.hi < 0 and hi >= 0:
            lo, hi = self.lo, self.hi
        return [str(lo), str(hi)]


def interval_expression(x):
    if x.is_Rational: return I(F(int(x.p), int(x.q)))
    if x.is_Add: return sum((interval_expression(v) for v in x.args), I())
    if x.is_Mul:
        out = I(1)
        for v in x.args: out = out * interval_expression(v)
        return out
    if x.is_Pow:
        a, p = x.args
        if p == sp.Rational(1, 2): return interval_expression(a).sqrt()
        if p == -1: return I(1) / interval_expression(a)
        if p.is_Integer and p > 0:
            out = I(1)
            for _ in range(int(p)): out = out * interval_expression(a)
            return out
    raise AssertionError(x)


def idot(a, b): return sum((x * y for x, y in zip(a, b)), I())
def imv(a, v): return [idot(row, v) for row in a]
def ids(mask): return [i for i in range(mask.bit_length()) if mask >> i & 1]


def rank_minor(matrix):
    """Certified lower rank, via complete pivoting with interval exclusions."""
    rows = [list(row) for row in matrix]
    row_ids, col_ids, pivots = list(range(len(rows))), list(range(10)), []
    for k in range(min(8, len(rows))):
        options = [(min(abs(rows[i][j].lo), abs(rows[i][j].hi)), i, j)
            for i in range(k, len(rows)) for j in range(k, 10)
            if rows[i][j].lo > 0 or rows[i][j].hi < 0]
        if not options: break
        _, i, j = max(options)
        rows[k], rows[i] = rows[i], rows[k]
        row_ids[k], row_ids[i] = row_ids[i], row_ids[k]
        for row in rows: row[k], row[j] = row[j], row[k]
        col_ids[k], col_ids[j] = col_ids[j], col_ids[k]
        pivot = rows[k][k]
        pivots.append(pivot)
        for i in range(k + 1, len(rows)):
            multiplier = rows[i][k] / pivot
            for j in range(k + 1, 10): rows[i][j] = rows[i][j] - multiplier * rows[k][j]
            rows[i][k] = I()
    return dict(rank_lower_bound=len(pivots), original_row_indices=row_ids[:len(pivots)],
        original_columns=col_ids[:len(pivots)], nonzero_pivot_bounds=[p.saved() for p in pivots])


def integrity(paths):
    queue, seen, count = list(paths), set(), 0
    while queue:
        p = queue.pop()
        if not p.endswith('.json') or p in seen: continue
        seen.add(p)
        for source, digest in json.loads((ROOT / p).read_text()).get('sources', {}).items():
            assert sha256((ROOT / source).read_bytes()).hexdigest() == digest, (p, source)
            count += 1
            if source.endswith('.json'): queue.append(source)
    return dict(manifests_checked=len(seen), source_references_checked=count)


def run():
    target, wall_doc, local, start, inherited = [json.loads((ROOT / p).read_text()) for p in (PRIMARY, WALLS, LOCAL_REVIEW, START, START_REVIEW)]
    assert inherited['status'] == 'PASS_INDEPENDENT_TWO_CAPS_AND_COMPLETE_32_WALL_CONE'
    assert inherited['sources'][START] == sha256((ROOT / START).read_bytes()).hexdigest()
    assert inherited['complete_census'] == start['complete_census']
    assert local['status'] == 'PASS_INDEPENDENT_LAMBERT_CAP_FINITE_CENSUS_AND_RENEWED_ESCAPE'
    assert local['sources'][WALLS] == sha256((ROOT / WALLS).read_bytes()).hexdigest()
    # A separate quartic field avoids the producer's octic reconstruction and
    # avoids importing any old module that mutates a shared arithmetic field.
    sq5 = sp.sqrt(5)
    rad = 6 + sp.Rational(16, 5) * sq5
    K4 = sp.QQ.algebraic_field(sq5, sp.sqrt(rad))
    zero, one = K4.zero, K4.one
    basis = [K4.from_sympy(v) for v in (sp.S.One, sq5, sp.sqrt(rad), sq5 * sp.sqrt(rad))]
    def decode4(raw): return sum((K4.convert(F(x)) * v for x, v in zip(raw[0] + raw[1], basis)), zero)
    def bdot(a, b): return sum((x * y for x, y in zip(a, b)), zero)
    def bmv(a, v): return [bdot(row, v) for row in a]
    def rref(matrix):
        rows, pivots = [list(row) for row in matrix], []
        # Reverse column order differs from both prior cone eliminators.
        for j in reversed(range(10)):
            k = len(pivots)
            at = next((i for i in range(k, len(rows)) if rows[i][j]), None)
            if at is None: continue
            rows[k], rows[at] = rows[at], rows[k]
            pivot = rows[k][j]
            rows[k] = [x / pivot for x in rows[k]]
            for i in range(len(rows)):
                if i != k and rows[i][j]:
                    scale = rows[i][j]
                    rows[i] = [x - scale * y for x, y in zip(rows[i], rows[k])]
            pivots.append(j)
            if len(pivots) == len(rows): break
        return rows, pivots
    theta4, theta8 = interval_expression(K4.ext.as_expr()), interval_expression(R.K.ext.as_expr())
    @lru_cache(None)
    def bound4(x):
        v = I()
        for c in x.to_list(): v = v * theta4 + I(F(str(c)))
        return v
    @lru_cache(None)
    def bound8(x):
        v = I()
        for c in x.to_list(): v = v * theta8 + I(F(str(c)))
        return v
    def decode8(raw): return R.K.new([R.K.dom.convert(F(x)) for x in raw])
    Q4 = [[decode4(x) for x in row] for row in start['form']]
    oldN = [[decode4(x) for x in row] for row in start['normals']]
    H4 = [decode4(x) for x in start['center']]
    oldcov = [bmv(Q4, n) for n in oldN]
    assert len(oldN) == 32 and len(rref(oldcov)[1]) == 10
    q8 = [[decode8(x) for x in row] for row in wall_doc['form']]
    assert q8 == R.m4(start['form'])
    assert [decode8(x) for x in wall_doc['center']] == R.v4(start['center'])
    QI = [[bound4(x) for x in row] for row in Q4]
    HI = [bound4(x) for x in H4]
    hcov = imv(QI, HI)
    assert idot(HI, hcov).strict_sign() < 0
    D = [R.O + R.BASIS[1], 24 * R.O + 72 * R.BASIS[1] / 5, 10 * R.O + 68 * R.BASIS[1] / 15]
    for p, delta in zip(('results/fresh_nonarithmetic_wall3_cap_20261016.json', 'results/fresh_nonarithmetic_escape34_20261016.json', WALLS), D):
        doc = json.loads((ROOT / p).read_text())
        assert decode8(doc['field']['extra_radical_square']) == delta
    radicals = [bound8(x).sqrt() for x in D]
    monomials = []
    for mask in range(8):
        v = I(1)
        for j in range(3):
            if mask >> j & 1: v = v * radicals[j]
        monomials.append(v)
    scalar_rows, N = [], []
    for row in wall_doc['normals']:
        exact = [[decode8(x[b2][b1][b0]) for b2 in range(2) for b1 in range(2) for b0 in range(2)] for x in row]
        scalar_rows.append(exact)
        N.append([sum((bound8(c) * m for c, m in zip(coefficients, monomials)), I()) for coefficients in exact])
    assert len(N) == 36
    assert [[x[0] for x in row] for row in scalar_rows[:32]] == R.m4(start['normals'])
    assert all(not any(x[1:]) for row in scalar_rows[:32] for x in row)
    cov = [imv(QI, n) for n in N]
    for i in range(32):
        # Use direct quartic bounds for old rows, avoiding needless octic
        # coordinate cancellation in rank certificates and ray signs.
        cov[i] = [bound4(x) for x in oldcov[i]]
    center_pairs = [idot(HI, c) for c in cov]
    assert all(p.strict_sign() < 0 for p in center_pairs)
    fixed_functional = [-sum((row[j] for row in oldcov), zero) for j in range(10)]
    rays, initial_counts = [], Counter()
    starting_digest = sha256()
    for kind, supports in start['complete_census'].items():
        for support in supports:
            rr, piv = rref([oldcov[i] for i in support])
            assert len(piv) == 9
            free = next(j for j in range(10) if j not in piv)
            v = [zero] * 10
            v[free] = one
            for k, j in enumerate(piv): v[j] = -rr[k][free]
            assert all(not bdot(v, oldcov[i]) for i in support)
            scale = bdot(v, fixed_functional)
            assert scale
            v = [x / scale for x in v]
            values = [bdot(v, row) for row in oldcov]
            assert [i for i, x in enumerate(values) if not x] == support
            assert all(not x or bound4(x).strict_sign() < 0 for x in values)
            assert bdot(v, fixed_functional) == one
            vi = [bound4(x) for x in v]
            norm, time = idot(vi, imv(QI, vi)), idot(vi, hcov)
            actual = 'hyperideal_rays' if norm.strict_sign() > 0 else 'finite_vertices' if time.strict_sign() < 0 else 'opposite_sheet_rays'
            assert actual == kind
            mask = sum(1 << i for i in support)
            rays.append(dict(v=vi, mask=mask, kind=actual, norm=norm, time=time))
            initial_counts[kind] += 1
            starting_digest.update(json.dumps([support, [list(map(str, x.to_list())) for x in v]], separators=(',', ':')).encode())
            if len(rays) % 200 == 0: print('independent quartic initial rays', len(rays), flush=True)
    assert len(rays) == len({r['mask'] for r in rays}) == 972
    stages = []
    for wall in range(32, 36):
        values = [idot(r['v'], cov[wall]) for r in rays]
        signs = [v.strict_sign() for v in values]
        negative = [i for i, s in enumerate(signs) if s < 0]
        positive = [i for i, s in enumerate(signs) if s > 0]
        new = {rays[i]['mask']: rays[i] for i in negative}
        cardinalities, rejected, crossings = Counter(), [], []
        for i in negative:
            for j in positive:
                common = rays[i]['mask'] & rays[j]['mask']
                support = ids(common)
                cardinalities[len(support)] += 1
                if len(support) < 8: continue
                minor = rank_minor([cov[k] for k in support])
                if minor['rank_lower_bound'] < 8:
                    assert wall == 32 and max(support) < 32
                    _, piv = rref([oldcov[k] for k in support])
                    assert len(piv) < 8
                    rejected.append(dict(negative_support=ids(rays[i]['mask']), positive_support=ids(rays[j]['mask']),
                        common_support=support, exact_quartic_rank=len(piv)))
                    continue
                # Common rows kill two distinct normalized rays; rank cannot
                # exceed8. The certified rank8 minor therefore proves a2-face.
                a, b = values[i], values[j]
                denominator = b - a
                alpha, beta = b / denominator, -a / denominator
                assert alpha.lo > 0 and beta.lo > 0
                vector = [alpha * x + beta * y for x, y in zip(rays[i]['v'], rays[j]['v'])]
                value = idot(vector, cov[wall])
                assert value.lo <= 0 <= value.hi
                norm, time = idot(vector, imv(QI, vector)), idot(vector, hcov)
                kind = 'hyperideal_rays' if norm.strict_sign() > 0 else 'finite_vertices' if time.strict_sign() < 0 else 'opposite_sheet_rays'
                mask = common | (1 << wall)
                # This is a uniqueness assertion, never silent deduplication.
                assert mask not in new
                new[mask] = dict(v=vector, mask=mask, kind=kind, norm=norm, time=time)
                crossings.append(dict(negative_support=ids(rays[i]['mask']), positive_support=ids(rays[j]['mask']),
                    new_support=ids(mask), common_rank8_minor=minor, positive_convex_coefficients=[alpha.saved(), beta.saved()],
                    kind=kind, norm_bounds=norm.saved(), center_pair_bounds=time.saved()))
        rays = sorted(new.values(), key=lambda r: ids(r['mask']))
        counts = Counter(r['kind'] for r in rays)
        stage = dict(wall=wall, negative=len(negative), positive=len(positive), zero=0,
            all_crossing_pair_common_cardinalities={str(k): v for k, v in sorted(cardinalities.items())},
            exact_rank_nonadjacency=rejected, crossings=crossings, counts=dict(counts), total_rays=len(rays))
        producer = target['clipping_stages'][wall - 32]
        assert {k: stage[k] for k in ('wall', 'negative', 'positive', 'zero', 'counts', 'total_rays')} == {k: producer[k] for k in ('wall', 'negative', 'positive', 'zero', 'counts', 'total_rays')}
        assert len(crossings) == producer['crossing_edges']
        assert sum(v for k, v in cardinalities.items() if k >= 8) == producer['adjacency_candidates']
        assert [{k: r[k] for k in ('negative_support', 'positive_support', 'new_support', 'kind')} for r in crossings] == producer['crossings']
        stages.append(stage)
        print('independent rank clip', wall, dict(counts), 'crossings', len(crossings), 'rank rejections', len(rejected), flush=True)
    final = []
    for index, ray in enumerate(rays):
        values = [idot(ray['v'], c) for c in cov]
        active = ids(ray['mask'])
        for i, value in enumerate(values):
            if i in active: assert value.lo <= 0 <= value.hi
            else: assert value.strict_sign() < 0
        final.append(dict(active=active, kind=ray['kind'], norm_bounds=ray['norm'].saved(),
            center_pair_bounds=ray['time'].saved(), normalized_coordinate_bounds=[x.saved() for x in ray['v']]))
        if (index + 1) % 300 == 0: print('final full36 inequality rows', index + 1, flush=True)
    census = {kind: sorted(row['active'] for row in final if row['kind'] == kind) for kind in ('finite_vertices', 'cusps', 'hyperideal_rays', 'opposite_sheet_rays')}
    assert census == target['complete_census']
    assert census['finite_vertices'] == sorted(wall_doc['finite_vertex_incidence'])
    counts = {k: len(v) for k, v in census.items()}
    assert counts == dict(finite_vertices=100, cusps=0, hyperideal_rays=1068, opposite_sheet_rays=0)
    assert sum(len(s['exact_rank_nonadjacency']) for s in stages) == 64
    paths = [f'scripts/{TAG}.py', PRIMARY, WALLS, LOCAL_REVIEW, START, START_REVIEW,
        'results/fresh_nonarithmetic_wall3_cap_20261016.json', 'results/fresh_nonarithmetic_escape34_20261016.json',
        'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
        'scripts/fresh_h9_h8_all_high_20261014_review.py',
        'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_COMPLETE_36_WALL_CONE_BY_RANK_CLIPPING',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        recursive_input_integrity=integrity(paths), interval_bits=BITS,
        arithmetic='Separate quartic primitive field for starting nullspaces and rank failures; independently rounded rational intervals for the octic tower and all subsequent strict signs. Exact positive-combination expression DAGs establish new coordinates and all zero identities. No producer arithmetic imports.',
        inherited_complete_cone=START_REVIEW, initial_counts=dict(initial_counts), initial_ray_count=972,
        initial_coordinate_digest=starting_digest.hexdigest(), normalization='f(v)=1, where f(v)=-sum_{i<32}(v.n_i). This functional is strictly positive on every nonzero old-cone vector because all old inequalities are nonpositive and their covectors have rank10.',
        full_dimension_and_pointedness=dict(old_covector_rank=10, all36_strict_center_pair_bounds=[x.saved() for x in center_pairs],
            argument='Rank10 implies zero lineality. The independently checked strict center gives full dimension10 after every clip.'),
        adjacency_proof='For two distinct normalized extreme rays u,v, their sum is relatively interior in the common-zero face: every other inequality is strict. Thus its dimension is10-rank(common active covectors). Adjacency is exactly rank8. Fewer than8 common walls cannot have rank8. Each accepted pair has a certified nonzero8x8 minor; two independent common-null rays give the matching upper rank8. Every rejected remaining candidate has exact quartic rank<8.',
        positive_combination_proof='For current values a<0<b, w=(b*u-a*v)/(b-a) has strictly positive coefficients summing to1. Its old zero set is exactly the intersection of the endpoint zero sets; its new-wall contact vanishes identically. The new equation is independent of the rank8 old support, so w has active rank9 and is an extreme ray.',
        mask_uniqueness_proof='Every active mask has rank9 in dimension10, hence a one-dimensional kernel; pointedness selects one ray orientation and f=1 selects one vector. The code asserts distinct masks and never silently discards duplicates.',
        producer_face_criterion_audit='For a pointed polyhedral cone with a complete ray list, a common-zero face contains exactly two distinct rays iff it is2-dimensional. Thus the producer criterion is equivalent to the independently verified rank8 condition. Completeness is inherited once at32 walls and propagated by the standard halfspace-clipping theorem at all four steps.',
        clipping_stages=stages, counts=counts, complete_census=census,
        final_rays=final, final_full_inequality_checks=36 * len(final),
        scope='Complete polyhedral cone census for the specified36 walls, conditional only on the frozen independently complete32-wall starting cone. Exactly100 finite rays,0 ideal rays and1068 hyperideal rays; no opposite-sheet rays. This does not establish finite volume; the separately reviewed strict escape proves infinite volume.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify: assert out == json.loads(path.read_text())
    else: path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], out['counts'], flush=True)
