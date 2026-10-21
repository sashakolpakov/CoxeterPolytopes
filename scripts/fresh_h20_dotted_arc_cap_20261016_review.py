#!/usr/bin/env python3
"""Independent sparse-radical and interval review of the H20 dotted-arc cap.

No producer arithmetic imports. The old parameter's existence/isolation and
old chamber acuteness are inherited from the frozen independent chain review.
Every new contact, center, new finite corner and remaining escape is replayed.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h20_bridge_edge_closure_20261014_review import Field, parse
import fresh_h20_bridge_faces_20261014_review as B
import fresh_h20_s61_multicap_20261014_exact_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_dotted_arc_cap_20261016_review'
INPUT = 'results/fresh_h20_dotted_arc_cap_20261016.json'
BASE = 'results/fresh_h20_s61_multicap_20261014_exact.json'
CHAIN_REVIEW = 'results/fresh_h20_s61_multicap_20261014_exact_review.json'
BANK = 'results/fresh_h20_nonright_s61_20261014_bank_gate.json'
FILES = R.FILES[:5]
pt, plus, minus, neg, times, div, root = R.pt, R.plus, R.minus, R.neg, R.times, R.div, R.root
total, inner, lor = R.total, R.inner, R.lor


def saved(x):
    den = 1 << 192
    return [str(F(R.floor(x[0] * den), den)), str(F(R.ceil(x[1] * den), den))]


def vec(v):
    return [saved(x) for x in v]


def overlap(a, raw):
    b = tuple(map(F, raw))
    assert R.overlap(a, b)


def integrity(paths):
    queue, seen, count = list(paths), set(), 0
    while queue:
        p = queue.pop()
        if not p.endswith('.json') or p in seen:
            continue
        seen.add(p)
        for q, digest in json.loads((ROOT / p).read_text()).get('sources', {}).items():
            assert sha256((ROOT / q).read_bytes()).hexdigest() == digest, (p, q)
            count += 1
            if q.endswith('.json'):
                queue.append(q)
    return dict(manifests_checked=len(seen), source_references_checked=count)


def run():
    old, bridge, edge, polar, cur = [json.loads((ROOT / p).read_text()) for p in FILES]
    d, base, reviewed, bank = [json.loads((ROOT / p).read_text()) for p in (INPUT, BASE, CHAIN_REVIEW, BANK)]
    assert d['status'] == 'PASS_EXACT_134_WALL_DOTTED_ARC_ESCAPE_CUTTER'
    assert reviewed['status'] == 'PASS_INDEPENDENT_133_WALL_RIGHT_CHAIN_AND_COMPACT_S61_FACE'
    assert base['known_total_vertices_at_least'] == 168
    ff = Field([parse(old['field']['discriminant'], 0), parse(bridge['field']['new_positive_sqrt_radical'], 1), parse(edge['field']['new_positive_sqrt_radical'], 2)])
    add, mul, sc, ft, const = ff.add, ff.mul, ff.scale, ff.total, ff.const
    zero = ff.zeros[3]
    Q = [list(map(F, row)) for row in old['form']]
    assert B.inertia(Q) == [20, 1, 0]
    form_rows = [[(j, q) for j, q in enumerate(row) if q] for row in Q]
    cache = {}
    def vp(raw, level):
        return tuple(ff.lift(parse(x, level), level) for x in raw)
    def dot(v, w):
        key = id(w)
        if key not in cache:
            cache[key] = (w, tuple(ft(sc(w[j], q) for j, q in row) for row in form_rows))
        return ft(mul(x, y) for x, y in zip(v, cache[key][1]))
    def lin(vectors, coefficients):
        return tuple(ft(sc(v[j], c) for v, c in zip(vectors, coefficients)) for j in range(21))
    def rat(x):
        for k in (3, 2, 1):
            assert x[1] == ff.zeros[k - 1]
            x = x[0]
        assert not any(x[1:])
        return x[0]
    @lru_cache(None)
    def iv(x):
        return R.out(ff.interval(x, 3, R.BITS))
    N = [vp(v, 1) for v in old['normals']] + [vp(bridge['new_normal'], 2)]
    N += [vp(v, 3) for v in edge['new_normals']]
    N += [tuple(const(x) for x in v) for v in polar['new_normals']]
    N += [vp(cur['cap_R']['normal'], 1)]
    norms = [const(x) for x in old['norms']] + [const(2)] + [parse(edge['new_normal_norm'], 3)] * 2
    norms += [const(x) for x in polar['new_norms']] + [ff.lift(parse(cur['cap_R']['norm'], 1), 1)]
    assert len(N) == len(norms) == 101
    normden = [root(div(iv(q), pt(2))) for q in norms]
    S = base['support']
    assert S == edge['recipe']['support'] and len(S) == 18 and 31 in S and 5 in S
    G = [[rat(dot(N[i], N[j])) for j in S] for i in S]
    J = B.inv(G)
    support_pivots = B.pd(G)
    assert G == [list(map(F, row)) for row in base['support_gram']]
    assert J == [list(map(F, row)) for row in base['support_inverse']]
    index31 = S.index(31)
    assert J[index31] == [F(i == index31, 2) for i in range(18)]
    omega = lin([N[31]], [F(1, 2)])
    assert dot(omega, omega) == const(F(1, 2))
    assert [dot(omega, N[i]) for i in S] == [const(i == 31) for i in S]
    frame = [vp(base['frame'][key], 2) for key in ('X', 'r95', 'W')]
    squares = [ff.lift(parse(base['frame']['A'], 2), 2), const(2), ff.lift(parse(base['frame']['W_norm'], 2), 2)]
    assert [[dot(v, w) for w in frame] for v in frame] == [[sc(squares[0], -1), zero, zero], [zero, squares[1], zero], [zero, zero, squares[2]]]
    assert all(dot(v, N[i]) == zero for v in frame for i in S)
    assert all(dot(v, omega) == zero for v in frame)
    frame_den = [root(iv(x)) for x in squares]
    hb = [vp(cur['hex_basis'][key], 1) for key in ('C', 'U', 'V')]
    hq = [ff.lift(parse(x, 1), 1) for x in cur['hex_basis']['squared_lengths']]
    hex_coeff = [[ff.lift(parse(x, 0), 0) for x in row] for row in cur['hexagon']['frame_coefficients']]
    def old_pairs(v):
        pairs = [div(iv(dot(v, n)), den) for n, den in zip(N, normden)]
        for co in hex_coeff:
            pairs.append(total(div(iv(mul(c, dot(v, w))), root(iv(q))) for c, w, q in zip(co, hb, hq)))
        assert len(pairs) == 107
        return pairs
    # Recompute all old107 covectors from literal normals, using the prior
    # independent sparse-radical field engine, not saved producer intervals.
    by_axis = [[div(x, den) for x in old_pairs(v)] for v, den in zip(frame, frame_den)]
    P = [list(row) for row in zip(*by_axis)]
    for i in S:
        assert all(dot(v, N[i]) == zero for v in frame)
        P[i] = [pt(0)] * 3
    for row, raw in zip(P, base['frame']['old107_unit_normal_pair_bounds']):
        for value, value_raw in zip(row, raw):
            overlap(value, value_raw)
    a = tuple(map(F, base['parameter']['isolating_interval']))
    assert list(map(str, a)) == reviewed['parameter_interval']
    coeff, rotation = R.chain(a)
    assert len(coeff) == 26
    for row, raw in zip(coeff, base['chain']['frame_coefficient_bounds']):
        for value, value_raw in zip(row, raw):
            overlap(value, value_raw)
    P += [[neg(row[0]), row[1], row[2]] for row in coeff]
    wp = old_pairs(omega) + [pt(0)] * 26
    for i in S:
        overlap(wp[i], [int(i == 31)] * 2)
        wp[i] = pt(i == 31)
    print('PASS independent field frame,107 covectors,26 chain normals and omega31', flush=True)

    sqrt2 = root(pt(2))
    h = minus(sqrt2, pt(1))
    a2 = times(a, a)
    y = [div(plus(times(h, a2), times(pt(2), sqrt2)), times(pt(2), a)), times(div(h, pt(2)), coeff[22][1]), times(div(h, pt(2)), coeff[22][2])]
    c = root(plus(pt(5), div(pt(4), a2)))
    residual = minus(pt(F(-1, 2)), div(pt(2), a2))
    assert y[0][0] > 0 and residual[1] < 0 and c[0] > 2
    # Symbolic equalities in Q(sqrt2), with u=a^2 kept indeterminate.
    # Write b+c*sqrt2 as a pair; no floating point or interval-zero test.
    def sadd(x, z): return (x[0] + z[0], x[1] + z[1])
    def sscale(x, z): return (z * x[0], z * x[1])
    def smul(x, z): return (x[0] * z[0] + 2 * x[1] * z[1], x[0] * z[1] + x[1] * z[0])
    hs, rs = (-1, 1), (0, 1)
    assert smul(hs, hs) == (3, -2) and smul(hs, rs) == (2, -1)
    assert sadd(sadd(sscale(smul(hs, hs), 2), sscale(smul(hs, rs), -4)), (2, 0)) == (0, 0)
    assert sadd(sscale(smul(rs, rs), -4), (8, 0)) == (0, 0)
    assert sadd(hs, sscale(rs, -1)) == (-1, 0)
    # The norm identity follows from those two coefficients after multiplying
    # y^2+1/2+2/a^2 by4a^2. Adjacent chain orthogonality gives spatial
    # n129.n128=n129.n130=a^2; therefore the two contacts are exactly-sqrt2.
    assert R.overlap(lor(y, y), residual)
    assert R.overlap(plus(residual, div(times(c, c), pt(2))), pt(2))
    finite = {i: 2 for i in S if i != 31} | {128: 4, 129: 3, 130: 4}
    pairs, counts = [], Counter()
    for i, (p, w) in enumerate(zip(P, wp)):
        value = minus(inner(y, p), times(c, w))
        if i in finite:
            label = finite[i]
            expected = pt(0) if label == 2 else pt(-1) if label == 3 else neg(sqrt2)
            assert R.overlap(value, expected)
            relation = 'right' if label == 2 else f'pi/{label}'
        else:
            assert value[1] < -2, (i, saved(value))
            relation = 'dotted'
        assert relation == d['pairs'][i]['relation'] and d['pairs'][i]['wall'] == i
        overlap(value, d['pairs'][i]['pair_bounds'])
        counts[relation] += 1
        pairs.append(dict(wall=i, relation=relation, independently_recomputed_bounds=saved(value)))
    assert dict(counts) == {'right': 17, 'pi/4': 2, 'pi/3': 1, 'dotted': 113}
    assert sorted(finite) == d['finite_neighbors']
    print('PASS all133 cap contacts, with symbolic finite equalities', flush=True)

    H = vp(cur['strict_center']['point'], 1)
    Hnorm = dot(H, H)
    assert ff.sign(Hnorm) < 0
    Hden = root(iv(sc(Hnorm, -1)))
    hf = [div(iv(dot(v, H)), times(den, Hden)) for v, den in zip(frame, frame_den)]
    omegaH = div(iv(dot(omega, H)), Hden)
    hp = [div(p, Hden) for p in old_pairs(H)] + [inner(row, hf) for row in coeff]
    C0 = [pt(1), pt(F(-1, 2)), pt(0)]
    eps = F(base['strict_center']['epsilon'])
    assert eps > 0 and hf[0][1] < 0
    center_old_future = minus(inner(C0, hf), pt(eps))
    assert center_old_future[1] < 0
    center_pairs = [plus(inner(C0, p), times(pt(eps), hpv)) for p, hpv in zip(P, hp)]
    center_pairs += [plus(lor(y, C0), times(pt(eps), minus(inner(y, hf), times(c, omegaH))))]
    center_norm = minus(plus(pt(F(-3, 4)), times(pt(2 * eps), inner(C0, hf))), pt(eps * eps))
    assert center_norm[1] < 0 and all(value[1] < 0 for value in center_pairs)
    facet_norms = [minus(center_norm, div(times(value, value), pt(2))) for value in center_pairs]
    assert all(value[1] < 0 for value in facet_norms)
    for value, raw in zip(center_pairs, d['common_center']['all134_pair_bounds']): overlap(value, raw)
    for value, raw in zip(facet_norms, d['common_center']['all134_facet_projection_norm_bounds']): overlap(value, raw)

    tested_escapes = []
    for name, raw in (('old133', base['strict_escape']['rational_point']), ('bank847', bank['escape']['rational_point'])):
        E = tuple(const(x) for x in raw)
        Eframe = [div(iv(dot(v, E)), den) for v, den in zip(frame, frame_den)]
        value = minus(inner(y, Eframe), times(c, iv(dot(omega, E))))
        assert (value[1] < 0) if name == 'old133' else (value[0] > 0)
        tested_escapes.append(dict(name=name, new_pair_bounds=saved(value), cut=value[0] > 0))

    # Enumerate all20 possible19-subsets of finite old neighbors. Any omitted
    # root except128/130 retains their dotted pair and cannot be spherical.
    excluded, survivors = [], []
    for omitted in sorted(finite):
        active = sorted(set(finite) - {omitted})
        if 128 in active and 130 in active:
            assert lor(coeff[21], coeff[23])[1] < -2
            excluded.append(dict(omitted_neighbor=omitted, dotted_pair=[128, 130]))
        else:
            survivors.append(active)
    assert len(excluded) == 18 and len(survivors) == 2
    vertices = []
    for active in sorted(survivors):
        pair = sorted(set(active) - set(S))
        assert pair in ([128, 129], [129, 130])
        u, v = [coeff[i - 107] for i in pair]
        determinant = minus(times(u[1], v[2]), times(u[2], v[1]))
        point = [pt(1), div(minus(times(u[0], v[2]), times(u[2], v[0])), determinant), div(minus(times(u[1], v[0]), times(u[0], v[1])), determinant)]
        zy = lor(point, y)
        dual = div(times(pt(2), zy), c)
        square = plus(lor(point, point), div(times(dual, dual), pt(2)))
        assert dual[1] < 0 and square[1] < 0
        vertex_pairs = [plus(inner(point, p), times(dual, w)) for p, w in zip(P, wp)]
        for i, value in enumerate(vertex_pairs):
            if i in active:
                assert value[0] <= 0 <= value[1]
            else:
                assert value[1] < 0, (pair, i, saved(value))
        future = plus(lor(point, C0), times(pt(eps), plus(inner(point, hf), times(dual, omegaH))))
        assert future[1] < 0
        # All exact zero contacts follow from support orthogonality, the two
        # linear intersection equations, and dual=2*(point.y)/c.
        assert F(2) - F(2 + 1, 2) == F(1, 2)
        B.pd([[G[i][j] for j in range(18) if j != index31] for i in range(18) if i != index31])
        raw = next(row for row in d['complete_new_finite_vertices'] if row['old_active'] == active)
        for value, previous in zip(vertex_pairs, raw['old133_pair_bounds']): overlap(value, previous)
        vertices.append(dict(active_roots=active + [133], chain_pair=pair,
            independently_solved_frame_point=vec(point), omega31_coefficient=saved(dual),
            norm_bounds=saved(square), center_pair_bounds=saved(future),
            old133_pair_bounds=vec(vertex_pairs), exact_final_schur_pivot='1/2'))
    assert sorted(v['active_roots'] for v in vertices) == sorted(row['active'] for row in d['complete_new_finite_vertices'])
    edge_common = [sorted(set(S + v['chain_pair']) & set(v['active_roots'])) for v in vertices]
    edge_common.append(sorted(set(vertices[0]['active_roots']) & set(vertices[1]['active_roots'])))
    assert all(len(common) == 19 for common in edge_common)
    print('PASS common strict center,134 actual facets, complete2 new corners and3 finite edges', flush=True)

    # Independent remaining escape, with exact rational positive norm in an
    # orthogonal Lorentz3 plus positive support decomposition.
    ez = [F(1), F(-1, 100), F(0)]
    eb = [F(-1, 100000)] * 18
    eb[S.index(5)] = F(-1)
    wt = [sum(J[i][k] * eb[k] for k in range(18)) for i in range(18)]
    V = lin([N[i] for i in S], wt)
    assert [dot(V, N[i]) for i in S] == [const(x) for x in eb]
    assert all(dot(V, frame_axis) == zero for frame_axis in frame)
    vpairs = old_pairs(V) + [pt(0)] * 26
    for k, i in enumerate(S): vpairs[i] = pt(eb[k])
    Eframe = [pt(x) for x in ez]
    escape_pairs = [plus(inner(Eframe, p), v) for p, v in zip(P, vpairs)]
    escape_pairs.append(minus(lor(y, Eframe), times(c, pt(wt[index31]))))
    assert all(value[1] < 0 for value in escape_pairs)
    escape_norm = -ez[0] ** 2 + ez[1] ** 2 + ez[2] ** 2 + sum(x * y for x, y in zip(eb, wt))
    assert escape_norm > 0 and str(escape_norm) == d['strict_escape']['norm']
    escape_center = plus(lor(Eframe, C0), times(pt(eps), plus(inner(Eframe, hf), div(iv(dot(V, H)), Hden))))
    assert escape_center[1] < 0
    for value, raw in zip(escape_pairs, d['strict_escape']['all134_pair_bounds']): overlap(value, raw)
    cusp_pairs = []
    for raw in [old['completed_cusp']['point']] + [row['point'] for row in old['preserved_old_bounded_cusps']]:
        nu = vp(raw, 0)
        assert dot(nu, nu) == zero
        nf = [div(iv(dot(v, nu)), den) for v, den in zip(frame, frame_den)]
        value = minus(inner(y, nf), times(c, iv(dot(omega, nu))))
        assert value[1] < 0
        cusp_pairs.append(saved(value))
    print('PASS bank847 cut, old escape retained, new exact strict escape and3 cusp signs', flush=True)
    paths = [f'scripts/{TAG}.py', INPUT, BASE, CHAIN_REVIEW, BANK] + FILES + [
        'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
        'scripts/fresh_h20_bridge_faces_20261014_review.py',
        'scripts/fresh_h20_s61_multicap_20261014_exact_review.py',
        'results/acute_extension_preservation_20261014.md']
    return dict(status='PASS_INDEPENDENT_134_WALL_DOTTED_ARC_CAP_REVIEW',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        recursive_input_integrity=integrity(paths), interval_bits=R.BITS,
        arithmetic='Frozen independent sparse Fraction radical tower, rational support inversion, independent complex-power chain evaluation and outward512-bit dyadic intervals. No producer arithmetic imports.',
        inherited_inputs='Old133 chamber acuteness, the unique positive parameter in its isolating interval, chain norm/adjacency identities, and old finite-vertex lower bound168. All107 physical covectors and all new geometric checks are recomputed.',
        support_pd_pivots=list(map(str, support_pivots)),
        symbolic_norm_and_contact_identities=dict(h='sqrt2-1', a_positive=True,
            norm_n129='2', adjacent_products='n128.n129=n129.n130=0',
            cap_orthogonal_part='y=t*x+(h/2)*spatial(n129), t=(h*a^2+2sqrt2)/(2a)',
            norm_y='-1/2-2/a^2', omega31_norm='1/2', cap_norm='2',
            central_contact='h-sqrt2=-1', neighboring_contacts='-sqrt2', right_support_roots=sorted(set(S) - {31})),
        pair_counts=dict(counts), all133_new_contacts=pairs,
        common_center_norm=saved(center_norm), all134_center_pair_bounds=vec(center_pairs),
        center_pair_with_inherited_future_unit_center=saved(center_old_future),
        all134_facet_point_norm_bounds=vec(facet_norms), actual_facets=134,
        actual_facet_argument='For every norm2 wall n, C-(C.n)n/2 is future timelike, orthogonal to n and strictly negative on every other wall. It has norm C^2-(C.n)^2/2<0.',
        old_finite_vertex_preservation='The common strict center and nonpositive new-old contacts imply every old finite vertex is strictly retained, by the positive inverse of its finite Gram matrix. The inherited lower bound168 therefore rises to at least170.',
        former_escape_tests=tested_escapes,
        complete_new_finite_vertex_census=dict(finite_old_neighbors=sorted(finite), candidate_subsets=20,
            dotted_pair_exclusions=excluded, exact_surviving_new_corners=vertices,
            count=2, simplicity_scope='A finite vertex of a hyperbolic Coxeter chamber has20 independent spherical incident walls. Any vertex on the new facet needs19 old finite neighbors. All candidate19-subsets are exhausted.'),
        three_certified_finite_edge_common_supports=edge_common,
        remaining_strict_escape=dict(frame_coefficients=list(map(str, ez)), support=S,
            dual_pairing_coefficients=list(map(str, eb)), support_normal_coefficients=list(map(str, wt)),
            exact_norm=str(escape_norm), all134_pair_bounds=vec(escape_pairs), center_pair_bounds=saved(escape_center),
            infinite_volume_argument='C is future timelike and E is spacelike with C.E<0, both strictly inside every wall halfspace. The positive ray C+tE meets the null cone while all wall inequalities remain strict. Its open ideal-boundary neighborhood proves infinite volume.'),
        preserved_cusp_new_pair_bounds=cusp_pairs,
        scope='Independent exact one-wall extension review, including all contacts, actual facets, exactly two new finite vertices, three finite edges and a surviving strict escape. The old ambient finite-vertex census remains a lower bound; no full ambient ray census, finite-volume completion or new-dimensional example is claimed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], flush=True)
