#!/usr/bin/env python3
"""Independent sparse-radical and interval review of the complementary H20 cap.

No producer arithmetic imports. The old parameter's existence/isolation and
old chamber acuteness are inherited from the frozen independent chain review.
Every new contact, center, finite-corner exclusion and remaining escape is replayed.
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
TAG = 'fresh_h20_complementary_cap_20261016_review'
INPUT = 'results/fresh_h20_complementary_cap_20261016.json'
PRIOR = 'results/fresh_h20_dotted_arc_cap_20261016.json'
PRIOR_REVIEW = 'results/fresh_h20_dotted_arc_cap_20261016_review.json'
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
    d, base, reviewed, prior, prior_review = [json.loads((ROOT / p).read_text()) for p in (INPUT, BASE, CHAIN_REVIEW, PRIOR, PRIOR_REVIEW)]
    assert d['status'] == 'PASS_EXACT_135_WALL_COMPLEMENTARY_E0_CUTTER'
    assert prior_review['sources'][PRIOR] == sha256((ROOT / PRIOR).read_bytes()).hexdigest()
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
    index5 = S.index(5)
    assert J[index5][index5] == F(30,11)
    omega = lin([N[i] for i in S], J[index5])
    assert dot(omega,omega) == const(F(30,11))
    assert [dot(omega,N[i]) for i in S] == [const(i==5) for i in S]
    assert dot(omega,N[0]) == const(F(-6,11))
    assert norms[0] == const(2)
    assert dot(omega,N[95]) == zero
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
        overlap(wp[i], [int(i == 5)] * 2)
        wp[i] = pt(i == 5)
    print('PASS independent field frame,107 covectors,26 chain normals and omega5', flush=True)

    # Anchor95 has zero e-coordinate contact exactly; anchor0 has pi/6.
    assert frame[1] == N[95] and dot(N[95],N[95]) == const(2)
    P[95] = [pt(0),root(pt(2)),pt(0)]
    wp[0] = pt(F(-6,11)); wp[95] = pt(0)
    p,q = P[0][0],P[0][2]
    assert p[1]<0 and q[1]<0
    h = minus(pt(F(-6,11)),root(pt(3)))
    alpha,beta = div(h,p),neg(div(q,p))
    qa = minus(pt(1),times(beta,beta))
    qb = neg(times(pt(2),times(alpha,beta)))
    qc = minus(pt(F(8,11)),times(alpha,alpha))
    disc = minus(times(qb,qb),times(pt(4),times(qa,qc)))
    assert qa[0]>0 and disc[0]>0
    z = div(plus(neg(qb),root(disc)),times(pt(2),qa))
    zminus = div(minus(neg(qb),root(disc)),times(pt(2),qa))
    t = plus(alpha,times(beta,z))
    y = [t,pt(0),z]
    assert t[0]>0 and z[0]>0 and zminus[1]<0
    for value,raw in zip(y,d['definition']['frame_bounds']):overlap(value,raw)
    # Symbolic equations certify exact norm/finite contacts, beyond intervals.
    import sympy as sy
    aa,bb,zz,pp,qq,hh=sy.symbols('aa bb zz pp qq hh')
    polynomial=(1-bb**2)*zz**2-2*aa*bb*zz+sy.Rational(8,11)-aa**2
    assert sy.expand((aa+bb*zz)**2-zz**2-sy.Rational(8,11)+polynomial)==0
    assert sy.simplify(pp*(hh/pp-qq*zz/pp)+qq*zz-hh)==0
    assert R.overlap(lor(y,y),pt(F(-8,11)))
    assert F(-8,11)+F(30,11)==2
    finite={i:2 for i in S}|{5:3,0:6,95:2}
    pairs=[];counts=Counter()
    for i,(cov,w) in enumerate(zip(P,wp)):
        value=minus(inner(y,cov),w)
        if i in finite:
            m=finite[i]
            expected=pt(0) if m==2 else pt(-1) if m==3 else neg(root(pt(3)))
            assert R.overlap(value,expected)
            relation='right' if m==2 else f'pi/{m}'
        else:
            assert value[1]<-2,(i,saved(value))
            relation='dotted'
        assert relation==d['pairs'][i]['relation']
        overlap(value,d['pairs'][i]['pair_bounds'])
        pairs.append(dict(wall=i,relation=relation,bounds=saved(value)))
        counts[relation]+=1
    # Reconstruct the preceding arc cap frame from its independently reviewed recipe.
    sqrt2=root(pt(2));h2=minus(sqrt2,pt(1));a2=times(a,a)
    y133=[div(plus(times(h2,a2),times(pt(2),sqrt2)),times(pt(2),a)),
          times(div(h2,pt(2)),coeff[22][1]),times(div(h2,pt(2)),coeff[22][2])]
    c133=root(plus(pt(5),div(pt(4),a2)))
    assert J[index5][S.index(31)]==0
    mutual=lor(y,y133)
    assert mutual[1]<-2
    overlap(mutual,d['mutual_arc_cap_pair_bounds'])
    pairs.append(dict(wall=133,relation='dotted',bounds=saved(mutual)));counts['dotted']+=1
    assert counts==dict(right=18,**{'pi/3':1,'pi/6':1,'dotted':114})
    assert sorted(finite)==d['finite_neighbors']
    print('PASS all134 independent contacts and symbolic finite equalities',flush=True)

    H=vp(cur['strict_center']['point'],1)
    Hnorm=dot(H,H);assert ff.sign(Hnorm)<0
    Hden=root(iv(sc(Hnorm,-1)))
    hf=[div(iv(dot(v,H)),times(den,Hden)) for v,den in zip(frame,frame_den)]
    omegaH=div(iv(dot(omega,H)),Hden)
    C0=[pt(1),pt(F(-1,2)),pt(0)];eps=F(base['strict_center']['epsilon'])
    newcenter=plus(lor(y,C0),times(pt(eps),minus(inner(y,hf),omegaH)))
    assert newcenter[1]<0
    center_norm=minus(plus(pt(F(-3,4)),times(pt(2*eps),inner(C0,hf))),pt(eps*eps))
    center_pairs=[tuple(map(F,row)) for row in prior['common_center']['all134_pair_bounds']]+[newcenter]
    assert center_norm[1]<0 and all(v[1]<0 for v in center_pairs)
    for value,raw in zip(center_pairs,d['common_center']['all135_pair_bounds']):overlap(value,raw)
    facets=[minus(center_norm,div(times(v,v),pt(2))) for v in center_pairs]
    assert all(v[1]<0 for v in facets)
    # A projection to a norm2 wall stays timelike; acuteness makes all its
    # other pairings strict. All135 walls are therefore actual facets.
    e0=plus(lor(y,[pt(1),pt(F(-1,100)),pt(0)]),pt(J[index5][index5]))
    assert e0[0]>0;overlap(e0,d['weak_escape_cut']['pair_bounds'])
    prior_eb=list(map(F,prior['strict_escape']['dual_coefficients']))
    oldpay=sum(J[index5][k]*prior_eb[k] for k in range(18))
    old_escape=minus(lor(y,[pt(1),pt(F(-1,100)),pt(0)]),pt(oldpay))
    assert old_escape[0]>0;overlap(old_escape,d['previous_strict_escape_cut'])

    # Exhaust the20 possible19-old-neighbor subsets using explicit witnesses.
    assert dot(N[0],N[1])==const(-1)
    assert ff.sign(add(dot(N[0],N[87]),const(2)))<0
    affine=sy.Matrix([[2,-sy.sqrt(3),0],[-sy.sqrt(3),2,-1],[0,-1,2]])
    kernel=sy.Matrix([sy.sqrt(3),2,1])
    assert affine*kernel==sy.zeros(3,1) and affine.rank()==2
    reasons=[]
    for subset in combinations(sorted(finite),19):
        active=set(subset)
        if {0,87}<=active:reason='dotted0-87'
        elif set(S)<=active:
            assert F(2)-J[index5][index5]==F(-8,11)
            reason='negative_Schur_S_plus_cap'
        else:
            assert {0,1}<=active
            reason='affine_cap_0_1'
        reasons.append(dict(omitted=next(iter(set(finite)-active)),reason=reason))
    assert len(reasons)==20
    assert Counter(r['reason'] for r in reasons)==dict({'dotted0-87':18,'negative_Schur_S_plus_cap':1,'affine_cap_0_1':1})

    # Recompute the replacement strict escape from literal physical covectors.
    ez=list(map(F,d['strict_escape']['frame_coefficients']))
    eb=list(map(F,d['strict_escape']['dual_coefficients']))
    assert ez==[1,F(-1,100),F(-49,50)]
    assert eb==[F(-1) if wall==17 else F(-1,100000) for wall in S]
    weights=[sum(J[i][k]*eb[k] for k in range(18)) for i in range(18)]
    V=lin([N[i] for i in S],weights)
    assert all(dot(v,V)==zero for v in frame)
    ve=old_pairs(V)+[pt(0)]*26
    for k,i in enumerate(S):
        assert dot(V,N[i])==const(eb[k]);ve[i]=pt(eb[k])
    escape_pairs=[plus(inner([pt(x) for x in ez],cov),w) for cov,w in zip(P,ve)]
    escape_pairs += [minus(lor(y133,[pt(x) for x in ez]),times(c133,pt(weights[S.index(31)]))),
                     minus(lor(y,[pt(x) for x in ez]),pt(weights[index5]))]
    assert len(escape_pairs)==135 and all(v[1]<0 for v in escape_pairs)
    for value,raw in zip(escape_pairs,d['strict_escape']['all135_pair_bounds']):overlap(value,raw)
    norm=-ez[0]**2+ez[1]**2+ez[2]**2+sum(eb[i]*weights[i] for i in range(18))
    assert norm==F(1842000049,4000000000)>0
    eh=plus(inner([pt(x) for x in ez],hf),div(iv(dot(V,H)),Hden))
    ec=plus(lor([pt(x) for x in ez],C0),times(pt(eps),eh))
    assert ec[1]<0;overlap(ec,d['strict_escape']['strict_center_pair_bounds'])
    cusp_pairs=[]
    for raw in [old['completed_cusp']['point']]+[r['point'] for r in old['preserved_old_bounded_cusps']]:
        nu=vp(raw,0)
        components=[div(iv(dot(v,nu)),den) for v,den in zip(frame,frame_den)]
        value=minus(inner(y,components),iv(dot(omega,nu)))
        assert value[1]<0;cusp_pairs.append(saved(value))
    paths=[f'scripts/{TAG}.py',INPUT,BASE,CHAIN_REVIEW,PRIOR,PRIOR_REVIEW]+FILES+[
        'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
        'scripts/fresh_h20_bridge_faces_20261014_review.py',
        'scripts/fresh_h20_s61_multicap_20261014_exact_review.py',
        'results/acute_extension_preservation_20261014.md']
    audited=integrity(paths)
    return dict(status='PASS_INDEPENDENT_135_WALL_COMPLEMENTARY_CAP_AND_REPLACEMENT_ESCAPE',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                recursive_manifest_audit=audited,frame=vec(y),pairs=pairs,
                actual_facets=135,pair_counts=dict(counts),center_pair=saved(newcenter),
                E0_cut=saved(e0),complete_new_finite_vertices=0,
                all20_finite_neighbor_exclusions=reasons,three_preserved_cusp_pairings=cusp_pairs,
                strict_escape_norm=str(norm),strict_escape_pairs=vec(escape_pairs),
                scope='Independent sparse-radical physical normals and covectors, chain recurrence, symbolic defining quadratic and finite contacts, all134 new Coxeter gates, coherent135actualfacets, all20 finite-corner exclusions and remaining strict escape. Inherits old134geometry/vertices/cusps and acute-extension preservation. Atleast170finitevertices; no ambientcompletecensus or finite-volume example.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)

