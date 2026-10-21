#!/usr/bin/env python3
"""Independent grounded-rational and symbolic-coordinate review of charge4/3.

No producer arithmetic is imported. Published geometry remains a mathematical
dependency, explicitly bound through the frozen181 case binder.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import sympy as sy
import fresh_nikulin_energy_charge_20261015_review as H

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nikulin_unit_charge_20261015_review'
PRIMARY = 'fresh_nikulin_unit_charge_20261015'
ALPHA = Q(4, 3)


def model(A, d2):
    n = len(A)
    ids = list(range(1, n))
    J = H.inv(H.sub(A, ids))
    ker = [Q(1)] + H.mv(J, [-A[i][0] for i in ids])
    assert H.mv(A, ker) == [0] * n and min(ker) > 0
    scale = min(d*x*x for d, x in zip(d2, ker))
    delta2 = [d*x*x/scale for d, x in zip(d2, ker)]
    assert min(delta2) == 1
    inverses, diagonals = [], []
    for omitted in range(n):
        retained = [i for i in range(n) if i != omitted]
        B = H.sub(A, retained)
        H.pd(B)
        K = H.inv(B)
        assert all(x >= 0 for row in K for x in row)
        inverses.append((retained, K))
        diagonals.append({i: d2[i]*K[k][k] for k, i in enumerate(retained)})

    @lru_cache(None)
    def resistance(i, j):
        b = [Q(0)]*n
        b[i] += 1/ker[i]
        b[j] -= 1/ker[j]
        potential = [Q(0)] + H.mv(J, b[1:])
        assert H.mv(A, potential) == b and H.dot(b, ker) == 0
        return scale*H.dot(b, potential)

    def cross(omitted, i, j):
        ids, K = inverses[omitted]
        square = d2[i]*d2[j]
        numerator, denominator = isqrt(square.numerator), isqrt(square.denominator)
        assert Q(numerator*numerator, denominator*denominator) == square
        return Q(numerator, denominator)*K[ids.index(i)][ids.index(j)]

    return delta2, diagonals, resistance, cross


def physical_gram(A, d2, saved):
    for i, row in enumerate(A):
        for j, entry in enumerate(row):
            text = saved[i][j]
            if text in ('-sqrt(2)', '-sqrt(3)'):
                sign, square = -1, Q(2 if text == '-sqrt(2)' else 3)
            else:
                value = Q(text)
                sign, square = (value > 0)-(value < 0), value*value
            assert sign == (entry > 0)-(entry < 0)
            assert entry*entry/(d2[i]*d2[j]) == square


def minimum_quadratic(coeff):
    a, b, c = coeff
    assert c >= 0
    candidates = [a, a+b+c]
    if c and 0 <= -b/(2*c) <= 1:
        candidates.append(a-b*b/(4*c))
    return min(candidates)


def exceptions(saved):
    literals = {
        'A1': ([[2, -2], [-2, 2]], [Q(1)]*2),
        'E6': (H.star([2, 2, 2]), [Q(1)]*7),
        'E7': (H.star([1, 3, 3]), [Q(1)]*8),
        'E8': (H.star([1, 2, 5]), [Q(1)]*9),
        'F4': ([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-2,0],
                [0,0,-2,4,-2],[0,0,0,-2,4]], [Q(1)]*3+[Q(2)]*2),
        'G2': ([[2,-1,0],[-1,2,-3],[0,-3,6]], [Q(1),Q(1),Q(3)])}
    ordinary = {'A1':(0,1),'E6':(2,4),'E7':(3,6),'E8':(1,6),'F4':(0,4),'G2':(0,2)}
    rows, censuses = [], []
    single_count = two_count = 0
    for name, (A, d2) in literals.items():
        delta2, diagonal, R, cross = model(A, d2)
        r = len(A)-1
        for forbidden in [None]+list(range(r+1)):
            if forbidden is None:
                omits = ordinary[name]
            elif name == 'A1':
                omits = (1-forbidden,)*2
            elif name == 'E6':
                omits = tuple(x for x in (2,4,6) if x != forbidden)[:2]
            elif name == 'E7':
                omits = (2,5) if forbidden in (4,7) else (4,7)
            elif name == 'E8':
                omits = (1,6) if forbidden in (3,7) else (3,7)
            elif name == 'F4':
                omits = (0,4) if forbidden in (1,3) else (1,3)
            else:
                omits = tuple(x for x in range(3) if x != forbidden)
            assert forbidden not in omits
            target = Q(2*(r-(forbidden is not None)),9)
            record = next(x for x in saved['exceptional_cases'] if x['type'] == name and x['forbidden'] == forbidden)
            assert record['omissions'] == list(omits) and Q(record['target']) == target
            physical_gram(A, d2, record['gram'])
            assert sy.Matrix([sy.sympify(x)**2 for x in record['kernel']]) == sy.Matrix(delta2)
            diags = [diagonal[i] for i in omits]
            assert [{str(k):str(v) for k,v in d.items()} for d in diags] == record['inverse_diagonals']
            L = min(delta2[i] for i in omits)*R(*omits)
            assert L == Q(record['point_energy']) and L >= target
            limit = target/ALPHA
            cheap = [(side,k,c) for side,dd in enumerate(diags) for k,c in dd.items() if c < limit]
            assert [(x['side'],x['node'],Q(x['cost'])) for x in record['single_contact_cases']] == cheap
            checks = []
            for side, k, cost in cheap:
                i,j = omits[side], omits[1-side]
                triple = [R(i,j),R(k,j),R(i,k)]
                coeff = [triple[0],-triple[0]+triple[1]-triple[2],triple[2]]
                minimum = minimum_quadratic(coeff)
                zero = min(delta2[k],delta2[j])*triple[1]+ALPHA*cost
                positive = min(delta2[i],delta2[j])*minimum+ALPHA*cost
                assert min(zero,positive) >= target
                other = next(x for x in record['single_contact_cases'] if (x['side'],x['node']) == (side,k))
                assert list(map(Q,other['resistances'])) == triple
                assert Q(other['minimum_energy']) == minimum
                assert Q(other['zero_omitted_total_lower']) == zero
                assert Q(other['positive_omitted_total_lower']) == positive
                checks.append(dict(side=side,node=k,cost=str(cost),energy_coefficients=list(map(str,coeff)),
                                   exact_quadratic_minimum=str(minimum),zero_bound=str(zero),positive_bound=str(positive)))
            # Enumerate every possible support, not just every pair. Monotonicity
            # in nonzero magnitudes follows from entrywise nonnegative inverse.
            surviving = []
            for size in range(2,len(cheap)+1):
                for support in combinations(cheap,size):
                    cost = sum((c for _,_,c in support),Q(0))
                    for (s,i,_),(t,j,_) in combinations(support,2):
                        if s == t:cost += 2*cross(omits[s],i,j)
                    if cost < limit:surviving.append((support,cost))
            if surviving:
                assert name == 'E7' and forbidden is None
                assert surviving == [(((0,4,Q(1,2)),(1,7,Q(1,2))),Q(1))]
                # Polarize the grounded resistance matrix to recover all six
                # coefficients of the two-mixture energy independently.
                nodes = [3,4,6,7]
                u,v = sy.symbols('u v')
                weights = [1-u,u,-(1-v),-v]
                E = sy.expand(-sum(weights[a]*weights[b]*R(nodes[a],nodes[b]) for a,b in combinations(range(4),2)))
                assert E == (1+u*u+v*v)/2 and Q(1,2)+ALPHA >= target
                assert record['two_contact_case'] is not None
                two_count += 1
            else:
                assert record['two_contact_case'] is None
            single_count += len(checks)
            rows.append(dict(type=name,forbidden=forbidden,omissions=list(omits),target=str(target),
                             point_energy=str(L),single_contacts=checks,nontrivial_multiple_supports=len(surviving)))

        if name in ('E6','E7','E8'):
            source = next(x for x in saved['complete_point_row_omission_censuses'] if x['type'] == name)
            eligible, excluded = [], []
            target = Q(2*r,9)
            for i,j in combinations(range(r+1),2):
                L = min(delta2[i],delta2[j])*R(i,j)
                if L < target:
                    excluded.append((i,j));continue
                costs = [{i:Q(0),**diagonal[i]},{j:Q(0),**diagonal[j]}]
                worst = Q(0)
                for a,b in product(range(r+1),repeat=2):
                    c,d = costs[0][a],costs[1][b]
                    if max(c,d)>2 or not c+d:continue
                    W = min(delta2[a],delta2[b])*R(a,b)
                    worst = max(worst,(target-W)/(c+d))
                primary = next(x for x in source['eligible_pairs'] if x['omissions'] == [i,j])
                assert Q(primary['required_charge_lower']) == worst
                eligible.append(dict(omissions=[i,j],necessary_coefficient=str(worst)))
            assert len(eligible)+len(excluded) == r*(r+1)//2
            assert [list(x) for x in excluded] == [x['omissions'] for x in source['ineligible_zero_cost_pairs']]
            assert len(eligible) == len(source['eligible_pairs'])
            censuses.append(dict(type=name,all_pairs=r*(r+1)//2,eligible=eligible,
                                 best_necessary_coefficient=str(min(Q(x['necessary_coefficient']) for x in eligible))))
    assert len(rows) == 40 and single_count == 46 and two_count == 1
    assert [x['best_necessary_coefficient'] for x in censuses] == ['1','10/9','4/3']
    return rows,censuses


t,u,v = sy.symbols('t u v',nonnegative=True,integer=True)
r,k = sy.symbols('r k',integer=True,positive=True)


def nonnegative(expr):
    num,den = sy.fraction(sy.cancel(expr))
    return all(x>=0 for x in sy.Poly(num,t,u,v).coeffs()) and all(x>=0 for x in sy.Poly(den,t,u,v).coeffs()) and den.subs({t:0,u:0,v:0})>0


def compare(a,b):
    if nonnegative(a-b):return 1
    assert nonnegative(b-a),(a,b)
    return -1


def potential(typ,node):
    if node == 0:return []
    if typ == 'C':return [(1,node,sy.S(1))]
    if node == 1:return [(1,1,sy.S(1))]
    if typ == 'D' and node == r-1:
        return [(1,r-1,sy.Rational(1,2)),(r,r,-sy.Rational(1,2))]
    return [(1,node,sy.Rational(1,2))]


def coordinate_distance(typ,a,b,subs):
    segments = potential(typ,a)+[(lo,hi,-c) for lo,hi,c in potential(typ,b)]
    segments = [(sy.sympify(lo).subs(subs),sy.sympify(hi).subs(subs),c) for lo,hi,c in segments]
    value = 0
    for lo,hi,c in segments:
        for low,high,d in segments:
            lower = lo if compare(lo,low)>0 else low
            upper = hi if compare(hi,high)<0 else high
            length = sy.expand(upper-lower+1)
            if nonnegative(length):value += c*d*length
            else:assert nonnegative(-length)
    return sy.factor(value/(2 if typ == 'C' else 1))


def classical():
    rows = []
    def check(typ,modified,i,j,node,d2,subs):
        Rij = coordinate_distance(typ,i,j,subs)
        Rik = coordinate_distance(typ,i,node,subs)
        Rkj = coordinate_distance(typ,node,j,subs)
        cost,loss = sy.factor(d2*Rik),sy.factor(Rij+Rik-Rkj)
        target = sy.Rational(2,9)*(r-int(modified)).subs(subs)
        assert nonnegative(Rij-target)
        margin = sy.factor(sy.Rational(4,3)*cost-loss)
        excluded = sy.factor(cost-sy.Rational(3,4)*target)
        if nonnegative(margin):reason='charge covers loss';proof=margin
        else:
            assert nonnegative(excluded) and excluded.subs({t:0,u:0,v:0})>0
            reason='contact exceeds trivial-charge threshold';proof=excluded
        rows.append(dict(type=typ,modified=modified,omissions=list(map(str,(i,j))),contact=str(node),
                         rank_domain={str(a):str(b) for a,b in subs.items()},point_distance=str(Rij),
                         diagonal=str(cost),linear_loss=str(loss),proof=reason,nonnegative_margin=str(proof)))
    # The indices below partition both ordinary B rows and both modified B rows.
    for i,j in ((0,r),(r,0)):
        for node,d2,subs in [(0,1,{r:2+t}),(1,1,{r:2+t}),(k,4,{k:2+u,r:3+u+v}),(r,2,{r:2+t})]:
            if node != i:check('B',False,i,j,node,d2,subs)
    for i,j in ((r-1,0),(0,r-1)):
        items = [(0,1,{r:3+t}),(1,1,{r:3+t}),(r,2,{r:3+t})]
        items += [(k,4,{k:2+u,r:4+u+v})] if i == r-1 else [(k,4,{k:2+u,r:3+u+v})]
        for node,d2,subs in items:
            if node != i:check('B',True,i,j,node,d2,subs)
    # Ordinary C's other row follows by path reversal.
    check('C',False,0,r,k,2,{k:1+u,r:2+u+v})
    check('C',False,0,r,r,1,{r:2+t})
    for i,j in ((1,r),(r,1)):
        items = [(0,1,{r:2+t})]
        if i == 1:items += [(k,2,{k:2+u,r:3+u+v}),(r,1,{r:2+t})]
        else:items += [(k,2,{k:1+u,r:2+u+v})]
        for node,d2,subs in items:check('C',True,i,j,node,d2,subs)
    # D's second row and all avoid-node leaf choices follow by fork symmetries.
    for node,d2,subs in [(1,1,{r:4+t}),(k,4,{k:2+u,r:4+u+v}),
                         (r-1,1,{r:4+t}),(r,1,{r:4+t})]:
        check('D',False,0,r,node,d2,subs)
    assert len(rows) == 24
    # Separate B2 modification avoids the coincident special prefix indices.
    A = [[2,0,-2],[0,2,-2],[-2,-2,4]]
    ds,diags,R,_ = model(A,[Q(1),Q(1),Q(2)])
    assert R(0,1) == 1 and min(x for i in (0,1) for x in diags[i].values()) >= Q(1,2)
    assert ALPHA*Q(1,2) >= Q(2,9)

    # Independent cycle derivation: the two directions yield losses d(m±eps)/m.
    m,d,e = sy.symbols('m d e',positive=True)
    h = (m-e)/2
    resist = lambda x:x*(m-x)/m
    assert sy.expand(resist(h)+resist(d)-resist(h-d)-(m+e)*d/m) == 0
    assert sy.expand(resist(h)+resist(d)-resist(m-h-d)-(m-e)*d/m) == 0
    assert sy.cancel(sy.Rational(4,3)*resist(d)-(m+e)*d/m-d*(m-4*d-3*e)/(3*m)) == 0
    small = []
    for rank in range(2,6):
        assert Q(rank,rank+1) >= Q(rank,6)
        small.append(dict(rank=rank,retained_contacts_below_threshold=0))
    for order in (7,8,9):
        assert Q(2*(order-2),order) >= Q(order-1,6)
        assert order-4-3*(order%2) >= 0
        small.append(dict(rank=order-1,allowed_distances=[1]))
    cycle = []
    for order,epsilon in ((10+2*t,0),(11+2*t,1)):
        assert nonnegative(order-9) and (order-9).subs(t,0)>0
        for distance in (1,2):
            margin = sy.expand(distance*(order-4*distance-3*epsilon))
            assert nonnegative(margin)
            cycle.append(dict(order=str(order),distance=distance,cleared_margin=str(margin)))
    # Point energy exceeds 2r/9, for both parities and all ranks.
    for epsilon,order in ((0,2+2*t),(1,3+2*t)):
        margin = sy.factor((order*order-epsilon)/(4*order)-sy.Rational(2,9)*(order-1))
        assert nonnegative(margin)
    return dict(prefix_coordinate_classes=rows,cycle_small_ranks=small,cycle_unbounded=cycle,
                avoid_node_symmetries='Cycle rotation, B left-fork switch, C reversal and D fork/end automorphisms; modified B/C rows explicit.')


def run():
    saved = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    rows,censuses = exceptions(saved)
    allranks = classical()
    # Verify the geometric identity's algebra without assuming equal heights.
    Rij,Rik,Rkj,Rjl,Rli,Rkl,w,z = sy.symbols('Rij Rik Rkj Rjl Rli Rkl w z')
    # |a-b|² = w²Rik+z²Rjl-wz(Rkj+Rli-Rkl-Rij).
    remainder = w*w*Rik+z*z*Rjl-w*z*(Rkj+Rli-Rkl-Rij)
    linear = Rij-w*(Rij+Rik-Rkj)-z*(Rij+Rjl-Rli)+remainder
    mixture = ((1-w)*(1-z)*Rij+w*(1-z)*Rkj+(1-w)*z*Rli+w*z*Rkl
               -w*(1-w)*Rik-z*(1-z)*Rjl)
    assert sy.expand(linear-mixture) == 0
    assert (Q(2*169,9)-ALPHA*4-Q(2,9))/2-2 == 14
    assert (Q(2*168,9)-ALPHA*4-Q(2,9))/2-2 < 14
    assert saved['combined_candidate_dimension_at_most'] == 169
    files = [f'scripts/{TAG}.py','scripts/fresh_nikulin_energy_charge_20261015_review.py']
    for tag in (PRIMARY,'fresh_nikulin_dimension_181_candidate_20261015'):
        files += [f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md']
    digest = lambda p:sha256((ROOT/p).read_bytes()).hexdigest()
    for p in files:
        if p.endswith('.json'):
            for f,h in json.loads((ROOT/p).read_text())['sources'].items():assert digest(f) == h,(p,f)
    return dict(status='PASS_INDEPENDENT_GROUNDED_ENERGY_ALL_RANK_AND_FOUR_THIRDS_REVIEW',
                sources={p:digest(p) for p in files},exceptional_cases=rows,omission_censuses=censuses,
                all_rank_coordinate_proof=allranks,charge_coefficient='4/3',
                complete_exceptional_cases=40,single_contact_cases=46,two_contact_supports=1,
                total_schur_cost_bound=4,avoid_node_penalty='2/9',first_excluded_dimension=170,
                ideal_branch_dimension_at_most=169,compact_facet_branch_dimension_at_most=143,
                combined_candidate_dimension_at_most=169,
                scope='Independent Fraction grounded solves, all exceptional support subsets, complete E omission census, symbolic interval-coordinate classes and analytic proof review. Inherits the181 binder for actual framing walls, common-height identity, published narrowness and compact143 branch. Internal research review, no external validation, formal proof-assistant certification or novelty claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    data = run()
    destination = ROOT/f'results/{TAG}.json'
    if args.verify:assert data == json.loads(destination.read_text())
    else:destination.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['combined_candidate_dimension_at_most'])
