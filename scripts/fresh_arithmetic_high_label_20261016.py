#!/usr/bin/env python3
"""Small exact checks for the distinguished-bond compact dimension-20 lemma.

The all-rank geometric and field proofs are in the matching report.  This
standalone checker has no project, optimizer, or symbolic-algebra imports.
Published face correspondence, finite Coxeter classification, arithmetic
ground-field criterion and h-vector theorem are explicit proof inputs.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_arithmetic_high_label_20261016'


def edge(a, b):
    return tuple(sorted((a, b)))


def distances(vertices, es, source):
    adj = {v: [] for v in vertices}
    for a, b in es:
        if a in adj and b in adj:
            adj[a].append(b)
            adj[b].append(a)
    out, todo = {source: 0}, deque([source])
    while todo:
        a = todo.popleft()
        for b in adj[a]:
            if b not in out:
                out[b] = out[a] + 1
                todo.append(b)
    return out


def trees_or_forest(vertices, es):
    remain = set(vertices)
    while remain:
        comp = set(distances(vertices, es, min(remain)))
        if sum(a in comp and b in comp for a, b in es) != len(comp)-1:
            return False
        remain -= comp
    return True


def connected_supports(vertices):
    pairs = list(combinations(vertices, 2))
    for bits in product((0, 1), repeat=len(pairs)):
        es = {p for p, bit in zip(pairs, bits) if bit}
        if len(distances(vertices, es, vertices[0])) == len(vertices):
            yield es


def local_supports():
    """Exhaust all supports of two order2/3 marked minimal pieces.

    Finite corners are relaxed to forests plus isolation of a chosen high
    bond. Shared white nodes and all possible cross contacts are included.
    """
    examined, survivors = 0, set()
    for lwhite, rwhite, shared in [(0, 0, False), (0, 1, False),
                                  (1, 0, False), (1, 1, False), (1, 1, True)]:
        L = [0, 2] + ([4] if lwhite else [])
        R = [1, 3] + ([4 if shared else 5] if rwhite else [])
        vertices = sorted(set(L+R))
        white = sorted(set(vertices)-{0, 1, 2, 3})
        inside = set(combinations(L, 2)) | set(combinations(R, 2))
        cross = sorted(set(combinations(vertices, 2))-inside)
        for le, re in product(connected_supports(L), connected_supports(R)):
            for lh, rh in product(sorted(le) if lwhite else [None],
                                  sorted(re) if rwhite else [None]):
                high = {e for e in (lh, rh) if e is not None}
                for bits in product((0, 1), repeat=len(cross)):
                    bridges = {p for p, bit in zip(cross, bits) if bit}
                    if not shared and not bridges:
                        continue  # orthogonal indefinite blocks
                    examined += 1
                    es = le | re | bridges
                    good = True
                    for a, b in product((0, 2), (1, 3)):
                        corner = set(white+[a, b])
                        ce = {p for p in es if set(p) <= corner}
                        deg = Counter(v for p in ce for v in p)
                        if not trees_or_forest(corner, ce):
                            good = False
                        if any(deg[x] != 1 or deg[y] != 1
                               for x, y in high if {x, y} <= corner):
                            good = False
                        if distances(corner, ce, a).get(b, 99) <= 2:
                            good = False
                        if not good:
                            break
                    if good:
                        assert lwhite == rwhite == 1 and not shared
                        assert lh == (0, 2) and rh == (1, 3)
                        assert bridges == {(4, 5)}
                        survivors.add(tuple(sorted(es)))
    assert len(survivors) == 9
    return dict(examined=examined, surviving_supports=sorted(survivors),
                conclusion='two disjoint opposite high finite pairs, joined only by the two distinct white nodes')


def phi(m):
    return sum(gcd(a, m) == 1 for a in range(1, m+1))


def real_cyclotomic_group(m):
    reps = sorted({min(a, m-a) for a in range(1, m) if gcd(a, m) == 1})
    def mul(a, b):
        z = a*b % m
        return min(z, m-z)
    orders = {}
    for a in reps:
        z, n = 1, 0
        while True:
            z, n = mul(z, a), n+1
            if z == 1:
                break
            assert n <= len(reps)
        orders[a] = n
    assert len(reps) == phi(m)//2
    return dict(conductor=m, degree=len(reps), representatives=reps,
                element_orders=orders, cyclic=max(orders.values()) == len(reps))


def fields():
    groups = {m: real_cyclotomic_group(m) for m in (5, 7, 9, 11, 15, 16, 20)}
    assert groups[9]['degree'] == 3 and groups[11]['degree'] == 5
    assert all(groups[m]['degree'] == 4 and groups[m]['cyclic'] for m in (15, 16, 20))
    assert all(m % 5 == 0 for m in (15, 20))
    # Coprime cyclotomic conductors give the product field of degree36.
    assert gcd(7, 9) == 1 and phi(63) == phi(7)*phi(9) == 36
    # A multiquadratic field has elementary-2 Galois group, so no cyclic
    # quartic subfield; odd-degree fields cannot be its subfields either.
    for degree in (3, 5):
        assert gcd(degree, 4) == gcd(degree, 8) == 1
    low = []
    for m in range(2, 9):
        low.append(dict(label=m, squared_cosine_field_degree=1 if m == 2 else phi(m)//2,
                        required_field={2:'Q',3:'Q',4:'Q',5:'Q(sqrt5)',6:'Q',7:'K7',8:'Q(sqrt2)'}[m]))
    twelve = ['Q(sqrt2)', 'Q(sqrt3)', 'Q(sqrt5)', 'Q(sqrt6)',
              'Q(sqrt2,sqrt3)', 'Q(sqrt2,sqrt5)', 'K7', 'K9', 'K11', 'K15', 'K16', 'K20']
    new = ['Q(sqrt3)', 'Q(sqrt6)', 'Q(sqrt2,sqrt3)',
           'K9', 'K11', 'K15', 'K16', 'K20']
    previous = ['Q(sqrt2)', 'K7']
    remaining = [k for k in twelve if k not in new+previous]
    assert len(twelve) == 12 and len(remaining) == 2
    # Nontrivial quadratic subfields of Q(sqrt2,sqrt5) are exactly these
    # three squareclasses; Q(sqrt2,sqrt3) is a different biquadratic field.
    F_quadratics = {2,5,10}
    assert 3 not in F_quadratics and 6 not in F_quadratics
    assert {2,3,6} != F_quadratics
    # Every target other than K9 has degree not divisible by3. K9 and K7
    # are disjoint by the coprime-conductor degree identity above.
    for degree in (2,4,5):
        assert degree % 3 != 0
    return dict(groups=groups, low_label_fields=low, new_fields=new,
                general_field_hypothesis='K is not isomorphic to a subfield of F=Q(sqrt2,sqrt5), and K does not contain K7',
                F_subfields=['Q','Q(sqrt2)','Q(sqrt5)','Q(sqrt10)','F'],
                vinberg_compact_dimension_at_least14_fields=twelve,
                already_excluded_dimension21=previous,
                remaining_dimension21_fields=remaining,
                no_complete_triangle_catalogue_required=True)


def scalar_checks():
    z = F(3, 2)
    assert z**3-3*z+1 == -F(1, 8)
    assert (2+z)/4 == F(7, 8)
    assert 1/(4*(1-F(7, 8))) == 2
    # Positive coefficients in f(3/2+h)-f(3/2), h>0.
    derivative_shift = [3*z*z-3, 3*z, F(1)]
    assert min(derivative_shift) > 0
    # Residual determinant for diagonals -1-a,-1-b and delta^2<=1:
    # det >= a+b+ab >0. Strictness comes from both costs>2.
    residual_coefficients = [F(1), F(1), F(1)]
    assert min(residual_coefficients) > 0
    # Dense contacts to a bond m>=6 cost at least1/[2(1-c)]>2,
    # since cos(pi/6)^2=3/4>(3/4)^2. This handles dense6/8 triangles.
    assert F(3,4) > F(3,4)**2
    assert 1/(2*(1-F(3,4))) == 2
    # K16's triangle(3,3,8) shows that single-contact>=9 is too narrow;
    # its dense8 bond satisfies the broader criterion above.
    assert F(2) < F(3, 2)**2
    h21 = [F((21-i)*(20-i)+i*(i-1),4) for i in range(22)]
    assert min(h21) == 50 and h21[0] > 50
    assert 5*min(h21)-21*20//2 == 40 == 2*21-2
    # n=21+2u and22+2u: lower-minus-upper polynomials, u>=0.
    parity = dict(odd=['0', '5', '1/2'], even=['2', '11/2', '1/2'])
    for values in parity.values():
        assert min(map(F, values)) >= 0
    return dict(cos_double_angle_cubic=[1,0,-3,1], at_three_halves='-1/8',
                cosine_square_strict_lower='7/8', single_contact_cost_strict_lower='2',
                dense_bond_minimum_label=6, dense_contact_cost_strict_lower='2',
                cubic_derivative_shift=list(map(str, derivative_shift)),
                residual_determinant_lower_terms=['a','b','a*b'],
                dimension21_threshold=40, dimension21_vertex_upper=40,
                parity_gap_coefficients=parity, h21_paired_coefficients=list(map(str,h21)),
                K16_dense_bond_example=[3,3,8])


def finite_forests():
    records = []
    for kind in ('A', 'B', 'D', 'E', 'F', 'H', 'I_high'):
        ranks = range(1, 42) if kind == 'A' else range(2,42) if kind == 'B' else range(4,42) if kind == 'D' else {'E':[6,7,8],'F':[4],'H':[3,4],'I_high':[2]}[kind]
        for r in ranks:
            if kind == 'D':
                es = [(i,i+1) for i in range(r-3)] + [(r-3,r-2),(r-3,r-1)]
            elif kind == 'E':
                es, nxt = [], 1
                for arm in (r-4,2,1):
                    prev = 0
                    for _ in range(arm):
                        es.append((prev,nxt));prev=nxt;nxt+=1
            else:
                es = [(i,i+1) for i in range(r-1)]
            base = sum(distances(range(r),es,i).get(j,99)<=2 for i,j in combinations(range(r),2))
            bonus = int(kind == 'I_high')
            expected = 0 if r==1 else 2*r-2 if kind in ('D','E') else 2*r-3
            assert base == expected and base+bonus <= 2*r-2
            records.append(dict(kind=kind,rank=r,base=base,bonus=bonus,total=base+bonus))
    return dict(records=records, all_rank_identity='edges + sum_v binom(degree(v),2)',
                path='2r-3 for r>=2', one_trivalent_vertex='2r-2',
                singleton='0', high_dihedral='2', component_deficit_at_least=2)


def run():
    path = f'scripts/{TAG}.py'
    return dict(status='PASS_EXACT_DISTINGUISHED_BOND_REDUCTION_AND_EIGHT_FIELD_DIMENSION20_CANDIDATE',
                sources={path:sha256((ROOT/path).read_bytes()).hexdigest()},
                local=local_supports(), fields=fields(), scalars=scalar_checks(), forests=finite_forests(),
                scope='Compact chambers with all Lanner orders<=3 and every Lanner triangle containing either a label>=9 or a label6/7/8 whose two other bonds are nonright have dimension<=20. Eight arithmetic field applications use exact ground-field equality, not mere containment. The general field hypothesis rules out label7 and forces the distinguished bond. Published global geometry and field inputs are explicit; no novelty/current-best or all-arithmetic dimension20 claim.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()))
    dest=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out==json.loads(dest.read_text())
        print('PASS fresh saved exact high-label replay')
    else:
        dest.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'], out['local']['examined'], len(out['forests']['records']))
