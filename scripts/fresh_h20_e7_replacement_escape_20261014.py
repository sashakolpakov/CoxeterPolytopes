#!/usr/bin/env python3
"""Exact common escape for the 42 saved mixed-contact cap proposals.

The proposed walls need not be mutually Coxeter. This is an obstruction
even after allowing all of them at once. No numerical acceptance is used.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_replacement_escape_20261014'
BASE = 'results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
PARENT = 'results/fresh_h20_e7_intrinsic_caps_20261014.json'
PROPOSALS = 'results/fresh_h20_e7_causal_residuals_20261014_interior_scout.json'
Z = K.Z


def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def ext(x): return (x, Z)
def add(x, y): return (K.add(x[0], y[0]), K.add(x[1], y[1]))
def neg(x): return (K.neg(x[0]), K.neg(x[1]))
def sub(x, y): return add(x, neg(y))
def mul(x, y):
    return (K.add(K.multiply(x[0], y[0]), K.scale(K.multiply(x[1], y[1]), 3)),
            K.add(K.multiply(x[0], y[1]), K.multiply(x[1], y[0])))
def scale(x, q): return (K.scale(x[0], q), K.scale(x[1], q))
def sign(x):
    a, b = x; sa, sb = K.sign(a), K.sign(b)
    if not sb: return sa
    if not sa or sa == sb: return sb
    v = K.sub(K.multiply(a, a), K.scale(K.multiply(b, b), 3))
    return sa*K.sign(v)
def radical_sign(a, b, q):
    """Sign(a+b sqrt(q)), with q positive in K(sqrt3)."""
    assert sign(q) > 0
    sa, sb = sign(a), sign(b)
    if not sb: return sa
    if not sa or sa == sb: return sb
    return sa*sign(sub(mul(a, a), mul(mul(b, b), q)))
def enc(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, (list, tuple)): return [enc(y) for y in x]
    if isinstance(x, dict): return {str(k): enc(v) for k, v in x.items()}
    return x


def run():
    base, parent, proposals = [json.loads((ROOT/p).read_text()) for p in [BASE, PARENT, PROPOSALS]]
    for doc in [base, parent, proposals]:
        for p, h in doc['sources'].items(): assert sha(p) == h, p
    assert base['normals'][:85] == parent['normals']
    Q = [[F(x) for x in row] for row in base['form']]
    sparse = [[(j, a) for j, a in enumerate(row) if a] for row in Q]
    cache = {}
    def dot(x, y):
        if y not in cache:
            cache[y] = tuple(tuple(sum(a*y[j][k] for j, a in row) for k in range(8)) for row in sparse)
        v = Z
        for a, b in zip(x, cache[y]):
            if any(a) and any(b): v = K.add(v, K.multiply(a, b))
        return v
    def linear(vectors, coeff):
        v = (Z,)*21
        for w, c in zip(vectors, coeff): v = K.vadd(v, K.vscale(w, c))
        return v
    N = list(map(K.vector, base['normals']))
    H = K.vector(base['strict_center'])
    oldE = K.vector(base['strict_escape']['point'])
    S = base['E7_support']; A = parent['face_facet_walls'][:13]
    I = K.inverse([[K.rational(dot(N[i], N[j])) for j in S] for i in S])
    J = K.inverse([[K.rational(dot(N[i], N[j])) for j in A] for i in A])
    W = [linear([N[i] for i in S], list(map(K.scalar, row))) for row in I]
    L = [linear([N[i] for i in A], list(map(K.scalar, row))) for row in J]
    sq2 = (F(0), F(1))+(F(0),)*6
    T = K.vscale(K.vector(parent['face_basis'][-1]), K.scale(sq2, F(1, 2)))
    assert dot(T, T) == K.scalar(-2)
    assert all(dot(T, w) == Z for w in W+L)
    assert all(dot(w, l) == Z for w in W for l in L)
    # Rationalized numerical proposal; its validity is checked from scratch.
    numerators = [1856660247,768568635,1211012212,1047964683,884927155,
                  721899626,558882098,395874569,232877041,139685024,
                  46503008,-46669008,-139831024,-232983041,-326125057,
                  -419257073,-512379090,-605491106,605501106,186248033,-93102016]
    E = tuple(K.scalar(F(x, 10**9)) for x in numerators)
    qE, EH = dot(E, E), dot(E, H)
    assert K.sign(qE) > 0 and K.sign(EH) < 0 and K.sign(dot(H, H)) < 0
    retained = [i for i in range(87) if i not in [61, 75]]
    old_pairs = [dot(E, N[i]) for i in retained]
    assert all(K.sign(p) < 0 for p in old_pairs)
    assert all(K.sign(dot(H, N[i])) < 0 for i in retained)
    values = [ext(K.scalar(1)), ext(sq2), (Z, K.scalar(1)), ext(K.scalar(2)), ext(K.scalar(3)), ext(K.scalar(4))]
    records = []
    candidates = proposals['delete61_75_candidates']; assert len(candidates) == 42
    for index, row in enumerate(candidates):
        support = row['E7_support']; assert row['E7_labels'] == [3]*len(support)
        assert len(support) == 3 and len(set(support)) == 3
        indices, choices = row['source_recipe']
        d = [ext(Z) for _ in range(13)]
        for i, j in zip(indices, choices): d[i] = values[j]
        mass = sum(I[i][j] for i in support for j in support)-2
        spaceq = ext(Z)
        for i in indices:
            for j in indices: spaceq = add(spaceq, scale(mul(d[i], d[j]), J[i][j]))
        t2 = scale(add(spaceq, ext(K.scalar(mass))), F(1, 2))
        assert sign(t2) > 0
        def pairing(v):
            a = ext(K.neg(dot(linear([W[i] for i in support], [K.scalar(1)]*len(support)), v)))
            for i in indices: a = sub(a, mul(d[i], ext(dot(L[i], v))))
            b = ext(dot(T, v))
            sg = radical_sign(a, b, t2)
            return dict(constant=a,radical_coefficient=b,sign=sg)
        ep, hp = pairing(E), pairing(H)
        assert ep['sign'] < 0 and hp['sign'] < 0
        records.append(dict(index=index,E7_support=support,source_indices=indices,source_value_indices=choices,
                            time_coefficient_squared=t2,escape_pairing=ep,center_pairing=hp))
    # Independently rebuild the particularly simple cap selected for exact construction.
    cap = linear([W[2],W[4],W[6],L[0],L[12],T],
                 [K.scalar(-1)]*3+[K.neg(sq2),K.scalar(-2),K.scalar(5)])
    assert dot(cap, cap) == K.scalar(2)
    assert K.sign(dot(cap, oldE)) > 0 and K.sign(dot(cap, E)) < 0
    assert candidates[34]['E7_support'] == [2,4,6] and candidates[34]['source_recipe'] == [[0,12],[1,3]]
    paths = [BASE,PARENT,PROPOSALS,f'scripts/{TAG}.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    return enc(dict(status='PASS_EXACT_COMMON_ESCAPE_FOR_ALL_42_MIXED_CONTACT_PROPOSALS',
        sources={p:sha(p) for p in paths},retained_base_walls=retained,removed_base_walls=[61,75],
        point=E,norm=qE,center_pairing=EH,retained_wall_pairings=old_pairs,records=records,
        selected_cap=dict(proposal_index=34,normal=cap,norm=dot(cap,cap),old_escape_pairing=dot(cap,oldE),
                          new_escape_pairing=dot(cap,E)),
        scope='One strict spacelike point survives all85 retained base walls and all42 listed cap proposals simultaneously. The inherited strict timelike point also survives all127 inequalities. Hence every subset of these proposals with this core has infinite volume. Mutual Coxeter compatibility, facet actuality, and unsaved continuous families are not asserted.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify',action='store_true'); args = ap.parse_args()
    result = run(); path = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['norm'])
