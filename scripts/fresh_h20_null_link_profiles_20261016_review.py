#!/usr/bin/env python3
"""Independent exact tower-arithmetic census of the through-cusp domain.

Gram inversion is direct; comparisons use recursive ordered quadratic
towers, not the producer's square-root interval comparisons.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_null_link_profiles_20261016_review'
PRIMARY = 'results/fresh_h20_null_link_profiles_20261016.json'
GEOMETRY = 'results/fresh_h20_affine_a17_completion_20261016_review.json'
PATHS = [[3,4], [6,7,8,9,10], [12,13,14,15], [17], [31]]
RADICALS = (1,2,3,6,5,10,15,30)
TOWER = (2,3,5)


def add(a, b):
    return tuple(x+y for x,y in zip(a,b))


def scale(a, q):
    return tuple(q*x for x in a)


@lru_cache(None)
def mul(a, b):
    if len(a) == 1:
        return (a[0]*b[0],)
    half = len(a)//2
    d = TOWER[half.bit_length()-1]
    aa, ab, ba, bb = a[:half], a[half:], b[:half], b[half:]
    return add(mul(aa,ba),scale(mul(ab,bb),d))+add(mul(aa,bb),mul(ab,ba))


@lru_cache(None)
def sign(a):
    if len(a) == 1:
        return (a[0]>0)-(a[0]<0)
    half = len(a)//2
    x,y = a[:half],a[half:]
    sx,sy_ = sign(x),sign(y)
    if not sx:
        return sy_
    if not sy_ or sx == sy_:
        return sx
    d = TOWER[half.bit_length()-1]
    cmp = sign(add(mul(x,x),scale(mul(y,y),-d)))
    assert cmp, 'Nonzero opposite-sign tower halves cannot cancel'
    return sx*cmp


def run():
    primary = json.loads((ROOT/PRIMARY).read_text())
    geometry = json.loads((ROOT/GEOMETRY).read_text())
    assert geometry['status']=='PASS_INDEPENDENT_140_WALL_G2_A17_CUSP_COMPLETION'
    white = sum(PATHS,[])
    assert geometry['old_null_incidence']==sorted([0,1,134]+white)
    G = sy.zeros(13)
    offset = 0
    for path in PATHS:
        for i in range(len(path)):
            G[offset+i,offset+i]=2
            if i:
                G[offset+i,offset+i-1]=G[offset+i-1,offset+i]=-1
        offset += len(path)
    J = [[F(x) for x in row] for row in G.inv().tolist()]
    assert all(x>=0 for row in J for x in row)
    assert G.is_positive_definite
    assert primary['white']==white and primary['components']==PATHS
    assert primary['inverse_Gram']==[[str(x) for x in row] for row in J]
    diagonal = sorted(J[i][i] for i in range(13))
    assert sum(diagonal[:4])==F(7,3)>2
    assert [white[i] for i in range(13) if J[i][i]<F(2,3)]==[17,31]
    assert min(diagonal)==F(1,2)

    zero = (F(0),)*8
    def element(terms):
        return tuple(F(terms.get(i,0)) for i in RADICALS)
    magnitudes = {3:element({1:1}),4:element({2:1}),
                  5:element({1:F(1,2),5:F(1,2)}),6:element({3:1})}
    pair_values = {(m,n):mul(magnitudes[m],magnitudes[n]) for m in magnitudes for n in magnitudes}
    two = element({1:2})
    accepted = {((),()):zero}
    trials = 0
    rejected = 0
    equal_budget = []
    for count in range(1,4):
        for indices in combinations(range(13),count):
            for labels in product((3,4,5,6),repeat=count):
                trials += 1
                cost = zero
                for a,i in enumerate(indices):
                    cost = add(cost,scale(pair_values[labels[a],labels[a]],J[i][i]))
                    for b in range(a):
                        cost = add(cost,scale(pair_values[labels[a],labels[b]],2*J[i][indices[b]]))
                comparison = sign(add(cost,scale(two,-1)))
                if comparison>0:
                    rejected += 1
                    continue
                key = tuple(white[i] for i in indices),labels
                accepted[key]=cost
                if not comparison:
                    equal_budget.append(dict(walls=list(key[0]),labels=list(labels)))
    recorded = {}
    for row in primary['profiles']:
        key=tuple(row['walls']),tuple(row['labels'])
        assert key not in recorded
        by_radical = dict(zip(primary['squarefree_basis'],map(F,row['cost'])))
        recorded[key] = tuple(by_radical[d] for d in RADICALS)
    assert accepted==recorded
    assert len(accepted)==primary['finite_profiles']==125
    assert trials==primary['finite_trials']==19604
    assert rejected==primary['rejected_over_budget']
    assert primary['continuous_tails']==[
        dict(wall=i,other_incident_contacts='right',magnitude_interval=['sqrt3','2'],
             overlaps_finite_label=6,includes_parallel=True) for i in (17,31)]
    # Exact sufficient inequalities for the completeness reduction:
    # d>2 costs >2 everywhere; four nonzero contacts cost >=7/3;
    # nonisolated d^2<=3; a second nonzero contact makes isolated d^2<=3.
    assert F(2,3)*3==2
    assert F(1,2)*3+F(1,2)==2
    assert F(1,2)*4==2
    paths=[f'scripts/{TAG}.py',PRIMARY,GEOMETRY]
    return dict(status='PASS_INDEPENDENT_COMPLETE_125_PLUS_TWO_INCIDENT_DOMAIN',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        independent_Gram=[[str(x) for x in row] for row in G.tolist()],
        independent_inverse=[[str(x) for x in row] for row in J],
        finite_trials=trials,finite_profiles=len(accepted),rejected_over_budget=rejected,
        support_size_counts=dict(sorted(Counter(len(k[0]) for k in accepted).items())),
        equal_budget_profiles=equal_budget,continuous_tail_roots=[17,31],
        comparison_method='Recursive exact signs in Q(sqrt2)(sqrt3)(sqrt5); no floating point or interval sign decisions',
        scope='Complete necessary all-label contact domain for norm2 acute walls through the specified nu; 125 finite profiles plus two continuous relaxations. No escape exclusion or wall existence claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args=ap.parse_args()
    result=run()
    # JSON object keys are strings in a saved certificate.
    result=json.loads(json.dumps(result))
    out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['finite_profiles'],result['finite_trials'],flush=True)
