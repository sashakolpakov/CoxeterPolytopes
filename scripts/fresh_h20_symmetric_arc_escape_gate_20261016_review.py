#!/usr/bin/env python3
"""Independent symbolic and rational review of symmetric-arc common escape."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_symmetric_arc_escape_gate_20261016_review'
PRIMARY='results/fresh_h20_symmetric_arc_escape_gate_20261016.json'
SINGLE='results/fresh_h20_dotted_arc_family_gate_20261016.json'
CHAIN='results/fresh_h20_s61_multicap_20261014_exact.json'
CAP='results/fresh_h20_dotted_arc_cap_20261016.json'
CAP_REVIEW='results/fresh_h20_dotted_arc_cap_20261016_review.json'
ABEL_REVIEW='results/fresh_h20_chain_abelian_gate_20261016_review.json'


def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for name,h in d.get('sources',{}).items():assert sha(name)==h,name
    return d
def rat(q):return S.Rational(str(q))
def expression(v):return rat(v[0])+rat(v[1])*S.sqrt(2)
def positive(q):return S.simplify(q).is_positive is True


def run():
    primary,single,chain,cap,cr,ar=map(read,(PRIMARY,SINGLE,CHAIN,CAP,CAP_REVIEW,ABEL_REVIEW))
    assert cr['sources'][CAP]==sha(CAP)
    assert ar['minimum_dotted_chain_contacts']==5
    support=chain['support']
    G=S.Matrix([[rat(x) for x in row] for row in chain['support_gram']])
    assert G.det()>0
    J=G.inv()
    assert J==S.Matrix([[rat(x) for x in row] for row in chain['support_inverse']])
    assert all(q>=0 for q in J)
    inside=[support.index(i) for i in range(1,11)]
    outside=sorted(set(range(18))-set(inside))
    assert J.extract(inside,inside)==S.Matrix([[S.Rational(min(i,j)*(11-max(i,j)),11) for j in range(1,11)] for i in range(1,11)])
    assert all(J[i,j]==0 for i in inside for j in outside)
    assert min(J[i,i] for i in range(18))==S.Rational(1,2)
    assert sum(J[i,j] for i in outside for j in outside)==13
    L,A,B,c=S.symbols('L A B c',positive=True)
    t=((B-A)*L+2*B)/(2*S.sqrt(L))
    ysquare=S.expand((B-A)**2*(L+2)/4-t*t)
    assert S.simplify(ysquare-(A*A/2-B*B*(L+2)/(2*L)))==0
    triple=S.Matrix([[2,0,-4*c-2],[0,2,0],[-4*c-2,0,2]])
    assert S.factor(triple.det())==-32*c*(c+1)
    lo,hi=map(rat,chain['parameter']['isolating_interval'])
    assert lo>0 and lo*lo>280
    total,small=S.Rational(281,70),S.Rational(141,70)
    assert 4+4/(lo*lo)<total and total-2==small
    assert total-S.Rational(40,11)<S.Rational(1,2)
    assert 4*S.Rational(18,11)>total and S.Rational(24,11)>small
    # Label>=5 has squared magnitude >=phi²>5/2.
    assert S.Rational(10,11)*S.Rational(5,2)>small
    assert S.Rational(18,11)*2>small
    ef=list(map(rat,cap['strict_escape']['frame_coefficients']))
    eb=S.Matrix(list(map(rat,cap['strict_escape']['dual_coefficients'])))
    assert ef==[1,-S.Rational(1,100),0]
    assert list(eb)==[-1 if i==5 else -S.Rational(1,100000) for i in support]
    z=-J*eb
    K=ef[0]**2-ef[1]**2-ef[2]**2
    assert K==S.Rational(9999,10000) and all(v>0 for v in z)
    tiny=S.Rational(1,1000)
    assert S.Rational(1,100000)**2*13*total<tiny*tiny
    # Independently check the broader theorem's two continuous endpoint cases,
    # and all eight continuous cases of the earlier single-support theorem.
    single_cases=[]
    ceiling=(1+(lo*lo+2)/(lo*lo))/2
    allowed=[j for j in range(18) if J[j,j]<ceiling]
    assert [support[j] for j in allowed]==[1,10,12,15,17,31,61,87]
    for j in allowed:
        u=max(4*J[j,j]-2,1/(hi*hi))
        lead=J[j,j]*K-z[j]*z[j]
        margin=u*lead-2*z[j]*z[j]
        assert lead>0 and margin>0
        row=next(r for r in single['remaining_exact_margins'] if r['wall']==support[j])
        assert rat(row['CS_strict_margin'])==margin
        single_cases.append(dict(wall=support[j],strict_margin=str(margin)))
    for row in primary['dotted_A10_endpoint_cases']:
        j=support.index(row['wall'])
        u=4*J[j,j]-2
        margin=u*(J[j,j]*K-z[j]*z[j])-2*z[j]*z[j]
        assert margin>0 and margin==rat(row['strict_margin'])
    # Enumerate in the opposite slot order with a symbolic radical algebra,
    # rather than the producer's two-component arithmetic.
    slots=[10,9,2,1]
    choices=[(0,1,S.sqrt(2)),(0,1),(0,1),(0,1,S.sqrt(2))]
    kept=[];excluded=0
    original={tuple(row['labels']):row for row in primary['A10_profiles']}
    rejected={tuple(row['labels']):row for row in primary['A10_excluded_profiles']}
    for vals in product(*choices):
        d=S.zeros(18,1)
        for wall,value in zip(slots,vals):d[support.index(wall)]=value
        cost=S.expand((d.T*J*d)[0])
        labels=tuple(reversed([0 if x==0 else 3 if x==1 else 4 for x in vals]))
        if positive(cost-small) or cost==small:
            excluded+=1
            assert S.simplify(cost-expression(rejected[labels]['cost']))==0
            continue
        assert positive(small-cost) and cost.is_Rational
        value=S.expand((d.T*z)[0])
        u=max(cost,1/(hi*hi))
        margin=S.expand(K*u-(value+tiny)**2)
        assert positive(margin)
        row=original[labels]
        assert cost==rat(row['A10_cost']) and S.simplify(value-expression(row['escape_A10_pair_magnitude']))==0
        assert S.simplify(margin-expression(row['strict_squared_margin']))==0
        kept.append(dict(labels=list(labels),cost=str(cost),margin=str(margin)))
    assert len(kept)==8 and excluded==28
    assert {tuple(row['labels']) for row in kept}==set(original)
    paths=[f'scripts/{TAG}.py',PRIMARY,SINGLE,CHAIN,CAP,CAP_REVIEW,ABEL_REVIEW]
    return dict(status='PASS_INDEPENDENT_ALL_SUPPORT_SYMMETRIC_ARC_ESCAPE_OBSTRUCTION',
                sources={p:sha(p) for p in paths},single_support_cases=single_cases,
                all_support_profiles=sorted(kept,key=lambda r:r['labels']),excluded_finite_profiles=excluded,
                all_label_support_budget=str(total),outside_inverse_sum='13',
                scope='Independent rational Gram inversion, symbolic triple/norm identities, all-label budget audit, eight single-support margins and complete36-profile symbolic census with eight retained margins. Inherits reviewed nonabelian chain theorem and independently verified strict escape. No finite-volume or general H20 completion exclusion; asymmetric or sparser nondotted contacts remain outside.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
