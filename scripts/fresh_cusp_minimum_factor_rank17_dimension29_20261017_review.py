#!/usr/bin/env python3
"""Independent Fraction-only all-rank algebra; geometric proof read in MD."""
import argparse,ast,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_cusp_minimum_factor_rank17_dimension29_20261017'
TAG=PRIMARY+'_review'

# Sparse polynomials in independent formal variables k,j. No producer or CAS.
def clean(p):return {e:F(c)for e,c in p.items()if c}
def add(a,b):
    p=dict(a)
    for e,c in b.items():p[e]=p.get(e,F())+c
    return clean(p)
def scale(a,c):return clean({e:v*c for e,v in a.items()})
def mul(a,b):
    p={}
    for e,c in a.items():
        for f,d in b.items():
            h=(e[0]+f[0],e[1]+f[1]);p[h]=p.get(h,F())+c*d
    return clean(p)
def C(c):return clean({(0,0):F(c)})
def sub(a,b):return add(a,scale(b,-1))
def enc(p):return [[*e,str(c)]for e,c in sorted(p.items())]
K={(1,0):F(1)};J={(0,1):F(1)}

def arithmetic():
    one=C(1);two=C(2);k1=add(K,one);k2=add(K,two)
    records=[]
    for parity in(0,1):
        r=add(scale(K,2),C(parity));r1=add(r,one)
        for half in(0,1):
            y=add(J,C(F(half,2)))
            n0=sub(add(K,C(1 if parity==1 and half==0 else 0)),J)
            n1=sub(add(K,C(1 if parity==1 or half==0 else 0)),J)
            # Sum |y+l-r/2| = (r+1)y+N(r-2y-N+1).
            z0=mul(n0,add(sub(sub(r,scale(y,2)),n0),one))
            z1=mul(n1,add(sub(sub(r1,scale(y,2)),n1),one))
            numerator=sub(mul(z1,add(r,one)),mul(z0,add(r1,one)))
            if half==0:expected=sub(mul(k1,k1),mul(J,J))
            elif parity==0:expected=sub(mul(K,k1),mul(J,add(J,one)))
            else:expected=sub(mul(k1,k2),mul(J,add(J,one)))
            assert numerator==expected,(parity,half,enc(numerator),enc(expected))
            records.append(dict(rank_parity=parity,argument_half=half,cleared_numerator=enc(numerator)))
    # Independent literal check of the identities, including zero tails.
    literals=[]
    for r in range(17,35):
        vals=[]
        def mean_abs(rank,y):return sum((abs(y+F(i)-F(rank,2))for i in range(rank+1)),F())/(rank+1)
        for q in range(0,r+5):
            y=F(q,2);delta=mean_abs(r+1,y)-mean_abs(r,y);assert delta>=0;vals.append(delta)
        maximum=F(r+2,4*(r+1))if r%2==0 else F(r+3,4*(r+2))
        assert max(vals)==maximum<=F(5,19)
        literals.append(dict(rank=r,breakpoints=len(vals),maximum=str(maximum)))
    # For r=18+2t or17+2t, 5/19 minus either maximum is
    # t/[38(19+2t)], so it is nonnegative for every t>=0.
    r_even=add(scale(K,2),C(18));r_odd=add(scale(K,2),C(17))
    cap_num_even=sub(scale(add(r_even,one),20),scale(add(r_even,two),19))
    cap_num_odd=sub(scale(add(r_odd,two),20),scale(add(r_odd,C(3)),19))
    assert cap_num_even==cap_num_odd==scale(K,2)
    increment=F(3,4)-F(5,2)*F(5,19);assert increment==F(7,76)>0
    # Single-factor / finite-vertex coefficient inequalities.
    finite=[]
    for parity,base in[(0,30),(1,31)]:
        n=add(scale(K,2),C(base))
        num=sub(mul(n,sub(n,C(6)))if parity==0 else mul(sub(n,one),sub(n,C(5))),scale(n,24))
        expected={(2,0):F(4),(1,0):F(60 if parity==0 else 64)}
        if parity:expected[(0,0)]=F(36)
        assert num==expected and all(c>=0 for c in num.values())
        finite.append(dict(parity=parity,eight_times_beta_minus_3n=enc(num)))
    # For s>=2, subtracting3s leaves39s/4 above the variance term.
    squared=sub(scale(mul(K,K),F(39,4)**2*48),scale(K,F(25,4)*F(323,12)*48))
    assert squared=={(2,0):F(4563),(1,0):F(-8075)}
    assert 4563*2-8075==1051>0
    assert F(17*19,12)==F(323,12)
    return dict(pointwise_breakpoint_identities=records,literal_breakpoint_checks=literals,
        pointwise_cap='5/19',cap_slack='t/(38*(19+2*t))',curvature_increment=str(increment),
        finite_and_single_factor_parities=finite,multiple_factor_squared_gap_times48=enc(squared),
        shifted_positive_linear_factor=[1051,4563],global_additive_slack='5*(n-1)>0 for n>=30',
        all_rank_scope='Polynomial identities and nonnegative coefficients establish the unbounded scalar claims; the18 literal ranks are corroboration only.')

def run():
    algebra=arithmetic()
    producer=ROOT/'results'/f'{PRIMARY}.json';p=json.loads(producer.read_text())
    assert p['status']=='PASS_ALL_RANK_MINIMUM_CUSP_FACTOR_RANK17_DIMENSION29_COMPOSITION'
    assert p['dimension_bound']==29 and p['affine_order_exclusion_through']==17
    assert p['arithmetic_hypothesis']is False and p['factor_count_limit']is None
    assert p['enumeration_cutoff_used']is False and p['unrestricted163_tail_used']is False
    assert F(p['curvature_increment'])==F(7,76)
    for name,h in p['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    paths=[Path(__file__),ROOT/'scripts'/f'{PRIMARY}.py',producer,ROOT/'results'/f'{PRIMARY}.md',
        ROOT/'results/fresh_affine_cycle_carrier_exclusion_20261017.md',
        ROOT/'scripts/fresh_affine_tree_carrier_exclusion_20261017.py',
        ROOT/'results/fresh_affine_tree_carrier_exclusion_20261017.json',
        ROOT/'results/fresh_affine_tree_carrier_exclusion_20261017.md']
    for path in paths:
        if path.suffix=='.py':ast.parse(path.read_text())
    return dict(status='PASS_INDEPENDENT_RATIONAL_SCALARS_AND_MATHEMATICAL_COMPOSITION_REVIEW',
        sources={str(q.relative_to(ROOT)):sha256(q.read_bytes()).hexdigest()for q in paths},
        independent_algebra=algebra,producer_direct_bindings_checked=len(p['sources']),
        conclusion='dimension<=29 under the minimum affine cusp factor rank17 hypothesis',
        scope='Independent Fraction-only scalar derivation, producer metadata/direct hashes, and separately written hostile mathematical/source review. No producer arithmetic imports, no SymPy, no duplicate large graph enumerations or published finite-diagram tables, no formal verification or external refereeing.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--arithmetic-only',action='store_true');a=ap.parse_args()
    if a.arithmetic_only:print(json.dumps(arithmetic(),indent=2))
    else:
        d=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
        if a.verify:assert d==json.loads(path.read_text())
        else:path.write_text(json.dumps(d,indent=2)+'\n')
        print(d['status'],flush=True)
