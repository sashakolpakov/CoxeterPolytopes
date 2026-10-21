#!/usr/bin/env python3
"""Independent all-label exterior-polar classification.

Uses certified old ray coordinates, independent nested-radical geometry,
and a different arctangent identity and alternating cosine enclosures.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import fresh_nonarithmetic_joint_regions_20261016_batch_review as B

ROOT,I=B.ROOT,B.I
TAG='fresh_nonarithmetic_joint_regions_20261016_polar_gate_review'
PRIMARY=B.POLAR
BATCH_REVIEW='results/fresh_nonarithmetic_joint_regions_20261016_batch_review.json'


def atan_recip(n,terms=240):
    total=F(0)
    for k in range(terms):total+=F((-1)**k,(2*k+1)*n**(2*k+1))
    remainder=F(1,(2*terms+1)*n**(2*terms+1))
    return I(total,total+remainder) if terms%2==0 else I(total-remainder,total)


# tan(atan(1/2)+atan(1/3))=1, with the sum in(0,pi/2).
PI=4*(atan_recip(2)+atan_recip(3))


def rational_cos(x,terms=90):
    assert 0<x<2
    value=F(1);term=F(1)
    for k in range(1,terms+1):
        term*=x*x/F((2*k-1)*(2*k))
        value+=(-1)**k*term
    remainder=term*x*x/F((2*terms+1)*(2*terms+2))
    return (value-remainder,value) if terms%2==0 else (value,value+remainder)


@lru_cache(None)
def threshold(m):
    assert m>=2
    if m==2:return I()
    # cos decreases on[0,pi/2]; enclose each rational endpoint separately.
    lo=rational_cos(PI.hi/m)[0]
    hi=rational_cos(PI.lo/m)[1]
    assert 0<lo<hi<1
    return I.out(4*lo*lo,4*hi*hi)


def run():
    target=B.read(PRIMARY);batchreview=B.read(BATCH_REVIEW);G=B.geometry()
    assert batchreview['status']=='PASS_INDEPENDENT_SIXTEEN_POLARS_AND_SIMULTANEOUS_48_FACETS'
    assert batchreview['sources'][PRIMARY]==B.digest(PRIMARY)
    hyper=G['inherited']['complete_census']['hyperideal_rays']
    assert len(hyper)==1068
    candidates=G['batch']['candidates']+target['additional_exact_polars']
    survivors=[]
    for row in candidates:
        index=row['ray_index'];active=row.get('right_support',row.get('active'))
        assert active==hyper[index]
        survivors.append(index)
    assert len(survivors)==len(set(survivors))==16
    assert sorted(survivors)==target['compatible_ray_indices']
    norms=[B.iv(B.pair(n,n,G['Q'])) for n in G['N'][:36]]
    assert all(q.lo>0 for q in norms)
    records=[];seen=set()
    for row in target['gap_records']:
        index=row['ray_index'];wall=row['wall'];lo,hi=row['neighbor_labels']
        assert index not in seen and index not in survivors and hi==lo+1 and lo>=2
        seen.add(index);assert row['active']==hyper[index]
        inherited=G['old_rays'][tuple(hyper[index])]
        v=[I(*map(F,x)) for x in inherited['normalized_coordinate_bounds']]
        q=B.idot(v,B.imv(G['QI'],v));g=B.idot(v,G['covI'][wall])
        assert q.lo>0 and g.hi<0
        cycle=4*g*g/(q*norms[wall])
        lower,upper=threshold(lo),threshold(hi)
        assert lower.hi<cycle.lo<=cycle.hi<upper.lo,(index,wall,lo)
        records.append(dict(ray_index=index,wall=wall,neighbor_labels=[lo,hi],
                            independent_cycle_bounds=cycle.saved(),
                            positive_lower_gap=str(cycle.lo-lower.hi),
                            positive_upper_gap=str(upper.lo-cycle.hi)))
        if len(records)%200==0:print('independent all-label polar gaps',len(records),flush=True)
    assert len(records)==1052 and set(survivors)|seen==set(range(1068))
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',PRIMARY,BATCH_REVIEW,B.BATCH,
           B.CONE_REVIEW,B.OLD,B.OLD_REVIEW]
    return dict(status='PASS_INDEPENDENT_COMPLETE_ALL_LABEL_1068_POLAR_CLASSIFICATION',
                sources={p:B.digest(p) for p in paths},recursive_integrity=B.R.integrity(paths),
                complete_polar_count=1068,compatible_polar_count=16,rejected_polar_count=1052,
                compatible_ray_indices=sorted(survivors),gap_records=records,
                independent_pi_interval=PI.saved(),
                independent_threshold_intervals={str(m):threshold(m).saved() for m in sorted({v for r in records for v in r['neighbor_labels']})},
                proof='Pi=4*(atan(1/2)+atan(1/3)), with240-term rational alternating-series brackets.90-term rational cosine alternating brackets at the opposite endpoints enclose each4cos²(pi/m). Independent normalized ray coordinate intervals and independently decoded normals give each squared product. Every rejection lies strictly between adjacent monotone label thresholds; all16accepted rays were exact-contact checked in the independent batch review.',
                scope='Complete one-generation polar screen for the fixed independently complete36wall cone. Exactly16Coxeter-compatible polars, all sharing the independently verified strict spacelike escape. No conclusion about nonpolar normals or polars of subsequently clipped cones; no finite-volume example.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==B.read(str(p.relative_to(ROOT)))
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
