#!/usr/bin/env python3
"""Tiny exact coefficients and equal-sign value for the sharp cusp bound.

The all-moment balancing/global-maximizer proof is written separately.
"""
import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import product
import json
from math import ceil,comb
from pathlib import Path
import resource
import signal
import sys
from run_lean_batches_20261018 import resources

ROOT=Path(__file__).resolve().parents[1]
TAG='ideal_ra6_sharp_reciprocal_20261018'


def encode(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return{k:encode(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return[encode(v)for v in x]
    return x


def run():
    coefficients=[]
    for j in range(17):
        row=[]
        for k in range(j+1):
            lhs=comb(2*j,2*k)
            rhs=sum(comb(j,2*l)*4**l*comb(j-2*l,k-l)
                for l in range(j//2+1)if l<=k<=j-l)
            assert lhs==rhs;row.append(lhs)
        coefficients.append(dict(j=j,monomial_coefficients=row))
    mu=Q(3,2);d2={}
    for signs in product((-1,1),repeat=5):
        square=Q(sum(signs)**2,20);d2[square]=d2.get(square,0)+1
    assert d2=={Q(5,4):2,Q(9,20):10,Q(1,20):20}
    terms=[dict(square=x,multiplicity=n,contribution=mu*Q(n,32)/(mu*mu-x))for x,n in sorted(d2.items())]
    rho=sum(t['contribution']for t in terms)
    assert rho==Q(103,132)<Q(25,32)<Q(71,90)
    degree=32*rho;corners=64+degree/2;residual=(corners-60)/12;ratio=27/residual
    assert degree==Q(824,33)and corners==Q(2524,33)>72
    assert residual==Q(136,99)and ratio==Q(2673,136)
    assert ceil(227*ratio)==4462 and 227*2673==606771
    paths=[Path(__file__),ROOT/'scripts/fresh_ideal_ra6_cusp_transport_20261018.py',
        ROOT/'results/ideal_ra6_cusp_transport_research_20261018.json',
        ROOT/'results/ideal_ra6_cusp_transport_independent_review_20261018.md',
        ROOT/'results/right_angled_ideal6_carrier_propagation_20261018.md']
    hashes={}
    for p in paths:
        with p.open('rb')as stream:hashes[str(p.relative_to(ROOT))]=file_digest(stream,'sha256').hexdigest()
    return encode(dict(status='PASS_IDEAL_RA6_SHARP_RECIPROCAL',sources=hashes,
        checked_pair_identity_degrees=list(range(17)),pair_identity_coefficients=coefficients,
        equal_square_grouped_terms=terms,sharp_mean_reciprocal=rho,
        cubic_good_degree_upper=degree,total_triangular_corner_average_upper=corners,
        ridge_minus_facet_to_vertex_upper=residual,vertex_to_facet_lower=ratio,
        inherited_facet_lower=227,new_vertex_lower=4462,
        scope=dict(all_moment_balancing_proved_in_written_report=True,
            finite_coefficient_checks_not_universal_proof=True,
            original_cusp_geometry_independently_reviewed=True,
            lean_verified=False,new_dimension_ceiling=False,finite_volume_example=False)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();signal.alarm(59)
    with(ROOT/'formal/.proof-batches/exclusive.lock').open('a')as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);resources();out=run()
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin'else 1024);assert rss<1024**3
        path=ROOT/f'results/{TAG}.json'
        if args.verify:
            with path.open()as stream:assert out==json.load(stream)
        else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
        resources(start=False)
    print(out['status'],out['sharp_mean_reciprocal'],out['new_vertex_lower'],rss,flush=True)
