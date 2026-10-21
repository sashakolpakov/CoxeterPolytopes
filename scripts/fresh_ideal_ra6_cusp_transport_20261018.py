#!/usr/bin/env python3
"""Small exact cusp-frame/moment checks; universal geometry is written separately."""
import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import combinations, product
import json
from math import ceil, comb
from pathlib import Path
import resource
import signal
import sys
from run_lean_batches_20261018 import resources

ROOT = Path(__file__).resolve().parents[1]
TAG = 'ideal_ra6_cusp_transport_research_20261018'


def dot(a, b):
    return -a[0]*b[1]-a[1]*b[0]+sum(x*y for x,y in zip(a[2:],b[2:]))


def plus(a, b):
    return tuple(x+y for x,y in zip(a,b))


def times(t, v):
    return tuple(t*x for x in v)


def frame(widths, signs):
    m = len(widths)
    widths = list(map(Q, widths))
    total = sum(w*w for w in widths)
    v = (Q(1), Q(0))+(Q(0),)*m
    z = (Q(0), Q(1))+(Q(0),)*m
    normals = [(widths[i]/2, Q(0))+tuple(Q(signs[i]*(i==j)) for j in range(m))
               for i in range(m)]
    thirds = [((total-2*widths[i]**2)/(4*widths[i]), 2/widths[i])+
              tuple(Q(0) if i==j else signs[j]*widths[j]/widths[i] for j in range(m))
              for i in range(m)]
    endpoint = (total/4, Q(2))+tuple(signs[i]*widths[i] for i in range(m))
    assert dot(endpoint,endpoint)==0 and dot(endpoint,v)==-2
    assert all(dot(n,n)==1 for n in normals+thirds)
    assert all(dot(thirds[i],thirds[j])==Q(i==j) for i in range(m) for j in range(m))
    for i in range(m):
        assert plus(normals[i],thirds[i])==times(1/widths[i],endpoint)
        assert all(dot(normals[j],thirds[i])==(-1 if i==j else 0) for j in range(m))
    if all(w==1 for w in widths):
        h = plus(z,times(Q(4-m,8),v))
        assert dot(h,h)==Q(m-4,4) and dot(h,v)==dot(h,endpoint)==-1
        new_z = plus(h,times(Q(m-4,8),endpoint))
        basis = [plus(n,times(Q(-1,2),endpoint)) for n in normals]
        assert dot(new_z,new_z)==0 and dot(new_z,endpoint)==-1
        assert all(dot(new_z,e)==0 for e in basis)
        assert all(dot(basis[i],basis[j])==Q(i==j) for i in range(m) for j in range(m))
    return thirds


def encode(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def run():
    frames = 0
    for widths in ((1,1,1),(1,1,1,1),(1,1,1,1,1),(1,2,3,4,5)):
        for signs in product((-1,1),repeat=len(widths)):
            frame(widths,signs)
            frames += 1
    # Rademacher square/fourth moments as exact multivariate polynomials.
    square, fourth = {}, {}
    for signs in product((-1,1),repeat=5):
        for degree, target in ((2,square),(4,fourth)):
            for word in product(range(5),repeat=degree):
                exponent=tuple(word.count(i) for i in range(5))
                coefficient=Q(1,32)
                for i in word:coefficient*=signs[i]
                target[exponent]=target.get(exponent,Q(0))+coefficient
    square={k:v for k,v in square.items() if v}
    fourth={k:v for k,v in fourth.items() if v}
    assert len(square)==5 and all(v==1 for v in square.values())
    assert len(fourth)==15
    for exponent,value in fourth.items():
        assert value==(1 if 4 in exponent else 6)
    mu=Q(3,2);variance=Q(1,4)
    fourth_upper=3*variance**2-2*variance**2/5
    assert fourth_upper==Q(13,80)
    mean_reciprocal=mu*(Q(4,9)+Q(16,81)*(variance+fourth_upper))
    assert mean_reciprocal==Q(71,90)
    degree_cap=32*mean_reciprocal
    triangle_corner_cap=64+degree_cap/2
    assert degree_cap==Q(1136,45) and triangle_corner_cap==Q(3448,45)
    residual_cap=(triangle_corner_cap-60)/12
    ratio=27/residual_cap
    assert residual_cap==Q(187,135) and ratio==Q(3645,187)
    assert ceil(227*ratio)==4425
    assert 5*2**4==80 and 2**5==32 and 2**4==16
    # A concrete rational p example checks the reciprocal bound on all32 signs.
    p=[Q(1,2),Q(0),Q(0),Q(0),Q(0)]
    ratios=[mu-sum(s*x for s,x in zip(signs,p)) for signs in product((-1,1),repeat=5)]
    assert sum(1/r for r in ratios)/32==Q(3,4)<mean_reciprocal
    paths=[Path(__file__),ROOT/'results/right_angled_ideal6_carrier_propagation_20261018.md',
           ROOT/'results/right_angled_five_face_carrier_20261018.md']
    pins={}
    for path in paths:
        with path.open('rb') as stream:pins[str(path.relative_to(ROOT))]=file_digest(stream,'sha256').hexdigest()
    return encode(dict(status='PASS_IDEAL_RA6_CUSP_TRANSPORT_MODEL_AND_MOMENTS',
        sources=pins,frames_checked=frames,moment_monomials=dict(second=len(square),fourth=len(fourth)),
        formal_cubic_neighbors=32,cusp_triangle_edges=80,noncubic_good_degree_upper=16,
        noncubic_triangle_corner_upper=72,reciprocal_average_upper=mean_reciprocal,
        cubic_good_degree_average_upper=degree_cap,total_triangle_corner_average_upper=triangle_corner_cap,
        ridge_minus_facet_to_vertex_upper=residual_cap,vertex_to_facet_lower=ratio,
        inherited_facet_lower=227,new_vertex_lower=4425,
        scope=dict(universal_geometry_in_written_report=True,independent_review_complete=False,
            model_checks_not_all_widths_classification=True,lean_verified=False,
            new_dimension_ceiling=False,finite_volume_example=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    signal.alarm(59)
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        before=resources()
        record=run()
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        assert rss<1024**3
        path=ROOT/f'results/{TAG}.json'
        if args.verify:
            with path.open() as stream:assert json.load(stream)==record
        else:
            assert not path.exists()
            path.write_text(json.dumps(record,indent=2)+'\n')
        after=resources(start=False)
    print(json.dumps(dict(status=record['status'],frames=record['frames_checked'],
        vertex_lower=4425,peak_rss_bytes=rss,resources=after)))
