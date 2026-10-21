#!/usr/bin/env python3
"""Exact local obstruction to the escaping single-inner-weight root plane."""
import argparse
from itertools import product
import json
from fresh_h9_glued_host_20261013_independent import *

OUT=ROOT/'results/fresh_h9_quartic_mixed_20261013_local_obstruction.json'


def residue(x):
    assert integral(x)
    return tuple(int(z)%3 for z in x)


def times(x,y):
    c=[0]*7
    for i in range(4):
        for j in range(4):c[i+j]+=x[i]*y[j]
    for i in range(6,3,-1):c[i-2]+=4*c[i];c[i-4]-=2*c[i]
    return tuple(z%3 for z in c[:4])


def run():
    rows=[]
    for q in POOL:
        valuation=0;part=q
        while integral(divide(part,T)):
            valuation+=1;part=divide(part,T)
        solutions=[list(x)for x in product(range(3),repeat=4)
                   if times(residue(neg(ALPHA)),times(x,x))==residue(q)]
        rows.append({'norm':enc(q),'T_adic_valuation':valuation,
                     'mod3_time_solutions':solutions,
                     'excluded_at_2':valuation in(1,3),
                     'excluded_at_3':not solutions})
    assert [r['T_adic_valuation']for r in rows]==[0,1,2,3,4]
    assert all(r['excluded_at_2']or r['excluded_at_3']for r in rows)
    assert all(not r['mod3_time_solutions']for r in rows if r['T_adic_valuation']in(0,2,4))
    return {'status':'PASS_GLOBAL_SINGLE_WEIGHT_PLANE_ROOT_OBSTRUCTION',
            'form':'6*u^2-(T-1)*x^2','norm_pool':rows,
            'dyadic_argument':'6u^2 is divisible by T^4 because (2)=(T)^4 up to a unit. The coefficient T-1 is a unit. If v_T(x)<2, the norm has even valuation 0 or 2; otherwise it has valuation at least 4. Normalized norm classes of valuations 1 and 3 are impossible.',
            'mod3_argument':'For the remaining normalized classes, reducing modulo 3 requires -(T-1)*x^2=q. All 81 residues in Z[T]/(3,T^4-4T^2+2) were enumerated; none solves this equation.',
            'scope':'No primitive integral reflection normal of the full quartic host lies in span(time, a single inner F4 fundamental weight), with all other coordinates zero. This excludes only that binary plane, not neighboring-weight slices or the full host.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();result=run()
    if a.verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
