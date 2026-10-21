#!/usr/bin/env python3
"""Exact composition of the reviewed order10 payments and rank7 curvature.

This checker binds written geometric inputs. It does not formalize them or
repeat historical rank enumerations, geometric probes, or optimization.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cusp_minimum_factor_rank9_dimension29_20261017'
PACKAGES=(
 'fresh_cusp_factor_payment_framework_20261017',
 'fresh_cusp_spectator_role_payments_20261017',
 'fresh_simple_cusp_vinberg_weight_frontier_20261017',
 'fresh_cusp_link_shape_curvature_20261017',
 'fresh_cusp_rank7_fine_curvature_20261017',
 'fresh_affine_carriers_order10_20261017',
 'fresh_affine_order10_quad_payment_20261017',
)
WRITTEN=(
 'fresh_simple_finite_quad_vinberg_transfer_20261017',
 'fresh_vinberg_finite_quad_affine_cutoff9_20261017',
)

def digest(p):
    h=sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def mul(a,b):
    c=[F()for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def sub(a,b):
    return [(a[i]if i<len(a)else F())-(b[i]if i<len(b)else F())for i in range(max(len(a),len(b)))]

def run():
    paths=[Path(__file__)];statuses={};docs={}
    for tag in PACKAGES:
        p=ROOT/f'results/{tag}.json';d=json.loads(p.read_text());docs[tag]=d
        assert d['status'].startswith('PASS');statuses[tag]=d['status']
        for source,h in d.get('sources',{}).items():assert digest(ROOT/source)==h,(tag,source)
        if 'source_sha256'in d:assert digest(ROOT/f'scripts/{tag}.py')==d['source_sha256']
        paths+=[p,ROOT/f'results/{tag}.md',ROOT/f'scripts/{tag}.py']
    paths += [ROOT/f'results/{tag}.md'for tag in WRITTEN]
    assert all(p.is_file()for p in paths)
    # Coefficientwise identities, ascending order in t>=0, for n=30/31+2t.
    even=sub(mul([30,2],[24,2]),[720,48]);assert even==[0,60,4]
    odd=sub(mul([30,2],[26,2]),[744,48]);assert odd==[36,64,4]
    assert all(a>=0 for a in even+odd)and odd[0]>0
    curve=docs['fresh_cusp_rank7_fine_curvature_20261017']
    assert curve['minimum_factor_rank']==7 and curve['ambient_dimension_at_least']==30
    assert F(curve['curvature_increment'])==F(1,18)>0
    assert sub([180,40],[180,36])==[0,4] # cleared cap slack5/18
    assert 16*F(9,4)**2==81 and 16*F(25,4)*F(21,4)==525
    assert 81*7-525==42>0
    assert curve['coefficients_8_beta_minus_3n_even']==even
    assert curve['coefficients_8_beta_minus_3n_odd']==odd
    profile_checks=[]
    for row in curve['minimal_profiles']:
        rr=row['ranks'];h=[1]
        for r in rr:
            hh=[0]*(len(h)+r)
            for i,a in enumerate(h):
                for j in range(r+1):hh[i+j]+=a
            h=hh
        d=sum(rr);V=prod(r+1 for r in rr)
        assert h==row['h_coefficients']and sum(h)==V
        kap=F(5*sum(a*min(j,d-j)for j,a in enumerate(h))-d*V,2*V)
        assert kap==F(row['normalized_curvature'])and kap-3*len(rr)>0
        profile_checks.append(dict(ranks=rr,margin=str(kap-3*len(rr))))
    assert len(profile_checks)==21 and min(F(q['margin'])for q in profile_checks)==F(41,256)
    carrier=docs['fresh_affine_carriers_order10_20261017']
    ideal=docs['fresh_affine_order10_quad_payment_20261017']
    assert carrier['minimum_affine_order']==ideal['minimum_affine_order']==10
    assert carrier['forbidden_affine_order_through']==ideal['forbidden_affine_order_through']==9
    assert 9>=curve['minimum_factor_rank']
    # Costs are linear in the affine factor's order m. End corrections are
    # precisely the inherited rule's extra fork/end-segment pairs.
    costs={'cycle':([0,3],0),'path':([-6,3],6),
           'single_fork':([-5,3],5),'double_fork':([-4,3],4)}
    for coeff,end in costs.values():assert [coeff[0]+end,coeff[1]]==[0,3]
    return dict(status='PASS_ALL_RANK_MINIMUM_CUSP_FACTOR_RANK9_DIMENSION29_COMPOSITION',
       sources={str(p.relative_to(ROOT)):digest(p)for p in paths},input_statuses=statuses,
       hypothesis='Every irreducible affine cusp factor has rank at least9.',
       minimum_cusp_factor_rank=9,dimension_bound=29,arithmetic_hypothesis=False,factor_count_limit=None,
       affine_order_exclusion_through=9,finite_cost_coefficient=3,ideal_cost_per_factor_per_link_vertex=3,
       even_finite_slack_times8_coefficients=list(map(str,even)),odd_finite_slack_times8_coefficients=list(map(str,odd)),
       curvature_minimum_rank_bound_used=7,curvature_increment='1/18',curvature_minimal_profile_checks=profile_checks,
       affine_cost_coefficients={k:dict(distance_coefficients=c,end_correction=e,total_coefficients=[0,3])for k,(c,e)in costs.items()},
       global_strict_slack='5*(n-1)>0 for every n>=30',
       geometric_steps=['actual cusp-factor affine exclusion through9','finite and varying-factor triangle payment',
         'ideal quadrilateral payment from order10','order10 complement-carrier exclusions',
         'published finite-quadrilateral cases with affine cutoff9'],
       enumeration_cutoff_used=False,unrestricted163_tail_used=False,optimizer_invocations=0,new_example=None,
       scope='Exact dependency/scalar composition with written internally reviewed geometry; no rederived published case tables, formal verification, or historical large census replay.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],len(data['sources']),'bindings',flush=True)
