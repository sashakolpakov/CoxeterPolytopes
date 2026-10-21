#!/usr/bin/env python3
"""All-real residual cutter gate from a three-term rational dual certificate.

This verifier uses only independent145 geometry and standard-library exact
rational intervals. The numerical residual-polytope scout is not a premise.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
from fresh_h20_h4_edge_escape_gate_20261016_review import Interval as I

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_residual_cutter_domain_20261016'
GEOMETRY='results/fresh_h20_residual_h4_corner_20261016_review.json'


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def sqrt_interval(x):
    assert x.lo>=0
    den=1<<192
    def lower(v):return isqrt((v.numerator*den*den)//v.denominator)
    a,b=lower(x.lo),lower(x.hi)
    if b*b*x.hi.denominator<x.hi.numerator*den*den:b+=1
    result=I(F(a,den),F(b,den))
    assert result.lo*result.lo<=x.lo and result.hi*result.hi>=x.hi
    return result


def run():
    g=json.loads((ROOT/GEOMETRY).read_text())
    assert g['status']=='PASS_INDEPENDENT_145_WALL_RESIDUAL_FINITE_CORNER_AND_ESCAPE_RECTANGLE'
    for p,h in g['sources'].items():assert digest(p)==h,p
    q=list(map(F,g['four_norms']));C=list(map(F,g['corner']['chart']))
    metric=[F(-1),F(1),F(1),F(1,2),F(1,2)]
    dot=lambda v,w:sum(a*b*c for a,b,c in zip(v,w,metric))
    cq=dot(C,C);assert str(cq)==g['corner']['norm'] and cq<0
    d=-cq
    assert all(v>0 for v in q)
    assert g['four_gram']==[[str(q[i] if i==j else F(0)) for j in range(4)] for i in range(4)]
    assert len(g['corner']['support16_pivots'])==16 and all(F(x)>0 for x in g['corner']['support16_pivots'])
    scales=[sqrt_interval(I(d/qi)) for qi in q]
    pairs={tuple(row['walls']):I(*row['bounds']) for row in g['new_pairs']}
    constraints=[]
    for old in (133,140):
        cp=I(*g['corner']['all145_pairs'][old]);assert cp.hi<0
        products=[pairs[(old,141+i)] for i in range(4)]
        assert all(x.hi<0 for x in products)
        coefficients=[p/cp*s for p,s in zip(products,scales)]
        assert all(x.lo>0 for x in coefficients)
        constraints.append(coefficients)
    e=[F(1),F(-1,100),F(-49,50),F(-3,10),F(-1,20)]
    ce=dot(C,e);assert ce<0
    escape=[I(*g['strict_escape']['all145_pairs'][141+i])/ce*scales[i] for i in range(4)]
    assert all(x.lo>0 for x in escape)
    multipliers=[F(529,100),F(268,25),F(759,50)]
    combined=[multipliers[0]*constraints[0][i]+multipliers[1]*constraints[1][i]-multipliers[2]*escape[i] for i in range(4)]
    assert all(combined[i].lo>1 for i in range(3)) and combined[3].lo>0
    bound=multipliers[0]+multipliers[1]-multipliers[2]
    assert bound==F(83,100)<1
    norm_upper=-d*(1-bound*bound)
    assert norm_upper<0 and 1-bound*bound==F(3111,10000)
    paths=[f'scripts/{TAG}.py',GEOMETRY,'scripts/fresh_h20_h4_edge_escape_gate_20261016_review.py']
    return dict(status='PASS_ALL_REAL_RESIDUAL_E4_CUTTERS_DOTTED_TO_WALL144',
                sources={p:digest(p) for p in paths},support='S16',caps=[141,142,143,144],
                coordinate_definition='r=C-sum(u_i*n_i), z_i=u_i*sqrt(q_i/(-C^2))>=0; positive C coefficient normalized to1',
                old_walls=[133,140],old_constraint_coefficients=[[x.saved() for x in row] for row in constraints],
                E4_constraint_coefficients=[x.saved() for x in escape],multipliers=list(map(str,multipliers)),
                combined_coefficients=[x.saved() for x in combined],sum_first_three_z_upper=str(bound),
                sum_first_three_squares_upper=str(bound*bound),corner_norm=str(cq),
                wall144_projection_norm_upper=str(norm_upper),projection_margin_fraction='3111/10000',
                proof='S16perpendicular has orthogonal basis n141..n144,C. Acuteness to the caps gives negative cap coefficients; acuteness to old133 forces a positive C coefficient for every spacelike candidate, permitting normalization. Two old acute inequalities and E4nonnegative pairing give the displayed three linear inequalities. Their positive rational combination has coefficients at least(1,1,1,0) and rightside83/100. Thus z0+z1+z2<=83/100, so sum of their squares<1. Projectionto n144perpendicular has strictlynegative square. Together with r.n144<=0 and r^2>0 this forces a strict dotted contact.',
                scope='Every real spacelike normal righttoS16, acute to the four residual caps and old133/140, with r.E4>=0, is strictlydotted towall144. Applies in particular to every145-compatible residual cutter. No field, angle, norm, height or retained-center constraint. Doesnot exclude cutters dottedto144 or cutters leavingS16; nofinitevolume conclusion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
