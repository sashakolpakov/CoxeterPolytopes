#!/usr/bin/env python3
"""Exact composition of the two-type rank-13 short-face incidence bound."""
import argparse,hashlib,json,math
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_irreducible_joint_incidence_20261017'
PARAM='fresh_sqrt5_rank13_irreducible_quad_parameters_20261017'
THIRD='fresh_sqrt5_rank13_triangle_third_costs_20261017'
ENDPOINT='fresh_sqrt5_rank13_high_ade_endpoints_20261017'

def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d['sources'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d

def run():
    par,third,end=map(read,[PARAM,THIRD,ENDPOINT]);n=13
    assert third['dimension']==n and len(end['targets'])==215
    targets=[]
    for t in par['targets']:
        assert t['forest'] in [['D13'],['A13']]
        tr=next(q for q in third['targets']if q['index']==t['index'])
        assert tr['forest']==t['forest']
        quad={tuple(p)for p in t['possible_quad_corners']}
        tri={tuple(p)for p in tr['half_charge_corners']}
        assert tri<=quad
        charges=[]
        for pair in combinations(range(n),2):
            a=F(1,2)if pair in tri else F(0)
            b=F(1,4)if pair in quad else F(0)
            charges.append(dict(ports=list(pair),triangle=str(a),quadrilateral=str(b),maximum=str(max(a,b))))
        alpha=sum((F(x['maximum'])for x in charges),F(0))
        assert alpha==F(len(tri)+len(quad),4)
        targets.append(dict(forest=t['forest'],index=t['index'],triangle_corners=len(tri),quad_corners=len(quad),
            corner_charges=charges,special_vertex_upper=str(alpha)))
    costs={t['forest'][0]:F(t['special_vertex_upper'])for t in targets}
    assert costs=={'D13':F(10),'A13':F(27,4)}
    corners=math.comb(n,2);max_h4=n//4;other=corners+max_h4
    threshold=F(5*((n-1)**2//4),2)-corners
    assert (corners,max_h4,other,threshold)==(78,3,81,12)
    # Every actual corner belongs to one face, so triangle/quad charges take
    # their maximum at that pair; they must not be added at the same pair.
    triangle_cases=[]
    for selected in range(3):
        for h4 in (False,True):
            pay=F(3-selected)
            if selected==2:pay+=1  # one extra at H4, or halves at selected ends
            assert pay>=2
            triangle_cases.append(dict(selected=selected,third_H4=h4,payment=str(pay)))
    quad_cases=[]
    for selected in range(5):
        pay=F(1)if selected==4 else F(4-selected)
        assert pay>=1;quad_cases.append(dict(selected=selected,payment=str(pay)))
    coefficients=[threshold-costs['D13'],threshold-costs['A13'],F(other)-threshold]
    scale=math.lcm(*(x.denominator for x in coefficients))
    integers=[int(scale*x)for x in coefficients]
    div=math.gcd(*integers);integers=[x//div for x in integers]
    assert integers==[8,21,276]
    paths=[f'scripts/{TAG}.py']+[f'results/{t}.json'for t in[PARAM,THIRD,ENDPOINT]]
    return dict(status='PASS_EXACT_JOINT_D13_A13_INCIDENCE',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=n,
        special_set=[['D13'],['A13']],targets=targets,other_vertex_upper=other,
        other_vertex_accounting=dict(all_root_pairs=corners,maximum_H4_components=max_h4,extra_per_H4_framing_5_pair=1),
        triangle_cases=triangle_cases,quadrilateral_cases=quad_cases,strict_h_threshold=str(threshold),
        frequency_coefficients=integers,frequency_inequality='8*V_D13 +21*V_A13 <276*V_other',
        hypotheses=['P is compact, hyperbolic Coxeter, arithmetic of simplest type over Q(sqrt5), of dimension13.',
            'Special means the entire vertex forest is D13 or A13; every other vertex is counted in V_other.',
            'The complete 215-type endpoint/triangle lists are used as upper lists for this smaller special set.',
            'The frozen finite-anchor parameter gate and its stated opposite-incidence and arithmetic lemmas apply.'],
        proof='No triangle has three special vertices by the rational-normalized-Gram/conjugate-positivity lemma. Give each nonspecial corner one unit, plus one extra on its H4 framing5 pair. A two-special triangle with an H4 third component is then paid2; every other two-special triangle gets another half at each special corner. An all-special quadrilateral gets one quarter at each corner. Other triangles and quadrilaterals are paid by their nonspecial corners. The exact root-pair lists bound special vertex totals by10 and27/4, and the nonspecial total by78+floor(13/4)=81. The strict h-vector inequality12*V<2*T+Q implies the displayed frequency inequality.',
        scope='A joint frequency restriction on two vertex types. It excludes a compact arithmetic Qsqrt5 polytope of dimension13 with only these two types, but does not exclude either type from a general polytope, control all other costly types, or prove a dimension bound.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['frequency_inequality'],flush=True)
