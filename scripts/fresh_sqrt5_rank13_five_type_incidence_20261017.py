#!/usr/bin/env python3
"""Exact joint short-face accounting for five explicitly classified ADE types."""
import argparse,hashlib,json,math
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_five_type_incidence_20261017'
PARAMS=['fresh_sqrt5_rank13_irreducible_quad_parameters_20261017','fresh_sqrt5_rank13_small_d_quad_parameters_20261017']
QUADS=['fresh_sqrt5_rank13_irreducible_joint_quads_20261017','fresh_sqrt5_rank13_small_d_joint_quads_20261017']
TRI='fresh_sqrt5_rank13_high_ade_triangles_20261017'
THIRD='fresh_sqrt5_rank13_triangle_third_costs_20261017'
END='fresh_sqrt5_rank13_high_ade_endpoints_20261017'
IDS=[0,1,35,53,66]

def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d

def charges(triangle,quad):
    rows=[]
    for pair in combinations(range(13),2):
        a=F(1,2)if pair in triangle else F(0);b=F(1,4)if pair in quad else F(0)
        rows.append(dict(ports=list(pair),triangle=str(a),quadrilateral=str(b),maximum=str(max(a,b))))
    return rows,sum((F(x['maximum'])for x in rows),F(0))

def run():
    parameters=[read(t)for t in PARAMS];quads=[read(t)for t in QUADS]
    tri,third,end=map(read,[TRI,THIRD,END]);allowed=set(IDS)
    tlookup={x['index']:x for x in tri['targets']};refined={x['index']:x for x in third['targets']}
    elookup={x['index']:x for x in end['targets']};targets=[]
    for pi,doc in enumerate(parameters):
        qlookup={x['index']:x for x in quads[pi]['targets']}
        for p in doc['targets']:
            index=p['index'];assert index in allowed;q=qlookup[index]
            assert p['forest']==q['forest']==elookup[index]['forest']
            accepted=p['accepted_map_indices']
            assert accepted==[i for i,m in enumerate(q['maps'])if doc['records'][m['record']]['accepted']]
            unrestricted_quad={tuple(q['maps'][i]['ports'])for i in accepted}
            assert unrestricted_quad=={tuple(x)for x in p['possible_quad_corners']}
            occurrences=tlookup[index]['accepted_occurrences']
            unrestricted_tri={tuple(sorted((a['port'],a['other'])))for a in occurrences if a['third_exceptional_component']!='H4'}
            assert unrestricted_tri=={tuple(x)for x in refined[index]['half_charge_corners']}
            filtered=[i for i in accepted if set(q['maps'][i]['adjacent_destination_types']+[q['maps'][i]['opposite_destination_type']])<=allowed]
            quad={tuple(q['maps'][i]['ports'])for i in filtered}
            triangle_occurrences=[i for i,a in enumerate(occurrences)if a['destination']in allowed and a['third_exceptional_component']!='H4']
            triangle={tuple(sorted((occurrences[i]['port'],occurrences[i]['other'])))for i in triangle_occurrences}
            rows,alpha=charges(triangle,quad);_,unrestricted=charges(unrestricted_tri,unrestricted_quad)
            same=[i for i in accepted if set(q['maps'][i]['adjacent_destination_types']+[q['maps'][i]['opposite_destination_type']])=={index}]
            destination_counts=Counter(d for i in accepted for d in q['maps'][i]['adjacent_destination_types']+[q['maps'][i]['opposite_destination_type']])
            targets.append(dict(index=index,forest=p['forest'],baseline=elookup[index]['baseline'],
                parameter_source=PARAMS[pi],quad_source=QUADS[pi],
                retained_quad_map_indices=filtered,retained_triangle_occurrence_indices=triangle_occurrences,
                quad_corners=[list(x)for x in sorted(quad)],triangle_corners=[list(x)for x in sorted(triangle)],
                quad_corner_count=len(quad),triangle_corner_count=len(triangle),corner_charges=rows,
                special_vertex_upper=str(alpha),margin_below12=str(12-alpha),
                all215_upper_lists=dict(quad_corner_count=len(unrestricted_quad),triangle_corner_count=len(unrestricted_tri),
                    corner_charge_upper=str(unrestricted),retained_literal_maps=len(accepted)),
                same_type_quad_corner_count=len({tuple(q['maps'][i]['ports'])for i in same}),
                destination_catalogue_counts=[dict(index=d,forest=elookup[d]['forest'],baseline=elookup[d]['baseline'],
                    map_corner_occurrences=count)for d,count in sorted(destination_counts.items())]))
    targets.sort(key=lambda t:IDS.index(t['index']))
    costs=[F(t['special_vertex_upper'])for t in targets]
    assert costs==[F(9,4),F(17,4),F(39,4),F(13,2),F(21,2)]
    n=13;corners=math.comb(n,2);threshold=F(5*((n-1)**2//4),2)-corners
    other=corners+n//4;assert (threshold,other)==(12,81)
    coefficients=[threshold-a for a in costs]+[other-threshold]
    assert all(x>0 for x in coefficients)
    scale=math.lcm(*(x.denominator for x in coefficients));ints=[int(scale*x)for x in coefficients]
    gcd=math.gcd(*ints);ints=[x//gcd for x in ints];assert ints==[39,31,9,22,6,276]
    face_payments=[]
    for selected in range(3):
        pay=F(3-selected)+(1 if selected==2 else 0);assert pay>=2
        face_payments.append(dict(face='triangle',selected=selected,payment=str(pay)))
    for selected in range(5):
        pay=F(1)if selected==4 else F(4-selected);assert pay>=1
        face_payments.append(dict(face='quadrilateral',selected=selected,payment=str(pay)))
    paths=[f'scripts/{TAG}.py']+[f'results/{t}.json'for t in PARAMS+QUADS+[TRI,THIRD,END]]
    return dict(status='PASS_EXACT_FIVE_TYPE_RANK13_JOINT_INCIDENCE',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=n,
        special_indices=IDS,targets=targets,other_vertex_upper=other,strict_h_threshold=str(threshold),
        frequency_coefficients=ints,
        frequency_inequality='39*V_D13 +31*V_A13 +9*V_D5_2D4 +22*V_2D5_A3 +6*V_3D4_A1 <276*V_other',
        face_payments=face_payments,
        hypotheses=['P is compact, hyperbolic Coxeter, arithmetic of simplest type over Q(sqrt5), of dimension13.',
            'Special is exactly the five whole vertex forests shown; every other type is included in V_other.',
            'The frozen complete finite-corner, arithmetic-anchor and refined triangle certificates are used with their stated hypotheses.'],
        proof='A triangle cannot have three special vertices because its normalized Gram would be rational, contradicting nonphysical positivity. Give one unit to every nonspecial corner and an extra one on each H4 framing5 pair. For a two-special triangle the extra H4 unit, when present, pays the missing unit; otherwise give one half to each special corner. Give one quarter at every corner of an all-special quadrilateral. Restrict both triangle destinations and all three quad destinations to the actual five-type set. Each root pair determines a single incident face, so take the maximum of its possible triangle and quad charges. The exact totals are9/4,17/4,39/4,13/2,21/2, while every other vertex costs at most78+3=81. Combining with12*V<2*T+Q yields the displayed strict inequality.',
        diagnostic_scope='The all215 lists and destination catalogue counts illustrate possible local migration; their counts are enumeration multiplicities, not frequencies in an actual polytope. Their large uniform corner bounds do not rule out improved joint payments.',
        scope='A simultaneous frequency restriction on five whole vertex types. In particular no compact arithmetic Qsqrt5 polytope of dimension13 has only these types. Other high-cost ADE and non-ADE types remain possible. No new dimension bound or global realization of retained local patterns is claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['frequency_inequality'],flush=True)
