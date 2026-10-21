#!/usr/bin/env python3
"""Independent Fraction checks of the three/four/five-port repair union.

No producer imports, optimizer, or graph-orbit enumeration. Literal feature
construction and the geometric inequalities retain their reviewed sources.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='review_sqrt5_joint_port_repairs_20261017'
BASE='fresh_sqrt5_restricted_flag_local_cones_20261017'
STAGES=('fresh_sqrt5_joint_triangle_square_ports_20261017',
        'fresh_sqrt5_joint_four_port_diagnostic_20261017',
        'fresh_sqrt5_joint_five_port_diagnostic_20261017')


@lru_cache(None)
def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def read(tag):
    doc=json.loads((ROOT/f'results/{tag}.json').read_text())
    for path,expected in doc['sources'].items():
        assert digest(path)==expected,path
    return doc


def masses(rows):
    out={i:F(v) for i,v in rows}
    assert len(out)==len(rows) and all(i>=0 and v>=0 for i,v in out.items())
    return out


def run():
    base=read(BASE)
    first=next(x for x in base['tests'] if x['name']=='legacy_export_seed')
    alpha=masses(first['face_multipliers'])
    initial=masses(first['vertex_multipliers'])
    numerator=F(first['objective_numerator'])
    assert len(alpha)==2219 and len(initial)==593
    assert numerator/sum(initial.values())==F(first['bound'])>12
    results=[]
    prior_bound=F(first['bound'])
    all_critical=[]
    for stage,tag in enumerate(STAGES,3):
        doc=read(tag)
        assert masses(doc['face_multipliers'])==alpha
        assert F(doc['objective_numerator'])==numerator
        assert F(doc['initial_bound'])==prior_bound
        if 'initial_vertex_multipliers' in doc:
            assert masses(doc['initial_vertex_multipliers'])==initial
        expected=dict(initial)
        loads={}
        if stage<5:
            for v,orbit,s,q,slack in doc['all_orbit_loads_and_slacks']:
                s,q=F(s),F(q)
                assert s>=0 and q>=0 and F(slack)==s+q-3*initial[v]
                key=(v,orbit,0)
                assert key not in loads
                loads[key]=(s+q,3)
                expected[v]=max(expected[v],(s+q)/3)
        else:
            facets=doc['facets']
            assert len(facets)==15
            for v,orbit,values in doc['all_orbit_histogram_loads']:
                hist=list(map(F,values))
                assert len(hist)==8 and min(hist)>=0
                for f in facets:
                    co=f['coefficients'];cap=f['capacity']
                    assert len(co)==8 and min(co)>=0 and cap>0
                    load=sum((a*b for a,b in zip(co,hist)),F())
                    key=(v,orbit,f['index'])
                    assert key not in loads
                    loads[key]=(load,cap)
                    expected[v]=max(expected[v],load/cap)
        actual=masses(doc['vertex_multipliers'])
        assert actual==expected and set(actual)==set(initial)
        critical=set()
        for col in doc['critical_columns']:
            orbit=col[{3:'three_port_orbit',4:'four_port_orbit',5:'five_port_orbit'}[stage]]
            key=(col['vertex'],orbit,col.get('facet',0))
            assert key not in critical
            critical.add(key)
            load,cap=loads[key]
            assert cap==col['vertex_cost'] and load==F(col['load'])
            assert F(col['old_vertex_mass'])==initial[key[0]] and load>cap*initial[key[0]]
            terms=masses(col['positive_face_coefficients'])
            assert set(terms)<=set(alpha)
            assert sum((alpha[fi]*a for fi,a in terms.items()),F())==load
            all_critical.append((key[0],load,cap))
        assert critical=={key for key,(load,cap) in loads.items() if load>cap*initial[key[0]]}
        den=sum(actual.values())
        assert den==F(doc['normalizer']) and numerator/den==F(doc['bound'])>12
        assert F(doc['bound'])<=prior_bound
        assert len(critical)==doc['counts']['violated_columns']
        assert sum(actual[v]>initial[v] for v in actual)==doc['counts']['increased_vertex_masses']
        for v,load,cap in all_critical:
            assert load<=cap*actual[v]
        results.append(dict(stage=stage,tested_inequalities=len(loads),
            verified_critical_columns=len(critical),vertex_mass_increases=doc['counts']['increased_vertex_masses'],
            initial_bound=str(prior_bound),bound=doc['bound'],normalizer=str(den)))
        initial=actual;prior_bound=F(doc['bound'])
    paths=[f'scripts/{TAG}.py',f'results/{BASE}.json',
           'results/review_sqrt5_restricted_flag_local_cover_20261017.json']
    for tag in STAGES:
        paths.extend((f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md'))
    return dict(status='PASS_INDEPENDENT_FRACTION_JOINT_PORT_REPAIR_UNION',
        sources={p:digest(p) for p in paths},stages=results,
        common_positive_faces=len(alpha),common_vertex_forests=len(initial),
        objective_numerator=str(numerator),final_bound=str(prior_bound),
        cumulative_verified_critical_columns=len(all_critical),
        scope='Standalone Fraction recount of all saved loads, critical positive face coefficients, per-vertex maxima, unchanged face masses, inherited repair preservation and exact ratios. No producer imports or optimizer. This does not independently reconstruct forest orbits, raw prism Grams or geometric facet validity; those retain explicit source/proof reviews and producer exact enumerations.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],result['final_bound'],flush=True)
