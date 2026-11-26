#!/usr/bin/env python3
"""Check homogeneous narrow intervals preserve shared populations and exact cones."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if out.exists() or not out.is_relative_to(ROOT/'results'):raise ValueError('Use fresh results output')
    out.mkdir();sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,conditioning,equilibration
    paths=[Path(__file__),*[ROOT/'research/global_formulation'/(n+'.py') for n in
                            ('conditioning','equilibration','conic','feasibility','cone_reduction')]]
    def bindings():return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    before=bindings();delta=F(1,10**30)
    P=conic.model(['m','u'],[
        dict(name='fixture_root_inverse_lower',entries=[[conic.form(m=-2,u=1)]]),
        dict(name='fixture_root_inverse_upper',entries=[[conic.form(m=2+delta,u=-1)]]),
        dict(name='mass_lower',entries=[[conic.form(-1,m=1)]]),
        dict(name='mass_upper',entries=[[conic.form(3,m=-1)]]),
        dict(name='shared_PSD',entries=[[conic.form(m=1),conic.form(1)],
                                        [conic.form(1),conic.form(m=1)]])])
    low=P['blocks'][0]['entries'][0][0];high=P['blocks'][1]['entries'][0][0]
    old=[]
    for v in low.keys()&high.keys():
        ratio=-low[v]/high[v];width=conic.plus(low,conic.times(ratio,high))
        size=max(abs(c) for c in width.values())
        magnitude=max(max(abs(c) for c in low.values()),max(abs(ratio*c) for c in high.values()))
        old.append((size/magnitude,v))
    if min(old)[1]!='m':raise ValueError('The regression no longer exhibits the legacy bad pivot')
    Q,change,_,report=conditioning.thin_intervals(P)
    if len(report['pairs'])!=1 or report['pairs'][0]['pivot']!='u' or 'm' not in Q['variables']:
        raise ValueError('Interval conditioning replaced the shared population')
    z=report['pairs'][0]['coordinate'];values={'m':F(2),z:F(1)}
    original={v:conic.value(q,values) for v,q in change.items()}
    if original!={'m':F(2),'u':F(4)+delta}:raise ValueError('Wrong homogeneous interval inverse')
    conic.check_point(P,original);conic.check_point(Q,values)
    try:conic.check_point(Q,dict(values,**{z:F(3)}))
    except ValueError:pass
    else:raise ValueError('An out-of-interval point was admitted')
    # Two dependent selected pivot columns still define an invertible change
    # of the full affine coordinates; one auxiliary coordinate is eliminated.
    R=conic.model(['m','u','v'],[
        dict(name='a_root_inverse_lower',entries=[[conic.form(m=-2,u=1,v=1)]]),
        dict(name='a_root_inverse_upper',entries=[[conic.form(m=2+delta,u=-1,v=-1)]]),
        dict(name='b_root_inverse_lower',entries=[[conic.form(m=-4,u=2,v=2)]]),
        dict(name='b_root_inverse_upper',entries=[[conic.form(m=4+2*delta,u=-2,v=-2)]]),
        dict(name='mass',entries=[[conic.form(m=1)]]),
        dict(name='observe_remaining_coordinate',entries=[[conic.form(m=1),conic.form(v=1)],
                                                         [conic.form(v=1),conic.form(m=10)]])])
    S,change,_,report=conditioning.thin_intervals(R)
    old={'m':F(2),'u':F(3)+delta,'v':F(1)}
    point={v:old[v] for v in S['variables'] if v in old}
    for pair in report['pairs']:
        if pair['coordinate'] in S['variables']:
            point[pair['coordinate']]=conic.value(pair['defining_row'],old)/pair['delta']
    conic.check_point(R,old);conic.check_point(S,point)
    if {v:conic.value(q,point) for v,q in change.items()}!=old:
        raise ValueError('Dependent selected pivots lost the exact inverse')
    if len(S['variables'])!=len(R['variables']) or not report['all_affine_entry_inverse_identities_checked']:
        raise ValueError('The sparse transform was not verified entry by entry')
    other=conditioning.self_check();units=equilibration.self_check()
    if bindings()!=before:raise ValueError('A conditioning source changed')
    result=dict(status='SPARSE_INTERVAL_CONDITIONING_REGRESSIONS_PASSED',
        shared_population_preserved=True,legacy_bad_pivot_exhibited=True,
        dependent_selected_pivots_checked=True,all_affine_entries_verified=True,exact_original_and_transformed_points_checked=True,
        out_of_interval_point_rejected=True,affine_interval_regression=other,equilibration_regression=units,
        source_bindings=before,scientific_feasibility_claimed=False)
    conic.write_json(out/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
