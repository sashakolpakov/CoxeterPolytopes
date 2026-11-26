#!/usr/bin/env python3
"""Regressions for cascading faces, proof replay, and exact tiny coefficients."""
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'): raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if out.exists() or not out.is_relative_to(ROOT/'results'): raise ValueError('Use fresh results output')
    out.mkdir();sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,cone_reduction as cr
    def scalar(name,q): return dict(name=name,entries=[[q]])
    def reduce(P):
        source=P;maps={v:{v:F(1)} for v in P['variables']};stages=[]
        while True:
            if len(stages)>12: raise ValueError('Fixture did not converge')
            Q,proof=cr.pass_once(source)
            maps={v:cr.substitute(q,proof['substitutions']) for v,q in maps.items()}
            stages.append((source,Q,proof))
            if cr.audit(Q)['passed']: return Q,maps,stages
            source=Q
    tiny=F(1,10**40)
    P=conic.model(['x','y','z','w','unused'],[
        scalar('x',conic.form(x=1)),scalar('y',conic.form(y=1)),
        scalar('sum_zero',conic.form(x=-1,y=-1)),
        dict(name='cascade',entries=[[conic.form(y=1),conic.form(z=1)],
                                    [conic.form(z=1),conic.form(1)]]),
        scalar('w',conic.form(w=1)),scalar('duplicate',conic.form(w=2)),
        scalar('upper',conic.form(2,w=-1)),scalar('tiny_is_real',conic.form(1,w=-tiny))])
    Q,maps,stages=reduce(P)
    if Q['variables']!=['w'] or len(stages[0][2]['face_rounds'])<2:
        raise ValueError('Internal exact face propagation failed')
    try:cr.verify_source_faces(P,stages[0][2]['face_rounds'][1]['faces'])
    except ValueError:pass
    else:raise ValueError('A later face was accepted without its preceding equations')
    point={v:conic.value(q,{'w':F(1)}) for v,q in maps.items()}
    conic.check_point(P,point)
    if not any(abs(c)==tiny for b in Q['blocks'] for row in b['entries'] for q in row for c in q.values()):
        raise ValueError('An exact small coefficient was erased')
    first,second,proof=stages[-1]
    damaged=copy.deepcopy(second)
    next(b for b in damaged['blocks'] if b['name']=='upper')['entries'][0][0]['']=F(4)
    try: cr.verify_pass(first,damaged,proof)
    except ValueError: pass
    else: raise ValueError('Changed constraint accepted by the equivalence verifier')
    bad_proof=copy.deepcopy(stages[0][2])
    bad_proof['inferred_equations'].append(conic.form(w=1))
    try: cr.verify_pass(P,stages[0][1],bad_proof)
    except ValueError: pass
    else: raise ValueError('An unproved equality was accepted')
    R=conic.model(['x','y'],[
        scalar('affine',conic.form(-1,x=1,y=-1)),
        scalar('opposite',conic.form(1,x=-1,y=1)),scalar('y',conic.form(y=1)),
        scalar('upper',conic.form(3,y=-1)),
        dict(name='dependent_columns',entries=[[conic.form(x=1),conic.form(x=1)],
                                               [conic.form(x=1),conic.form(x=1)]])])
    S,subs,other=reduce(R)
    free={v:F(0) for v in S['variables']}
    conic.check_point(R,{v:conic.value(q,free) for v,q in subs.items()})
    if max(len(b['entries']) for b in S['blocks'])!=1:
        raise ValueError('Symbolic column kernel was not removed')
    impossible=conic.model(['z'],[
        dict(name='zero_diagonal',entries=[[{},conic.form(z=1)],[conic.form(z=1),conic.form(1)]]),
        scalar('must_positive',conic.form(-1,z=1))])
    try: cr.pass_once(impossible)
    except cr.Contradiction as error:
        if error.proof['kind'] not in ('negative_constant_after_elimination','negative_diagonal_constant'): raise
    else: raise ValueError('The forced-zero contradiction was missed')
    paths=[Path(__file__),ROOT/'research/global_formulation/cone_reduction.py',
           ROOT/'research/global_formulation/conic.py',ROOT/'research/global_formulation/feasibility.py']
    report=dict(status='EXACT_CONE_REDUCTION_REGRESSIONS_PASSED',cascade_passes=len(stages),
        internal_face_rounds=len(stages[0][2]['face_rounds']),
        affine_opposite_passes=len(other),unproved_equation_rejected=True,tampered_constraint_rejected=True,
        tiny_exact_coefficients_retained=True,original_witnesses_checked=True,
        later_face_requires_preceding_equalities=True,
        source_bindings={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    conic.write_json(out/'result.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':main()
