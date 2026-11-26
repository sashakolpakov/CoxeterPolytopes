#!/usr/bin/env python3
"""Check a reduced-certificate infeasibility ray as a compressed primal point.

An exact compressed point proves the certificate search infeasible. It does
not by itself prove feasibility of the original primal before dual facial
reduction. That requires a separate full-model lift/check.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, solver_gate
from scripts.compare_conic_solvers_20261126 import checked_dual_face_model


def projected_model(problem,layout,embedding,columns):
    import numpy as np
    active=np.zeros(layout['scalar_cones']+sum(n*(n+1)//2 for n in layout['PSD_sizes']),dtype=bool)
    active[embedding]=True
    scalar=0
    base=layout['scalar_cones']
    blocks=[]
    for block in problem['blocks']:
        entries=block['entries']
        n=len(entries)
        if n==1:
            if active[scalar]: blocks.append(block)
            scalar+=1
        else:
            i,j=np.triu_indices(n)
            kept=np.flatnonzero(active[base+np.flatnonzero(i==j)])
            if len(kept):
                blocks.append(dict(name=block['name'],entries=[[entries[i][j] for j in kept] for i in kept]))
            base+=n*(n+1)//2
    names=[problem['variables'][int(k)] for k in columns]
    allowed=set(names)|{''}
    if any(v not in allowed for b in blocks for row in b['entries'] for q in row for v in q):
        raise ValueError('A removed primal variable still occurs on the retained certificate face')
    return conic.model(names,blocks)


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-input',type=Path,required=True)
    parser.add_argument('--preflight',type=Path,required=True)
    parser.add_argument('--encoding',type=Path,required=True)
    parser.add_argument('--exact-model',type=Path,required=True)
    parser.add_argument('--dual-face-model',type=Path,required=True)
    parser.add_argument('--dual-face-replay',type=Path,required=True)
    parser.add_argument('--solver-output',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    started=time.monotonic()
    gate=solver_gate.verify_native(args.native_input,args.preflight)
    checked_dual_face_model(args.dual_face_model,args.dual_face_replay,gate)
    source_paths=[Path(__file__),args.solver_output/'result.json',args.solver_output/'solution.npz',
                  args.exact_model,args.dual_face_replay,args.dual_face_model/'embedding.npz',
                  ROOT/'scripts/compare_conic_solvers_20261126.py',ROOT/'research/global_formulation/conic.py']
    def bindings():
        return {str(p):{'bytes':p.stat().st_size,'sha256':solver_gate.file_digest(p)} for p in source_paths}
    before=bindings()
    args.output.mkdir()
    conic.write_json(args.output/'source_bindings.json',before)
    result=json.loads((args.solver_output/'result.json').read_text())
    if result['native']['native_status']!='PrimalInfeasible' or result['mode']!='certificate':
        raise ValueError('Expected a certificate-feasibility infeasibility ray')
    arrays=np.load(args.solver_output/'solution.npz',allow_pickle=False)
    layout=json.loads((args.dual_face_model/'layout.json').read_text())
    maps=np.load(args.dual_face_model/'embedding.npz',allow_pickle=False)
    columns=maps['retained_ordinary_equation_columns']
    embedding=maps['dual_coordinate_embedding']
    v=arrays['y']
    z=layout['zero_cones']
    if z!=len(columns)+1 or not np.isfinite(v).all() or v[z-1]<=0:
        raise ValueError('The returned ray has no valid positive normalizing coordinate')
    eta=float(v[z-1])
    x=-v[:z-1]/eta
    data,original_layout=backend.load_encoding(args.encoding)
    data,original_layout=backend.original_feasibility(data,original_layout)
    A=data['A'][embedding][:,columns].tocsc()
    b=data['b'][embedding]
    cone_layout=dict(layout,zero_cones=0)
    affine=b-A@x
    slack=v[z:]/eta
    numerical={'normalizing_coordinate':eta,
        'compressed_primal_minimum_affine_cone_eigenvalue':backend.cone_minimum(affine,cone_layout),
        'proposed_cone_slack_minimum_eigenvalue':backend.cone_minimum(slack,cone_layout),
        'maximum_affine_residual':float(np.max(np.abs(affine-slack)))}
    np.savez_compressed(args.output/'compressed_numerical_point.npz',x=x,original_columns=columns)
    proof=json.loads((args.dual_face_model/'face_proof.json').read_text())
    if solver_gate.file_digest(args.exact_model)!=proof['exact_model_sha256']:
        raise ValueError('Different exact source model')
    problem=conic.read_model(args.exact_model)
    projected=projected_model(problem,original_layout,embedding,columns)
    del problem
    raw=[F(str(float(c))) for c in x]
    identities={v:{v:F(1)} for v in projected['variables']}
    report={'status':'NO_EXACT_CERTIFICATE_FEASIBILITY_DECISION_FROM_NUMERICAL_RAY',
        'native_status':'PrimalInfeasible','numerical_ray_diagnostics':numerical,
        'projected_primal_variables':len(projected['variables']),
        'projected_primal_cones':len(projected['blocks']),
        'original_primal_feasibility_status':'UNDECIDED',
        'original_primal_point_lift_verified':False,
        'accepted_bounds_changed':False,'geometric_constraints_changed':False}
    try:
        point,checked=conic.reconstruct(projected,projected,identities,raw)
    except ValueError as error:
        report['exact_recovery_failure']=str(error)
    else:
        conic.write_json(args.output/'compressed_feasible_point.exact.json',point)
        conic.write_json(args.output/'compressed_point_check.json',checked)
        report.update(status='EXACT_REDUCED_CERTIFICATE_INFEASIBILITY_VERIFIED',
            normalized_full_certificate_search_infeasible_by_verified_face_equivalence=True,
            exact_compressed_point_sha256=checked['check']['candidate_sha256'])
    if bindings()!=before: raise ValueError('A candidate or exact checker input changed')
    report.update(elapsed_seconds=time.monotonic()-started,source_bindings=before)
    conic.write_json(args.output/'result.json',report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k!='source_bindings'})),flush=True)


if __name__=='__main__':
    main()
