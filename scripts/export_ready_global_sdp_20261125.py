#!/usr/bin/env python3
"""Condition a verified fixed point and export an input bound to a mandatory launch gate."""
import argparse
from fractions import Fraction as F
import gzip
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared',required=True);parser.add_argument('--verified',required=True)
    parser.add_argument('--reference',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();prepared=(ROOT/args.prepared).resolve();verified=(ROOT/args.verified).resolve()
    reference=(ROOT/args.reference).resolve();out=(ROOT/args.output).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in (prepared,verified,reference,out)):
        raise ValueError('Use retained inputs and fresh output')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,cone_reduction as cr,conditioning,equilibration,slater_margin,solver_gate
    executed_sources=[Path(__file__),*[ROOT/'research/global_formulation'/(n+'.py') for n in
        ('conic','cone_reduction','feasibility','conditioning','equilibration','slater_margin','solver_gate',
         'sdpa_cli','checkpoint_solver')],ROOT/'scripts/run_sdpa_checkpoint_20261125.py']
    before={str(p.relative_to(ROOT)):solver_gate.file_digest(p) for p in executed_sources}
    started=time.monotonic();out.mkdir()
    def progress(stage,**details):
        item=dict(stage=stage,elapsed_seconds=time.monotonic()-started,**details)
        with (out/'progress.jsonl').open('a') as stream:stream.write(json.dumps(item)+'\n')
        print(json.dumps(item),flush=True)
    check=json.loads((verified/'result.json').read_text())
    if check.get('status')!='EXACT_CONIC_REDUCTION_CHAIN_VERIFIED':raise ValueError('The reduction chain has not been verified')
    solver_gate.require_clean_audit(check['fixed_point_audit'])
    artifacts=dict(check['bound_artifacts'])
    for name,sha in artifacts.items():
        if solver_gate.file_digest(ROOT/name)!=sha:raise ValueError('A verified artifact changed')
    P=conic.read_model(prepared/'compiled.exact.json')
    if conic.digest_bytes(conic.canonical(P))!=check['prepared_matrix_sha256']:raise ValueError('Changed prepared model')
    substitutions={v:{w:F(c) for w,c in q.items()} for v,q in json.loads((prepared/'substitutions.exact.json').read_text()).items()}
    with (gzip.open(reference,'rt') if reference.suffix=='.gz' else reference.open()) as stream:
        supplied={v:F(c) for v,c in json.load(stream).items()}
    progress('verified_input_loaded',variables=len(P['variables']),cones=len(P['blocks']))
    P,change,_,intervals=conditioning.thin_intervals(P,progress=lambda **d:progress(**d))
    progress('interval_inverse_identities_verified',pairs=len(intervals['pairs']))
    substitutions={v:cr.substitute(q,change) for v,q in substitutions.items()}
    for item in intervals['pairs']:
        supplied[item['coordinate']]=conic.value(item['defining_row'],supplied)/item['delta']
    # A failed candidate only selects units. It adds no feasible-point restrictions.
    reference_units={v:max(F(1,10**6),min(F(10**6),abs(supplied.get(v,F(1))))) for v in P['variables']}
    P,change,_,_,units=equilibration.equilibrate(P,reference_units)
    progress('unit_inverse_identities_verified')
    substitutions={v:cr.substitute(q,change) for v,q in substitutions.items()}
    conic.write_json(out/'equilibrated.exact.json',P)
    conic.write_json(out/'equilibrated_substitutions.exact.json',substitutions)
    conic.write_json(out/'interval_conditioning.exact.json',intervals)
    conic.write_json(out/'equilibration.exact.json',units)
    checkpoint=dict(status='EQUIVALENT_REPRESENTATION_AWAITING_FINAL_AUDIT',solver_ready=False,
        frozen_input_sha256=check['frozen_input_sha256'],verified_reduction_receipt=str((verified/'result.json').relative_to(ROOT)),
        all_scaling_inverse_identities_verified=True,source_bindings=before,
        bound_artifacts={str((out/name).relative_to(ROOT)):solver_gate.file_digest(out/name) for name in
            ('equilibrated.exact.json','equilibrated_substitutions.exact.json','interval_conditioning.exact.json','equilibration.exact.json')})
    conic.write_json(out/'representation_checkpoint.json',checkpoint)
    progress('equivalent_representation_checkpoint_written')
    audit=cr.audit(P)
    progress('post_conditioning_audit_completed',passed=audit['passed'],**audit['counts'])
    conic.write_json(out/'fixed_point_audit.json',audit)
    if not audit['passed']:
        conic.write_json(out/'result.json',dict(status='POST_CONDITIONING_REDUCTION_REQUIRED',solver_ready=False,audit=audit))
        raise ValueError('Conditioning exposed another exact reduction; no solver input is authorized')
    if not P['variables']:
        point={v:conic.value(q,{}) for v,q in substitutions.items()}
        original=conic.read_model(ROOT/check['frozen_input'])
        checked=conic.check_point(original,point)
        conic.write_json(out/'original_feasible_point.exact.json',point)
        conic.write_json(out/'result.json',dict(status='EXACT_ORIGINAL_FEASIBLE_BY_ELIMINATION',check=checked))
        return
    smallest,largest,coefficient_count=None,F(0),0
    for block in P['blocks']:
        for row in block['entries']:
            for q in row:
                for c in q.values():
                    if c:
                        magnitude=abs(c);coefficient_count+=1
                        smallest=magnitude if smallest is None else min(smallest,magnitude)
                        largest=max(largest,magnitude)
    scale_report=dict(full_matrix_affine_coefficients=coefficient_count,
        minimum_nonzero_absolute_coefficient=smallest,maximum_absolute_coefficient=largest,
        coefficient_range_is_not_a_condition_number=True,
        scalar_cones=sum(len(b['entries'])==1 for b in P['blocks']),
        PSD_blocks=sum(len(b['entries'])>1 for b in P['blocks']),
        largest_PSD_block=max(len(b['entries']) for b in P['blocks']))
    conic.write_json(out/'representation_statistics.json',scale_report)
    Q=slater_margin.augmented(P);slater_margin.write_input(out/'input.dat-s',Q)
    seed=slater_margin.strict_initial(P,out/'initial.dat')
    conic.write_json(out/'initialization.json',seed)
    for path in (verified/'result.json',reference,prepared/'compiled.exact.json',prepared/'substitutions.exact.json',
                 out/'equilibrated.exact.json',out/'equilibrated_substitutions.exact.json',
                 out/'interval_conditioning.exact.json',out/'equilibration.exact.json',out/'fixed_point_audit.json',
                 Path(__file__),ROOT/'research/global_formulation/cone_reduction.py',
                 ROOT/'research/global_formulation/solver_gate.py',ROOT/'research/global_formulation/constraint_policy.json',
                 ROOT/'research/global_formulation/conditioning.py',ROOT/'research/global_formulation/equilibration.py',
                 ROOT/'research/global_formulation/slater_margin.py'):
        artifacts[str(path.relative_to(ROOT))]=solver_gate.file_digest(path)
    if any(solver_gate.file_digest(ROOT/name)!=sha for name,sha in before.items()):
        raise ValueError('An export helper changed during conditioning')
    artifacts.update(before)
    artifacts.update(check['scientific_registration']['bound_artifacts'])
    receipt=dict(schema='solver-ready-v1',status='SOLVER_READY',presolve_fixed_point=True,
        all_exact_reduction_proofs_verified=True,all_scaling_inverse_identities_verified=True,
        fixed_point_audit=audit,scientific_contract=check['scientific_contract'],
        frozen_input_sha256=check['frozen_input_sha256'],bound_artifacts=artifacts,
        native_input_sha256=solver_gate.file_digest(out/'input.dat-s'),
        native_free_variables=len(Q['variables']),native_block_sizes=solver_gate.native_shape(out/'input.dat-s')[1],
        auxiliary_problem='common PSD margin',accepted_bounds_changed=False)
    conic.write_json(out/'solver_ready.json',receipt)
    gate=solver_gate.verify_native(out/'input.dat-s',out/'solver_ready.json')
    report=dict(status='VERIFIED_FIXED_POINT_SDP_EXPORTED',solver_ready=True,
        frozen_input_sha256=check['frozen_input_sha256'],prepared_matrix_sha256=conic.digest_bytes(conic.canonical(P)),
        prepared_substitutions_sha256=conic.digest_bytes(conic.canonical(substitutions)),
        variables=len(P['variables']),cones=len(P['blocks']),strict_auxiliary_seed=seed,preflight=gate,
        numerical_representation=scale_report,
        constraints_added_or_dropped_without_equivalence_proof=False,accepted_bounds_changed=False,
        elapsed_seconds=time.monotonic()-started)
    conic.write_json(out/'result.json',report);print(json.dumps(conic.plain(report)),flush=True)


if __name__=='__main__':main()
