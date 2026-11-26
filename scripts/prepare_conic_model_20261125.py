#!/usr/bin/env python3
"""Make one fully verified exact reduction checkpoint; never confuse it with readiness."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original',required=True)
    parser.add_argument('--previous')
    parser.add_argument('--registration',help='Required after changing the frozen scientific constraints')
    parser.add_argument('--output',required=True)
    args=parser.parse_args();original=(ROOT/args.original).resolve();out=(ROOT/args.output).resolve()
    previous=(ROOT/args.previous).resolve() if args.previous else None
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in (original,out,*([previous] if previous else []))):
        raise ValueError('Use retained results and fresh output')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,cone_reduction as cr,solver_gate
    started=time.monotonic();out.mkdir()
    source_path=previous/'compiled.exact.json' if previous else original
    paths=[original,source_path,Path(__file__),ROOT/'research/global_formulation/cone_reduction.py',
           ROOT/'research/global_formulation/conic.py',ROOT/'research/global_formulation/feasibility.py',
           ROOT/'research/global_formulation/solver_gate.py',ROOT/'research/global_formulation/constraint_policy.json']
    if previous:paths.extend(previous/name for name in ('result.json','substitutions.exact.json','reduction_proof.exact.json'))
    def digest(path):
        h=hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda:stream.read(1024**2),b''):h.update(chunk)
        return h.hexdigest()
    def bindings():return {str(p.relative_to(ROOT)):digest(p) for p in paths}
    before=bindings();conic.write_json(out/'source_bindings.json',before)
    registration=solver_gate.check_registered_input(before[str(original.relative_to(ROOT))],
        (ROOT/args.registration).resolve() if args.registration else None)
    def progress(stage,**details):
        item=dict(stage=stage,elapsed_seconds=time.monotonic()-started,**details)
        with (out/'progress.jsonl').open('a') as stream:stream.write(json.dumps(item)+'\n')
        print(json.dumps(item),flush=True)
    P=conic.read_model(source_path)
    source_sha=conic.digest_bytes(conic.canonical(P))
    if previous:
        prior=json.loads((previous/'result.json').read_text())
        if prior['frozen_input_sha256']!=before[str(original.relative_to(ROOT))] or prior['prepared_matrix_sha256']!=source_sha:
            raise ValueError('The checkpoint does not match the frozen input')
        old_subs={v:{w:F(c) for w,c in q.items()} for v,q in json.loads((previous/'substitutions.exact.json').read_text()).items()}
        if conic.digest_bytes(conic.canonical(old_subs))!=prior['prepared_substitutions_sha256']:
            raise ValueError('Changed original-coordinate reconstruction')
        chain=list(prior['reduction_chain'])
        contract=prior['scientific_contract']
    else:
        old_subs={v:{v:F(1)} for v in P['variables']};chain=[]
        contract=solver_gate.scientific_contract(P)
    progress('loaded_source',variables=len(P['variables']),cones=len(P['blocks']),pass_number=len(chain)+1)
    try: Q,proof=cr.pass_once(P)
    except cr.Contradiction as error:
        # Keep the exact obstruction and its source binding for independent
        # verification; a global geometric bound is never changed here.
        conic.write_json(out/'contradiction.exact.json',error.proof)
        report=dict(status='EXACT_SOURCE_OBSTRUCTION_REQUIRES_CHAIN_REVIEW',solver_ready=False,
            frozen_input_sha256=before[str(original.relative_to(ROOT))],source_matrix_sha256=source_sha,
            source_model=str(source_path.relative_to(ROOT)),reduction_chain=chain,
            accepted_bounds_changed=False,elapsed_seconds=time.monotonic()-started)
        conic.write_json(out/'result.json',report);print(json.dumps(report),flush=True);return
    progress('pass_equivalence_verified',variables=len(Q['variables']),cones=len(Q['blocks']),**proof['counts'])
    substitutions={v:cr.substitute(q,proof['substitutions']) for v,q in old_subs.items()}
    conic.write_json(out/'compiled.exact.json',Q)
    conic.write_json(out/'substitutions.exact.json',substitutions)
    conic.write_json(out/'reduction_proof.exact.json',proof)
    matrix_sha=conic.digest_bytes(conic.canonical(Q))
    chain.append(dict(source_matrix_sha256=source_sha,target_matrix_sha256=matrix_sha,
        source_model=str(source_path.relative_to(ROOT)),target_model=str((out/'compiled.exact.json').relative_to(ROOT)),
        proof=str((out/'reduction_proof.exact.json').relative_to(ROOT)),proof_sha256=digest(out/'reduction_proof.exact.json'),
        exact_equivalence_replayed=True))
    try: audit=cr.audit(Q)
    except cr.Contradiction as error:
        conic.write_json(out/'contradiction.exact.json',error.proof)
        conic.write_json(out/'result.json',dict(status='EXACT_POST_REDUCTION_OBSTRUCTION_REQUIRES_CHAIN_REVIEW',
            solver_ready=False,presolve_fixed_point=False,frozen_input=str(original.relative_to(ROOT)),
            frozen_input_sha256=before[str(original.relative_to(ROOT))],prepared_matrix_sha256=matrix_sha,
            prepared_substitutions_sha256=conic.digest_bytes(conic.canonical(substitutions)),
            reduction_chain=chain,scientific_contract=contract,scientific_registration=registration,
            accepted_bounds_changed=False))
        return
    conic.write_json(out/'fixed_point_audit.json',audit)
    progress('post_pass_audit',fixed_point=audit['passed'],**audit['counts'])
    report=dict(status='EXACT_REDUCTION_FIXED_POINT' if audit['passed'] else 'EXACT_REDUCTION_CHECKPOINT_NOT_SOLVER_READY',
        solver_ready=False,presolve_fixed_point=audit['passed'],conditioning_and_native_preflight_still_required=True,
        frozen_input=str(original.relative_to(ROOT)),frozen_input_sha256=before[str(original.relative_to(ROOT))],
        prepared_matrix_sha256=matrix_sha,prepared_substitutions_sha256=conic.digest_bytes(conic.canonical(substitutions)),
        compiled_variables=len(Q['variables']),compiled_cones=len(Q['blocks']),reduction_chain=chain,
        scientific_contract=contract,scientific_registration=registration,
        source_bindings=before,accepted_bounds_changed=False,geometric_constraints_added=False,
        elapsed_seconds=time.monotonic()-started)
    if bindings()!=before:raise ValueError('A reduction source changed')
    if solver_gate.check_registered_input(before[str(original.relative_to(ROOT))],
        (ROOT/args.registration).resolve() if args.registration else None)!=registration:
        raise ValueError('The scientific constraint registration changed')
    conic.write_json(out/'result.json',report)
    progress('checkpoint_written',status=report['status'])


if __name__=='__main__':main()
