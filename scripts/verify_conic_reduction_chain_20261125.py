#!/usr/bin/env python3
"""Replay exact reduction proofs in resumable, non-authorizing checkpoints."""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def fractions(value):
    if isinstance(value,list):return [fractions(v) for v in value]
    if isinstance(value,dict):return {k:fractions(v) for k,v in value.items()}
    if isinstance(value,str):
        try:return F(value)
        except ValueError:return value
    return value


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared',required=True);parser.add_argument('--output',required=True)
    parser.add_argument('--previous',help='A preceding partial replay of this same final preparation')
    parser.add_argument('--max-stages',type=int,default=2)
    args=parser.parse_args();prepared=(ROOT/args.prepared).resolve();out=(ROOT/args.output).resolve()
    previous=(ROOT/args.previous).resolve() if args.previous else None
    if not 1<=args.max_stages<=4:raise ValueError('Use one to four complete stages per replay batch')
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in (prepared,out,*([previous] if previous else []))):
        raise ValueError('Use retained preparation and fresh output')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,cone_reduction as cr,solver_gate
    executed=[Path(__file__),*[ROOT/'research/global_formulation'/(n+'.py') for n in
                              ('conic','cone_reduction','feasibility','solver_gate')]]
    before={str(p.relative_to(ROOT)):solver_gate.file_digest(p) for p in executed}
    started=time.monotonic();prior=json.loads((prepared/'result.json').read_text())
    if not prior.get('presolve_fixed_point'):raise ValueError('Incomplete preparation cannot pass verification')
    chain=prior['reduction_chain'];chain_sha=conic.digest_bytes(conic.canonical(chain))
    original=ROOT/prior['frozen_input']
    if solver_gate.file_digest(original)!=prior['frozen_input_sha256']:raise ValueError('Changed scientific input')
    prepared_sha=solver_gate.file_digest(prepared/'result.json')
    artifacts={str(original.relative_to(ROOT)):prior['frozen_input_sha256'],
               str((prepared/'result.json').relative_to(ROOT)):prepared_sha,**before}
    def unchanged(bindings):
        for name,sha in bindings.items():
            if solver_gate.file_digest(ROOT/name)!=sha:raise ValueError('Changed replay artifact: '+name)
    out.mkdir()
    if previous:
        old=json.loads((previous/'result.json').read_text())
        if old.get('status')!='EXACT_REDUCTION_CHAIN_REPLAY_CHECKPOINT' or old.get('solver_ready') is not False:
            raise ValueError('Use a non-authorizing partial replay')
        if old['reduction_chain_sha256']!=chain_sha or old['final_preparation_receipt_sha256']!=prepared_sha:
            raise ValueError('The replay belongs to a different preparation')
        unchanged(old['bound_artifacts']);unchanged(old['checkpoint_artifacts'])
        for name,sha in before.items():
            if old['bound_artifacts'].get(name)!=sha:raise ValueError('A replay helper changed between batches')
        artifacts.update(old['bound_artifacts'])
        artifacts[str((previous/'result.json').relative_to(ROOT))]=solver_gate.file_digest(previous/'result.json')
        done=old['verified_stages']
        if not 0<done<len(chain):raise ValueError('Invalid replay stage count')
        current=ROOT/old['checkpoint_model']
        if old['checkpoint_matrix_sha256']!=chain[done-1]['target_matrix_sha256']:
            raise ValueError('Wrong replay model')
        P=conic.read_model(current)
        if conic.digest_bytes(conic.canonical(P))!=chain[done-1]['target_matrix_sha256']:
            raise ValueError('The replay checkpoint is not the previously proved target')
        composed={v:{w:F(c) for w,c in q.items()} for v,q in
                  json.loads((previous/'composed_substitutions.exact.json').read_text()).items()}
    else:
        P=conic.read_model(original);done=0
        composed={v:{v:F(1)} for v in P['variables']}
    stop=min(len(chain),done+args.max_stages)
    for stage in range(done,stop):
        item=chain[stage]
        if conic.digest_bytes(conic.canonical(P))!=item['source_matrix_sha256']:
            raise ValueError('A reduction chain skips its source')
        path=ROOT/item['proof']
        if solver_gate.file_digest(path)!=item['proof_sha256']:raise ValueError('Changed exact reduction proof')
        proof=fractions(json.loads(path.read_text()));updates=proof['substitutions'];blocks={}
        for mapping in proof['block_maps']:
            if mapping['kind']!='principal_congruence' or mapping['target'] in blocks:continue
            source=P['blocks'][mapping['source']];indices=mapping['keep'];scale=mapping['scale']
            M=[[cr.substitute(source['entries'][i][j],updates) for j in indices] for i in indices]
            blocks[mapping['target']]=dict(name=source['name'],entries=[[
                conic.times(1/scale,q) for q in row] for row in M])
        Q=conic.model([v for v in P['variables'] if updates[v]=={v:F(1)}],
                      [blocks[k] for k in range(len(blocks))])
        if conic.digest_bytes(conic.canonical(Q))!=item['target_matrix_sha256']:
            raise ValueError('The proof reconstructs a different target')
        cr.verify_pass(P,Q,proof)
        composed={v:cr.substitute(q,updates) for v,q in composed.items()}
        artifacts[str(path.relative_to(ROOT))]=item['proof_sha256'];P=Q
        print(json.dumps(dict(stage=stage+1,variables=len(P['variables']),cones=len(P['blocks']),
                              all_source_constraints_verified=True)),flush=True)
    unchanged(artifacts)
    common=dict(frozen_input=prior['frozen_input'],frozen_input_sha256=prior['frozen_input_sha256'],
        reduction_chain_sha256=chain_sha,final_preparation_receipt_sha256=prepared_sha,
        verified_stages=stop,total_stages=len(chain),bound_artifacts=artifacts,
        accepted_bounds_changed=False,elapsed_seconds=time.monotonic()-started)
    if stop<len(chain):
        state=out/'composed_substitutions.exact.json';conic.write_json(state,composed)
        model_path=out/'verified_checkpoint_model.exact.json';conic.write_json(model_path,P)
        tip=str(model_path.relative_to(ROOT))
        checkpoint={tip:solver_gate.file_digest(model_path),str(state.relative_to(ROOT)):solver_gate.file_digest(state)}
        unchanged(artifacts)
        report=dict(common,status='EXACT_REDUCTION_CHAIN_REPLAY_CHECKPOINT',solver_ready=False,
            all_exact_reduction_proofs_verified=False,checkpoint_model=tip,checkpoint_artifacts=checkpoint,
            checkpoint_matrix_sha256=chain[stop-1]['target_matrix_sha256'])
    else:
        audit=cr.audit(P);solver_gate.require_clean_audit(audit)
        if conic.digest_bytes(conic.canonical(P))!=prior['prepared_matrix_sha256']:
            raise ValueError('The verified chain does not reach the prepared SDP')
        expected={v:{w:F(c) for w,c in q.items()} for v,q in
                  json.loads((prepared/'substitutions.exact.json').read_text()).items()}
        if composed!=expected:raise ValueError('The full original-coordinate reconstruction changed')
        for path in (prepared/'compiled.exact.json',prepared/'substitutions.exact.json'):
            artifacts[str(path.relative_to(ROOT))]=solver_gate.file_digest(path)
        unchanged(artifacts)
        report=dict(common,status='EXACT_CONIC_REDUCTION_CHAIN_VERIFIED',presolve_fixed_point=True,
            prepared_matrix_sha256=prior['prepared_matrix_sha256'],
            prepared_substitutions_sha256=prior['prepared_substitutions_sha256'],
            all_exact_reduction_proofs_verified=True,all_original_coordinates_reconstructed=True,
            scientific_contract=prior['scientific_contract'],scientific_registration=prior['scientific_registration'],
            fixed_point_audit=audit,stages=len(chain),elapsed_seconds=time.monotonic()-started)
    conic.write_json(out/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('bound_artifacts','checkpoint_artifacts')}),flush=True)


if __name__=='__main__':main()
