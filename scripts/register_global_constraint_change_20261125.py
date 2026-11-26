#!/usr/bin/env python3
"""Check a documented fixed-scope extension before admitting a changed model."""
import argparse
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent',required=True);parser.add_argument('--input',required=True)
    parser.add_argument('--family-record',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();parent=(ROOT/args.parent).resolve();new=(ROOT/args.input).resolve()
    family_path=(ROOT/args.family_record).resolve();out=(ROOT/args.output).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT) for p in (parent,new,family_path,out)):
        raise ValueError('Use project inputs and a fresh output')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,solver_gate,strengthening
    family=json.loads(family_path.read_text())
    required={'family_id','parameter_domain','uniform_lift_argument','source_references','moment_registry',
              'implementation_modules','local_lift_receipt','rank_conditions_treatment','independence_status'}
    if not required<=set(family) or any(not family[k] for k in required):
        raise ValueError('Document scope, uniform lift, sources, shared moments, rank treatment and regression evidence')
    document=ROOT/family['uniform_lift_argument'].split('#',1)[0]
    local=ROOT/family['local_lift_receipt'];regression=json.loads(local.read_text())
    if not regression.get('status','').endswith('PASSED'):
        raise ValueError('The local lift regression has not passed')
    source_bindings=regression.get('source_bindings',{})
    paths=[parent,new,family_path,document,local,Path(__file__),ROOT/'research/global_formulation/strengthening.py',
           ROOT/'research/global_formulation/solver_gate.py']
    for name in family['implementation_modules']:
        p=ROOT/name;old=source_bindings.get(name)
        expected=old.get('sha256') if isinstance(old,dict) else old
        if expected!=solver_gate.file_digest(p):
            raise ValueError('The lift regression is stale or does not cover '+name)
        paths.append(p)
    before={str(p.relative_to(ROOT)):solver_gate.file_digest(p) for p in paths}
    P,Q=conic.read_model(parent),conic.read_model(new)
    old_contract=solver_gate.scientific_contract(P);contract=solver_gate.scientific_contract(Q)
    if any(old_contract[k]!=contract[k] for k in ('threshold','parity','upper_dimension')):
        raise ValueError('Changing dimension scope is not fixed-scope constraint addition; validate a new uniform instantiation separately')
    retained=strengthening.compare_and_test(P,Q,None)
    report=dict(schema='verified-constraint-change-v1',status='CONSTRAINT_EXTENSION_CHECKED',
        new_input_sha256=before[str(new.relative_to(ROOT))],parent_input_sha256=before[str(parent.relative_to(ROOT))],
        family_id=family['family_id'],scientific_contract=contract,all_parent_constraints_retained=True,
        exact_retention=retained,independence_status=family['independence_status'],
        local_lift_regression_current=True,geometric_truth_is_not_certified_by_registration=True,
        independent_review_required_before_bound_update=True,bound_artifacts=before)
    if any(solver_gate.file_digest(ROOT/p)!=sha for p,sha in before.items()):raise ValueError('A constraint-change input changed')
    out.mkdir();conic.write_json(out/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='bound_artifacts'}),flush=True)


if __name__=='__main__':main()
