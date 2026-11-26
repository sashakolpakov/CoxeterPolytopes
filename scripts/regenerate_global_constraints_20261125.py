#!/usr/bin/env python3
"""Regenerate the complete level-7 constraints and compare the frozen model."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'): raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--frozen',default='results/global_edge_closure_current/necessary_primal.exact.json')
    args=parser.parse_args();out=(ROOT/args.output).resolve();frozen=(ROOT/args.frozen).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in (out,frozen)):
        raise ValueError('Use retained frozen input and fresh output')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,data,model
    started=time.monotonic();out.mkdir()
    paths=[Path(__file__),frozen,ROOT/data.LITERAL,*sorted((ROOT/'research/global_formulation').glob('*.py'))]
    def bindings():
        answer={}
        for p in paths:
            h=hashlib.sha256()
            with p.open('rb') as stream:
                for chunk in iter(lambda:stream.read(1024**2),b''):h.update(chunk)
            answer[str(p.relative_to(ROOT))]=h.hexdigest()
        return answer
    before=bindings();conic.write_json(out/'source_bindings.json',before)
    old=json.loads(frozen.read_text())
    expected=conic.digest_bytes(conic.canonical(old));del old
    regenerated=model.build(57,'odd',data.intrinsic_tables(),57,level=7)
    observed=conic.digest_bytes(conic.canonical(regenerated))
    conic.write_json(out/'necessary_primal.exact.json',regenerated)
    same=observed==expected
    report=dict(status='EXACT_GLOBAL_CONSTRAINTS_REGENERATED' if same else 'REGENERATION_MISMATCH_REQUIRES_REVIEW',
        exactly_same_scientific_model=same,source_model_sha256=expected,regenerated_model_sha256=observed,
        source_input_sha256=before[str(frozen.relative_to(ROOT))],threshold=57,parity='odd',level=7,
        upper_dimension=57,variables=len(regenerated['variables']),equalities=len(regenerated['equations']),
        cones=len(regenerated['blocks']),all_geometric_constraints_preserved=same,
        feature_coefficients_remain_free=26,accepted_bounds_changed=False,
        historical_successful_geometry_checks_not_rerun=True,elapsed_seconds=time.monotonic()-started)
    if bindings()!=before:raise ValueError('A regeneration source changed')
    conic.write_json(out/'result.json',report);print(json.dumps(report),flush=True)
    if not same:raise ValueError('Regeneration changed the mathematical model')


if __name__=='__main__':main()
