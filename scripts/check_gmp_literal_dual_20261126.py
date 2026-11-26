#!/usr/bin/env python3
"""Check a GMP xVec directly over Q in every exact G(u) cone.

All printed digits are preserved. A cone-feasible point is only a candidate
until the separate full adjoint verifier has replayed it.
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import re
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import balanced_dual,conic,field_rounding,solver_gate
from scripts import run_balanced_dual_gmp_20261126 as gmp
from scripts.search_balanced_dual_field_20261126 import checked_model,save_point


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('checkpoint','model','replay','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use a serial scientific watchdog')
    args.output=args.output.resolve()
    if args.output.exists() or not args.output.is_relative_to(ROOT/'results'):raise ValueError('Use fresh output')
    native=gmp.read(args.checkpoint/'result.json');path=args.checkpoint/'native.result'
    if (native['native']['returncode']!=0
            or not all(native['native']['solution_arrays_present'].values())
            or native['native']['native_result_sha256']!=solver_gate.file_digest(path)):
        raise ValueError('No complete unaltered GMP output')
    checked_model(args.model,args.replay)
    prepared=Path(native['prepared']);info=gmp.read(prepared/'result.json')
    if (Path(info['model']).resolve()!=args.model.resolve()
            or info['native_input_sha256']!=solver_gate.file_digest(prepared/'input.dat-s')):
        raise ValueError('This checkpoint concerns a different exact model')
    gmp.check_bindings(info['source_bindings'])
    meta=gmp.read(args.model/'coordinates.json')
    raw,phase=conic.native_vector(path)
    if len(raw)!=len(meta['variables']):raise ValueError('Incomplete free vector')
    field=field_rounding.Field(0);point=[field.element(v) for v in raw]
    forms=balanced_dual.read_forms(args.model/'forms.exact.jsonl.gz',meta['forms'])
    with gzip.open(args.model/'blocks.json.gz','rt') as stream:blocks=json.load(stream)
    args.output.mkdir();started=time.monotonic()
    def progress(**event):
        with (args.output/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(elapsed_seconds=time.monotonic()-started,**event))+'\n')
    checked=field_rounding.check_point(forms,blocks,point,field,progress)
    report=dict(status='GMP_LITERAL_VECTOR_NOT_FEASIBLE_IN_EXACT_G_CONES',
        checkpoint=str(args.checkpoint),native_phase=phase,
        method='All printed GMP xVec digits interpreted exactly as rationals; no binary64 conversion',
        exact_cone_check=checked,exact_certificate_obtained=False,
        original_primal_feasibility_status='UNDECIDED',separate_replay_still_required=True,
        model=str(args.model),separate_dual_parametrization_replay=str(args.replay),
        source_bindings=gmp.bind([Path(__file__),Path(gmp.__file__),
            ROOT/'research/global_formulation/field_rounding.py',
            ROOT/'scripts/search_balanced_dual_field_20261126.py',
            args.checkpoint/'result.json',path,args.model/'result.json',args.replay]))
    if checked['valid']:
        save_point(args.output/'field_point.exact.json.gz',point,field)
        report.update(status='EXACT_FIELD_DUAL_CERTIFICATE_CANDIDATE_REQUIRES_SEPARATE_REPLAY',
            field='Q',field_point_sha256=solver_gate.file_digest(args.output/'field_point.exact.json.gz'),
            all_cones_checked=checked)
    else:
        index=checked['block_index'];block=blocks[index]
        if block['size']==1:
            row=forms[block['start']]
            total=abs(row.get(-1,field.Q(0)))
            for col,c in row.items():
                if col>=0:total+=abs(c*point[col][0])
            value=field.Q(checked['value'][0])
            report['first_scalar_diagnostic']=dict(value_exact=str(value),
                absolute_term_sum_exact=str(total),
                relative_to_absolute_term_sum=None if not total else float(value/total))
    text=path.read_text()
    report['native_reported_metrics']={k:v for k,v in re.findall(
        r'^\s*(mu|objValPrimal|objValDual|p\.feas\.error|d\.feas\.error|relative eps|gap)\s*=\s*(\S+)',text,re.M)}
    conic.write_json(args.output/'result.json',report)
    print(json.dumps(dict(status=report['status'],cone_check=checked,output=str(args.output))),flush=True)


if __name__=='__main__':main()
