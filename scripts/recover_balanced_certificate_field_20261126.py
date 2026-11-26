#!/usr/bin/env python3
"""Reconstruct the opposite certificate if the explicit dual is infeasible.

This searches the original normalized certificate side, using its exact
homogeneous parametrization. It is complementary to the explicit-dual search.
"""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, field_rounding as rounding, solver_gate
from scripts.solve_balanced_certificate_20261126 import verified_balanced_model
from scripts.search_balanced_dual_field_20261126 import checked_model, proposals, save_point


def homogeneous_forms(rational,field):
    import numpy as np
    with gzip.open(rational/'substitutions.exact.json.gz','rt') as stream:subs=json.load(stream)
    free=np.load(rational/'free.npz',allow_pickle=False)['original_dual_free_coordinates']
    inverse={int(col):i for i,col in enumerate(free)}
    n=len(free)+len(subs);forms=[]
    for col in range(n):
        if col in inverse:forms.append({inverse[col]:field.Q(1)})
        else:
            q={}
            for old,c in subs[str(col)].items():
                if int(old)<0:raise ValueError('The certificate cone is not homogeneous')
                q[inverse[int(old)]]=field.Q(c)
            forms.append(q)
    return forms,free


def normalize(point,field,units,free):
    inverse={int(col):i for i,col in enumerate(free)}
    row={int(col):field.Q(c) for col,c in units['source_normalization'].items()}
    largest=max(abs(c) for col,c in row.items() if col>=0)
    value=field.zero
    for col,c in row.items():
        if col>=0:value=field.add(value,field.scale(c/largest,point[inverse[col]]))
    if field.sign(value)>=0:raise ValueError('The homogeneous point has no negative normalization')
    factor=field.div(field.neg(field.one),value)
    if field.sign(factor)<=0:raise ValueError('Normalization must preserve the PSD cones')
    result=[field.mul(factor,v) for v in point]
    checked=field.zero
    for col,c in row.items():
        if col>=0:checked=field.add(checked,field.scale(c/largest,result[inverse[col]]))
    if checked!=field.neg(field.one):raise ValueError('Exact normalization failed')
    return result,dict(positive_normalization_factor=field.text(factor),exact_normalization=-1)


def source_adjoint(forms,point,field,exact_model,face,encoding,tau):
    """Check the full prepared source Farkas identities, including zero faces."""
    import numpy as np
    from flint import fmpq
    model=conic.read_model(exact_model)
    layout=json.loads((encoding/'layout.json').read_text());scalar_total=layout['scalar_cones']-1
    embedding=np.load(face/'embedding.npz',allow_pickle=False)['dual_coordinate_embedding']
    inverse={int(old):i for i,old in enumerate(embedding)}
    position={v:i for i,v in enumerate(model['variables'])}
    a=[fmpq(0) for _ in model['variables']];b=[fmpq(0) for _ in model['variables']]
    constant=field.zero;scalar=0;base=scalar_total
    for block in model['blocks']:
        M=block['entries'];n=len(M)
        if n==1:locations=[(scalar,0,0)];scalar+=1
        else:
            locations=[]
            for i in range(n):
                for j in range(i,n):locations.append((base,i,j));base+=1
        for old,i,j in locations:
            if old not in inverse:continue
            value=field.scale(tau*(1 if i==j else 2),rounding.evaluate(forms[inverse[old]],point,field))
            if value==field.zero:continue
            for name,c in M[i][j].items():
                c=fmpq(str(c))
                if not name:constant=field.add(constant,field.scale(c,value))
                else:
                    k=position[name];a[k]+=c*value[0]
                    if value[1]:b[k]+=c*value[1]
    if constant!=field.neg(field.one):raise ValueError('Nonnegative source-cone combination has wrong exact constant')
    if any(x or y for x,y in zip(a,b,strict=True)):raise ValueError('A full prepared-source Farkas identity is nonzero')
    return dict(full_prepared_source_adjoint_equations_checked=len(a),exact_constant_pairing=-1,
        removed_dual_rows_reinserted_as_zero=True,prepared_necessary_data_model_infeasible=True)


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','attempt','verify'])
    parser.add_argument('--balanced',type=Path)
    parser.add_argument('--rational',type=Path)
    parser.add_argument('--dual-model',type=Path)
    parser.add_argument('--dual-verified',type=Path)
    parser.add_argument('--face',type=Path)
    parser.add_argument('--encoding',type=Path)
    parser.add_argument('--exact-model',type=Path)
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--solver-output',type=Path)
    parser.add_argument('--source-role',choices=['balanced-primal','explicit-dual-ray'],default='explicit-dual-ray')
    parser.add_argument('--candidate',type=Path)
    parser.add_argument('--fields',default='0,2,5')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    paths=[Path(__file__),ROOT/'research/global_formulation/field_rounding.py',
        ROOT/'scripts/search_balanced_dual_field_20261126.py',ROOT/'scripts/solve_balanced_certificate_20261126.py',
        ROOT/'research/global_formulation/conic.py',ROOT/'research/global_formulation/conic_backend.py']
    if args.action!='check':
        paths += [args.balanced/'result.json',args.balanced/'units.exact.json',
            args.rational/'substitutions.exact.json.gz',args.rational/'free.npz',
            args.dual_model/'blocks.json.gz',args.dual_verified]
        if args.action=='attempt':paths += [args.solver_output/'result.json',args.solver_output/'solution.npz']
        if args.action=='verify':paths += [args.candidate/'result.json',args.candidate/'field_point.exact.json.gz',
            args.exact_model,args.face/'embedding.npz',args.encoding/'layout.json']
    def bindings():return {str(p):dict(bytes=p.stat().st_size,sha256=solver_gate.file_digest(p)) for p in paths}
    before=bindings();args.output.mkdir();conic.write_json(args.output/'source_bindings.json',before)
    def progress(**event):
        event['elapsed_seconds']=time.monotonic()-started
        with (args.output/'progress.jsonl').open('a') as stream:stream.write(json.dumps(event)+'\n')
    if args.action=='check':
        f=rounding.Field(2);point=[f.element(2),f.element(0,2),f.element(4)]
        units=dict(source_normalization={'-1':'1','0':'-1'})
        result,receipt=normalize(point,f,units,np.arange(3))
        if result!=[f.one,f.element(0,1),f.element(2)]:raise ValueError('Field normalization changed the ray')
        forms=[{0:f.Q(1)},{1:f.Q(1)},{2:f.Q(1)}]
        if not rounding.check_point(forms,[dict(name='rank_one',start=0,size=2)],result,f)['valid']:
            raise ValueError('Positive field normalization lost PSD')
        try:normalize(point,f,dict(source_normalization={'-1':'1','0':'1'}),np.arange(3))
        except ValueError:pass
        else:raise ValueError('A negative cone scaling was accepted')
        # A genuinely nonrational Farkas witness for x>=0, -1-x>=0, 1>=0.
        source=args.output/'fixture_source.json'
        conic.write_json(source,conic.model(['x'],[
            dict(name='positive_x',entries=[[{'x':F(1)}]]),
            dict(name='negative_x',entries=[[{'':F(-1),'x':F(-1)}]]),
            dict(name='constant',entries=[[{'':F(1)}]])]))
        face=args.output/'fixture_face';face.mkdir()
        encoding=args.output/'fixture_encoding';encoding.mkdir()
        np.savez_compressed(face/'embedding.npz',dual_coordinate_embedding=np.arange(3))
        conic.write_json(encoding/'layout.json',dict(scalar_cones=4,PSD_sizes=[]))
        q=[{0:f.Q(1)},{0:f.Q(1)},{1:f.Q(1)}]
        witness=[f.element(0,1),f.element(-1,1)]
        audit=source_adjoint(q,witness,f,source,face,encoding,f.Q(1))
        if audit['exact_constant_pairing']!=-1:raise ValueError('Wrong exact field Farkas constant')
        try:source_adjoint(q,[witness[0],witness[0]],f,source,face,encoding,f.Q(1))
        except ValueError:pass
        else:raise ValueError('A wrong field certificate normalization was accepted')
        report=dict(status='HOMOGENEOUS_FIELD_CERTIFICATE_NORMALIZATION_FIXTURES_VERIFIED',
            nonrational_Farkas_source_identity_checked=True,tampered_normalization_rejected=True)
    else:
        data,layout,_,balanced=verified_balanced_model(args.balanced,args.native_input,args.preflight)
        checked_model(args.dual_model,args.dual_verified)
        units=balanced['exact_units'];field=rounding.Field(0)
        forms,free=homogeneous_forms(args.rational,field)
        with gzip.open(args.dual_model/'blocks.json.gz','rt') as stream:blocks=json.load(stream)
        if args.action=='verify':
            recorded=json.loads((args.candidate/'result.json').read_text());path=args.candidate/'field_point.exact.json.gz'
            if solver_gate.file_digest(path)!=recorded['field_point_sha256']:raise ValueError('Changed field candidate')
            with gzip.open(path,'rt') as stream:raw=json.load(stream)
            field=rounding.Field(raw['radicand']);point=[field.element(a,b) for a,b in raw['values']]
            point,norm=normalize(point,field,units,free)
            test=rounding.check_point(forms,blocks,point,field,progress)
            if not test['valid']:raise ValueError('Exact certificate cone check failed: '+str(test))
            identities=source_adjoint(forms,point,field,args.exact_model,args.face,args.encoding,
                field.Q(units['global_free_variable_scale']))
            report=dict(status='EXACT_FIELD_INFEASIBILITY_CERTIFICATE_FOR_PREPARED_SOURCE_VERIFIED',
                field=raw['radicand'],cone_check=test,source_identity_check=identities,exact_certificate_obtained=True)
        else:
            arrays=np.load(args.solver_output/'solution.npz',allow_pickle=False)
            native=json.loads((args.solver_output/'result.json').read_text())
            if args.source_role=='balanced-primal':approximate=arrays['x']
            else:
                values=arrays['y'].copy();base=layout['scalar_cones']
                for n in layout['PSD_sizes']:
                    i,j=np.triu_indices(n);values[base:base+len(i)][i!=j]/=2**.5;base+=len(i)
                approximate=values[free]
            if len(approximate)!=len(free) or not np.isfinite(approximate).all():raise ValueError('Incomplete field seed')
            cone_data=dict(A=data['A'][1:].tocsc(),b=data['b'][1:].copy(),c=data['c'])
            report=dict(status='SOURCE_FIELD_CERTIFICATE_ATTEMPTS_NO_EXACT_CERTIFICATE',
                source_role=args.source_role,native_status=native['native']['native_status'],attempts=[],exact_certificate_obtained=False)
            accepted=None
            for denominator in [2**32,2**44,2**50]:
                field=rounding.Field(0);point=rounding.rational_point(approximate,field,denominator)
                try:
                    point,norm=normalize(point,field,units,free)
                    test=rounding.check_point(forms,blocks,point,field,progress)
                except ValueError as error:test=dict(valid=False,reason=str(error))
                report['attempts'].append(dict(method='common_grid_and_positive_normalization',grid=denominator,check=test))
                if test['valid']:accepted=(field,point,test);break
            if accepted is None:
                for d in map(int,args.fields.split(',')):
                    field=rounding.Field(d)
                    for tolerance in [1e-10,1e-12]:
                        eqs,kernels,scalar,limited=proposals(forms,blocks,cone_data,approximate,field,tolerance,progress)
                        attempt=dict(method='field_kernel_projection',radicand=d,tolerance=tolerance,kernel_blocks=len(kernels),
                            scalar_equations=scalar,candidate_term_limit_reached=limited)
                        try:
                            point,info=rounding.project(eqs,approximate,field,denominator=2**48)
                            point,norm=normalize(point,field,units,free)
                            test=rounding.check_point(forms,blocks,point,field,progress)
                            attempt.update(projection=info,normalization=norm,check=test)
                        except ValueError as error:test={'valid':False};attempt['candidate_failure']=str(error)
                        report['attempts'].append(attempt);progress(stage='source_field_attempt',**attempt)
                        if test['valid']:accepted=(field,point,test);break
                    if accepted is not None:break
            if accepted is not None:
                field,point,test=accepted;save_point(args.output/'field_point.exact.json.gz',point,field)
                report.update(status='EXACT_SOURCE_FIELD_CERTIFICATE_CANDIDATE_REQUIRES_SEPARATE_REPLAY',
                    radicand=field.d,field_point_sha256=solver_gate.file_digest(args.output/'field_point.exact.json.gz'),
                    all_cones_checked=test,separate_replay_still_required=True)
    if before!=bindings():raise ValueError('A field input or source changed during reconstruction')
    report.update(source_bindings=before,elapsed_seconds=time.monotonic()-started,
        accepted_bounds_changed=False,geometric_constraints_changed=False,original_feasibility_status='UNDECIDED')
    conic.write_json(args.output/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in {'source_bindings','attempts'}}),flush=True)


if __name__=='__main__':main()
