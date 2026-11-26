"""Serial native SDPA, with warm starts and checks in the frozen original."""
from fractions import Fraction as F
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from . import conic,feasibility,sdpa_double,solver_gate

EXECUTABLE=Path.home()/'.local/opt/coxeter-sdpa-7.3.21/bin/sdpa'
GMP_EXECUTABLE=Path.home()/'.local/opt/sdpa-gmp/bin/sdpa_gmp'


def parameters(path,iterations,precision=None):
    accuracy='1e-24' if precision else '1e-10'
    values=[str(iterations),accuracy,'10','10000','-1e12','1e12',
        '0.1','0.3','0.9',accuracy]
    if precision:values.append(str(precision))
    # Parameter::readFile uses "%s %*[^\n]" for these four lines.  With
    # bare tokens its whitespace directive crosses the newline and the
    # scanset discards the NEXT format.  Keep a nonempty label on each
    # line, as in the native parameter-file examples.
    output_format='%+90.80Fe' if precision else '%+26.17e'
    values += [output_format+' = '+name for name in ('xPrint','XPrint','YPrint','infPrint')]
    path.write_text('\n'.join(values)+'\n')


def call(input_path,out,iterations,seconds,initial=None,precision=None,*,preflight=None,fixture=False):
    checked=solver_gate.verify_native(input_path,preflight,fixture=fixture)
    executable=GMP_EXECUTABLE if precision else EXECUTABLE
    parameters(out/'param.sdpa',iterations,precision)
    command=[str(executable),'-ds',str(input_path),'-o',str(out/'native.result'),
        '-p',str(out/'param.sdpa')]
    initial_source=initial
    if initial is not None:
        # IO::readInit scans numbers from the beginning of its input. Full
        # solver reports contain dates, parameters and iteration tables
        # BEFORE xVec; passing the report directly silently corrupts state.
        dense=out/('gmp_initial.dat' if precision else 'native_initial.dat')
        dense.write_text('\n'.join(feasibility.native_array(initial,field) for field in ('xVec','xMat','yMat'))+'\n')
        initial=dense
    if not precision:command+=['-numThreads','1']
    if initial is not None:command+=['-id',str(initial)]
    metadata=dict(command=command,input_sha256=checked['native_input_sha256'],preflight_check=checked,
        executable_sha256=conic.digest_bytes(executable.read_bytes()),precision_bits=precision,
        initial_path=str(initial) if initial else None,
        initial_sha256=conic.digest_bytes(initial.read_bytes()) if initial else None,
        initial_source_path=str(initial_source) if initial_source else None,
        initial_source_sha256=conic.digest_bytes(initial_source.read_bytes()) if initial_source else None,
        restart_uses_only_the_three_solution_arrays=initial is not None,
        parameters_sha256=conic.digest_bytes((out/'param.sdpa').read_bytes()),one_thread=True)
    conic.write_json(out/'native_metadata.json',metadata)
    timed_out=False
    with (out/'solver.stdout').open('x') as stream:
        try:code=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=seconds,check=False).returncode
        except subprocess.TimeoutExpired:code=None;timed_out=True
    actual_formats={}
    if not timed_out and code==0 and (out/'native.result').exists():
        with (out/'native.result').open() as stream:header=stream.read(32768)
        actual_formats=dict(re.findall(r'^\s*(xPrint|XPrint|YPrint|infPrint)\s*=\s*(\S+)',header,re.M))
        expected='%+90.80Fe' if precision else '%+26.17e'
        if actual_formats!={name:expected for name in ('xPrint','XPrint','YPrint','infPrint')}:
            raise ValueError('SDPA did not apply all four requested solution formats: '+str(actual_formats))
    return dict(timed_out=timed_out,returncode=code,actual_output_formats=actual_formats,**metadata)


def self_check(root,out,precision=None):
    out.mkdir()
    P=conic.model(['x'],[
        dict(name='PSD',entries=[[conic.form(2),conic.form(x=1)],[conic.form(x=1),conic.form(2)]]),
        dict(name='lower',entries=[[conic.form(-1,x=1)]]),
        dict(name='upper',entries=[[conic.form(1,x=-1)]])])
    conic.write_sdpa(out/'input.dat-s',P)
    first=out/'first';first.mkdir();call(out/'input.dat-s',first,4,20,precision=precision,fixture=True)
    resumed=out/'resumed';resumed.mkdir();call(out/'input.dat-s',resumed,60,20,first/'native.result',precision,fixture=True)
    raw,phase=conic.native_vector(resumed/'native.result')
    if abs(raw[0]-1)>F(1,10**7):raise ValueError('The native warm start changed the free-vector sign')
    matrices=feasibility.native_cones(resumed/'native.result',P,field='xMat')
    for block,M in zip(P['blocks'],matrices,strict=True):
        for i,row in enumerate(block['entries']):
            for j,q in enumerate(row):
                if abs(float(M[i][j]-conic.value(q,{'x':raw[0]})))>1e-6:
                    raise ValueError('The native cone order or off-diagonal convention changed')
    import re
    first_text=(first/'solver.stdout').read_text();second_text=(resumed/'solver.stdout').read_text()
    first_mu=re.search(r'^\s*0\s+(\S+)',first_text,re.M)
    second_mu=re.search(r'^\s*0\s+(\S+)',second_text,re.M)
    if first_mu is None or second_mu is None or abs(float(second_mu[1])-float(first_mu[1]))<1:
        raise ValueError('The second native call did not use its saved initial point')
    first_end=re.search(r'^\s*mu\s*=\s*(\S+)',(first/'native.result').read_text(),re.M)
    if first_end is None or abs(float(second_mu[1])-float(first_end[1]))>0.06*abs(float(first_end[1])):
        raise ValueError('The restart did not preserve the saved complementarity')
    central=out/'central';central.mkdir()
    conic.write_json(central/'reference.json',{'x':F(1)})
    initialization=proposed_initial(P,central,central/'reference.json',central/'initial.dat',1e-3,1e-3)
    if initialization['maximum_relative_centrality_error']>1e-8:raise ValueError('The proposed central seed failed its matrix identity')
    call(out/'input.dat-s',central,60,20,central/'initial.dat',precision,fixture=True)
    central_raw,_=conic.native_vector(central/'native.result')
    if abs(central_raw[0]-1)>F(1,10**7):raise ValueError('The centrally initialized fixture failed its known point')
    report=dict(status='NATIVE_SDPA_WARM_START_AND_ENCODING_CHECK_PASSED',phase=phase,
        free_coordinate_sign=True,original_scalar_PSD_order=True,offdiagonal_convention=True,
        saved_initial_point_used=True,central_inverse_slack_seed_checked=True,
        restart_complementarity_matches_saved_terminal_value=True,
        log_header_numbers_excluded_from_initial_state=True,
        all_four_native_solution_output_formats_verified=True,
        precision_bits=precision,
        executable_sha256=conic.digest_bytes((GMP_EXECUTABLE if precision else EXECUTABLE).read_bytes()))
    conic.write_json(out/'result.json',report);print(json.dumps(report),flush=True)


def proposed_initial(problem,source,reference_path,path,epsilon=1e-4,tau=1e-12):
    import numpy as np
    supplied={v:F(c) for v,c in json.loads(reference_path.read_text()).items()}
    if (source/'equilibration.exact.json').exists():
        units=json.loads((source/'equilibration.exact.json').read_text())
        intervals=json.loads((source/'interval_conditioning.exact.json').read_text())
    elif (source/'conditioning.exact.json').exists():
        combined=json.loads((source/'conditioning.exact.json').read_text())
        units=combined['units'];intervals=combined['intervals']
    else:units=dict(column_scales={v:'1' for v in problem['variables']});intervals=dict(pairs=[])
    columns={v:F(c) for v,c in units['column_scales'].items()}
    for pair in intervals['pairs']:
        width={v:F(c) for v,c in pair['width'].items()}
        supplied[pair['coordinate']]=abs(conic.value(width,supplied))/(2*F(pair['delta']))
    values={v:supplied[v]/columns[v] for v in problem['variables']}
    full=[];dual_full=[];scalar=[];maximum_repair=0.;centrality=0.
    for block in problem['blocks']:
        M=np.asarray([[float(conic.value(q,values)) for q in row] for row in block['entries']])
        scale=max(1.,float(np.max(np.abs(M))));floor=epsilon*scale
        if len(M)==1:
            scalar.append(max(float(M[0,0]),floor));maximum_repair=max(maximum_repair,floor-float(M[0,0]))
            continue
        eigen,U=np.linalg.eigh((M+M.T)/2)
        maximum_repair=max(maximum_repair,float(np.max(np.maximum(floor-eigen,0))))
        clipped=np.maximum(eigen,floor)
        M=(U*clipped[None,:])@U.T
        dual=(U*(tau/clipped)[None,:])@U.T
        centrality=max(centrality,float(np.max(np.abs(M@dual-tau*np.eye(len(M)))))/tau)
        full.append(M);dual_full.append(dual)
    def number(value):return format(float(value),'+.17e')
    def matrices(dual):
        blocks=[]
        for M,Adual in zip(full,dual_full,strict=True):
            A=Adual if dual else M
            blocks.append('{'+',\n'.join('{'+', '.join(number(v) for v in row)+'}' for row in A)+'}')
        if scalar:blocks.append('{'+', '.join(number(tau/v if dual else v) for v in scalar)+'}')
        return '{'+',\n'.join(blocks)+'}'
    with path.open('x') as stream:
        stream.write('xVec = {'+', '.join(number(values[v]) for v in problem['variables'])+'}\n')
        stream.write('xMat = '+matrices(False)+'\nyMat = '+matrices(True)+'\n')
    return dict(kind='candidate_only_strict_cone_initialization',reference=str(reference_path),
        projected_initial_cones_are_not_original_feasibility=True,minimum_relative_primal_eigenvalue=epsilon,
        initial_common_complementarity=tau,dual_initialization='tau times inverse proposed primal slack',
        maximum_relative_centrality_error=centrality,free_coordinates_not_clipped=True,
        maximum_cone_eigenvalue_repair=maximum_repair)


def run(root,out,source,reference_path,iterations,seconds,resume=None,epsilon=1e-3,tau=1e-3,precision=None):
    raise ValueError('Retired solver path: use the verified preparation pipeline and run_sdpa_checkpoint_20261125.py')
