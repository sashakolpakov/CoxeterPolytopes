"""A strictly feasible auxiliary SDP for feasibility and cone-face discovery.

Maximize e<=1 with F_b(x)-e*I PSD. Unlike a nonnegative violation
variable, this permits e>0 and has a nontrivial trace-one face dual at
optimum zero. Only exact original points or exact cone identities count.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import time

from . import conic,dual_completion,face_phase,feasibility,sdpa_cli,sdpa_double

E='_common_PSD_margin'


def augmented(problem):
    if problem['equations'] or E in problem['variables']:raise ValueError('Use an equation-free source')
    blocks=[dict(b,entries=[[conic.plus(q,conic.form(**{E:-1})) if i==j else q
        for j,q in enumerate(row)] for i,row in enumerate(b['entries'])]) for b in problem['blocks']]
    blocks.append(dict(name='margin_upper_bound',entries=[[conic.form(1,**{E:-1})]]))
    return conic.model([*problem['variables'],E],blocks)


def write_input(path,problem):
    conic.write_sdpa(path,problem);temporary=path.with_suffix('.partial')
    with path.open() as source,temporary.open('x') as target:
        for _ in range(4):target.write(source.readline())
        source.readline();target.write(' '.join('-1' if v==E else '0' for v in problem['variables'])+'\n')
        shutil.copyfileobj(source,target,length=1024**2)
    temporary.replace(path)


def strict_initial(problem,path):
    import numpy as np
    rho=F(1)
    for block in problem['blocks']:
        for i,row in enumerate(block['entries']):
            rho=max(rho,1+sum((abs(q.get('',F(0))) for j,q in enumerate(row) if j!=i),F(0))-row[i].get('',F(0)))
    full=[];inverses=[];scalar=[];trace=0.
    for block in problem['blocks']:
        exact=[[q.get('',F(0))+(rho if i==j else 0) for j,q in enumerate(row)]
            for i,row in enumerate(block['entries'])]
        for i,row in enumerate(exact):
            if row[i]-sum((abs(c) for j,c in enumerate(row) if j!=i),F(0))<1:
                raise ValueError('The exact Gershgorin start is not strictly positive')
        M=np.asarray([[float(c) for c in row] for row in exact])
        if len(M)==1:scalar.append(float(M[0,0]));trace+=1/float(M[0,0])
        else:
            inverse=np.linalg.inv(M);full.append(M);inverses.append(inverse);trace+=float(np.trace(inverse))
    scalar.append(float(1+rho));trace+=1/float(1+rho);tau=1/trace
    def number(c):return format(float(c),'+.17e')
    def matrix_text(dual):
        blocks=[]
        for M,inverse in zip(full,inverses,strict=True):
            A=tau*inverse if dual else M
            blocks.append('{'+',\n'.join('{'+', '.join(number(c) for c in row)+'}' for row in A)+'}')
        blocks.append('{'+', '.join(number(tau/c if dual else c) for c in scalar)+'}')
        return '{'+',\n'.join(blocks)+'}'
    with path.open('x') as stream:
        stream.write('xVec = {'+', '.join(['0']*len(problem['variables'])+[number(-rho)])+'}\n')
        stream.write('xMat = '+matrix_text(False)+'\nyMat = '+matrix_text(True)+'\n')
    return dict(exact_Gershgorin_rho=rho,all_original_augmented_blocks_strictly_diagonally_dominant=True,
        exact_minimum_Gershgorin_margin=1,initial_common_complementarity=tau,
        epsilon_dual_trace_equation_satisfied_numerically=True)


def sparse_exposure(problem,matrices):
    weights=[max(0.,sum(float(M[i][i]) for i in range(len(M)))) for M in matrices]
    selected=set(sorted(range(len(weights)),key=lambda b:weights[b],reverse=True)[:400])
    # Unit rows complete the adjoint without retaining all 90k coupled rows.
    for b,block in enumerate(problem['blocks']):
        if len(block['entries'])==1:
            q=block['entries'][0][0]
            if len(q)==1 and '' not in q:selected.add(b)
    indices=sorted(selected)
    restricted=conic.model(problem['variables'],[problem['blocks'][b] for b in indices])
    candidate,report=face_phase.exposure(restricted,[matrices[b] for b in indices],seconds=120)
    report.update(candidate_blocks=len(indices),all_rays_are_from_source_cones=True)
    if candidate is None:return None,report
    equations,proof=candidate
    proof=[dict(item,block=indices[item['block']]) for item in proof]
    if conic.plus(*(conic.times(item['weight'],face_phase.quadratic(
            problem['blocks'][item['block']],item['direction'])) for item in proof)):
        raise ValueError('The sparse exposure does not give an exact source identity')
    return (equations,proof),report


def self_check(root,out):
    out.mkdir()
    P=conic.model(['x'],[dict(name='PSD',entries=[
        [conic.form(1),conic.form(x=1)],[conic.form(x=1),conic.form(1)]]),
        dict(name='lower',entries=[[conic.form(-1,x=1)]]),dict(name='nonnegative_x',entries=[[conic.form(x=1)]])])
    Q=augmented(P);write_input(out/'input.dat-s',Q)
    initialization=strict_initial(P,out/'initial.dat')
    sdpa_cli.call(out/'input.dat-s',out,80,20,out/'initial.dat',fixture=True)
    raw,phase=conic.native_vector(out/'native.result')
    if abs(raw[-1])>F(1,10**6):raise ValueError('The known singular fixture has the wrong margin optimum')
    matrices=feasibility.native_cones(out/'native.result',Q,field='yMat')[:-1]
    face,checked=sparse_exposure(P,matrices)
    if face is None:raise ValueError('The margin dual did not recover the known exact face: '+str(checked))
    equations,proof=face
    if any(conic.value(q,{'x':F(1)}) for q in equations):raise ValueError('A false fixture face was exposed')
    report=dict(status='STRICT_MARGIN_FIXTURE_AND_EXACT_FACE_PASSED',phase=phase,
        initialization=initialization,forced_equations=equations,positive_zero_identity=proof,exposure_check=checked)
    conic.write_json(out/'result.json',report);print(json.dumps(conic.plain(report)),flush=True)


def run(root,out,source,seconds,input_path=None):
    raise ValueError('Retired solver path: use the verified preparation pipeline and run_sdpa_checkpoint_20261125.py')
