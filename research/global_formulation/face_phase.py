"""An auxiliary SDPA face search with exact exposure and original checks.

The auxiliary trace-one problem is not the original feasibility problem.
A numerical answer is useful only after an original-point check or an
exact nonnegative zero identity exposes a proper cone face.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import time

from . import conic, feasibility, sdpa_double


def face_input(path,problem):
    variables={v:i+1 for i,v in enumerate(problem["variables"])}
    const=len(variables)+1;trace=const+1
    full=[b for b in problem["blocks"] if len(b["entries"])>1]
    scalar=[b for b in problem["blocks"] if len(b["entries"])==1]
    sizes=[len(b["entries"]) for b in full]+([-len(scalar)] if scalar else [])
    with path.open("x") as stream:
        stream.write('"Auxiliary trace-one face problem; all acceptances checked exactly\n')
        stream.write(f'{trace}\n{len(sizes)}\n'+" ".join(map(str,sizes))+"\n")
        stream.write(" ".join(["0"]*const+["-1"])+"\n")
        def entry(block,i,j,q):
            for v,c in q.items():
                if c:stream.write(f'{variables[v] if v else const} {block} {i} {j} {conic.decimal_text(-c)}\n')
            if i==j:stream.write(f'{trace} {block} {i} {j} -1\n')
        for b,block in enumerate(full,1):
            for i,row in enumerate(block["entries"]):
                for j in range(i,len(row)):entry(b,i+1,j+1,row[j])
        for i,block in enumerate(scalar,1):entry(len(full)+1,i,i,block["entries"][0][0])


def quadratic(block,direction):
    return conic.plus(*(conic.times(a*b,block["entries"][i][j])
        for i,a in direction for j,b in direction))


def self_check(out):
    out.mkdir()
    P=conic.model(["x"],[dict(name="PSD",entries=[
        [conic.form(1),conic.form(x=1)],[conic.form(x=1),conic.form(1)]]),
        dict(name="lower",entries=[[conic.form(-1,x=1)]])])
    face_input(out/"input.dat-s",P)
    sdpa_double.native_call(out/"input.dat-s",out,80,.1,1e-9,fixture=True)
    values,phase=conic.native_vector(out/"native.result")
    matrices=feasibility.native_cones(out/"native.result",P,field="yMat")
    candidate,report=exposure(P,matrices,20)
    if candidate is None:raise ValueError("No exact face recovered in the known singular fixture: "+str(report))
    equations,proof=candidate
    if any(conic.value(q,{"x":F(1)}) for q in equations):raise ValueError("A false kernel equation was exposed")
    if not any(conic.value(q,{"x":F(0)}) for q in equations):raise ValueError("The fixture face is not proper")
    conic.write_json(out/"result.json",dict(status="EXACT_FACE_FIXTURE_PASSED",phase=phase,
        checked_nonnegative_identity=proof,forced_equations=equations,LP_report=report))
    print(json.dumps(dict(status="EXACT_FACE_FIXTURE_PASSED",forced_kernel_rows=len(equations))),flush=True)


def exposure(problem,matrices,seconds=120):
    """LP proposes positive combinations; a rational identity gates use."""
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix, vstack
    rows=[];keys=set()
    def add(b,d):
        q=quadratic(problem["blocks"][b],d)
        if not q:return
        scale=max(abs(c) for c in q.values());q=conic.times(1/scale,q)
        key=conic.canonical(q)
        if key in keys:return
        keys.add(key);rows.append(dict(block=b,direction=d,row=q,scale=scale))
    for b,(block,M) in enumerate(zip(problem["blocks"],matrices,strict=True)):
        r=len(M)
        if r==1:add(b,[(0,F(1))]);continue
        for i in range(r):add(b,[(i,F(1))])
        values,vectors=np.linalg.eigh(np.asarray([[float(c) for c in row] for row in M]))
        cutoff=max(1e-9,float(values[-1])*1e-5)
        for k,value in enumerate(values):
            if value<cutoff:continue
            v=vectors[:,k];v=v/max(abs(v))
            for denominator in (1,10,1000):
                d=[(i,F(str(float(c))).limit_denominator(denominator))
                    for i,c in enumerate(v) if abs(c)>1e-7]
                d=[(i,c) for i,c in d if c]
                if d:add(b,d)
    names=sorted({v for item in rows for v in item["row"]});index={v:i for i,v in enumerate(names)}
    ii=[];jj=[];vv=[]
    for j,item in enumerate(rows):
        for v,c in item["row"].items():ii.append(index[v]);jj.append(j);vv.append(float(c))
    B=coo_matrix((vv,(ii,jj)),shape=(len(names),len(rows))).tocsc()
    normalization=coo_matrix(np.ones((1,len(rows))))
    # Prefer a PSD direction over scalar-only exposures, with bounded units.
    objective=np.array([-1. if len(problem["blocks"][item["block"]]["entries"])>1 else 0. for item in rows])
    result=linprog(objective,A_eq=vstack([B,normalization]).tocsc(),
        b_eq=np.r_[np.zeros(len(names)),1.],bounds=(0,None),method="highs",
        options=dict(time_limit=seconds,threads=1,primal_feasibility_tolerance=1e-9,
            dual_feasibility_tolerance=1e-9))
    report=dict(status=int(result.status),message=str(result.message),candidate_quadratic_rows=len(rows))
    if result.x is None:return None,report
    support=[j for j,c in enumerate(result.x) if c>1e-8]
    report["proposed_support"]=len(support)
    exact=None
    for limit in (100,10000,10**8):
        proposed={j:F(str(float(result.x[j]))).limit_denominator(limit) for j in support}
        if proposed and all(c>=0 for c in proposed.values()) and not conic.plus(*(
                conic.times(c,rows[j]["row"]) for j,c in proposed.items())):
            exact=proposed;break
    if exact is None and 0<len(support)<=3500:
        variables={j:"weight_"+str(j) for j in support};equations={}
        for j in support:
            for v,c in rows[j]["row"].items():equations.setdefault(v,{})[variables[j]]=c
        normalization={v:F(1) for v in variables.values()};normalization[""]=-F(1)
        P=conic.model(list(variables.values()),[dict(name=v,entries=[[conic.form(**{v:1})]])
            for v in variables.values()],list(equations.values())+[normalization])
        try:
            Q,subs=conic.eliminate(P)
            numerical={variables[j]:F(str(float(result.x[j]))) for j in support}
            point,_=conic.reconstruct(P,Q,subs,[numerical[v] for v in Q["variables"]])
            exact={j:point[v] for j,v in variables.items() if point[v]>0}
            report["exact_nullspace_projection"]=True
        except ValueError as error:report["exact_projection_failure"]=str(error)
    if exact is None:
        report["exact_identity_verified"]=False;return None,report
    combination=[dict(block=rows[j]["block"],direction=rows[j]["direction"],weight=c/rows[j]["scale"])
        for j,c in exact.items() if c>0]
    if not combination or conic.plus(*(conic.times(item["weight"],quadratic(
            problem["blocks"][item["block"]],item["direction"])) for item in combination)):
        raise ValueError("The cone exposure failed its exact identity")
    equations={}
    for item in combination:
        for row in problem["blocks"][item["block"]]["entries"]:
            q=conic.plus(*(conic.times(c,row[i]) for i,c in item["direction"]))
            if q:equations[conic.canonical(q)]=q
    if not equations:
        report["exact_identity_verified_but_no_new_kernel_rows"]=True;return None,report
    report.update(exact_identity_verified=True,positive_summands=len(combination),forced_kernel_rows=len(equations))
    return (list(equations.values()),combination),report


def run(root,out,source,prepared,seconds=550):
    raise ValueError('Retired solver path: use the verified preparation pipeline and run_sdpa_checkpoint_20261125.py')
