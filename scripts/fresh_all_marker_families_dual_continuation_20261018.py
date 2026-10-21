#!/usr/bin/env python3
"""Equivalent normalized-dual numerical continuation on the identical matrix.

The first primal-form simplex run and its iterate remain immutable. This second
method uses IPM without crossover; exact recovery and full-matrix checks, rather
than a numerical status, certify every retained bound. No geometry is repeated.
"""
import argparse,json,time
from fractions import Fraction as F
from pathlib import Path
import fresh_all_marker_families_full_model_20261018 as A

ROOT=A.ROOT;P=A.P;M=A.M
TAG='fresh_all_marker_families_dual_continuation_20261018'
OUT=ROOT/f'results/{TAG}.json'


def fixed_dual(d):
    old=M.read(A.B.TAG)
    prior=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    r=M.read(A.INPUTS[-1][0])['exported']['repair']
    assert r['face_multipliers']==prior['face_multipliers']
    alpha={f:F(x)for f,x in r['face_multipliers']};mu={v:F(x)for v,x in r['vertex_multipliers']}
    numerator=sum(d['faces'][f]['required']*x for f,x in alpha.items());denominator=sum(mu.values())
    assert numerator/denominator==F(r['bound'])
    c=dict(kind='EXACT_PARTIAL_MODEL_DUAL',face_multipliers=r['face_multipliers'],
        vertex_multipliers=r['vertex_multipliers'],objective_numerator=str(numerator),
        normalizer=str(denominator),bound=str(numerator/denominator),
        excludes_target12=numerator>12*denominator,
        recovery='Frozen raw112 fixed-face-mass repair, verified against every column of this exact model.')
    P.verify(d,c);return c


def optimize(d):
    import highspy,numpy as np,scipy
    from scipy.sparse import coo_matrix
    nf=len(d['faces']);nb=len(d['budgets']);active=d['active'];index={j:i for i,j in enumerate(active)}
    rr=[];cc=[];vv=[];inactive=0
    for f in d['faces']:
        for j,a in f['coefficients']:
            if j in index:rr.append(index[j]);cc.append(f['index']);vv.append(a)
            else:assert a<=0;inactive+=1
    for b in d['budgets']:
        for j,a in b['coefficients']:
            assert a>0
            if j in index:rr.append(index[j]);cc.append(nf+b['index']);vv.append(-a)
        rr.append(len(active));cc.append(nf+b['index']);vv.append(1)
    matrix=coo_matrix((np.asarray(vv,dtype=float),(rr,cc)),shape=(len(active)+1,nf+nb)).tocsc()
    lp=highspy.HighsLp();lp.num_col_=nf+nb;lp.num_row_=len(active)+1
    lp.col_cost_=np.r_[-np.asarray([f['required']for f in d['faces']],dtype=float),np.zeros(nb)]
    lp.col_lower_=np.zeros(nf+nb);lp.col_upper_=np.full(nf+nb,highspy.kHighsInf)
    lp.row_lower_=np.full(len(active)+1,-highspy.kHighsInf);lp.row_upper_=np.r_[np.zeros(len(active)),1.]
    lp.a_matrix_.format_=highspy.MatrixFormat.kColwise;lp.a_matrix_.start_=matrix.indptr
    lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
    h=highspy.Highs()
    options=dict(solver='ipm',run_crossover='off',time_limit=900.,threads=4,presolve='on',random_seed=0)
    for k,v in options.items():assert h.setOptionValue(k,v)==highspy.HighsStatus.kOk
    assert h.passModel(lp)==highspy.HighsStatus.kOk
    print('EQUIVALENT_NORMALIZED_DUAL_IPM_BEGIN',matrix.shape,matrix.nnz,'inactive_nonpositive_entries',inactive,flush=True)
    start=time.monotonic();runstatus=h.run();elapsed=time.monotonic()-start
    sol=h.getSolution();info=h.getInfo();status=h.getModelStatus()
    solver=dict(interface='highspy',highs_version=h.version(),scipy_version=scipy.__version__,numpy_version=np.__version__,
        method='IPM on equivalent normalized dual',options=options,run_status=str(runstatus),model_status=str(status),
        model_status_string=h.modelStatusToString(status),success=status==highspy.HighsModelStatus.kOptimal,
        elapsed_seconds=elapsed,ipm_iterations=info.ipm_iteration_count,simplex_iterations=info.simplex_iteration_count,
        objective=info.objective_function_value,primal_solution_status=info.primal_solution_status,
        dual_solution_status=info.dual_solution_status,value_valid=sol.value_valid,dual_valid=sol.dual_valid,
        numerical_warm_start=False,feasible_exact_start='Frozen raw112 normalized dual remains an exact feasible candidate; IPM is started without a basis so presolve stays available.')
    path=ROOT/f'results/{TAG}_numerical.json';assert not path.exists()
    path.write_text(json.dumps(dict(solver=solver,col_value=list(sol.col_value),row_dual=list(sol.row_dual)))+'\n')
    print('EQUIVALENT_NORMALIZED_DUAL_IPM_END',solver,flush=True)
    candidates=[];errors=[]
    if sol.value_valid:
        try:candidates.append(P.recover_dual(d,[-float(x)for x in sol.col_value]))
        except (AssertionError,ValueError,OverflowError)as error:errors.append(['dual',repr(error)])
    if sol.dual_valid:
        try:candidates.append(P.ENGINE.recover_primal(d,[-float(x)for x in sol.row_dual]))
        except (AssertionError,ValueError,OverflowError)as error:errors.append(['primal',repr(error)])
    for c in candidates:P.verify(d,c)
    solver['exact_recovery_errors']=errors
    return solver,candidates,str(path.relative_to(ROOT)),A.M.digest(str(path.relative_to(ROOT)))


def run(verify):
    first=M.read(A.TAG)
    assert first['solver_invocations']==1 and first['solver']['model_status_string']=='Time limit reached'
    assert M.digest(first['numerical_snapshot'])==first['numerical_snapshot_sha256']
    d,prepared=A.build();assert prepared==first['prepared']
    sources={f'scripts/{TAG}.py':M.digest(f'scripts/{TAG}.py'),f'results/{A.TAG}.json':M.digest(f'results/{A.TAG}.json'),
        first['numerical_snapshot']:first['numerical_snapshot_sha256']}
    if verify:
        saved=M.read(TAG);assert saved['prepared']==prepared and saved['sources']==sources
        assert saved['solver_invocations']==2 and saved['solver_history'][0]==first['solver']
        assert M.digest(saved['numerical_snapshot'])==saved['numerical_snapshot_sha256']
        for c in saved['certificates']:P.verify(d,c)
        assert F(saved['exact_interval']['lower'])<=F(saved['exact_interval']['upper'])
        print('PASS_ALL_MARKER_FAMILIES_DUAL_CONTINUATION_FULL_MATRIX_REPLAY',saved['exact_interval'],flush=True);return
    assert not OUT.exists()
    out=dict(status='EQUIVALENT_DUAL_CONTINUATION_STARTED',prepared=prepared,sources=sources,solver_invocations=2,
        solver_history=[first['solver']],certificates=[])
    OUT.write_text(P.BASE.dump(out)+'\n')
    solver,candidates,path,h=optimize(d)
    candidates+=first['certificates']+[fixed_dual(d)]
    oldprimal=next(c for c in M.read(A.B.TAG)['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    candidates.append(oldprimal)
    for c in candidates:P.verify(d,c)
    primal=min((c for c in candidates if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL'),key=lambda c:F(c['maximum']))
    dual=max((c for c in candidates if c['kind']=='EXACT_PARTIAL_MODEL_DUAL'),key=lambda c:F(c['bound']))
    lo=F(dual['bound']);hi=F(primal['maximum']);assert lo<=hi
    out.update(status='PASS_EXACT_ALL_MARKER_FAMILIES_CONTINUATION_INTERVAL',solver=solver,
        solver_history=[first['solver'],solver],certificates=[primal,dual],numerical_snapshot=path,numerical_snapshot_sha256=h,
        candidate_certificate_summary=[dict(kind=c['kind'],value=c.get('maximum',c.get('bound')))for c in candidates],
        exact_interval=dict(lower=str(lo),upper=str(hi),gap=str(hi-lo),lower_exceeds12=lo>12,upper_at_most12=hi<=12),
        same_exact_matrix=True,numerical_method_changed=True,
        scope='Two chronological numerical methods on one identical finite exported matrix. Exact feasibility is decided by full original-matrix checks. Numerical optimality and exact gap are separate; no complete geometric coverage or H13 exclusion is inferred.')
    OUT.write_text(P.BASE.dump(out)+'\n');print(out['status'],out['exact_interval'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();run(args.verify)
