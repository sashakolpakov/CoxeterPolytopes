#!/usr/bin/env python3
"""Same H3-center family; fixed35 and a bounded-below secondary objective.

All bonuses are nonnegative and have objective coefficient1. Discounts are
bounded0..1 and T is fixed35, so the objective cannot be unbounded below.
Every iteration writes literal vertex cuts to an atomic temporary checkpoint.
"""
import argparse,json,os
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import numpy as np
import highspy
import flint
import fresh_sqrt5_marked_h3_center_cross_20261016 as M
BASE_TAG=M.TAG;TAG='fresh_sqrt5_marked_h3_center_sparse_recovery_20261016';ROOT=M.ROOT
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
RESUME=False

def matrix(rows,n):
    rr=[];cc=[];vv=[]
    for i,row in enumerate(rows):
        for j,z in row.items():rr.append(i);cc.append(j);vv.append(float(z))
    return M.coo_matrix((vv,(rr,cc)),shape=(len(rows),n)).tocsr()

def highs_model(rows,rhs,n):
    A=matrix(rows,n);lp=highspy.HighsLp();lp.num_col_=n;lp.num_row_=len(rows)
    lp.col_cost_=np.array([0,0]+[1]*(n-3)+[0],dtype=float)
    lp.col_lower_=np.array([0]*(n-1)+[35],dtype=float)
    lp.col_upper_=np.array([1,1]+[highspy.kHighsInf]*(n-3)+[35],dtype=float)
    lp.row_lower_=np.full(len(rows),-highspy.kHighsInf);lp.row_upper_=np.array(list(map(float,rhs)))
    lp.a_matrix_.format_=highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_=A.indptr;lp.a_matrix_.index_=A.indices;lp.a_matrix_.value_=A.data
    h=highspy.Highs();h.setOptionValue('output_flag',False);h.setOptionValue('solver','simplex')
    h.setOptionValue('threads',1);h.setOptionValue('random_seed',0)
    h.setOptionValue('primal_feasibility_tolerance',1e-9);h.setOptionValue('dual_feasibility_tolerance',1e-9)
    assert h.passModel(lp)==highspy.HighsStatus.kOk
    return h

def append_rows(h,rows,rhs,n):
    if not rows:return
    A=matrix(rows,n)
    assert h.addRows(len(rows),np.full(len(rows),-highspy.kHighsInf),np.array(list(map(float,rhs))),
        len(A.data),A.indptr,A.indices,A.data)==highspy.HighsStatus.kOk

def exact_basis(h,rows,rhs,n):
    import heapq
    basis=h.getBasis();B=highspy.HighsBasisStatus
    basic=[j for j,s in enumerate(basis.col_status)if s==B.kBasic]
    tight=[i for i,s in enumerate(basis.row_status)if s!=B.kBasic]
    assert len(basic)==len(tight),(len(basic),len(tight))
    x=[F(0)]*n
    for j,s in enumerate(basis.col_status):
        if s==B.kBasic:continue
        if j==n-1:x[j]=F(35)
        elif s==B.kUpper:assert j<2;x[j]=F(1)
        else:assert s in(B.kLower,B.kZero,B.kNonbasic)
    active=set(basic);eqs={};right={};colrows={j:set()for j in basic};heap=[];version={}
    q=lambda z:flint.fmpq(z.numerator,z.denominator)
    for k,i in enumerate(tight):
        assert basis.row_status[i]in(B.kUpper,B.kNonbasic)
        eqs[k]={j:q(z)for j,z in rows[i].items()if j in active}
        right[k]=q(rhs[i]-sum(z*x[j]for j,z in rows[i].items()if j not in active))
        assert eqs[k]
        for j in eqs[k]:colrows[j].add(k)
        version[k]=0;heapq.heappush(heap,(len(eqs[k]),0,k))
    pivots=[];initial=len(active)
    print('EXACT_BASIS_SPARSE_START',initial,'nonzeros',sum(map(len,eqs.values())),flush=True)
    while eqs:
        while heap:
            size,v,k=heapq.heappop(heap)
            if k in eqs and version[k]==v:break
        else:raise AssertionError('Lost sparse pivot rows')
        assert size==len(eqs[k])and size>0
        # Sparse rows make exact substitution cheap. Leave a small dense
        # remainder to FLINT rather than create large fill in Python.
        if len(active)<=400 or size>16:break
        j=min(eqs[k],key=lambda j:(len(colrows[j]),j));a=eqs[k][j]
        normalized={v:z/a for v,z in eqs[k].items()if v!=j};b=right[k]/a
        neighbors=sorted(colrows[j]-{k})
        for v in eqs[k]:colrows[v].remove(k)
        del eqs[k];del right[k];active.remove(j);pivots.append((j,b,normalized))
        for i in neighbors:
            coefficient=eqs[i].pop(j);colrows[j].remove(i);right[i]-=coefficient*b
            for v,z in normalized.items():
                old=eqs[i].get(v,flint.fmpq(0));new=old-coefficient*z
                if new:
                    eqs[i][v]=new
                    if not old:colrows[v].add(i)
                elif old:del eqs[i][v];colrows[v].remove(i)
            assert eqs[i],('Dependent basis row',i,right[i])
            version[i]+=1;heapq.heappush(heap,(len(eqs[i]),version[i],i))
        assert not colrows[j]
        if len(pivots)%1000==0:print('EXACT_SPARSE_ELIMINATED',len(pivots),'remaining',len(active),flush=True)
    active=sorted(active);position={j:k for k,j in enumerate(active)};ids=sorted(eqs)
    assert len(active)==len(ids)
    print('EXACT_DENSE_REMAINDER',len(active),'nonzeros',sum(map(len,eqs.values())),flush=True)
    rational={}
    if active:
        A=flint.fmpq_mat(len(active),len(active));b=flint.fmpq_mat(len(active),1)
        for k,i in enumerate(ids):
            b[k,0]=right[i]
            for j,z in eqs[i].items():A[k,position[j]]=z
        solution=A.solve(b,algorithm='dixon');assert A*solution==b
        rational.update({j:solution[k,0]for j,k in position.items()})
    for j,b,row in reversed(pivots):rational[j]=b-sum((z*rational[v]for v,z in row.items()),flint.fmpq(0))
    for j,v in rational.items():x[j]=F(str(v))
    assert x[-1]==35 and min(x)>=0 and max(x[:2])<=1
    assert all(sum(z*x[j]for j,z in row.items())<=bound for row,bound in zip(rows,rhs))
    return x,dict(basic_variables=initial,exact_fraction_equations=len(tight),sparse_pivots=len(pivots),dense_remainder=len(active))


def checkpoint(vertices,steps,weights):
    doc=dict(status='PARTIAL_CHECKPOINT_NOT_A_CERTIFICATE',sources={
        f'scripts/{BASE_TAG}.py':sha256((ROOT/f'scripts/{BASE_TAG}.py').read_bytes()).hexdigest(),
        f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()},
        vertex_cuts=vertices,steps=steps,latest_numeric_weights=list(map(float,weights))if weights is not None else None)
    temporary=CHECKPOINT.with_suffix('.tmp');temporary.write_text(json.dumps(doc,separators=(',',':'))+'\n');os.replace(temporary,CHECKPOINT)

def solve(obligations):
    n=M.N+1;rows=[{j:-z for j,z in p.items()if j>=0}for p,t,w in obligations]
    rhs=[p.get(-1,F(0))-t for p,t,w in obligations];meta=[w for p,t,w in obligations]
    cuts=set();vertices=[];steps=[]
    def cut(front):
        profile=front['profile'];k=M.key(profile)
        if k in cuts:return False
        cuts.add(k);row={j:z for j,z in profile.items()if j>=0};row[M.N]=F(-1)
        bound=-profile.get(-1,F(0));m=dict(kind='vertex',**{k:v for k,v in front.items()if k not in('profile','maximum','kind')})
        rows.append(row);rhs.append(bound);meta.append(m)
        vertices.append(dict(row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(bound),metadata=m))
        return True
    if RESUME and CHECKPOINT.is_file():
        prior=json.loads(CHECKPOINT.read_text())
        for p,wanted in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==wanted
        for rec in prior['vertex_cuts']:
            profile={int(j):F(z)for j,z in rec['row']if int(j)!=M.N};
            if F(rec['rhs']):profile[-1]=-F(rec['rhs'])
            cut(dict(**rec['metadata'],profile=profile))
        steps=prior['steps'];print('RESUMED',len(cuts),'literal cuts',flush=True)
    else:
        for front in M.oracle([0.0]*M.N):cut(front)
    checkpoint(vertices,steps,None)
    h=highs_model(rows,rhs,n)
    for iteration in range(len(steps),200):
        assert h.run()==highspy.HighsStatus.kOk;status=h.getModelStatus()
        if status==highspy.HighsModelStatus.kInfeasible:
            A=matrix(rows,n);objective=np.zeros(n);objective[M.N]=1
            lp=M.linprog(objective,A_ub=A,b_ub=np.array(list(map(float,rhs))),
                bounds=[(0,1),(0,1)]+[(0,None)]*(M.N-1),method='highs')
            assert lp.success and lp.fun>35+1e-7,(lp.message,lp.fun)
            cert=M.dual_certificate(rows,rhs,meta,lp);checkpoint(vertices,steps,lp.x)
            return cert,steps
        assert status==highspy.HighsModelStatus.kOptimal,h.modelStatusToString(status)
        weights=np.array(h.getSolution().col_value);fronts=M.oracle(weights[:M.N]);maximum=max(f['maximum']for f in fronts)
        steps.append(dict(iteration=iteration,fixed_target=35,secondary_objective=h.getObjectiveValue(),
            cuts=len(cuts),maximum=float(maximum),discounts=list(map(float,weights[:2]))))
        print('SECONDARY',iteration,'sum',h.getObjectiveValue(),'cuts',len(cuts),'maximum',maximum,'discounts',weights[:2],flush=True)
        checkpoint(vertices,steps,weights)
        if maximum<=35+1e-7:
            x,basis_record=exact_basis(h,rows,rhs,n);fronts=M.oracle(x[:M.N],True)
            upper=max(f['maximum']for f in fronts);print('EXACT_BASIS',basis_record,'maximum',float(upper),flush=True)
            if upper<=35:
                assert all(M.value(p,x)>=t for p,t,w in obligations)
                return dict(kind='exact_feasible_restricted_rule',primal=list(map(str,x[:M.N])),maximum=str(upper),
                    exact_basis_recovery=basis_record),steps
        before=len(rows)
        for front in fronts:
            if front['maximum']>35+(0 if isinstance(front['maximum'],F)else 1e-7):cut(front)
        assert len(rows)>before,'No new cut despite an unresolved exact feasibility check'
        append_rows(h,rows[before:],rhs[before:],n);checkpoint(vertices,steps,weights)
    raise AssertionError('Secondary iteration cap; literal cuts saved')

def run(saved=None):
    M.solve=solve;M.TAG=TAG
    result=M.run(saved)
    result['status']='PASS_EXACT_H3_CENTER_SECONDARY_MARKED_PROPOSAL'
    result['sources'][f'scripts/{BASE_TAG}.py']=sha256((ROOT/f'scripts/{BASE_TAG}.py').read_bytes()).hexdigest()
    for path in ['scripts/fresh_sqrt5_marked_h3_center_secondary_20261016.py','data/enumerations/fresh_sqrt5_marked_h3_center_secondary_20261016_checkpoint.json','results/fresh_sqrt5_marked_h3_center_secondary_numerical_20261016.json']:
        result['sources'][path]=sha256((ROOT/path).read_bytes()).hexdigest()
    result['optimizer_scope']='Same constraint family as the H3-center source; T fixed35, minimize sum of all nonnegative bonus coefficients. Discounts are bounded0..1 and have zero objective coefficient. Objective is bounded below by0, and any feasible objective level bounds every bonus variable. Atomic literal cut checkpoints are temporary research state, not certificates.'
    result['predecessor_handoff']='The pure min-T process stopped without a cut checkpoint. Its secondary replacement reached numerical full-oracle35 and checkpointed17972 literal cuts. This run resumes that literal checkpoint and recovers its rational basis by exact sparse elimination followed by a smaller FLINT solve; it changes no model constraints.'
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--resume',action='store_true');args=ap.parse_args();RESUME=args.resume
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
