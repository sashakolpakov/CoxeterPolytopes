#!/usr/bin/env python3
"""Five-minute zero-discount restriction of the existing global payment family."""
import argparse,json,os,time
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_global_discount_20261016 as G
M,X,ROOT=G.M,G.X,G.ROOT;N=G.N
TAG='fresh_sqrt5_global_zero_discount_20261017'
INITIAL=f'data/enumerations/{TAG}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
PRIOR_INITIAL='data/enumerations/fresh_sqrt5_global_b1_first_repair_20261016_initial.json'
OLD_HIGHS=G.highs_model;OLD_LINPROG=M.linprog
DEADLINE=None

class BudgetEnded(Exception):pass

def remaining():
    answer=DEADLINE-time.monotonic()
    if answer<=0:raise BudgetEnded('Five-minute optimization budget exhausted')
    return answer

class Bounded:
    def __init__(self,h):self.h=h
    def __getattr__(self,key):return getattr(self.h,key)
    def run(self):
        left=remaining();assert self.h.setOptionValue('time_limit',left)==G.highspy.HighsStatus.kOk
        print('ZERO_FIXED40_START',left,flush=True)
        status=self.h.run();model=self.h.getModelStatus()
        print('ZERO_FIXED40_END',str(status),self.h.modelStatusToString(model),flush=True)
        if model==G.highspy.HighsModelStatus.kTimeLimit:raise BudgetEnded('Solver time limit reached')
        return status

def highs(rows,rhs,n):
    h=OLD_HIGHS(rows,rhs,n);h.setOptionValue('output_flag',True)
    return Bounded(h)

def minimum_t(*args,**kw):
    left=remaining();kw['options']=dict(kw.get('options')or{},time_limit=left,disp=True)
    print('ZERO_MINIMUM_T_START',left,flush=True);answer=OLD_LINPROG(*args,**kw)
    print('ZERO_MINIMUM_T_END',answer.message,answer.fun,flush=True)
    if answer.status==1:raise BudgetEnded('Minimum-T solver time limit reached')
    return answer

def install():
    G.install();G.TAG=TAG;G.INITIAL=INITIAL;G.CHECKPOINT=CHECKPOINT
    G.highs_model=highs;M.linprog=minimum_t

def faces():
    obligations={};old=M.H.read('fresh_sqrt5_parameter_weight_stress_20261015')
    queue=M.H.read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    rows=[('prior',i,r['row'])for i,r in enumerate(old['faces'])if r['row']['nodes']<=23]
    rows += [('finite_queue',r['index'],r['row'])for r in queue['records']if r['result']['singular']and r['row']['nodes']<=23]
    cases=0
    for family,index,row in rows:
        es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        compiled=[X.OLD_COMPILE(row['white']+[a,b],es,[a,b],[None])for a,b in product((0,2),(1,3))]
        assert all(not any(v[3][None])for v in compiled)
        for ext in M.contexts(23-row['nodes']):
            p,tables=X.OLD_ADD(compiled,ext)
            assert not any(j in X.SPLIT for j in p)and all(not t[None]for t in tables)
            assert not any(j in p for j in(0,1))and all(z>=0 for j,z in p.items()if j>=0)
            obligations.setdefault((M.key(p),F(1)),dict(kind='quad',family=family,index=index,external_counts=list(ext)))
            cases+=1
    # Both variables stay in their original order but are explicitly zero.
    for j in(0,1):obligations[M.key({j:F(-1)}),F(0)]=dict(kind='model_bound',name='discount_zero',index=j)
    tri=json.loads((ROOT/f'results/{X.A.TRI}_rank22.json').read_text())
    assert len(tri['records'])==996 and not tri['pending']and tri['complete_through_bound']
    assert all(F(r['baseline'])>=2 for r in tri['records'])
    return [(dict(p),t,meta)for(p,t),meta in obligations.items()],dict(
        baseline_zero_rows=len(rows),common_contexts=cases,complete_triangle_baseline_minimum=str(min(F(r['baseline'])for r in tri['records'])),
        automatic_paid_quad_argument='Both discounts are zero and every added payment is nonnegative; every quadrilateral with baseline at least1 retains payment at least1, irrespective of white extensions.')

def prepare():
    p=ROOT/INITIAL
    if not p.exists():
        d=json.loads((ROOT/PRIOR_INITIAL).read_text())
        for q,h in d['sources'].items():assert sha256((ROOT/q).read_bytes()).hexdigest()==h
        assert len(d['vertex_cuts'])==4762
        paths=[f'scripts/{TAG}.py',PRIOR_INITIAL]
        out=dict(status='UNCHANGED_GLOBAL_LITERAL_CUTS_ZERO_DISCOUNT_RESTRICTION',
            sources={q:sha256((ROOT/q).read_bytes()).hexdigest()for q in paths},
            vertex_cuts=d['vertex_cuts'],steps=[],latest_numeric_weights=None)
        p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    if not CHECKPOINT.exists():CHECKPOINT.write_bytes(p.read_bytes())

def run(saved=None):
    global DEADLINE
    install();obligations,counts=faces();print('ZERO_MODEL',N,len(obligations),counts,flush=True)
    certificate=None;proof=None;outcome=None;steps=[]
    if saved is None:
        DEADLINE=time.monotonic()+300
        try:
            certificate,steps=G.solve(obligations)
            proof=G.verify(obligations,certificate)
            outcome='EXACT_RESTRICTED_FAMILY_RESULT'
        except BudgetEnded as error:
            outcome='UNRESOLVED_TIME_LIMIT';print('ZERO_UNRESOLVED',str(error),flush=True)
            state=json.loads(CHECKPOINT.read_text());steps=state['steps']
            archived=ROOT/f'data/enumerations/{TAG}_terminal_checkpoint.json';archived.write_bytes(CHECKPOINT.read_bytes())
    else:
        certificate=saved['certificate'];steps=saved['steps'];outcome=saved['outcome']
        if certificate is not None:proof=G.verify(obligations,certificate)
        else:
            state=json.loads((ROOT/f'data/enumerations/{TAG}_terminal_checkpoint.json').read_text())
            assert outcome=='UNRESOLVED_TIME_LIMIT'and steps==state['steps']
            for q,h in state['sources'].items():assert sha256((ROOT/q).read_bytes()).hexdigest()==h
    paths=[f'scripts/{TAG}.py',INITIAL,'scripts/fresh_sqrt5_global_discount_20261016.py',
           f'results/{X.A.TRI}_rank22.json','results/fresh_sqrt5_parameter_weight_stress_20261015.json',
           'results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json']
    if certificate is None:paths.append(f'data/enumerations/{TAG}_terminal_checkpoint.json')
    return dict(status='PASS_REPLAYABLE_ZERO_DISCOUNT_EXPERIMENT',outcome=outcome,
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=21,target='40',variables=N,
        optimization_wall_budget_seconds=300,face_counts=counts,unique_face_profiles=len(obligations),
        certificate=certificate,verification=proof,steps=steps,
        scope='Existing global payment family with both baseline discounts fixed to zero. All baseline-positive quadrilaterals and all996 triangles are automatic under nonnegative bonuses. The inherited complete baseline-zero obligations are retained. A feasible result still requires exact rational recovery and the full global finite-forest oracle; an unresolved solve implies neither feasibility nor infeasibility. No automatic dimension-bound promotion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    print('TASK_PID',os.getpid(),flush=True);install()
    if not a.verify:prepare()
    p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text())if a.verify else None
    d=json.loads(json.dumps(run(saved),default=str))
    if a.verify:assert d==saved
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['outcome'],flush=True)
