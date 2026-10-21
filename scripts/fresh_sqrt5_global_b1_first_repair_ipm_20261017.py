#!/usr/bin/env python3
"""Bounded IPM replacement for the unchanged first-stop global40 model.
No constraint, normalization, variable or vertex-oracle change is made.
"""
import argparse,json,os,time
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_global_b1_first_repair_20261016 as R
G,M,ROOT=R.G,R.M,R.ROOT
TAG='fresh_sqrt5_global_b1_first_repair_ipm_20261017';N=R.N
INITIAL=f'data/enumerations/{TAG}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
OLD_CHECKPOINT=R.CHECKPOINT
OLD_HIGHS=G.highs_model;OLD_LINPROG=M.linprog


class LoggedModel:
    def __init__(self,h):self.h=h;self.calls=0
    def __getattr__(self,name):return getattr(self.h,name)
    def run(self):
        self.calls+=1;start=time.monotonic();print('IPM_FIXED40_START',self.calls,flush=True)
        status=self.h.run()
        print('IPM_FIXED40_END',self.calls,str(status),self.h.modelStatusToString(self.h.getModelStatus()),
              'seconds',round(time.monotonic()-start,3),flush=True)
        return status


def highs(rows,rhs,n):
    h=OLD_HIGHS(rows,rhs,n)
    for key,value in [('solver','ipm'),('time_limit',300.0),('output_flag',True),('run_crossover','on')]:
        assert h.setOptionValue(key,value)==G.highspy.HighsStatus.kOk
    return LoggedModel(h)


def minimum_t(*args,**kw):
    kw['method']='highs-ipm';options=dict(kw.get('options')or{});options['time_limit']=300.0;options['disp']=True;kw['options']=options
    start=time.monotonic();print('IPM_MINIMUM_T_START',flush=True);answer=OLD_LINPROG(*args,**kw)
    print('IPM_MINIMUM_T_END',answer.message,answer.fun,'seconds',round(time.monotonic()-start,3),flush=True)
    return answer


def install():
    R.install();G.TAG=TAG;G.INITIAL=INITIAL;G.CHECKPOINT=CHECKPOINT
    G.highs_model=highs;M.linprog=minimum_t


def prepare_initial():
    p=ROOT/INITIAL
    if p.exists():return
    old=json.loads(OLD_CHECKPOINT.read_text())
    for source,h in old['sources'].items():assert sha256((ROOT/source).read_bytes()).hexdigest()==h
    # The previous wrapper already literally reconstructed all these cuts.
    # Preserve every row and its complete forest metadata without dropping cuts.
    assert len(old['vertex_cuts'])>=4762 and not old['steps']
    original=json.loads((ROOT/R.INITIAL).read_text())
    assert old['vertex_cuts']==original['vertex_cuts']
    paths=[f'scripts/{TAG}.py',f'scripts/{R.TAG}.py',R.INITIAL,R.FINITE,R.DOTTED]
    d=dict(status='PRESERVED_FIRST_REPAIR_LITERAL_CUTS_NO_SOLVER_RESULT',
        sources={q:sha256((ROOT/q).read_bytes()).hexdigest()for q in paths},
        vertex_cuts=old['vertex_cuts'],steps=[],latest_numeric_weights=None,
        prior_solver='Task-owned PID89604 stopped after an opaque initial solve; no completed optimization iterate or exact result was claimed.')
    p.write_text(json.dumps(d,separators=(',',':'))+'\n')
    print('IPM_PRESERVED_CUTS',len(d['vertex_cuts']),flush=True)


def run(saved=None):
    install();obligations,counts,failures=R.faces()
    print('UNCHANGED_FIRST_REPAIR_MODEL',N,'variables',len(obligations),'profiles',counts,flush=True)
    cert,steps=G.solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=G.verify(obligations,cert)
    paths=[f'scripts/{TAG}.py',INITIAL,f'scripts/{R.TAG}.py',R.FINITE,R.DOTTED,R.OLD_RESULT,
           'scripts/fresh_sqrt5_global_discount_20261016.py']
    old=json.loads((ROOT/R.OLD_RESULT).read_text())
    return dict(status='PASS_EXACT_GLOBAL_FIRST_STOP_REPAIR_IPM_MODEL',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=21,target='40',variables=N,
        variable_order=old['variable_order'],face_counts=counts,unique_face_profiles=len(obligations),
        original_rule_first_stop_failures=failures,certificate=cert,verification=proof,steps=steps,
        vertex_context_order=old['vertex_context_order'],triangle_external_context_order=old['triangle_external_context_order'],
        quadrilateral_external_context_order=old['quadrilateral_external_context_order'],
        scope='Unchanged global40 first-stop repair: all196 finite and14 unpaid-dotted admissible singular first stops, every common spectator context, all996 triangles,37 earlier paid stops and1979 zero-payment quad rows. Complete deeper B1 closure remains separate and open. Restricted model only; no dimension bound.',
        solver_change='Interior-point with crossover,300-second per-solve limit and explicit stages; fallback minimum-T also uses IPM with300-second limit. No overlapping optimization loop. Every inherited cut is preserved exactly; rational recovery/full oracle requirements are unchanged.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('TASK_PID',os.getpid(),flush=True);install()
    if not a.verify:
        prepare_initial()
        if not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text())if a.verify else None
    d=json.loads(json.dumps(run(saved),default=str))
    if a.verify:assert d==saved
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    value=d['verification'].get('maximum',d['verification'].get('lower_bound'))
    print(d['status'],d['certificate']['kind'],value,flush=True)
