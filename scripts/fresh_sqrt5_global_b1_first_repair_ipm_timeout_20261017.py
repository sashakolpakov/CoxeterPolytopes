#!/usr/bin/env python3
"""Replay the evidence of a bounded unresolved optimization; no LP claim."""
import argparse,json,re
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017'
PARENT='fresh_sqrt5_global_b1_first_repair_ipm_20261017'

def run():
    state_path=f'data/enumerations/{TAG}.json';log_path=f'data/enumerations/{TAG}.log'
    initial_path=f'data/enumerations/{PARENT}_initial.json'
    state=json.loads((ROOT/state_path).read_text());initial=json.loads((ROOT/initial_path).read_text())
    for doc in(state,initial):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert state['vertex_cuts']==initial['vertex_cuts'] and len(state['vertex_cuts'])==4762
    assert not state['steps'] and state['latest_numeric_weights'] is None
    log=(ROOT/log_path).read_text()
    assert log.count('IPM_FIXED40_START ')==1 and log.count('IPM_FIXED40_END ')==1
    assert 'UNCHANGED_FIRST_REPAIR_MODEL 32554 variables 70156 profiles 'in log
    end=re.findall(r'IPM_FIXED40_END 1 HighsStatus.kWarning Time limit reached seconds ([0-9.]+)',log)
    assert len(end)==1 and 300<=float(end[0])<310
    assert log.rstrip().endswith('AssertionError')
    assert 'IPM_MINIMUM_T_START'not in log and 'EXACT_BASIS 'not in log and 'EXACT_DUAL 'not in log
    paths=[f'scripts/{TAG}.py',f'scripts/{PARENT}.py',state_path,log_path,initial_path,
           'scripts/fresh_sqrt5_global_b1_first_repair_20261016.py',
           'results/fresh_sqrt5_global_b1_first_repair_stopped_20261017.json']
    return dict(status='PASS_PRESERVED_BOUNDED_UNRESOLVED_GLOBAL_REPAIR',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        dimension=21,target='40',variables=32554,unique_face_profiles=70156,
        preserved_literal_vertex_cuts=4762,completed_optimization_iterations=0,
        configured_solver_seconds=300,reported_solver_seconds=end[0],
        mathematical_solver_result=None,authoritative_exit_status=1,
        solver_outcome='Interior-point stopped making progress; automatic simplex cleanup reached the configured 300-second limit. The inherited exact-recovery driver rejected the non-success status, so no basis recovery, oracle iteration or fallback minimum-T solve occurred.',
        scope='Preserved unresolved numerical experiment only. No feasible rule, infeasibility, objective bound or geometric dimension conclusion follows. The model retains every initial cut and every selected face profile; complete deeper B1 face coverage remains a separate task.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
