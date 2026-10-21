#!/usr/bin/env python3
"""Preserve an explicitly unresolved first-stop optimization, not a certificate."""
import argparse,ast,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_global_b1_first_repair_stopped_20261017'
STATE=f'data/enumerations/{TAG}.json';LOG=f'data/enumerations/{TAG}.log'
INITIAL='data/enumerations/fresh_sqrt5_global_b1_first_repair_20261016_initial.json'

def run():
    state=json.loads((ROOT/STATE).read_text());initial=json.loads((ROOT/INITIAL).read_text())
    for doc in(state,initial):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert state['vertex_cuts']==initial['vertex_cuts']and len(state['vertex_cuts'])==4762
    assert not state['steps']and state['latest_numeric_weights']is None
    log=(ROOT/LOG).read_text();lines=[r for r in log.splitlines()if r.startswith('GLOBAL_FIRST_REPAIR_MODEL ')]
    assert len(lines)==1 and '32554 variables 70156 profiles 'in lines[0]
    counts=ast.literal_eval(lines[0][lines[0].index('{'):]);assert counts['finite_first_singular_rows']==196 and counts['dotted_first_singular_rows']==14
    assert counts['finite_first_common_contexts']==98877 and counts['dotted_first_common_contexts']==10304 and counts['old_rule_unpaid_stops']==52
    assert log.rstrip().endswith('RESUMED 4762 literal cuts')and not any(s in log for s in('SECONDARY ','EXACT_DUAL ','EXACT_BASIS '))
    paths=[f'scripts/{TAG}.py',STATE,LOG,INITIAL,'scripts/fresh_sqrt5_global_b1_first_repair_20261016.py',
           'scripts/fresh_sqrt5_global_discount_20261016.py','results/fresh_sqrt5_global_discount_20261016.json',
           'results/fresh_sqrt5_global_b1_finite_first_20261016.json',
           'results/fresh_sqrt5_global_b1_dotted_first_completions_20261016.json']
    return dict(status='PASS_PRESERVED_UNRESOLVED_GLOBAL_FIRST_STOP_REPAIR',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=21,target='40',
        variables=32554,unique_face_profiles=70156,face_counts=counts,preserved_literal_vertex_cuts=4762,
        completed_optimization_iterations=0,mathematical_solver_result=None,
        source_cut_validation='The original wrapper literally recounted and independently reassembled every inherited global forest cut before the preserved initial checkpoint; this checker binds that completed handoff without rerunning its arithmetic.',
        stop='Only the identified task-owned PID89604 was terminated after the prolonged initial optimization. Checkpoint and log were copied before termination; authoritative session94072 returned143. No infeasibility, feasible vector, optimal value or completed oracle iteration is inferred.',
        continuation='A separate bounded IPM wrapper uses the identical variables, constraints and oracle; its outcome is outside this preservation record.',
        scope='Reproducible unresolved optimization state only. The exact196 finite and14 dotted first-stop domains and all common contexts were rebuilt, revealing52 original-rule unpaid stops, but this run did not produce a repaired rule or an exact obstruction. Complete deeper quadrilateral closure is independently pending. No dimension bound or new polytope claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
