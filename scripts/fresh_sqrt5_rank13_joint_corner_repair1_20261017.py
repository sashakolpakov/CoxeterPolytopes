#!/usr/bin/env python3
"""First counterexample-guided successor of the rank13 corner-orbit LP."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_joint_corner_repair1_20261017'
SEP='fresh_sqrt5_rank13_joint_corner_separation_20261017'

def domain():
    d=B.domain();sep=B.read(SEP)
    known={(f['required'],tuple(map(tuple,f['coefficients'])))for f in d['faces']}
    for w in sep['witnesses']:
        coef=sorted(Counter(w['corner_variable_indices']).items());assert(1,tuple(coef))not in known
        for c in w['corners']:
            v=d['variables'][c['variable']]
            assert v['forest']==c['forest']and v['marked_forest_code']==c['marked_forest_code']
        d['faces'].append(dict(index=len(d['faces']),required=1,coefficients=[list(x)for x in coef],
            provenance=[dict(kind='exact_omitted_compact_local_quad',source=SEP,class_index=w['class_index'],
                spectator_forest=w['spectator_forest'])]))
    d['counts']['added_exact_compact_local_quads']=len(sep['witnesses'])
    d['counts']['literal_face_rows']+=len(sep['witnesses'])
    d['counts']['distinct_face_constraints']=len(d['faces'])
    for p in[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'results/{B.TAG}.json',f'results/{SEP}.json']:
        d['sources'][p]=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    return d

def run(verify):
    d=domain();path=ROOT/f'results/{TAG}.json'
    if verify:
        saved=json.loads(path.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:B.verify_certificate(d,c)
        return saved
    solver,certificates=B.optimize(d)
    feasible=any(c.get('reaches_target12',False)for c in certificates)
    excluded=any(c.get('excludes_target12',False)for c in certificates);assert not(feasible and excluded)
    status='PASS_EXACT_FIRST_REPAIRED_PARTIAL_RULE_AT12'if feasible else'PASS_EXACT_FIRST_REPAIRED_PARTIAL_MODEL_OBSTRUCTION'if excluded else'UNRESOLVED_FIRST_REPAIRED_PARTIAL_MODEL'
    d.update(status=status,dimension=13,target='12',solver=solver,certificates=certificates,
        change='Two independently checked actual compact local quadrilaterals, classes4266/7746 with A2+D5/D5 spectators, are added to the unchanged full original partial-domain constraints. Both had zero payment under the frozen predecessor. No old coefficient is fixed and no prior constraint is removed.',
        scope='Counterexample-guided partial-model successor only. It still lacks complete rank13 triangle and quadrilateral coverage and all absent finite vertex types. A feasible budget12 vector is not a dimension theorem; a restricted-model obstruction would not exclude more global geometric restrictions.')
    path.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound')))for c in d['certificates']],flush=True)
