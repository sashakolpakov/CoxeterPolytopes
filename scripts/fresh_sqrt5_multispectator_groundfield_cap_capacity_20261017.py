#!/usr/bin/env python3
"""Fixed-Q(sqrt5) strengthening of the frozen two-A1 cap family."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_sqrt5_multispectator_cap_capacity_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_multispectator_groundfield_cap_capacity_20261017'
MODEL='fresh_sqrt5_multispectator_cap_joint_ports_20261017'


def run():
    old=M.read(M.TAG);model=M.read(MODEL)
    assert old['status']=='PASS_EXACT_TWO_A1_MULTISPECTATOR_CAP_CAPACITIES'
    alpha={i:F(a)for i,a in model['face_multipliers']}
    tests=[];patterns={}
    for row in old['tau_pair_tests']:
        tau=M.K(*map(F,row['tau']));h=(3-tau)/2
        parameters=[(h,0),(h+F(1,2),0),(h,F(1,2)),(h,1)]
        pairs={}
        for name,(a,b) in zip(M.PATTERN_NAMES,parameters):
            rec=M.pair_gate(a,b,tau)
            assert rec==row['pairs'][name]
            rec['passes_acuteness_and_integrality']=rec['allowed_by_necessary_tests']
            rec['two_cycle_in_Qsqrt5']=rec['branch']!='relative_quadratic_norm'
            rec['allowed_by_necessary_tests']=rec['passes_acuteness_and_integrality'] and rec['two_cycle_in_Qsqrt5']
            pairs[name]=rec
        pattern=tuple(pairs[name]['allowed_by_necessary_tests']for name in M.PATTERN_NAMES)
        assert all(not b or a for a,b in zip(row['pattern'],pattern))
        patterns[tuple(row['tau'])]=pattern
        tests.append(dict(tau=row['tau'],h=h.out(),pairs=pairs,old_pattern=row['pattern'],pattern=list(pattern)))
    columns=[];required=defaultdict(F)
    for oldc in old['columns']:
        pattern=patterns[tuple(oldc['tau'])]
        one,kind=M.cap(pattern,oldc['A1_count'])
        capacity=one*oldc['orientation_count']
        assert 0<capacity<=oldc['capacity']
        c=dict(oldc,previous_capacity=oldc['capacity'],capacity=capacity,
               single_orientation_capacity=one,proof_kind=kind,pair_pattern=list(pattern))
        load=sum((m*alpha.get(i,F())for i,m in c['positive_face_coefficients']),F())
        mass=alpha.get(c['simplex_face'],F())
        c.update(load=str(load),initial_simplex_mass=str(mass),slack=str(capacity*mass-load))
        assert c['positive_face_coefficients']==oldc['positive_face_coefficients']
        assert c['provenance_coefficients']==oldc['provenance_coefficients']
        required[c['simplex_face']]=max(required[c['simplex_face']],load/capacity)
        columns.append(c)
    assert len(columns)==555 and len(tests)==11
    violations=[c['index']for c in columns if F(c['slack'])<0]
    paths=[Path(__file__),Path(M.__file__),ROOT/f'results/{M.TAG}.json',
           ROOT/f'results/{M.TAG}.md',ROOT/f'results/{MODEL}.json']
    return dict(status='PASS_EXACT_FIXED_GROUNDFIELD_TWO_A1_CAP_TIGHTENING',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        arithmetic_ground_field='Q(sqrt5)',inherited_raw_cap_package=M.TAG,
        pair_pattern_order=M.PATTERN_NAMES,tau_pair_tests=tests,columns=columns,
        counts=dict(columns=555,tau_values=11,pair_tests=44,
            tightened_columns=sum(c['capacity']<c['previous_capacity']for c in columns),
            positive_columns=sum(F(c['load'])>0 for c in columns),violations=len(violations)),
        violating_columns=violations,
        capacity_histogram=sorted(Counter(c['capacity']for c in columns).items()),
        required_simplex_masses=[[i,str(alpha.get(i,F())),str(a)]for i,a in sorted(required.items())],
        input_bound=model['bound'],new_repaired_bound=None,optimizer_invocations=0,
        scope='Fresh strengthening under fixed arithmetic ground field Q(sqrt5): every normalized norm2 two-cycle squared product must lie in K. The exact frozen555 raw maps/projection norms/colored orientation multiplicities are inherited unchanged. All44 pair tests are recomputed; outside-K squared products remove compatible-position edges. Existing proved matching capacities apply to the resulting patterns. No optimizer, repaired-witness claim, all-arithmetic-field extension or dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],flush=True)
