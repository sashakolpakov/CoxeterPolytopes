#!/usr/bin/env python3
"""One-coordinate monotone repair of an exact face-incidence witness.

All face masses and every inherited geometric load remain unchanged.  The
only edit is a nonnegative increment to one vertex mass. No LP or face scan.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_core23_25_incidence_repair_20261017'
BASE='fresh_sqrt5_a3_d5_cap_joint_ports_20261017'
GATE='fresh_sqrt5_core23_25_forced_quad_20261017'
LOCAL='fresh_sqrt5_a3_d5_cap_local_cover_update_20261017'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'


def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for obj in (d,d.get('prepared',{})):
        for p,h in obj.get('sources',{}).items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d


def run():
    old,gate,local,forests=map(read,[BASE,GATE,LOCAL,FORESTS])
    assert old['status'].startswith('PASS') and gate['status']=='PASS_EXACT_CORE23_25_FORCED_QUADRILATERAL_INCIDENCE'
    assert old['counts']['total_inequalities']==2711432 and old['counts']['violated_columns']==0
    assert old['face_multipliers']==local['face_multipliers']
    assert old['vertex_multipliers']==local['vertex_multipliers']
    alpha={f:F(a)for f,a in old['face_multipliers']};mu={v:F(a)for v,a in old['vertex_multipliers']}
    assert all(a>=0 for a in alpha.values()) and all(a>=0 for a in mu.values())
    row=gate['inequality'];assert row['positive_faces']==[[26410,2],[26880,1]]
    assert row['negative_faces']==[[580,2]] and row['vertex_capacity']==[[5572,3]]
    assert forests['vertices'][5572]['forest']==['A4','H3','H3','H3']
    assert len(forests['vertices'])==6558
    lhs=sum((c*alpha.get(f,F())for f,c in row['positive_faces']),F())-sum((c*alpha.get(f,F())for f,c in row['negative_faces']),F())
    before=mu.get(5572,F());assert before==0
    gap=lhs-3*before;assert gap==F(row['gap'])==F(3363761,10**9)>0
    delta=gap/3;updated=dict(mu);updated[5572]=before+delta
    assert lhs==3*updated[5572]
    assert all(updated.get(v,F())>=a for v,a in mu.items())
    assert {v for v in set(mu)|set(updated)if mu.get(v,F())!=updated.get(v,F())}=={5572}
    normalizer=sum(updated.values(),F());objective=F(old['objective_numerator']);bound=objective/normalizer
    assert F(old['normalizer'])==sum(mu.values(),F())
    assert F(old['bound'])==objective/F(old['normalizer']) and normalizer==F(old['normalizer'])+delta
    assert 12<bound<F(old['bound']) and F(gate['input_bound'])==F(old['bound'])
    paths=[Path(__file__)]+[ROOT/f'results/{tag}.json'for tag in (BASE,GATE,LOCAL,FORESTS)]
    paths += [ROOT/f'results/{GATE}.md',ROOT/f'results/{GATE}_review.md']
    return dict(status='PASS_EXACT_MONOTONE_FORCED_QUAD_VERTEX_MASS_REPAIR',dimension=13,
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        input_model=BASE,input_bound=old['bound'],bound=str(bound),objective_numerator=str(objective),normalizer=str(normalizer),
        face_multipliers=old['face_multipliers'],vertex_multipliers=[[v,str(a)]for v,a in sorted(updated.items())],
        face_masses_unchanged=True,vertex_updates=[[5572,str(before),str(updated[5572])]],
        added_vertex_mass=str(delta),new_incidence_row=row,repaired_incidence_load=str(lhs),repaired_vertex_capacity=str(3*updated[5572]),
        inherited_local_cover_source=LOCAL,inherited_joint_feature_source=BASE,
        preservation=dict(inherited_signed_and_cap_loads='identical because alpha is unchanged',
            complete_local_packing_marginals='identical; every saved disjoint packing mixture is retained with weakly larger vertex allowance',
            joint_3_4_5_family_loads='identical; all capacities are nonnegative multiples of mu',
            inherited_positive_cost_columns='unchanged face loads and weakly larger nonnegative-cost budget',
            old_inequality_count=2711432,new_incidence_count=1),
        optimizer_invocations=0,repeated_large_scans=0,excludes_target12=True,
        scope='Exact repair of one saved partial-model witness by increasing only mu5572. All inherited constraints persist by monotonicity; one new actual-face incidence inequality is added. No complete global incidence model, optimality or dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],float(F(data['bound'])),data['vertex_updates'],flush=True)
