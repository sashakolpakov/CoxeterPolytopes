#!/usr/bin/env python3
"""Shared prism-base capacity: exact zero-cost export and analytic dual repair."""
import argparse,json
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_prism_base_transfers_20261017 as T

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_prism_shared_transfers_20261017'
PROOF='results/fresh_prism_simplex_shared_capacity_20261017.md'


def run():
    prior=P.read(T.TAG);d=P.read(T.I.LP)
    old=next(c for c in d['certificates'] if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    B.dual_bound=D.sparse_dual_bound;B.verify_certificate(d,old)
    alpha={i:F(a) for i,a in old['face_multipliers']}
    groups=defaultdict(list)
    for col in prior['transfer_columns']:groups[col['base_type']].append(col)
    columns=[]
    for t,group in sorted(groups.items()):
        dimension=group[0]['simplex_dimension'];sf=group[0]['simplex_face']
        assert all(c['simplex_dimension']==dimension and c['simplex_face']==sf for c in group)
        capacity=(dimension+1)*(13-dimension);coeff=Counter()
        for col in group:
            k=col['prism_dimension'];assert 2<=k<=dimension+1
            assert k*col['capacity_per_simplex']>=capacity
            for fi,m in col['positive_quad_coefficients']:coeff[fi]+=k*m
        assert all(fi>=10185 and n>0 for fi,n in coeff.items())
        columns.append(dict(index=len(columns),base_type=t,old_star=group[0]['old_star'],
            simplex_face=sf,simplex_dimension=dimension,vertex_cost=0,
            shared_capacity_per_simplex=capacity,positive_quad_coefficients=[list(x) for x in sorted(coeff.items())],
            negative_simplex_coefficient=[sf,-capacity],separate_columns=[c['index'] for c in group]))
    repaired=dict(alpha);requirements=[];violations=[]
    for col in columns:
        load=sum((m*alpha.get(fi,F(0)) for fi,m in col['positive_quad_coefficients']),F(0))
        sf=col['simplex_face'];capacity=col['shared_capacity_per_simplex'];minimum=load/capacity
        if minimum>repaired.get(sf,F(0)):repaired[sf]=minimum
        slack=capacity*alpha.get(sf,F(0))-load
        if slack<0:violations.append(col['index'])
        if load:requirements.append([col['index'],str(load),str(minimum),str(slack)])
    repaired={i:a for i,a in repaired.items() if a>0}
    shared_loads=[];separate_loads=[]
    for col in columns:
        load=sum((m*repaired.get(fi,F(0)) for fi,m in col['positive_quad_coefficients']),F(0))
        load-=col['shared_capacity_per_simplex']*repaired.get(col['simplex_face'],F(0))
        assert load<=0;shared_loads.append(str(load))
    for col in prior['transfer_columns']:
        load=sum((m*repaired.get(fi,F(0)) for fi,m in col['positive_quad_coefficients']),F(0))
        load-=col['capacity_per_simplex']*repaired.get(col['simplex_face'],F(0))
        assert load<=0;separate_loads.append(str(load))
    bound,beta,tight,num,den=D.sparse_dual_bound(d,repaired)
    cert=dict(kind='EXACT_RESTRICTED_MODEL_DUAL',bound=str(bound),objective_numerator=str(num),normalizer=str(den),
        face_multipliers=[[i,str(a)] for i,a in sorted(repaired.items())],
        vertex_multipliers=[[i,str(a)] for i,a in sorted(beta.items())],
        tight_orbits={str(i):x for i,x in tight.items()},excludes_target12=bound>12)
    B.verify_certificate(d,cert)
    changes=[]
    for sf,a in sorted(repaired.items()):
        previous=alpha.get(sf,F(0))
        if a!=previous:
            assert sf<10185 and a>previous
            changes.append(dict(simplex_face=sf,old_mass=str(previous),new_mass=str(a),increase=str(a-previous),
                                star=d['faces'][sf]['provenance'][0]['star']))
    assert all(repaired.get(fi,F(0))==a for fi,a in alpha.items() if fi>=10185)
    separate_repair={i:F(a) for i,a in prior['repaired_dual']['face_multipliers']}
    assert all(repaired.get(i,F(0))>=a for i,a in separate_repair.items())
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{P.TAG}.py',
        f'scripts/{T.TAG}.py',f'results/{T.TAG}.json',f'results/{T.I.LP}.json',PROOF]
    return dict(status='PASS_EXACT_SHARED_PRISM_BASE_TRANSFERS_AND_ANALYTIC_REPAIR',sources={p:P.sha(p) for p in paths},
        dimension=13,target=12,counts=dict(shared_columns=len(columns),
            positive_quad_entries=sum(len(c['positive_quad_coefficients']) for c in columns),
            dominated_separate_columns=len(prior['transfer_columns']),old_violated_shared_columns=len(violations),
            increased_simplex_masses=len(changes),newly_positive_simplex_rows=sum(F(c['old_mass'])==0 for c in changes),
            repaired_positive_face_rows=len(repaired),repaired_positive_vertex_masses=len(beta)),
        transfer_columns=columns,old_bound=old['bound'],separate_repair_bound=prior['repaired_dual']['bound'],
        original_violated_shared_columns=violations,original_positive_load_requirements=requirements,
        simplex_mass_changes=changes,repaired_dual=cert,all_shared_dual_loads=shared_loads,
        all_separate_dual_loads=separate_loads,
        additional_face_demand=str(num-F(old['objective_numerator'])),
        additional_vertex_mass=str(den-F(old['normalizer'])),
        provenance_rule='For each quad row its prism dimension k is constant across every provenance. Thus multiplying the already verified minimum base multiplicity by k gives the minimum shared coefficient. Group by represented base type and sum across k; keep every provenance reference in the frozen input. Missing base types remain unknown and are not assigned zero mass.',
        dominance='k*C(d+1,k)=(d+1)*C(d,k-1)>=d+1. Shared load<= (d+1)(13-d)*alphaS implies every separate-k load<=C(d+1,k)(13-d)*alphaS. The exact repaired vector is explicitly checked against both exports.',
        scope='No optimizer. Original quad masses are retained and original simplex masses are increased only as necessary for represented shared capacities. Every inherited pair/clique/disjoint-port column is rechecked, and all new zero-cost signed loads are nonpositive. The separate mixed-three-port columns are not included, nor any missing base type or unenumerated face. A bound above12 is a partial-family obstruction; a bound below12 would only fail to retain this particular obstruction, not establish feasibility or a dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['repaired_dual']['bound'],flush=True)
