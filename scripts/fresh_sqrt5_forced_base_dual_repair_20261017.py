#!/usr/bin/env python3
"""Exact analytic dual repair on the prepared forced-base domain."""
import argparse,json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_forced_base_simplex_domain_20261017 as N
import fresh_sqrt5_prism_shared_transfers_20261017 as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_forced_base_dual_repair_20261017'


def domain():
    patch=P.read(N.TAG);shared=P.read(S.TAG);d=P.read(N.I.LP)
    assert patch['inherited_model']==f'results/{N.I.LP}.json'
    assert patch['inherited_shared_transfers']==f'results/{S.TAG}.json'
    original=dict(faces=len(d['faces']),variables=len(d['variables']),vertices=len(d['vertices']))
    assert original==patch['inherited_counts']
    for var in patch['new_variables']:
        assert var['index']==len(d['variables']);d['variables'].append(var)
    for v in patch['new_vertices']:
        assert v['index']==len(d['vertices']);d['vertices'].append(v)
    for v,entries in patch['existing_vertex_budget_additions']:
        assert v<original['vertices'];d['vertices'][v]['variables'].extend(entries)
    for face in patch['new_faces']:
        assert face['index']==len(d['faces']);d['faces'].append(face)
    seen=set()
    for v in d['vertices']:
        for j,m in v['variables']:
            assert j not in seen;seen.add(j)
            var=d['variables'][j]
            assert var['index']==j and var['vertex']==v['index'] and var['multiplicity']==m and m>0
    assert seen==set(range(len(d['variables'])))
    for face in d['faces']:
        assert face['required']>=0
        assert len({j for j,m in face['coefficients']})==len(face['coefficients'])
        assert all(j in seen and m>0 for j,m in face['coefficients'])
    for rec in patch['touched_vertex_budgets']:
        actual=dict(d['vertices'][rec['vertex']]['variables'])
        for j,m in rec['pair_variables']+rec['port_variables']:assert actual[j]==m
        assert sum(m for j,m in rec['pair_variables'])==78
        assert sum(m for j,m in rec['port_variables'])==13
    transfers=[]
    for source,columns in [('represented',shared['transfer_columns']),('newly_forced',patch['new_shared_transfer_columns'])]:
        for c in columns:
            sf,negative=c['negative_simplex_coefficient'];capacity=c['shared_capacity_per_simplex']
            assert negative==-capacity and c['vertex_cost']==0 and capacity>0
            assert d['faces'][sf]['kind'] in ('maximal_simplex_aggregate','maximal_edge_incidence')
            assert all(d['faces'][i]['kind']=='canonical_dotted_quad_aggregate' and m>0
                       for i,m in c['positive_quad_coefficients'])
            transfers.append(dict(index=len(transfers),source=source,source_column=c['index'],
                base_type=c['base_type'],simplex_face=sf,simplex_dimension=c['simplex_dimension'],
                capacity=capacity,positive_quad_coefficients=c['positive_quad_coefficients']))
    d['active_variables']=sorted({j for f in d['faces'] for j,m in f['coefficients']})
    d['counts']=dict(vertices=len(d['vertices']),variables=len(d['variables']),faces=len(d['faces']),
        active_variables=len(d['active_variables']),zero_only_variables=len(d['variables'])-len(d['active_variables']),
        represented_transfers=len(shared['transfer_columns']),new_transfers=len(patch['new_shared_transfer_columns']),
        signed_zero_budget_transfers=len(transfers),new_edge_rows=patch['counts']['maximal_edge_rows'])
    assert d['counts']['signed_zero_budget_transfers']==5693
    d['transfers']=transfers
    return d,patch,shared


def signed_loads(d,alpha):
    answer=[]
    for c in d['transfers']:
        value=sum((m*alpha.get(i,F(0)) for i,m in c['positive_quad_coefficients']),F(0))
        value-=c['capacity']*alpha.get(c['simplex_face'],F(0))
        answer.append(value)
    return answer


def repair(d,alpha):
    answer=dict(alpha)
    for c in d['transfers']:
        need=sum((m*alpha.get(i,F(0)) for i,m in c['positive_quad_coefficients']),F(0))/c['capacity']
        sf=c['simplex_face']
        if need>answer.get(sf,F(0)):answer[sf]=need
    answer={i:a for i,a in answer.items() if a>0}
    assert all(x<=0 for x in signed_loads(d,answer))
    return answer


def certificate(d,alpha):
    bound,beta,tight,num,den=D.sparse_dual_bound(d,alpha)
    c=dict(kind='EXACT_RESTRICTED_MODEL_DUAL',bound=str(bound),objective_numerator=str(num),normalizer=str(den),
        face_multipliers=[[i,str(a)] for i,a in sorted(alpha.items())],
        vertex_multipliers=[[i,str(a)] for i,a in sorted(beta.items())],
        tight_orbits={str(i):a for i,a in tight.items()},excludes_target12=bound>12,
        all_signed_transfer_loads=list(map(str,signed_loads(d,alpha))))
    verify_dual(d,c);return c


def verify_dual(d,c):
    B.dual_bound=D.sparse_dual_bound;B.verify_certificate(d,c)
    loads=signed_loads(d,{i:F(a) for i,a in c['face_multipliers']})
    assert all(x<=0 for x in loads)
    assert list(map(str,loads))==c['all_signed_transfer_loads']


def run():
    d,patch,shared=domain();old=shared['repaired_dual']
    alpha={i:F(a) for i,a in old['face_multipliers']};repaired=repair(d,alpha)
    old_faces=patch['inherited_counts']['faces']
    assert all(repaired.get(i,F(0))==a for i,a in alpha.items())
    assert all(i>=old_faces for i in set(repaired)-set(alpha))
    changes=[]
    for f in patch['new_faces']:
        need=F(f['shared_dual_forced_minimum_mass'])
        assert repaired.get(f['index'],F(0))==need
        if need:changes.append(dict(face=f['index'],base_type=f['base_type'],
            simplex_dimension=f['simplex_dimension'],required=f['required'],mass=str(need)))
    c=certificate(d,repaired);old_mu={i:F(a) for i,a in old['vertex_multipliers']}
    increases=[[i,str(F(a)-old_mu.get(i,F(0)))] for i,a in c['vertex_multipliers'] if F(a)>old_mu.get(i,F(0))]
    assert all(F(a)>=old_mu.get(i,F(0)) for i,a in c['vertex_multipliers'])
    positive=defaultdict(int)
    for r in changes:positive[r['simplex_dimension']]+=1
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{P.TAG}.py',
        f'scripts/{N.TAG}.py',f'results/{N.TAG}.json',f'scripts/{S.TAG}.py',f'results/{S.TAG}.json',
        f'results/{N.I.LP}.json']
    return dict(status='PASS_EXACT_FORCED_BASE_ANALYTIC_DUAL_BARRIER' if F(c['bound'])>12 else
        'PASS_EXACT_FORCED_BASE_ANALYTIC_REPAIR_AT_MOST12',sources={p:P.sha(p) for p in paths},
        dimension=13,target=12,counts=dict(d['counts'],new_positive_face_masses=len(changes),
            positive_face_masses=len(c['face_multipliers']),positive_vertex_masses=len(c['vertex_multipliers']),
            increased_vertex_masses=len(increases)),old_shared_bound=old['bound'],repaired_dual=c,
        positive_new_faces_by_dimension=dict(sorted(positive.items())),new_face_masses=changes,
        additional_face_demand=str(F(c['objective_numerator'])-F(old['objective_numerator'])),
        additional_vertex_mass=str(F(c['normalizer'])-F(old['normalizer'])),vertex_mass_increases=increases,
        all_new_base_rows=f'results/{N.TAG}.json',
        repair_definition='Retain every represented shared-repaired simplex mass and every old quad mass. Add each missing edge/simplex mass at its exact minimum shared quad load divided by capacity. Maximal edges have demand0. Recompute every original and added positive-cost column; recompute all5693 original/new signed zero-budget transfer loads. No numerical optimizer is used in this analytic step.',
        scope='Only the inherited pair,22clique,old/newone-port and represented/newshared-transfer domain is covered. The separate mixed-three-port columns are excluded because their features on newly activated stars have not been reconstructed here. Other absent bases and face classes remain outside. A bound>12 obstructs this weighting family only; no dimension theorem or globally realizable incidence distribution is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['repaired_dual']['bound'],flush=True)
