#!/usr/bin/env python3
"""Exact union of vertex-labelled base capacity and integer port packing."""
import argparse,json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_forced_base_dual_repair_20261017 as R
import fresh_sqrt5_actual_provenance_dual_20261017 as A
import fresh_sqrt5_vertex_refined_base_capacity_20261017 as V

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_vertex_refined_integer_union_20261017'
PACKING='fresh_sqrt5_integer_simplex_port_packing_20261017'


def run():
    vertex=R.P.read(V.TAG);packing=R.P.read(PACKING);actual=R.P.read(A.TAG)
    d,patch,shared=R.domain()
    alpha={i:F(a) for i,a in vertex['face_multipliers']}
    oldmu={v:F(a) for v,a in vertex['vertex_multipliers']}
    assert set(alpha)=={i for i,a in packing['face_multipliers']}=={i for i,a in actual['face_multipliers']}
    assert set(oldmu)=={v for v,a in packing['vertex_multipliers']}
    pfaces={r['face']:r for r in packing['positive_simplex_features']}
    quads={i for i in alpha if d['faces'][i]['kind']=='canonical_dotted_quad_aggregate'}
    assert set(pfaces).isdisjoint(quads) and set(pfaces)|quads==set(alpha)
    assert len(pfaces)==5462 and len(quads)==1582
    pc=packing['columns'];mixed=actual['triple_orbit_columns']
    for c,p in zip(mixed,pc):
        assert c['index']==p['index'] and c['vertex']==p['vertex']
        assert c['vertex_cost']==p['vertex_cost'] and c['marked_forest_code']==p['marked_forest_code']
    assert len(pc)==len(mixed)+3*len(oldmu)==62353
    loads=defaultdict(F)
    for fi,rec in pfaces.items():
        assert d['faces'][fi]['required'] in (0,2,8,20)
        for j,m in rec['coefficients']:
            assert m>0;loads[j]+=alpha[fi]*m
    mu=dict(oldmu);violations=[]
    for c in pc:
        j,v,cost=c['index'],c['vertex'],c['vertex_cost'];need=loads[j]/cost
        if need>oldmu[v]:violations.append(dict(column=j,vertex=v,kind=c['kind'],old_vertex_mass=str(oldmu[v]),
            load=str(loads[j]),capacity=cost,required_vertex_mass=str(need)))
        mu[v]=max(mu[v],need)
    # Recheck every old nonnegative column, not merely the newly active ones.
    positive=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in d['faces'][fi]['coefficients']:positive[j]+=a*m
    for c in d['variables']:assert positive[c['index']]<=c['multiplicity']*mu.get(c['vertex'],F(0))
    features={r['face']:r for r in actual['positive_face_mixed_features']}
    assert set(features)==set(alpha)
    mloads=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in features[fi]['coefficients']:mloads[j]+=a*m
    for c in mixed:assert mloads[c['index']]<=c['vertex_cost']*mu[c['vertex']]
    for c in pc:assert loads[c['index']]<=c['vertex_cost']*mu[c['vertex']]
    inherited=R.signed_loads(d,alpha);assert all(a<=0 for a in inherited)
    shared_loads=[]
    for c in actual['actual_shared_base_constraints']:
        load=sum((m*F(actual['provenance_masses'][j]['mass']) for j,m in c['positive_provenance_coefficients']),F(0))
        value=load-c['capacity']*alpha[c['simplex_face']]
        assert load==F(c['actual_positive_load']) and value<=0
        shared_loads.append([c['base_type'],str(value)])
    labelled=[]
    for c in vertex['labelled_transfer_columns']:
        load=sum((m*F(actual['provenance_masses'][j]['mass']) for j,m in c['provenance_coefficients']),F(0))
        value=load-c['capacity']*alpha[c['simplex_face']]
        assert load==F(c['actual_load']) and value==F(c['repaired_signed_load'])<=0
        labelled.append(str(value))
    # Explicit provenance masses are unchanged under both repairs.
    for row in actual['merged_row_mass_constraints']:
        mass=sum((F(actual['provenance_masses'][j]['mass']) for j in row['provenance_indices']),F(0))
        assert mass==alpha[row['face']]==F(row['mass'])
    num=sum((d['faces'][i]['required']*a for i,a in alpha.items()),F(0))
    assert num==F(vertex['objective_numerator'])
    den=sum(mu.values(),F(0));old_den=sum(oldmu.values(),F(0));bound=num/den
    changes=[[v,str(oldmu[v]),str(a),str(a-oldmu[v])] for v,a in sorted(mu.items()) if a>oldmu[v]]
    paths=[f'scripts/{TAG}.py',f'scripts/{V.TAG}.py',f'results/{V.TAG}.json',f'results/{V.TAG}.md',
        f'scripts/{PACKING}.py',f'results/{PACKING}.json',f'results/{PACKING}.md',
        f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',f'scripts/{R.TAG}.py']
    return dict(status='PASS_EXACT_LABELLED_BASE_INTEGER_PACKING_UNION',sources={p:R.P.sha(p) for p in paths},
        dimension=13,target=12,bound=str(bound),excludes_target12=bound>12,
        vertex_refined_bound=vertex['bound'],integer_packing_separate_bound=packing['bound'],
        objective_numerator=str(num),normalizer=str(den),additional_vertex_mass=str(den-old_den),
        counts=dict(positive_rows=len(alpha),simplex_rows=len(pfaces),quad_rows=len(quads),
            positive_forests=len(mu),positive_cost_columns=len(d['variables']),mixed_columns=len(mixed),
            integer_packing_columns=len(pc),actual_provenances=len(actual['provenance_masses']),
            labelled_signed_columns=len(labelled),actual_unlabelled_signed_columns=len(shared_loads),
            inherited_signed_columns=len(inherited),violating_packing_columns=len(violations),
            increased_vertex_masses=len(changes)),
        face_multipliers=vertex['face_multipliers'],vertex_multipliers=[[v,str(a)] for v,a in sorted(mu.items())],
        packing_column_loads=[[j,str(loads[j])] for j in range(len(pc))],packing_violations=violations,
        vertex_mass_changes=changes,mixed_column_loads=[[j,str(a)] for j,a in sorted(mloads.items()) if a],
        labelled_signed_loads=labelled,actual_shared_signed_loads=shared_loads,
        inherited_signed_loads=list(map(str,inherited)),
        scope='One exact union check at the vertex-refined alpha and frozen actual uniform provenance split. All integer-packing loads are recomputed at this alpha; the two separate normalizers are not combined heuristically. A ratio at most12 removes only this fixed-alpha obstruction and is not a feasible primal, optimum, complete face coverage or dimension bound. No optimizer is invoked.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(R.B.dump(d)+'\n')
    print(d['status'],d['counts'],d['bound'],flush=True)
