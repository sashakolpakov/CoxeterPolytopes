#!/usr/bin/env python3
"""One new combined model with the proved disjoint-simplex-port charges."""
import argparse,json
from collections import defaultdict
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
PORTS='fresh_simplex_disjoint_ports_20261017'


def domain():
    d=P.read(P.TAG);ports=P.read(PORTS)
    assert d['counts']['vertices']==6558 and len(d['faces'])==78335
    assert d['counts']['old_simplex_rows']==10185
    shared=f'results/{P.OLD}.json'
    assert d['sources'][shared]==ports['sources'][shared]
    original_counts=dict(d['counts']);old_n=len(d['variables'])
    assert old_n==203587
    models=dict(B.C.models());literal={}
    # Reconstruct all thirteen one-marked-root codes at each old forest. The
    # exporter and combined model use the same canonical component ordering.
    for old_v,new_v in enumerate(d['old_vertex_to_new']):
        forest=d['vertices'][new_v]['forest'];g={}
        for name in forest:
            offset=len(g)
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13
        pieces=B.C.pieces(g);orbits=defaultdict(list)
        for port in g:orbits[B.dump(sorted(B.C.code(p,(port,)) for p in pieces))].append(port)
        assert sum(map(len,orbits.values()))==13
        literal[old_v]=orbits
    columns=[];skipped=[];old_row_sums=[0]*10185
    for col in ports['successor_columns']:
        old_v=col['vertex'];v=d['old_vertex_to_new'][old_v]
        assert col['forest']==d['vertices'][v]['forest']
        code=B.dump(col['marked_forest_code'])
        assert literal[old_v][code]==col['literal_port_orbit']
        assert col['vertex_cost']==len(col['literal_port_orbit'])>0
        if not col['face_coefficients']:
            skipped.append(col['index']);continue
        j=len(d['variables']);cost=col['vertex_cost']
        d['variables'].append(dict(index=j,vertex=v,forest=col['forest'],
            kind='maximal_simplex_disjoint_port_charge',multiplicity=cost,
            source_column=col['index'],marked_forest_code=col['marked_forest_code'],
            literal_port_orbit=col['literal_port_orbit'],normalization='one unit per used port; vertex cost orbit size'))
        d['vertices'][v]['variables'].append([j,cost])
        coeff=[]
        for fi,n in col['face_coefficients']:
            assert 0<=fi<10185 and isinstance(n,int) and n>0
            face=d['faces'][fi]
            assert face['index']==face['original_face']==fi
            assert face['kind']=='maximal_simplex_aggregate'
            face['coefficients'].append([j,n]);coeff.append([fi,n]);old_row_sums[fi]+=n
        assert len({i for i,n in coeff})==len(coeff)
        columns.append(dict(variable=j,source_column=col['index'],old_vertex=old_v,
                            vertex=v,vertex_cost=cost,face_coefficients=coeff))
    assert len(columns)==ports['counts']['nonzero_columns']==6528
    assert len(columns)+len(skipped)==ports['counts']['columns']==15119
    assert sum(len(c['face_coefficients']) for c in columns)==54749
    for fi,face in enumerate(d['faces'][:10185]):
        assert len(face['provenance'])==1
        dimension=face['provenance'][0]['simplex_dimension']
        assert old_row_sums[fi]==dimension*(dimension+1)
    for face in d['faces'][10185:]:
        assert face['kind']=='canonical_dotted_quad_aggregate'
        assert all(j<old_n for j,n in face['coefficients'])
    pair_n=d['counts']['pair_variables']
    for vertex in d['vertices']:
        assert sum(m for j,m in vertex['variables'] if j<pair_n)==78
    d['active_variables']=sorted({j for face in d['faces'] for j,n in face['coefficients']})
    d['counts'].update(variables=len(d['variables']),active_variables=len(d['active_variables']),
        zero_only_variables=len(d['variables'])-len(d['active_variables']),
        disjoint_port_variables=len(columns),disjoint_port_face_coefficients=54749,
        omitted_zero_only_port_columns=len(skipped))
    d['disjoint_port_columns']=columns;d['omitted_zero_only_port_columns']=skipped
    d['inherited_counts']=original_counts
    d['sources'].update({p:P.sha(p) for p in [f'scripts/{TAG}.py',f'results/{P.TAG}.json',
        f'scripts/{PORTS}.py',f'results/{PORTS}.json']})
    for key in ('status','solver','certificates','scope'):d.pop(key,None)
    d['definition']+=f' Add all6528 active exported one-port charges, with literal port count coefficients and orbit-size vertex costs. Recheck all13 port codes and the old forest/row maps. All previous variables, face rows, budget entries and22 simplex packing charges remain; all quadrilateral coefficients of the new columns are zero. Omit only the{len(skipped)} exported columns absent from every face row.'
    return d


def optimize(d):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    active=d['active_variables'];lookup={i:j for j,i in enumerate(active)};N=len(active)
    rows=[];cols=[];values=[];rhs=[]
    for face in d['faces']:
        r=len(rhs);rhs.append(-face['required'])
        for j,m in face['coefficients']:rows.append(r);cols.append(lookup[j]);values.append(-m)
    for vertex in d['vertices']:
        r=len(rhs);rhs.append(0)
        for j,m in vertex['variables']:
            if j in lookup:rows.append(r);cols.append(lookup[j]);values.append(m)
        rows.append(r);cols.append(N);values.append(-1)
    matrix=coo_matrix((values,(rows,cols)),shape=(len(rhs),N+1)).tocsr()
    objective=np.zeros(N+1);objective[-1]=1
    print('DISJOINT_PORT_SINGLE_IPM_START',d['counts'],flush=True)
    result=linprog(objective,A_ub=matrix,b_ub=np.array(rhs,dtype=float),bounds=(0,None),
                   method='highs-ipm',options={'time_limit':300})
    solver=dict(status=int(result.status),message=str(result.message),success=bool(result.success),
                objective=None if result.fun is None else float(result.fun),
                method='highs-ipm',time_limit_seconds=300)
    print('DISJOINT_PORT_SINGLE_IPM_END',solver,flush=True)
    certificates=[]
    if result.success:
        x=np.zeros(len(d['variables']))
        for i,j in lookup.items():x[i]=result.x[j]
        primal=B.recover_primal(d,x)
        if primal is not None:certificates.append(primal)
        certificates.append(D.recover_unminimized(d,result.ineqlin.marginals[:len(d['faces'])]))
    B.dual_bound=D.sparse_dual_bound
    for c in certificates:B.verify_certificate(d,c)
    return solver,certificates


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=domain();out=ROOT/f'results/{TAG}.json';B.dual_bound=D.sparse_dual_bound
    if a.verify:
        saved=json.loads(out.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:B.verify_certificate(d,c)
        d=saved
    else:
        d.update(status='PREPARED_DISJOINT_PORT_COMBINED_DOMAIN',solver=None,certificates=[])
        out.write_text(B.dump(d)+'\n')
        solver,certificates=optimize(d)
        feasible=any(c.get('reaches_target12',False) for c in certificates)
        excluded=any(c.get('excludes_target12',False) for c in certificates)
        assert not(feasible and excluded)
        d.update(status='PASS_EXACT_COMBINED_DISJOINT_PORT_RULE_AT12' if feasible else
                 'PASS_EXACT_COMBINED_DISJOINT_PORT_OBSTRUCTION' if excluded else
                 'UNRESOLVED_COMBINED_DISJOINT_PORT_MODEL',solver=solver,certificates=certificates,
            scope='One300-second highs-ipm solve of a genuinely enlarged model with6528 proved one-port simplex charges. Every inherited combined row, variable and budget is preserved. Fixed-determinant stars, other connected cores, all-finite quad stars and general remaining faces are still absent;231 untouched finite forests remain outside the budget domain. An exact rule is partial feasibility only, and a dual only obstructs this charge family. No dimension theorem or simultaneous compact ambient realization is asserted.')
        out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound'))) for c in d['certificates']],flush=True)


if __name__=='__main__':main()
