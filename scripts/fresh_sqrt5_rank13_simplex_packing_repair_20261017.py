#!/usr/bin/env python3
"""One augmented rank13 test with22 proved maximal-simplex clique charges."""
import argparse,json
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_simplex_packing_repair_20261017'
INPUT='fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017'
PACK='fresh_simplex_port_incidence_20261017'


def sha(path):
    import hashlib
    h=hashlib.sha256()
    with (ROOT/path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def domain():
    old=B.read(INPUT);pack=B.read(PACK)
    assert len(old['vertices'])==1870 and len(old['faces'])==10185
    assert len(pack['successor_columns'])==22
    variables=old['variables'];vertices=old['vertices'];faces=old['faces']
    old_variables=len(variables)
    for v in vertices:assert sum(m for j,m in v['variables'])==78
    columns=[]
    for col in pack['successor_columns']:
        v=col['vertex'];assert col['forest']==vertices[v]['forest']
        assert F(col['vertex_cost'])==1
        coeff=[(i,F(a)) for i,a in col['face_coefficients']]
        assert len({i for i,a in coeff})==len(coeff) and all(a>0 for i,a in coeff)
        denominator=lcm(*(a.denominator for i,a in coeff))
        j=len(variables)
        variables.append(dict(index=j,vertex=v,forest=col['forest'],
            kind='maximal_simplex_clique_charge',multiplicity=denominator,
            source_column=col['index'],unscaled_vertex_cost='1',
            charge_scaling=denominator,subset=col['subset'],
            minimum_block_size=col['minimum_block_size'],subset_orbit=col['subset_orbit']))
        vertices[v]['variables'].append([j,denominator])
        scaled=[]
        for i,a in coeff:
            n=a*denominator;assert n.denominator==1
            faces[i]['coefficients'].append([j,n.numerator]);scaled.append([i,n.numerator])
        columns.append(dict(variable=j,source_column=col['index'],vertex=v,
            denominator=denominator,integer_vertex_cost=denominator,
            original_face_coefficients=col['face_coefficients'],integer_face_coefficients=scaled))
    for v in vertices:
        assert sum(m for j,m in v['variables'] if j<old_variables)==78
    assert len(variables)==old_variables+22
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',
           f'results/{INPUT}.json',f'scripts/{PACK}.py',f'results/{PACK}.json']
    return dict(sources={p:sha(p) for p in paths},dimension=13,target=12,
        vertices=vertices,variables=variables,faces=faces,packing_columns=columns,
        counts=dict(vertices=len(vertices),variables=len(variables),old_pair_variables=old_variables,
            packing_variables=len(columns),distinct_face_constraints=len(faces),
            added_face_rows=0,packing_nonzero_coefficients=sum(len(c['integer_face_coefficients']) for c in columns)),
        definition='Retain all1870 complete78-pair vertex budgets and all10185 flexible-star simplex aggregate rows. Add exactly22 nonnegative automorphism-averaged maximal-simplex clique charges, each with proved unscaled vertex cost1. Each new column is scaled independently by its exact common denominator; physical charge equals denominator times its LP variable. No old coefficient or constraint is removed.')


def run(verify):
    d=domain();out=ROOT/f'results/{TAG}.json'
    B.dual_bound=D.sparse_dual_bound
    if verify:
        saved=json.loads(out.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:B.verify_certificate(d,c)
        return saved
    B.recover_dual=D.recover_unminimized
    solver,certificates=B.optimize(d)
    feasible=any(c.get('reaches_target12',False) for c in certificates)
    excluded=any(c.get('excludes_target12',False) for c in certificates)
    assert not(feasible and excluded)
    d.update(status='PASS_EXACT_PARTIAL_SIMPLEX_PACKING_RULE_AT12' if feasible else 'PASS_EXACT_PARTIAL_SIMPLEX_PACKING_OBSTRUCTION' if excluded else 'UNRESOLVED_PARTIAL_SIMPLEX_PACKING_MODEL',
        solver=solver,certificates=certificates,
        scope='One300-second augmented partial model. The22 new charges use a proved global maximal-simplex port incidence restriction. Only the inherited flexible-star simplex domain is tested; fixed-class stars, other missing simplex stars, quadrilaterals and absent whole forests remain outside. No complete coverage, compact ambient realization or dimension bound is claimed.')
    out.write_text(json.dumps(d,indent=2)+'\n');return d


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['counts'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound'))) for c in d['certificates']],flush=True)
