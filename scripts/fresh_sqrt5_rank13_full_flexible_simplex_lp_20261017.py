#!/usr/bin/env python3
"""One enlarged generic-discriminant simplex model, before any large class batch."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017'
INPUT='fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'
PREVIOUS='fresh_sqrt5_rank13_discriminant_simplex_lp_20261017'

def domain():
    inp=B.read(INPUT);previous=B.read(PREVIOUS)
    models=dict(B.C.models());finite={B.C.code(g):name for name,g in models.items() if len(g)<=13}
    forests={tuple(v['forest']) for v in previous['vertices']};old_forests=set(forests);raw=[]
    for star in inp['stars']:
        if star['fixed_discriminant']:continue
        g={i:{} for i in range(14)}
        for i,j,m in star['edges']:g[i][j]=g[j][i]=m
        S=[s['omitted'] for s in star['spherical_deletions']];keys=[]
        for omitted in S:
            end={i:{j:m for j,m in row.items() if j!=omitted} for i,row in g.items() if i!=omitted}
            pieces=B.C.pieces(end);forest=tuple(sorted(finite[B.C.code(p)] for p in pieces));forests.add(forest)
            for pair in combinations([i for i in S if i!=omitted],2):keys.append((forest,B.dump(sorted(B.C.code(p,pair) for p in pieces))))
        assert len(keys)==3*comb(len(S),3)
        raw.append((2*comb(len(S),3),keys,dict(star=star['index'],simplex_dimension=len(S)-1,component_sizes=star['component_sizes'])))
    variables=[];vertices=[];lookup={}
    for forest in sorted(forests):
        g={}
        for name in forest:
            h=models[name];offset=len(g)
            for i,row in h.items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13;pieces=B.C.pieces(g);orbits=defaultdict(list)
        for pair in combinations(range(13),2):orbits[B.dump(sorted(B.C.code(p,pair) for p in pieces))].append(list(pair))
        vid=len(vertices);vv=[]
        for code,pairs in sorted(orbits.items()):
            index=len(variables);lookup[forest,code]=index;vv.append([index,len(pairs)])
            variables.append(dict(index=index,vertex=vid,forest=list(forest),marked_forest_code=json.loads(code),multiplicity=len(pairs),literal_pairs=pairs))
        assert sum(m for i,m in vv)==78
        vertices.append(dict(index=vid,forest=list(forest),variables=vv))
    faces=[];face_lookup={}
    for required,keys,provenance in raw:
        coeff=tuple(sorted(Counter(lookup[k] for k in keys).items()));key=(required,coeff)
        if key not in face_lookup:
            face_lookup[key]=len(faces);faces.append(dict(index=len(faces),required=required,coefficients=[list(x) for x in coeff],provenance=[]))
        faces[face_lookup[key]]['provenance'].append(provenance)
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'results/{INPUT}.json',f'results/{PREVIOUS}.json']
    return dict(sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},vertices=vertices,variables=variables,faces=faces,
        counts=dict(vertices=len(vertices),variables=len(variables),distinct_face_constraints=len(faces),literal_simplex_stars=len(raw),
            inherited_whole_forests=len(old_forests),new_whole_forests=len(forests-old_forests)),
        definition='All determinant-flexible maximal simplex stars in the complete197957-occurrence selected-endpoint triangle catalogue, with aggregate2-per-triangle demands and all78-pair budgets for the old605 forests plus every new literal corner forest.')

def run(verify):
    d=domain();B.dual_bound=D.sparse_dual_bound;out=ROOT/f'results/{TAG}.json'
    if verify:
        saved=json.loads(out.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:B.verify_certificate(d,c)
        return saved
    B.recover_dual=D.recover_unminimized;solver,certificates=B.optimize(d)
    feasible=any(c.get('reaches_target12',False) for c in certificates);excluded=any(c.get('excludes_target12',False) for c in certificates)
    assert not(feasible and excluded)
    d.update(status='PASS_EXACT_ENLARGED_FLEXIBLE_SIMPLEX_RULE_AT12' if feasible else 'PASS_EXACT_ENLARGED_FLEXIBLE_SIMPLEX_OBSTRUCTION' if excluded else 'UNRESOLVED_ENLARGED_FLEXIBLE_SIMPLEX_MODEL',
        solver=solver,certificates=certificates,
        scope='One partial generic-discriminant model only, with a300-second optimizer limit. Fixed-class stars, triangles outside the saved selected-endpoint domain, other simplex stars, quadrilaterals and absent vertex forests remain outside. Exact feasibility is not a dimension theorem; an exact dual only obstructs this whole-forest corner-payment family on the witnessed flexible stars.')
    out.write_text(json.dumps(d,indent=2)+'\n');return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['counts'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound'))) for c in d['certificates']],flush=True)
