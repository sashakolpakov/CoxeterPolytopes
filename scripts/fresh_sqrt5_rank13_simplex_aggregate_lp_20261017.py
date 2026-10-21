#!/usr/bin/env python3
"""Partial corner-orbit model with complete maximal-simplex triangle aggregates."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_simplex_aggregate_lp_20261017'
INPUT='fresh_sqrt5_rank13_triangle_simplex_flags_20261017'

def domain():
    inp=B.read(INPUT);models=dict(B.C.models());finite={B.C.code(g):name for name,g in models.items() if len(g)<=13}
    raw=[];forests=set()
    for record in inp['records']:
        g={i:{} for i in range(14)}
        for i,j,m in record['edges']:g[i][j]=g[j][i]=m
        S=[v['omitted'] for v in record['spherical_deletions']];keys=[]
        for omitted in S:
            end={i:{j:m for j,m in row.items() if j!=omitted} for i,row in g.items() if i!=omitted}
            pieces=B.C.pieces(end);forest=tuple(sorted(finite[B.C.code(p)] for p in pieces));forests.add(forest)
            for pair in combinations([i for i in S if i!=omitted],2):
                keys.append((forest,B.dump(sorted(B.C.code(p,pair) for p in pieces))))
        assert len(keys)==3*len(record['triangle_root_triples'])
        raw.append((record['aggregate_required_payment'],keys,dict(original_face=record['face'],simplex_dimension=len(S)-1,
            edges=record['edges'],spherical_deletions=S)))
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
            j=len(variables);lookup[forest,code]=j;vv.append([j,len(pairs)])
            variables.append(dict(index=j,vertex=vid,forest=list(forest),marked_forest_code=json.loads(code),
                multiplicity=len(pairs),literal_pairs=pairs))
        assert sum(m for j,m in vv)==78
        vertices.append(dict(index=vid,forest=list(forest),variables=vv))
    faces=[];face_lookup={}
    for required,keys,provenance in raw:
        coefficients=tuple(sorted(Counter(lookup[k] for k in keys).items()));key=(required,coefficients)
        if key not in face_lookup:
            face_lookup[key]=len(faces)
            faces.append(dict(index=len(faces),required=required,coefficients=[list(x) for x in coefficients],provenance=[]))
        faces[face_lookup[key]]['provenance'].append(provenance)
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'results/{INPUT}.json']
    return dict(sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        variables=variables,vertices=vertices,faces=faces,
        counts=dict(vertices=len(vertices),variables=len(variables),distinct_face_constraints=len(faces),literal_simplex_stars=len(raw)),
        definition='For every retained maximal simplex star, sum the payments at all corners of all its triangular2-faces, requiring2 per triangle only in aggregate. All78 pairs of each represented whole rank13 forest remain in its nonnegative vertex budget.')

def run(verify):
    d=domain();B.dual_bound=D.sparse_dual_bound;out=ROOT/f'results/{TAG}.json'
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
    d.update(status='PASS_EXACT_PARTIAL_SIMPLEX_AGGREGATE_RULE_AT12' if feasible else
        'PASS_EXACT_PARTIAL_SIMPLEX_AGGREGATE_OBSTRUCTION' if excluded else 'UNRESOLVED_PARTIAL_SIMPLEX_AGGREGATE_MODEL',
        solver=solver,certificates=certificates,
        scope='Only maximal simplex extensions of the saved positive-dual triangle representatives are imposed. General maximal-simplex domains, quadrilaterals, and all other finite vertex forests are not covered. Any feasible rule is a structural partial-model advance only; no dimension theorem or simultaneous realization is claimed.')
    out.write_text(json.dumps(d,indent=2)+'\n');return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['counts'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound'))) for c in d['certificates']],flush=True)
