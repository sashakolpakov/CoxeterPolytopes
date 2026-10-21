#!/usr/bin/env python3
"""Joint corner-orbit successor with every finite spectator context of five cores."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_joint_corner_repair1_20261017 as R

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017'
FAMILY='fresh_sqrt5_rank13_five_core_spectators_20261017'

def domain():
    d=R.domain();family=B.read(FAMILY);graphs=dict(B.C.models())
    lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index']for v in d['variables']}
    present={tuple(v['forest'])for v in d['vertices']}
    required={tuple(c['forest'])for r in family['contexts']for c in r['corners']}
    new_forests=sorted(required-present)
    for forest in new_forests:
        g={}
        for name in forest:
            h=graphs[name];off=len(g)
            for i,row in h.items():g[off+i]={off+j:m for j,m in row.items()}
        assert len(g)==13;pieces=B.C.pieces(g);orbits=defaultdict(list)
        for pair in combinations(range(13),2):
            code=B.dump(sorted(B.C.code(p,pair)for p in pieces));orbits[code].append(list(pair))
        vertex=len(d['vertices']);variables=[]
        for code,pairs in sorted(orbits.items()):
            index=len(d['variables']);lookup[forest,code]=index;variables.append([index,len(pairs)])
            d['variables'].append(dict(index=index,vertex=vertex,forest=list(forest),marked_forest_code=json.loads(code),
                multiplicity=len(pairs),literal_pairs=pairs))
        assert sum(m for _,m in variables)==78
        d['vertices'].append(dict(index=vertex,forest=list(forest),variables=variables))
    face_lookup={(f['required'],tuple(map(tuple,f['coefficients']))):f['index']for f in d['faces']}
    mappings=[]
    for context in family['contexts']:
        indices=[]
        for corner in context['corners']:
            index=lookup[tuple(corner['forest']),B.dump(corner['marked_forest_code'])]
            if corner['predecessor_variable']is not None:assert corner['predecessor_variable']==index
            indices.append(index)
        coeff=tuple(sorted(Counter(indices).items()));key=(1,coeff)
        provenance=dict(kind='complete_finite_spectator_family',source=FAMILY,context=context['index'])
        if key in face_lookup:d['faces'][face_lookup[key]]['provenance'].append(provenance)
        else:
            face_lookup[key]=len(d['faces'])
            d['faces'].append(dict(index=len(d['faces']),required=1,coefficients=[list(x)for x in coeff],provenance=[provenance]))
        mappings.append(dict(context=context['index'],face=face_lookup[key],corner_variables=indices))
    assert len(new_forests)==family['counts']['new_corner_forests']==1660
    d['counts'].update(vertices=len(d['vertices']),variables=len(d['variables']),
        literal_face_rows=d['counts']['literal_face_rows']+len(family['contexts']),distinct_face_constraints=len(d['faces']),
        complete_spectator_contexts=len(family['contexts']),added_whole_vertex_forests=len(new_forests))
    d['spectator_context_rows']=mappings
    for p in[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{R.TAG}.py',f'results/{R.TAG}.json',f'results/{FAMILY}.json']:
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
    status='PASS_EXACT_SPECTATOR_REPAIRED_PARTIAL_RULE_AT12'if feasible else'PASS_EXACT_SPECTATOR_REPAIRED_PARTIAL_MODEL_OBSTRUCTION'if excluded else'UNRESOLVED_SPECTATOR_REPAIRED_PARTIAL_MODEL'
    d.update(status=status,dimension=13,target='12',solver=solver,certificates=certificates,
        change='All682 finite orthogonal spectator contexts for five exactly certified compact local cores are included in one batch. All1660 previously absent whole vertex forests receive complete78-pair budgets and free nonnegative corner-orbit variables. All predecessor constraints are preserved.',
        scope='Exact partial-model successor only. The five fixed-core spectator families are complete, but general rank13 triangle/quad coverage and the remaining vertex types are not. A feasible payment vector is not a dimension bound; any exact dual applies only to this local model.')
    path.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound')))for c in d['certificates']],flush=True)
