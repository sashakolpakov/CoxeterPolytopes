#!/usr/bin/env python3
"""Joint corner-orbit successor with the complete all-selected non-ADE triangle batch."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_joint_corner_double_h4_repair_20261017 as R
import fresh_sqrt5_rank13_all_selected_nonade_separation_20261017 as E

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_joint_corner_all_nonade_repair_20261017'

def domain():
    d=R.domain();tri=B.read(E.TRI);graphs=dict(B.C.models())
    finite={B.C.code(g):name for name,g in graphs.items()if len(g)<=13}
    lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index']for v in d['variables']}
    present={tuple(v['forest'])for v in d['vertices']};required=set();raw=[]
    for ti,target in enumerate(tri['targets']):
        for oi,occurrence in enumerate(target['triangle_occurrences']):
            g,black=E.literal(target,occurrence);keys=[]
            for removed in [13,occurrence['port'],occurrence['other']]:
                pair=[i for i in black if i!=removed]
                end={i:{j:m for j,m in row.items()if j!=removed}for i,row in g.items()if i!=removed}
                pieces=B.C.pieces(end);forest=tuple(sorted(finite[B.C.code(p)]for p in pieces))
                code=B.dump(sorted(B.C.code(p,pair)for p in pieces));keys.append((forest,code));required.add(forest)
            assert keys[0][0]==tuple(target['forest'])and keys[1][0]==tuple(occurrence['destination'])
            assert keys[2]==(tuple(occurrence['third_forest']),B.dump(occurrence['third_marked_forest_code']))
            raw.append((keys,dict(kind='all_selected_nonade_triangle',source=E.TRI,target=ti,occurrence=oi)))
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
    for keys,provenance in raw:
        indices=[lookup[k]for k in keys];coeff=tuple(sorted(Counter(indices).items()));key=(2,coeff)
        if key in face_lookup:d['faces'][face_lookup[key]]['provenance'].append(provenance)
        else:
            face_lookup[key]=len(d['faces'])
            d['faces'].append(dict(index=len(d['faces']),required=2,coefficients=[list(x)for x in coeff],provenance=[provenance]))
        mappings.append(dict(target=provenance['target'],occurrence=provenance['occurrence'],face=face_lookup[key],corner_variables=indices))
    d['counts'].update(vertices=len(d['vertices']),variables=len(d['variables']),
        literal_face_rows=d['counts']['literal_face_rows']+len(raw),distinct_face_constraints=len(d['faces']),
        all_selected_nonade_triangle_occurrences=len(raw),all_selected_nonade_added_whole_vertex_forests=len(new_forests))
    d['all_selected_nonade_triangle_rows']=mappings
    for p in[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{R.TAG}.py',f'scripts/{E.TAG}.py',
             f'results/{R.TAG}.json',f'results/{E.TAG}.json',f'results/{E.TRI}.json']:
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
    status='PASS_EXACT_ALL_NONADE_REPAIRED_PARTIAL_RULE_AT12'if feasible else'PASS_EXACT_ALL_NONADE_REPAIRED_PARTIAL_MODEL_OBSTRUCTION'if excluded else'UNRESOLVED_ALL_NONADE_REPAIRED_PARTIAL_MODEL'
    d.update(status=status,dimension=13,target='12',solver=solver,certificates=certificates,
        change='All literal triangle occurrences at the363 selected non-ADE targets are included together. Every previously absent corner forest receives all78-pair budgets and free nonnegative corner-orbit variables. All predecessor constraints and variable indices are preserved.',
        scope='Exact partial-model successor only. The new triangle domain is complete for all363 non-ADE starting types and the stated578-type selected endpoint set. General rank13 face coverage and absent finite vertex types remain outstanding. A feasible payment vector is not a dimension bound; any exact dual applies only to this local model.')
    path.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound')))for c in d['certificates']],flush=True)
