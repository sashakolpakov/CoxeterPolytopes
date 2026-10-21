#!/usr/bin/env python3
"""One bounded combined flexible-simplex/dotted-quad aggregate model in rank13."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_simplex_quad_combined_20261017'
OLD='fresh_sqrt5_rank13_simplex_packing_repair_20261017'
QUADS='fresh_sqrt5_rank13_packing_quad_separation_20261017'
SPECTATORS='fresh_sqrt5_rank13_dotted_b1_spectators_20261017'


def sha(path):
    h=hashlib.sha256()
    with (ROOT/path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for path,h in d['sources'].items():assert sha(path)==h,path
    return d


def domain():
    old=read(OLD);quad=read(QUADS);spdata=read(SPECTATORS)
    assert len(old['vertices'])==1870 and len(old['faces'])==10185
    assert len(old['packing_columns'])==22
    models=[(name,g) for name,g in B.C.models() if len(g)<=13]
    model=dict(models);ranks={name:len(g) for name,g in models}
    codes={name:json.loads(B.dump(B.C.code(g))) for name,g in models}
    forests={tuple(v['forest']) for v in old['vertices']};old_forests=set(forests)
    raw=[]
    for source_context in quad['representative_context_indices']:
        ctx=spdata['contexts'][source_context];assert ctx['index']==source_context
        core=quad['canonical_cores'][quad['source_core_maps'][ctx['core']]['canonical_core']]
        spectator=ctx['spectator_forest']
        assert sum(ranks[n] for n in spectator)+core['nodes']==15
        if core['quotient_even'] and all(ranks[n]%2==0 for n in spectator):continue
        keys=[]
        for corner in core['aggregate_corner_types']:
            forest=tuple(sorted(corner['forest']+spectator))
            assert sum(ranks[n] for n in forest)==13
            code=B.dump(sorted(corner['marked_forest_code']+[codes[n] for n in spectator]))
            keys.append((forest,code,corner['multiplicity']));forests.add(forest)
        assert sum(m for f,c,m in keys)==4*core['required']
        raw.append((source_context,core['index'],core['required'],keys))
    assert len(raw)==quad['counts']['canonical_flexible']==81183
    del spdata,quad
    variables=[];vertices=[];lookup={};vertex_lookup={}
    for forest in sorted(forests):
        g={}
        for name in forest:
            h=model[name];offset=len(g)
            for i,row in h.items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13
        pieces=B.C.pieces(g);orbits=defaultdict(list)
        for pair in combinations(range(13),2):
            orbits[B.dump(sorted(B.C.code(p,pair) for p in pieces))].append(list(pair))
        vid=len(vertices);vertex_lookup[forest]=vid;vv=[]
        for code,pairs in sorted(orbits.items()):
            index=len(variables);lookup[forest,code]=index;vv.append([index,len(pairs)])
            variables.append(dict(index=index,vertex=vid,forest=list(forest),marked_forest_code=json.loads(code),
                                  multiplicity=len(pairs),literal_pairs=pairs))
        assert sum(m for i,m in vv)==78
        vertices.append(dict(index=vid,forest=list(forest),variables=vv))
        if len(vertices)%1000==0:print('COMBINED_VERTEX_BUDGETS',len(vertices),len(variables),flush=True)
    pair_variables=len(variables)
    old_vertex_map=[vertex_lookup[tuple(v['forest'])] for v in old['vertices']]
    old_variable_map=[]
    for var in old['variables']:
        v=old_vertex_map[var['vertex']]
        if var.get('kind')=='maximal_simplex_clique_charge':
            j=len(variables);new=dict(var);new.update(index=j,vertex=v)
            variables.append(new);vertices[v]['variables'].append([j,var['multiplicity']])
        else:
            j=lookup[tuple(var['forest']),B.dump(var['marked_forest_code'])]
            new=variables[j]
            assert new['vertex']==v and new['multiplicity']==var['multiplicity']
            assert new['literal_pairs']==var['literal_pairs']
        old_variable_map.append(j)
    assert len(set(old_variable_map))==len(old['variables'])
    assert len(variables)-pair_variables==22
    for vertex in vertices:assert sum(m for j,m in vertex['variables'] if j<pair_variables)==78
    faces=[]
    for face in old['faces']:
        assert face['index']==len(faces)
        new=dict(face);new['coefficients']=sorted([[old_variable_map[j],m] for j,m in face['coefficients']])
        new['kind']='maximal_simplex_aggregate';new['original_face']=face['index'];faces.append(new)
    # No quad receives a simplex packing charge. Deduplicate only the new quad rows;
    # keep every old row, its demand and every old packing coefficient verbatim.
    quad_lookup={}
    for ctx,ci,required,keys in raw:
        cc=Counter()
        for forest,code,m in keys:cc[lookup[forest,code]]+=m
        coeff=tuple(sorted(cc.items()));assert sum(m for j,m in coeff)==4*required
        assert all(j<pair_variables for j,m in coeff)
        key=required,coeff
        if key not in quad_lookup:
            quad_lookup[key]=len(faces)
            faces.append(dict(index=len(faces),kind='canonical_dotted_quad_aggregate',required=required,
                              coefficients=[list(x) for x in coeff],source_contexts=[]))
        faces[quad_lookup[key]]['source_contexts'].append([ctx,ci])
    assert sum(len(f.get('source_contexts',[])) for f in faces)==81183
    packing=[]
    for col in old['packing_columns']:
        new=dict(col);new['variable']=old_variable_map[col['variable']]
        new['vertex']=old_vertex_map[col['vertex']];packing.append(new)
    @lru_cache(None)
    def finite_forests(rank,first=0):
        if rank==0:return ((),)
        return tuple(tuple(sorted((name,)+tail)) for j,(name,g) in enumerate(models[first:],first)
                     if len(g)<=rank for tail in finite_forests(rank-len(g),j))
    universe=set(finite_forests(13));assert len(universe)==6789 and forests<=universe
    active=sorted({j for face in faces for j,m in face['coefficients']})
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{B.C.TAG}.py',f'scripts/{D.TAG}.py',
           f'results/{OLD}.json',f'scripts/{QUADS}.py',f'results/{QUADS}.json',f'results/{SPECTATORS}.json']
    return dict(sources={p:sha(p) for p in paths},dimension=13,target=12,variables=variables,vertices=vertices,
        faces=faces,active_variables=active,packing_columns=packing,
        old_vertex_to_new=old_vertex_map,old_variable_to_new=old_variable_map,
        absent_vertex_forests=[list(f) for f in sorted(universe-forests)],
        counts=dict(vertices=len(vertices),all_finite_rank13_forests=6789,absent_forests=len(universe-forests),
            inherited_vertices=len(old_forests),new_vertices=len(forests-old_forests),variables=len(variables),
            pair_variables=pair_variables,packing_variables=22,active_variables=len(active),
            zero_only_variables=len(variables)-len(active),old_simplex_rows=10185,
            literal_canonical_flexible_quad_stars=81183,distinct_quad_rows=len(quad_lookup),
            distinct_face_constraints=len(faces)),
        definition='Retain all10185 old flexible simplex aggregate rows and all22 proved simplex packing charges. Reconstruct all81183 canonical flexible dotted full-star quadrilateral aggregates, giving them zero simplex-charge coefficients. Rebuild every78-pair orbit budget for the union of all touched finite rank13 forests, with explicit injective old-index maps. Variables absent from every face row are fixed to zero numerically but retained in the exact domain. No mixed simplex/quad packing columns are included.')


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
    print('COMBINED_SINGLE_LP_START',d['counts'],flush=True)
    result=linprog(objective,A_ub=matrix,b_ub=np.array(rhs,dtype=float),bounds=(0,None),
                   method='highs',options={'time_limit':300})
    solver=dict(status=int(result.status),message=str(result.message),success=bool(result.success),
                objective=None if result.fun is None else float(result.fun),time_limit_seconds=300)
    print('COMBINED_SINGLE_LP_END',solver,flush=True)
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
        # Save one concise authoritative domain before the sole optimization;
        # overwrite this same artifact at completion, without a duplicate export.
        d.update(status='PREPARED_COMBINED_PARTIAL_DOMAIN',solver=None,certificates=[])
        out.write_text(B.dump(d)+'\n')
        solver,certificates=optimize(d)
        feasible=any(c.get('reaches_target12',False) for c in certificates)
        excluded=any(c.get('excludes_target12',False) for c in certificates)
        assert not(feasible and excluded)
        d.update(status='PASS_EXACT_COMBINED_PARTIAL_RULE_AT12' if feasible else
                 'PASS_EXACT_COMBINED_PARTIAL_MODEL_OBSTRUCTION' if excluded else 'UNRESOLVED_COMBINED_PARTIAL_MODEL',
            solver=solver,certificates=certificates,
            scope='One300-second optimization only. The inherited flexible simplex and canonical flexible dotted-quad domains are complete as saved inputs, but fixed-determinant stars, other connected cores, all-finite quad stars and general remaining faces are not covered. Only touched finite forests have budget variables. A feasible rule is not a dimension theorem; an exact dual obstructs this particular nonnegative pair-orbit plus22-simplex-charge family on these local configurations. Mixed aggregate packing charges are intentionally absent.')
        out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound'))) for c in d['certificates']],flush=True)


if __name__=='__main__':main()
