#!/usr/bin/env python3
"""Lift every positive merged quad mass to actual saved provenances.

The deterministic uniform split is exact. No optimizer is invoked. Shared
base loads use literal provenance coefficients, never coordinatewise minima.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_forced_base_dual_repair_20261017 as R
import fresh_sqrt5_mixed_forced_base_transfers_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_actual_provenance_dual_20261017'


def code(pieces,marks=()):return tuple(sorted(B.C.code(p,marks) for p in pieces))


def freeze(x):return tuple(freeze(v) for v in x) if isinstance(x,(list,tuple)) else x


def run():
    saved=P.read(M.TAG);incidence=P.read(R.N.I.TAG)
    d,patch,shared=R.domain()
    old_alpha={i:F(a) for i,a in saved['face_multipliers']}
    old_mu={i:F(a) for i,a in saved['augmented_vertex_multipliers']}
    byface=defaultdict(list)
    for row in incidence['contexts']:
        if old_alpha.get(row[1],F(0)):byface[row[1]].append(row)
    assert len(byface)==1582 and sum(map(len,byface.values()))==1858
    assert set(byface)=={i for i in old_alpha if d['faces'][i]['kind']=='canonical_dotted_quad_aggregate'}
    typemap={t['index']:t['old_face'] for t in incidence['base_types'] if t['represented_in_lp']}
    for t,fi in patch['base_to_new_face']:
        assert t not in typemap;typemap[t]=fi
    needed=sorted({t for rows in byface.values() for row in rows for t,m in row[4]})
    assert set(needed)<=set(typemap)
    assert len({typemap[t] for t in needed})==len(needed)
    loads=defaultdict(F);origins=defaultdict(list);provenances=[];row_constraints=[];inconsistent=[]
    for fi,rows in sorted(byface.items()):
        face=d['faces'][fi]
        assert sorted((r[0],r[2]) for r in rows)==sorted(map(tuple,face['source_contexts']))
        ks={r[3] for r in rows};assert len(ks)==1;k=next(iter(ks))
        total=4 if k==2 else 2
        assert all(sum(m for t,m in r[4])==total for r in rows)
        vectors=[dict(r[4]) for r in rows];keys=sorted(set().union(*(set(v) for v in vectors)))
        minimum={t:min(v.get(t,0) for v in vectors) for t in keys}
        minimum={t:m for t,m in minimum.items() if m}
        varying=any(v!=vectors[0] for v in vectors[1:])
        assert (sum(minimum.values())<total)==varying
        ids=[];a=old_alpha[fi]/len(rows)
        for ctx,check,ci,kk,entries in rows:
            assert check==fi and kk==k
            j=len(provenances);ids.append(j)
            coeff=[[t,k*m] for t,m in entries]
            provenances.append(dict(index=j,source_context=ctx,merged_face=fi,canonical_prism_core=ci,
                prism_dimension=k,mass=str(a),base_coefficients=coeff,
                positive_cost_and_mixed_feature_source_face=fi))
            for t,m in coeff:loads[t]+=a*m;origins[t].append([j,m])
        assert sum((F(provenances[j]['mass']) for j in ids),F(0))==old_alpha[fi]
        row_constraints.append(dict(face=fi,mass=str(old_alpha[fi]),provenance_indices=ids,
            prism_dimension=k,constant_base_total=total,coordinatewise_minimum_base_total=sum(minimum.values()),
            minimum_coefficients=[[t,k*m] for t,m in sorted(minimum.items())]))
        if varying:
            inconsistent.append(dict(face=fi,provenance_indices=ids,prism_dimension=k,
                fixed_weighted_total=k*total,minimum_weighted_total=k*sum(minimum.values()),
                missing_weighted_total=k*(total-sum(minimum.values()))))
    assert len(inconsistent)==243
    alpha=dict(old_alpha);constraints=[];increases=[]
    for t in needed:
        typ=incidence['base_types'][t];fi=typemap[t];dim=typ['simplex_dimension'];cap=(dim+1)*(13-dim)
        assert d['faces'][fi]['required']==2*comb(dim+1,3)
        need=loads[t]/cap;old=alpha.get(fi,F(0));alpha[fi]=max(old,need)
        if alpha[fi]>old:increases.append([fi,t,str(old),str(alpha[fi]),str(alpha[fi]-old)])
        constraints.append(dict(base_type=t,simplex_face=fi,simplex_dimension=dim,capacity=cap,
            positive_provenance_coefficients=origins[t],actual_positive_load=str(loads[t]),
            old_simplex_mass=str(old),repaired_simplex_mass=str(alpha[fi]),
            signed_load=str(loads[t]-cap*alpha[fi])))
    alpha={i:a for i,a in alpha.items() if a}
    assert all(alpha[i]==a for i,a in old_alpha.items() if i in byface)
    assert all(alpha.get(i,F(0))>=a for i,a in old_alpha.items())
    print('ACTUAL_SPLIT',len(provenances),'BASES',len(needed),'RAISED',len(increases),flush=True)

    # Old pair/clique/one-port vectors do not depend on quad provenance by
    # the defining exact row merger. The frozen mixed certificate separately
    # compared every one of these1858 actual provenance footprints.
    _,positive_mu,tight,num,_=R.D.sparse_dual_bound(d,alpha)
    assert all(positive_mu.get(v,F(0))>=a for v,a in old_mu.items())
    forests=set(positive_mu)
    models=dict(B.C.models());finite={B.C.code(g):name for name,g in models.items()}
    vlookup={tuple(v['forest']):v['index'] for v in d['vertices']}
    pairlookup={(v['vertex'],B.dump(v['marked_forest_code'])):v['index'] for v in d['variables'][:203565]}
    portlookup={(v['vertex'],B.dump(v['marked_forest_code'])):v['index'] for v in d['variables']
                if v.get('kind')=='maximal_simplex_disjoint_port_charge'}
    columns=[];orbit={}
    for v in sorted(forests):
        g={}
        for name in d['vertices'][v]['forest']:
            offset=len(g)
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13;pieces=B.C.pieces(g);groups=Counter(code(pieces,U) for U in combinations(range(13),3))
        assert sum(groups.values())==286
        for key,cost in sorted(groups.items()):
            j=len(columns);orbit[v,key]=j
            columns.append(dict(index=j,vertex=v,marked_forest_code=key,vertex_cost=cost))
    saved_columns=saved['triple_orbit_columns'];features={};inherited_provenances=0
    for rec in saved['face_features']:
        fi=rec['face']
        if fi not in alpha:continue
        assert rec['provenance_dependent'] is False
        coeff=[]
        for j,m in rec['coefficients']:
            oldc=saved_columns[j];key=freeze(oldc['marked_forest_code']);v=oldc['vertex']
            dest=orbit[v,key];assert columns[dest]['vertex_cost']==oldc['vertex_cost']
            coeff.append([dest,m])
        if fi in byface:
            assert rec['provenances']==len(byface[fi]);inherited_provenances+=rec['provenances']
        features[fi]=dict(face=fi,source='frozen_complete_provenance_features',coefficients=coeff,
            provenance_count=rec['provenances'])
    assert inherited_provenances==1858
    missing=sorted(set(alpha)-set(features));reverse={fi:t for t,fi in typemap.items()}
    corner_count=0
    for fi in missing:
        assert fi not in byface;t=reverse[fi];typ=incidence['base_types'][t];core=incidence['cores'][typ['core']]
        g={i:{} for i in range(core['nodes'])}
        for i,j,(m,a,b) in core['edges']:g[i][j]=g[j][i]=m
        offset=core['nodes']
        for name in typ['spectator_forest']:
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
            offset+=len(models[name])
        assert offset==len(g)==14
        S=set(core['spherical_deletions']);dim=typ['simplex_dimension'];assert len(S)==dim+1
        cc=Counter();pp=Counter();uu=Counter();literal=[]
        for omitted in sorted(S):
            h={i:{j:m for j,m in row.items() if j!=omitted} for i,row in g.items() if i!=omitted}
            assert len(h)==13 and all(m>0 for row in h.values() for m in row.values())
            pieces=B.C.pieces(h)
            assert all(sum(map(len,p.values()))==2*len(p)-2 for p in pieces)
            forest=tuple(sorted(finite[B.C.code(p)] for p in pieces));v=vlookup[forest]
            block=sorted(S-{omitted});footprint=list(combinations(block,2))
            cp=Counter(pairlookup[v,B.dump(code(pieces,p))] for p in footprint)
            cu=Counter(portlookup[v,B.dump(code(pieces,(p,)))] for p in block)
            cm=Counter(orbit[v,code(pieces,U)] for U in M.M.triple_indicators(footprint))
            assert sum(cm.values())==comb(dim,3)
            pp.update(cp);uu.update(cu);cc.update(cm)
            literal.append(dict(omitted=omitted,vertex=v,forest=forest,port_block=block,
                                mixed_coefficients=[list(x) for x in sorted(cm.items())]))
        assert dict(pp)=={j:m for j,m in d['faces'][fi]['coefficients'] if j<203565}
        assert dict(uu)=={j:m for j,m in d['faces'][fi]['coefficients'] if d['variables'][j].get('kind')=='maximal_simplex_disjoint_port_charge'}
        assert sum(cc.values())==(dim+1)*comb(dim,3)
        features[fi]=dict(face=fi,source='reconstructed_newly_positive_represented_base',base_type=t,
            simplex_dimension=dim,coefficients=[list(x) for x in sorted(cc.items())],literal_corners=literal,
            provenance_count=1)
        corner_count+=len(S)
    assert set(features)==set(alpha)
    mixed_loads=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in features[fi]['coefficients']:mixed_loads[j]+=m*a
    mu=dict(positive_mu);mixed_increases=[]
    for col in columns:
        j,v,cost=col['index'],col['vertex'],col['vertex_cost'];need=mixed_loads[j]/cost
        if need>mu.get(v,F(0)):
            mixed_increases.append([j,v,str(mu.get(v,F(0))),str(need)]);mu[v]=need
    # Fresh full arithmetic checks: every positive-cost column, every mixed
    # triple orbit, all inherited conservative signed columns, and all actual
    # provenance-shared columns. Other zero-mass provenances have zero load.
    column_loads=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in d['faces'][fi]['coefficients']:column_loads[j]+=m*a
    for var in d['variables']:assert column_loads[var['index']]<=var['multiplicity']*mu.get(var['vertex'],F(0))
    for col in columns:assert mixed_loads[col['index']]<=col['vertex_cost']*mu[col['vertex']]
    conservative=R.signed_loads(d,alpha);assert all(v<=0 for v in conservative)
    for rec in constraints:
        actual=sum((m*F(provenances[j]['mass']) for j,m in rec['positive_provenance_coefficients']),F(0))
        assert actual==F(rec['actual_positive_load'])
        assert actual-rec['capacity']*alpha[rec['simplex_face']]==F(rec['signed_load'])<=0
    assert num==sum((d['faces'][fi]['required']*a for fi,a in alpha.items()),F(0))
    den=sum(mu.values(),F(0));bound=num/den
    weighted_deficit=sum((old_alpha[r['face']]*r['missing_weighted_total'] for r in inconsistent),F(0))
    paths=[f'scripts/{TAG}.py',f'scripts/{R.TAG}.py',f'scripts/{R.D.TAG}.py',f'scripts/{B.C.TAG}.py',
        f'scripts/{M.TAG}.py',f'results/{M.TAG}.json',f'results/{R.N.I.TAG}.json',
        f'results/{R.N.TAG}.json',f'results/{R.N.I.LP}.json',f'results/{R.S.TAG}.json']
    return dict(status='PASS_EXACT_ACTUAL_PROVENANCE_DUAL_BARRIER' if bound>12 else 'PASS_EXACT_UNIFORM_PROVENANCE_REPAIR_AT_MOST12',
        sources={p:P.sha(p) for p in paths},dimension=13,target=12,bound=str(bound),excludes_target12=bound>12,
        old_relaxed_bound=saved['bound'],objective_numerator=str(num),normalizer=str(den),
        counts=dict(positive_merged_quad_rows=len(byface),actual_positive_provenances=len(provenances),
            incompatible_minimum_vectors=len(inconsistent),actual_shared_base_constraints=len(constraints),
            inherited_conservative_transfer_columns=len(conservative),raised_simplex_rows=len(increases),
            newly_positive_rows=len(set(alpha)-set(old_alpha)),positive_rows=len(alpha),positive_forests=len(mu),
            positive_cost_columns=len(d['variables']),mixed_triple_orbit_columns=len(columns),
            newly_reconstructed_represented_simplex_features=len(missing),newly_reconstructed_corners=corner_count,
            mixed_vertex_increases=len({r[1] for r in mixed_increases})),
        provenance_masses=provenances,merged_row_mass_constraints=row_constraints,
        unattainable_minimum_vectors=inconsistent,weighted_minimum_incidence_deficit=str(weighted_deficit),
        actual_shared_base_constraints=constraints,simplex_mass_increases=increases,
        face_multipliers=[[fi,str(a)] for fi,a in sorted(alpha.items())],
        vertex_multipliers=[[v,str(a)] for v,a in sorted(mu.items())],
        positive_cost_tight_columns={str(v):j for v,j in tight.items()},
        triple_orbit_columns=columns,positive_face_mixed_features=[features[fi] for fi in sorted(features)],
        mixed_column_loads=[[j,str(a)] for j,a in sorted(mixed_loads.items()) if a],mixed_mass_increases=mixed_increases,
        inherited_signed_loads=list(map(str,conservative)),
        split_definition='Each merged quad mass is divided equally among ALL its saved provenances. The original merged masses are preserved exactly. Actual k-times-base multiplicities give each shared load. Retain every old simplex mass and raise it only as required by the actual load/capacity; then recount every vertex and mixed column.',
        minimum_diagnosis='For fixed prism dimension every actual base vector has the same sum (four for squares, two otherwise). A varying row has a coordinatewise-minimum vector with strictly smaller sum, hence that minimum is outside the convex hull of its actual provenances. This proves the old relaxation is strict, not that its dual cannot be repaired.',
        scope='Exact dual test on the previously saved positive support, split by actual local full-star provenance; all pair/clique/one-port/mixed and represented/forced shared-transfer constraints are checked. Uniform splitting is one feasible allocation, not an optimum. The new distribution is not asserted globally realizable in one compact polytope or arithmetic form. Missing face/base classes remain outside. A bound above12 obstructs only this partial weighting family; no dimension bound follows.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['bound'],flush=True)
