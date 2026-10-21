#!/usr/bin/env python3
"""Prepare, without optimization, the positively forced missing base stars.

The export is an additive patch to a frozen model, not a duplicated model.
It includes maximal-edge incidence rows with zero curvature demand.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_prism_simplex_base_incidence_20261017 as I
import fresh_sqrt5_prism_shared_transfers_20261017 as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_forced_base_simplex_domain_20261017'


def run():
    incidence=P.read(I.TAG);shared=P.read(S.TAG);old=P.read(I.LP)
    alpha={i:F(a) for i,a in shared['repaired_dual']['face_multipliers']}
    nfaces=len(old['faces']);nvars=len(old['variables']);nvertices=len(old['vertices'])
    assert (nfaces,nvars,nvertices)==(78335,210115,6558)
    byface=defaultdict(list)
    for record in incidence['contexts']:byface[record[1]].append(record)
    positive=[r for r in incidence['contexts'] if alpha.get(r[1],F(0))>0]
    selected=sorted({t for r in positive for t,m in r[4]
                     if not incidence['base_types'][t]['represented_in_lp']})
    assert len(positive)==1858 and len(selected)==4195
    selected_set=set(selected);origins=defaultdict(list)
    for ctx,fi,ci,k,entries in positive:
        for t,m in entries:
            if t in selected_set:origins[t].append([ctx,fi,ci,k,m])
    assert set(origins)==selected_set
    # The complete minimum is over EVERY provenance, even when a particular
    # base only appears in one of them. Zero minima are explicitly retained
    # in the provenance/source maps and do not become false forced loads.
    transfers=defaultdict(list);varying=0
    for fi,records in sorted(byface.items()):
        face=old['faces'][fi]
        assert sorted((r[0],r[2]) for r in records)==sorted(map(tuple,face['source_contexts']))
        ks={r[3] for r in records};assert len(ks)==1;k=next(iter(ks))
        alternatives=[dict(r[4]) for r in records]
        varying+=any(x!=alternatives[0] for x in alternatives[1:])
        keys=set().union(*(set(x) for x in alternatives)) & selected_set
        for t in sorted(keys):
            m=min(x.get(t,0) for x in alternatives)
            if m:transfers[t].append([fi,k*m])
    models={name:g for name,g in B.C.models() if len(g)<=13}
    ranks={name:len(g) for name,g in models.items()}
    finite={B.C.code(g):name for name,g in models.items()}
    def restrict(g,ids):
        ids=set(ids)
        return {i:{j:m for j,m in g[i].items() if j in ids} for i in sorted(ids)}
    def finite_name(g):
        assert sum(map(len,g.values()))==2*len(g)-2
        assert all(m!=0 for row in g.values() for m in row.values())
        return finite[B.C.code(g)]
    def forest_code(pieces,marked=()):
        return B.dump(sorted(B.C.code(p,marked) for p in pieces))
    def forest_graph(forest):
        g={}
        for name in forest:
            offset=len(g)
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13
        return g
    @lru_cache(None)
    def core_corners(ci):
        c=incidence['cores'][ci];g={i:{} for i in range(c['nodes'])}
        for i,j,(m,a,b) in c['edges']:g[i][j]=g[j][i]=m
        spherical=[];corners=[]
        for p in g:
            h=restrict(g,set(g)-{p});pieces=B.C.pieces(h)
            try:forest=tuple(sorted(finite_name(piece) for piece in pieces))
            except (AssertionError,KeyError):continue
            spherical.append(p)
            corners.append((p,h,pieces,forest))
        assert spherical==c['spherical_deletions']
        assert len(spherical)==c['simplex_dimension']+1
        return tuple(corners)
    raw=[];touched=set();dimensions=Counter();core_usage=Counter()
    for t in selected:
        typ=incidence['base_types'][t];ci=typ['core'];sp=typ['spectator_forest']
        c=incidence['cores'][ci];dimension=typ['simplex_dimension'];assert dimension==c['simplex_dimension']
        assert c['nodes']+sum(ranks[name] for name in sp)==14
        assert not(c['nodes']%2==0 and all(ranks[name]%2==0 for name in sp))
        Sset=set(c['spherical_deletions']);rows=[]
        for omitted,h0,pieces0,forest0 in core_corners(ci):
            h={i:dict(row) for i,row in h0.items()};offset=c['nodes']
            for name in sp:
                for i,row in models[name].items():h[offset+i]={offset+j:m for j,m in row.items()}
                offset+=ranks[name]
            assert offset==14 and len(h)==13
            pieces=B.C.pieces(h)
            forest=tuple(sorted(finite_name(piece) for piece in pieces))
            assert forest==tuple(sorted(list(forest0)+sp));touched.add(forest)
            block=sorted(Sset-{omitted});assert len(block)==dimension and set(block)<=set(h)
            pair=Counter(forest_code(pieces,p) for p in combinations(block,2))
            port=Counter(forest_code(pieces,(p,)) for p in block)
            rows.append(dict(omitted=omitted,forest=forest,literal_port_block=block,
                             pair_codes=pair,port_codes=port))
        assert len(rows)==dimension+1
        raw.append((t,rows));dimensions[dimension]+=1;core_usage[ci]+=1
        if len(raw)%1000==0:print('FORCED_BASE_CORNERS',len(raw),len(touched),flush=True)
    old_forest={tuple(v['forest']):v['index'] for v in old['vertices']}
    pair_lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index']
                 for v in old['variables'][:old['counts']['pair_variables']]}
    assert len(pair_lookup)==203565
    port_lookup={}
    for v in old['variables']:
        if v.get('kind')=='maximal_simplex_disjoint_port_charge':
            key=tuple(v['forest']),B.dump(v['marked_forest_code']);assert key not in port_lookup
            port_lookup[key]=v['index']
    assert len(port_lookup)==6528
    additions=[];new_vertices=[];budget_additions=defaultdict(list);touched_budgets=[]
    vertex_lookup=dict(old_forest);new_pair_indices=[];new_port_indices=[];used_port_keys=set()
    def append_variable(v,forest,code,cost,kind,**extra):
        j=nvars+len(additions)
        rec=dict(index=j,vertex=v,forest=list(forest),marked_forest_code=json.loads(code),
                 multiplicity=cost,kind=kind,**extra)
        additions.append(rec);budget_additions[v].append([j,cost]);return j
    for forest in sorted(touched):
        if forest not in vertex_lookup:
            v=nvertices+len(new_vertices);vertex_lookup[forest]=v
            new_vertices.append(dict(index=v,forest=list(forest)))
        v=vertex_lookup[forest];g=forest_graph(forest);pieces=B.C.pieces(g)
        pairs=defaultdict(list);ports=defaultdict(list)
        for pair in combinations(range(13),2):pairs[forest_code(pieces,pair)].append(list(pair))
        for p in range(13):ports[forest_code(pieces,(p,))].append(p)
        assert sum(map(len,pairs.values()))==78 and sum(map(len,ports.values()))==13
        pv=[];uv=[]
        for code,literal in sorted(pairs.items()):
            key=forest,code
            if key not in pair_lookup:
                assert forest not in old_forest
                j=append_variable(v,forest,code,len(literal),'ordinary_pair_charge',literal_pairs=literal)
                pair_lookup[key]=j;new_pair_indices.append(j)
            j=pair_lookup[key]
            saved=old['variables'][j] if j<nvars else additions[j-nvars]
            assert saved['vertex']==v and saved['literal_pairs']==literal and saved['multiplicity']==len(literal)
            pv.append([j,len(literal)])
        for code,literal in sorted(ports.items()):
            key=forest,code
            if key not in port_lookup:
                j=append_variable(v,forest,code,len(literal),'maximal_simplex_disjoint_port_charge',
                    literal_port_orbit=literal,normalization='one unit per used port; maximal edges included')
                port_lookup[key]=j;new_port_indices.append(j)
            j=port_lookup[key]
            saved=old['variables'][j] if j<nvars else additions[j-nvars]
            assert saved['vertex']==v and saved['literal_port_orbit']==literal and saved['multiplicity']==len(literal)
            uv.append([j,len(literal)])
        touched_budgets.append(dict(vertex=v,forest=list(forest),pair_variables=pv,port_variables=uv))
    for v in new_vertices:
        v['variables']=budget_additions[v['index']]
        assert sum(m for j,m in v['variables'] if j in set(new_pair_indices))==78
    faces=[];base_to_face={};columns=[];robust_mass={};literal_corners=0
    for t,rows in raw:
        typ=incidence['base_types'][t];dimension=typ['simplex_dimension'];fi=nfaces+len(faces)
        base_to_face[t]=fi;pair_coeff=Counter();port_coeff=Counter();corners=[]
        for row in rows:
            forest=row['forest'];v=vertex_lookup[forest]
            pc=sorted((pair_lookup[forest,code],m) for code,m in row['pair_codes'].items())
            uc=sorted((port_lookup[forest,code],m) for code,m in row['port_codes'].items())
            for j,m in pc:pair_coeff[j]+=m
            for j,m in uc:port_coeff[j]+=m
            used_port_keys.update((forest,code) for code in row['port_codes'])
            corners.append(dict(omitted=row['omitted'],vertex=v,forest=list(forest),
                literal_port_block=row['literal_port_block'],pair_coefficients=pc,port_coefficients=uc))
        assert sum(pair_coeff.values())==3*comb(dimension+1,3)
        assert sum(port_coeff.values())==dimension*(dimension+1)
        required=2*comb(dimension+1,3)
        if dimension==1:assert required==0 and not pair_coeff and len(corners)==2
        coeff=pair_coeff+port_coeff;assert len(coeff)==len(pair_coeff)+len(port_coeff)
        load=sum((m*alpha.get(q,F(0)) for q,m in transfers[t]),F(0))
        capacity=(dimension+1)*(13-dimension);mass=load/capacity;robust_mass[t]=mass
        faces.append(dict(index=fi,kind='maximal_edge_incidence' if dimension==1 else 'maximal_simplex_aggregate',
            required=required,base_type=t,core=typ['core'],spectator_forest=typ['spectator_forest'],
            simplex_dimension=dimension,coefficients=[list(x) for x in sorted(coeff.items())],
            pair_coefficients=[list(x) for x in sorted(pair_coeff.items())],
            port_coefficients=[list(x) for x in sorted(port_coeff.items())],corners=corners,
            positive_quad_origins=origins[t],shared_dual_forced_minimum_mass=str(mass)))
        literal_corners+=len(corners)
        if transfers[t]:
            columns.append(dict(index=len(columns),base_type=t,simplex_face=fi,simplex_dimension=dimension,
                vertex_cost=0,shared_capacity_per_simplex=capacity,
                positive_quad_coefficients=transfers[t],negative_simplex_coefficient=[fi,-capacity],
                current_shared_dual_load=str(load),current_forced_simplex_mass=str(mass)))
    # The original disjoint columns exhaust their nonzero old-domain entries.
    # Newly activated ports therefore have zero on every old simplex row;
    # quads have zero on all port columns. Old clique columns can safely be
    # assigned zero on the new rows (a restriction of their flag domain).
    old_used_ports={j for f in old['faces'][:10185] for j,m in f['coefficients'] if j>=203587}
    assert old_used_ports=={v['index'] for v in old['variables'] if v.get('kind')=='maximal_simplex_disjoint_port_charge'}
    for t in selected:
        assert not incidence['base_types'][t]['represented_in_lp']
    assert len(faces)==len(base_to_face)==4195
    active_additions={j for f in faces for j,m in f['coefficients'] if j>=nvars}
    new_forests={tuple(v['forest']) for v in new_vertices}
    assert new_forests<=set(map(tuple,old['absent_vertex_forests']))
    # Exact comparison diagnostics need no new solve. These are requirements
    # on the CURRENT quad masses, not a repaired dual or optimized bound.
    mass_by_dimension=defaultdict(F)
    for t,mass in robust_mass.items():mass_by_dimension[incidence['base_types'][t]['simplex_dimension']]+=mass
    core_profiles=[]
    for ci in sorted(core_usage):
        c=incidence['cores'][ci]
        core_profiles.append(dict(core=ci,nodes=c['nodes'],simplex_dimension=c['simplex_dimension'],
            stars=core_usage[ci],finite_labels=sorted({m for i,j,(m,a,b) in c['edges'] if m}),
            dotted_squared_parameters=sorted({(a,b) for i,j,(m,a,b) in c['edges'] if not m}),
            inherited_Qsqrt5_admissible_gauge=True))
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{B.C.TAG}.py',f'scripts/{P.TAG}.py',
        f'scripts/{I.TAG}.py',f'results/{I.TAG}.json',f'results/{I.LP}.json',
        f'scripts/{S.TAG}.py',f'results/{S.TAG}.json',
        'results/fresh_simplex_disjoint_ports_20261017.md',
        'results/fresh_prism_simplex_shared_capacity_20261017.md']
    return dict(status='PASS_EXACT_FORCED_MISSING_BASE_DOMAIN_WITH_MAXIMAL_EDGES',
        sources={p:P.sha(p) for p in paths},dimension=13,target=12,
        inherited_model=f'results/{I.LP}.json',inherited_shared_transfers=f'results/{S.TAG}.json',
        inherited_counts=dict(faces=nfaces,variables=nvars,vertices=nvertices),
        counts=dict(positive_quad_provenances=len(positive),selected_missing_stars=len(selected),
            maximal_edge_rows=dimensions[1],triangle_rows=dimensions[2],tetrahedron_rows=dimensions[3],
            four_simplex_rows=dimensions[4],connected_base_cores=len(core_usage),literal_corners=literal_corners,
            robustly_positive_missing_stars=sum(m>0 for m in robust_mass.values()),
            provenance_only_stars=sum(m==0 for m in robust_mass.values()),touched_forests=len(touched),
            new_forests=len(new_vertices),inherited_touched_forests=len(touched)-len(new_vertices),
            new_pair_variables=len(new_pair_indices),new_port_variables=len(new_port_indices),
            new_variables=len(additions),active_new_variables=len(active_additions),
            zero_only_new_variables=len(additions)-len(active_additions),new_shared_transfer_columns=len(columns),
            positive_quad_transfer_entries=sum(len(c['positive_quad_coefficients']) for c in columns),
            all_quad_provenances_checked=len(incidence['contexts']),all_quad_rows_checked=len(byface),
            varying_quad_provenance_rows=varying,combined_forests=nvertices+len(new_vertices),
            remaining_absent_forests=6789-nvertices-len(new_vertices)),
        selected_base_types=selected,base_to_new_face=[[t,fi] for t,fi in sorted(base_to_face.items())],
        core_profiles=core_profiles,new_variables=additions,new_pair_variable_indices=new_pair_indices,
        new_port_variable_indices=new_port_indices,new_vertices=new_vertices,
        existing_vertex_budget_additions=[[v,entries] for v,entries in sorted(budget_additions.items()) if v<nvertices],
        touched_vertex_budgets=touched_budgets,new_faces=faces,new_shared_transfer_columns=columns,
        forced_mass_totals_by_dimension={str(k):str(v) for k,v in sorted(mass_by_dimension.items())},
        first_domain='All missing base stars appearing in ANY provenance of a positive quad row in the frozen shared-repaired dual. This includes provenance-only stars whose coefficientwise minimum is zero; no chosen representative replaces the full provenance set.',
        maximal_edge_semantics='The triangular-port equivalence classes at a simple finite vertex include isolated singleton ports. These are maximal edges, and every port belongs to exactly one maximal simplex, including dimension1. A maximal-edge row has demand0, no pair-corner entries, and two endpoint one-port charges. It is an incidence bookkeeping row, not positive curvature demand. Its shared prism capacity is(1+1)*(13-1)=24; squares use all four edges.',
        patch_semantics='Keep every inherited variable,row,budget and all represented shared transfers unchanged. Append new variable indices, new forest budgets and listed additions to existing budgets, then append the new rows and zero-vertex-cost shared transfers. Reused one-port variables receive the new row coefficients, preventing separate double budgets. New one-port variables have zero old-domain coefficients; all quads have zero one-port coefficients. The22 old clique charges are deliberately zero on new rows. Ordinary pair variables no longer form one initial index interval after patching; use the explicit new_pair_variable_indices with the inherited pair interval.',
        arithmetic_scope='Every base is obtained by deleting one normal from a frozen exact locally admissible Q(sqrt5) dotted prism star. All complete finite one-deletions and corner forests are reconstructed here; the ambient admissible gauge is inherited from the parent star. No new ground field is introduced. Every selected full rank14 star is determinant-flexible by the same even-component criterion. This does not classify all possible arithmetic forms or assert a common global chamber.',
        coverage='This is a prepared partial-domain patch, not a solve or a dimension bound. It covers these4195 forced missing stars only; remaining missing base stars, fixed-discriminant classes, other connected simplex cores, all-finite quadrilateral aggregates and general unenumerated faces remain outside. The newly activated simplex rows require any future mixed-port or clique coefficients to be reconstructed rather than inherited from a prior dual diagnostic.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],flush=True)
