#!/usr/bin/env python3
"""Saved-only support/provenance migration to the new5121 dual masses.

No matrix reconstruction, geometry, clique enumeration or optimization.
Old local charges are reevaluated; old optimal masses are never transported.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cross_vertex_support_migration_20261018'
MODEL='fresh_cross_vertex_current_full_model_20261018'
OLD='fresh_weighted_raw110_full_model_20261017'
H='fresh_current_historical_weighted_corners_20261017'
HD='fresh_historical22_corner_domain_20261017'
HC='fresh_historical22_corner_cones_20261017'
PATCH='fresh_unowned_cover_pair_geometry_20261017'
EIGHT='fresh_historical22_eight_pair_geometry_20261017'
RAW88='fresh_historical22_raw88_raw25_pair_geometry_20261017'
RAW125='fresh_raw125_raw205_actual_corner_frontier_20261017'
SHORT='fresh_290_class_shortlist_pair_geometry_20261017'
P2599='fresh_2599_raw110_raw24_pair_geometry_20261017'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
FORCED='fresh_sqrt5_forced_base_simplex_domain_20261017'
ACTUAL='fresh_sqrt5_actual_provenance_dual_20261017'
QUAD='fresh_sqrt5_rank13_packing_quad_separation_20261017'
LOCAL='fresh_mixed_raw39_full_model_20261017_local_cones'
EXTRA=('fresh_mixed_raw39_actual_corners_20261017','fresh_mixed_raw23_raw25_corners_20261017')
DOMAIN_STEMS=(H,PATCH,HD,HC,EIGHT,RAW88,RAW125,SHORT,P2599)


@lru_cache(None)
def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())


@lru_cache(None)
def digest(path):
    h=sha256()
    with (ROOT/path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def dual(stem):return next(c for c in read(stem)['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')


def current_domains():
    """Latest complete cone rows, with the matching immutable orbit domain."""
    rows={}
    def put(tag,row,groups=None):
        v=row['vertex'];previous=rows.get(v)
        if groups is None:groups=row.get('corner_orbits') or (row.get('union') or {}).get('corner_orbits')
        if groups is None:
            assert previous is not None;groups=previous['groups'];group_source=previous['group_source']
        else:group_source=tag
        assert [g['index'] for g in groups]==list(range(len(groups)))
        assert all(len(a)==len(groups) for a in row['maximal_count_vectors'])
        assert len(row['maximal_count_vectors'])==len(row['literal_representatives'])
        assert len(row['certificates'][0]['weights'])==len(groups)
        if previous and len(previous['groups'])==len(groups):
            # Patches at fixed domains retain exact coefficient coordinates.
            assert [g['positive_face_coefficients'] for g in previous['groups']]==[g['positive_face_coefficients'] for g in groups]
        rows[v]=dict(vertex=v,source=tag,group_source=group_source,groups=groups,row=row)
    for row in read(H)['forest_records']:put(H,row)
    for row in read(PATCH)['forest_records']:put(PATCH,row)
    domain={r['vertex']:r for r in read(HD)['forest_records']}
    for row in read(HC)['forest_records']:put(HC,row,domain[row['vertex']]['corner_orbits']);rows[row['vertex']]['group_source']=HD
    put(EIGHT,read(EIGHT))
    for tag in (RAW88,RAW125,SHORT,P2599):
        for row in read(tag)['forest_records']:put(tag,row)
    assert len(rows)==571
    return rows


def exact_charge_load(vectors,demands,weights):
    w=list(map(F,weights));assert len(w)==len(demands) and all(a>=0 for a in w)
    assert all(sum((F(n)*a for n,a in zip(vector,w)),F())<=1 for vector in vectors)
    return sum((a*b for a,b in zip(w,demands)),F())


def bindings(tags):
    checked={};out={f'scripts/{TAG}.py':digest(f'scripts/{TAG}.py')}
    for tag in tags:
        d=read(tag);assert d['status'].startswith('PASS'),tag
        for p,h in d.get('sources',d.get('prepared',{}).get('sources',{})).items():
            if p not in checked:checked[p]=digest(p)
            assert checked[p]==h,(tag,p)
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{tag}.{ext}'
            if (ROOT/p).exists():out[p]=digest(p)
    return out


def run():
    model=read(MODEL);cert=dual(MODEL);old=dual(OLD)
    alpha={i:F(a) for i,a in cert['face_multipliers']};oa={i:F(a) for i,a in old['face_multipliers']}
    mu={v:F(a) for v,a in cert['vertex_multipliers']};assert len(alpha)==5121 and len(oa)==5053
    assert len(mu)==1050 and all(a>0 for a in alpha.values())
    assert model['prepared']['provenance_to_face']==read(OLD)['prepared']['provenance_to_face']
    pmap={f:j for j,f in model['prepared']['provenance_to_face']};assert len(pmap)==len(model['prepared']['provenance_to_face'])
    faces=list(read(FORESTS)['faces'])
    assert len(faces)==read(FORCED)['inherited_counts']['faces']==78335
    for row in read(FORCED)['new_faces']:
        assert row['index']==len(faces);faces.append(row)
    assert len(faces)==82530
    pr=read(ACTUAL)['provenance_masses'];cores=read(QUAD)['canonical_cores']
    selected={r['face']:r for r in read(H)['selected_faces']}
    for z in read(HD)['added_faces']:selected[z['face']]=z['selected_face']
    for z in read(RAW125)['selected_faces']:selected[z['face']]=z
    assert len(selected)==290
    historical_order=read(H)['source_precedence']+list(EXTRA)
    historical={}
    for tag in historical_order:
        for row in read(tag)['selected_faces']:historical.setdefault(row['face'],[]).append(dict(source=tag,record=row))
    qrows=[];srows=[]
    for fi,mass in sorted(alpha.items()):
        if fi in pmap:
            pi=pmap[fi];r=pr[pi];ci=r['canonical_prism_core'];ctx=r['source_context'];core=cores[ci]
            q=dict(face=fi,provenance=pi,provenance_kind='explicit',source_context=ctx,
                   canonical_core=ci,raw_family=core['source_core'],required=core['required'])
        elif fi<len(faces) and faces[fi]['kind']=='canonical_dotted_quad_aggregate':
            ff=faces[fi];contexts=ff['source_contexts'];q=dict(face=fi,required=ff['required'],source_contexts=contexts)
            if len(contexts)==1:
                ctx,ci=contexts[0];q.update(provenance=None,provenance_kind='unique_saved_context',source_context=ctx,
                                           canonical_core=ci,raw_family=cores[ci]['source_core'])
            else:q.update(provenance_kind='merged_unexpanded',raw_families=sorted({cores[ci]['source_core'] for ctx,ci in contexts}))
        else:
            assert fi<len(faces) and ('simplex' in faces[fi]['kind'] or 'edge' in faces[fi]['kind']),(fi,len(faces),faces[fi]['kind'] if fi<len(faces) else None)
            srows.append(dict(face=fi,mass=mass,required=faces[fi]['required'],kind=faces[fi]['kind']));continue
        q.update(mass=mass,demand_mass=mass*q['required'],in_current290=fi in selected,
                 saved_historical_sources=[r['source'] for r in historical.get(fi,[])])
        if fi in selected:
            z=selected[fi];assert q['provenance_kind']!='merged_unexpanded'
            assert (z['source_context'],z['canonical_prism_core'],z['raw_family'],z['provenance'])==(q['source_context'],q['canonical_core'],q['raw_family'],q['provenance'])
        for h in historical.get(fi,[]):
            z=h['record'];assert q['provenance_kind']!='merged_unexpanded'
            assert (z['source_context'],z['canonical_prism_core'],z['raw_family'],z['provenance'])==(q['source_context'],q['canonical_core'],q['raw_family'],q['provenance'])
        qrows.append(q)
    assert set(alpha)=={r['face'] for r in qrows+srows}
    numerator=sum((r['mass']*r['required'] for r in qrows+srows),F())
    assert numerator==F(cert['objective_numerator']) and sum(mu.values(),F())==F(cert['normalizer'])
    assert numerator/sum(mu.values(),F())==F(cert['bound'])
    rows=current_domains();local=[];faceunion=set();repaired=dict(mu);origins=set()
    for v,entry in sorted(rows.items()):
        groups=entry['groups'];row=entry['row']
        demands=[sum((alpha.get(f,F())*F(m) for f,m in g['positive_face_coefficients']),F()) for g in groups]
        olddemands=[sum((oa.get(f,F())*F(m) for f,m in g['positive_face_coefficients']),F()) for g in groups]
        assert all(x>0 for x in olddemands)
        loads=[]
        for c in row['certificates']:
            loads.append(exact_charge_load(row['maximal_count_vectors'],demands,c['weights']))
        oldfull=exact_charge_load(row['maximal_count_vectors'],olddemands,row['certificates'][0]['weights'])
        assert oldfull==F(row['certificates'][0]['optimal_mass'])
        lower=max(loads);repaired[v]=max(mu.get(v,F()),lower)
        for g in groups:
            faceunion.update(f for f,m in g['positive_face_coefficients'])
            for o in g['origins']:
                key=(o['face'],o['provenance'],tuple(o['omitted']),o['vertex'],o['raw_family'])
                assert key not in origins;origins.add(key)
        local.append(dict(vertex=v,source=entry['source'],group_source=entry['group_source'],
            forest=row.get('forest',(row.get('union') or {}).get('forest')),
            orbit_count=len(groups),complete_count_vectors=len(row['maximal_count_vectors']),
            demands=demands,old_demands=olddemands,positive_orbits=sum(bool(x) for x in demands),
            old_charge_loads_at_new_alpha=loads,known_valid_charge_lower_bound=lower,
            baseline_vertex_mass=mu.get(v,F()),known_required_increase=max(F(),lower-mu.get(v,F())),
            old_optimum_not_reused=True))
    assert faceunion==set(selected)
    expected={(o['face'],o['provenance'],tuple(o['omitted']),o['vertex'],o['raw_family']) for r in selected.values() for o in r['corners']}
    assert origins==expected
    columns=[];two=0
    for c in model['prepared']['new_columns']:
        load=sum((alpha.get(f,F())*F(m) for f,m in c['integer_face_coefficients']),F())
        cost=sum((mu.get(v,F())*F(m) for v,m in c['integer_vertex_capacities']),F());assert load<=cost
        two+=len(c['integer_vertex_capacities'])>1
        columns.append(dict(variable=c['variable'],vertex_capacities=c['integer_vertex_capacities'],load=load,cost=cost,slack=cost-load))
    assert len(columns)==209 and two==14
    sf={r['face'] for r in read(LOCAL)['positive_simplex_features']}
    newly_historical=[r for r in qrows if r['saved_historical_sources'] and not r['in_current290']]
    counts=dict(positive_faces=len(alpha),newly_positive=len(set(alpha)-set(oa)),dropped=len(set(oa)-set(alpha)),
        changed_common=sum(alpha[i]!=oa[i] for i in set(alpha)&set(oa)),positive_S=len(srows),positive_Q=len(qrows),
        Q_explicit=sum(r['provenance_kind']=='explicit' for r in qrows),
        Q_unique=sum(r['provenance_kind']=='unique_saved_context' for r in qrows),
        Q_merged=sum(r['provenance_kind']=='merged_unexpanded' for r in qrows),
        current_actual_classes=290,positive_current_actual_classes=len(set(selected)&set(alpha)),
        inactive_current_actual_classes=len(set(selected)-set(alpha)),all_historical_classes=len(historical),
        positive_historical_classes=len(set(historical)&set(alpha)),positive_historical_outside290=len(newly_historical),
        forests=len(local),positive_forests=sum(r['positive_orbits']>0 for r in local),
        changed_demands=sum(r['demands']!=r['old_demands'] for r in local),
        known_charge_violations=sum(r['known_required_increase']>0 for r in local),
        S_outside_last_simplex_cover=sum(r['face'] not in sf for r in srows),preserved_two_owner_columns=two)
    tags=tuple(dict.fromkeys((MODEL,OLD,*DOMAIN_STEMS,FORESTS,FORCED,ACTUAL,QUAD,LOCAL,*historical_order)))
    result=dict(status='PASS_NEW5121_SAVED_SUPPORT_AND_PROVENANCE_MIGRATION',sources=bindings(tags),
        model_sha256=model['prepared']['model_sha256'],face_multipliers=cert['face_multipliers'],vertex_multipliers=cert['vertex_multipliers'],
        objective_numerator=numerator,normalizer=F(cert['normalizer']),baseline_bound=F(cert['bound']),counts=counts,
        simplex_rows=srows,quad_rows=qrows,current290_classes=sorted(selected),
        newly_positive_historical_outside290=sorted(newly_historical,key=lambda r:(-r['demand_mass'],r['face'])),
        current_local_domains=local,ranked_known_charge_violations=sorted((r for r in local if r['known_required_increase']),key=lambda r:(-r['known_required_increase'],r['vertex'])),
        new_exported_column_loads=columns,
        known_charge_only_repair=dict(normalizer=sum(repaired.values(),F()),bound=numerator/sum(repaired.values(),F()),
            interpretation='Feasible monotone dual repair for the inherited exported model and these reevaluated old charges only; not a full new-alpha local-cone upper cover.'),
        global_optimizer_invocations=0,local_optimizer_invocations=0,clique_enumerations=0,
        old_alpha_optima_reused=0,new_geometry_asserted=False,
        scope='Saved-only support/provenance migration. Every current290 corner occurrence and complete saved count-vector domain is tracked. Old valid charges are reevaluated at new demands, never treated as new local optima. New positive historical classes and missing S/Q coverage stay explicit. All14 multi-owner exports retain every vertex cost. No complete new-alpha cover or dimension theorem.')
    print(result['status'],counts,'known-charge-only',float(result['known_charge_only_repair']['bound']),flush=True)
    return enc(result)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:
        assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
