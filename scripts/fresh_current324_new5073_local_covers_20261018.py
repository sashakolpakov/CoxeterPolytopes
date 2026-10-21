#!/usr/bin/env python3
"""Separate all600 current324 graph cones at the newly optimized5073 masses.

Only small local LPs propose rational certificates. Saved replay invokes no
optimizer. Immutable complete count-vector domains are transported, not
reenumerated; no physical pair or global-model optimization is performed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import time
import fresh_cross_vertex_support_migration_20261018 as M
import fresh_current_historical_weighted_corners_20261017 as H

ROOT=M.ROOT; W=H.W
TAG='fresh_current324_new5073_local_covers_20261018'
MODEL='fresh_current324_full_model_20261018'
L='fresh_cross_vertex_current290_local_covers_20261018'
U='fresh_cross_vertex_historical32_local_covers_20261018'
C='fresh_new5121_raw4_two_face_local_covers_20261018'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.jsonl')


def jhash(x):
    return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def current_domains():
    """Chronological saved-domain transport; all old pair outcomes inherited."""
    physical=M.current_domains();local=M.read(L);historical=M.read(U);raw4=M.read(C)
    assert (len(physical),len(local['forest_records']),len(historical['forest_records']),len(raw4['forest_records']))==(571,571,161,6)
    latest={}
    for r in local['forest_records']:
        v=r['vertex'];d=physical[v]
        assert r['count_vectors_sha256']==jhash(d['row']['maximal_count_vectors'])
        assert r['literal_representatives_sha256']==jhash(d['row']['literal_representatives'])
        assert r['orbit_coefficients']==[g['positive_face_coefficients'] for g in d['groups']]
        latest[v]=dict(source=L,record=r,geometry_source=d['source'],group_source=d['group_source'],
                       groups=d['groups'],domain=d['row'])
    for tag,doc in ((U,historical),(C,raw4)):
        for r in doc['forest_records']:
            v=r['vertex'];groups=r['domain']['corner_orbits']
            if v in latest:
                old=latest[v];n=r['domain']['old_literal_stars']
                assert n==sum(len(g['literal_stars']) for g in old['groups'])
                assert len(groups)>=len(old['groups'])
                if tag==U:
                    assert r['domain']['baseline_count_vectors_sha256']==jhash(old['domain']['maximal_count_vectors'])
                else:
                    assert r['domain']['baseline']['record_sha256']==jhash(old['record'])
                # Historical32 adds complete new orbits. Raw4 can add a new
                # origin to an identical complete orbit, but cannot remove one.
                for a,b in zip(old['groups'],groups):
                    assert a['literal_stars']==b['literal_stars']
                    assert b['origins'][:len(a['origins'])]==a['origins']
            latest[v]=dict(source=tag,record=r,geometry_source=tag,group_source=tag,groups=groups,domain=r)
    assert len(latest)==600 and len(raw4['inherited_forests'])==594
    for ref in raw4['inherited_forests']:
        d=latest[ref['vertex']]
        assert ref['source']==d['source'] and ref['record_sha256']==jhash(d['record'])
    selected={r['face']:r for r in M.read(M.H)['selected_faces']}
    for r in M.read(M.HD)['added_faces']:selected[r['face']]=r['selected_face']
    for r in M.read(M.RAW125)['selected_faces']:selected[r['face']]=r
    assert len(selected)==290
    for r in historical['added_faces']:selected[r['selected_face']['face']]=r['selected_face']
    for r in raw4['selected_faces']:selected[r['face']]=r
    assert len(selected)==324
    expected={H.origin_key(o) for r in selected.values() for o in r['corners']}
    assert len(expected)==sum(len(r['corners']) for r in selected.values())
    actual=set();faces=set()
    for v,d in latest.items():
        groups=d['groups'];vectors=d['domain']['maximal_count_vectors'];reps=d['domain']['literal_representatives']
        assert [g['index'] for g in groups]==list(range(len(groups)))
        assert len(vectors)==len(reps) and vectors and all(len(x)==len(groups) for x in vectors)
        owner={s:g['index'] for g in groups for s in g['literal_stars']}
        assert len(owner)==sum(len(g['literal_stars']) for g in groups) and set(owner)==set(range(len(owner)))
        for g in groups:
            assert Counter(o['face'] for o in g['origins'])==Counter(dict(g['positive_face_coefficients']))
            faces.update(f for f,m in g['positive_face_coefficients'])
            for o in g['origins']:
                key=H.origin_key(o);assert key not in actual and o['vertex']==v;actual.add(key)
                assert o['literal_orientations']==len(g['literal_stars'])
        for vector,clique in zip(vectors,reps):
            assert len(clique)==len(set(clique)) and all(i in owner for i in clique)
            assert vector==[Counter(owner[i] for i in clique).get(k,0) for k in range(len(groups))]
    assert actual==expected and faces==set(selected)
    return latest,selected


def run(saved=None):
    started=time.monotonic();model=M.read(MODEL);cert=M.dual(MODEL)
    assert model['status']=='PASS_EXACT_CURRENT324_MODEL_CERTIFICATES'
    alpha={i:F(x) for i,x in cert['face_multipliers']};mu={v:F(x) for v,x in cert['vertex_multipliers']}
    assert len(alpha)==5073 and len(mu)==1095
    domains,selected=current_domains();assert len(set(alpha)&set(selected))==240
    tags=tuple(dict.fromkeys((MODEL,L,U,C,*M.DOMAIN_STEMS)))
    sources=M.bindings(tags)
    for mod in (M,H,W):
        p=str(Path(mod.__file__).relative_to(ROOT));sources[p]=M.digest(p)
    sources[f'scripts/{TAG}.py']=M.digest(f'scripts/{TAG}.py')
    key=jhash([sources,cert['face_multipliers']]);cached={}
    if saved is not None:
        assert saved['sources']==sources and saved['face_multipliers']==cert['face_multipliers']
        cached={r['vertex']:r for r in saved['forest_records']};assert len(cached)==len(saved['forest_records'])==600
    elif PROGRESS.exists():
        for line in PROGRESS.read_text().splitlines():
            r=json.loads(line);assert r['key']==key and r['row']['vertex'] not in cached
            cached[r['row']['vertex']]=r['row']
    assert set(cached)<=set(domains)
    print('NEW5073_DOMAIN',len(domains),'forests',len(selected),'classes','cached',len(cached),flush=True)
    records=[];columns=[];newmu=dict(mu);zero=0;solves=0
    order=sorted(domains,key=lambda v:(-mu.get(v,F()),v))
    for step,v in enumerate(order,1):
        d=domains[v];groups=d['groups'];vectors=d['domain']['maximal_count_vectors']
        demand=[sum((alpha.get(f,F())*F(m) for f,m in g['positive_face_coefficients']),F()) for g in groups]
        if v in cached:c=cached[v]['certificate']
        else:
            assert saved is None;c=H.exact_profile(vectors,demand);solves+=bool(any(demand))
        value=W.check_certificate(vectors,demand,c);newmu[v]=max(mu.get(v,F()),value)
        terms=Counter()
        for w,g in zip(map(F,c['weights']),groups):
            for f,m in g['positive_face_coefficients']:terms[f]+=w*F(m)
        coeffs=[[f,str(x)] for f,x in sorted(terms.items()) if x]
        assert sum((alpha.get(f,F())*F(x) for f,x in coeffs),F())==value
        row=M.enc(dict(vertex=v,source=d['source'],source_record_sha256=jhash(d['record']),
            geometry_source=d['geometry_source'],group_source=d['group_source'],
            orbit_coefficients=[g['positive_face_coefficients'] for g in groups],
            complete_count_vectors=len(vectors),count_vectors_sha256=jhash(vectors),
            literal_representatives_sha256=jhash(d['domain']['literal_representatives']),
            demands=demand,certificate=c,baseline_vertex_mass=mu.get(v,F()),
            exact_local_optimum=value,repaired_vertex_mass=newmu[v],required_increase=newmu[v]-mu.get(v,F())))
        if v in cached:assert row==cached[v]
        records.append(row);zero+=not any(demand)
        if value>mu.get(v,F()):
            columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',
                positive_face_coefficients=coeffs,orbit_weights=c['weights'],exact_load=str(value),
                baseline_vertex_mass=str(mu.get(v,F())),source_geometry=d['geometry_source']))
        if saved is None and v not in cached:
            with PROGRESS.open('a') as f:f.write(json.dumps(dict(key=key,row=row),separators=(',',':'))+'\n')
        if step%50==0 or value>mu.get(v,F()):
            print('NEW5073_COVER',step,'/600','violations',len(columns),'seconds',round(time.monotonic()-started,2),flush=True)
    assert all(newmu.get(v,F())>=x for v,x in mu.items())
    inherited=model['prepared']['inherited_several_vertex_columns'];added=model['prepared']['new_columns']
    assert len(inherited)==209 and len(added)==203
    multi=[];loads=[]
    for col in inherited+added:
        load=sum((alpha.get(i,F())*m for i,m in col['integer_face_coefficients']),F())
        before=sum((mu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        after=sum((newmu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        assert load<=before<=after
        row=dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after))
        loads.append(row)
        if len(col['integer_vertex_capacities'])>1:multi.append(row)
    assert len(multi)==14
    num=F(cert['objective_numerator']);den=sum(newmu.values(),F())
    assert num/F(cert['normalizer'])==F(cert['bound']) and sum(mu.values(),F())==F(cert['normalizer'])
    result=dict(status='PASS_NEW5073_COMPLETE_CURRENT324_EXACT_LOCAL_COVERS',sources=sources,
        model_sha256=model['prepared']['model_sha256'],face_multipliers=cert['face_multipliers'],
        initial_vertex_multipliers=cert['vertex_multipliers'],vertex_multipliers=[[v,str(x)] for v,x in sorted(newmu.items()) if x],
        initial_bound=cert['bound'],objective_numerator=str(num),normalizer=str(den),bound=str(num/den),
        represented_classes=sorted(selected),forest_records=records,columns=columns,
        preserved_multi_owner_columns=multi,checked_recent_model_column_loads=loads,
        counts=dict(represented_classes=324,positive_represented_classes=240,forests=600,
            exact_new_demand_local_optima=600,zero_demand_forests=zero,nonzero_demand_forests=600-zero,
            new_charge_columns=len(columns),raised_vertex_masses=len(columns),
            complete_count_vectors=sum(r['complete_count_vectors'] for r in records),
            orbit_coordinates=sum(len(r['demands']) for r in records),
            mixture_terms=sum(len(r['certificate']['mixture']) for r in records)),
        generation_local_LPs=600-zero,replay_local_LPs=0,global_optimizer_invocations=0,
        new_physical_pair_tests=0,clique_enumerations=0,old_mass_optima_reused=0,
        preservation='All5073 face masses are fixed. Every one of600 complete current324 conservative graph cones receives a new-demand exact charge/cover certificate. Only vertex multipliers increase, preserving all439822 inherited full-model columns, including signed/zero-budget records; all14 two-owner costs are explicitly checked. Graphs, native domains, owners and unknown outcomes are unchanged.',
        scope='Complete existing324-class/600-forest local graph-cone separation at current5073 masses, with chronological frozen count-vector/origin transport. No physical pair or clique recensus, new class geometry, global optimization, target12 primal, physical realization or dimension theorem.')
    print(result['status'],result['counts'],'bound',float(num/den),'new_local_solves_this_run',solves,
          'seconds',round(time.monotonic()-started,2),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    if not args.verify:assert not path.exists(),'Output exists; use --verify.'
    result=run(saved)
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
