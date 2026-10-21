#!/usr/bin/env python3
"""Screen all missing raw76 masses, then test the best two-class batch.

Preserve every old pair outcome; only new-node pairs are tested. Complete
changed graph cones receive exact new-demand charges/covers. No global LP.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time
import fresh_new5073_raw76_support_frontier_20261018 as Q
import fresh_new5121_raw4_two_face_local_covers_20261018 as Z

L=Q.L;M=L.M;A=Q.A;R=Q.R;H=L.H;W=L.W;ROOT=Q.ROOT
TAG='fresh_new5073_raw76_high_mass_local_covers_20261018'
FACES=(42980,43028)
PROGRESS=Path('/private/tmp')/(TAG+'_progress.jsonl')
normal=Q.normal


def screen(inp,local):
    alpha={f:F(x) for f,x in local['face_multipliers']};mu={v:F(x) for v,x in local['vertex_multipliers']}
    old={r['vertex']:F(r['exact_local_optimum']) for r in local['forest_records']}
    represented=set(local['represented_classes']);forests=inp.read(M.FORESTS);quads=inp.read(M.QUAD)
    spectators=inp.read(A.B.Q.SPECTATORS);model=inp.read(L.MODEL)
    provenance=inp.read(M.ACTUAL)['provenance_masses'];pmap={f:i for i,f in model['prepared']['provenance_to_face']}
    lookup={tuple(r['forest']):r['index'] for r in forests['vertices']};core=quads['canonical_cores'][65]
    assert core['source_core']==76 and len(core['aggregate_corner_occurrences'])==4
    rows=[];missing=0
    for fi,mass in sorted(alpha.items()):
        if fi in represented:continue
        if fi in pmap:
            pi=pmap[fi];r=provenance[pi];ci=r['canonical_prism_core'];ctx=r['source_context'];pk='explicit'
        elif fi<len(forests['faces']):
            f=forests['faces'][fi]
            if f['kind']!='canonical_dotted_quad_aggregate' or len(f['source_contexts'])!=1:continue
            ctx,ci=f['source_contexts'][0];pi=None;pk='unique_saved_context'
        else:continue
        missing+=1
        if ci!=65:continue
        assert pi is None and pk=='unique_saved_context'
        context=spectators['contexts'][ctx];assert quads['source_core_maps'][context['core']]['canonical_core']==65
        spectator=context['spectator_forest'];multiplicities=Counter()
        for c in core['aggregate_corner_occurrences']:
            v=lookup[tuple(sorted(spectator+c['forest']))];multiplicities[v]+=1
        assert sum(multiplicities.values())==4
        rows.append(dict(face=fi,mass=mass,required=1,demand_mass=mass,raw_family=76,canonical_core=65,
            source_context=ctx,provenance=None,provenance_kind=pk,saved_historical_sources=[],
            spectator_forest=spectator,corner_multiplicities=sorted(multiplicities.items())))
    assert missing==1519 and len(rows)==22
    num,den=F(local['objective_numerator']),F(local['normalizer']);need=num/12-den
    def potential(faces):
        extra=Counter()
        for r in faces:
            for v,n in r['corner_multiplicities']:extra[v]+=r['mass']*n
        increases={v:max(mu.get(v,F()),old.get(v,F())+e)-mu.get(v,F()) for v,e in extra.items()}
        gain=sum(increases.values(),F())
        return dict(faces=[r['face'] for r in faces],extra_demands=sorted(extra.items()),
            maximum_mu_increase=gain,ratio_lower_bound=num/(den+gain),
            passes_necessary_gain_threshold=gain>=need)
    singles=[potential([r]) for r in rows];pairs=[potential(list(rs)) for rs in combinations(rows,2)]
    assert not any(r['passes_necessary_gain_threshold'] for r in singles)
    feasible=[r for r in pairs if r['passes_necessary_gain_threshold']];assert len(feasible)==4
    best=max(pairs,key=lambda r:(r['maximum_mu_increase'],[-f for f in r['faces']]))
    assert tuple(best['faces'])==FACES
    selected=[next(r for r in rows if r['face']==f) for f in FACES]
    return dict(positive_unambiguous_outside324=missing,raw76_classes=len(rows),rows=rows,
        single_class_ceilings=singles,two_class_ceilings=pairs,passing_two_class_batches=feasible,
        best_two_class_batch=best,all22_ceiling=potential(rows),required_mu_increment=need,
        minimal_potential_batch_size_within_raw76=2,
        selection='Among all22 currently positive unambiguous raw76 classes, choose the two-class batch with largest exact singleton/slack upper gain. Passing this necessary threshold is not a prediction of actual pair exclusions or local optimum.'),selected,quads,alpha


def current(v,p,local,domains,patches):
    nodes,groups,g,owner,forest,unused=Z.current(v,p)
    nodes=normal(nodes);groups=normal(groups)
    if v in patches:
        row=patches[v];d=row['domain'];oldn=len(nodes)
        assert d['literal_stars'][:oldn]==nodes and d['corner_orbits'][:len(groups)]==groups
        nodes=normal(d['literal_stars']);groups=normal(d['corner_orbits'])
        g.add_nodes_from(range(oldn,len(nodes)))
        pairs=[(i,j) for i,j in combinations(range(len(nodes)),2) if j>=oldn]
        assert [tuple(z[:2]) for z in row['newly_owned_pair_domain']]==pairs
        for i,j,code in row['newly_owned_pair_domain']:
            if code is None:g.add_edge(i,j)
        assert L.jhash(sorted(tuple(sorted(e)) for e in g.edges()))==row['graph_edges_sha256']
    d=domains[v];row=local[v]
    assert groups==normal(d['groups']) and len(nodes)==sum(len(z['literal_stars']) for z in groups)
    assert row['source_record_sha256']==L.jhash(d['record'])
    assert row['count_vectors_sha256']==L.jhash(d['domain']['maximal_count_vectors'])
    assert row['orbit_coefficients']==[z['positive_face_coefficients'] for z in groups]
    assert W.check_certificate(d['domain']['maximal_count_vectors'],list(map(F,row['demands'])),row['certificate'])==F(row['exact_local_optimum'])
    owner={i:z['index'] for z in groups for i in z['literal_stars']}
    assert set(owner)==set(range(len(nodes)))
    assert all(all(g.has_edge(i,j) for i,j in combinations(cl,2)) for cl in d['domain']['literal_representatives'])
    return nodes,groups,g,owner,forest,dict(source=L.TAG,record_sha256=L.jhash(row),
        exact_local_optimum=row['exact_local_optimum'],graph_edges_sha256=L.jhash(sorted(tuple(sorted(e)) for e in g.edges())))


def merge(raw,p,local,domains,patches,alpha):
    v=raw['vertex'];nodes,groups,g,owner,forest,baseline=current(v,p,local,domains,patches)
    assert forest==raw['forest'] and normal(R.P.forest_graph(forest,dict(R.C.models())))==raw['target']
    oldn=len(nodes);oldedges={tuple(sorted(e)) for e in g.edges()};oldgroups=normal(groups)
    keys={H.frozen(s['contacts']):i for i,s in enumerate(nodes)};assert len(keys)==oldn
    mapping=[]
    for s in raw['literal_stars']:
        native=s['native'];k=H.frozen(native['contacts'])
        if k in keys:
            i=keys[k];assert all(nodes[i][f]==native[f] for f in H.NATIVE_FIELDS)
        else:
            i=len(nodes);keys[k]=i;nodes.append(dict(index=i,**normal(native)));g.add_node(i)
        mapping.append(i)
    assert len(mapping)==len(set(mapping))
    orbitkeys={frozenset(z['literal_stars']):z['index'] for z in groups}
    origins={H.origin_key(o) for z in groups for o in z['origins']};added=[]
    for occurrence in raw['occurrences']:
        o=occurrence['origin'];ok=H.origin_key(o);assert ok not in origins;origins.add(ok)
        ids=sorted(mapping[i] for i in occurrence['literal_star_indices']);assert len(ids)==len(set(ids))==o['literal_orientations']
        sk=frozenset(ids)
        if sk in orbitkeys:gi=orbitkeys[sk]
        else:
            assert not any(i in owner for i in ids),'partial native-orbit overlap'
            gi=len(groups);orbitkeys[sk]=gi
            for i in ids:owner[i]=gi
            groups.append(dict(index=gi,face=o['face'],literal_stars=ids,origins=[],positive_face_coefficients=[]))
        z=groups[gi];z['origins'].append(normal(o));coeff=Counter(dict(z['positive_face_coefficients']));coeff[o['face']]+=1
        z['positive_face_coefficients']=[[f,m] for f,m in sorted(coeff.items())]
        added.append(dict(origin=normal(o),orbit=gi,literal_star_indices=ids))
    assert set(owner)==set(range(len(nodes))) and sum(len(z['literal_stars']) for z in groups)==len(nodes)
    for old,new in zip(oldgroups,groups):
        assert old['literal_stars']==new['literal_stars'] and new['origins'][:len(old['origins'])]==old['origins']
    assert all(Counter(o['face'] for o in z['origins'])==Counter(dict(z['positive_face_coefficients'])) for z in groups)
    pairs=[(i,j) for i,j in combinations(range(len(nodes)),2) if j>=oldn];g.add_edges_from(pairs)
    demand=[sum((alpha.get(f,F())*m for f,m in z['positive_face_coefficients']),F()) for z in groups]
    return dict(vertex=v,forest=forest,baseline=baseline,literal_stars=nodes,corner_orbits=groups,
        raw76_to_merged_star_indices=mapping,added_occurrences=added,demands=list(map(str,demand)),
        old_literal_stars=oldn,new_literal_stars=len(nodes)-oldn,new_pair_positions=len(pairs)),g,owner,pairs,oldedges,demand


def run(saved=None,inventory_only=False):
    started=time.monotonic();inp=Q.Inputs();base=inp.read(L.TAG)
    screened,selected,quads,alpha=screen(inp,base)
    raw=normal(M.enc(Q.geometry(inp,[{k:v for k,v in r.items() if k!='spectator_forest'} for r in selected],quads,alpha)))
    for r in raw['selected_faces']:
        assert Counter(c['vertex'] for c in r['corners'])==Counter(dict(r['corner_multiplicities']))
    domains,represented=L.current_domains();p=Z.prepare();patches={r['vertex']:r for r in M.read(Z.TAG)['forest_records']}
    records0={r['vertex']:r for r in base['forest_records']};mu={v:F(x) for v,x in base['vertex_multipliers']}
    prepared=[merge(r,p,records0,domains,patches,alpha) for r in raw['forest_records']]
    assert len(prepared)==6 and all(u['new_pair_positions']<=60000 for u,*_ in prepared)
    totalpairs=sum(u['new_pair_positions'] for u,*_ in prepared)
    print('RAW76_BATCH_PREPARED',len(selected),'classes',len(prepared),'forests',
          sum(u['new_literal_stars'] for u,*_ in prepared),'newstars',totalpairs,'newpairs',
          'seconds',round(time.monotonic()-started,2),flush=True)
    for u,*_ in prepared:print('RAW76_TARGET',u['vertex'],u['old_literal_stars'],u['new_literal_stars'],u['new_pair_positions'],flush=True)
    if inventory_only:return None
    paths=set(inp.paths)
    for s in (L.TAG,Z.TAG,Z.U.TAG,L.MODEL):
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):paths.add(f'{folder}/{s}.{ext}')
    for mod in (Q,L,M,A,R,H,W,Z,Z.U,Z.T):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(str(Path(__file__).relative_to(ROOT)));sources={f:M.digest(f) for f in sorted(paths)}
    key=L.jhash(sources);cache={}
    if saved:
        assert saved['sources']==sources;cache={r['vertex']:r for r in saved['forest_records']}
        assert len(cache)==len(saved['forest_records'])==6
    elif PROGRESS.exists():
        for line in PROGRESS.read_text().splitlines():
            z=json.loads(line);assert z['key']==key and z['row']['vertex'] not in cache;cache[z['row']['vertex']]=z['row']
    out=[];newmu=dict(mu);columns=[];counts=Counter()
    for u,g,owner,pairs,oldedges,demand in prepared:
        v=u['vertex'];print('RAW76_PHYSICAL_START',v,len(pairs),flush=True)
        physical=Z.physical(u,g,pairs)
        assert {tuple(sorted(e)) for e in g.edges() if max(e)<u['old_literal_stars']}==oldedges
        enum,error=H.vectors_for(g,owner,len(demand));assert error is None,(v,error)
        vectors,reps,ncliques,npatterns=enum;vectors=list(map(list,vectors));reps=list(map(list,reps))
        c=cache[v]['certificate'] if v in cache else H.exact_profile(vectors,demand)
        value=W.check_certificate(vectors,demand,c);assert value>=F(u['baseline']['exact_local_optimum'])
        newmu[v]=max(mu.get(v,F()),value);co=Counter()
        for w,z in zip(map(F,c['weights']),u['corner_orbits']):
            for f,m in z['positive_face_coefficients']:co[f]+=w*m
        coefficients=[[f,str(x)] for f,x in sorted(co.items()) if x]
        assert sum(alpha.get(f,F())*F(x) for f,x in coefficients)==value
        row=normal(M.enc(dict(vertex=v,domain=u,physical=physical,graph_edges_sha256=L.jhash(sorted(tuple(sorted(e)) for e in g.edges())),
            maximal_cliques=ncliques,count_patterns=npatterns,maximal_count_vectors=vectors,literal_representatives=reps,
            certificate=c,exact_local_optimum=value,baseline_vertex_mass=mu.get(v,F()),repaired_vertex_mass=newmu[v],
            multiplier_increase=newmu[v]-mu.get(v,F()),positive_face_coefficients=coefficients)))
        if v in cache:assert row==cache[v]
        elif not saved:
            with PROGRESS.open('a') as f:f.write(json.dumps(dict(key=key,row=row),separators=(',',':'))+'\n')
        out.append(row)
        if value>mu.get(v,F()):columns.append(dict(vertex=v,vertex_cost='1',positive_face_coefficients=coefficients,exact_load=str(value)))
        counts.update(new_pairs=len(pairs),new_exclusions=physical['new_exclusions'],new_unknowns=physical['new_unknown_pairs'],
            new_stars=u['new_literal_stars'],actual_corners=len(u['added_occurrences']),maximal_cliques=ncliques,
            count_vectors=len(vectors),mixture_terms=len(c['mixture']))
        print('RAW76_DONE',v,'excluded',physical['new_exclusions'],'unknown',physical['new_unknown_pairs'],
            'vectors',len(vectors),'mu_increase',str(newmu[v]-mu.get(v,F())),'seconds',round(time.monotonic()-started,2),flush=True)
    model=M.read(L.MODEL);loads=[];two=0
    for col in model['prepared']['inherited_several_vertex_columns']+model['prepared']['new_columns']:
        load=sum(alpha.get(f,F())*m for f,m in col['integer_face_coefficients'])
        before=sum(mu.get(v,F())*m for v,m in col['integer_vertex_capacities'])
        after=sum(newmu.get(v,F())*m for v,m in col['integer_vertex_capacities']);assert load<=before<=after
        two+=len(col['integer_vertex_capacities'])>1
        loads.append(dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after)))
    assert len(loads)==412 and two==14
    inherited=[dict(vertex=v,record_sha256=L.jhash(r)) for v,r in sorted(records0.items()) if v not in {u['vertex'] for u,*_ in prepared}]
    assert len(inherited)==594 and len(represented)==324 and not set(FACES)&set(represented)
    num,den=F(base['objective_numerator']),sum(newmu.values(),F())
    result=M.enc(dict(status='PASS_NEW5073_RAW76_HIGH_MASS_EXACT_LOCAL_COVERS',sources=sources,screen=screened,
        selected_faces=raw['selected_faces'],raw_domain=raw,forest_records=out,inherited_forests=inherited,
        face_multipliers=base['face_multipliers'],vertex_multipliers=sorted((v,x) for v,x in newmu.items() if x),
        baseline_ratio=base['bound'],objective_numerator=num,normalizer=den,bound=num/den,columns=columns,
        counts=dict(counts),represented_classes=326,represented_forests=600,changed_forests=6,
        checked_recent_model_loads=loads,preserved_multi_owner_columns=14,
        global_optimizer_invocations=0,new_dimension_bound=None,
        scope='Two complete raw76 actual classes added after all22 exact potential screening. All old-old graph outcomes and full multi-owner costs preserved; only new-node pairs tested. Six complete graph cones have exact charge/cover optima at fixed5073 masses;594 unchanged optima inherited. Monotone repair preserves439822 old model columns plus prior25 local exports. No global optimization, target12 primal, physical realization or dimension theorem.'))
    print(result['status'],result['counts'],'ratio',float(num/den),flush=True)
    return normal(result)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--inventory',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    if not args.verify and not args.inventory:assert not path.exists()
    out=run(saved,args.inventory)
    if args.inventory:pass
    elif args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
