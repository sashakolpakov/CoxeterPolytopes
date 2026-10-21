#!/usr/bin/env python3
"""Merge two exact raw4 domains into current322 cones at fixed5121 masses.

All old graph outcomes are immutable. Only pairs involving newly added native
stars receive physical tests. Six complete local cones, no global optimizer.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_cross_vertex_historical32_local_covers_20261018 as U
import fresh_new5121_raw4_two_face_domain_20261018 as Q

ROOT=U.ROOT;M=U.M;L=U.L;T=U.T;N=U.N;H=U.H;W=U.W;R=N.R
TAG='fresh_new5121_raw4_two_face_local_covers_20261018'
normal=U.normal


def prepare():
    p=U.prepare();base=M.read(U.TAG);raw=M.read(Q.TAG)
    assert base['represented_classes']==322 and base['represented_forests']==600
    assert base['face_multipliers']==p['base']['face_multipliers']
    assert raw['counts']==dict(faces=2,actual_corner_occurrences=12,vertex_types=6,
        active_embeddings=56,literal_orientation_occurrences=40,distinct_native_stars=28)
    assert [r['face'] for r in raw['selected_faces']]==[14792,14927]
    represented=set(M.read(M.TAG)['current290_classes'])|{r['selected_face']['face'] for r in base['added_faces']}
    assert len(represented)==322 and not represented&{14792,14927}
    alpha={i:F(a) for i,a in base['face_multipliers']}
    for r in raw['selected_faces']:
        assert alpha[r['face']]==F(r['mass']) and F(r['demand_mass'])==3*alpha[r['face']]
        assert r['required']==3 and len(r['corners'])==6
        assert len({tuple(c['omitted']) for c in r['corners']})==6
    sources=M.bindings((U.TAG,Q.TAG,L.TAG,M.MODEL))
    for mod in (U,Q,T,N,H,W,R):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=M.digest(path)
    sources[f'scripts/{TAG}.py']=M.digest(f'scripts/{TAG}.py')
    return dict(previous=p,base=base,raw=raw,alpha=alpha,
        mu={i:F(x) for i,x in base['vertex_multipliers']},sources=sources,
        rows32={r['vertex']:r for r in base['forest_records']},
        rows290={r['vertex']:r for r in p['base']['forest_records']})


def current(v,p):
    """Recover old native stars/graph without retesting any old physical pair."""
    up=p['previous']
    old,g,owner,unowned,demands,nodes,groups,forest,refs,source=T.current(v,up['p'],up['changed'])
    nodes=normal(nodes);groups=normal(groups)
    if v in up['patches']:
        patch=up['patches'][v]
        assert patch['literal_stars']==nodes and patch['corner_orbits']==groups
        for i,j,code in patch['newly_owned_pair_domain']:
            assert g.has_edge(i,j)
            if code is not None:g.remove_edge(i,j)
    assert groups==up['latest'][v]['groups']
    if v in p['rows32']:
        row=p['rows32'][v];d=row['domain'];oldn=len(nodes)
        assert d['old_literal_stars']==oldn and d['old_corner_orbits']==len(groups)
        assert d['old_edges_sha256']==L.jhash(sorted(tuple(sorted(e)) for e in g.edges()))
        assert d['corner_orbits'][:len(groups)]==groups and d['forest']==forest
        assert d['baseline_count_vectors_sha256']==L.jhash(up['latest'][v]['row']['maximal_count_vectors'])
        for tag,oi in d['new_literal_references']:
            s=normal(W.forests(tag)[v]['literal_stars'][oi]);assert s['index']==oi
            nodes.append(s)
        assert len(nodes)==d['counts']['literal_stars']
        assert L.jhash([{k:s[k] for k in H.NATIVE_FIELDS} for s in nodes])==d['native_stars_sha256']
        g.add_nodes_from(range(oldn,len(nodes)))
        for i,j in combinations(range(len(nodes)),2):
            if j>=oldn:g.add_edge(i,j)
        for i,j in d['new_excluded_pairs']:
            assert j>=oldn and g.has_edge(i,j);g.remove_edge(i,j)
        groups=d['corner_orbits'];vectors=row['maximal_count_vectors'];reps=row['literal_representatives']
        certificate=row['certificate'];origin=U.TAG
    else:
        row=p['rows290'][v];vectors=up['latest'][v]['row']['maximal_count_vectors']
        reps=up['latest'][v]['row']['literal_representatives'];certificate=row['certificate'];origin=L.TAG
        assert row['count_vectors_sha256']==L.jhash(vectors)
        assert row['literal_representatives_sha256']==L.jhash(reps)
    groups=normal(groups);owner={i:z['index'] for z in groups for i in z['literal_stars']}
    assert set(owner)==set(range(len(nodes))) and sum(len(z['literal_stars']) for z in groups)==len(nodes)
    demand=[sum((p['alpha'].get(f,F())*F(m) for f,m in z['positive_face_coefficients']),F()) for z in groups]
    assert list(map(str,demand))==(row['domain']['demands'] if origin==U.TAG else row['demands'])
    value=W.check_certificate(vectors,demand,certificate)
    assert value==F(row['exact_local_optimum'])<=p['mu'].get(v,F())
    assert len(vectors)==len(reps)
    for vec,cl in zip(vectors,reps):
        assert all(g.has_edge(i,j) for i,j in combinations(cl,2))
        assert vec==[Counter(owner[i] for i in cl).get(j,0) for j in range(len(groups))]
    return nodes,groups,g,owner,forest,dict(source=origin,record_sha256=L.jhash(row),
        graph_edges_sha256=L.jhash(sorted(tuple(sorted(e)) for e in g.edges())),
        exact_local_optimum=str(value),literal_stars=len(nodes),corner_orbits=len(groups))


def merge(raw,p):
    v=raw['vertex'];nodes,groups,g,owner,forest,baseline=current(v,p)
    assert raw['forest']==forest
    assert raw['target']==normal(R.P.forest_graph(forest,dict(R.C.models())))
    oldn=len(nodes);oldedges={tuple(sorted(e)) for e in g.edges()};oldgroups=normal(groups)
    keys={H.frozen(s['contacts']):i for i,s in enumerate(nodes)};assert len(keys)==oldn
    mapping=[]
    for s in raw['literal_stars']:
        assert s['index']==len(mapping);native=s['native'];key=H.frozen(native['contacts'])
        if key in keys:
            i=keys[key];assert all(nodes[i][k]==native[k] for k in H.NATIVE_FIELDS)
        else:
            i=len(nodes);keys[key]=i;nodes.append(dict(index=i,**normal(native)));g.add_node(i)
        mapping.append(i)
    assert len(mapping)==len(set(mapping))
    orbitkeys={frozenset(z['literal_stars']):z['index'] for z in groups}
    origins={H.origin_key(o) for z in groups for o in z['origins']};added=[]
    for occurrence in raw['occurrences']:
        o=occurrence['origin'];okey=H.origin_key(o);assert okey not in origins;origins.add(okey)
        ids=sorted(mapping[i] for i in occurrence['literal_star_indices'])
        assert len(ids)==len(set(ids))==o['literal_orientations']
        os=frozenset(ids)
        if os in orbitkeys:gi=orbitkeys[os]
        else:
            assert not any(i in owner for i in ids),'partial native orbit overlap'
            gi=len(groups);orbitkeys[os]=gi
            for i in ids:owner[i]=gi
            groups.append(dict(index=gi,face=o['face'],literal_stars=ids,origins=[],positive_face_coefficients=[]))
        group=groups[gi];group['origins'].append(normal(o))
        co=Counter(dict(group['positive_face_coefficients']));co[o['face']]+=1
        group['positive_face_coefficients']=[[f,m] for f,m in sorted(co.items())]
        added.append(dict(origin=normal(o),orbit=gi,literal_star_indices=ids))
    assert set(owner)==set(range(len(nodes))) and sum(len(z['literal_stars']) for z in groups)==len(nodes)
    for z in groups:
        assert z['index']<len(groups)
        assert Counter(o['face'] for o in z['origins'])==Counter(dict(z['positive_face_coefficients']))
    for i,z in enumerate(oldgroups):
        assert groups[i]['literal_stars']==z['literal_stars']
        assert groups[i]['origins'][:len(z['origins'])]==z['origins']
    pairs=[(i,j) for i,j in combinations(range(len(nodes)),2) if j>=oldn]
    g.add_edges_from(pairs)
    demand=[sum((p['alpha'].get(f,F())*F(m) for f,m in z['positive_face_coefficients']),F()) for z in groups]
    return dict(vertex=v,forest=forest,baseline=baseline,literal_stars=nodes,corner_orbits=groups,
        raw4_to_merged_star_indices=mapping,added_occurrences=added,demands=list(map(str,demand)),
        old_literal_stars=oldn,new_literal_stars=len(nodes)-oldn,new_pair_positions=len(pairs)),g,owner,pairs,oldedges,demand


def physical(u,g,pairs):
    nodes=u['literal_stars'];pp=[N.native(s) for s in nodes]
    target=R.P.forest_graph(u['forest'],dict(R.C.models()));J,active,components=R.projection_data(target)
    checked=sorted({i for pair in pairs for i in pair});heights={};native_checks=[]
    for i in checked:
        star=nodes[i];cov=pp[i];edge=tuple(star['native_omitted_product']);field=N.pair_field(cov,(edge,));hs=[]
        assert all(j in active for c in cov for j,e in c)
        for c in cov:
            z=R.projection(c,c,J,field);assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;hs.append(h)
        assert [h.out() for h in hs]==star['squared_timelike_heights']
        cross=R.projection(cov[0],cov[1],J,field);native=field.negative_contact(edge)
        assert field.signs_difference(tuple(a-b for a,b in zip(cross,native)),hs[0]*hs[1])==0
        heights[i]=hs;native_checks.append(dict(index=i,projection=R.qout(cross),
            pair_field_radicals=[field.tau.out(),field.sigma.out()],squared_timelike_heights=[h.out() for h in hs]))
    cache={};expressions=[];exprids={};tested=[];bad=[];reasons=Counter()
    for i,j in pairs:
        shared=set(map(tuple,nodes[i]['occupied_pairs']))&set(map(tuple,nodes[j]['occupied_pairs']))
        why=[i,j,0] if shared else None
        if why is None:
            for a,b in ((0,0),(0,1),(1,0),(1,1)):
                pi,pj=pp[i][a],pp[j][b]
                if pi==pj:continue
                ck=tuple(sorted((pi,pj)))
                if ck not in cache:
                    field=N.pair_field((pi,pj));z=R.projection(pi,pj,J,field);h=heights[i][a]*heights[j][b]
                    sg=field.signs_difference(z,h);sg1=field.signs_difference(R.qaddk(z,1),h)
                    code=1 if sg>0 else 2 if sg<0 and sg1>0 else None
                    expr=dict(pair_field_radicals=[field.tau.out(),field.sigma.out()],projection=R.qout(z),
                        height_product=h.out(),sign_product=sg,sign_product_plus_one=sg1)
                    if code is None:
                        code,extra=field.arithmetic_wall_gate(tuple((x.a,x.b) for x in z),h.a,h.b)
                        if code is not None:expr['arithmetic_gate']=extra
                    cache[ck]=(code,expr)
                code,expr=cache[ck]
                if code is not None:
                    key=json.dumps(expr,sort_keys=True,separators=(',',':'))
                    if key not in exprids:exprids[key]=len(expressions);expressions.append(expr)
                    why=[i,j,code,a,b,exprids[key]];break
        if why is None:tested.append([i,j,None])
        else:
            bad.append(why);tested.append([i,j,why[2]]);reasons[why[2]]+=1;g.remove_edge(i,j)
    return dict(newly_owned_pair_domain=tested,incompatibility_certificates=bad,
        wall_expressions=expressions,reason_counts=dict(reasons),native_checks=native_checks,
        projection_components=components,new_exclusions=len(bad),new_unknown_pairs=len(pairs)-len(bad))


def run(saved=None,inventory=False):
    start=time.monotonic();p=prepare();cached={r['vertex']:r for r in saved['forest_records']} if saved else {}
    if saved:assert saved['sources']==p['sources']
    records=[];columns=[];newmu=dict(p['mu']);counts=Counter()
    for raw in p['raw']['forest_records']:
        u,g,owner,pairs,oldedges,demand=merge(raw,p);v=u['vertex']
        print('RAW4_DOMAIN',v,'old',u['old_literal_stars'],'new',u['new_literal_stars'],
            'orbits',len(u['corner_orbits']),'newpairs',len(pairs),'seconds',round(time.monotonic()-start,2),flush=True)
        counts.update(new_literal_stars=u['new_literal_stars'],new_actual_corners=len(u['added_occurrences']),new_pair_tests=len(pairs))
        if inventory:continue
        physical_data=physical(u,g,pairs)
        assert {tuple(sorted(e)) for e in g.edges() if max(e)<u['old_literal_stars']}==oldedges
        assert len(pairs)-sum(g.has_edge(i,j) for i,j in pairs)==physical_data['new_exclusions']
        enum,error=H.vectors_for(g,owner,len(demand));assert error is None,(v,error)
        vectors,reps,ncliques,npatterns=enum;vectors=list(map(list,vectors));reps=list(map(list,reps))
        for vec,cl in zip(vectors,reps):
            assert all(g.has_edge(i,j) for i,j in combinations(cl,2))
            assert vec==[Counter(owner[i] for i in cl).get(j,0) for j in range(len(demand))]
        c=cached[v]['certificate'] if saved else H.exact_profile(vectors,demand)
        value=W.check_certificate(vectors,demand,c);assert value>=F(u['baseline']['exact_local_optimum'])
        newmu[v]=max(newmu.get(v,F()),value)
        co=Counter()
        for w,z in zip(map(F,c['weights']),u['corner_orbits']):
            for f,m in z['positive_face_coefficients']:co[f]+=w*F(m)
        coeffs=[[f,str(x)] for f,x in sorted(co.items()) if x]
        assert sum((p['alpha'].get(f,F())*F(x) for f,x in coeffs),F())==value
        if value>p['mu'].get(v,F()):columns.append(dict(vertex=v,vertex_cost='1',positive_face_coefficients=coeffs,
            orbit_weights=c['weights'],exact_load=str(value)))
        row=normal(dict(vertex=v,domain=u,**physical_data,graph_edges_sha256=L.jhash(sorted(tuple(sorted(e)) for e in g.edges())),
            maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,maximal_count_vectors=vectors,
            literal_representatives=reps,certificate=c,baseline_vertex_mass=str(p['mu'].get(v,F())),
            exact_local_optimum=str(value),required_increase=str(newmu[v]-p['mu'].get(v,F()))))
        if saved:assert row==cached[v]
        records.append(row);counts.update(new_exclusions=physical_data['new_exclusions'],new_unknown_pairs=physical_data['new_unknown_pairs'])
        print('RAW4_COVER',v,'excluded',physical_data['new_exclusions'],'unknown',physical_data['new_unknown_pairs'],
            'cliques',ncliques,'vectors',len(vectors),'increase',row['required_increase'],'seconds',round(time.monotonic()-start,2),flush=True)
    assert counts['new_actual_corners']==12
    if inventory:print('INVENTORY_COMPLETE',dict(counts),flush=True);return
    assert len(records)==6
    affected={r['vertex'] for r in records};inherited=[]
    for v in sorted(set(p['previous']['latest'])|set(p['rows32'])):
        if v in affected:continue
        row=p['rows32'].get(v,p['rows290'].get(v));assert row is not None
        inherited.append(dict(vertex=v,source=U.TAG if v in p['rows32'] else L.TAG,record_sha256=L.jhash(row)))
    assert len(inherited)==594
    multi=[]
    for col in M.read(M.MODEL)['prepared']['new_columns']:
        load=sum((p['alpha'].get(i,F())*m for i,m in col['integer_face_coefficients']),F())
        before=sum((p['mu'].get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        after=sum((newmu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        assert load<=before<=after
        if len(col['integer_vertex_capacities'])>1:multi.append(dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after)))
    assert len(multi)==14 and all(newmu.get(v,F())>=x for v,x in p['mu'].items())
    num=F(p['base']['objective_numerator']);den=sum(newmu.values(),F())
    out=dict(status='PASS_NEW5121_RAW4_TWO_FACE_EXACT_LOCAL_COVERS',sources=p['sources'],
        face_multipliers=p['base']['face_multipliers'],initial_vertex_multipliers=p['base']['vertex_multipliers'],
        vertex_multipliers=[[v,str(x)] for v,x in sorted(newmu.items()) if x],objective_numerator=str(num),
        normalizer=str(den),bound=str(num/den),baseline_bound=p['base']['bound'],selected_faces=p['raw']['selected_faces'],
        forest_records=records,inherited_forests=inherited,columns=columns,preserved_multi_owner_columns=multi,counts=dict(counts),
        represented_classes=324,positive_represented_classes=261,represented_forests=600,
        changed_forests=6,new_forests=0,inherited_forest_count=594,new_charge_columns=len(columns),
        complete_count_vectors=sum(len(r['maximal_count_vectors']) for r in records),
        complete_maximal_cliques=sum(r['maximal_literal_cliques'] for r in records),
        generation_local_LPs=6,replay_local_LPs=0,global_optimizer_invocations=0,
        scope='Complete current324-class/600-forest conservative local graph cones at fixed5121 masses. Only six affected graphs are rebuilt; every old-old edge and exclusion, including all owned unknowns, is preserved literally. New physical tests use immutable native covectors, per-pair fields, positive future heights, port uniqueness and exact necessary wall gates; unknowns remain allowed. Actual corner multiplicities, all14 two-owner costs and prior model columns survive monotone vertex-mass repair. No complete missing S/Q coverage, global optimum, target12 feasibility, physical realization or dimension theorem.')
    print(out['status'],dict(counts),'newcolumns',len(columns),'bound',float(num/den),'seconds',round(time.monotonic()-start,2),flush=True)
    return normal(out)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--inventory',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    out=run(saved,args.inventory)
    if args.inventory:pass
    elif args.verify:assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,separators=(',',':'))+'\n')
