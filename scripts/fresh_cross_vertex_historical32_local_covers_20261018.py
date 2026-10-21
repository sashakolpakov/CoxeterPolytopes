#!/usr/bin/env python3
"""Merge32 newly positive historical classes at the frozen5121 face masses.

No physical pair arithmetic or global solve. Preserve the complete old graph;
new pairs take their first complete historical owner, including unknowns.
Bounded clique enumeration and local LPs produce exact charge/cover proofs.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json,time
from pathlib import Path
import networkx as nx
import fresh_cross_vertex_support_migration_20261018 as M
import fresh_cross_vertex_current290_local_covers_20261018 as L
import fresh_290_class_shortlist_pair_geometry_20261017 as T

ROOT=M.ROOT;H=T.H;W=T.W;N=T.N
TAG='fresh_cross_vertex_historical32_local_covers_20261018'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.jsonl')


def normal(x):return json.loads(json.dumps(x,sort_keys=True))


def prepare():
    migration=M.read(M.TAG);base=M.read(L.TAG);p=N.prepare()
    order=M.read(M.H)['source_precedence']+list(M.EXTRA)
    added={r['face']:r for r in migration['newly_positive_historical_outside290']}
    selected={};presentations={}
    for tag in order:
        for r in W.read(tag)['selected_faces']:
            fi=r['face']
            if fi not in added:continue
            q=added[fi]
            assert (r['source_context'],r['canonical_prism_core'],r['provenance'],r['raw_family'])==(q['source_context'],q['canonical_core'],q['provenance'],q['raw_family'])
            if fi in selected:assert selected[fi]['corners']==r['corners']
            else:selected[fi]=r
            presentations.setdefault(fi,[]).append(tag)
    assert set(selected)==set(added) and len(selected)==32
    affected=sorted({c['vertex'] for r in selected.values() for c in r['corners']})
    assert len(affected)==161 and sum(len(r['corners']) for r in selected.values())==248
    latest=M.current_domains();assert len(set(affected)-set(latest))==29
    changed={r['vertex']:r for r in W.read(N.TAG)['forest_records']}
    patches={r['vertex']:r for tag in (M.SHORT,M.P2599) for r in W.read(tag)['forest_records']}
    sources=M.bindings(tuple(dict.fromkeys((M.MODEL,M.TAG,L.TAG,*M.DOMAIN_STEMS,*order))))
    for mod in (M,L,T,T.U,N,N.O,N.O.D,H,W):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=M.digest(path)
    sources[f'scripts/{TAG}.py']=M.digest(f'scripts/{TAG}.py')
    return dict(p=p,order=order,selected=selected,presentations=presentations,affected=affected,
                latest=latest,changed=changed,patches=patches,base=base,sources=sources,
                alpha={i:F(a) for i,a in base['face_multipliers']},
                mu={v:F(a) for v,a in base['vertex_multipliers']})


def assemble(v,p):
    """Transport old graph exactly; conservative first-owner merge for new nodes."""
    if v in p['latest']:
        old,g,owner,unowned,demands,nodes,groups,forest,refs,source=T.current(v,p['p'],p['changed'])
        nodes=normal(nodes);groups=normal(groups)
        if v in p['patches']:
            patch=p['patches'][v]
            assert patch['literal_stars']==nodes and patch['corner_orbits']==groups
            for i,j,code in patch['newly_owned_pair_domain']:
                assert g.has_edge(i,j)
                if code is not None:g.remove_edge(i,j)
        assert groups==p['latest'][v]['groups']
        baseline=p['latest'][v]['row'];oldvectors=baseline['maximal_count_vectors']
    else:
        g=nx.Graph();owner={};nodes=[];groups=[];oldvectors=[]
        forest=next(W.forests(tag)[v]['forest'] for tag in p['order'] if v in W.forests(tag))
    oldn=len(nodes);oldgroups=len(groups);oldedges={tuple(sorted(e)) for e in g.edges()}
    keys={H.frozen(s['contacts']):i for i,s in enumerate(nodes)};assert len(keys)==oldn
    orbitkeys={frozenset(z['literal_stars']):j for j,z in enumerate(groups)}
    origins={};newrefs=[]
    rows=[(tag,W.forests(tag)[v]) for tag in p['order'] if v in W.forests(tag)]
    for tag,row in rows:
        assert row['forest']==forest
        for oldg in row['corner_orbits']:
            ff={f for f,m in oldg['positive_face_coefficients']};assert len(ff)==1
            fi=next(iter(ff))
            if fi not in p['selected']:continue
            oo={H.origin_key(c):c for c in oldg['origins']}
            expected={H.origin_key(c):c for c in p['selected'][fi]['corners'] if c['vertex']==v}
            assert len(oo)==len(oldg['origins']) and set(oo)<=set(expected)
            assert all(c==expected[k] for k,c in oo.items())
            assert oldg['positive_face_coefficients']==[[fi,len(oo)]]
            indices=[]
            for oi in oldg['literal_stars']:
                star=row['literal_stars'][oi];assert star['index']==oi
                key=H.frozen(star['contacts'])
                if key in keys:
                    i=keys[key];assert all(star[k]==nodes[i][k] for k in H.NATIVE_FIELDS)
                else:
                    i=len(nodes);keys[key]=i;nodes.append(normal(star));g.add_node(i);newrefs.append([tag,oi])
                assert i>=oldn,'a new historical class coincides with an old complete orbit'
                indices.append(i)
            assert len(indices)==len(set(indices)) and indices
            assert all(c['literal_orientations']==len(indices) for c in oo.values())
            ok=frozenset(indices)
            if ok in orbitkeys:
                gi=orbitkeys[ok];assert groups[gi]['positive_face_coefficients']==[[fi,len(oo)]]
                assert {H.origin_key(c):c for c in groups[gi]['origins']}==oo
            else:
                assert not any(i in owner for i in indices),'partial native orbit overlap'
                gi=len(groups);orbitkeys[ok]=gi
                for i in indices:owner[i]=gi
                groups.append(dict(index=gi,face=fi,literal_stars=sorted(indices),
                    origins=[oo[k] for k in sorted(oo,key=str)],positive_face_coefficients=[[fi,len(oo)]]))
                for k,c in oo.items():assert k not in origins;origins[k]=c
    expected={H.origin_key(c):c for r in p['selected'].values() for c in r['corners'] if c['vertex']==v}
    assert origins==expected and set(owner)==set(range(len(nodes)))
    domains=[]
    for tag,row in rows:
        mp={}
        for star in row['literal_stars']:
            key=H.frozen(star['contacts'])
            if key in keys:
                i=keys[key];assert i not in mp
                assert all(star[k]==nodes[i][k] for k in H.NATIVE_FIELDS)
                mp[i]=star['index']
        if mp:domains.append((tag,mp,W.relation(tag,v)))
    owned=Counter();excluded=Counter();newbad=[];unowned=0
    for i,j in combinations(range(len(nodes)),2):
        if j<oldn:continue  # even old unowned/allowed pairs are kept unchanged.
        chosen=None;bad=False
        for tag,mp,relation in domains:
            if i in mp and j in mp:
                chosen=tag;bad=tuple(sorted((mp[i],mp[j]))) in relation;break
        if chosen is None:unowned+=1
        else:owned[chosen]+=1
        if bad:newbad.append([i,j]);excluded[chosen]+=1
        else:g.add_edge(i,j)
    assert {tuple(sorted(e)) for e in g.edges() if max(e)<oldn}==oldedges
    assert sum(owned.values())+unowned==(len(nodes)*(len(nodes)-1)-oldn*(oldn-1))//2
    W.relation.cache_clear()
    demand=[sum((p['alpha'].get(f,F())*F(m) for f,m in z['positive_face_coefficients']),F()) for z in groups]
    u=dict(vertex=v,forest=forest,old_literal_stars=oldn,old_corner_orbits=oldgroups,
        new_literal_references=newrefs,corner_orbits=groups,
        native_stars_sha256=L.jhash([{k:s[k] for k in H.NATIVE_FIELDS} for s in nodes]),
        old_edges_sha256=L.jhash(sorted(oldedges)),new_excluded_pairs=newbad,
        source_domains=[dict(source=tag,new_to_old_star_indices=sorted(mp.items()),first_owned_pairs=owned[tag],first_owned_exclusions=excluded[tag]) for tag,mp,rel in domains],
        counts=dict(literal_stars=len(nodes),new_literal_stars=len(nodes)-oldn,corner_orbits=len(groups),
            new_corner_orbits=len(groups)-oldgroups,new_actual_corners=len(origins),
            new_pairs=sum(owned.values())+unowned,new_exclusions=len(newbad),
            owned_unknown_pairs=sum(owned.values())-len(newbad),unowned_allowed_pairs=unowned),
        demands=list(map(str,demand)),baseline_count_vectors_sha256=L.jhash(oldvectors))
    return u,g,owner,demand


def run(saved=None,inventory=False):
    start=time.monotonic();p=prepare();key=L.jhash([p['sources'],p['base']['face_multipliers']]);cached={}
    if saved:
        assert saved['sources']==p['sources'];cached={r['vertex']:r for r in saved['forest_records']}
    elif PROGRESS.exists() and not inventory:
        for line in PROGRESS.read_text().splitlines():
            r=json.loads(line);assert r['key']==key and r['row']['vertex'] not in cached;cached[r['row']['vertex']]=r['row']
    records=[];counts=Counter();newmu=dict(p['mu']);columns=[];opened=[]
    for step,v in enumerate(p['affected'],1):
        u,g,owner,demand=assemble(v,p);counts.update(u['counts'])
        if inventory:
            print('INVENTORY',step,v,u['counts'],flush=True);continue
        old=cached.get(v)
        if old and not saved:
            assert old['domain']==normal(u);row=old
        else:
            enum,error=H.vectors_for(g,owner,len(demand))
            if error:
                assert saved is None;opened.append(dict(vertex=v,domain=u,reason=error));continue
            vectors,reps,ncliques,npatterns=enum;vectors=list(map(list,vectors));reps=list(map(list,reps))
            for vec,cl in zip(vectors,reps):
                assert all(g.has_edge(i,j) for i,j in combinations(cl,2))
                assert vec==[Counter(owner[i] for i in cl).get(j,0) for j in range(len(demand))]
            c=old['certificate'] if old else H.exact_profile(vectors,demand)
            value=W.check_certificate(vectors,demand,c)
            row=dict(vertex=v,domain=u,maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,
                maximal_count_vectors=vectors,literal_representatives=reps,certificate=c,
                baseline_vertex_mass=str(p['mu'].get(v,F())),exact_local_optimum=str(value),
                required_increase=str(max(F(),value-p['mu'].get(v,F()))))
            row=normal(row)
            if old:assert row==old
            elif not saved:
                with PROGRESS.open('a') as f:f.write(json.dumps(dict(key=key,row=row),separators=(',',':'))+'\n')
        value=W.check_certificate(row['maximal_count_vectors'],demand,row['certificate']);newmu[v]=max(newmu.get(v,F()),value)
        if value>p['mu'].get(v,F()):
            terms=Counter()
            for w,z in zip(map(F,row['certificate']['weights']),u['corner_orbits']):
                for f,m in z['positive_face_coefficients']:terms[f]+=w*F(m)
            coeffs=[[f,str(x)] for f,x in sorted(terms.items()) if x]
            assert sum((p['alpha'].get(f,F())*F(x) for f,x in coeffs),F())==value
            columns.append(dict(vertex=v,vertex_cost='1',positive_face_coefficients=coeffs,
                orbit_weights=row['certificate']['weights'],exact_load=str(value)))
        records.append(row)
        if step%10==0 or value>p['mu'].get(v,F()):print('MERGED32',step,'/161',v,'columns',len(columns),'seconds',round(time.monotonic()-start,2),flush=True)
    assert counts['new_actual_corners']==248
    if inventory:print('INVENTORY_COMPLETE',dict(counts),flush=True);return
    assert not opened,opened
    assert len(records)==161 and len(p['latest'])-len(set(p['affected'])&set(p['latest']))==439
    inherited=[dict(vertex=r['vertex'],record_sha256=L.jhash(r)) for r in p['base']['forest_records'] if r['vertex'] not in p['affected']]
    assert len(inherited)==439
    multi=[]
    for col in M.read(M.MODEL)['prepared']['new_columns']:
        load=sum((p['alpha'].get(i,F())*m for i,m in col['integer_face_coefficients']),F())
        before=sum((p['mu'].get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        after=sum((newmu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        assert load<=before<=after
        if len(col['integer_vertex_capacities'])>1:multi.append(dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after)))
    assert len(multi)==14 and all(newmu.get(v,F())>=x for v,x in p['mu'].items())
    num=F(p['base']['objective_numerator']);den=sum(newmu.values(),F())
    out=dict(status='PASS_NEW5121_HISTORICAL32_MERGED_EXACT_LOCAL_COVERS',sources=p['sources'],
        face_multipliers=p['base']['face_multipliers'],initial_vertex_multipliers=p['base']['vertex_multipliers'],
        vertex_multipliers=[[v,str(x)] for v,x in sorted(newmu.items()) if x],
        objective_numerator=str(num),normalizer=str(den),bound=str(num/den),baseline_bound=p['base']['bound'],
        added_faces=[dict(selected_face=r,historical_sources=p['presentations'][f]) for f,r in sorted(p['selected'].items())],
        owner_precedence=p['order'],forest_records=records,inherited_forests=inherited,
        columns=columns,preserved_multi_owner_columns=multi,counts=dict(counts),
        represented_classes=322,positive_represented_classes=259,represented_forests=600,
        changed_forests=132,new_forests=29,inherited_forest_count=439,new_charge_columns=len(columns),
        complete_count_vectors=sum(len(r['maximal_count_vectors']) for r in records),
        complete_maximal_cliques=sum(r['maximal_literal_cliques'] for r in records),
        generation_local_LPs=161,replay_local_LPs=0,global_optimizer_invocations=0,new_physical_pair_tests=0,
        scope='Complete conservative merged322-class/600-forest graph cones at fixed5121 alpha. Every old pair outcome, including old unowned allowed pairs, is preserved; new pairs use their first full historical owner, with unknowns allowed. Native compatibility facts are frozen inputs; new complete clique domains and exact local charge/cover certificates are checked. All14 two-owner costs retained. Missing classes, merged Q provenance, simplex/joint coverage and physical realization remain open. No global optimum, target12 primal feasibility or dimension theorem.')
    print(out['status'],'newcolumns',len(columns),'bound',float(num/den),'seconds',round(time.monotonic()-start,2),flush=True)
    return normal(out)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--inventory',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    out=run(saved,args.inventory)
    if args.inventory:pass
    elif args.verify:assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,separators=(',',':'))+'\n')
