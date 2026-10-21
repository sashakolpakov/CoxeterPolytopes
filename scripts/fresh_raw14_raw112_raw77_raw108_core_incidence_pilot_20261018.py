#!/usr/bin/env python3
"""Bounded exact core/face preflight for four uncovered arithmetic families.

Complete selected core supports and whole spectator inventory; no new pair
census, local cover, column export or global optimization.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_raw56_raw57_extra_vertex_family_20261018 as T
D=T.D;A=T.A;M=T.M;R=T.R;V=T.V;K=T.K;P=T.P;ROOT=T.ROOT
TAG='fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018'
PARENT='fresh_raw56_raw57_raw76_joint_marker_family_20261018'
SUPPORT='fresh_raw76_complete_family_support_20261018'
CORES={14:12,112:92,77:66,108:88}
normal=D.normal


def scaled(g):
    n=len(g);sc=[None]*n
    for root in range(n):
        if sc[root] is not None:continue
        sc[root]=K(1);todo=[root]
        while todo:
            i=todo.pop()
            for j in range(n):
                if i==j or g[i][j]==A.RIGHT or sc[j] is not None:continue
                sc[j]=sc[i]*V.squared(g[i][j]);todo.append(j)
    G=[[2*sc[i]if i==j else K()for j in range(n)]for i in range(n)]
    for i,j in combinations(range(n),2):
        if g[i][j]==A.RIGHT:continue
        z=A.M.ksqrt(V.squared(g[i][j])*sc[i]*sc[j]);assert z is not None and z.sign()>0
        G[i][j]=G[j][i]=-z
    return G,sc


def active_graph(g,ids):
    return {i:{j:g[i][j][0]for j in ids if i!=j and g[i][j]!=A.RIGHT}for i in ids}


def forest_type(g,ids):
    names={R.C.code(x):n for n,x in R.C.models()}
    return sorted(names[R.C.code(x)]for x in R.C.pieces(active_graph(g,ids)))


def core_geometry(raw,canonical,ri):
    c=canonical[CORES[ri]];g=A.V.raw_graph(raw['records'][ri]);n=len(g);rank=n-2;G,sc=scaled(g)
    assert c['source_core']==ri and c['nodes']==n
    assert P.congruence(G)['inertia']==[rank,1,1]
    assert P.congruence([[K(x.a,-x.b)for x in row]for row in G])['inertia']==[rank+1,0,1]
    records=[];pd=set()
    for om in combinations(range(n),2):
        ids=sorted(set(range(n))-set(om));ii=P.congruence(P.sub(G,ids))['inertia']
        row=dict(omitted=list(om),support=ids,inertia=ii)
        if ii==[rank,0,0]:
            pd.add(om);row['active_forest']=forest_type(g,ids)
        records.append(row)
    assert sorted(map(list,pd))==c['pd_deletions']
    occurrences=c['aggregate_corner_occurrences'];used={tuple(sorted(x['vertex_omitted']))for x in occurrences}
    assert len(occurrences)==4*c['required']and used<=pd
    extras=sorted(pd-used);assert extras==([(2,4),(2,5)]if ri in(112,108)else[])
    autos=V.automorphisms(g);assert autos==[list(range(n))]
    # Every maximal star in the omission graph defines a candidate simplex.
    simplexes=[]
    for e in range(n):
        sides=sorted(next(i for i in om if i!=e)for om in pd if e in om)
        if len(sides)<3:continue
        fixed=sorted(set(range(n))-set(sides)-{e})
        assert P.congruence(P.sub(G,fixed))['inertia']==[len(fixed),0,0]
        sg=P.schur(G,fixed,sides);inv,_=V.inverse(sg)
        assert P.congruence(sg)['inertia']==[len(sides)-1,1,0]
        assert all(x.sign()<0 for row in inv for x in row)
        supports=[sorted(fixed+list(set(sides)-{j}))for j in sides]
        assert all(tuple(sorted(set(range(n))-set(s)))in pd for s in supports)
        simplexes.append(dict(excluded_wall=e,fixed=fixed,sides=sides,dimension=len(sides)-1,
            vertex_supports=supports,residual_gram=P.out(sg),inverse=P.out(inv)))
    extra=[]
    for om in extras:
        ids=sorted(set(range(n))-set(om));J,_=V.inverse(P.sub(G,ids));cvs=[[G[i][j]for j in ids]for i in om]
        pr=[[sum((cvs[i][a]*J[a][b]*cvs[j][b]for a in range(rank)for b in range(rank)),K())for j in range(2)]for i in range(2)]
        hs=[pr[i][i]-G[om[i]][om[i]]for i in range(2)];res=pr[0][1]-G[om[0]][om[1]]
        assert all(h.sign()>0 for h in hs)and res.sign()>0 and res*res==hs[0]*hs[1]
        triangles=[sorted(s)for sx in simplexes if set(ids)==set(sx['fixed'])|set(sx['sides'])-{om[1]}
                   for s in combinations(sorted(set(sx['sides'])&set(ids)),2)]
        extra.append(dict(omitted=list(om),support=ids,active_forest=forest_type(g,ids),
            scaled_covectors=[[x.out()for x in cv]for cv in cvs],scaled_projection=P.out(pr),
            scaled_heights=[h.out()for h in hs],positive_height_product_root=res.out(),
            forced_triangle_port_pairs=sorted(set(map(tuple,triangles)))))
    return g,dict(raw_family=ri,canonical_core=CORES[ri],nodes=n,spectator_rank=13-rank,
        native_gram=g,squared_scales=[x.out()for x in sc],scaled_gram=P.out(G),
        inertia=[rank,1,1],conjugate_inertia=[rank+1,0,1],support_records=records,
        required_quads=c['required'],quad_partitions=c['quad_partitions'],counted_vertex_omissions=sorted(used),
        ordinary_corner_occurrences=occurrences,extra_vertices=extra,compact_simplex_faces=simplexes,
        native_automorphisms=autos)


def expected_maps(active_types,target_forest):
    ac=Counter(active_types);tc=Counter(target_forest);models=dict(R.C.models());out=1
    for name,k in ac.items():
        assert tc[name]>=k
        # Compute label-preserving finite-component maps directly.
        import networkx as nx
        gg=nx.Graph();gg.add_nodes_from(models[name]);gg.add_edges_from((i,j,{'label':m})for i,row in models[name].items()for j,m in row.items()if i<j)
        aut=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(gg,gg,edge_match=lambda a,b:a['label']==b['label']).isomorphisms_iter())
        for j in range(k):out*=tc[name]-j
        out*=aut**k
    return out


def inventory(inp,native,geom,positive):
    forests=inp.read(M.FORESTS);spectators=inp.read(A.B.Q.SPECTATORS);quads=inp.read(M.QUAD)
    lookup={tuple(r['forest']):r['index']for r in forests['vertices']};models=dict(R.C.models());rows=[];markers=[]
    for fi,ff in enumerate(forests['faces']):
        if ff['kind']!='canonical_dotted_quad_aggregate':continue
        relevant=[(ctx,ci)for ctx,ci in ff['source_contexts']if ci in CORES.values()]
        if not relevant:continue
        assert len(relevant)==len(ff['source_contexts'])==1
        ctx,ci=relevant[0];ri=next(r for r,c in CORES.items()if c==ci);sp=spectators['contexts'][ctx]['spectator_forest']
        assert sum(len(models[n])for n in sp)==geom[ri]['spectator_rank']
        assert quads['source_core_maps'][spectators['contexts'][ctx]['core']]['canonical_core']==ci
        mass=F(positive.get(fi,{}).get('mass','0'))
        if mass:
            pp=positive[fi];assert pp['raw_family']==ri and pp['source_context']==ctx and pp['required']==geom[ri]['required_quads']
        rows.append(dict(face=fi,raw_family=ri,source_context=ctx,canonical_core=ci,spectator_forest=sp,
            required_quads=geom[ri]['required_quads'],current_mass=str(mass),extra_vertex_occurrences=len(geom[ri]['extra_vertices'])))
        for ex in geom[ri]['extra_vertices']:
            vf=sorted(sp+ex['active_forest']);v=lookup[tuple(vf)];m=expected_maps(ex['active_forest'],vf)
            markers.append(dict(face=fi,raw_family=ri,source_context=ctx,omitted=ex['omitted'],
                spectator_forest=sp,active_forest=ex['active_forest'],target_vertex=v,target_forest=vf,
                active_role_orientations=m,actual_extra_occurrences=1,current_mass=str(mass)))
    # One highest-positive class per family: complete actual active embeddings
    # and native recovery, with orientation keys separate from occurrence mass.
    pilots=[]
    for ri in(112,108):
        chosen=min((r for r in rows if r['raw_family']==ri and F(r['current_mass'])>0),key=lambda r:(-F(r['current_mass']),r['face']))
        for row in[r for r in markers if r['face']==chosen['face']]:
            ex=next(r for r in geom[ri]['extra_vertices']if r['omitted']==row['omitted']);g=native[ri];ids=ex['support']
            target=R.P.forest_graph(row['target_forest'],models);J,_,_=R.projection_data(target)
            maps=R.active_embeddings(active_graph(g,ids),target);assert len(maps)==row['active_role_orientations']
            saved=[];seen=set()
            for mp in maps:
                cov=[tuple(sorted((mp[j],tuple(g[i][j]))for j in ids if g[i][j]!=A.RIGHT))for i in row['omitted']]
                assert len(set(cov))==2;key=tuple(sorted(cov));assert key not in seen;seen.add(key)
                field=A.pair_field(tuple(cov),(g[row['omitted'][0]][row['omitted'][1]],));hs=[]
                for c in cov:
                    z=R.projection(c,c,J,field);assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;hs.append(h)
                z=R.projection(cov[0],cov[1],J,field);prescribed=field.negative_contact(g[row['omitted'][0]][row['omitted'][1]])
                assert field.signs_difference(tuple(x-y for x,y in zip(z,prescribed)),hs[0]*hs[1])==0
                saved.append(dict(active_embedding=sorted(mp.items()),full_frame_covectors=cov,
                    positive_heights=[h.out()for h in hs],mutual_projection=R.qout(z)))
            pilots.append(dict(**row,finite_frame=target,complete_active_maps=saved))
    return rows,markers,pilots


def workload(inp,markers):
    # Saved-domain reconstruction only: no native or pair census.
    domains,_=D.L.current_domains()
    ordinary={v:sum(len(g['literal_stars'])for g in d['groups'])for v,d in domains.items()}
    for tag in(D.B.TAG,D.THIRD):
        for r in inp.read(tag)['forest_records']:ordinary[r['vertex']]=len(r['domain']['literal_stars'])
    assert len(ordinary)==600
    previous=Counter()
    for r in inp.read(T.U.TAG)['forest_records']:previous[r['vertex']]+=r['counts']['markers']
    for r in inp.read(T.TAG)['forest_records']:previous[r['vertex']]+=r['counts']['markers']
    newer=Counter()
    for r in markers:newer[r['target_vertex']]+=r['active_role_orientations']
    records=[dict(vertex=v,new_markers=n,old_ordinary_stars=ordinary.get(v,0),
        old_family_markers=previous[v],new_ordinary_pairs=n*ordinary.get(v,0),
        new_old_marker_pairs=n*previous[v],new_marker_pairs=n*(n-1)//2)for v,n in sorted(newer.items())]
    baseline=inp.read(D.L.TAG)
    for p,h in baseline['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return dict(records=records,target_forests=len(records),ordinary_target_forests=sum(r['old_ordinary_stars']>0 for r in records),
        old_marker_target_forests=sum(r['old_family_markers']>0 for r in records),
        new_ordinary_pairs=sum(r['new_ordinary_pairs']for r in records),
        new_old_marker_pairs=sum(r['new_old_marker_pairs']for r in records),
        new_marker_pairs=sum(r['new_marker_pairs']for r in records),
        scope='Unreduced pair-position inventory only. No compatibility test; symmetry reductions and retained old outcomes remain prerequisites for a future extension.')


def run():
    inp=D.Q.Inputs();raw=inp.read(A.M.TAGS['raw']);quads=inp.read(M.QUAD);support=inp.read(SUPPORT);parent=inp.read(PARENT)
    assert parent['remaining_support_diagnostic']['highest_saved_demand_families'][:4]==[next(r for r in support['ranked_unambiguous_missing_families']if r['raw_family']==i)for i in CORES]
    native={};geom={}
    for ri in CORES:native[ri],geom[ri]=core_geometry(raw,quads['canonical_cores'],ri)
    positive={r['face']:r for r in support['positive_faces']};rows,markers,pilots=inventory(inp,native,geom,positive)
    work=workload(inp,markers)
    summaries=[]
    for ri in CORES:
        rr=[r for r in rows if r['raw_family']==ri];mm=[r for r in markers if r['raw_family']==ri]
        summaries.append(dict(raw_family=ri,complete_classes=len(rr),positive_classes=sum(F(r['current_mass'])>0 for r in rr),
            positive_mass=str(sum((F(r['current_mass'])for r in rr),F())),required_quads=geom[ri]['required_quads'],
            actual_extra_vertices=len(geom[ri]['extra_vertices']),marker_domains=len(mm),
            role_orientation_count=sum(r['active_role_orientations']for r in mm),
            positive_extra_occurrence_mass=str(sum((F(r['current_mass'])for r in mm),F()))))
    paths=set(inp.paths)|{str(Path(__file__).relative_to(ROOT))}
    for mod in(T,T.U,D,D.L,R,V,P,A,A.M,R.C):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(f'scripts/{SUPPORT}.py')
    for tag in(PARENT,'fresh_raw76_complete_family_full_model_20261018','fresh_quad_product_aggregation_20261017','fresh_forest1325_cross_vertex_incidence_20261017'):
        paths.add(f'results/{tag}.md')
    sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)}
    result=dict(status='PASS_BOUNDED_FOUR_RAW_CORE_INCIDENCE_PILOT',sources=sources,geometries=list(geom.values()),
        complete_family_inventory=rows,extra_marker_inventory=markers,native_recovery_pilots=pilots,summaries=summaries,
        future_pair_workload=work,
        counts=dict(core_support_tests=sum(len(g['support_records'])for g in geom.values()),
            exact_compact_simplex_faces=sum(len(g['compact_simplex_faces'])for g in geom.values()),
            family_classes=len(rows),extra_marker_domains=len(markers),
            extra_role_orientations=sum(r['active_role_orientations']for r in markers),
            exact_native_pilot_maps=sum(len(r['complete_active_maps'])for r in pilots)),
        current_joint_bound=parent['exported']['repair']['bound'],
        normalizer_gap_to12=parent['remaining_support_diagnostic']['added_normalizer_needed_to_reach12'],
        new_pair_tests=0,exported_columns=[],old_outcomes_changed=False,optimizer_invocations=0,
        scope='Exact four-core support and compact-simplex diagnostic, complete selected spectator inventory, and four bounded full native recovery domains. No marker capacity, full-family pair compatibility, new charge, cover, ratio improvement or ambient realization is claimed.')
    return normal(result)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(out.read_text())==result
    else:assert not out.exists();out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['counts'],flush=True)
    for r in result['summaries']:print(r,flush=True)
