#!/usr/bin/env python3
"""Complete raw76 spectator family by exact relative-role and orbit transport.

No old Q-Q outcome is changed, and no global optimizer is called. The frozen
three-domain pilot is reused. Only one marker's new ordinary cross pairs are
tested in each remaining represented forest.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_raw76_forced_extra_vertex_family_20261018 as D

ROOT=D.ROOT;M=D.M;H=D.H;R=D.R;A=D.A;L=D.L;Z=D.Z;B=D.B;W=D.W
TAG='fresh_raw76_extra_vertex_complete_family_20261018'
normal=D.normal


def role(marker):
    p=dict(marker['active_embedding']);return p[0],p[4],tuple(p[i]for i in(1,3,6))


def relative(a,b):
    x,y,h=role(a);u,v,k=role(b);names={x:0,y:1}
    for i in(u,v):
        if i not in names:names[i]=len(names)
    return names[u],names[v],int(h!=k)


def marker_rule(a,b):
    x,y,h=role(a);u,v,k=role(b)
    return x==u and y!=v and h!=k


def marker_theorem(family,g,active):
    host=next(r for r in family if Counter(r['target_forest'])=={'A1':7,'H3':2})
    target,J,markers=D.marker_maps(host,active,g);anchor=markers[0]
    reps={}
    for i,m in enumerate(markers[1:],1):reps.setdefault(relative(anchor,m),i)
    assert len(reps)==13
    out=[]
    for pattern,i in sorted(reps.items()):
        graph=nx.Graph();dec,expr=D.physical(target,J,[markers[i]['native']],[anchor],graph)
        assert len(dec)==1
        allowed=dec[0]['exclusion']is None
        assert allowed==marker_rule(anchor,markers[i])
        out.append(dict(relative=list(pattern),other_marker=i,allowed=allowed,
            decision=dec[0],wall_expressions=expr))
    return dict(host=host,finite_frame=target,markers=markers,anchor=0,relative_patterns=out,
        capacity_formula='min(number_of_A1_components-1,number_of_H3_components)',
        interpretation='Exact conservative pair-gate capacity, not simultaneous physical realization.')


def complete_current(v,row,p,domains,local,patches):
    if v not in domains:
        return [],[],nx.Graph(),{},row['target_forest'],[]
    if v in p['rows32']and p['rows32'][v]['domain']['old_literal_stars']==0:
        r=p['rows32'][v];u=r['domain'];assert u['old_corner_orbits']==0
        nodes=[normal(W.forests(tag)[v]['literal_stars'][i])for tag,i in u['new_literal_references']]
        assert L.jhash([{k:s[k]for k in H.NATIVE_FIELDS}for s in nodes])==u['native_stars_sha256']
        groups=normal(u['corner_orbits']);graph=nx.complete_graph(len(nodes));graph.remove_edges_from(u['new_excluded_pairs'])
        assert groups==normal(domains[v]['groups'])and v not in patches
        assert not any(rr['vertex']==v for tag in(B.TAG,D.THIRD)for rr in M.read(tag)['forest_records'])
        owner={i:z['index']for z in groups for i in z['literal_stars']}
        assert set(owner)==set(range(len(nodes)))
        for cl in r['literal_representatives']:assert all(graph.has_edge(i,j)for i,j in combinations(cl,2))
        return nodes,groups,graph,owner,u['forest'],[Z.U.TAG]
    return D.current(v,p,local,domains,patches)


def graph_components(target):
    graph=nx.Graph();graph.add_nodes_from(target)
    for i,row in target.items():
        for j,m in row.items():graph.add_edge(i,j,label=m)
    return graph,[tuple(sorted(c))for c in nx.connected_components(graph)]


def automorphism(target,base,other):
    graph,components=graph_components(target)
    p=dict(base['active_embedding']);q=dict(other['active_embedding']);perm={p[i]:q[i]for i in p}
    old_used=set(perm);new_used=set(perm.values())
    left=[c for c in components if not set(c)&old_used];right=[c for c in components if not set(c)&new_used]
    assert sum(map(len,left))==sum(map(len,right))==8
    for c in left:
        for k,e in enumerate(right):
            if len(c)!=len(e):continue
            iso=nx.algorithms.isomorphism.GraphMatcher(graph.subgraph(c),graph.subgraph(e),edge_match=lambda a,b:a['label']==b['label'])
            if iso.is_isomorphic():perm.update(next(iso.isomorphisms_iter()));right.pop(k);break
        else:raise AssertionError('unmatched passive component')
    assert not right and set(perm)==set(perm.values())==set(range(13))
    assert all(target[i].get(j,2)==target[perm[i]].get(perm[j],2)for i,j in combinations(range(13),2))
    return perm


def transport_native(star,perm):
    co=[sorted([[perm[i],e]for i,e in c])for c in star['contacts']]
    order=sorted(range(2),key=lambda i:co[i])
    return dict(contacts=[co[i]for i in order],
        occupied_pairs=sorted(sorted(perm[i]for i in pair)for pair in star['occupied_pairs']),
        squared_timelike_heights=[star['squared_timelike_heights'][i]for i in order],
        native_omitted_product=star['native_omitted_product'])


def process(row,active,g,p,domains,local,patches,frozen):
    v=row['target_vertex'];target,J,markers=D.marker_maps(row,active,g)
    nodes,groups,oldgraph,owner,forest,history=complete_current(v,row,p,domains,local,patches)
    assert forest==row['target_forest'];n=len(nodes);m=len(markers)
    active_coordinates={i for i,j in J if i==j}
    assert all(i in active_coordinates for s in nodes for cov in s['contacts']for i,e in cov)
    canonical_nodes=[transport_native(s,dict(zip(range(13),range(13))))for s in nodes]
    oldkeys={H.frozen(s['contacts']):i for i,s in enumerate(canonical_nodes)}
    markkeys={H.frozen(s['native']['contacts']):i for i,s in enumerate(markers)}
    assert len(oldkeys)==n and len(markkeys)==m and not(set(oldkeys)&set(markkeys))
    transports=[]
    for marker in markers:
        perm=automorphism(target,markers[0],marker);on=[];mn=[]
        for i,s in enumerate(nodes):
            t=transport_native(s,perm);j=oldkeys[H.frozen(t['contacts'])]
            assert all(normal(t[k])==normal(canonical_nodes[j][k])for k in H.NATIVE_FIELDS)
            assert owner[i]==owner[j]  # preserves every actual origin/coefficient.
            on.append(j)
        for s in markers:
            t=transport_native(s['native'],perm);j=markkeys[H.frozen(t['contacts'])]
            assert all(normal(t[k])==normal(markers[j]['native'][k])for k in H.NATIVE_FIELDS)
            mn.append(j)
        assert len(set(on))==n and len(set(mn))==m and mn[0]==marker['index']
        transports.append(dict(frame_permutation=sorted(perm.items()),ordinary_permutation=on,marker_permutation=mn))
    if v in frozen:
        prior=frozen[v];assert prior['old_nodes']==nodes and prior['old_groups']==groups and prior['markers']==normal(M.enc(markers))
        assert prior['old_graph_edges']==sorted(tuple(sorted(e))for e in oldgraph.edges())or prior['old_graph_edges']==normal(sorted(tuple(sorted(e))for e in oldgraph.edges()))
        dec=[x for x in prior['new_pair_decisions']if x['pair'][1]==n]
        assert len(dec)==n;expr=prior['wall_expressions'];reuse=True
    else:
        dec,expr=D.physical(target,J,nodes,[markers[0]],nx.Graph());assert len(dec)==n;reuse=False
    anchor={x['pair'][0]for x in dec if x['exclusion']is None}
    assert [x['pair']for x in dec]==[[i,n]for i in range(n)]
    cross=[]
    for k,t in enumerate(transports):
        cross.extend([i,k]for i in sorted(t['ordinary_permutation'][i]for i in anchor))
    medges=[list(pair)for pair in combinations(range(m),2)if marker_rule(*(markers[i]for i in pair))]
    if v in frozen:
        assert cross==sorted([i,j-n]for i,j in frozen[v]['graph_edges']if i<n<=j)
        assert medges==[[i-n,j-n]for i,j in frozen[v]['graph_edges']if i>=n]
    c=Counter(forest);capacity=min(c['A1']-1,c['H3']);assert capacity>=1
    # Explicit abstract attaining marker clique; upper bound is the proved
    # injectivity into distinct second A1 and H3 roles at one first A1 role.
    mg=nx.Graph();mg.add_nodes_from(range(m));mg.add_edges_from(medges)
    mc=max(nx.find_cliques(mg),key=len);assert len(mc)==capacity
    return normal(M.enc(dict(face=row['face'],vertex=v,forest=forest,finite_frame=target,old_nodes=nodes,
        old_groups=groups,old_graph_edges=sorted(tuple(sorted(e))for e in oldgraph.edges()),prior_history=history,
        markers=markers,transports=transports,anchor_decisions=dec,wall_expressions=expr,
        reused_frozen_pilot=reuse,marker_ordinary_edges=cross,marker_edges=medges,
        marker_capacity=capacity,attaining_marker_clique=sorted(mc),
        counts=dict(ordinary_stars=n,markers=m,cross_pairs=n*m,marker_pairs=m*(m-1)//2,
            direct_anchor_tests=0 if reuse else n,reused_anchor_tests=n if reuse else 0,
            retained_cross_pairs=len(cross),retained_marker_pairs=len(medges)))))


def export_charges(rows,domains,model):
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};mu={v:F(x)for v,x in dual['vertex_multipliers']}
    certificates=defaultdict(list);vectors={v:d['domain']['maximal_count_vectors']for v,d in domains.items()}
    for tag in(L.TAG,B.TAG,D.THIRD):
        for r in M.read(tag)['forest_records']:
            certificates[r['vertex']].append((tag,r['certificate']['weights']))
            if tag!=L.TAG:vectors[r['vertex']]=r['maximal_count_vectors']
    columns=[];newmu=dict(mu);beforeviol=0;markerinduced=0;quotas=[]
    for r in rows:
        v=r['vertex'];groups=r['old_groups'];owner={i:z['index']for z in groups for i in z['literal_stars']}
        omitted=sorted({owner[i]for i,k in r['marker_ordinary_edges']})
        assert not omitted or(v==1304 and omitted==[0]and groups[0]['positive_face_coefficients']==[[60808,1]])
        cap=r['marker_capacity'];candidates=certificates[v]if groups else [('marker_only',[])]
        assert candidates
        unique={};oldvalues=[]
        for tag,weights in candidates:
            original=list(map(F,weights))+[F()]*(len(groups)-len(weights));assert len(original)==len(groups)and all(w>=0 for w in original)
            vv=vectors[v]if groups else [[]]
            assert all(len(x)==len(original)for x in vv)
            assert all(sum((w*n for w,n in zip(original,x)),F())<=1 for x in vv)
            w=[F()if i in omitted else z for i,z in enumerate(original)]
            assert all(w[owner[i]]==0 for i,k in r['marker_ordinary_edges'])
            cf=Counter()
            for weight,z in zip(w,groups):
                for face,n in z['positive_face_coefficients']:cf[face]+=weight*n
            oldload=sum((alpha.get(f,F())*x for f,x in cf.items()),F());oldvalues.append(oldload)
            cf[r['face']]+=F(1,cap);load=sum((alpha.get(f,F())*x for f,x in cf.items()),F())
            key=tuple((f,x)for f,x in sorted(cf.items())if x)
            record=dict(vertex=v,source_face=r['face'],source_charge=tag,old_orbit_weights=list(map(str,original)),
                omitted_orbits=omitted,retained_orbit_weights=list(map(str,w)),marker_weight=str(F(1,cap)),
                positive_face_coefficients=[[f,str(x)]for f,x in key],vertex_cost='1',old_charge_load=str(oldload),
                load=str(load),baseline_cost=str(mu.get(v,F())),violation=str(load-mu.get(v,F())),
                old_vector_count=len(vv),old_vectors_sha256=L.jhash(vv))
            if key in unique:unique[key]['other_charge_sources'].append(tag)
            else:record['other_charge_sources']=[];unique[key]=record
        cols=list(unique.values());columns.extend(cols)
        newmu[v]=max([mu.get(v,F())]+[F(c['load'])for c in cols])
        beforeviol+=any(x>mu.get(v,F())for x in oldvalues)
        markerinduced+=any(F(c['old_charge_load'])<=mu.get(v,F())<F(c['load'])for c in cols)
        quotas.append(dict(vertex=v,face=r['face'],marker_capacity=cap,ordinary_omitted_orbits=omitted,
            old_vertex_mass=str(mu.get(v,F())),new_vertex_mass=str(newmu[v]),increase=str(newmu[v]-mu.get(v,F()))))
    multi=[]
    assert len(model['prepared']['inherited_multi_vertex_columns'])==14
    for c in model['prepared']['inherited_multi_vertex_columns']:
        load=sum((alpha.get(f,F())*n for f,n in c['integer_face_coefficients']),F())
        old=sum((mu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F())
        new=sum((newmu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F())
        assert load<=old<=new
        multi.append(dict(variable=c['variable'],all_vertex_costs=c['integer_vertex_capacities'],load=str(load),old_cost=str(old),new_cost=str(new)))
    assert all(newmu.get(v,F())>=x for v,x in mu.items())
    num=F(dual['objective_numerator']);den=sum(newmu.values(),F());assert sum(mu.values(),F())==F(dual['normalizer'])
    return dict(columns=columns,target_allowances=quotas,preserved_multi_owner_costs=multi,
        repair=dict(face_multipliers=dual['face_multipliers'],vertex_multipliers=[[v,str(x)]for v,x in sorted(newmu.items())if x],
            objective_numerator=str(num),normalizer=str(den),bound=str(num/den),initial_bound=dual['bound'],
            scope='Least monotone fixed-alpha repair for the exported charge columns only; not the full enlarged local cones or a global optimum.'),
        counts=dict(exported_columns=len(columns),strictly_violated_columns=sum(F(c['violation'])>0 for c in columns),
            raised_target_allowances=sum(F(r['increase'])>0 for r in quotas),
            old_charge_already_violated_targets=beforeviol,marker_induced_violated_targets=markerinduced))


def run():
    start=time.monotonic();inp=D.Q.Inputs();pilot=inp.read(D.TAG);model=inp.read(D.MODEL)
    for path,h in pilot['sources'].items():assert M.digest(path)==h,path
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};g,geo=D.core_geometry(inp)
    family,active=D.family(inp,g,alpha,model);assert normal(M.enc(family))==pilot['family']and normal(M.enc(geo))==pilot['core_geometry']
    theorem=marker_theorem(family,g,active)
    domains,represented=L.current_domains();p=Z.prepare();local={r['vertex']:r for r in M.read(L.TAG)['forest_records']}
    patches={r['vertex']:r for r in M.read(Z.TAG)['forest_records']};frozen={r['vertex']:r for r in pilot['pilots']}
    rows=[]
    for k,row in enumerate(family,1):
        r=process(row,active,g,p,domains,local,patches,frozen);rows.append(r)
        if k%25==0 or r['marker_ordinary_edges']:
            print('RAW76_FAMILY',k,'/235',r['vertex'],r['counts'],'seconds',round(time.monotonic()-start,2),flush=True)
    assert sum(len(r['marker_ordinary_edges'])for r in rows)==6
    exported=export_charges(rows,domains,model)
    paths=set(inp.paths)|set(p['sources'])
    for tag in(L.TAG,B.TAG,D.THIRD):
        paths.update(f'{f}/{tag}.{e}'for f,e in(('scripts','py'),('results','json'),('results','md')))
    for mod in(D,R,A,H,W,M,L,Z,Z.U):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    return normal(M.enc(dict(status='PASS_RAW76_COMPLETE_FAMILY_COUPLED_CHARGES',sources={p:M.digest(p)for p in sorted(paths)},
        core_geometry=geo,family=family,marker_theorem=theorem,forest_records=rows,
        counts=dict(classes=len(rows),existing_ordinary_forests=sum(bool(r['old_nodes'])for r in rows),
            **{key:sum(r['counts'][key]for r in rows)for key in rows[0]['counts']},**exported['counts']),
        exported=exported,optimizer_invocations=0,
        scope='Complete235-class marker family relative to the fixed327-class ordinary domain. All old Q-Q outcomes and columns are preserved. Full finite-frame marker recovery and one occurrence per source quad are inherited from the frozen pilot. No new global solve or completeness claim for enlarged local charge facets.')))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
