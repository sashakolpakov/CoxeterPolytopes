#!/usr/bin/env python3
"""Raw76 complete spectator-family inventory and a bounded extra-vertex pilot.

Synthetic markers denote actual raw76 quads at their unique fifth PD support.
No new ordinary quad corner is invented. Old Q-Q graphs are kept unchanged.
This staged diagnostic tests three marker domains; no global solve.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_new5073_raw76_high_mass_local_covers_20261018 as B
import fresh_raw39_shared_edge_target_incidence_20261018 as E

Q=B.Q;A=B.A;R=B.R;L=B.L;M=B.M;Z=B.Z;H=B.H;W=B.W;ROOT=B.ROOT
P=E.P;K=E.K;V=E.V
TAG='fresh_raw76_forced_extra_vertex_family_20261018'
MODEL='fresh_current327_full_model_20261018'
THIRD='fresh_new5073_raw76_third_class_local_covers_20261018'
PILOT=(43092,43012,43010)
IDS=(0,1,3,4,6)
normal=B.normal


def core_geometry(inp):
    raw=inp.read(A.M.TAGS['raw']);quads=inp.read(M.QUAD);g=A.V.raw_graph(raw['records'][76])
    core=quads['canonical_cores'][65];assert core['source_core']==76 and core['quad_partitions']==[[[1,3],[0,2]]]
    scale=[K(F(7,2),F(1,2))]+[K(1)]*6
    gram=[[2*scale[i]if i==j else K()for j in range(7)]for i in range(7)]
    for i,j in combinations(range(7),2):
        sq=V.squared(g[i][j]);z=A.M.ksqrt(sq*scale[i]*scale[j])if sq else K()
        assert z is not None;gram[i][j]=gram[j][i]=-z
    assert P.congruence(gram)['inertia']==[5,1,1]
    supports=[]
    for omitted in combinations(range(7),2):
        ids=sorted(set(range(7))-set(omitted));inertia=P.congruence(P.sub(gram,ids))['inertia']
        supports.append(dict(omitted=list(omitted),support=ids,inertia=inertia))
    pd={tuple(r['omitted'])for r in supports if r['inertia']==[5,0,0]}
    quad={tuple(sorted(c['vertex_omitted']))for c in core['aggregate_corner_occurrences']}
    assert pd==set(map(tuple,core['pd_deletions'])) and pd-quad=={(2,5)}
    assert quad=={(0,1),(0,3),(1,2),(2,3)}
    assert all({4,5,6}<=set(range(7))-set(q)for q in quad)
    assert sorted(set(range(7))-{2,5})==list(IDS)
    triangle=P.schur(gram,[0,4,6],[1,3,5]);inverse,_=V.inverse(triangle)
    assert P.congruence(triangle)['inertia']==[2,1,0]and all(x.sign()<0 for row in inverse for x in row)
    tri_supports=[sorted([0,4,6]+list(s))for s in combinations([1,3,5],2)]
    assert all(tuple(sorted(set(range(7))-set(s)))in pd for s in tri_supports)
    ng=nx.Graph();ng.add_nodes_from(range(7));ng.add_edges_from((i,j,{'native':g[i][j]})for i,j in combinations(range(7),2))
    autos=list(nx.algorithms.isomorphism.GraphMatcher(ng,ng,edge_match=lambda a,b:a['native']==b['native']).isomorphisms_iter())
    assert autos==[dict(zip(range(7),range(7)))]
    return g,dict(native_gram=[[list(e)for e in row]for row in g],squared_scales=[x.out()for x in scale],
        scaled_gram=P.out(gram),inertia=[5,1,1],support_records=supports,quad_omissions=sorted(map(list,quad)),
        quad_fixed=[4,5,6],extra_omission=[2,5],extra_support=list(IDS),triangle_fixed=[0,4,6],triangle_sides=[1,3,5],
        triangle_supports=tri_supports,triangle_residual=P.out(triangle),triangle_inverse=P.out(inverse),
        native_automorphisms=[sorted(a.items())for a in autos])


def family(inp,g,alpha,model):
    forests=inp.read(M.FORESTS);spectators=inp.read(A.B.Q.SPECTATORS);quads=inp.read(M.QUAD)
    actual=inp.read(M.ACTUAL)['provenance_masses'];pm={f:p for p,f in model['prepared']['provenance_to_face']}
    lookup={tuple(r['forest']):r['index']for r in forests['vertices']};models=dict(R.C.models())
    names={R.C.code(h):n for n,h in models.items()}
    active={i:{j:g[i][j][0]for j in IDS if i!=j and g[i][j]!=A.RIGHT}for i in IDS}
    assert sorted(names[R.C.code(h)]for h in R.C.pieces(active))==['A1','A1','H3']
    rows=[]
    for fi,ff in enumerate(forests['faces']):
        if ff['kind']!='canonical_dotted_quad_aggregate' or not any(ci==65 for _,ci in ff['source_contexts']):continue
        assert len(ff['source_contexts'])==1;ctx,ci=ff['source_contexts'][0];assert ci==65
        pi=pm.get(fi)
        if pi is not None:assert (actual[pi]['source_context'],actual[pi]['canonical_prism_core'])==(ctx,65)
        sp=spectators['contexts'][ctx]['spectator_forest'];assert sum(len(models[n])for n in sp)==8
        assert quads['source_core_maps'][spectators['contexts'][ctx]['core']]['canonical_core']==65
        vf=sorted(sp+['A1','A1','H3']);v=lookup[tuple(vf)];co=Counter(vf)
        expected=co['A1']*(co['A1']-1)*co['H3']
        rows.append(dict(face=fi,source_context=ctx,provenance=pi,spectator_forest=sp,target_vertex=v,
            target_forest=vf,marker_orientations=expected,current_mass=str(alpha.get(fi,F())),actual_extra_occurrences=1))
    assert len(rows)==len({r['target_vertex']for r in rows})==235
    assert sum(r['marker_orientations']for r in rows)==3204
    return rows,active


def marker_maps(row,active,g):
    target=R.P.forest_graph(row['target_forest'],dict(R.C.models()));J,ac,blocks=R.projection_data(target)
    maps=R.active_embeddings(active,target);assert len(maps)==row['marker_orientations']
    markers=[];seen=set()
    for mp in maps:
        assert set(mp)==set(IDS) and len(set(mp.values()))==5
        cov=[tuple((mp[j],tuple(g[i][j]))for j in IDS if g[i][j]!=A.RIGHT)for i in(2,5)]
        cov=[tuple(sorted(c))for c in cov];assert len(set(cov))==2
        hs=[];selfs=[]
        for c in cov:
            field=A.pair_field((c,));p=R.projection(c,c,J,field);assert not any(p[1:]);h=p[0]-2;assert h.sign()>0
            hs.append(h);selfs.append(R.qout(p))
        field=A.pair_field(tuple(cov),(g[2][5],));cross=R.projection(cov[0],cov[1],J,field)
        prescribed=field.negative_contact(g[2][5]);assert field.signs_difference(tuple(a-b for a,b in zip(cross,prescribed)),hs[0]*hs[1])==0
        order=sorted(range(2),key=lambda i:cov[i]);native=dict(contacts=[[[j,list(e)]for j,e in cov[i]]for i in order],
            squared_timelike_heights=[hs[i].out()for i in order],native_omitted_product=list(g[2][5]),occupied_pairs=[])
        key=H.frozen(native['contacts']);assert key not in seen;seen.add(key)
        markers.append(dict(index=len(markers),native=native,active_embedding=sorted(mp.items()),
            triangle_ports=sorted([mp[1],mp[3]]),self_projections=[selfs[i]for i in order],
            mutual_projection=R.qout(cross),height_product=(hs[0]*hs[1]).out()))
    return target,J,markers


def current(v,p,local,domains,patches):
    nodes,groups,graph,owner,forest,baseline=B.current(v,p,local,domains,patches)
    # Complete chronological newest raw76 patches, including explicit unknowns.
    applied=[]
    for tag in(B.TAG,THIRD):
        doc=M.read(tag)
        for r in doc['forest_records']:
            if r['vertex']!=v:continue
            u=r['domain'];n=len(nodes);assert u['literal_stars'][:n]==nodes
            assert u['corner_orbits'][:len(groups)]==groups
            nodes=normal(u['literal_stars']);groups=normal(u['corner_orbits']);graph.add_nodes_from(range(n,len(nodes)))
            for i,j,code in r['physical']['newly_owned_pair_domain']:
                assert j>=n
                if code is None:graph.add_edge(i,j)
            assert L.jhash(sorted(tuple(sorted(e))for e in graph.edges()))==r['graph_edges_sha256']
            applied.append(tag)
    owner={i:z['index']for z in groups for i in z['literal_stars']}
    return nodes,groups,graph,owner,forest,applied


def physical(target,J,nodes,markers,graph):
    oldn=len(nodes);allnodes=nodes+[m['native']for m in markers];graph.add_nodes_from(range(oldn,len(allnodes)))
    oldkeys={H.frozen(n['contacts'])for n in nodes}
    assert not oldkeys&{H.frozen(m['native']['contacts'])for m in markers}
    decisions=[];expressions=[];cache={}
    for i,j in combinations(range(len(allnodes)),2):
        if j<oldn:continue
        bad=None
        if i<oldn and markers[j-oldn]['triangle_ports']in nodes[i]['occupied_pairs']:
            bad=dict(kind='actual_triangle_port',ports=markers[j-oldn]['triangle_ports'])
        else:
            for ai,ca in enumerate(A.native(allnodes[i])):
                for bj,cb in enumerate(A.native(allnodes[j])):
                    if ca==cb:continue
                    key=tuple(sorted((ca,cb)))
                    if key not in cache:
                        field=A.pair_field(key);hs=[]
                        for c in key:
                            q=R.projection(c,c,J,field);assert not any(q[1:]);h=q[0]-2;assert h.sign()>0;hs.append(h)
                        proj=R.projection(key[0],key[1],J,field);hp=hs[0]*hs[1]
                        code,cert=field.arithmetic_wall_gate(tuple((x.a,x.b)for x in proj),hp.a,hp.b)
                        sign=field.signs_difference(proj,hp)
                        # Acute Coxeter product must be nonpositive and cannot
                        # lie strictly between the norm2 values-1 and0.
                        if sign>0:code='positive'
                        elif sign<0 and field.signs_difference(tuple(x+(K(1)if k==0 else K())for k,x in enumerate(proj)),hp)>0:code='gap'
                        index=len(expressions);expressions.append(dict(contacts=normal(key),projection=R.qout(proj),
                            heights=[h.out()for h in hs],height_product=hp.out(),code=code,gate_certificate=cert))
                        cache[key]=(code,index)
                    code,ei=cache[key]
                    if code is not None:bad=dict(kind='physical_product',wall_indices=[ai,bj],expression=ei,code=code);break
                if bad:break
        if bad is None:graph.add_edge(i,j)
        decisions.append(dict(pair=[i,j],exclusion=bad))
    return decisions,expressions


def pilot(row,active,g,p,local,domains,patches,alpha,mu):
    target,J,markers=marker_maps(row,active,g);v=row['target_vertex']
    nodes,groups,graph,owner,forest,applied=current(v,p,local,domains,patches);assert forest==row['target_forest']
    oldedges=sorted(tuple(sorted(e))for e in graph.edges());oldn=len(nodes);oldclique_count=0
    decisions,expressions=physical(target,J,nodes,markers,graph)
    assert sorted(tuple(sorted(e))for e in graph.edges()if max(e)<oldn)==oldedges
    # Exact mass-independent marker-only capacity and the strongest count
    # envelope for each already proved nonnegative local charge.
    candidates=[]
    for tag in(L.TAG,B.TAG,THIRD):
        for r in M.read(tag)['forest_records']:
            if r['vertex']==v:
                w=list(map(F,r['certificate']['weights']))
                if len(w)<len(groups):w +=[F()]*(len(groups)-len(w))
                assert len(w)==len(groups);candidates.append((tag,w))
    pairs_by_count=[{}for _ in candidates];cliques=[];start=time.monotonic();maxk=0
    for cl in nx.find_cliques(graph):
        assert len(cliques)<3000000 and time.monotonic()-start<30
        old=[i for i in cl if i<oldn];k=sum(i>=oldn for i in cl);maxk=max(maxk,k)
        counts=Counter(owner[i]for i in old);vec=[counts.get(j,0)for j in range(len(groups))]
        cliques.append(sorted(cl))
        for dd,(_,w)in zip(pairs_by_count,candidates):
            load=sum((x*n for x,n in zip(w,vec)),F());assert load<=1
            if k not in dd or load>dd[k][0]:dd[k]=(load,vec,sorted(cl))
    assert maxk>0
    charges=[]
    for (source,w),maxima in zip(candidates,pairs_by_count):
        # Vertices of lambda>=0,t>=0,lambda*x_k+t*k<=1.
        constraints=[(x,k)for k,(x,_,_)in maxima.items()]
        points={(F(),F()),(F(),F(1,maxk))}
        positive=[x for x,k in constraints if x]
        if positive:points.add((1/max(positive),F()))
        for (a,b),(c,d)in combinations(constraints,2):
            det=a*d-b*c
            if det:points.add((F(d-b)/det,(a-c)/det))
        for lam,t in sorted(points):
            if lam<0 or t<=0 or any(lam*x+t*k>1 for x,k in constraints):continue
            coeff=Counter({row['face']:t})
            for x,z in zip(w,groups):
                for f,n in z['positive_face_coefficients']:coeff[f]+=lam*x*n
            coeffs=[[f,str(x)]for f,x in sorted(coeff.items())if x]
            load=sum((alpha.get(f,F())*F(x)for f,x in coeffs),F())
            charges.append(dict(old_charge_source=source,lambda_=str(lam),marker_weight=str(t),
                positive_face_coefficients=coeffs,vertex=v,vertex_cost='1',load=str(load),
                baseline_cost=str(mu.get(v,F())),violation=str(load-mu.get(v,F()))))
    best=max(charges,key=lambda z:F(z['violation']))
    return normal(M.enc(dict(face=row['face'],vertex=v,forest=forest,finite_frame=target,old_nodes=nodes,old_groups=groups,
        prior_raw76_patches=applied,old_graph_edges=oldedges,markers=markers,new_pair_decisions=decisions,wall_expressions=expressions,
        graph_edges=sorted(tuple(sorted(e))for e in graph.edges()),maximal_cliques=cliques,marker_capacity=maxk,
        charge_envelopes=[dict(source=s,weights=list(map(str,w)),count_maxima=[dict(markers=k,old_charge=str(x),vector=vec,clique=cl)
            for k,(x,vec,cl)in sorted(dd.items())])for(s,w),dd in zip(candidates,pairs_by_count)],
        charges=charges,best_charge=best,counts=dict(old_stars=oldn,markers=len(markers),new_pair_positions=len(decisions),
            exclusions=sum(z['exclusion']is not None for z in decisions),unknowns=sum(z['exclusion']is None for z in decisions),maximal_cliques=len(cliques)))))


def run():
    inp=Q.Inputs();model=inp.read(MODEL);dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};mu={v:F(x)for v,x in dual['vertex_multipliers']}
    assert len(alpha)==4624
    for p,h in model['prepared']['sources'].items():assert M.digest(p)==h,p
    g,geo=core_geometry(inp);family_rows,active=family(inp,g,alpha,model)
    selected=sorted(family_rows,key=lambda r:(-F(r['current_mass']),r['face']))[:3];assert tuple(r['face']for r in selected)==PILOT
    domains,represented=L.current_domains();prepared=Z.prepare();patches={r['vertex']:r for r in M.read(Z.TAG)['forest_records']}
    local={r['vertex']:r for r in M.read(L.TAG)['forest_records']};rows=[]
    for r in selected:
        out=pilot(r,active,g,prepared,local,domains,patches,alpha,mu);rows.append(out)
        print('RAW76_EXTRA_PILOT',out['face'],out['vertex'],out['counts'],'best_violation',out['best_charge']['violation'],flush=True)
    paths=set(inp.paths)
    for tag in(L.TAG,B.TAG,THIRD,Z.TAG):
        for f,e in(('scripts','py'),('results','json'),('results','md')):paths.add(f'{f}/{tag}.{e}')
    for tag in('fresh_forest1325_cross_vertex_incidence_20261017','fresh_simplex_disjoint_ports_20261017',
        'fresh_quad_product_aggregation_20261017','fresh_sqrt5_prism_simplex_base_incidence_20261017'):
        paths.add(f'results/{tag}.md')
    for mod in(B,Q,A,R,L,M,Z,H,W,E,P,V):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    return normal(M.enc(dict(status='PASS_RAW76_FORCED_EXTRA_VERTEX_BOUNDED_FAMILY_PILOT',sources={p:M.digest(p)for p in sorted(paths)},
        core_geometry=geo,family=family_rows,pilots=rows,counts=dict(family_classes=235,extra_marker_orientations=3204,
            family_marker_pair_positions=39966,positive_family_classes=sum(F(r['current_mass'])>0 for r in family_rows),
            positive_family_mass=sum((F(r['current_mass'])for r in family_rows),F()),pilot_classes=3,
            physical_markers=sum(r['counts']['markers']for r in rows),
            new_pair_positions=sum(r['counts']['new_pair_positions']for r in rows),
            exclusions=sum(r['counts']['exclusions']for r in rows),
            retained_unknowns=sum(r['counts']['unknowns']for r in rows),
            maximal_cliques=sum(r['counts']['maximal_cliques']for r in rows),
            strictly_violated_targets=sum(F(r['best_charge']['violation'])>0 for r in rows)),
        current_global_model=MODEL,current_represented_classes=327,current_graphs_changed=False,
        optimizer_invocations=0,new_dimension_bound=None,
        scope='Complete235-class spectator metadata plus three exact synthetic-extra-vertex marker domains. Actual pointed-quad incidence uses the accompanying strict-center and complete-frame recovery proof. Only saved old graph charges and their exact2D count envelopes are considered. No new global solve, no all-family physical census or completeness of mixed charge facets.')))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not a.verify:assert not path.exists()
    out=run()
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
