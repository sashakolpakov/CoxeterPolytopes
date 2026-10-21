#!/usr/bin/env python3
"""Exact raw56/57 fifth-vertex geometry and complete spectator families.

Synthetic events are actual source quadrilaterals at a uniquely recovered
extra finite vertex. No old pair outcome or optimizer is changed here.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
import networkx as nx
import fresh_raw76_extra_vertex_complete_family_20261018 as U

D=U.D;ROOT=D.ROOT;M=D.M;A=D.A;R=D.R;P=D.P;V=D.V;K=D.K;H=D.H
TAG='fresh_raw56_raw57_extra_vertex_family_20261018'
MODEL='fresh_raw76_complete_family_full_model_20261018'
IDS=(0,1,3,5)
CORES={56:45,57:46}
TAU=K(F(5,2),F(1,2))
normal=D.normal


def geometry(inp,raw_index):
    raw=inp.read(A.M.TAGS['raw']);quads=inp.read(M.QUAD)
    g=A.V.raw_graph(raw['records'][raw_index]);core=quads['canonical_cores'][CORES[raw_index]]
    assert core['source_core']==raw_index and core['required']==1
    assert core['quad_partitions']==[[[1,3],[0,2]]]
    scale=[K(1)]*6;scale[1]=TAU
    if raw_index==56:scale[0]=K(5,2)
    gram=[[2*scale[i]if i==j else K()for j in range(6)]for i in range(6)]
    for i,j in combinations(range(6),2):
        sq=V.squared(g[i][j]);z=A.M.ksqrt(sq*scale[i]*scale[j])if sq else K()
        assert z is not None;gram[i][j]=gram[j][i]=-z
    assert all(s.sign()>0 for s in scale)and P.congruence(gram)['inertia']==[4,1,1]
    records=[]
    for omitted in combinations(range(6),2):
        ids=sorted(set(range(6))-set(omitted))
        records.append(dict(omitted=list(omitted),support=ids,inertia=P.congruence(P.sub(gram,ids))['inertia']))
    pd={tuple(r['omitted'])for r in records if r['inertia']==[4,0,0]}
    quad={tuple(sorted(c['vertex_omitted']))for c in core['aggregate_corner_occurrences']}
    assert pd==set(map(tuple,core['pd_deletions']))and quad=={(0,1),(0,3),(1,2),(2,3)}
    assert pd-quad=={(2,4)}and sorted(set(range(6))-{2,4})==list(IDS)
    # The triangle fixes0,5: its sides1,3,4 have all three vertices actual.
    fixed=[0,5];sides=[1,3,4]
    supports=[sorted(fixed+list(s))for s in combinations(sides,2)]
    assert all(tuple(sorted(set(range(6))-set(s)))in pd for s in supports)
    triangle=P.schur(gram,fixed,sides);inverse,_=V.inverse(triangle)
    assert P.congruence(triangle)['inertia']==[2,1,0]
    assert all(x.sign()<0 for row in inverse for x in row)
    autos=V.automorphisms(g);assert autos==[list(range(6))]
    ng=nx.Graph();ng.add_nodes_from(range(6));ng.add_edges_from((i,j)for i,j in combinations(range(6),2)if g[i][j]!=A.RIGHT)
    assert nx.is_connected(ng)
    active={i:{j:g[i][j][0]for j in IDS if i!=j and g[i][j]!=A.RIGHT}for i in IDS}
    names={R.C.code(h):n for n,h in R.C.models()}
    types=sorted(names[R.C.code(h)]for h in R.C.pieces(active))
    assert types==(['A1','A1','I2(10)']if raw_index==56 else['A2','I2(10)'])
    return g,active,dict(raw_family=raw_index,canonical_core=CORES[raw_index],native_gram=g,
        squared_scales=[x.out()for x in scale],scaled_gram=P.out(gram),inertia=[4,1,1],
        support_records=records,quad_omissions=sorted(map(list,quad)),quad_fixed=[4,5],
        extra_omission=[2,4],extra_support=list(IDS),extra_active_types=types,
        triangle_fixed=fixed,triangle_sides=sides,triangle_supports=supports,
        triangle_residual=P.out(triangle),triangle_inverse=P.out(inverse),native_automorphisms=autos)


def families(inp,alpha,model,geometries):
    forests=inp.read(M.FORESTS);spectators=inp.read(A.B.Q.SPECTATORS);quads=inp.read(M.QUAD)
    actual=inp.read(M.ACTUAL)['provenance_masses'];pm={f:p for p,f in model['prepared']['provenance_to_face']}
    lookup={tuple(r['forest']):r['index']for r in forests['vertices']};models=dict(R.C.models());rows=[]
    for fi,ff in enumerate(forests['faces']):
        if ff['kind']!='canonical_dotted_quad_aggregate':continue
        relevant=[(ctx,ci)for ctx,ci in ff['source_contexts']if ci in CORES.values()]
        if not relevant:continue
        assert len(relevant)==len(ff['source_contexts'])==1
        ctx,ci=relevant[0];raw_index=next(k for k,v in CORES.items()if v==ci)
        pi=pm.get(fi)
        if pi is not None:assert(actual[pi]['source_context'],actual[pi]['canonical_prism_core'])==(ctx,ci)
        sp=spectators['contexts'][ctx]['spectator_forest']
        assert sum(len(models[n])for n in sp)==9
        assert quads['source_core_maps'][spectators['contexts'][ctx]['core']]['canonical_core']==ci
        vf=sorted(sp+geometries[raw_index]['extra_active_types']);v=lookup[tuple(vf)];c=Counter(vf)
        expected=2*c['I2(10)']*(c['A1']*(c['A1']-1)if raw_index==56 else 2*c['A2'])
        rows.append(dict(face=fi,source_context=ctx,provenance=pi,raw_family=raw_index,
            canonical_core=ci,spectator_forest=sp,target_vertex=v,target_forest=vf,
            marker_orientations=expected,actual_extra_occurrences=1,current_mass=str(alpha.get(fi,F()))))
    for raw_index in CORES:
        rr=[r for r in rows if r['raw_family']==raw_index]
        assert len(rr)==len({r['target_vertex']for r in rr})
    return rows


def projection_data(target):
    J,active,components=R.projection_data(target);tens=[]
    for c in R.C.pieces(target):
        if len(c)==2 and list(next(iter(c.values())).values())==[10]:
            i,j=sorted(c);tens.append((i,j));active.update((i,j))
    return J,active,components,tens


def field_for(covectors,extra=()):
    # The inherited squareclass decoder names non-K native magnitudes as
    # dotted encodings. This is only a numeric encoding, not a Coxeter label.
    return A.pair_field(tuple(covectors),tuple(extra)+((0,str(TAU.a),str(TAU.b)),))


def projection(p,q,J,tens,field):
    z=R.projection(p,q,J,field);p=dict(p);q=dict(q)
    den=K(4)-TAU;assert den.sign()>0
    off=tuple(-x/den for x in field.negative_contact((0,str(TAU.a),str(TAU.b))))
    diag=(K(2)/den,K(),K(),K())
    for i,j in tens:
        for a,b,c in((i,i,diag),(j,j,diag),(i,j,off),(j,i,off)):
            if a not in p or b not in q:continue
            term=field.qm(field.qm(field.negative_contact(p[a]),field.negative_contact(q[b])),c)
            z=tuple(x+y for x,y in zip(z,term))
    return z


def marker_maps(row,active,g):
    target=R.P.forest_graph(row['target_forest'],dict(R.C.models()))
    J,ac,blocks,tens=projection_data(target);maps=R.active_embeddings(active,target)
    assert len(maps)==row['marker_orientations'];markers=[];seen=set()
    for mp in maps:
        assert set(mp)==set(IDS)and len(set(mp.values()))==4
        cov=[tuple(sorted((mp[j],tuple(g[i][j]))for j in IDS if g[i][j]!=A.RIGHT))for i in(2,4)]
        assert len(set(cov))==2 and all(i in ac for c in cov for i,e in c)
        heights=[];selfs=[]
        for c in cov:
            field=field_for((c,));p=projection(c,c,J,tens,field)
            assert not any(p[1:]);h=p[0]-2;assert h.sign()>0;heights.append(h);selfs.append(R.qout(p))
        field=field_for(cov,(g[2][4],));cross=projection(cov[0],cov[1],J,tens,field)
        prescribed=field.negative_contact(g[2][4])
        assert field.signs_difference(tuple(a-b for a,b in zip(cross,prescribed)),heights[0]*heights[1])==0
        order=sorted(range(2),key=lambda i:cov[i]);native=dict(contacts=[[[j,list(e)]for j,e in cov[i]]for i in order],
            squared_timelike_heights=[heights[i].out()for i in order],native_omitted_product=list(g[2][4]),occupied_pairs=[])
        key=H.frozen(native['contacts']);assert key not in seen;seen.add(key)
        markers.append(dict(index=len(markers),native=native,active_embedding=sorted(mp.items()),
            triangle_ports=sorted([mp[1],mp[3]]),self_projections=[selfs[i]for i in order],
            mutual_projection=R.qout(cross),height_product=(heights[0]*heights[1]).out()))
    return target,J,tens,markers


def pair_test(J,tens,a,b):
    expressions=[]
    for ai,ca in enumerate(A.native(a)):
        for bj,cb in enumerate(A.native(b)):
            if ca==cb:continue
            field=field_for((ca,cb));heights=[]
            for c in(ca,cb):
                p=projection(c,c,J,tens,field);assert not any(p[1:]);h=p[0]-2
                assert h.sign()>0;heights.append(h)
            pr=projection(ca,cb,J,tens,field);hp=heights[0]*heights[1]
            code,cert=field.arithmetic_wall_gate(tuple((x.a,x.b)for x in pr),hp.a,hp.b)
            sign=field.signs_difference(pr,hp)
            if sign>0:code='positive'
            elif sign<0 and field.signs_difference(tuple(x+(K(1)if k==0 else K())for k,x in enumerate(pr)),hp)>0:code='gap'
            expressions.append(dict(walls=[ai,bj],contacts=normal((ca,cb)),projection=R.qout(pr),
                heights=[h.out()for h in heights],height_product=hp.out(),code=code,gate_certificate=cert))
            if code is not None:return dict(allowed=False,expressions=expressions)
    return dict(allowed=True,expressions=expressions)


def relative(raw_a,ma,raw_b,mb):
    a=dict(ma['active_embedding']);b=dict(mb['active_embedding'])
    ten=0 if(a[1],a[3])==(b[1],b[3])else 1 if(a[1],a[3])==(b[3],b[1])else 2
    if raw_a==raw_b==56:
        names={a[0]:0,a[5]:1}
        for i in(b[0],b[5]):
            if i not in names:names[i]=len(names)
        return(names[b[0]],names[b[5]],ten)
    if raw_a==raw_b==57:
        two=0 if(a[0],a[5])==(b[0],b[5])else 1 if(a[0],a[5])==(b[5],b[0])else 2
        return(two,ten)
    assert(raw_a,raw_b)==(56,57)
    return(ten,)


def allowed_pattern(raw_a,raw_b,pattern):
    if raw_a!=raw_b:return False
    if raw_a==57:return tuple(pattern)==(0,2)
    return tuple(pattern)in((0,1,2),(2,3,0),(2,3,2))


def marker_theorem(family,native,active):
    host_specs=[(56,56,['A1']*7+['I2(10)']*3,20),
        (57,57,['A1']+['A2']*3+['I2(10)']*3,8),
        (56,57,['A1']*7+['A2']+['I2(10)']*2,3)]
    tests=[]
    for ra,rb,forest,count in host_specs:
        row=next(r for r in family if r['raw_family']==ra and Counter(r['target_forest'])==Counter(forest))
        target,J,tens,ma=marker_maps(row,active[ra],native[ra]);anchor=ma[0]
        rowb=next(r for r in family if r['raw_family']==rb and r['target_vertex']==row['target_vertex'])
        _,_,_,mb=marker_maps(rowb,active[rb],native[rb]);reps={}
        for m in mb:
            if ra==rb and m['index']==0:continue
            reps.setdefault(relative(ra,anchor,rb,m),m)
        assert len(reps)==count
        records=[]
        for pattern,other in sorted(reps.items()):
            result=pair_test(J,tens,anchor['native'],other['native'])
            assert result['allowed']==allowed_pattern(ra,rb,pattern)
            records.append(dict(pattern=list(pattern),other_marker=other,result=result))
        tests.append(dict(raw_families=[ra,rb],host_forest=forest,finite_frame=target,anchor=anchor,records=records))
    assert sum(len(r['records'])for r in tests)==31
    return dict(tests=tests,raw56_capacity='floor(a1/2)*i210',raw57_capacity='i210',
        mixed_families_mutually_excluded=True,
        scope='Complete relative-role physical pair-gate classification. Retained cases are unknown, not realized configurations; capacities are exact for this conservative marker graph.')


def automorphism(target,base,other):
    graph,components=U.graph_components(target)
    a=dict(base['active_embedding']);b=dict(other['active_embedding']);perm={a[i]:b[i]for i in a}
    used=set(perm);dest=set(perm.values())
    left=[c for c in components if not set(c)&used];right=[c for c in components if not set(c)&dest]
    assert sum(map(len,left))==sum(map(len,right))==9
    for c in left:
        for k,e in enumerate(right):
            iso=nx.algorithms.isomorphism.GraphMatcher(graph.subgraph(c),graph.subgraph(e),edge_match=lambda a,b:a['label']==b['label'])
            if iso.is_isomorphic():perm.update(next(iso.isomorphisms_iter()));right.pop(k);break
        else:raise AssertionError('unmatched passive component')
    assert not right and set(perm)==set(perm.values())==set(range(13))
    assert all(target[i].get(j,2)==target[perm[i]].get(perm[j],2)for i,j in combinations(range(13),2))
    return perm


def process(rows,native,active,p,domains,local,patches):
    row=rows[0];v=row['target_vertex']
    target=R.P.forest_graph(row['target_forest'],dict(R.C.models()));J,ac,blocks,tens=projection_data(target)
    nodes,groups,graph,owner,forest,history=U.complete_current(v,row,p,domains,local,patches)
    assert forest==row['target_forest'];n=len(nodes)
    tenids=set(sum((list(t)for t in tens),[]))
    assert all(i in ac and i not in tenids for s in nodes for c in s['contacts']for i,e in c)
    canonical=[U.transport_native(s,dict(zip(range(13),range(13))))for s in nodes]
    oldkeys={H.frozen(s['contacts']):i for i,s in enumerate(canonical)};assert len(oldkeys)==n
    family_records=[];allmarkers=[]
    for row in rows:
        ri=row['raw_family'];_,_,_,markers=marker_maps(row,active[ri],native[ri])
        markkeys={H.frozen(s['native']['contacts']):i for i,s in enumerate(markers)}
        assert len(markkeys)==len(markers)and not set(oldkeys)&set(markkeys)
        anchor=[]
        for i,s in enumerate(nodes):
            if markers[0]['triangle_ports']in s['occupied_pairs']:
                result=dict(allowed=False,triangle_port_exclusion=markers[0]['triangle_ports'],expressions=[])
            else:result=pair_test(J,tens,s,markers[0]['native'])
            anchor.append(dict(ordinary=i,result=result))
        allowed={r['ordinary']for r in anchor if r['result']['allowed']};transports=[];cross=[]
        for marker in markers:
            perm=automorphism(target,markers[0],marker);on=[];mn=[]
            for i,s in enumerate(nodes):
                t=U.transport_native(s,perm);j=oldkeys[H.frozen(t['contacts'])]
                assert all(normal(t[k])==normal(canonical[j][k])for k in H.NATIVE_FIELDS)
                assert owner[i]==owner[j];on.append(j)
            for s in markers:
                t=U.transport_native(s['native'],perm);j=markkeys[H.frozen(t['contacts'])]
                assert all(normal(t[k])==normal(markers[j]['native'][k])for k in H.NATIVE_FIELDS);mn.append(j)
            assert len(set(on))==n and len(set(mn))==len(markers)and mn[0]==marker['index']
            cross.extend([on[i],marker['index']]for i in sorted(allowed))
            transports.append(dict(frame_permutation=sorted(perm.items()),ordinary_permutation=on,marker_permutation=mn))
        c=Counter(forest);cap=(c['A1']//2)*c['I2(10)']if ri==56 else c['I2(10)'];assert cap>=1
        # Build an explicit largest compatible marker set. Its members are
        # graph witnesses only; the proof gives the matching/order upper bound.
        if ri==56:
            ones=sorted(i for i in target if not target[i]);pairs=list(zip(ones[::2],ones[1::2]))
            chosen=[m['index']for m in markers if(dict(m['active_embedding'])[0],dict(m['active_embedding'])[5])in pairs
                and dict(m['active_embedding'])[1]<dict(m['active_embedding'])[3]]
        else:
            base=dict(markers[0]['active_embedding'])
            chosen=[m['index']for m in markers if(dict(m['active_embedding'])[0],dict(m['active_embedding'])[5])==(base[0],base[5])
                and dict(m['active_embedding'])[1]<dict(m['active_embedding'])[3]]
        assert len(chosen)==cap
        assert all(allowed_pattern(ri,ri,relative(ri,markers[i],ri,markers[j]))for i,j in combinations(chosen,2))
        family_records.append(dict(face=row['face'],raw_family=ri,markers=markers,anchor_decisions=anchor,
            transports=transports,marker_ordinary_edges=sorted(cross),marker_capacity=cap,attaining_marker_clique=chosen))
        allmarkers.extend((ri,m)for m in markers)
    keys=[H.frozen(m['native']['contacts'])for ri,m in allmarkers];assert len(keys)==len(set(keys))
    counts=dict(ordinary_stars=n,markers=len(allmarkers),anchor_tests=n*len(rows),cross_pairs=n*len(allmarkers),
        marker_pairs=len(allmarkers)*(len(allmarkers)-1)//2,
        retained_cross_pairs=sum(len(r['marker_ordinary_edges'])for r in family_records),
        anchor_unknowns=sum(z['result']['allowed']for r in family_records for z in r['anchor_decisions']))
    return normal(M.enc(dict(vertex=v,forest=forest,finite_frame=target,old_nodes=nodes,old_groups=groups,
        old_graph_edges=sorted(tuple(sorted(e))for e in graph.edges()),prior_history=history,
        marker_families=family_records,counts=counts)))


def export_charges(rows,domains,model):
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};mu={v:F(x)for v,x in dual['vertex_multipliers']}
    certificates=defaultdict(list);vectors={v:d['domain']['maximal_count_vectors']for v,d in domains.items()}
    for tag in(D.L.TAG,D.B.TAG,D.THIRD):
        for r in M.read(tag)['forest_records']:
            certificates[r['vertex']].append((tag,r['certificate']['weights']))
            if tag!=D.L.TAG:vectors[r['vertex']]=r['maximal_count_vectors']
    columns=[];newmu=dict(mu);quotas=[]
    for r in rows:
        v=r['vertex'];groups=r['old_groups'];owner={i:z['index']for z in groups for i in z['literal_stars']}
        omitted=sorted({owner[i]for family in r['marker_families']for i,k in family['marker_ordinary_edges']})
        candidates=certificates[v]if groups else [('marker_only',[])];assert candidates
        unique={};vv=vectors[v]if groups else [[]]
        for tag,weights in candidates:
            original=list(map(F,weights))+[F()]*(len(groups)-len(weights))
            assert len(original)==len(groups)and all(w>=0 for w in original)
            assert all(len(x)==len(original)and sum((w*n for w,n in zip(original,x)),F())<=1 for x in vv)
            w=[F()if i in omitted else x for i,x in enumerate(original)];cf=Counter()
            for weight,z in zip(w,groups):
                for face,n in z['positive_face_coefficients']:cf[face]+=weight*n
            for family in r['marker_families']:
                assert all(w[owner[i]]==0 for i,k in family['marker_ordinary_edges'])
                cf[family['face']]+=F(1,family['marker_capacity'])
            key=tuple((f,x)for f,x in sorted(cf.items())if x);load=sum((alpha.get(f,F())*x for f,x in key),F())
            record=dict(vertex=v,source_charge=tag,old_orbit_weights=list(map(str,original)),
                omitted_orbits=omitted,retained_orbit_weights=list(map(str,w)),
                marker_terms=[[z['raw_family'],z['face'],str(F(1,z['marker_capacity']))]for z in r['marker_families']],
                positive_face_coefficients=[[f,str(x)]for f,x in key],vertex_cost='1',load=str(load),
                baseline_cost=str(mu.get(v,F())),violation=str(load-mu.get(v,F())),
                old_vector_count=len(vv),old_vectors_sha256=D.L.jhash(vv))
            if key in unique:unique[key]['other_charge_sources'].append(tag)
            else:record['other_charge_sources']=[];unique[key]=record
        cols=list(unique.values());columns.extend(cols);newmu[v]=max([mu.get(v,F())]+[F(c['load'])for c in cols])
        quotas.append(dict(vertex=v,old_vertex_mass=str(mu.get(v,F())),new_vertex_mass=str(newmu[v]),
            increase=str(newmu[v]-mu.get(v,F())),ordinary_omitted_orbits=omitted))
    multi=[]
    assert len(model['prepared']['inherited_multi_vertex_columns'])==14
    for c in model['prepared']['inherited_multi_vertex_columns']:
        load=sum((alpha.get(f,F())*n for f,n in c['integer_face_coefficients']),F())
        old=sum((mu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F())
        new=sum((newmu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F());assert load<=old<=new
        multi.append(dict(variable=c['variable'],all_vertex_costs=c['integer_vertex_capacities'],load=str(load),old_cost=str(old),new_cost=str(new)))
    assert all(newmu.get(v,F())>=x for v,x in mu.items())
    prior=[]
    for c in M.read(U.TAG)['exported']['columns']:
        load=sum((alpha.get(f,F())*F(n)for f,n in c['positive_face_coefficients']),F());v=c['vertex']
        assert c['vertex_cost']=='1'and load<=mu.get(v,F())<=newmu.get(v,F())
        prior.append(dict(vertex=v,source_face=c['source_face'],load=str(load),old_cost=str(mu.get(v,F())),new_cost=str(newmu.get(v,F()))))
    assert len(prior)==242
    num=F(dual['objective_numerator']);den=sum(newmu.values(),F());assert sum(mu.values(),F())==F(dual['normalizer'])
    return dict(columns=columns,target_allowances=quotas,preserved_multi_owner_costs=multi,preserved_raw76_columns=prior,
        repair=dict(face_multipliers=dual['face_multipliers'],vertex_multipliers=[[v,str(x)]for v,x in sorted(newmu.items())if x],
            objective_numerator=str(num),normalizer=str(den),bound=str(num/den),initial_bound=dual['bound'],
            scope='Least monotone fixed-alpha repair for exported columns only; neither complete enlarged local cones nor a global optimum.'),
        counts=dict(exported_columns=len(columns),strictly_violated_columns=sum(F(c['violation'])>0 for c in columns),
            raised_target_allowances=sum(F(r['increase'])>0 for r in quotas)))


def run():
    start=time.monotonic();inp=D.Q.Inputs();model=inp.read(MODEL)
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};geometries={};native={};active={}
    for ri in CORES:native[ri],active[ri],geometries[ri]=geometry(inp,ri)
    assert not any(all(native[56][i][j]==native[57][p[i]][p[j]]for i,j in combinations(range(6),2))for p in permutations(range(6)))
    rows=families(inp,alpha,model,geometries);theorem=marker_theorem(rows,native,active)
    domains,represented=D.L.current_domains();p=D.Z.prepare()
    local={r['vertex']:r for r in M.read(D.L.TAG)['forest_records']};patches={r['vertex']:r for r in M.read(D.Z.TAG)['forest_records']}
    grouped=defaultdict(list)
    for r in rows:grouped[r['target_vertex']].append(r)
    preflight=dict(targets=len(grouped),ordinary_targets=sum(v in domains for v in grouped),
        ordinary_stars=sum(sum(len(g['literal_stars'])for g in domains[v]['groups'])for v in grouped if v in domains),
        anchor_tests=sum(len(rr)*sum(len(g['literal_stars'])for g in domains[v]['groups'])for v,rr in grouped.items()if v in domains),
        cross_pairs=sum(sum(r['marker_orientations']for r in rr)*sum(len(g['literal_stars'])for g in domains[v]['groups'])for v,rr in grouped.items()if v in domains),
        marker_pairs=sum((m:=sum(r['marker_orientations']for r in rr))*(m-1)//2 for rr in grouped.values()))
    assert preflight==dict(targets=1212,ordinary_targets=91,ordinary_stars=2883,anchor_tests=4115,cross_pairs=108604,marker_pairs=577768)
    print('RAW56_RAW57_PREFLIGHT',preflight,flush=True);mapped=[]
    for k,(v,rr)in enumerate(sorted(grouped.items()),1):
        out=process(rr,native,active,p,domains,local,patches);mapped.append(out)
        if k%100==0 or out['counts']['anchor_unknowns']:
            print('RAW56_RAW57_FAMILY',k,'/1212',v,out['counts'],'seconds',round(time.monotonic()-start,2),flush=True)
    exported=export_charges(mapped,domains,model)
    counts=dict(family_classes=len(rows),marker_orientations=sum(r['marker_orientations']for r in rows),
        target_forests=len({r['target_vertex']for r in rows}),positive_classes=sum(F(r['current_mass'])>0 for r in rows))
    byraw={str(k):dict(classes=sum(r['raw_family']==k for r in rows),
        orientations=sum(r['marker_orientations']for r in rows if r['raw_family']==k),
        positive_classes=sum(r['raw_family']==k and F(r['current_mass'])>0 for r in rows),
        positive_mass=str(sum((F(r['current_mass'])for r in rows if r['raw_family']==k),F())))for k in CORES}
    paths=set(inp.paths)|set(p['sources'])
    for tag in(D.L.TAG,D.B.TAG,D.THIRD,U.TAG):
        paths.update(f'{f}/{tag}.{e}'for f,e in(('scripts','py'),('results','json'),('results','md')))
    for mod in(U,D,A,R,P,V,D.H,D.W,D.M,D.L,D.Z,D.Z.U):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    for tag in('fresh_forest1325_cross_vertex_incidence_20261017','fresh_quad_product_aggregation_20261017'):
        paths.add(f'results/{tag}.md')
    paths.add(str(Path(__file__).relative_to(ROOT)))
    counts.update(exported['counts']);counts.update({key:sum(r['counts'][key]for r in mapped)for key in mapped[0]['counts']})
    assert counts==dict(family_classes=1408,marker_orientations=21584,target_forests=1212,positive_classes=75,
        exported_columns=1212,strictly_violated_columns=60,raised_target_allowances=60,
        ordinary_stars=2883,markers=21584,anchor_tests=4115,cross_pairs=108604,marker_pairs=577768,
        retained_cross_pairs=13104,anchor_unknowns=414)
    return normal(M.enc(dict(status='PASS_RAW56_RAW57_COMPLETE_EXTRA_VERTEX_FAMILY',
        sources={p:M.digest(p)for p in sorted(paths)},geometries=geometries,family=rows,
        marker_theorem=theorem,preflight=preflight,forest_records=mapped,exported=exported,
        counts=counts,by_raw_family=byraw,optimizer_invocations=0,current_graphs_changed=False,
        scope='Complete1408-class extra-vertex families relative to the fixed327 ordinary-class domain. Every oldQ-Q outcome and prior column is preserved. Capacities concern conservative marker graphs; exported one-vertex coupled charges are valid independently of masses. No complete enlarged-cone or global-optimum claim.')))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:
        assert not path.exists()
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['by_raw_family'],flush=True)
