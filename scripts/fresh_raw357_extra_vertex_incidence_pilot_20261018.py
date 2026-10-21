#!/usr/bin/env python3
"""Bounded exact raw357 extra-vertex geometry and complete spectator inventory.

No pair test, marker capacity, local charge, covering LP or global solve.
Complete native recovery is tested at three highest-mass target domains;
all47 source classes and141 target-role domains are inventoried separately.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018 as C

ROOT=C.ROOT;A=C.A;R=C.R;V=C.V;P=C.P;K=C.K;D=C.D;M=C.M
TAG='fresh_raw357_extra_vertex_incidence_pilot_20261018'
MODEL='fresh_all_marker_families_dual_continuation_20261018'
SUPPORT='fresh_all_marker_families_model_support_20261018'
PRIOR=('fresh_raw76_extra_vertex_complete_family_20261018','fresh_raw56_raw57_extra_vertex_family_20261018',
       'fresh_raw108_complete_marker_family_20261018','fresh_raw112_complete_marker_family_20261018')


def geometry(raw,canonical):
    c=canonical[318];assert c['source_core']==357 and c['required']==1
    g=A.V.raw_graph(raw['records'][357]);n=len(g);assert n==10
    G,sc=C.scaled(g)
    assert P.congruence(G)['inertia']==[8,1,1]
    conjugate=[[K(x.a,-x.b)for x in row]for row in G]
    assert P.congruence(conjugate)['inertia']==[9,0,1]
    records=[];pd=set()
    for om in combinations(range(n),2):
        ids=sorted(set(range(n))-set(om));inertia=P.congruence(P.sub(G,ids))['inertia']
        row=dict(omitted=list(om),support=ids,inertia=inertia)
        if inertia==[8,0,0]:pd.add(om);row['active_forest']=C.forest_type(g,ids)
        records.append(row)
    assert sorted(map(list,pd))==c['pd_deletions']
    used={tuple(sorted(x['vertex_omitted']))for x in c['aggregate_corner_occurrences']}
    assert len(c['aggregate_corner_occurrences'])==len(used)==4
    assert pd-used=={(2,4),(2,5),(2,6)}
    assert V.automorphisms(g)==[list(range(n))]
    simplex=[]
    for e in range(n):
        sides=sorted(next(i for i in om if i!=e)for om in pd if e in om)
        if len(sides)<3:continue
        fixed=sorted(set(range(n))-set(sides)-{e})
        assert P.congruence(P.sub(G,fixed))['inertia']==[len(fixed),0,0]
        sg=P.schur(G,fixed,sides);inv,_=V.inverse(sg)
        assert P.congruence(sg)['inertia']==[len(sides)-1,1,0]and all(x.sign()<0 for row in inv for x in row)
        supports=[sorted(fixed+list(set(sides)-{j}))for j in sides]
        assert all(tuple(sorted(set(range(n))-set(s)))in pd for s in supports)
        simplex.append(dict(excluded_wall=e,fixed=fixed,sides=sides,dimension=len(sides)-1,
            vertex_supports=supports,residual_gram=P.out(sg),inverse=P.out(inv)))
    assert len(simplex)==1 and simplex[0]['fixed']==[0,7,8,9]and simplex[0]['sides']==[1,3,4,5,6]
    extra=[]
    for om in sorted(pd-used):
        ids=sorted(set(range(n))-set(om));J,_=V.inverse(P.sub(G,ids));cov=[[G[i][j]for j in ids]for i in om]
        pr=[[sum((cov[i][a]*J[a][b]*cov[j][b]for a in range(8)for b in range(8)),K())for j in range(2)]for i in range(2)]
        hs=[pr[i][i]-G[om[i]][om[i]]for i in range(2)];z=pr[0][1]-G[om[0]][om[1]]
        assert all(h.sign()>0 for h in hs)and z.sign()>0 and z*z==hs[0]*hs[1]
        ports=sorted(set(simplex[0]['sides'])&set(ids));assert len(ports)==4
        extra.append(dict(omitted=list(om),support=ids,active_forest=C.forest_type(g,ids),
            scaled_covectors=[[x.out()for x in cv]for cv in cov],scaled_projection=P.out(pr),
            scaled_heights=[h.out()for h in hs],positive_height_product_root=z.out(),
            forced_triangle_port_pairs=list(combinations(ports,2))))
    assert[e['active_forest']for e in extra]==[['A2','A2','H4'],['A1','A2','A2','H3'],['A2','A2','A2','I2(5)']]
    return g,dict(raw_family=357,canonical_core=318,nodes=n,spectator_rank=5,required_quads=1,
        native_gram=g,scaled_gram=P.out(G),squared_scales=[x.out()for x in sc],
        inertia=[8,1,1],conjugate_inertia=[9,0,1],support_records=records,
        counted_vertex_omissions=sorted(used),ordinary_corner_occurrences=c['aggregate_corner_occurrences'],
        quad_partitions=c['quad_partitions'],extra_vertices=extra,compact_simplex_faces=simplex,
        native_automorphisms=[list(range(n))])


def inventory(inp,g,geom,support):
    forests=inp.read(M.FORESTS);spectators=inp.read(A.B.Q.SPECTATORS);quads=inp.read(M.QUAD)
    positive={r['face']:r for r in support['positive_faces']};models=dict(R.C.models())
    lookup={tuple(r['forest']):r['index']for r in forests['vertices']};rows=[];markers=[]
    for ff in forests['faces']:
        relevant=[(ctx,ci)for ctx,ci in ff.get('source_contexts',[])if ci==318]
        if not relevant:continue
        assert ff['kind']=='canonical_dotted_quad_aggregate'and len(relevant)==len(ff['source_contexts'])==1
        ctx,ci=relevant[0];sp=spectators['contexts'][ctx]['spectator_forest'];fi=ff['index']
        assert sum(len(models[n])for n in sp)==5 and quads['source_core_maps'][spectators['contexts'][ctx]['core']]['canonical_core']==318
        mass=F(positive.get(fi,{}).get('mass','0'))
        if mass:assert positive[fi]['raw_family']==357 and positive[fi]['source_context']==ctx and positive[fi]['required']==1
        rows.append(dict(face=fi,source_context=ctx,canonical_core=318,raw_family=357,spectator_forest=sp,
            required_quads=1,current_mass=str(mass),extra_vertex_occurrences=3))
        for ex in geom['extra_vertices']:
            vf=sorted(sp+ex['active_forest']);v=lookup[tuple(vf)];count=C.expected_maps(ex['active_forest'],vf)
            cc=Counter(vf);b=cc['A2']
            formula=4*b*(b-1)*cc['H4']if ex['omitted']==[2,4]else 4*cc['A1']*b*(b-1)*cc['H3']if ex['omitted']==[2,5]else 16*b*(b-1)*(b-2)*cc['I2(5)']
            assert count==formula
            markers.append(dict(face=fi,source_context=ctx,omitted=ex['omitted'],spectator_forest=sp,
                active_forest=ex['active_forest'],target_vertex=v,target_forest=vf,
                active_role_orientations=count,actual_extra_occurrences=1,current_mass=str(mass)))
    assert len(rows)==47 and len(markers)==141
    positive_rows=[r for r in rows if F(r['current_mass'])>0]
    assert len(positive_rows)==9 and sum(F(r['current_mass'])for r in positive_rows)==F(83284683,10**10)
    chosen=min(positive_rows,key=lambda r:(-F(r['current_mass']),r['face']));assert chosen['face']==70408
    pilots=[]
    for row in[r for r in markers if r['face']==chosen['face']]:
        ex=next(x for x in geom['extra_vertices']if x['omitted']==row['omitted']);ids=ex['support']
        target=R.P.forest_graph(row['target_forest'],models);J,_,_=R.projection_data(target)
        maps=R.active_embeddings(C.active_graph(g,ids),target);assert len(maps)==row['active_role_orientations']
        saved=[];seen=set()
        for mp in maps:
            cov=[tuple(sorted((mp[j],tuple(g[i][j]))for j in ids if g[i][j]!=A.RIGHT))for i in row['omitted']]
            assert len(set(cov))==2;key=tuple(sorted(cov));assert key not in seen;seen.add(key)
            field=A.pair_field(tuple(cov),(g[row['omitted'][0]][row['omitted'][1]],));hs=[]
            for cv in cov:
                z=R.projection(cv,cv,J,field);assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;hs.append(h)
            z=R.projection(cov[0],cov[1],J,field);prescribed=field.negative_contact(g[row['omitted'][0]][row['omitted'][1]])
            assert field.signs_difference(tuple(x-y for x,y in zip(z,prescribed)),hs[0]*hs[1])==0
            saved.append(dict(active_embedding=sorted(mp.items()),full_frame_covectors=cov,
                positive_heights=[h.out()for h in hs],mutual_projection=R.qout(z)))
        pilots.append(dict(**row,finite_frame=target,complete_active_maps=saved))
    return rows,markers,pilots


def workload(inp,markers):
    # Current ordinary327/600 chronology, plus all five marker source families.
    ordinary={};cache={}
    for tag in('fresh_current324_new5073_local_covers_20261018','fresh_new5073_raw76_high_mass_local_covers_20261018','fresh_new5073_raw76_third_class_local_covers_20261018'):
        doc=inp.read(tag)
        for r in doc['forest_records']:
            if 'domain'in r:q=r['domain']
            else:
                name=r['group_source']
                if name not in cache:cache[name]={q['vertex']:q for q in inp.read(name)['forest_records']}
                q=cache[name][r['vertex']];q=q.get('domain',q.get('union',q))
            groups=q['corner_orbits'];keys=[i for z in groups for i in z['literal_stars']]
            assert len(keys)==len(set(keys))and set(keys)==set(range(len(keys)))
            ordinary[r['vertex']]=len(keys)
    assert len(ordinary)==600
    old=Counter()
    for tag in PRIOR:
        for r in inp.read(tag)['forest_records']:old[r['vertex']]+=r['counts']['markers']
    new=Counter()
    for r in markers:new[r['target_vertex']]+=r['active_role_orientations']
    records=[dict(vertex=v,new_markers=n,old_ordinary_stars=ordinary.get(v,0),old_family_markers=old[v],
        new_ordinary_pairs=n*ordinary.get(v,0),new_old_marker_pairs=n*old[v],new_marker_pairs=n*(n-1)//2)for v,n in sorted(new.items())]
    return dict(records=records,target_forests=len(records),
        new_ordinary_pairs=sum(r['new_ordinary_pairs']for r in records),
        new_old_marker_pairs=sum(r['new_old_marker_pairs']for r in records),
        new_marker_pairs=sum(r['new_marker_pairs']for r in records),
        scope='Unreduced pair positions only. All old outcomes/unknowns stay unchanged. Relative-role reduction and capacity proofs are still required.')


def run():
    inp=D.Q.Inputs();raw=inp.read(A.M.TAGS['raw']);quad=inp.read(M.QUAD)
    support=inp.read(SUPPORT);model=inp.read(MODEL)
    assert support['sources'][f'results/{MODEL}.json']==sha256((ROOT/f'results/{MODEL}.json').read_bytes()).hexdigest()
    assert F(model['exact_interval']['lower'])>12 and support['face_masses_changed']
    g,geom=geometry(raw,quad['canonical_cores']);rows,markers,pilots=inventory(inp,g,geom,support)
    work=workload(inp,markers)
    paths=set(inp.paths)|{str(Path(__file__).relative_to(ROOT))}
    for mod in(C,C.T,C.T.U,D,R,V,P,A,A.M,R.C):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    for tag in('fresh_all_marker_families_full_model_20261018',C.TAG,PRIOR[-1]):paths.add(f'results/{tag}.md')
    sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)}
    return D.normal(dict(status='PASS_BOUNDED_RAW357_EXTRA_VERTEX_INCIDENCE_PILOT',sources=sources,geometry=geom,
        complete_family_inventory=rows,extra_marker_inventory=markers,native_recovery_pilots=pilots,future_pair_workload=work,
        counts=dict(core_support_tests=45,PD_supports=7,ordinary_quad_vertices=4,extra_vertices=3,compact_four_simplexes=1,
            source_classes=len(rows),positive_classes=9,marker_domains=len(markers),
            role_orientations=sum(r['active_role_orientations']for r in markers),
            native_pilot_domains=len(pilots),native_pilot_maps=sum(len(r['complete_active_maps'])for r in pilots)),
        current_positive_mass='83284683/10000000000',source_quad_demand=1,extra_occurrences_per_source=3,
        source_model_interval=model['exact_interval'],old_outcomes_changed=False,new_pair_tests=0,exported_columns=[],optimizer_invocations=0,
        scope='Complete raw357 core PD-support and whole spectator inventory; exact actual-vertex/simplex and three bounded full native recovery domains. No compatibility capacity, exported charge, fixed-mass gain, global feasibility or realized chamber claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'],out['counts'],out['future_pair_workload'],flush=True)
