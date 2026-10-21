#!/usr/bin/env python3
"""All remaining raw4 actual-corner demands and exact pair-only benefit bounds.

No new wall pairs, graph cliques or covering LPs are computed. A small exact
cardinality convolution selects subsets of the singleton-cover upper bound.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
import fresh_new5121_raw4_next_four_domain_20261018 as Q

ROOT=Q.ROOT;Z=Q.Z;D=Z.M;R=Q.R;C=Q.C
TAG='fresh_new5121_raw4_all89_cover_barrier_20261018'
normal=Q.normalize


def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def components(rows):
    byv=defaultdict(list);adj={r['face']:set() for r in rows}
    for r in rows:
        for v,n in r['corner_coefficients']:byv[v].append(r['face'])
    for fs in byv.values():
        for a,b in combinations(fs,2):adj[a].add(b);adj[b].add(a)
    left=set(adj);out=[]
    while left:
        stack=[min(left)];seen=set(stack)
        while stack:
            for x in sorted(adj[stack.pop()]-seen):seen.add(x);stack.append(x)
        left-=seen;out.append(sorted(seen))
    return out


def build():
    inputs={}
    def read(tag):
        if tag not in inputs:inputs[tag]=D.read(tag)
        return inputs[tag]
    migration,u,base,l=map(read,(Q.MIGRATION,Q.CURRENT,Q.BASE,Z.L.TAG))
    represented=set(migration['current290_classes'])|{r['selected_face']['face'] for r in u['added_faces']}|{r['face'] for r in base['selected_faces']}
    assert len(represented)==base['represented_classes']==324
    family=sorted((r for r in migration['quad_rows'] if r['face'] not in represented and r.get('raw_family')==4),key=lambda r:(-F(r['demand_mass']),r['face']))
    assert len(family)==89
    alpha={i:F(x) for i,x in base['face_multipliers']};mu={v:F(x) for v,x in base['vertex_multipliers']}
    assert len(alpha)==5121 and base['face_multipliers']==u['face_multipliers']==l['face_multipliers']
    latest={r['vertex']:(Z.L.TAG,r) for r in l['forest_records']}
    latest.update({r['vertex']:(Q.CURRENT,r) for r in u['forest_records']})
    latest.update({r['vertex']:(Q.BASE,r) for r in base['forest_records']});assert len(latest)==600
    for ref in base['inherited_forests']:
        tag,row=latest[ref['vertex']];assert tag==ref['source'] and jhash(row)==ref['record_sha256']
    olddomains=D.current_domains()
    for tag in D.DOMAIN_STEMS:read(tag)
    groups={v:r['groups'] for v,r in olddomains.items()}
    groups.update({r['vertex']:r['domain']['corner_orbits'] for r in u['forest_records']})
    groups.update({r['vertex']:r['domain']['corner_orbits'] for r in base['forest_records']})
    assert set(groups)==set(latest)
    opt={v:F(r['exact_local_optimum']) for v,(tag,r) in latest.items()}
    assert all(x<=mu.get(v,F()) for v,x in opt.items())
    oldn={v:sum(len(g['literal_stars']) for g in gs) for v,gs in groups.items()}
    for v,gs in groups.items():assert oldn[v]==len({i for g in gs for i in g['literal_stars']})

    raw,quads,spectators,forests=map(read,(Q.M.TAGS['raw'],Q.M.TAGS['quad'],Q.B.Q.SPECTATORS,Q.P.FORESTS))
    models=dict(C.models());finite={C.code(g):name for name,g in models.items()}
    lookup={tuple(v['forest']):v['index'] for v in forests['vertices']};names={v['index']:v['forest'] for v in forests['vertices']}
    core=quads['canonical_cores'][4]
    assert (core['index'],core['source_core'],core['nodes'],core['required'])==(4,4,5,3)
    graph=Q.V.raw_graph(raw['records'][4])
    corners=Q.X.quad_features(dict(quadrilateral_black_partitions=core['quad_partitions']))
    assert len(corners)==len({tuple(c['vertex_omitted']) for c in corners})==6
    templates=[]
    for corner in corners:
        omitted=tuple(corner['vertex_omitted']);ids=[i for i in range(5) if i not in omitted]
        h={i:{j:graph[i][j][0] for j in ids if i!=j and graph[i][j]!=Q.A.RIGHT} for i in ids}
        parts=C.pieces(h);typ=[finite[C.code(a)] for a in parts]
        automorphisms={}
        for a,name in zip(parts,typ):
            kk=sorted(a);count=0
            for perm in permutations(kk):
                mp=dict(zip(kk,perm))
                count+=all(a[i].get(j,2)==a[mp[i]].get(mp[j],2) for i in kk for j in kk)
            assert count>0;automorphisms[name]=count
        cov=tuple(tuple((j,graph[x][j]) for j in ids if graph[x][j]!=Q.A.RIGHT) for x in omitted)
        field=Q.A.pair_field(cov,(graph[omitted[0]][omitted[1]],));J,active,blocks=R.projection_data(h)
        heights=[]
        for p in cov:
            z=R.projection(p,p,J,field);assert not any(z[1:]);height=z[0]-2;assert height.sign()>0;heights.append(height)
        proj=R.projection(cov[0],cov[1],J,field);native=field.negative_contact(graph[omitted[0]][omitted[1]])
        assert field.signs_difference(tuple(a-b for a,b in zip(proj,native)),heights[0]*heights[1])==0
        templates.append(dict(omitted=omitted,active_graph=h,finite_components=sorted(typ),
            component_automorphisms=automorphisms,projection_blocks=blocks,contacts=cov,
            squared_timelike_heights=[x.out() for x in heights],native_omitted_product=graph[omitted[0]][omitted[1]],
            pair_field_radicals=[field.tau.out(),field.sigma.out()],projection=R.qout(proj)))

    provenance_data,model=read(D.ACTUAL),read(D.MODEL)
    explicit={f:j for j,f in model['prepared']['provenance_to_face']}
    assert len(explicit)==len(model['prepared']['provenance_to_face'])
    rows=[];allvertices=set()
    for q in family:
        fi=q['face']
        if q['provenance_kind']=='explicit':
            j=q['provenance'];assert explicit[fi]==j
            origin=provenance_data['provenance_masses'][j]
            assert origin['index']==j and origin['merged_face']==fi
            assert origin['source_context']==q['source_context'] and origin['canonical_prism_core']==4 and origin['prism_dimension']==3
        else:
            assert q['provenance_kind']=='unique_saved_context' and q['provenance'] is None and fi not in explicit
        assert q['canonical_core']==4 and q['required']==3 and alpha[fi]==F(q['mass'])>0
        face=forests['faces'][fi];assert face['kind']=='canonical_dotted_quad_aggregate' and face['source_contexts']==[[q['source_context'],4]]
        ctx=spectators['contexts'][q['source_context']];assert quads['source_core_maps'][ctx['core']]['canonical_core']==4
        spectator=ctx['spectator_forest'];assert sum(len(models[n]) for n in spectator)==10
        co=Counter();maps=Counter();actual=[]
        for template in templates:
            vf=tuple(sorted(spectator+template['finite_components']));v=lookup[vf];co[v]+=1;allvertices.add(v)
            available=Counter(vf);required=Counter(template['finite_components']);count=1
            for typ,n in required.items():
                assert available[typ]>=n
                for j in range(n):count*=available[typ]-j
                count*=template['component_automorphisms'][typ]**n
            assert count>0;maps[v]+=count
            actual.append(dict(omitted=template['omitted'],vertex=v,whole_component_embedding_count=count))
        assert sum(co.values())==6 and len(actual)==6
        rows.append(dict(face=fi,source_context=q['source_context'],provenance=q['provenance'],provenance_kind=q['provenance_kind'],raw_family=4,
            spectator_forest=spectator,required=3,mass=str(alpha[fi]),demand_mass=str(3*alpha[fi]),
            actual_corners=actual,corner_coefficients=sorted(co.items()),embedding_count_upper_by_vertex=sorted(maps.items())))

    rowmap={r['face']:r for r in rows}
    def extension(fs):
        e=Counter();emb=Counter()
        for fi in fs:
            for v,n in rowmap[fi]['corner_coefficients']:e[v]+=n*alpha[fi]
            for v,n in rowmap[fi]['embedding_count_upper_by_vertex']:emb[v]+=n
        records=[];gain=F()
        for v,x in sorted(e.items()):
            upper=max(mu.get(v,F()),opt.get(v,F())+x);inc=upper-mu.get(v,F());gain+=inc
            records.append(dict(vertex=v,forest=names[v],actual_corner_demand=str(x),
                baseline_local_cover=str(opt.get(v,F())),baseline_vertex_mass=str(mu.get(v,F())),
                possible_vertex_mass_upper=str(upper),maximum_possible_mu_increase=str(inc),
                baseline_record=None if v not in latest else dict(source=latest[v][0],record_sha256=jhash(latest[v][1])),
                old_literal_stars=oldn.get(v,0),candidate_new_stars_upper=emb[v],
                candidate_new_pair_tests_upper=oldn.get(v,0)*emb[v]+emb[v]*(emb[v]-1)//2))
        return gain,records
    def gain(fs):return extension(fs)[0]
    comps=components(rows);local=[];allsubsets=0;previous={0:(F(),())};stages=[]
    for ci,cc in enumerate(comps):
        bycount={};subsets=[]
        for k in range(len(cc)+1):
            for fs in combinations(cc,k):
                z=gain(fs);subsets.append(dict(faces=fs,gain=str(z)));allsubsets+=1
                if k not in bycount or z>bycount[k][0] or(z==bycount[k][0] and fs<bycount[k][1]):bycount[k]=(z,fs)
        local.append(dict(index=ci,faces=cc,complete_subsets=subsets,
            cardinality_maxima=[dict(count=k,gain=str(z),faces=fs) for k,(z,fs) in sorted(bycount.items())]))
        current={}
        for a,(x,fs) in previous.items():
            for b,(y,gs) in bycount.items():
                z=x+y;combo=tuple(sorted(fs+gs))
                if a+b not in current or z>current[a+b][0] or(z==current[a+b][0] and combo<current[a+b][1]):current[a+b]=(z,combo)
        stages.append([dict(count=k,gain=str(z),faces=fs) for k,(z,fs) in sorted(current.items())]);previous=current
    num,den=F(base['objective_numerator']),F(base['normalizer']);assert den==sum(mu.values(),F())
    needed=num/12-den;fullgain,fullrecords=extension(sorted(rowmap));assert fullgain==previous[89][0]
    crossing=[k for k,(z,fs) in previous.items() if z>needed];minimum=min(crossing) if crossing else None
    assert len(previous)==90 and minimum==42
    bestfaces=previous[minimum][1];bestgain,bestrecords=extension(bestfaces);assert bestgain==previous[minimum][0]
    prefix=[]
    for k in range(1,90):
        fs=[r['face'] for r in rows[:k]];z=gain(fs)
        prefix.append(dict(count=k,gain=str(z),ratio_lower_bound=str(num/(den+z))))
    next4=read(Q.TAG)
    assert gain([r['face'] for r in rows[:4]])==F(next4['cover_inventory']['maximum_possible_mu_increase'])
    sources={}
    for tag,d in inputs.items():
        assert d['status'].startswith('PASS'),(tag,d['status'])
        for p,h in d.get('sources',d.get('prepared',{}).get('sources',{})).items():assert D.digest(p)==h,(tag,p)
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{tag}.{ext}'
            if (ROOT/p).exists():sources[p]=D.digest(p)
    for mod in (Q,Z,D,Q.A,R,C,Q.V,Q.M,Q.X,Q.B,Q.P):
        p=str(Path(mod.__file__).relative_to(ROOT));sources[p]=D.digest(p)
    sources[f'scripts/{TAG}.py']=D.digest(f'scripts/{TAG}.py')
    result=dict(status='PASS_ALL89_RAW4_CORNER_MASS_AND_PAIR_ONLY_CARDINALITY_BARRIER',sources=sources,
        baseline_classes=324,baseline_forests=600,baseline_ratio=base['bound'],face_multipliers=base['face_multipliers'],
        vertex_multipliers=base['vertex_multipliers'],templates=templates,faces=rows,
        counts=dict(classes=89,actual_corners=534,provenance_kinds=dict(Counter(r['provenance_kind'] for r in rows)),target_forests=len(allvertices),
            existing_forests=len(allvertices&set(latest)),additional_forests=len(allvertices-set(latest)),
            independent_face_components=len(comps),largest_face_component=max(map(len,comps)),enumerated_component_subsets=allsubsets),
        objective_numerator=str(num),current_normalizer=str(den),needed_mu_increase_for_ratio12=str(needed),
        all89=dict(forests=fullrecords,maximum_possible_mu_increase=str(fullgain),ratio_lower_bound=str(num/(den+fullgain)),
            remaining_saving_to_restore_barrier=str(fullgain-needed),candidate_pair_tests_upper=sum(r['candidate_new_pair_tests_upper'] for r in fullrecords)),
        component_certificates=local,cardinality_stages=stages,
        cardinality_frontier=[dict(count=k,gain=str(z),faces=fs,ratio_lower_bound=str(num/(den+z)))for k,(z,fs)in sorted(previous.items())],
        minimum_cardinality_permitting_ratio_below12=minimum,
        best_minimum_batch=dict(faces=bestfaces,forests=bestrecords,gain=str(bestgain),ratio_lower_bound=str(num/(den+bestgain)),
            candidate_pair_tests_upper=sum(r['candidate_new_pair_tests_upper']for r in bestrecords)),
        exact_mass_prefix=prefix,new_physical_pair_tests=0,new_clique_enumerations=0,covering_LPs=0,global_model_optimizations=0,
        actual_represented_classes=324,actual_represented_forests=600,new_dimension_bound=None,
        scope='Exact complete provenance and six actual corner demands for all89 remaining unambiguous raw4 classes. Six three-root native projection templates prove singleton feasibility under any orthogonal spectator extension. Frozen current local cover masses and slack give an upper bound on possible monotone repair after new-node pair exclusions only. Small318-subset/cardinality convolution certifies the minimum batch passing that bound, not actual exclusion or realizability. Embedding counts bound new-star/pair workload; no old/native-overlap or clique census is claimed. Current324 graph outcomes, all owned unknowns, face masses and two-owner costs remain unchanged.')
    print(result['status'],result['counts'],'minimum_batch',minimum,'all89_lower',float(num/(den+fullgain)),flush=True)
    return normal(result)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';r=build()
    if args.verify:assert r==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(r,separators=(',',':'))+'\n')
