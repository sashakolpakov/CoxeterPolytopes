#!/usr/bin/env python3
"""Exact joint actual-corner capacities for raw87/88/89 at fixed new masses.

External contacts lie in K(sqrt(tau87),sqrt(tau88)); tau89 has a K root.
Complete order8 cores and spectator rank7 are retained in one vertex frame.
No optimizer. Unknown pairs remain compatible.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from functools import lru_cache
from pathlib import Path
import json,time
import networkx as nx
import fresh_dominant_cap_compatibility_20261017 as P
import fresh_mixed_raw24_actual_corners_20261017 as T
ROOT=P.ROOT
TAG='fresh_raw87_88_89_actual_corners_20261017'
FULL='fresh_raw178_296_297_full_model_20261017'
SUPPORT=FULL+'_support'
M,V,K,C,X,B=P.M,P.V,P.K,P.C,P.X,P.B
RIGHT=V.RIGHT
DIMENSIONS={87:4,88:4,89:4}
CORE_IDS={87:76,88:77,89:78}
CORE_ORDERS={87:8,88:8,89:8}
TAUS={87:K(F(9,2),F(3,2)),88:K(5,2),89:K(6,2)}
ROOTS={89:K(1,1)}
FIELD=T.PairField(TAUS[87],TAUS[88])
active_embeddings=T.active_embeddings
projection_data=T.projection_data
sparse_covector=T.sparse_covector
negative_contact=FIELD.negative_contact
signs_difference=FIELD.signs_difference
qaddk=T.qaddk
qout=T.qout

def projection(p,q,J):
    return T.projection(p,q,J,FIELD)

def run(domain_only=False):
    start=time.monotonic()
    full=P.read(FULL);support=P.read(SUPPORT)
    witness=next(c for c in full['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    actual=P.read(M.TAGS['actual']);quads=P.read(M.TAGS['quad'])
    raw=P.read(M.TAGS['raw']);spectators=P.read(B.Q.SPECTATORS)
    old=P.read(P.FORESTS);vertices=old['vertices'];oldfaces=old['faces']
    lookup={tuple(v['forest']):v['index']for v in vertices}
    alpha={i:F(a)for i,a in witness['face_multipliers']}
    mu={v:F(a)for v,a in witness['vertex_multipliers']}
    assert len(alpha)==4410
    assert all(a>0 for a in alpha.values())and all(a>=0 for a in mu.values())
    assert sum(mu.values())==F(witness['normalizer'])
    assert F(witness['objective_numerator'])/sum(mu.values())==F(witness['bound'])
    assert M.ksqrt(TAUS[89])==ROOTS[89] and ROOTS[89].sign()>0
    assert all(M.ksqrt(t)is None for t in (TAUS[87],TAUS[88],TAUS[87]/TAUS[88]))
    assert witness['face_multipliers']==support['face_multipliers']
    assert witness['vertex_multipliers']==support['vertex_multipliers']
    pmap=dict(full['prepared']['provenance_to_face'])
    assert len(set(pmap.values()))==len(pmap)
    inverse={fi:j for j,fi in pmap.items()}
    models=dict(C.models());finite={C.code(g):n for n,g in models.items()}
    selected=[];explicit=unique=0;ambiguous_selected=[]
    for fi in sorted(alpha):
        if fi in inverse:
            j=inverse[fi];pr=dict(actual['provenance_masses'][j])
            cc=quads['canonical_cores'][pr['canonical_prism_core']]
            if cc['source_core'] not in DIMENSIONS:continue
            assert pr['prism_dimension']==DIMENSIONS[cc['source_core']]
            explicit+=1;pr['provenance_kind']='existing_explicit_map'
        elif fi<len(oldfaces) and oldfaces[fi]['kind']=='canonical_dotted_quad_aggregate':
            source=oldfaces[fi]
            # No arbitrary representative of a merged class is permitted.
            # Select only if EVERY saved realization is this single context.
            if len(source['source_contexts'])!=1:
                if any(quads['canonical_cores'][ci]['source_core'] in DIMENSIONS for ctx,ci in source['source_contexts']):
                    ambiguous_selected.append(fi)
                continue
            ctx,ci=source['source_contexts'][0]
            cc=quads['canonical_cores'][ci]
            if cc['source_core'] not in DIMENSIONS:continue
            assert source['required']==cc['required']
            assert quads['source_core_maps'][spectators['contexts'][ctx]['core']]['canonical_core']==ci
            j=None;unique+=1
            pr=dict(source_context=ctx,canonical_prism_core=ci,prism_dimension=DIMENSIONS[cc['source_core']],
                provenance_kind='unique_saved_context',source_contexts=source['source_contexts'])
        else:continue
        assert cc['index']==CORE_IDS[cc['source_core']]
        assert cc['required']==DIMENSIONS[cc['source_core']]*(DIMENSIONS[cc['source_core']]-1)//2
        ctx=spectators['contexts'][pr['source_context']]
        assert quads['source_core_maps'][ctx['core']]['canonical_core']==cc['index']
        assert sum(len(models[name])for name in ctx['spectator_forest'])==15-cc['nodes']
        selected.append((fi,j,pr,cc))
    occurrences=defaultdict(list);face_records=[];target_graphs={}
    for fi,j,pr,cc in selected:
        family=cc['source_core'];dimension=DIMENSIONS[family]
        assert pr['prism_dimension']==dimension
        rg=V.raw_graph(raw['records'][family]);assert len(rg)==CORE_ORDERS[family]==cc['nodes']
        assert all(V.squared(e)==TAUS[family]for row in rg for e in row if e[0]==0)
        assert all(e[0]in(0,2,3,5)for row in rg for e in row)
        sp=spectators['contexts'][pr['source_context']]['spectator_forest']
        corners=X.quad_features(dict(quadrilateral_black_partitions=cc['quad_partitions']))
        assert len(corners)==2*dimension
        assert len({tuple(sorted(c['vertex_omitted']))for c in corners})==2*dimension
        records=[]
        for c in corners:
            omitted=tuple(c['vertex_omitted']);ids=[i for i in range(len(rg))if i not in omitted]
            assert len(ids)==len(rg)-2
            h={i:{j:rg[i][j][0]for j in ids if j!=i and rg[i][j]!=RIGHT}for i in ids}
            assert all(m in(3,5)for row in h.values()for m in row.values())
            vf=tuple(sorted(sp+[finite[C.code(a)]for a in C.pieces(h)]));vi=lookup[vf]
            if vi not in target_graphs:target_graphs[vi]=P.forest_graph(vf,models)
            target=target_graphs[vi]
            assert len(target)==13
            maps=active_embeddings(h,target)
            orientations={}
            for mp in maps:
                es=[]
                for x in omitted:
                    e=[RIGHT]*13
                    for i,z in mp.items():e[z]=rg[x][i]
                    es.append(tuple(e))
                key=tuple(sorted(es))
                assert key[0]!=key[1]
                occ=tuple(sorted(tuple(sorted(mp[z]for z in pair))for pair in c['occupied_pairs']))
                assert len(occ)==dimension-1 and len(set(occ))==dimension-1
                if key in orientations:assert orientations[key][0]==occ
                else:orientations[key]=(occ,sorted(mp.items()))
            record=dict(face=fi,provenance=j,raw_family=family,omitted=omitted,vertex=vi,
                        active_embeddings=len(maps),literal_orientations=len(orientations))
            records.append(record)
            occurrences[vi].append((record,orientations,rg[omitted[0]][omitted[1]]))
        face_records.append(dict(face=fi,provenance=j,raw_family=family,canonical_prism_core=pr['canonical_prism_core'],
            source_context=pr['source_context'],provenance_kind=pr['provenance_kind'],
            source_contexts=pr.get('source_contexts'),raw_core_order=len(rg),
            spectator_rank=sum(len(models[name])for name in sp),spectator_forest=sp,corners=records))
    print('DOMAIN',len(selected),'faces',sum(map(len,occurrences.values())),'corners',
          len(occurrences),'forests',sum(r[0]['literal_orientations']for rs in occurrences.values()for r in rs),
          'orientation occurrences',flush=True)
    family_counts={}
    for family in DIMENSIONS:
        rr=[r for rs in occurrences.values()for r in rs if r[0]['raw_family']==family]
        family_counts[family]=dict(classes=sum(cc['source_core']==family for fi,j,pr,cc in selected),
            corners=len(rr),forests=len({r[0]['vertex']for r in rr}),orientations=sum(r[0]['literal_orientations']for r in rr))
    common_forests=sum(len({r[0]['raw_family']for r in rs})>1 for rs in occurrences.values())
    expected_cross_pairs=sum(sum(a[0]['literal_orientations']*b[0]['literal_orientations']
        for a,b in combinations(rs,2)if a[0]['raw_family']!=b[0]['raw_family'])for rs in occurrences.values())
    expected={87:(36,288),88:(27,216),89:(26,208)}
    assert {f:(r['classes'],r['corners'])for f,r in family_counts.items()}==expected
    supported={r['face']for r in support['quad_rows']if r.get('raw_family')in DIMENSIONS and r['provenance_kind']!='merged_unexpanded'}
    assert {r['face']for r in face_records}==supported and len(selected)==89
    assert (explicit,unique)==(51,38)
    print('FAMILIES',family_counts,'common forests',common_forests,'cross pairs',expected_cross_pairs,flush=True)
    if domain_only:
        print('PASS_COMPLETE_RAW87_RAW88_RAW89_DOMAIN',len(selected),sum(map(len,occurrences.values())),flush=True)
        return None
    face_family={fi:cc['source_core']for fi,j,pr,cc in selected}
    forest_records=[];columns=[];updates=[];single_updates=[];counts=Counter()
    for step,(vi,orecs)in enumerate(sorted(occurrences.items()),1):
        target=target_graphs[vi];J,active,components=projection_data(target)
        stars={};orbits={}
        for rec,oris,native in orecs:
            orbit=tuple(sorted(oris))
            if orbit not in orbits:orbits[orbit]=dict(origins=[],coefficients=Counter())
            orbits[orbit]['origins'].append(rec);orbits[orbit]['coefficients'][rec['face']]+=1
            for key,(occupied,mp)in oris.items():
                if key in stars:
                    assert stars[key]['occupied_pairs']==occupied
                    assert stars[key]['native']==native
                else:
                    p=[sparse_covector(e)for e in key]
                    assert all(i in active for q in p for i,z in q)
                    heights=[]
                    for q in p:
                        val=projection(q,q,J)
                        assert not any(val[1:])
                        height=val[0]-2;assert height.sign()>0;heights.append(height)
                    val=projection(p[0],p[1],J)
                    n=negative_contact(native)
                    # All external roots use the same future sign, agreeing
                    # with their prescribed normalized product in the core.
                    assert signs_difference(tuple(a-b for a,b in zip(val,n)),heights[0]*heights[1])==0
                    stars[key]=dict(occupied_pairs=occupied,embedding=mp,covectors=p,
                        heights=heights,native=native,origin=rec)
        # Distinct complete automorphism orbits must be disjoint. Overlap
        # would expose an incomplete map domain or duplicate class coding.
        orbit_keys=list(sorted(orbits))
        for a,b in combinations(orbit_keys,2):assert not(set(a)&set(b))
        for o in orbits.values():assert len(o['coefficients'])==1
        keys=sorted(stars);ix={key:i for i,key in enumerate(keys)}
        literal=[]
        for key in keys:
            s=stars[key]
            literal.append(dict(index=ix[key],contacts=[[[i,list(e)]for i,e in enumerate(q)if e!=RIGHT]for q in key],
                occupied_pairs=s['occupied_pairs'],active_embedding=s['embedding'],origin=s['origin'],
                squared_timelike_heights=[h.out()for h in s['heights']],native_omitted_product=list(s['native'])))
        pair_cache={};pair_cert=[];expressions=[];expression_ids={};reason_counts=Counter();wall_cache={}
        cross_tested=cross_excluded=0
        def incompatible(a,b):
            nonlocal cross_tested,cross_excluded
            if a==b:return False
            ij=tuple(sorted((ix[a],ix[b])))
            if ij in pair_cache:return pair_cache[ij]
            sa,sb=stars[a],stars[b]
            cross=sa['origin']['raw_family']!=sb['origin']['raw_family']
            cross_tested+=int(cross)
            shared=sorted(set(sa['occupied_pairs'])&set(sb['occupied_pairs']))
            reason=None
            if shared:
                reason=[*ij,0];reason_counts['distinct_full_stars_share_actual_quad_port_pair']+=1
            else:
                for u,v in ((0,0),(0,1),(1,0),(1,1)):
                    if a[u]==b[v]:continue # identical actual wall, not a forbidden pair
                    wk=tuple(sorted((a[u],b[v])))
                    if wk not in wall_cache:
                        p=projection(sa['covectors'][u],sb['covectors'][v],J)
                        h=sa['heights'][u]*sb['heights'][v]
                        sg=signs_difference(p,h);sg1=signs_difference(qaddk(p,1),h)
                        code=1 if sg>0 else 2 if sg<0 and sg1>0 else None
                        expr=dict(pair_field_radicals=[FIELD.tau.out(),FIELD.sigma.out()],projection=qout(p),height_product=h.out(),
                            sign_product=sg,sign_product_plus_one=sg1)
                        if code is None:
                            code,extra=FIELD.arithmetic_wall_gate(tuple((z.a,z.b)for z in p),h.a,h.b)
                            if code is not None:expr['arithmetic_gate']=extra
                        wall_cache[wk]=(code,expr)
                    code,expr=wall_cache[wk]
                    if code is not None:
                        ek=P.dump(expr)
                        if ek not in expression_ids:
                            expression_ids[ek]=len(expressions);expressions.append(expr)
                        reason=[*ij,code,u,v,expression_ids[ek]]
                        reason_counts[{1:'positive_actual_wall_product',2:'actual_wall_product_in_open_minus1_zero',3:'nonintegral_native_two_cycle',4:'native_two_cycle_outside_fixed_K',5:'forbidden_finite_Coxeter_product'}[code]]+=1
                        break
            pair_cache[ij]=reason is not None
            if reason is not None:
                pair_cert.append(reason);cross_excluded+=int(cross)
            return reason is not None
        groups=[];graph=nx.Graph()
        for oi,key in enumerate(orbit_keys):
            oo=orbits[key];internal=all(incompatible(a,b)for a,b in combinations(key,2))
            group=dict(index=oi,literal_stars=[ix[k]for k in key],positive_face_coefficients=sorted(oo['coefficients'].items()),
                origins=oo['origins'],capacity_one_proved=internal)
            groups.append(group)
            if internal:graph.add_node(oi)
        for a,b in combinations(sorted(graph),2):
            if all(incompatible(x,y)for x in orbit_keys[a]for y in orbit_keys[b]):graph.add_edge(a,b)
        local_columns={}
        for clique in nx.find_cliques(graph):
            clique=sorted(clique);coeff=Counter()
            for g in clique:coeff.update(orbits[orbit_keys[g]]['coefficients'])
            key=tuple(sorted(coeff.items()))
            if key not in local_columns:local_columns[key]=clique
        # Actual simultaneous stars form a clique in the conservative
        # compatibility graph: every unexcluded pair is retained as an edge.
        # A proper coloring therefore gives an exact upper capacity, without
        # needing a maximum-clique or optimization claim.
        compatibility=nx.Graph();compatibility.add_nodes_from(range(len(keys)))
        for ia,ib in combinations(range(len(keys)),2):
            if not incompatible(keys[ia],keys[ib]):compatibility.add_edge(ia,ib)
        assert len(pair_cache)==len(keys)*(len(keys)-1)//2
        candidates={}
        for coeff,clique in local_columns.items():
            candidates[coeff]=dict(cost=1,orbit_clique=clique,kind='capacity_one_orbit_clique')
        subsets=[(oi,)for oi in range(len(orbit_keys))]
        if len(orbit_keys)>1:subsets.append(tuple(range(len(orbit_keys))))
        for subset in subsets:
            nodes=sorted(ix[k]for oi in subset for k in orbit_keys[oi])
            assert len(nodes)==len(set(nodes))
            subgraph=compatibility.subgraph(nodes)
            partitions=[]
            for strategy in ('saturation_largest_first','largest_first','smallest_last'):
                coloring=nx.coloring.greedy_color(subgraph,strategy=strategy)
                parts=tuple(sorted(tuple(sorted(i for i in nodes if coloring[i]==color))
                    for color in set(coloring.values())))
                assert sorted(i for part in parts for i in part)==nodes
                assert all(incompatible(keys[i],keys[j])for part in parts for i,j in combinations(part,2))
                partitions.append(parts)
            parts=min(partitions,key=lambda x:(len(x),x))
            coeff=Counter()
            for oi in subset:coeff.update(orbits[orbit_keys[oi]]['coefficients'])
            coeff=tuple(sorted(coeff.items()));cap=len(parts)
            if coeff not in candidates or cap<candidates[coeff]['cost']:
                candidates[coeff]=dict(cost=cap,kind='proper_compatibility_coloring',
                    orbit_indices=list(subset),literal_coloring=parts)
        best=best_single=mu.get(vi,F());local_ids=[]
        for coeff,cert in sorted(candidates.items()):
            load=sum((m*alpha[fi]for fi,m in coeff),F());idx=len(columns);cost=cert.pop('cost')
            columns.append(dict(index=idx,vertex=vi,vertex_cost=cost,positive_face_coefficients=coeff,
                raw_families=sorted({face_family[fi]for fi,m in coeff}),**cert,load=str(load),baseline_vertex_mass=str(mu.get(vi,F())),
                violated=load>cost*mu.get(vi,F())))
            local_ids.append(idx);best=max(best,load/cost)
            if len(columns[-1]['raw_families'])==1:best_single=max(best_single,load/cost)
        if best>mu.get(vi,F()):updates.append([vi,str(mu.get(vi,F())),str(best)])
        if best_single>mu.get(vi,F()):single_updates.append([vi,str(mu.get(vi,F())),str(best_single)])
        forest_records.append(dict(vertex=vi,forest=vertices[vi]['forest'],active_components=components,
            corner_occurrences=len(orecs),literal_stars=literal,corner_orbits=groups,
            incompatibility_certificates=sorted(pair_cert),radical_expressions=expressions,
            group_incompatibility_edges=sorted(tuple(sorted(e))for e in graph.edges),columns=local_ids,
            tested_literal_pairs=len(pair_cache),unknown_literal_pairs=sum(not b for b in pair_cache.values()),
            cross_family_pairs=cross_tested,cross_family_exclusions=cross_excluded))
        counts.update(dict(literal_orientations=len(keys),corner_orbits=len(groups),
            capacity_one_orbits=len(graph),group_edges=graph.number_of_edges(),literal_pairs_tested=len(pair_cache),
            exact_literal_exclusions=len(pair_cert),cross_family_pairs=cross_tested,cross_family_exclusions=cross_excluded,columns=len(candidates)))
        counts.update(reason_counts)
        if step%25==0 or step==len(occurrences):
            print('FORESTS',step,'/',len(occurrences),'stars',counts['literal_orientations'],
                  'columns',len(columns),'violated forests',len(updates),flush=True)
    new_mu=dict(mu)
    for v,old,new in updates:assert new_mu.get(v,F())==F(old);new_mu[v]=F(new)
    for col in columns:
        assert sum((m*alpha[fi]for fi,m in col['positive_face_coefficients']),F())<=col['vertex_cost']*new_mu[col['vertex']]
    assert all(new_mu[v]>=x for v,x in mu.items())
    single_mu=dict(mu)
    for v,old,new in single_updates:single_mu[v]=F(new)
    normalizer=sum(new_mu.values());objective=F(witness['objective_numerator']);ratio=objective/normalizer
    counts.update(actual_Q_classes=len(selected),existing_explicit_classes=explicit,new_unique_context_classes=unique,corner_occurrences=sum(map(len,occurrences.values())),
        forests=len(occurrences),violated_columns=sum(c['violated']for c in columns),vertex_mass_increases=len(updates))
    counts.update(mixed_columns=sum(len(c['raw_families'])>1 for c in columns),
        violated_mixed_columns=sum(c['violated']and len(c['raw_families'])>1 for c in columns),
        mixed_additional_vertex_raises=sum(new_mu.get(v,F())>single_mu.get(v,F())for v in new_mu))
    assert counts['cross_family_pairs']==expected_cross_pairs
    counts['cross_family_unknown']=expected_cross_pairs-counts['cross_family_exclusions']
    paths=[Path(__file__),Path(P.__file__),Path(M.__file__),Path(V.__file__),Path(C.__file__),Path(X.__file__),Path(B.__file__),Path(B.Q.__file__)]
    paths +=[ROOT/f'results/{t}.json'for t in(FULL,SUPPORT,P.FORESTS,M.TAGS['actual'],M.TAGS['quad'],M.TAGS['raw'],B.Q.SPECTATORS)]
    paths +=[ROOT/f'scripts/{FULL}.py',ROOT/f'results/{FULL}.md',ROOT/f'scripts/{SUPPORT}.py',ROOT/f'results/{SUPPORT}.md',Path(T.__file__)]
    paths +=[ROOT/'scripts/fresh_mixed_raw39_actual_corners_20261017.py',ROOT/'results/fresh_mixed_raw39_actual_corners_20261017.md']
    paths +=[ROOT/'scripts/fresh_sqrt5_labelled_packing_full_model_20261017.py',ROOT/'scripts/fresh_actual_corner_enlarged_model_20261017.py']
    paths +=[ROOT/'results/fresh_quad_product_aggregation_20261017.md',ROOT/f'results/{P.TAG}.md',ROOT/'results/fresh_sqrt5_spectator_groundfield_gate_20261017.md']
    result=dict(status='PASS_EXACT_RAW87_RAW88_RAW89_JOINT_FINITE_VERTEX_CORNER_CAPACITIES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        raw_sources=sorted(DIMENSIONS),canonical_cores=list(CORE_IDS.values()),raw_core_orders=CORE_ORDERS,
        dotted_squares={f:t.out()for f,t in TAUS.items()},positive_dotted_roots={f:t.out()for f,t in ROOTS.items()},contact_field="K(sqrt((9+3sqrt5)/2),sqrt(5+2sqrt5)), K=Q(sqrt5)",
        pair_field_basis_squares=[z.out()for z in FIELD.bases],fixed_ambient_ground_field="Q(sqrt5)",
        family_counts=family_counts,common_forests=common_forests,
        omitted_ambiguous_selected_rows=ambiguous_selected,selected_faces=face_records,forest_records=forest_records,
        columns=columns,counts=dict(counts),
        fixed_witness_repair=dict(face_masses_unchanged=True,baseline=FULL,vertex_updates=updates,
            objective_numerator=str(objective),normalizer=str(normalizer),bound=str(ratio),bound_below12=ratio<12,
            baseline_bound=witness['bound'],optimizer_invocations=0,
            without_mixed_family_updates=single_updates,without_mixed_family_bound=str(objective/sum(single_mu.values())),
            preservation='All face masses fixed and only nonnegative vertex allowances increase; every436262 exported inherited column, including96176 zero-budget and52 positive-cost signed columns, remains valid by monotonicity. Old complete implicit covers at previous masses are not claimed.'),
        certificate_schema=dict(pair_reason_codes={'0':'shared occupied square pair','1':'positive actual-wall product','2':'actual-wall product in (-1,0)','3':'nonintegral native two-cycle','4':'native two-cycle outside fixed K','5':'forbidden finite Coxeter product'},pair_records='[star_a,star_b,0] or [star_a,star_b,reason,wall_a,wall_b,radical_expression_index]'),
        inherited_dual_convention='Every active primal charge variable has bounds(0,None). Signed denotes face coefficients, not a free variable. Its dual condition is load_j<=cost_j*mu_v, with positive vertex cost, or load_j<=0 for a zero-budget transfer.',
        face_multipliers=witness['face_multipliers'],initial_vertex_multipliers=witness['vertex_multipliers'],
        vertex_multipliers=[[v,str(a)]for v,a in sorted(new_mu.items())],
        scope='The positive raw87/raw88/raw89 classes with explicit or uniquely saved provenance in the stated fixed witness, all their finite prism corners in common complete vertex frames and exact active orientation orbits. Exact pair-exclusion cliques and proper colorings of the conservative compatibility graph give the stated capacities. Every unknown pair remains allowed; coloring need not be optimal. No optimizer is invoked. The implicit full local packing/joint families have not been rebuilt at these4410 face masses; no complete-family coverage, target12 feasibility or dimension theorem is asserted.')
    print('COMPLETE',dict(counts),'bound',float(ratio),'seconds',round(time.monotonic()-start,2),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--domain-only',action='store_true');args=ap.parse_args()
    result=run(args.domain_only)
    if result is None:raise SystemExit(0)
    d=json.loads(json.dumps(result,sort_keys=True,separators=(',',':')));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n')
    print(d['status'],flush=True)
