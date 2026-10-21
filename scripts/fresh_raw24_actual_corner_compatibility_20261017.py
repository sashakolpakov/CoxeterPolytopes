#!/usr/bin/env python3
"""Exact finite-vertex corner-orbit clique/coloring capacities for 45 positive raw24 prism classes.

No optimizer. Unknown pairs are left compatible. Inactive spectator
automorphisms are quotiented only after retaining every active embedding.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from functools import lru_cache
import json
from pathlib import Path
import time
import networkx as nx
import fresh_dominant_cap_compatibility_20261017 as P

ROOT=P.ROOT
TAG='fresh_raw24_actual_corner_compatibility_20261017'
WITNESS='fresh_mixed_raw39_full_model_20261017_joint_ports'
FULL='fresh_mixed_raw39_full_model_20261017'
M,V,K,C,X,B=P.M,P.V,P.K,P.C,P.X,P.B
TAU=K(F(23,2),F(9,2))
RIGHT=V.RIGHT


def qout(q):return [q[0].out(),q[1].out()]
def qm(q,r):return M.mulrad(q,r,TAU)
def qa(q,r):return (q[0]+r[0],q[1]+r[1])
def negative_contact(e):
    if e==RIGHT:return K(),K()
    if e[0]==0:
        assert V.squared(e)==TAU
        return K(),K(-1)
    assert e[0] in (3,5)
    return -M.ksqrt(V.squared(e)),K()


def signs_difference(a,b,h):
    """Sign of a+b*sqrt(TAU)-sqrt(h), with both square roots positive."""
    assert h.sign()>0
    s=M.signrad(a,b,TAU)
    if s<=0:return -1
    return M.signrad(a*a+b*b*TAU-h,2*a*b,TAU)


@lru_cache(None)
def arithmetic_wall_gate(xa,xb,ya,yb,ha,hb):
    """Exact fixed-K membership of (a+b sqrt(TAU)-sqrt(h))^2.

    TAU is nonsquare. The three branches are the two quadratic cases
    and the independent biquadratic case; no floating field decision.
    """
    a,b,h=K(xa,xb),K(ya,yb),K(ha,hb)
    r=M.ksqrt(h);r_tau=M.ksqrt(h/TAU)
    cert=dict(height_product=h.out())
    if r is not None:
        assert r.sign()>0
        A,B=a-r,b;cert.update(branch='height_square_in_K',root=r.out())
    elif r_tau is not None:
        assert r_tau.sign()>0
        A,B=a,b-r_tau;cert.update(branch='height_over_tau_square_in_K',root=r_tau.out())
    else:
        cert.update(branch='independent_biquadratic')
        if a or b:
            cert['two_cycle_in_fixed_K']=False
            return 4,cert
        A,B=None,None;square=h
    if A is not None:
        cert['product_coefficients']=[A.out(),B.out()]
        if A and B:
            cert['two_cycle_in_fixed_K']=False
            return 4,cert
        square=A*A+B*B*TAU
    cert.update(two_cycle_in_fixed_K=True,cross_two_cycle=square.out())
    if not M.integral(square):
        cert['two_cycle_integral']=False
        return 3,cert
    # Complete finite labels with squared product in K are 2,3,4,5,6,10.
    allowed=(K(),K(1),K(2),K(F(3,2),F(1,2)),K(3),K(F(5,2),F(1,2)))
    if signs_difference(a,b,h)<0 and signs_difference(a+2,b,h)>0 and square not in allowed:
        cert['finite_angle_in_alphabet']=False
        return 5,cert
    return None,None


def active_embeddings(raw_corner,target):
    """All isomorphisms of each whole raw component into a whole target one.

    An inactive component permutation has zero external contacts and no
    occupied pair; it changes none of the exported literal corner data.
    """
    left=C.pieces(raw_corner);right=C.pieces(target)
    domains=[]
    for a in left:
        choices=[]
        for j,b in enumerate(right):
            if C.code(a)!=C.code(b):continue
            gm=nx.algorithms.isomorphism.GraphMatcher(P.ngraph(a),P.ngraph(b),
                edge_match=lambda x,y:x['label']==y['label'])
            maps=list(gm.isomorphisms_iter());assert maps
            choices += [(j,mp)for mp in maps]
        assert choices
        domains.append(choices)
    out=[]
    def visit(k,used,mp):
        if k==len(domains):out.append(dict(mp));return
        for j,part in domains[k]:
            if j not in used:visit(k+1,used|{j},mp|part)
    visit(0,set(),{})
    assert len({tuple(sorted(x.items()))for x in out})==len(out)
    for mp in out:
        assert len(mp)==len(raw_corner) and len(set(mp.values()))==len(mp)
        assert all(raw_corner[i].get(j,2)==target[mp[i]].get(mp[j],2)
                   for i in raw_corner for j in raw_corner)
        # Every mapped component is complete, so unused target roots have
        # zero incident products with all active roots.
        assert all(not (set(target[j])-set(mp.values()))for j in mp.values())
    return out


def projection_data(target):
    """Only 3/5 components can carry raw24 covectors.

    Passive B/F/I2(6)/I2(10) components need no normalized K inverse:
    they are orthogonal to every active external covector.
    """
    J={};active=set();components=[]
    for c in C.pieces(target):
        ids=sorted(c)
        if any(m not in (3,5)for row in c.values()for m in row.values()):continue
        G=[[K(2)if i==j else -M.ksqrt(V.squared(V.edge(c[i][j])))
            if j in c[i]else K()for j in ids]for i in ids]
        assert V.pd(G)
        Ji,det=V.inverse(G);assert det.sign()>0
        for a,i in enumerate(ids):
            for b,j in enumerate(ids):J[i,j]=Ji[a][b]
        active.update(ids)
        components.append(dict(roots=ids,gram=[[x.out()for x in r]for r in G],
            inverse=[[x.out()for x in r]for r in Ji]))
    return J,active,components


def sparse_covector(es):
    return [(i,negative_contact(e))for i,e in enumerate(es)if e!=RIGHT]


def projection(p,q,J):
    z=(K(),K())
    for i,a in p:
        for j,b in q:
            if (i,j)in J:
                x=qm(a,b);z=qa(z,(J[i,j]*x[0],J[i,j]*x[1]))
    return z


def run():
    start=time.monotonic()
    witness=P.read(WITNESS);full=P.read(FULL)
    assert M.ksqrt(TAU) is None and TAU.sign()>0
    actual=P.read(M.TAGS['actual']);quads=P.read(M.TAGS['quad'])
    raw=P.read(M.TAGS['raw']);spectators=P.read(B.Q.SPECTATORS)
    old=P.read(P.FORESTS);vertices=old['vertices'];oldfaces=old['faces']
    lookup={tuple(v['forest']):v['index']for v in vertices}
    alpha={i:F(a)for i,a in witness['face_multipliers']}
    mu={v:F(a)for v,a in witness['vertex_multipliers']}
    assert all(a>0 for a in alpha.values())and all(a>=0 for a in mu.values())
    assert sum(mu.values())==F(witness['normalizer'])
    assert F(witness['objective_numerator'])/sum(mu.values())==F(witness['bound'])
    pmap=dict(full['prepared']['provenance_to_face'])
    assert len(set(pmap.values()))==len(pmap)
    inverse={fi:j for j,fi in pmap.items()}
    models=dict(C.models());finite={C.code(g):n for n,g in models.items()}
    selected=[];explicit=unique=0;ambiguous_raw24=[]
    for fi in sorted(alpha):
        if fi in inverse:
            j=inverse[fi];pr=dict(actual['provenance_masses'][j])
            cc=quads['canonical_cores'][pr['canonical_prism_core']]
            if cc['source_core']!=24:continue
            assert pr['prism_dimension']==4
            explicit+=1;pr['provenance_kind']='existing_explicit_map'
        elif fi<len(oldfaces) and oldfaces[fi]['kind']=='canonical_dotted_quad_aggregate':
            source=oldfaces[fi]
            # No arbitrary representative of a merged class is permitted.
            # Select only if EVERY saved realization is this single context.
            if len(source['source_contexts'])!=1:
                if any(quads['canonical_cores'][ci]['source_core']==24 for ctx,ci in source['source_contexts']):
                    ambiguous_raw24.append(fi)
                continue
            ctx,ci=source['source_contexts'][0]
            cc=quads['canonical_cores'][ci]
            if cc['source_core']!=24:continue
            assert source['kind']=='canonical_dotted_quad_aggregate'and source['required']==6
            assert quads['source_core_maps'][spectators['contexts'][ctx]['core']]['canonical_core']==ci
            j=None;unique+=1
            pr=dict(source_context=ctx,canonical_prism_core=ci,prism_dimension=4,
                provenance_kind='unique_saved_context',source_contexts=source['source_contexts'])
        else:continue
        assert cc['index']==20 and cc['required']==6
        ctx=spectators['contexts'][pr['source_context']]
        assert quads['source_core_maps'][ctx['core']]['canonical_core']==cc['index']
        assert sum(len(models[name])for name in ctx['spectator_forest'])==8
        selected.append((fi,j,pr,cc))
    assert (explicit,unique,len(selected))==(15,30,45)
    rg=V.raw_graph(raw['records'][24]);assert len(rg)==7
    assert all(e[0]in(0,2,3,5)for row in rg for e in row)
    occurrences=defaultdict(list);face_records=[];target_graphs={}
    for fi,j,pr,cc in selected:
        assert pr['prism_dimension']==4 and cc['source_core']==24
        sp=spectators['contexts'][pr['source_context']]['spectator_forest']
        corners=X.quad_features(dict(quadrilateral_black_partitions=cc['quad_partitions']))
        assert len(corners)==8
        assert len({tuple(sorted(c['vertex_omitted']))for c in corners})==8
        records=[]
        for c in corners:
            omitted=tuple(c['vertex_omitted']);ids=[i for i in range(7)if i not in omitted]
            assert len(ids)==5
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
                assert len(occ)==3 and len(set(occ))==3
                if key in orientations:assert orientations[key][0]==occ
                else:orientations[key]=(occ,sorted(mp.items()))
            record=dict(face=fi,provenance=j,omitted=omitted,vertex=vi,
                        active_embeddings=len(maps),literal_orientations=len(orientations))
            records.append(record)
            occurrences[vi].append((record,orientations,rg[omitted[0]][omitted[1]]))
        face_records.append(dict(face=fi,provenance=j,canonical_prism_core=pr['canonical_prism_core'],
            source_context=pr['source_context'],provenance_kind=pr['provenance_kind'],
            source_contexts=pr.get('source_contexts'),spectator_forest=sp,corners=records))
    print('DOMAIN',len(selected),'faces',sum(map(len,occurrences.values())),'corners',
          len(occurrences),'forests',sum(r[0]['literal_orientations']for rs in occurrences.values()for r in rs),
          'orientation occurrences',flush=True)
    forest_records=[];columns=[];updates=[];counts=Counter()
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
                        a,b=projection(q,q,J);assert not b
                        height=a-2;assert height.sign()>0;heights.append(height)
                    a,b=projection(p[0],p[1],J)
                    n=negative_contact(native)
                    # All external roots use the same future sign, agreeing
                    # with their prescribed normalized product in the core.
                    assert signs_difference(a-n[0],b-n[1],heights[0]*heights[1])==0
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
        pair_cache={};pair_cert=[];expressions=[];expression_ids={};reason_counts=Counter()
        def incompatible(a,b):
            if a==b:return False
            ij=tuple(sorted((ix[a],ix[b])))
            if ij in pair_cache:return pair_cache[ij]
            sa,sb=stars[a],stars[b];shared=sorted(set(sa['occupied_pairs'])&set(sb['occupied_pairs']))
            reason=None
            if shared:
                reason=[*ij,0];reason_counts['distinct_full_stars_share_actual_quad_port_pair']+=1
            else:
                for u,v in ((0,0),(0,1),(1,0),(1,1)):
                    if a[u]==b[v]:continue # identical actual wall, not a forbidden pair
                    x,y=projection(sa['covectors'][u],sb['covectors'][v],J)
                    h=sa['heights'][u]*sb['heights'][v]
                    sg=signs_difference(x,y,h)
                    sg1=signs_difference(x+1,y,h)
                    code=1 if sg>0 else 2 if sg<0 and sg1>0 else None
                    expr=dict(projection=qout((x,y)),height_product=h.out(),
                        sign_product=sg,sign_product_plus_one=sg1)
                    if code is None:
                        code,extra=arithmetic_wall_gate(x.a,x.b,y.a,y.b,h.a,h.b)
                        if code is not None:expr['arithmetic_gate']=extra
                    if code is not None:
                        ek=P.dump(expr)
                        if ek not in expression_ids:
                            expression_ids[ek]=len(expressions);expressions.append(expr)
                        reason=[*ij,code,u,v,expression_ids[ek]]
                        reason_counts[{1:'positive_actual_wall_product',2:'actual_wall_product_in_open_minus1_zero',3:'nonintegral_native_two_cycle',4:'native_two_cycle_outside_fixed_K',5:'forbidden_finite_Coxeter_product'}[code]]+=1
                        break
            pair_cache[ij]=reason is not None
            if reason is not None:pair_cert.append(reason)
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
        best=mu.get(vi,F());local_ids=[]
        for coeff,cert in sorted(candidates.items()):
            load=sum((m*alpha[fi]for fi,m in coeff),F());idx=len(columns);cost=cert.pop('cost')
            columns.append(dict(index=idx,vertex=vi,vertex_cost=cost,positive_face_coefficients=coeff,
                **cert,load=str(load),baseline_vertex_mass=str(mu.get(vi,F())),
                violated=load>cost*mu.get(vi,F())))
            local_ids.append(idx);best=max(best,load/cost)
        if best>mu.get(vi,F()):updates.append([vi,str(mu.get(vi,F())),str(best)])
        forest_records.append(dict(vertex=vi,forest=vertices[vi]['forest'],active_components=components,
            corner_occurrences=len(orecs),literal_stars=literal,corner_orbits=groups,
            incompatibility_certificates=sorted(pair_cert),radical_expressions=expressions,
            group_incompatibility_edges=sorted(tuple(sorted(e))for e in graph.edges),columns=local_ids,
            tested_literal_pairs=len(pair_cache),unknown_literal_pairs=sum(not b for b in pair_cache.values())))
        counts.update(dict(literal_orientations=len(keys),corner_orbits=len(groups),
            capacity_one_orbits=len(graph),group_edges=graph.number_of_edges(),literal_pairs_tested=len(pair_cache),
            exact_literal_exclusions=len(pair_cert),columns=len(candidates)))
        counts.update(reason_counts)
        if step%25==0 or step==len(occurrences):
            print('FORESTS',step,'/',len(occurrences),'stars',counts['literal_orientations'],
                  'columns',len(columns),'violated forests',len(updates),flush=True)
    new_mu=dict(mu)
    for v,old,new in updates:assert new_mu.get(v,F())==F(old);new_mu[v]=F(new)
    for col in columns:
        assert sum((m*alpha[fi]for fi,m in col['positive_face_coefficients']),F())<=col['vertex_cost']*new_mu[col['vertex']]
    assert all(new_mu[v]>=x for v,x in mu.items())
    normalizer=sum(new_mu.values());objective=F(witness['objective_numerator']);ratio=objective/normalizer
    counts.update(actual_Q_classes=len(selected),existing_explicit_classes=explicit,new_unique_context_classes=unique,corner_occurrences=sum(map(len,occurrences.values())),
        forests=len(occurrences),violated_columns=sum(c['violated']for c in columns),vertex_mass_increases=len(updates))
    paths=[Path(__file__),Path(P.__file__),Path(M.__file__),Path(V.__file__),Path(C.__file__),Path(X.__file__),Path(B.__file__),Path(B.Q.__file__)]
    paths +=[ROOT/f'results/{t}.json'for t in(WITNESS,FULL,P.FORESTS,M.TAGS['actual'],M.TAGS['quad'],M.TAGS['raw'],B.Q.SPECTATORS)]
    paths +=[ROOT/'results/fresh_quad_product_aggregation_20261017.md',ROOT/f'results/{P.TAG}.md',ROOT/'results/fresh_sqrt5_spectator_groundfield_gate_20261017.md']
    result=dict(status='PASS_EXACT_RAW24_UNIQUE_PROVENANCE_FINITE_VERTEX_CORNER_CAPACITIES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        raw_source=24,canonical_core=20,dotted_square=TAU.out(),
        omitted_ambiguous_raw24_rows=ambiguous_raw24,selected_faces=face_records,forest_records=forest_records,
        columns=columns,counts=dict(counts),
        fixed_witness_repair=dict(face_masses_unchanged=True,baseline=WITNESS,vertex_updates=updates,
            objective_numerator=str(objective),normalizer=str(normalizer),bound=str(ratio),bound_below12=ratio<12,
            baseline_bound=witness['bound'],optimizer_invocations=0,
            preservation='All face masses fixed and only nonnegative vertex allowances increase; every inherited signed/packing/joint/cap/incidence inequality remains valid by monotonicity.'),
        certificate_schema=dict(pair_reason_codes={'0':'shared occupied square pair','1':'positive actual-wall product','2':'actual-wall product in (-1,0)','3':'nonintegral native two-cycle','4':'native two-cycle outside fixed K','5':'forbidden finite Coxeter product'},pair_records='[star_a,star_b,0] or [star_a,star_b,reason,wall_a,wall_b,radical_expression_index]'),
        face_multipliers=witness['face_multipliers'],initial_vertex_multipliers=witness['vertex_multipliers'],
        vertex_multipliers=[[v,str(a)]for v,a in sorted(new_mu.items())],
        scope='The 45 positive raw24 classes with already explicit or uniquely saved provenance in the stated fixed witness, all their eight finite prism corners and exact active orientation orbits. Exact pair-exclusion cliques and proper colorings of the conservative compatibility graph give the stated capacities. Every unknown pair remains allowed; coloring need not be optimal. No optimum, feasibility of a target12 rule, full joint-configuration classification, or dimension theorem is asserted.')
    print('COMPLETE',dict(counts),'bound',float(ratio),'seconds',round(time.monotonic()-start,2),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(P.dump(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(P.dump(d)+'\n')
    print(d['status'],flush=True)
