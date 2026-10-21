#!/usr/bin/env python3
"""Exact finite-vertex corner-orbit cliques for every positive raw23 prism.

No optimizer. Unknown pairs are left compatible. Inactive spectator
automorphisms are quotiented only after retaining every active embedding.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time
import networkx as nx
import fresh_dominant_cap_compatibility_20261017 as P

ROOT=P.ROOT
TAG='fresh_dominant_cap_compatibility_family_20261017'
M,V,K,C,X,B=P.M,P.V,P.K,P.C,P.X,P.B
TAU=K(F(7,2),F(1,2))
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
    """Only 3/5 components can carry raw23 covectors.

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
    witness=P.read(P.WITNESS);prior=P.read(P.TAG);full=P.read(P.FULL)
    actual=P.read(M.TAGS['actual']);quads=P.read(M.TAGS['quad'])
    raw=P.read(M.TAGS['raw']);spectators=P.read(B.Q.SPECTATORS)
    vertices=P.read(P.FORESTS)['vertices']
    lookup={tuple(v['forest']):v['index']for v in vertices}
    alpha={i:F(a)for i,a in witness['face_multipliers']}
    mu={v:F(a)for v,a in witness['vertex_multipliers']}
    assert all(a>=0 for a in alpha.values())and all(a>=0 for a in mu.values())
    for v,old,new in prior['fixed_witness_repair']['vertex_updates']:
        assert mu[v]==F(old);mu[v]=F(new)
    assert sum(mu.values())==F(prior['fixed_witness_repair']['normalizer'])
    pmap=dict(full['prepared']['provenance_to_face'])
    assert len(set(pmap.values()))==len(pmap)
    models=dict(C.models());finite={C.code(g):n for n,g in models.items()}
    selected=[]
    for pr in actual['provenance_masses']:
        j=pr['index'];fi=pmap[j]
        cc=quads['canonical_cores'][pr['canonical_prism_core']]
        if alpha.get(fi,0)>0 and cc['source_core']==23:selected.append((fi,j,pr,cc))
    assert len(selected)==57
    rg=V.raw_graph(raw['records'][23]);assert len(rg)==7
    assert all(e[0]in(0,2,3,5)for row in rg for e in row)
    occurrences=defaultdict(list);face_records=[];target_graphs={}
    for fi,j,pr,cc in selected:
        assert pr['prism_dimension']==5 and cc['source_core']==23
        sp=spectators['contexts'][pr['source_context']]['spectator_forest']
        corners=X.quad_features(dict(quadrilateral_black_partitions=cc['quad_partitions']))
        assert len(corners)==10
        assert len({tuple(sorted(c['vertex_omitted']))for c in corners})==10
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
                assert len(occ)==4 and len(set(occ))==4
                if key in orientations:assert orientations[key][0]==occ
                else:orientations[key]=(occ,sorted(mp.items()))
            record=dict(face=fi,provenance=j,omitted=omitted,vertex=vi,
                        active_embeddings=len(maps),literal_orientations=len(orientations))
            records.append(record)
            occurrences[vi].append((record,orientations,rg[omitted[0]][omitted[1]]))
        face_records.append(dict(face=fi,provenance=j,canonical_prism_core=pr['canonical_prism_core'],
            source_context=pr['source_context'],spectator_forest=sp,corners=records))
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
                if key in stars:assert stars[key]['occupied_pairs']==occupied
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
        pair_cache={};pair_cert=[]
        def incompatible(a,b):
            if a==b:return False
            ij=tuple(sorted((ix[a],ix[b])))
            if ij in pair_cache:return pair_cache[ij]
            sa,sb=stars[a],stars[b];shared=sorted(set(sa['occupied_pairs'])&set(sb['occupied_pairs']))
            reason=None
            if shared:
                reason=dict(stars=ij,reason='distinct_full_stars_share_actual_quad_port_pair',shared_pairs=shared)
            else:
                for u,v in ((0,0),(0,1),(1,0),(1,1)):
                    if a[u]==b[v]:continue # identical actual wall, not a forbidden pair
                    x,y=projection(sa['covectors'][u],sb['covectors'][v],J)
                    h=sa['heights'][u]*sb['heights'][v]
                    sg=signs_difference(x,y,h)
                    sg1=signs_difference(x+1,y,h)
                    if sg>0 or (sg<0 and sg1>0):
                        reason=dict(stars=ij,reason='positive_actual_wall_product'if sg>0 else 'actual_wall_product_in_open_minus1_zero',
                            walls=[u,v],projection=qout((x,y)),height_product=h.out(),
                            sign_product=sg,sign_product_plus_one=sg1)
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
        best=mu.get(vi,F())
        local_ids=[]
        for coeff,clique in sorted(local_columns.items()):
            load=sum((m*alpha[fi]for fi,m in coeff),F());idx=len(columns)
            columns.append(dict(index=idx,vertex=vi,vertex_cost=1,positive_face_coefficients=coeff,
                orbit_clique=clique,load=str(load),baseline_vertex_mass=str(mu.get(vi,F())),
                violated=load>mu.get(vi,F())))
            local_ids.append(idx);best=max(best,load)
        if best>mu.get(vi,F()):updates.append([vi,str(mu.get(vi,F())),str(best)])
        forest_records.append(dict(vertex=vi,forest=vertices[vi]['forest'],active_components=components,
            corner_occurrences=len(orecs),literal_stars=literal,corner_orbits=groups,
            incompatibility_certificates=sorted(pair_cert,key=lambda x:x['stars']),
            group_incompatibility_edges=sorted(tuple(sorted(e))for e in graph.edges),columns=local_ids,
            tested_literal_pairs=len(pair_cache),unknown_literal_pairs=sum(not b for b in pair_cache.values())))
        counts.update(dict(literal_orientations=len(keys),corner_orbits=len(groups),
            capacity_one_orbits=len(graph),group_edges=graph.number_of_edges(),literal_pairs_tested=len(pair_cache),
            exact_literal_exclusions=len(pair_cert),columns=len(local_columns)))
        counts.update(r['reason']for r in pair_cert)
        if step%25==0 or step==len(occurrences):
            print('FORESTS',step,'/',len(occurrences),'stars',counts['literal_orientations'],
                  'columns',len(columns),'violated forests',len(updates),flush=True)
    new_mu=dict(mu)
    for v,old,new in updates:assert new_mu.get(v,F())==F(old);new_mu[v]=F(new)
    for col in columns:
        assert sum((m*alpha[fi]for fi,m in col['positive_face_coefficients']),F())<=new_mu[col['vertex']]
    # The first six-orientation certificate is reproduced as one member of
    # this larger family, with the same actual-corner multiplicities.
    assert any(c['vertex']==prior['target_vertex'] and
               list(map(list,c['positive_face_coefficients']))==prior['inequality']['positive_face_coefficients']
               for c in columns)
    assert all(new_mu[v]>=x for v,x in mu.items())
    normalizer=sum(new_mu.values());objective=F(witness['objective_numerator']);ratio=objective/normalizer
    counts.update(actual_Q_classes=len(selected),corner_occurrences=sum(map(len,occurrences.values())),
        forests=len(occurrences),violated_columns=sum(c['violated']for c in columns),vertex_mass_increases=len(updates))
    paths=[Path(__file__),Path(P.__file__),Path(M.__file__),Path(V.__file__),Path(C.__file__),Path(X.__file__),Path(B.__file__),Path(B.Q.__file__)]
    paths +=[ROOT/f'results/{t}.json'for t in(P.TAG,P.WITNESS,P.FULL,P.FORESTS,M.TAGS['actual'],M.TAGS['quad'],M.TAGS['raw'],B.Q.SPECTATORS)]
    paths +=[ROOT/'results/fresh_quad_product_aggregation_20261017.md',ROOT/f'results/{P.TAG}.md']
    result=dict(status='PASS_EXACT_POSITIVE_RAW23_FINITE_VERTEX_CORNER_CLIQUE_FAMILY',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        raw_source=23,dotted_square=TAU.out(),selected_faces=face_records,forest_records=forest_records,
        columns=columns,counts=dict(counts),
        fixed_witness_repair=dict(face_masses_unchanged=True,baseline=P.TAG,vertex_updates=updates,
            objective_numerator=str(objective),normalizer=str(normalizer),bound=str(ratio),bound_below12=ratio<12,
            baseline_bound=prior['fixed_witness_repair']['bound'],optimizer_invocations=0,
            preservation='All face masses fixed and only nonnegative vertex allowances increase; every inherited signed/packing/joint/cap/incidence inequality remains valid by monotonicity.'),
        scope='Complete positive raw23 face cohort of the stated fixed witness, all its finite corners and exact active orientation orbits. Only proven pair-exclusion clique charges are imposed; unknown pairs remain allowed. No optimum, feasibility of a target12 rule, full joint-configuration classification, or dimension theorem is asserted.')
    print('COMPLETE',dict(counts),'bound',float(ratio),'seconds',round(time.monotonic()-start,2),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(P.dump(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(P.dump(d)+'\n')
    print(d['status'],flush=True)
