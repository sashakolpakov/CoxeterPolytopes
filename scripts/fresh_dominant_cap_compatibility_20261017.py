#!/usr/bin/env python3
"""Exact rooted-corner incompatibility at one dominant finite vertex.

No optimizer. Seven pairs share an actual quadrilateral port pair; eight
remaining pairs force a native product strictly between minus one and zero.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import networkx as nx
import fresh_sqrt5_multispectator_cap_capacity_20261017 as M
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C
import fresh_mixed_face_port_incidence_20261017 as X
import fresh_sqrt5_prism_simplex_base_incidence_20261017 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_dominant_cap_compatibility_20261017'
WITNESS='fresh_sqrt5_scaled_core_cap_edge_repair_20261017'
FULL='fresh_sqrt5_cross_coface_cap_full_model_20261017'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
V,K=M.V,M.K
TARGET=3137
SELECTED={26448:490,26474:497}


def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def dump(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def dot(p,J,q):return sum((p[i]*J[i][j]*q[j]for i in range(len(p))for j in range(len(q))),K())
def ngraph(g):
    h=nx.Graph();h.add_nodes_from(g)
    for i,row in g.items():
        for j,m in row.items():h.add_edge(i,j,label=m)
    return h

def forest_graph(forest,models):
    g={}
    for name in forest:
        off=len(g)
        for i,row in models[name].items():g[off+i]={off+j:m for j,m in row.items()}
    return g


def run():
    witness=read(WITNESS);full=read(FULL);actual=read(M.TAGS['actual'])
    quads=read(M.TAGS['quad']);raw=read(M.TAGS['raw']);spectators=read(B.Q.SPECTATORS)
    forests=read(FORESTS)['vertices'];lookup={tuple(v['forest']):v['index']for v in forests}
    alpha={i:F(a)for i,a in witness['face_multipliers']};mu={v:F(a)for v,a in witness['vertex_multipliers']}
    assert len(alpha)==2524 and all(a>=0 for a in alpha.values())and all(a>=0 for a in mu.values())
    pmap=dict(full['prepared']['provenance_to_face'])
    assert len(set(pmap.values()))==len(pmap) and all(pmap[j]==fi for fi,j in SELECTED.items())
    forest=forests[TARGET]['forest'];assert forest==['A1','A3','H3','H4','I2(5)']
    models=dict(C.models());finite={C.code(g):n for n,g in models.items()}
    g=forest_graph(forest,models);assert sorted(g)==list(range(13))
    G=[[K(2)if i==j else -M.ksqrt(V.squared(V.edge(g[i][j])))if j in g[i]else K()for j in range(13)]for i in range(13)]
    assert V.pd(G);J,det=V.inverse(G);assert det.sign()>0
    # All components have distinct types. Only A3 reversal and I2(5) reversal
    # are nontrivial weighted automorphisms, yielding these four permutations.
    aut=[]
    for a,b in product((False,True),repeat=2):
        p=list(range(13))
        if a:p[1],p[3]=p[3],p[1]
        if b:p[11],p[12]=p[12],p[11]
        assert all(g[i].get(j,2)==g[p[i]].get(p[j],2)for i in g for j in g)
        aut.append(tuple(p))
    graph_aut=[tuple(p[i]for i in range(13))for p in nx.algorithms.isomorphism.GraphMatcher(ngraph(g),ngraph(g),edge_match=lambda x,y:x['label']==y['label']).isomorphisms_iter()]
    assert set(graph_aut)==set(aut) and len(graph_aut)==4
    tau=K(F(7,2),F(1,2));stars={};corner_records=[];counts=Counter();full_records=[]
    for fi,j in sorted(SELECTED.items()):
        pr=actual['provenance_masses'][j];cc=quads['canonical_cores'][pr['canonical_prism_core']]
        assert pr['prism_dimension']==5 and cc['source_core']==23
        rg=V.raw_graph(raw['records'][23]);sp=spectators['contexts'][pr['source_context']]['spectator_forest']
        gg=[list(row)for row in rg]
        for name in sp:
            h=models[name];off=len(gg);z=len(h)
            gg=[row+[V.RIGHT]*z for row in gg]+[[V.RIGHT]*(off+z)for _ in range(z)]
            for x,row in h.items():
                for y,m in row.items():gg[off+x][off+y]=V.edge(m)
        assert len(gg)==15
        cs=X.quad_features(dict(quadrilateral_black_partitions=cc['quad_partitions']));assert len(cs)==10
        all_corners=[]
        for c in cs:
            omitted=c['vertex_omitted'];ids=[i for i in range(15)if i not in omitted]
            h={i:{z:gg[i][z][0]for z in ids if z!=i and gg[i][z]!=V.RIGHT}for i in ids}
            assert all(m>0 for row in h.values()for m in row.values())
            vf=tuple(sorted(finite[C.code(p)]for p in C.pieces(h)));vi=lookup[vf]
            all_corners.append(dict(omitted=omitted,vertex=vi,forest=vf))
            if vi!=TARGET:continue
            counts[fi]+=1;maps=list(nx.algorithms.isomorphism.GraphMatcher(ngraph(h),ngraph(g),edge_match=lambda x,y:x['label']==y['label']).isomorphisms_iter())
            assert len(maps)==4
            first=maps[0];expected={tuple(p[first[i]]for i in ids)for p in aut}
            assert {tuple(mp[i]for i in ids)for mp in maps}==expected
            ckeys=[]
            for mp in maps:
                ex=[]
                for x in omitted:
                    es=[V.RIGHT]*13
                    for old,new in mp.items():es[new]=gg[x][old]
                    ex.append(tuple(es))
                key=tuple(sorted(ex));occupied=tuple(sorted(tuple(sorted(mp[z]for z in e))for e in c['occupied_pairs']))
                assert len(occupied)==4
                if key in stars:assert stars[key]['occupied_pairs']==occupied
                else:stars[key]=dict(contacts=key,occupied_pairs=occupied,origins=[])
                origin=dict(face=fi,provenance=j,omitted=omitted,isomorphism=sorted(mp.items()))
                if origin not in stars[key]['origins']:stars[key]['origins'].append(origin)
                ckeys.append(key)
            corner_records.append(dict(face=fi,provenance=j,omitted=omitted,isomorphisms=[sorted(mp.items())for mp in maps],literal_star_keys=sorted(set(ckeys))))
        full_records.append(dict(face=fi,provenance=j,canonical_prism_core=pr['canonical_prism_core'],raw_source=23,
            source_context=pr['source_context'],spectator_forest=sp,all_corners=all_corners))
    assert counts=={26448:1,26474:2} and len(stars)==6 and len(corner_records)==3
    ordered=[];index={}
    for key,r in sorted(stars.items()):
        ix=len(ordered);index[key]=ix
        candidates=[]
        for k,es in enumerate(key):
            ports=[i for i,e in enumerate(es)if e!=V.RIGHT]
            if len(ports)==1 and es[ports[0]][0]==0:
                q=ports[0];assert V.squared(es[q])==tau;candidates.append((k,q))
        assert len(candidates)==1
        k,q=candidates[0];other=key[1-k]
        p=[-M.ksqrt(V.squared(e))if e!=V.RIGHT else K()for e in other]
        h0=tau*J[q][q]-2;h1=dot(p,J,p)-2;assert h0.sign()>0 and h1.sign()>0
        b=-sum((J[q][z]*p[z]for z in range(13)),K());assert b.sign()>=0
        # Select the common future side and verify the omitted-root product.
        origin=r['origins'][0];u,v=origin['omitted'];z=-M.ksqrt(V.squared(rg[u][v]))if rg[u][v]!=V.RIGHT else K()
        assert not b*z and h0*h1==b*b*tau+z*z and (b.sign()>0 or z.sign()<0)
        if q==0:assert h0==K(F(-1,4),F(1,4))
        else:assert q in(1,3) and h0==K(F(5,8),F(3,8))
        ordered.append(dict(index=ix,contacts=key,occupied_pairs=r['occupied_pairs'],origins=r['origins'],
            forced_wall=k,single_contact_root=q,squared_timelike_height=h0.out(),
            other_covector=[x.out()for x in p],other_squared_height=h1.out(),
            omitted_product=z.out(),projection_coefficient_of_sqrt_tau=b.out()))
    pairs=[]
    for a,b in combinations(ordered,2):
        shared=sorted(set(a['occupied_pairs'])&set(b['occupied_pairs']))
        if shared:
            # Their external contact pairs differ, hence they are distinct
            # full wall stars; one actual quadrilateral has only one such star.
            assert a['contacts']!=b['contacts']
            rec=dict(stars=[a['index'],b['index']],reason='distinct_full_stars_share_actual_quad_port_pair',shared_pairs=shared)
        else:
            i,j=a['single_contact_root'],b['single_contact_root'];assert J[i][j]==0
            h=K(*map(F,a['squared_timelike_height']))*K(*map(F,b['squared_timelike_height']))
            assert h==K(F(5,16),F(1,16)) and h.sign()>0 and (1-h).sign()>0
            rec=dict(stars=[a['index'],b['index']],reason='forced_distinct_wall_product_in_open_minus1_zero_gap',
                single_contact_roots=[i,j],projection_product='0',negative_product_squared=h.out())
        pairs.append(rec)
    assert len(pairs)==15 and Counter(r['reason']for r in pairs)=={
        'distinct_full_stars_share_actual_quad_port_pair':7,'forced_distinct_wall_product_in_open_minus1_zero_gap':8}
    for c in corner_records:
        c['literal_stars']=[index[key]for key in c.pop('literal_star_keys')]
    positive=sorted(counts.items());load=sum((m*alpha[fi]for fi,m in positive),F())
    assert load==F(310256599,10000000000) and mu[TARGET]==F(200076613,10000000000)
    gap=load-mu[TARGET];assert gap==F(55089993,5000000000)>0
    normalizer=F(witness['normalizer'])+gap;ratio=F(witness['objective_numerator'])/normalizer
    assert ratio<12
    paths=[Path(__file__),Path(M.__file__),Path(V.__file__),Path(C.__file__),Path(X.__file__),Path(B.__file__),Path(B.Q.__file__)]
    paths +=[ROOT/f'results/{t}.json'for t in(WITNESS,FULL,FORESTS,M.TAGS['actual'],M.TAGS['quad'],M.TAGS['raw'],B.Q.SPECTATORS)]
    paths +=[ROOT/'results/fresh_quad_product_aggregation_20261017.md']
    return dict(status='PASS_EXACT_DOMINANT_FINITE_VERTEX_Q_CORNER_INCOMPATIBILITY',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        target_vertex=TARGET,target_forest=forest,incident_gram=[[x.out()for x in row]for row in G],
        incident_inverse=[[x.out()for x in row]for row in J],automorphisms=aut,dotted_square=tau.out(),
        full_prism_records=full_records,corner_occurrences=corner_records,literal_stars=ordered,incompatibilities=pairs,
        inequality=dict(vertex=TARGET,vertex_cost=1,positive_face_coefficients=positive,load=str(load),
            old_vertex_mass=str(mu[TARGET]),gap=str(gap)),
        fixed_witness_repair=dict(face_masses_unchanged=True,vertex_updates=[[TARGET,str(mu[TARGET]),str(load)]],
            objective_numerator=witness['objective_numerator'],normalizer=str(normalizer),bound=str(ratio),
            bound_below12=ratio<12,old_bound=witness['bound'],optimizer_invocations=0,
            preservation='All old face masses fixed; one nonnegative vertex allowance increases. Every old signed constraint and local/joint/cap/incidence charge remains valid by monotonicity.'),
        counts=dict(actual_Q_classes=2,target_corner_occurrences=3,literal_orientations=6,pair_checks=15,
            shared_quad_exclusions=7,open_angle_gap_exclusions=8),
        scope='One necessary actual finite-vertex incidence inequality, for the two stated full-star classes. It invalidates the current fixed dual; its monotone repair below12 is neither optimized feasibility nor a dimension theorem. No global optimizer or historical large scan is run.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(dump(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['inequality'],float(F(d['fixed_witness_repair']['bound'])),flush=True)
