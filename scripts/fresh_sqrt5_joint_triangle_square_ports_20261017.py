#!/usr/bin/env python3
"""Exact shared three-port charge 2*s_U+t_U+q_U <=3."""
import argparse
import json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path

import fresh_sqrt5_prism_simplex_base_incidence_20261017 as B
import fresh_mixed_face_port_incidence_20261017 as X

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_joint_triangle_square_ports_20261017'
LOCAL='fresh_sqrt5_restricted_flag_local_cones_20261017'
MODEL='fresh_sqrt5_flag_local_packing_full_model_20261017'
ACTUAL='fresh_sqrt5_actual_provenance_dual_20261017'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
C=B.B.C


def dump(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def read(tag):return B.P.read(tag)
def freeze(x):return tuple(map(freeze,x)) if isinstance(x,list) else x


def triples(g):
    pieces=C.pieces(g)
    groups=defaultdict(list)
    for U in combinations(sorted(g),3):
        code=tuple(sorted(C.code(p,U) for p in pieces))
        groups[code].append(U)
    assert sum(map(len,groups.values()))==286
    return groups


def small_witness():
    # All three mutually intersecting two-port blocks: their complete
    # disjoint-packing polytope is simply sum of their masses <=1.
    pairs=list(combinations(range(3),2))
    s=q=F(0)
    loads=[]
    for pair in pairs:
        s+=F(1,3);q+=F(2,3)
        loads.append(F(1,3)+F(2,3))
    assert s==1 and loads==[F(1)]*3 and 2*s+q==4>3
    return dict(simplex_blocks=[dict(ports=list(p),mass='1/3')for p in pairs],
        quadrilateral_single_pair_footprints=[dict(ports=list(p),mass='2/3')for p in pairs],
        pair_loads=list(map(str,loads)),simplex_packing_mass='1',old_mixed_triple_load='0',
        new_joint_load='4',capacity='3',
        scope='Fractional relaxation only; no common arithmetic geometry or global flag realization asserted.')


def run():
    local=read(LOCAL)
    t=next(t for t in local['tests'] if t['name']=='legacy_export_seed')
    assert t['complete_positive_simplex_forest_coverage'] and not t['open_forests']
    alpha={i:F(a) for i,a in t['face_multipliers']}
    mu={v:F(a) for v,a in t['vertex_multipliers']}
    assert len(alpha)==2219
    model=json.loads((ROOT/f'results/{MODEL}.json').read_text())
    for p,h in model['prepared']['sources'].items():assert B.P.sha(p)==h,p
    pmap=dict(model['prepared']['provenance_to_face']);inverse={f:j for j,f in pmap.items()}
    selected_q=sorted(set(alpha)&set(inverse));assert len(selected_q)==355
    actual=read(ACTUAL);quad=read(B.Q.TAG);raw=read(B.Q.CORES);spectators=read(B.Q.SPECTATORS)
    forests=read(FORESTS)['vertices']
    assert len(forests)==6558
    forest_lookup={tuple(v['forest']):v['index']for v in forests};assert len(forest_lookup)==6558
    models=dict(C.models());finite={C.code(g):n for n,g in models.items()}
    graphs={};orbit_data={};canonical_u={};code_to_oi={}

    def ensure(v):
        if v in graphs:return
        g={}
        for name in forests[v]['forest']:
            h=models[name];offset=len(g)
            for i,row in h.items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert set(g)==set(range(13))
        graphs[v]=g
        groups=triples(g)
        rows=[]
        for oi,(code,us) in enumerate(sorted(groups.items())):
            rows.append(dict(orbit=oi,code=json.loads(dump(code)),representative=list(us[0]),
                             multiplicity=len(us)))
            canonical_u[v,oi]=sum(1<<i for i in us[0])
            code_to_oi[v,code]=oi
        orbit_data[v]=rows

    blockcols={c['index']:c for c in t['literal_block_orbits']}
    blockweights=defaultdict(F)
    for f in t['positive_simplex_features']:
        for j,m in f['block_orbits']:blockweights[j]+=m*alpha[f['face']]
    assert set(blockweights)==set(blockcols)
    block_coeff={};sload=defaultdict(F)
    for j,c in sorted(blockcols.items()):
        v=c['vertex'];ensure(v);masks=c['literal_masks']
        assert masks and len(set(masks))==len(masks)
        assert all(0<m<1<<13 and m.bit_count()==c['block_size']for m in masks)
        for row in orbit_data[v]:
            oi=row['orbit'];u=canonical_u[v,oi]
            a=F(sum(2*((m&u).bit_count()>=2)+((m&u).bit_count()==3) for m in masks),len(masks))
            if a:
                block_coeff[j,oi]=a;sload[v,oi]+=blockweights[j]*a
        k=c['block_size']
        assert sum((block_coeff.get((j,r['orbit']),F(0))*r['multiplicity']
                    for r in orbit_data[v]),F(0))==2*comb(k,2)*(13-k)+3*comb(k,3)
    print('JOINT_SIMPLEX_FEATURES',len(t['positive_simplex_features']),len(blockcols),len(graphs),flush=True)
    sfacecoeff=defaultdict(Counter)
    for f in t['positive_simplex_features']:
        for j,m in f['block_orbits']:
            v=blockcols[j]['vertex']
            for row in orbit_data[v]:
                oi=row['orbit'];a=block_coeff.get((j,oi),F(0))
                if a:sfacecoeff[v,oi][f['face']]+=m*a

    qload=defaultdict(F);qcoeff=defaultdict(Counter);qrecords=[]
    for fi in selected_q:
        j=inverse[fi];rec=actual['provenance_masses'][j]
        ci=rec['canonical_prism_core'];core=quad['canonical_cores'][ci]
        g0,_=B.Q.graph(raw['records'][core['source_core']])
        g={i:{k:m[0]for k,m in row.items()}for i,row in g0.items()}
        sp=spectators['contexts'][rec['source_context']]['spectator_forest']
        for name in sp:
            h=models[name];off=len(g)
            for i,row in h.items():g[off+i]={off+k:m for k,m in row.items()}
        assert len(g)==15
        cs=X.quad_features(dict(quadrilateral_black_partitions=core['quad_partitions']))
        assert len(cs)==2*rec['prism_dimension']
        corner_rows=[]
        for c in cs:
            omitted=set(c['vertex_omitted'])
            h={i:{k:m for k,m in row.items()if k not in omitted}
               for i,row in g.items()if i not in omitted}
            assert len(h)==13 and all(m>0 for row in h.values() for m in row.values())
            ps=C.pieces(h);forest=tuple(sorted(finite[C.code(p)]for p in ps));v=forest_lookup[forest]
            ensure(v)
            edges=set(map(tuple,c['occupied_pairs']));assert len(edges)==rec['prism_dimension']-1
            groups=triples(h);hits=Counter()
            for code,us in groups.items():
                oi=code_to_oi[v,code]
                assert len(us)==orbit_data[v][oi]['multiplicity']
                hit=sum(any(tuple(e)in edges for e in combinations(U,2))for U in us)
                if hit:
                    a=F(hit,len(us));hits[oi]=a;qcoeff[v,oi][fi]+=a;qload[v,oi]+=alpha[fi]*a
            e=len(edges)
            assert sum((a*orbit_data[v][oi]['multiplicity'] for oi,a in hits.items()),F(0))==11*e-comb(e,2)
            corner_rows.append(dict(omitted=sorted(omitted),vertex=v,occupied_pairs=c['occupied_pairs'],
                orbit_coefficients=[[oi,str(a)]for oi,a in sorted(hits.items())]))
        qrecords.append(dict(face=fi,provenance=j,canonical_core=ci,
            prism_dimension=rec['prism_dimension'],spectator_forest=sp,corners=corner_rows))
    print('JOINT_QUAD_FEATURES',len(qrecords),sum(len(r['corners'])for r in qrecords),len(graphs),flush=True)
    # Export only genuinely violated columns, but check every orbit load.
    updates=dict(mu);checks=[];critical=[]
    for v,rows in sorted(orbit_data.items()):
        for row in rows:
            oi=row['orbit'];s,q=sload[v,oi],qload[v,oi];load=s+q
            need=load/3
            checks.append([v,oi,str(s),str(q),str(load-3*mu.get(v,F(0)))])
            updates[v]=max(updates.get(v,F(0)),need)
            if need>mu.get(v,F(0)):
                terms=Counter(qcoeff[v,oi])
                terms.update(sfacecoeff[v,oi])
                assert sum((alpha[fi]*a for fi,a in terms.items()),F(0))==load
                critical.append(dict(index=len(critical),vertex=v,three_port_orbit=oi,
                    representative=row['representative'],orbit_multiplicity=row['multiplicity'],
                    vertex_cost=3,positive_face_coefficients=[[fi,str(a)]for fi,a in sorted(terms.items())],
                    simplex_load=str(s),quad_load=str(q),load=str(load),old_vertex_mass=str(mu.get(v,F(0)))))
    numerator=F(t['objective_numerator']);den=sum(updates.values(),F(0));bound=numerator/den
    assert all(updates[v]>=a for v,a in mu.items())
    for v,oi,s,q,_ in checks:assert F(s)+F(q)<=3*updates[v]
    paths=[f'scripts/{TAG}.py',f'scripts/{LOCAL}.py',f'results/{LOCAL}.json',
        f'results/{MODEL}.json',f'results/{ACTUAL}.json',f'results/{FORESTS}.json',
        f'scripts/{X.TAG}.py',f'results/{X.TAG}.json',f'results/{B.Q.TAG}.json',
        f'results/{B.Q.CORES}.json',f'results/{B.Q.SPECTATORS}.json',
        'results/fresh_prism_simplex_shared_capacity_20261017.md',
        'results/fresh_quad_product_aggregation_20261017.md',
        'results/fresh_simplex_disjoint_ports_20261017.md']
    return dict(status='PASS_EXACT_JOINT_TRIANGLE_SQUARE_THREE_PORT_DIAGNOSTIC',
        sources={p:B.P.sha(p)for p in paths},dimension=13,target=12,optimizer_invocations=0,
        initial_bound=t['bound'],bound=str(bound),objective_numerator=str(numerator),normalizer=str(den),
        counts=dict(simplex_faces=len(t['positive_simplex_features']),simplex_block_orbits=len(blockcols),
            actual_quad_provenances=len(qrecords),quad_vertex_occurrences=sum(len(r['corners'])for r in qrecords),
            tested_forests=len(graphs),three_port_orbits=len(checks),violated_columns=len(critical),
            increased_vertex_masses=sum(updates[v]>mu.get(v,F(0))for v in updates)),
        finite_forest_triple_orbits=[dict(vertex=v,forest=forests[v]['forest'],orbits=rows)
                                   for v,rows in sorted(orbit_data.items())],
        actual_quad_features=qrecords,all_orbit_loads_and_slacks=checks,critical_columns=critical,
        fractional_relaxation_witness=small_witness(),face_multipliers=t['face_multipliers'],
        vertex_multipliers=[[v,str(a)]for v,a in sorted(updates.items())],
        scope='Compact actual faces only. Triangle-port transitivity and triangle-plus-two-squares completion imply 2*s_U+t_U+q_U<=3, where t_U indicates all3ports in one S block; each canonical quad aggregate counts once even when it occupies two U-pairs. Exact full-forest automorphism averages and one common actual-provenance mass are used. All positive legacy-restricted S and355 actual Q are covered; other positive coefficients are set to zero, giving valid weaker charge columns. Old full local-simplex covers survive because face masses are unchanged and only mu increases. No optimizer, joint realizability, complete joint local polytope or dimension claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(path.read_text())
    else:path.write_text(dump(result)+'\n')
    print(result['status'],result['counts'],result['bound'],float(F(result['bound'])),flush=True)
