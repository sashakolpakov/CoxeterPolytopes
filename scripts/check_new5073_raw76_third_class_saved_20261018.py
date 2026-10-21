#!/usr/bin/env python3
"""Independent stdlib saved accounting for the third raw76 class.

Checks actual occurrence accounting, complete new-pair index domains,
charge/cover certificates, unchanged record transport and multiplier repair.
Inherits native geometry, historical graph ownership and clique completeness.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_new5073_raw76_third_class_local_covers_20261018'
BASE='fresh_new5073_raw76_high_mass_local_covers_20261018'
LOCAL='fresh_current324_new5073_local_covers_20261018'


def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def run():
    d=read(TAG);assert d['status']=='PASS_NEW5073_RAW76_THIRD_CLASS_EXACT_LOCAL_COVERS'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    base=read(BASE);local=read(LOCAL)
    assert base['status']=='PASS_NEW5073_RAW76_HIGH_MASS_EXACT_LOCAL_COVERS'
    assert base['represented_classes']==326 and base['represented_forests']==600
    assert d['face_multipliers']==base['face_multipliers']==local['face_multipliers']
    alpha={f:F(x) for f,x in base['face_multipliers']};mu={v:F(x) for v,x in base['vertex_multipliers']}
    oldmu={v:F(x) for v,x in local['vertex_multipliers']}
    old={r['vertex']:r for r in local['forest_records']};assert len(old)==600
    prior={r['vertex']:r for r in base['forest_records']};assert len(prior)==6
    assert d['selected_face']==42946 and 42946 not in local['represented_classes']
    assert 42946 not in {r['face'] for r in base['selected_faces']}
    selected={r['face']:r for r in d['raw_domain']['selected_faces']}
    assert set(selected)=={42946} and len(selected[42946]['corners'])==4
    sr=next(r for r in base['screen']['rows'] if r['face']==42946)
    assert F(sr['mass'])==alpha[42946]==F(31942181,10000000000)
    assert sr['source_context']==59136 and sr['spectator_forest']==['A1','A1','A3','H3']
    forests=read('fresh_sqrt5_rank13_combined_disjoint_ports_20261017')
    quads=read('fresh_sqrt5_rank13_packing_quad_separation_20261017')
    spectators=read('fresh_sqrt5_rank13_dotted_b1_spectators_20261017')
    assert forests['faces'][42946]['source_contexts']==[[59136,65]]
    assert spectators['contexts'][59136]['spectator_forest']==sr['spectator_forest']
    lookup={tuple(r['forest']):r['index'] for r in forests['vertices']}
    cm=Counter(lookup[tuple(sorted(sr['spectator_forest']+r['forest']))]
        for r in quads['canonical_cores'][65]['aggregate_corner_occurrences'])
    assert cm==Counter({534:2,788:1,1665:1}) and dict(sr['corner_multiplicities'])==dict(cm)
    expected_origins={(42946,tuple(r['omitted']),r['vertex']) for r in selected[42946]['corners']}
    assert len(expected_origins)==4 and Counter(v for _,_,v in expected_origins)==cm
    rawmap={r['vertex']:r for r in d['raw_domain']['forest_records']}
    changed={r['vertex']:r for r in d['forest_records']}
    assert set(changed)==set(rawmap)==set(cm) and len(d['forest_records'])==3
    assert not set(changed)&set(prior) and all(mu[v]==oldmu[v] for v in changed)
    num=F(base['objective_numerator']);den=F(base['normalizer'])
    actual=set();counts=Counter();newmu=dict(mu);columns=[];ineq=terms=0
    geometry={}
    for v,r in changed.items():
        u=r['domain'];groups=u['corner_orbits'];nodes=u['literal_stars'];n=len(nodes);oldn=u['old_literal_stars']
        assert u['baseline']['record_sha256']==digest(old[v]) and u['baseline']['exact_local_optimum']==old[v]['exact_local_optimum']
        assert n-oldn==u['new_literal_stars'] and [z['index'] for z in groups]==list(range(len(groups)))
        owner={i:z['index'] for z in groups for i in z['literal_stars']}
        assert len(owner)==sum(len(z['literal_stars']) for z in groups)==n and set(owner)==set(range(n))
        assert [z['positive_face_coefficients'] for z in groups[:len(old[v]['demands'])]]==old[v]['orbit_coefficients']
        for z in groups:assert Counter(o['face'] for o in z['origins'])==Counter(dict(z['positive_face_coefficients']))
        for o in u['added_occurrences']:
            origin=o['origin'];k=(origin['face'],tuple(origin['omitted']),v)
            assert k not in actual;actual.add(k)
            assert o['literal_star_indices']==groups[o['orbit']]['literal_stars']
            assert len(o['literal_star_indices'])==origin['literal_orientations']
            assert origin in selected[origin['face']]['corners']
        pp=r['physical'];expected_pairs=[(i,j) for i,j in combinations(range(n),2) if j>=oldn]
        assert [tuple(z[:2]) for z in pp['newly_owned_pair_domain']]==expected_pairs
        rejected={tuple(z[:2]):z[2] for z in pp['newly_owned_pair_domain'] if z[2] is not None}
        assert len(rejected)==pp['new_exclusions'] and len(expected_pairs)-len(rejected)==pp['new_unknown_pairs']
        assert {tuple(z[:2]):z[2] for z in pp['incompatibility_certificates']}==rejected
        for i,j,code,*rest in pp['incompatibility_certificates']:
            if code==0:assert set(map(tuple,nodes[i]['occupied_pairs']))&set(map(tuple,nodes[j]['occupied_pairs']))
            else:assert len(rest)==3 and rest[0] in (0,1) and rest[1] in (0,1) and 0<=rest[2]<len(pp['wall_expressions'])
        vectors=r['maximal_count_vectors'];reps=r['literal_representatives'];assert len(vectors)==len(reps)
        for a,cl in zip(vectors,reps):
            assert len(cl)==len(set(cl)) and all(0<=i<n for i in cl)
            assert not any(tuple(sorted(p)) in rejected for p in combinations(cl,2))
            assert a==[Counter(owner[i] for i in cl).get(k,0) for k in range(len(groups))]
        demand=[sum((alpha.get(f,F())*m for f,m in z['positive_face_coefficients']),F()) for z in groups]
        assert list(map(str,demand))==u['demands']
        cert=r['certificate'];w=list(map(F,cert['weights']));assert len(w)==len(groups) and all(x>=0 for x in w)
        assert all(len(a)==len(w) and sum(n*x for n,x in zip(a,w))<=1 for a in vectors);ineq+=len(vectors)
        cover=[F() for _ in w];mass=F()
        for ci,weight in cert['mixture']:
            assert 0<=ci<len(vectors)
            a=vectors[ci];x=F(weight);assert x>=0
            mass+=x;terms+=1
            for i,n in enumerate(a):cover[i]+=n*x
        assert all(x>=y for x,y in zip(cover,demand))
        value=sum(x*y for x,y in zip(w,demand));assert value==mass==F(cert['optimal_mass'])==F(r['exact_local_optimum'])
        assert value>=F(old[v]['exact_local_optimum'])
        newmu[v]=max(mu.get(v,F()),value)
        assert newmu[v]==F(r['repaired_vertex_mass']) and newmu[v]-mu.get(v,F())==F(r['multiplier_increase'])
        coeff=Counter()
        for x,z in zip(w,groups):
            for f,m in z['positive_face_coefficients']:coeff[f]+=m*x
        coeffs=[[f,str(x)] for f,x in sorted(coeff.items()) if x]
        assert coeffs==r['positive_face_coefficients']
        if value>mu.get(v,F()):columns.append(dict(vertex=v,vertex_cost='1',positive_face_coefficients=coeffs,exact_load=str(value)))
        counts.update(new_pairs=len(expected_pairs),new_exclusions=len(rejected),new_unknowns=len(expected_pairs)-len(rejected),
            new_stars=len(nodes)-oldn,actual_corners=len(u['added_occurrences']),maximal_cliques=r['maximal_cliques'],
            count_vectors=len(vectors),mixture_terms=len(cert['mixture']))
    assert actual==expected_origins and d['columns']==columns and d['counts']==dict(counts)
    latest={v:(LOCAL,r) for v,r in old.items()}
    latest.update({v:(BASE,r) for v,r in prior.items()})
    expected_inherited=[dict(vertex=v,source=s,record_sha256=digest(r))
        for v,(s,r) in sorted(latest.items()) if v not in changed]
    assert d['inherited_forests']==expected_inherited and len(expected_inherited)==597
    assert len([r for r in expected_inherited if r['source']==BASE])==6
    assert [[v,str(x)] for v,x in sorted(newmu.items()) if x]==d['vertex_multipliers']
    model=read('fresh_current324_full_model_20261018');loads=[];two=0
    for c in model['prepared']['inherited_several_vertex_columns']+model['prepared']['new_columns']:
        load=sum(alpha.get(f,F())*m for f,m in c['integer_face_coefficients'])
        before=sum(mu.get(v,F())*m for v,m in c['integer_vertex_capacities'])
        after=sum(newmu.get(v,F())*m for v,m in c['integer_vertex_capacities']);assert load<=before<=after
        two+=len(c['integer_vertex_capacities'])>1
        loads.append(dict(variable=c['variable'],all_vertex_costs=c['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after)))
    assert loads==d['checked_recent_model_loads'] and len(loads)==412 and two==14
    newden=sum(newmu.values(),F());assert newden==F(d['normalizer']) and num/newden==F(d['bound'])
    assert num==F(d['objective_numerator']) and num/den==F(d['baseline_ratio'])
    assert d['represented_classes']==327 and d['represented_forests']==600 and d['global_optimizer_invocations']==0
    assert d['changed_forests']==3 and len(changed)==3
    assert d['previous_patch_vertices']==sorted(prior)
    assert d['previous_unknown_pairs_preserved']==base['counts']['new_unknowns']==312
    assert d['new_dimension_bound'] is None
    print('PASS_SAVED_RAW76_THIRD_CLASS',counts,'inequalities',ineq,'terms',terms,'ratio',float(num/newden),flush=True)


if __name__=='__main__':run()
