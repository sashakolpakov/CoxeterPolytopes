#!/usr/bin/env python3
"""Independent stdlib saved arithmetic check; no project imports or LPs.

Reconstructs chronological saved cone coordinates, new exact demands, every
charge/cover equality and monotone repair. Native geometry, owner chronology
inside each frozen graph and maximal-clique completeness remain bound inputs.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_current324_new5073_local_covers_20261018'
MODEL='fresh_current324_full_model_20261018'
L='fresh_cross_vertex_current290_local_covers_20261018'
U='fresh_cross_vertex_historical32_local_covers_20261018'
C='fresh_new5121_raw4_two_face_local_covers_20261018'


def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def domains():
    raw={}
    def put(tag,row,groups=None,group_source=None):
        v=row['vertex']
        if groups is None:groups=row.get('corner_orbits') or row.get('union',{}).get('corner_orbits')
        if groups is None:
            groups=raw[v]['groups'];group_source=raw[v]['group_source']
        raw[v]=dict(groups=groups,domain=row,geometry_source=tag,group_source=group_source or tag)
    for tag in ('fresh_current_historical_weighted_corners_20261017','fresh_unowned_cover_pair_geometry_20261017'):
        for r in read(tag)['forest_records']:put(tag,r)
    domain_tag='fresh_historical22_corner_domain_20261017'
    groups={r['vertex']:r['corner_orbits'] for r in read(domain_tag)['forest_records']}
    tag='fresh_historical22_corner_cones_20261017'
    for r in read(tag)['forest_records']:put(tag,r,groups[r['vertex']],domain_tag)
    tag='fresh_historical22_eight_pair_geometry_20261017';put(tag,read(tag))
    for tag in ('fresh_historical22_raw88_raw25_pair_geometry_20261017',
                'fresh_raw125_raw205_actual_corner_frontier_20261017',
                'fresh_290_class_shortlist_pair_geometry_20261017',
                'fresh_2599_raw110_raw24_pair_geometry_20261017'):
        for r in read(tag)['forest_records']:put(tag,r)
    assert len(raw)==571
    for r in read(L)['forest_records']:
        d=raw[r['vertex']];d.update(source=L,record=r)
        assert r['count_vectors_sha256']==jhash(d['domain']['maximal_count_vectors'])
        assert r['orbit_coefficients']==[g['positive_face_coefficients'] for g in d['groups']]
    for tag in (U,C):
        for r in read(tag)['forest_records']:
            raw[r['vertex']]=dict(source=tag,record=r,domain=r,groups=r['domain']['corner_orbits'],
                                 geometry_source=tag,group_source=tag)
    assert len(raw)==600
    for ref in read(C)['inherited_forests']:
        d=raw[ref['vertex']]
        assert d['source']==ref['source'] and jhash(d['record'])==ref['record_sha256']
    return raw


def main():
    saved=read(TAG);assert saved['status']=='PASS_NEW5073_COMPLETE_CURRENT324_EXACT_LOCAL_COVERS'
    for p,h in saved['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    model=read(MODEL);dual=next(c for c in model['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert saved['face_multipliers']==dual['face_multipliers'] and saved['initial_vertex_multipliers']==dual['vertex_multipliers']
    alpha={i:F(x) for i,x in dual['face_multipliers']};mu={i:F(x) for i,x in dual['vertex_multipliers']};newmu=dict(mu)
    assert len(alpha)==5073 and len(mu)==1095
    dd=domains();records={r['vertex']:r for r in saved['forest_records']}
    assert len(records)==len(saved['forest_records'])==600 and set(records)==set(dd)
    inequalities=terms=zero=orbits=0;exports={};faces=set()
    for v,d in sorted(dd.items()):
        r=records[v];groups=d['groups'];vectors=d['domain']['maximal_count_vectors'];reps=d['domain']['literal_representatives']
        assert (r['source'],r['geometry_source'],r['group_source'])==(d['source'],d['geometry_source'],d['group_source'])
        assert r['source_record_sha256']==jhash(d['record']) and r['count_vectors_sha256']==jhash(vectors)
        assert r['literal_representatives_sha256']==jhash(reps) and r['complete_count_vectors']==len(vectors)
        coefficients=[g['positive_face_coefficients'] for g in groups]
        assert r['orbit_coefficients']==coefficients
        demand=[sum((alpha.get(f,F())*m for f,m in co),F()) for co in coefficients]
        assert r['demands']==list(map(str,demand))
        owner={i:g['index'] for g in groups for i in g['literal_stars']}
        assert len(owner)==sum(len(g['literal_stars']) for g in groups) and len(vectors)==len(reps)
        for vec,clique in zip(vectors,reps):
            assert len(clique)==len(set(clique)) and vec==[Counter(owner[i] for i in clique).get(k,0) for k in range(len(groups))]
        for g in groups:
            assert Counter(o['face'] for o in g['origins'])==Counter(dict(g['positive_face_coefficients']))
            faces.update(f for f,m in g['positive_face_coefficients'])
        c=r['certificate'];w=list(map(F,c['weights']));mix=[(j,F(x)) for j,x in c['mixture']]
        assert len(w)==len(groups) and all(x>=0 for x in w)
        assert len({j for j,x in mix})==len(mix) and all(0<=j<len(vectors) and x>0 for j,x in mix)
        assert all(sum((a*b for a,b in zip(vec,w)),F())<=1 for vec in vectors)
        assert all(sum((vectors[j][i]*x for j,x in mix),F())>=dem for i,dem in enumerate(demand))
        mass=sum((x for j,x in mix),F());load=sum((a*b for a,b in zip(w,demand)),F())
        assert mass==load==F(c['optimal_mass'])==F(r['exact_local_optimum'])
        old=mu.get(v,F());newmu[v]=max(old,mass)
        assert F(r['baseline_vertex_mass'])==old and F(r['repaired_vertex_mass'])==newmu[v]
        assert F(r['required_increase'])==newmu[v]-old
        co=Counter()
        for weight,ls in zip(w,coefficients):
            for f,m in ls:co[f]+=weight*m
        if mass>old:
            exports[v]=dict(vertex=v,vertex_cost='1',positive_face_coefficients=[[f,str(x)] for f,x in sorted(co.items()) if x],
                orbit_weights=c['weights'],exact_load=str(mass),baseline_vertex_mass=str(old),source_geometry=d['geometry_source'])
        inequalities+=len(vectors);terms+=len(mix);orbits+=len(groups);zero+=not any(demand)
    assert saved['represented_classes']==sorted(faces) and len(faces)==324 and len(faces&set(alpha))==240
    assert len(saved['columns'])==len(exports)
    for index,col in enumerate(saved['columns']):
        assert col['index']==index and {k:v for k,v in col.items() if k!='index'}==exports.pop(col['vertex'])
    assert not exports
    newpairs=[[v,str(x)] for v,x in sorted(newmu.items()) if x]
    assert saved['vertex_multipliers']==newpairs and all(newmu.get(v,F())>=x for v,x in mu.items())
    loads=[];multi=[]
    for col in model['prepared']['inherited_several_vertex_columns']+model['prepared']['new_columns']:
        load=sum((alpha.get(i,F())*m for i,m in col['integer_face_coefficients']),F())
        before=sum((mu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        after=sum((newmu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        assert load<=before<=after
        row=dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after))
        loads.append(row)
        if len(col['integer_vertex_capacities'])>1:multi.append(row)
    assert len(loads)==412 and len(multi)==14
    assert saved['checked_recent_model_column_loads']==loads and saved['preserved_multi_owner_columns']==multi
    num=F(dual['objective_numerator']);den=sum(newmu.values(),F())
    assert saved['objective_numerator']==str(num) and saved['normalizer']==str(den) and saved['bound']==str(num/den)
    assert saved['initial_bound']==dual['bound'] and num/sum(mu.values(),F())==F(dual['bound'])
    assert saved['counts']==dict(represented_classes=324,positive_represented_classes=240,forests=600,
        exact_new_demand_local_optima=600,zero_demand_forests=zero,nonzero_demand_forests=600-zero,
        new_charge_columns=len(saved['columns']),raised_vertex_masses=len(saved['columns']),
        complete_count_vectors=inequalities,orbit_coordinates=orbits,mixture_terms=terms)
    assert saved['generation_local_LPs']==600-zero
    assert all(saved[k]==0 for k in ('replay_local_LPs','global_optimizer_invocations','new_physical_pair_tests','clique_enumerations','old_mass_optima_reused'))
    print('PASS_STDLIB_CURRENT324_NEW5073_SAVED',dict(forests=600,inequalities=inequalities,mixture_terms=terms,
        violated_columns=len(saved['columns']),recent_column_loads=len(loads),multi_owner=len(multi),bound=saved['bound']))


if __name__=='__main__':main()
