#!/usr/bin/env python3
"""Independent stdlib saved arithmetic checker for the changed family model.

Checks complete export/origin accounting, arbitrary-mass charge coefficients,
primitive normalization, saved dual loads on new and14 multi-owner columns,
and exact interval arithmetic. Frozen geometry/envelope and old graph/vector
completeness are inherited. The producer's separate full-matrix replay checks
all old signed inequalities, weighted objective and every primal obligation.
"""
import json
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from math import gcd,lcm
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_all_marker_families_dual_continuation_20261018'
FIRST='fresh_all_marker_families_full_model_20261018'
FAMILIES=(
 ('fresh_raw76_extra_vertex_complete_family_20261018',242),
 ('fresh_raw56_raw57_extra_vertex_family_20261018',1212),
 ('fresh_raw56_raw57_raw76_joint_marker_family_20261018',1430),
 ('fresh_raw108_complete_marker_family_20261018',415),
 ('fresh_raw112_complete_marker_family_20261018',529))
LOCAL='fresh_current324_new5073_local_covers_20261018'
PATCH='fresh_new5073_raw76_high_mass_local_covers_20261018'
THIRD='fresh_new5073_raw76_third_class_local_covers_20261018'


def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def hashfile(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def hashdata(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def primitive(faces,vertex):
    co=Counter()
    for f,a in faces:assert Q(a)>0;co[f]+=Q(a)
    denomin=lcm(*(x.denominator for x in co.values()))
    numbers=[int(co[f]*denomin)for f in sorted(co)]+[denomin]
    common=gcd(*numbers)
    return(tuple((f,int(co[f]*denomin)//common)for f in sorted(co)),((vertex,denomin//common),)),denomin,common


def run():
    data=read(TAG);p=data['prepared'];assert data['solver_invocations']==2
    assert data['status'].startswith(('PASS_EXACT_','PARTIAL_'))
    for path,h in p['sources'].items():assert hashfile(path)==h,path
    for path,h in data['sources'].items():assert hashfile(path)==h,path
    first=read(FIRST);assert first['prepared']==p and first['solver_invocations']==1
    assert first['solver']['model_status_string']=='Time limit reached'
    assert data['solver_history']==[first['solver'],data['solver']]
    assert hashfile(first['numerical_snapshot'])==first['numerical_snapshot_sha256']
    assert hashfile(data['numerical_snapshot'])==data['numerical_snapshot_sha256']
    docs={s:read(s)for s,n in FAMILIES};rows={s:{r['vertex']:r for r in d['forest_records']}for s,d in docs.items()}
    assert all(len(rows[s])==len(d['forest_records'])for s,d in docs.items())
    vectors={};cache={}
    for tag in(LOCAL,PATCH,THIRD):
        for r in read(tag)['forest_records']:
            if 'maximal_count_vectors'in r:vv=r['maximal_count_vectors']
            else:
                geo=r['geometry_source']
                if geo not in cache:cache[geo]={z['vertex']:z for z in read(geo)['forest_records']}
                vv=cache[geo][r['vertex']]['maximal_count_vectors']
                assert hashdata(vv)==r['count_vectors_sha256']and len(vv)==r['complete_count_vectors']
            vectors[r['vertex']]=vv
    assert len(vectors)==600
    exports={};inequalities=marker_checks=0
    for stage,(s,n)in enumerate(FAMILIES):
        d=docs[s];assert len(d['exported']['columns'])==n
        quotas={x['vertex']:x for x in d['exported'].get('target_allowances',[])}
        for index,c in enumerate(d['exported']['columns']):
            v=c['vertex'];r=rows[s][v]
            if stage==2:
                g=rows[FAMILIES[1][0]].get(v,rows[FAMILIES[0][0]].get(v));assert g is not None
                assert hashdata(g['old_groups'])==r['old_groups_sha256']
            else:g=r
            groups=g['old_groups'];original=list(map(Q,c['old_orbit_weights']));retained=list(map(Q,c['retained_orbit_weights']))
            assert len(original)==len(retained)==len(groups)and all(x>=0 for x in original)
            omitted=set(c['omitted_orbits']);assert len(omitted)==len(c['omitted_orbits'])and all(0<=i<len(groups)for i in omitted)
            assert retained==[Q()if i in omitted else x for i,x in enumerate(original)]
            vv=vectors[v]if groups else[[]]
            assert hashdata(vv)==c['old_vectors_sha256']and len(vv)==c['old_vector_count']
            for q in vv:
                assert len(q)==len(original)and all(isinstance(z,int)and z>=0 for z in q)
                assert sum(z*x for z,x in zip(q,original))<=1;inequalities+=1
            terms=[[76,c['source_face'],c['marker_weight']]]if stage==0 else c['marker_terms']
            weights={k:Q(x)for k,f,x in terms};assert len(weights)==len(terms)and all(x>0 for x in weights.values())
            if stage==0:order=[76];boundary=[[r['marker_capacity']]]
            elif stage==1:
                cap={f['raw_family']:f['marker_capacity']for f in r['marker_families']};order=sorted(cap)
                boundary=[[cap[k]if k==j else 0 for k in order]for j in order]
            else:
                env=r['envelope']if stage==2 else quotas[v]['envelope']
                order=env.get('raw_order',env.get('kind_order'));boundary=env['boundary_counts']
            assert set(weights)<=set(order)
            for q in boundary:
                assert len(q)==len(order)and all(isinstance(z,int)and z>=0 for z in q)
                assert sum(z*weights.get(k,Q())for z,k in zip(q,order))<=1;marker_checks+=1
            cf=Counter()
            for w,z in zip(retained,groups):
                assert Counter(o['face']for o in z['origins'])==Counter(dict(z['positive_face_coefficients']))
                for f,a in z['positive_face_coefficients']:cf[f]+=w*a
            for k,f,x in terms:cf[f]+=Q(x)
            assert [[f,str(x)]for f,x in sorted(cf.items())if x]==c['positive_face_coefficients']and c['vertex_cost']=='1'
            key,den,div=primitive(c['positive_face_coefficients'],v)
            exports[s,index]=dict(key=key,column_sha256=hashdata(c),scaling=dict(denominator=den,integer_gcd=div))
    assert len(exports)==3828
    allcols=p['new_columns']+p['reused_columns'];seen=set();keys=set()
    for c in allcols:
        key=tuple(map(tuple,c['integer_face_coefficients'])),tuple(map(tuple,c['integer_vertex_capacities']))
        assert key not in keys;keys.add(key)
        assert len(key[1])==1 and all(a>0 for f,a in key[0])and key[1][0][1]>0
        for o in c['origins']:
            origin=o['package'],o['column_index'];assert origin not in seen;seen.add(origin)
            expected=exports[origin]
            assert expected==dict(key=key,column_sha256=o['column_sha256'],scaling=o['scaling'])
        if 'matching_variables'in c:
            assert c['matching_variables']and all(0<=j<439856 for j in c['matching_variables'])
            assert all(Q(o['old_dual_violation'])<=0 for o in c['origins'])
    assert seen==set(exports)
    n=p['counts'];assert n['inherited_variables']==439856 and n['zero_budget_signed_variables']==96176
    assert n['inherited_vertex_cost_signed_variables']==52 and n['preserved_multi_owner_columns']==14
    assert n['export_origins']==3828 and n['distinct_family_charges']==len(keys)==len(allcols)
    assert n['appended_family_columns']==len(p['new_columns'])and n['reused_family_columns']==len(p['reused_columns'])
    assert n['variables']==439856+len(p['new_columns'])
    assert[j['variable']for j in p['new_columns']]==list(range(439856,n['variables']))
    audit=p['family_audit'];assert audit['ordinary_vector_inequalities']==inequalities and audit['marker_boundary_inequalities']==marker_checks
    prior=read('fresh_current327_full_model_20261018')['prepared']
    multi=p['inherited_multi_vertex_columns'];assert multi==prior['inherited_multi_vertex_columns']and len(multi)==14
    assert p['provenance_to_face']==prior['provenance_to_face']and p['inherited_model_sha256']==prior['model_sha256']
    duals=[];primals=[]
    for c in data['certificates']:
        if c['kind']=='EXACT_PARTIAL_MODEL_DUAL':
            alpha={f:Q(x)for f,x in c['face_multipliers']};mu={v:Q(x)for v,x in c['vertex_multipliers']}
            assert len(alpha)==len(c['face_multipliers'])and len(mu)==len(c['vertex_multipliers'])
            assert all(0<=f<82806 and x>0 for f,x in alpha.items())and all(0<=v<6558 and x>0 for v,x in mu.items())
            for z in allcols+multi:
                left=sum(alpha.get(f,Q())*a for f,a in z['integer_face_coefficients'])
                right=sum(mu.get(v,Q())*a for v,a in z['integer_vertex_capacities']);assert left<=right
            normal=sum(mu.values(),Q());num=Q(c['objective_numerator']);bound=Q(c['bound'])
            assert normal==Q(c['normalizer'])>0 and num/normal==bound and c['excludes_target12']==(bound>12)
            duals.append(bound)
        else:
            assert c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL'
            values=dict(c['rounded_variable_numerators']);assert len(values)==len(c['rounded_variable_numerators'])
            assert c['integer_denominator']>0 and Q(c['repair_delta'])>=0
            assert all(0<=j<n['variables']and isinstance(x,int)and x>=0 for j,x in values.items())
            maximum=Q(c['maximum']);assert max(Q(x)for v,x in c['vertex_budgets'])==maximum
            assert c['meets_target12']==(maximum<=12);primals.append(maximum)
    assert primals
    if duals:assert max(duals)<=min(primals)
    e=data['exact_interval'];assert Q(e['lower'])==max(duals)and Q(e['upper'])==min(primals)
    assert Q(e['gap'])==min(primals)-max(duals)
    assert e['lower_exceeds12']==(max(duals)>12)and e['upper_at_most12']==(min(primals)<=12)
    print('PASS_ALL_MARKER_FAMILIES_SAVED_CHECK',dict(origins=len(exports),primitive_charges=len(keys),ordinary_inequalities=inequalities,
        marker_inequalities=marker_checks,full_old_matrix_replayed=False,primal_face_obligations_replayed=False),
        'lower',str(max(duals))if duals else None,'upper',str(min(primals)),flush=True)


if __name__=='__main__':run()
