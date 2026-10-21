#!/usr/bin/env python3
"""Stdlib-only accounting for all242 marker-family exports and global bounds.

No producer/solver imports, native geometry or full old matrix reconstruction.
The separate full-matrix replay proves old-column matching/dual conditions,
all primal face requirements and the face-weighted dual numerator.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import gcd,lcm
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw76_complete_family_full_model_20261018'
PARENT='fresh_current327_full_model_20261018'
FAMILY='fresh_raw76_extra_vertex_complete_family_20261018'
LOCAL='fresh_current324_new5073_local_covers_20261018'
PATCHES=('fresh_new5073_raw76_high_mass_local_covers_20261018','fresh_new5073_raw76_third_class_local_covers_20261018')

@lru_cache(None)
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def pairs(xs):
    out=Counter()
    for i,x in xs:out[i]+=F(x)
    return {i:x for i,x in out.items()if x}
def canonical(face,cost):
    den=lcm(*(x.denominator for x in list(face.values())+list(cost.values())))
    div=gcd(*(int(x*den)for x in list(face.values())+list(cost.values())))
    factor=F(den,div)
    return (tuple(sorted((i,int(x*factor))for i,x in face.items())),
            tuple(sorted((i,int(x*factor))for i,x in cost.items()))),factor
@lru_cache(None)
def rows(tag):return {r['vertex']:r for r in read(tag)['forest_records']}


def main():
    d=read(TAG);p=d['prepared'];old=read(PARENT);fam=read(FAMILY)
    assert d['status'].startswith('PASS_EXACT')and d['solver_invocations']==1
    assert p['inherited_model_sha256']==old['prepared']['model_sha256']
    assert p['provenance_to_face']==old['prepared']['provenance_to_face']
    inherited=old['prepared']['inherited_multi_vertex_columns']
    assert p['inherited_multi_vertex_columns']==inherited and len(inherited)==14
    checked={}
    for bindings in(p['sources'],fam['sources'],old['prepared']['sources'],d.get('certification_sources',{})):
        for path,h in bindings.items():
            if path not in checked:checked[path]=sha256((ROOT/path).read_bytes()).hexdigest()
            assert checked[path]==h,path
    dual0=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha=pairs(dual0['face_multipliers']);mu0=pairs(dual0['vertex_multipliers'])
    fixed=fam['exported']['repair'];mu=pairs(fixed['vertex_multipliers'])
    assert pairs(fixed['face_multipliers'])==alpha and all(mu.get(v,F())>=x for v,x in mu0.items())
    assert F(fixed['bound'])==F(fixed['objective_numerator'])/sum(mu.values(),F())>12
    vectors={}
    for r in read(LOCAL)['forest_records']:
        vv=rows(r['geometry_source'])[r['vertex']]['maximal_count_vectors']
        assert jhash(vv)==r['count_vectors_sha256'];vectors[r['vertex']]=vv
    for tag in PATCHES:
        for r in read(tag)['forest_records']:vectors[r['vertex']]=r['maximal_count_vectors']
    records={r['vertex']:r for r in fam['forest_records']}
    assert len(records)==235 and len(fam['exported']['columns'])==242
    raw={};keys={};ineq=0
    for i,col in enumerate(fam['exported']['columns']):
        v=col['vertex'];r=records[v];groups=r['old_groups']
        owner={s:g['index']for g in groups for s in g['literal_stars']}
        assert set(owner)==set(range(len(r['old_nodes'])))
        assert all(Counter(o['face']for o in g['origins'])==Counter(dict(g['positive_face_coefficients']))for g in groups)
        w0=list(map(F,col['old_orbit_weights']));w=list(map(F,col['retained_orbit_weights']))
        assert len(w0)==len(w)==len(groups)and all(x>=0 for x in w0)
        vv=vectors[v]if groups else [[]]
        assert len(vv)==col['old_vector_count']and jhash(vv)==col['old_vectors_sha256']
        for a in vv:
            assert len(a)==len(w0)and all(isinstance(n,int)and n>=0 for n in a)
            assert sum((x*n for x,n in zip(w0,a)),F())<=1;ineq+=1
        omit=sorted({owner[j]for j,k in r['marker_ordinary_edges']})
        assert omit==col['omitted_orbits']and w==[F()if j in omit else x for j,x in enumerate(w0)]
        assert not omit or(v==1304 and omit==[0]and groups[0]['positive_face_coefficients']==[[60808,1]])
        types=Counter(r['forest']);cap=min(types['A1']-1,types['H3'])
        assert cap==r['marker_capacity']and F(col['marker_weight'])==F(1,cap)
        assert len(r['markers'])==types['A1']*(types['A1']-1)*types['H3']
        face=Counter()
        for x,g in zip(w,groups):
            for f,n in g['positive_face_coefficients']:face[f]+=x*n
        face[r['face']]+=F(1,cap);face={f:x for f,x in face.items()if x}
        assert pairs(col['positive_face_coefficients'])==face and col['vertex_cost']=='1'
        key,factor=canonical(face,{v:F(1)})
        assert key not in keys;keys[key]=i;raw[i]=(face,{v:F(1)},r,col,factor)
    assert ineq==1048 and len(raw)==len(keys)==242
    seen=set();usedkeys=set();exports=[]
    assert len(p['new_columns'])+len(p['reused_columns'])==242
    for status,recs in(('new',p['new_columns']),('reused',p['reused_columns'])):
        for i,r in enumerate(recs):
            face=pairs(r['integer_face_coefficients']);cost=pairs(r['integer_vertex_capacities'])
            assert len(face)==len(r['integer_face_coefficients'])and len(cost)==len(r['integer_vertex_capacities'])==1
            assert all(0<=f<82806 and x>0 and x.denominator==1 for f,x in face.items())
            assert all(0<=v<6558 and x>0 and x.denominator==1 for v,x in cost.items())
            key=tuple(sorted((f,int(x))for f,x in face.items())),tuple(sorted((v,int(x))for v,x in cost.items()))
            assert key in keys and key not in usedkeys;usedkeys.add(key)
            assert gcd(*(int(x)for x in list(face.values())+list(cost.values())))==1
            load=sum((alpha.get(f,F())*x for f,x in face.items()),F())
            oldcost=sum((mu0.get(v,F())*x for v,x in cost.items()),F())
            assert load<=sum((mu.get(v,F())*x for v,x in cost.items()),F())
            if status=='new':assert r['variable']==439856+i
            else:
                assert r['matching_variables']and all(0<=j<439856 for j in r['matching_variables'])
                assert load<=oldcost
            assert len(r['origins'])==1
            o=r['origins'][0];ci=o['column_index'];assert ci not in seen;seen.add(ci)
            of,oc,rr,col,factor=raw[ci]
            assert o['package']==FAMILY and o['vertex']==rr['vertex']and keys[key]==ci
            assert o['record_sha256']==jhash(rr)and o['column_sha256']==jhash(col)
            assert {f:x*factor for f,x in of.items()}==face and {v:x*factor for v,x in oc.items()}==cost
            assert factor==F(o['scaling']['denominator'],o['scaling']['integer_gcd'])
            assert F(o['old_dual_load'])==load and F(o['old_dual_cost'])==oldcost and F(o['old_dual_violation'])==load-oldcost
            exports.append((r.get('variable'),face,cost,load>oldcost,status))
    assert seen==set(raw)and len(usedkeys)==242
    assert sum(z[3]for z in exports)==10 and all(z[4]=='new'for z in exports if z[3])
    counts=p['counts'];assert counts['inherited_variables']==439856
    assert counts['variables']==439856+len(p['new_columns'])
    assert counts['exported_family_columns']==242 and counts['appended_family_columns']==len(p['new_columns'])and counts['reused_family_columns']==len(p['reused_columns'])
    assert counts['zero_budget_signed_variables']==96176 and counts['inherited_vertex_cost_signed_variables']==52
    assert counts['faces']==82806 and counts['vertex_forests']==6558 and counts['preserved_multi_owner_columns']==14
    assert p['family_audit']['old_charge_vector_inequalities']==1048 and p['family_audit']['old_dual_violations']==10
    assert p['family_audit']['fixed_mass_exported_only_ratio']==fixed['bound']
    primal=next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    dual=next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    a,m=pairs(dual['face_multipliers']),pairs(dual['vertex_multipliers'])
    assert len(a)==len(dual['face_multipliers'])and len(m)==len(dual['vertex_multipliers'])
    assert all(0<=f<82806 and x>0 for f,x in a.items())and all(0<=v<6558 and x>0 for v,x in m.items())
    for j,face,cost,_,_ in exports:
        assert sum((a.get(f,F())*x for f,x in face.items()),F())<=sum((m.get(v,F())*x for v,x in cost.items()),F())
    for c in inherited:
        assert sum((a.get(f,F())*x for f,x in c['integer_face_coefficients']),F())<=sum((m.get(v,F())*x for v,x in c['integer_vertex_capacities']),F())
    x=pairs(primal['rounded_variable_numerators']);den=primal['integer_denominator']
    assert den>0 and len(x)==len(primal['rounded_variable_numerators'])
    assert all(0<=j<counts['variables']and z>=0 and z.denominator==1 for j,z in x.items())
    assert F(primal['repair_delta'])>=0 and all(0<=j<439856 for j in primal['repair_variable_indices'])
    budgets=pairs(primal['vertex_budgets']);contribution=Counter()
    for j,face,cost,_,status in exports:
        if status=='new':
            for v,z in cost.items():contribution[v]+=z*x.get(j,F())/den
    assert all(z<=budgets.get(v,F())for v,z in contribution.items())and max(budgets.values())==F(primal['maximum'])
    lower=F(dual['objective_numerator'])/F(dual['normalizer']);upper=F(primal['maximum'])
    assert sum(m.values(),F())==F(dual['normalizer'])>0 and lower==F(dual['bound'])<=upper
    assert dual['excludes_target12']==(lower>12)and primal['meets_target12']==(upper<=12)
    assert d['status']=='PASS_EXACT_RAW76_FAMILY_FALLBACK_INTERVAL_OPTIMUM_UNRESOLVED'
    assert d['optimization_resolved']is False and d['changed_face_masses']is False
    proposal=read(TAG+'_proposal')
    assert proposal['status']=='UNRESOLVED_SINGLE_RAW76_COMPLETE_FAMILY_SOLVE'
    assert proposal['solver_invocations']==1 and proposal['certificates']==[]
    assert proposal['prepared']==p and proposal['solver']==d['solver']
    assert d['solver']['status']==1 and not d['solver']['success']and d['solver']['objective']is None
    assert primal==next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    assert a==alpha and m==mu and lower==F(fixed['bound'])
    assert d['exact_interval']==dict(lower=str(lower),upper=str(upper),gap=str(upper-lower))
    represented=set(read(LOCAL)['represented_classes'])
    represented.update(r['face']for r in read(PATCHES[0])['selected_faces']);represented.add(read(PATCHES[1])['selected_face']);assert len(represented)==327
    markerfaces={r['face']for r in fam['family']};assert len(markerfaces)==235
    migration=dict(old_positive_faces=len(alpha),new_positive_faces=len(a),
        newly_positive=len(set(a)-set(alpha)),newly_zero=len(set(alpha)-set(a)),
        newly_positive_outside327=len(set(a)-set(alpha)-represented),ordinary327_positive=len(set(a)&represented),
        marker235_positive=len(set(a)&markerfaces),marker235_mass=str(sum((a.get(f,F())for f in markerfaces),F())),
        new_columns_used_by_primal=sum(x.get(j,F())>0 for j,_,_,_,s in exports if s=='new'))
    print('PASS_RAW76_COMPLETE_FAMILY_MODEL_SAVED_CHECK',dict(charge_columns=242,old_vector_inequalities=ineq,
        appended=len(p['new_columns']),reused=len(p['reused_columns']),two_owner_columns=14,
        lower=str(lower),upper=str(upper),exact_gap=str(upper-lower),
        lower_decimal=float(lower),upper_decimal=float(upper),support_migration=migration),flush=True)


if __name__=='__main__':main()
