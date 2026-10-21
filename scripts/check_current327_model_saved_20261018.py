#!/usr/bin/env python3
"""Stdlib saved34-charge and new global-certificate accounting.

No producer/solver imports, native geometry, clique census or full old matrix
reconstruction. Old-column/primal-face obligations and full matching remain
obligations of the separate source-bound complete matrix replay.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_current327_full_model_20261018';PARENT='fresh_current324_full_model_20261018'
INPUTS=('fresh_current324_new5073_local_covers_20261018',
        'fresh_new5073_raw76_high_mass_local_covers_20261018',
        'fresh_new5073_raw76_third_class_local_covers_20261018')

@lru_cache(None)
def read(t):return json.loads((ROOT/f'results/{t}.json').read_text())
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def pairs(xs):
    c=Counter()
    for i,x in xs:c[i]+=F(x)
    return {i:x for i,x in c.items() if x}
def row_at(tag,v):
    d=read(tag);rs=d.get('forest_records',[d]);r=[r for r in rs if r.get('vertex')==v]
    assert len(r)==1;return r[0]


def main():
    d=read(TAG);p=d['prepared'];parent=read(PARENT)
    assert d['status'].startswith('PASS_EXACT') and d['solver_invocations']==1
    assert p['inherited_model_sha256']==parent['prepared']['model_sha256']
    assert p['provenance_to_face']==parent['prepared']['provenance_to_face']
    inherited=[r for r in parent['prepared']['inherited_several_vertex_columns'] if len(r['integer_vertex_capacities'])>1]
    assert p['inherited_multi_vertex_columns']==inherited and len(inherited)==14
    for path,h in p['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    final=read(INPUTS[-1]);olddual=next(c for c in parent['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha=pairs(final['face_multipliers']);mu=pairs(final['vertex_multipliers']);oldmu=pairs(olddual['vertex_multipliers'])
    assert alpha==pairs(olddual['face_multipliers']) and F(final['bound'])<12
    assert all(mu.get(v,F())>=x for v,x in oldmu.items())
    raw={};ineq=terms=0;binding_cache={};counts=[]
    for tag,expected_n,expected_c in zip(INPUTS,(600,6,3),(25,6,3)):
        doc=read(tag);assert doc['status'].startswith('PASS') and pairs(doc['face_multipliers'])==alpha
        assert len(doc['forest_records'])==expected_n and len(doc['columns'])==expected_c
        for path,h in doc['sources'].items():
            if path not in binding_cache:binding_cache[path]=sha256((ROOT/path).read_bytes()).hexdigest()
            assert binding_cache[path]==h,path
        for ci,col in enumerate(doc['columns']):
            v=col['vertex'];r=row_at(tag,v);cert=r['certificate']
            if tag==INPUTS[0]:
                rr=row_at(r['geometry_source'],v);vectors=rr['maximal_count_vectors'];co=r['orbit_coefficients'];dem=r['demands']
                assert len(vectors)==r['complete_count_vectors'] and jhash(vectors)==r['count_vectors_sha256']
                assert jhash(row_at(r['source'],v))==r['source_record_sha256']
            else:vectors=r['maximal_count_vectors'];co=[g['positive_face_coefficients'] for g in r['domain']['corner_orbits']];dem=r['domain']['demands']
            w=list(map(F,cert['weights']));demands=[sum((alpha.get(f,F())*F(n) for f,n in cc),F()) for cc in co]
            assert len(w)==len(co) and all(x>=0 for x in w) and list(map(str,demands))==dem
            for a in vectors:
                assert len(a)==len(w) and all(isinstance(n,int) and n>=0 for n in a)
                assert sum((x*n for x,n in zip(w,a)),F())<=1;ineq+=1
            coverage=[F() for _ in w];mass=F()
            for vi,amount in cert['mixture']:
                amount=F(amount);assert 0<=vi<len(vectors) and amount>=0;mass+=amount;terms+=1
                for i,n in enumerate(vectors[vi]):coverage[i]+=n*amount
            assert all(x>=y for x,y in zip(coverage,demands))
            load=sum((x*y for x,y in zip(w,demands)),F())
            assert load==mass==F(cert['optimal_mass'])==F(r['exact_local_optimum'])<=mu.get(v,F())
            face=Counter()
            for x,cc in zip(w,co):
                for f,n in cc:face[f]+=x*F(n)
            face={f:x for f,x in face.items() if x}
            assert face==pairs(col['positive_face_coefficients']) and F(col['vertex_cost'])==1
            raw[tag,ci]=(face,{v:F(1)},r)
        counts.append(dict(package=tag,available_records=expected_n,exported_charges=expected_c))
    assert len(raw)==34 and ineq==530 and terms==244
    assert p['local_audit']['input_counts']==counts
    assert p['local_audit']['charge_inequalities']==ineq and p['local_audit']['mixture_terms']==terms
    assert p['local_audit']['fixed_mass_ratio']==final['bound']
    seen=set();keys=set();added=[]
    assert len(p['new_columns'])==34
    for i,r in enumerate(p['new_columns']):
        face=pairs(r['integer_face_coefficients']);cost=pairs(r['integer_vertex_capacities'])
        assert len(face)==len(r['integer_face_coefficients']) and len(cost)==len(r['integer_vertex_capacities'])==1
        assert all(0<=f<82806 and x>0 and x.denominator==1 for f,x in face.items())
        assert all(0<=v<6558 and x>0 and x.denominator==1 for v,x in cost.items())
        assert gcd(*(int(x) for x in list(face.values())+list(cost.values())))==1
        key=tuple(sorted(face.items())),tuple(sorted(cost.items()));assert key not in keys;keys.add(key)
        assert r['variable']==439822+i and r['origins']
        load=sum((alpha.get(f,F())*x for f,x in face.items()),F())
        oldcost=sum((oldmu.get(v,F())*x for v,x in cost.items()),F())
        assert load>oldcost and load<=sum((mu.get(v,F())*x for v,x in cost.items()),F())
        for o in r['origins']:
            q=o['package'],o['column_index'];assert q not in seen;seen.add(q);of,oc,rr=raw[q]
            assert o['vertex']==rr['vertex'] and o['record_sha256']==jhash(rr)
            assert set(face)==set(of) and set(cost)==set(oc)
            ratio=face[min(face)]/of[min(face)]
            assert ratio>0 and all(face[f]==ratio*x for f,x in of.items()) and all(cost[v]==ratio*x for v,x in oc.items())
            assert ratio==F(o['scaling']['denominator'],o['scaling']['integer_gcd'])
            assert F(o['old_dual_load'])==load and F(o['old_dual_cost'])==oldcost
        added.append((r['variable'],face,cost))
    assert seen==set(raw) and len(keys)==34
    c=p['counts'];assert c['inherited_variables']==439822 and c['variables']==439856
    assert c['appended_local_columns']==c['old_dual_violated_new_columns']==34 and c['reused_local_columns']==0
    assert c['zero_budget_signed_variables']==96176 and c['inherited_vertex_cost_signed_variables']==52
    assert c['preserved_multi_owner_columns']==14 and c['faces']==82806 and c['vertex_forests']==6558
    primal=next(x for x in d['certificates'] if x['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    dual=next(x for x in d['certificates'] if x['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    a,m=pairs(dual['face_multipliers']),pairs(dual['vertex_multipliers'])
    assert len(a)==len(dual['face_multipliers']) and len(m)==len(dual['vertex_multipliers'])
    assert all(0<=f<82806 and x>0 for f,x in a.items()) and all(0<=v<6558 and x>0 for v,x in m.items())
    for j,face,cost in added:
        assert sum((a.get(f,F())*x for f,x in face.items()),F())<=sum((m.get(v,F())*x for v,x in cost.items()),F())
    for r in inherited:
        assert sum((a.get(f,F())*x for f,x in r['integer_face_coefficients']),F())<=sum((m.get(v,F())*x for v,x in r['integer_vertex_capacities']),F())
    x=pairs(primal['rounded_variable_numerators']);den=primal['integer_denominator']
    assert den>0 and len(x)==len(primal['rounded_variable_numerators'])
    assert all(0<=j<c['variables'] and z>=0 and z.denominator==1 for j,z in x.items())
    assert F(primal['repair_delta'])>=0 and all(0<=j<439822 for j in primal['repair_variable_indices'])
    budgets=pairs(primal['vertex_budgets']);contribution=Counter()
    for j,face,cost in added:
        for v,z in cost.items():contribution[v]+=z*x.get(j,F())/den
    assert all(z<=budgets.get(v,F()) for v,z in contribution.items())
    assert max(budgets.values())==F(primal['maximum'])
    lower=F(dual['objective_numerator'])/F(dual['normalizer']);upper=F(primal['maximum'])
    assert sum(m.values(),F())==F(dual['normalizer'])>0 and lower==F(dual['bound'])<=upper
    assert dual['excludes_target12']==(lower>12) and primal['meets_target12']==(upper<=12)
    represented=set(read(INPUTS[0])['represented_classes'])
    represented.update(r['face'] for r in read(INPUTS[1])['selected_faces']);represented.add(final['selected_face']);assert len(represented)==327
    migration=dict(old_positive_faces=len(alpha),new_positive_faces=len(a),
        added_positive_faces=len(set(a)-set(alpha)),removed_positive_faces=len(set(alpha)-set(a)),
        newly_positive_outside327=len(set(a)-set(alpha)-represented),represented327_positive=len(set(a)&represented),
        appended_columns_used_by_primal=sum(x.get(j,F())>0 for j,_,_ in added))
    print('PASS_CURRENT327_SAVED_CHARGES_AND_CERTIFICATES',dict(local_certificates=len(raw),charge_inequalities=ineq,
        mixture_terms=terms,appended=len(added),two_owner_columns=len(inherited),lower=str(lower),upper=str(upper),
        lower_decimal=float(lower),upper_decimal=float(upper),support_migration=migration),flush=True)


if __name__=='__main__':main()
