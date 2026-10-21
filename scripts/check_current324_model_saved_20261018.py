#!/usr/bin/env python3
"""Stdlib check of saved local charges, new columns and certificate accounting.

No producer imports, inherited full-matrix reconstruction, solver, graph or
native geometry replay. The separate full-model replay verifies old columns,
primal face obligations, and matching to reused inherited variables.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_current324_full_model_20261018'
PARENT='fresh_cross_vertex_current_full_model_20261018'
INPUTS=('fresh_cross_vertex_current290_local_covers_20261018',
        'fresh_cross_vertex_historical32_local_covers_20261018',
        'fresh_new5121_raw4_two_face_local_covers_20261018')

@lru_cache(None)
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def pairs(xs):
    co=Counter()
    for i,x in xs:co[i]+=F(x)
    return {i:x for i,x in co.items() if x}

def main():
    d=read(TAG);p=d['prepared'];parent=read(PARENT)
    assert d['status'].startswith('PASS_EXACT') and d['solver_invocations']==1
    assert p['inherited_model_sha256']==parent['prepared']['model_sha256'] and p['provenance_to_face']==parent['prepared']['provenance_to_face']
    assert p['inherited_several_vertex_columns']==parent['prepared']['new_columns']
    for path,h in p['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    alpha=pairs(read(INPUTS[-1])['face_multipliers']);mu=pairs(read(INPUTS[-1])['vertex_multipliers'])
    previous_dual=next(c for c in parent['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert alpha==pairs(previous_dual['face_multipliers']);oldmu=pairs(previous_dual['vertex_multipliers'])
    assert all(mu.get(v,F())>=x for v,x in oldmu.items())
    raw={};zero=inequalities=mixtureterms=0
    for tag in INPUTS:
        doc=read(tag);assert pairs(doc['face_multipliers'])==alpha
        for index,row in enumerate(doc['forest_records']):
            v=row['vertex'];cert=row['certificate']
            if tag==INPUTS[0]:
                src=read(row['source']);rr=next(x for x in src.get('forest_records',[src]) if x['vertex']==v)
                vectors=rr['maximal_count_vectors'];coeffs=row['orbit_coefficients']
                assert jhash(vectors)==row['count_vectors_sha256']
            else:vectors=row['maximal_count_vectors'];coeffs=[g['positive_face_coefficients'] for g in row['domain']['corner_orbits']]
            w=list(map(F,cert['weights']));dem=[sum((alpha.get(f,F())*F(m) for f,m in co),F()) for co in coeffs]
            assert len(w)==len(coeffs) and all(x>=0 for x in w)
            assert list(map(str,dem))==(row['demands'] if tag==INPUTS[0] else row['domain']['demands'])
            assert all(len(vector)==len(w) for vector in vectors)
            for vector in vectors:assert sum((x*y for x,y in zip(w,vector)),F())<=1;inequalities+=1
            coverage=[F() for _ in w];mass=F()
            for vi,amount in cert['mixture']:
                amount=F(amount);assert 0<=vi<len(vectors) and amount>=0;mass+=amount;mixtureterms+=1
                for i,n in enumerate(vectors[vi]):coverage[i]+=n*amount
            assert all(x>=y for x,y in zip(coverage,dem))
            load=sum((x*y for x,y in zip(w,dem)),F())
            assert mass==load==F(cert['optimal_mass'])==F(row['exact_local_optimum'])<=mu.get(v,F())
            face=Counter()
            for x,co in zip(w,coeffs):
                for f,n in co:face[f]+=x*F(n)
            face={f:x for f,x in face.items() if x}
            if not face:zero+=1;continue
            raw[tag,index]=(face,{v:F(1)},row)
    assert len(raw)==679 and zero==59 and inequalities==3273
    declared={};seen=set();keys=set();added=[];reused=[]
    for kind,records in [('new',p['new_columns']),('reused',p['reused_local_columns'])]:
        for offset,row in enumerate(records):
            face=pairs(row['integer_face_coefficients']);cost=pairs(row['integer_vertex_capacities'])
            assert len(face)==len(row['integer_face_coefficients']) and len(cost)==len(row['integer_vertex_capacities'])==1
            assert all(0<=f<82806 and x>0 and x.denominator==1 for f,x in face.items())
            assert all(0<=v<6558 and x>0 and x.denominator==1 for v,x in cost.items())
            assert gcd(*(int(x) for x in list(face.values())+list(cost.values())))==1
            key=tuple(sorted(face.items())),tuple(sorted(cost.items()));assert key not in keys;keys.add(key)
            assert row['origins']
            for origin in row['origins']:
                q=origin['package'],origin['record_index'];assert q not in seen;seen.add(q)
                of,oc,rr=raw[q];assert origin['vertex']==rr['vertex'] and origin['record_sha256']==jhash(rr)
                assert set(of)==set(face) and set(oc)==set(cost)
                f=min(face);ratio=face[f]/of[f]
                assert ratio>0 and all(face[i]==ratio*x for i,x in of.items()) and all(cost[i]==ratio*x for i,x in oc.items())
                assert ratio==F(origin['scaling']['denominator'],origin['scaling']['integer_gcd'])
                declared[q]=key
            assert sum((alpha.get(f,F())*x for f,x in face.items()),F())<=sum((mu.get(v,F())*x for v,x in cost.items()),F())
            if kind=='new':assert row['variable']==439619+offset;added.append((row['variable'],face,cost,row))
            else:
                assert row['inherited_variables'] and len(set(row['inherited_variables']))==len(row['inherited_variables'])
                assert all(0<=j<439619 for j in row['inherited_variables']);reused.append(row)
    assert seen==set(raw) and len(keys)==647
    exported=set()
    for tag in INPUTS:
        for col in read(tag)['columns']:
            matches=[q for q in raw if q[0]==tag and raw[q][2]['vertex']==col['vertex']];assert len(matches)==1
            key=declared[matches[0]];exported.add(key)
            face,cost=map(dict,key)
            assert sum((alpha.get(f,F())*x for f,x in face.items()),F())>sum((oldmu.get(v,F())*x for v,x in cost.items()),F())
    assert len(exported)==47
    assert sum(row['violates_previous_optimized_dual'] for j,f,c,row in added)==47
    for j,face,cost,row in added:
        actual_violation=sum((alpha.get(f,F())*z for f,z in face.items()),F())>sum((oldmu.get(v,F())*z for v,z in cost.items()),F())
        assert row['violates_previous_optimized_dual']==actual_violation==((tuple(sorted(face.items())),tuple(sorted(cost.items()))) in exported)
    counts=p['counts'];assert counts['variables']==439619+len(added) and counts['appended_local_columns']==len(added)
    assert counts['reused_local_columns']==len(reused) and counts['preserved_multi_owner_columns']==14
    assert counts['zero_budget_signed_variables']==96176 and counts['inherited_vertex_cost_signed_variables']==52
    primal=next(c for c in d['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    dual=next(c for c in d['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    a,m=pairs(dual['face_multipliers']),pairs(dual['vertex_multipliers'])
    assert len(a)==len(dual['face_multipliers']) and len(m)==len(dual['vertex_multipliers'])
    assert all(0<=f<82806 and x>0 for f,x in a.items()) and all(0<=v<6558 and x>0 for v,x in m.items())
    for j,face,cost,row in added:
        assert sum((a.get(f,F())*x for f,x in face.items()),F())<=sum((m.get(v,F())*x for v,x in cost.items()),F())
    two=0
    for row in parent['prepared']['new_columns']:
        if len(row['integer_vertex_capacities'])>1:
            two+=1
            assert sum((a.get(f,F())*x for f,x in row['integer_face_coefficients']),F())<=sum((m.get(v,F())*x for v,x in row['integer_vertex_capacities']),F())
    assert two==14
    x=pairs(primal['rounded_variable_numerators']);den=primal['integer_denominator']
    assert den>0 and len(x)==len(primal['rounded_variable_numerators'])
    assert all(0<=j<counts['variables'] and z>=0 and z.denominator==1 for j,z in x.items())
    assert F(primal['repair_delta'])>=0 and all(0<=j<439619 for j in primal['repair_variable_indices'])
    budgets=pairs(primal['vertex_budgets']);contribution=Counter()
    for j,face,cost,row in added:
        for v,z in cost.items():contribution[v]+=z*x.get(j,F())/den
    assert all(y<=budgets.get(v,F()) for v,y in contribution.items())
    assert max(budgets.values())==F(primal['maximum'])
    assert sum(m.values(),F())==F(dual['normalizer'])>0
    lower=F(dual['objective_numerator'])/F(dual['normalizer']);upper=F(primal['maximum'])
    assert lower==F(dual['bound'])<=upper
    assert dual['excludes_target12']==(lower>12) and primal['meets_target12']==(upper<=12)
    migration=dict(old_positive_faces=len(alpha),new_positive_faces=len(a),
        added_positive_faces=len(set(a)-set(alpha)),removed_positive_faces=len(set(alpha)-set(a)),
        new_columns_used_by_primal=sum(x.get(j,F())>0 for j,face,cost,row in added))
    print('PASS_CURRENT324_SAVED_LOCAL_COLUMNS_AND_CERTIFICATES',dict(local_certificates=738,charge_inequalities=inequalities,mixture_terms=mixtureterms,origins=len(seen),distinct_charges=len(keys),appended=len(added),reused=len(reused),two_owner_columns=two,lower=float(lower),upper=float(upper),support_migration=migration))

if __name__=='__main__':main()
