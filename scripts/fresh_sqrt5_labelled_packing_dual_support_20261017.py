#!/usr/bin/env python3
"""Read-only exact support diagnosis of the one-shot enlarged model dual."""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_labelled_packing_full_model_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_labelled_packing_dual_support_20261017'


def run():
    saved=json.loads((ROOT/f'results/{M.TAG}.json').read_text());d,prepared=M.build()
    assert saved['prepared']==json.loads(json.dumps(prepared))
    cert=next(c for c in saved['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    M.verify_dual(d,cert)
    actual=M.R.P.read(M.A.TAG);packing=M.R.P.read(M.U.PACKING)
    alpha={i:F(a) for i,a in cert['face_multipliers']};mu={i:F(a) for i,a in cert['vertex_multipliers']}
    assert all(a>=0 for a in alpha.values())
    previous={i for i,a in actual['face_multipliers']};previous_vertices={v for v,a in actual['vertex_multipliers']}
    mixed={r['face'] for r in actual['positive_face_mixed_features']};integer={r['face'] for r in packing['positive_simplex_features']}
    inverse={dest:old for old,copies in prepared['old_quad_to_actual_faces'] for dest in copies}
    groups=Counter();demands=defaultdict(F);rows=[];outside={}
    for fi,a in sorted(alpha.items()):
        f=d['faces'][fi];old=inverse.get(fi,fi);quad=f['kind']=='canonical_dotted_quad_aggregate'
        key=('quad' if quad else 'simplex')+('_exported' if old in mixed else '_outside_exported')
        groups[key]+=1;demands[key]+=a*f['required']
        vertices=sorted({d['variables'][j]['vertex'] for j,m in f['coefficients'] if m>0 and d['variables'][j]['vertex'] is not None})
        rows.append(dict(face=fi,inherited_face=old,kind=f['kind'],required=f['required'],mass=str(a),
            objective_contribution=str(a*f['required']),old_positive_class=old in previous,
            has_exported_mixed_feature=old in mixed,has_exported_integer_feature=old in integer,
            vertex_forests=vertices))
        if old not in mixed:outside[fi]=a
    assert sum(demands.values(),F(0))==F(cert['objective_numerator'])
    top=sorted(rows,key=lambda r:(-F(r['objective_contribution']),r['face']))[:30]
    topvertices=[dict(vertex=v,forest=d['budgets'][v]['forest'],mass=str(a),
        normalized_mass=str(a/F(cert['normalizer'])),old_positive_forest=v in previous_vertices)
        for v,a in sorted(mu.items(),key=lambda x:(-x[1],x[0]))[:30]]
    restricted=M.dual_certificate(d,dict(outside));M.verify_dual(d,restricted)
    exported=M.dual_certificate(d,{i:a for i,a in alpha.items() if i not in outside});M.verify_dual(d,exported)
    exported_added=[i for i,a in exported['face_multipliers'] if i not in alpha or i in outside]
    assert all(d['faces'][i]['kind']!='canonical_dotted_quad_aggregate' for i in exported_added)
    restricted_support={i:F(a) for i,a in restricted['face_multipliers']}
    added=[i for i in restricted_support if i not in outside]
    assert all(d['faces'][i]['kind']!='canonical_dotted_quad_aggregate' for i in added)
    paths=[f'scripts/{TAG}.py',f'scripts/{M.TAG}.py',f'results/{M.TAG}.json',
        f'results/{M.A.TAG}.json',f'results/{M.U.PACKING}.json']
    return dict(status='PASS_EXACT_ENLARGED_MODEL_DUAL_SUPPORT',sources={p:M.R.P.sha(p) for p in paths},
        full_dual_bound=cert['bound'],counts=dict(positive_faces=len(alpha),positive_vertex_forests=len(mu),
            newly_positive_forests=len(set(mu)-previous_vertices),
            outside_exported_positive_faces=len(outside),old_exported_positive_faces=len(alpha)-len(outside),
            outside_repair_added_simplex_rows=len(added)),
        face_class_counts=dict(sorted(groups.items())),face_class_objective_contributions={k:str(a) for k,a in sorted(demands.items())},
        positive_face_support=rows,top_objective_faces=top,top_vertex_masses=topvertices,
        outside_exported_repaired_dual=restricted,outside_repair_added_simplex_rows=added,
        exported_only_repaired_dual=exported,exported_repair_added_simplex_rows=exported_added,
        scope='No optimization or face search. The saved exact dual is recounted by its inherited full-star feature class. Two exact restrictions retain respectively the masses outside or inside the previous exported classes, then monotonically restore any containing-simplex masses forced by the existing signed columns. These ratios are witnessed lower bounds, not optima or global realizations. A class with zero new indicator coefficients is intentionally underpaid in this restricted model; the missing features have not been inferred zero geometrically.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(M.R.B.dump(d)+'\n')
    print(d['status'],d['counts'],d['face_class_counts'],'outside',d['outside_exported_repaired_dual']['bound'],'exported',d['exported_only_repaired_dual']['bound'],flush=True)
