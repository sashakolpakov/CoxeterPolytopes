#!/usr/bin/env python3
"""Exact support split of the new actual-corner model dual; no optimizer."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import fresh_actual_corner_enlarged_model_20261017 as E
import fresh_dominant_cap_compatibility_20261017 as P

ROOT=E.ROOT
TAG='fresh_actual_corner_enlarged_model_20261017_support'


def run():
    model=E.read(E.TAG);d,prepared=E.build()
    assert json.loads(E.dump(prepared))==model['prepared']
    full=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    E.verify(d,full)
    alpha={fi:F(a)for fi,a in full['face_multipliers']};mu={v:F(a)for v,a in full['vertex_multipliers']}
    pmap=dict(prepared['provenance_to_face']);inverse={fi:j for j,fi in pmap.items()}
    kinds=Counter(d['faces'][fi]['kind']for fi in alpha)
    qrows={fi for fi in alpha if d['faces'][fi]['kind']=='canonical_dotted_quad_aggregate'}
    selected=qrows&set(inverse);unselected=qrows-set(inverse);simplex=set(alpha)-qrows
    assert not set(inverse)&simplex
    old=E.read(E.WITNESS);oldS={f['face']for f in old['positive_simplex_features']}
    input_alpha={fi:a for fi,a in alpha.items()if fi not in unselected}
    restricted=E.ENGINE.dual_certificate(d,dict(input_alpha));E.verify(d,restricted)
    newa={fi:F(a)for fi,a in restricted['face_multipliers']}
    assert not unselected&set(newa) and selected<=set(newa)
    assert all(newa.get(fi,F())==a for fi,a in alpha.items()if fi in selected)
    restored=sorted(set(newa)-set(input_alpha))
    assert not restored # Removing Q masses cannot require a larger S mass.
    assert newa==input_alpha
    actual=P.read(P.M.TAGS['actual']);quads=P.read(P.M.TAGS['quad'])
    byraw=defaultdict(lambda:dict(faces=[],mass=F(),objective=F()))
    for fi in sorted(selected):
        pr=actual['provenance_masses'][inverse[fi]]
        raw=quads['canonical_cores'][pr['canonical_prism_core']]['source_core']
        byraw[raw]['faces'].append(fi);byraw[raw]['mass']+=alpha[fi]
        byraw[raw]['objective']+=d['faces'][fi]['required']*alpha[fi]
    rawrows=[dict(raw_source=raw,faces=x['faces'],total_face_mass=str(x['mass']),
        objective_contribution=str(x['objective']))for raw,x in sorted(byraw.items(),key=lambda z:(-z[1]['objective'],z[0]))]
    loads=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in d['faces'][fi]['coefficients']:loads[j]+=a*m
    added=[]
    for b in prepared['added_column_blocks']:
        tight=positive=0;maximum=None
        for j in range(b['start'],b['start']+b['count']):
            c=d['variables'][j];cap=F()if c['vertex']is None else c['cost']*mu.get(c['vertex'],F())
            assert loads[j]<=cap
            positive+=loads[j]>0;tight+=loads[j]==cap
            excess=loads[j]-cap;maximum=excess if maximum is None else max(maximum,excess)
        added.append(dict(kind=b['kind'],columns=b['count'],positive_load_columns=positive,
            tight_columns=tight,largest_signed_slack=str(maximum)))
    old23={c['face']for c in E.read(E.FAMILY)['selected_faces']}
    new23=set(byraw.get(23,{}).get('faces',[]))
    paths=[Path(__file__),Path(E.__file__),ROOT/f'results/{E.TAG}.json',
        ROOT/f'results/{E.WITNESS}.json',ROOT/f'results/{E.FAMILY}.json',
        ROOT/f'results/{P.M.TAGS["actual"]}.json',ROOT/f'results/{P.M.TAGS["quad"]}.json']
    result=dict(status='PASS_EXACT_ACTUAL_CORNER_NEW_DUAL_SUPPORT_RESTRICTION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        input_model_sha256=prepared['model_sha256'],optimizer_invocations=0,
        full_dual_bound=full['bound'],counts=dict(positive_faces=len(alpha),positive_vertex_masses=len(mu),
            face_kinds=dict(kinds),simplex_faces=len(simplex),selected_actual_Q=len(selected),
            unselected_merged_Q=len(unselected),new_S_outside_prior_fixed_cover=len(simplex-oldS),
            old_raw23_positive_classes=len(old23),new_raw23_positive_classes=len(new23)),
        unselected_Q_rows=sorted(unselected),new_simplex_rows_outside_old_cover=sorted(simplex-oldS),
        selected_Q_by_raw_core=rawrows,new_column_loads=added,
        restricted_actual_Q_test=dict(name='all_positive_S_and_selected_actual_Q',
            dropped_unselected_Q_rows=sorted(unselected),restored_S_outside_seed_domain=restored,
            face_masses_preserved_on_retained_rows=True,repaired_dual=restricted),
        scope='Exact restriction of the new finite-model dual to all its S rows and only common, explicitly split actual Q provenances. No arbitrary provenance is assigned to a discarded merged row. Per-vertex loads are recomputed against every finite exported column. Neither old local covers nor old complete joint-family scans are asserted at these new masses.')
    print(result['status'],result['counts'],'full',float(F(full['bound'])),
          'restricted',float(F(restricted['bound'])),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(E.dump(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(E.dump(d)+'\n')
