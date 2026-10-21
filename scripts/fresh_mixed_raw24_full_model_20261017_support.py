#!/usr/bin/env python3
"""Exact saved-dual support/coverage diagnostic; no optimizer or new geometry."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_mixed_raw24_full_model_20261017 as N
import fresh_dominant_cap_compatibility_20261017 as P
ROOT=N.ROOT
TAG='fresh_mixed_raw24_full_model_20261017_support'

def run():
    saved=N.read(N.TAG);d,prepared=N.build()
    assert json.loads(N.dump(prepared))==saved['prepared']
    for c in saved['certificates']:N.verify(d,c)
    cert=next(c for c in saved['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={i:F(a)for i,a in cert['face_multipliers']}
    old=N.read('fresh_mixed_raw24_actual_corners_20261017')
    oldalpha={i:F(a)for i,a in old['face_multipliers']}
    local=N.read('fresh_mixed_raw39_full_model_20261017_local_cones')
    joint=N.read('fresh_mixed_raw39_full_model_20261017_joint_ports')
    sf={r['face']for r in local['positive_simplex_features']}
    jq={r['face']for r in joint['actual_quad_features']}
    mf={r['face']for r in old['selected_faces']}
    previous=N.read(P.FORESTS);faces=previous['faces']
    actual=N.read(P.M.TAGS['actual']);quad=N.read(P.M.TAGS['quad']);raw=N.read(P.M.TAGS['raw'])
    pmap=dict(prepared['provenance_to_face']);assert len(pmap)==len(set(pmap.values()))
    inverse={fi:j for j,fi in pmap.items()}
    Q={fi for fi in alpha if fi in inverse or (fi<len(faces)and faces[fi]['kind']=='canonical_dotted_quad_aggregate')}
    S=set(alpha)-Q
    assert all('simplex'in d['faces'][fi]['kind']or'edge'in d['faces'][fi]['kind']for fi in S)
    groups=defaultdict(list);qrows=[]
    for fi in sorted(Q):
        row=dict(face=fi,mass=str(alpha[fi]),required=d['faces'][fi]['required'],
            demand_mass=str(alpha[fi]*d['faces'][fi]['required']),in_mixed_export=fi in mf,in_joint_export=fi in jq)
        if fi in inverse:
            j=inverse[fi];pr=actual['provenance_masses'][j];ci=pr['canonical_prism_core']
            family=quad['canonical_cores'][ci]['source_core']
            row.update(provenance_kind='explicit',provenance=j,source_context=pr['source_context'],
                canonical_core=ci,raw_family=family,source_context_count=1)
        else:
            contexts=faces[fi]['source_contexts']
            fs=sorted({quad['canonical_cores'][ci]['source_core']for ctx,ci in contexts})
            row.update(source_context_count=len(contexts),raw_families=fs)
            if len(contexts)==1:
                ctx,ci=contexts[0];row.update(provenance_kind='unique_saved_context',
                    source_context=ctx,canonical_core=ci,raw_family=fs[0])
            else:row.update(provenance_kind='merged_unexpanded')
        qrows.append(row)
        if not row['in_mixed_export'] and row['provenance_kind']!='merged_unexpanded':groups[row['raw_family']].append(row)
    cohorts=[]
    for family,rows in groups.items():
        cohorts.append(dict(raw_family=family,classes=len(rows),explicit=sum(r['provenance_kind']=='explicit'for r in rows),
            unique=sum(r['provenance_kind']=='unique_saved_context'for r in rows),
            mass=str(sum((F(r['mass'])for r in rows),F())),
            demand_mass=str(sum((F(r['demand_mass'])for r in rows),F())),
            top_faces=sorted(rows,key=lambda r:(-F(r['demand_mass']),r['face']))[:10]))
    cohorts.sort(key=lambda r:(-F(r['demand_mass']),r['raw_family']))
    for cohort in cohorts:
        family=cohort['raw_family'];g=P.V.raw_graph(raw['records'][family])
        labels=sorted({e[0]for row in g for e in row})
        dotted=sorted({e for row in g for e in row if e[0]==0})
        cc=next(c for c in quad['canonical_cores']if c['source_core']==family)
        corners=P.X.quad_features(dict(quadrilateral_black_partitions=cc['quad_partitions']))
        roots=[]
        for e in dotted:
            tau=P.V.squared(e);r=P.M.ksqrt(tau)
            roots.append(dict(square=tau.out(),positive_root_in_K=None if r is None else r.out()))
        cohort.update(core_order=len(g),finite_labels=[m for m in labels if m],dotted_parameters=roots,
            canonical_core=cc['index'],corners_per_class=len(corners),actual_corner_occurrences=len(corners)*cohort['classes'])
    mu={i:F(a)for i,a in cert['vertex_multipliers']}
    capstats=[]
    for tag,field,coeffield,expected in N.INPUTS:
        package=N.read(tag);loads=[]
        for c in package[field]:
            load=sum((F(m)*alpha.get(fi,F())for fi,m in c[coeffield]),F())
            capacity=F(c['vertex_cost'])*mu.get(c['vertex'],F())
            assert load<=capacity
            loads.append(dict(index=c['index'],vertex=c['vertex'],load=str(load),capacity=str(capacity),slack=str(capacity-load)))
        capstats.append(dict(source=tag,columns=len(loads),positive_loads=sum(F(r['load'])>0 for r in loads),
            exactly_tight_positive=sum(F(r['load'])>0 and F(r['slack'])==0 for r in loads),
            loads=loads))
    restrictions=[]
    exactQ={r['face']for r in qrows if r['provenance_kind']!='merged_unexpanded'}
    for name,keepQ in [('already_mixed_exported_Q',Q&mf),('all_explicit_or_unique_Q',exactQ)]:
        keep=S|keepQ;seed={fi:alpha[fi]for fi in keep}
        restricted=N.ENGINE.dual_certificate(d,dict(seed));N.verify(d,restricted)
        assert {fi:F(a)for fi,a in restricted['face_multipliers']}==seed
        assert all(F(a)<=mu.get(v,F())for v,a in restricted['vertex_multipliers'])
        restrictions.append(dict(name=name,kept_S=len(S),kept_Q=len(keepQ),dropped_Q=len(Q-keepQ),
            outside_S_restored=False,alpha_unchanged_on_retained_rows=True,certificate=restricted))
    counts=dict(positive_faces=len(alpha),positive_S=len(S),positive_Q=len(Q),positive_vertices=len(mu),
        newly_positive_faces=len(set(alpha)-set(oldalpha)),dropped_faces=len(set(oldalpha)-set(alpha)),
        changed_common_faces=sum(alpha[i]!=oldalpha[i]for i in set(alpha)&set(oldalpha)),
        S_inside_old_complete_feature_domain=len(S&sf),S_outside_old_complete_feature_domain=len(S-sf),
        Q_inside_mixed_export=len(Q&mf),Q_outside_mixed_export=len(Q-mf),
        Q_inside_old_joint_export=len(Q&jq),Q_outside_old_joint_export=len(Q-jq),
        Q_explicit=sum(r['provenance_kind']=='explicit'for r in qrows),
        Q_unique_saved_context=sum(r['provenance_kind']=='unique_saved_context'for r in qrows),
        Q_merged_unexpanded=sum(r['provenance_kind']=='merged_unexpanded'for r in qrows))
    paths=[Path(__file__),Path(N.__file__),ROOT/f'results/{N.TAG}.json',
        ROOT/f'results/{P.FORESTS}.json',ROOT/f'results/{P.M.TAGS["actual"]}.json',ROOT/f'results/{P.M.TAGS["quad"]}.json',ROOT/f'results/{P.M.TAGS["raw"]}.json',Path(P.__file__),Path(P.V.__file__),Path(P.M.__file__),Path(P.X.__file__)]
    paths +=[ROOT/f'results/{tag}.json'for tag,field,coeffield,n in N.INPUTS]
    result=dict(status='PASS_NEW_MIXED_RAW24_DUAL_SUPPORT_AND_EXPORTED_LOADS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        model_sha256=prepared['model_sha256'],bound=cert['bound'],counts=counts,
        face_multipliers=cert['face_multipliers'],vertex_multipliers=cert['vertex_multipliers'],
        face_kind_counts=dict(Counter(d['faces'][fi]['kind']for fi in alpha)),
        simplex_rows_outside_old_features=sorted(S-sf),quad_rows=qrows,uncovered_exact_provenance_cohorts=cohorts,
        appended_export_loads=capstats,exact_restrictions=restrictions,optimizer_invocations=0,
        scope='The exact finite-model certificates and every appended exported load are rechecked. This diagnostic classifies new support only; no old-alpha full local cover or implicit joint scan is reused, no merged row is given an arbitrary actual representative, and no global geometry or dimension theorem is asserted.')
    print('SUPPORT',counts,'BOUND',float(F(cert['bound'])),flush=True)
    print('TOP_UNCOVERED',[(r['raw_family'],r['classes'],float(F(r['demand_mass'])))for r in cohorts[:8]],flush=True)
    print('EXACT_RESTRICTIONS',[(r['name'],r['kept_Q'],float(F(r['certificate']['bound'])))for r in restrictions],flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(N.dump(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(N.dump(result)+'\n')
    print(result['status'],flush=True)
