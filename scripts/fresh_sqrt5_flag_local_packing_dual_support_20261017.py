#!/usr/bin/env python3
"""Exact support restrictions of the flag/local successor; no optimizer."""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import fresh_sqrt5_flag_local_packing_full_model_20261017 as N

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_flag_local_packing_dual_support_20261017'
SUBSETS=N.SUBSETS
M=N.M


def run():
    saved=json.loads((ROOT/f'results/{N.TAG}.json').read_text())
    for p,h in saved['prepared']['sources'].items():assert M.R.P.sha(p)==h,p
    d,prepared=N.build()
    assert json.loads(json.dumps(prepared))==saved['prepared']
    original=next(c for c in saved['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    N.verify(d,original)
    old=json.loads((ROOT/f'results/{M.TAG}.json').read_text())
    oldcert=next(c for c in old['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    actual=M.R.P.read(M.A.TAG)
    integer=M.R.P.read(M.U.PACKING)
    subset=M.R.P.read(SUBSETS)
    alpha={i:F(a) for i,a in original['face_multipliers']}
    oldalpha={i:F(a) for i,a in oldcert['face_multipliers']}
    pmap=dict(prepared['provenance_to_face'])
    flag_q=set(pmap.values())
    assert len(flag_q)==len(pmap)==1858
    copies=dict(old['prepared']['old_quad_to_actual_faces'])
    def expand(features):
        return {dest for r in features for dest in copies.get(r['face'],[r['face']])}
    mixed=expand(actual['positive_face_mixed_features'])
    integer_s=expand(integer['positive_simplex_features'])
    block_s={r['face'] for r in subset['positive_simplex_features']}
    assert all(d['faces'][i]['kind']!='canonical_dotted_quad_aggregate' for i in block_s|integer_s)
    assert flag_q<=mixed
    all_q={f['index'] for f in d['faces'] if f['kind']=='canonical_dotted_quad_aggregate'}
    all_s=set(range(len(d['faces'])))-all_q
    assert flag_q<=all_q
    legacy_s=mixed&all_s
    assert legacy_s==integer_s
    complete_export_s=block_s&legacy_s
    latest_export=flag_q|block_s
    legacy_export=flag_q|legacy_s
    complete_export=flag_q|complete_export_s
    groups=Counter();demands=defaultdict(F);rows=[]
    for fi,a in sorted(alpha.items()):
        f=d['faces'][fi]
        if fi in all_q:
            group='Q_selected_actual_flag' if fi in flag_q else 'Q_unselected_merged_no_flag'
        elif fi in complete_export_s:
            group='S_exported_mixed_integer_and_full_block'
        elif fi in block_s:
            group='S_full_block_only'
        elif fi in legacy_s:
            group='S_legacy_mixed_integer_only'
        else:
            group='S_outside_both_exports'
        groups[group]+=1;demands[group]+=a*f['required']
        rows.append(dict(face=fi,kind=f['kind'],group=group,mass=str(a),
            objective_contribution=str(a*f['required']),new_positive_index=fi not in oldalpha,
            has_selected_actual_flag=fi in flag_q,has_exported_mixed=fi in mixed,
            has_exported_legacy_integer=fi in integer_s,has_full_block_feature=fi in block_s))
    assert sum(demands.values(),F(0))==F(original['objective_numerator'])
    cases=[
        ('remove_unselected_Q_keep_all_S',all_s|flag_q),
        ('legacy_export_seed',legacy_export),
        ('latest_block_export_seed',latest_export),
        ('all_feature_export_seed',complete_export),
        ('outside_latest_export_seed',set(range(len(d['faces'])))-latest_export),
    ]
    tests=[]
    for name,domain in cases:
        initial={i:a for i,a in alpha.items() if i in domain}
        cert=M.dual_certificate(d,dict(initial))
        N.verify(d,cert)
        repaired={i:F(a) for i,a in cert['face_multipliers']}
        assert all(repaired.get(i,F(0))>=a for i,a in initial.items())
        assert all(repaired.get(i,F(0))==initial.get(i,F(0)) for i in all_q)
        changed=[i for i,a in repaired.items() if a>initial.get(i,F(0))]
        assert set(changed)<=all_s
        outside=[i for i in repaired if i not in domain]
        assert set(outside)<=all_s
        rec=dict(name=name,seed_positive_faces=len(initial),
            seed_positive_Q=sum(i in all_q for i in initial),seed_positive_S=sum(i in all_s for i in initial),
            deleted_positive_Q=sum(i in all_q and i not in domain for i in alpha),
            increased_or_restored_S=len(changed),restored_S_outside_seed_domain=outside,
            restored_S_outside_latest_block_export=[i for i in repaired if i in all_s and i not in block_s],
            simplex_mass_changes=[[i,str(initial.get(i,F(0))),str(repaired[i])] for i in changed],
            complete_signed_columns_checked=len(d['zeros']),repaired_dual=cert)
        tests.append(rec)
        print('SUPPORT_RESTRICTION',name,len(initial),'OUTSIDE_S',len(outside),
              'BOUND',cert['bound'],float(F(cert['bound'])),flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{N.TAG}.py',f'results/{N.TAG}.json',
        f'scripts/{M.TAG}.py',f'results/{M.TAG}.json',f'results/{M.A.TAG}.json',
        f'results/{M.U.PACKING}.json',f'results/{SUBSETS}.json']
    return dict(status='PASS_EXACT_FLAG_LOCAL_MODEL_SUPPORT_RESTRICTIONS',
        sources={p:M.R.P.sha(p) for p in paths},optimizer_invocations=0,
        model_sha256=prepared['model_sha256'],full_bound=original['bound'],
        counts=dict(positive_faces=len(alpha),positive_vertex_forests=len(original['vertex_multipliers']),
            new_positive_face_indices=len(set(alpha)-set(oldalpha)),
            selected_actual_flag_Q_domain=len(flag_q),full_block_S_domain=len(block_s),
            legacy_export_S_domain=len(legacy_s),common_export_S_domain=len(complete_export_s),
            all_face_obligations=len(d['faces']),signed_columns=len(d['zeros'])),
        positive_support_counts=dict(sorted(groups.items())),
        objective_by_support_group={k:str(a) for k,a in sorted(demands.items())},
        positive_face_support=rows,tests=tests,
        scope='No optimizer or changes to the geometric domain. Every restriction starts with exact masses on explicitly named classes and then uses the full successor dual repair, preserving all94481 signed inequalities including77391 flags. Any forced containing-simplex mass outside the seed domain is listed, never hidden as an exported-only certificate. The unselected Q rows are inherited merged rows without selected actual-provenance splitting; assigning them flag coefficients would require consistent full provenance expansion, not independent minima. Below12 restrictions do not establish primal feasibility or an ambient dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(result['status'],result['counts'],result['positive_support_counts'],flush=True)
