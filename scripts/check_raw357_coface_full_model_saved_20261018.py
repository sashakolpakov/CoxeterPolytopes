#!/usr/bin/env python3
"""Independent saved arithmetic for the16-column coface extension.

The producer's separate full-matrix replay carries old signed inequalities,
matrix preservation and primal feasibility. This stdlib checker independently
checks all16 saved new loads, exact multiplier deletion and objective arithmetic.
It does not reconstruct old faces, budgets, geometry or a numerical optimizer.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw357_coface_full_model_20261018'
COFACE='fresh_raw357_simplex_coface_transfer_20261018'
MODEL='fresh_all_marker_families_dual_continuation_20261018'


def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())


def main():
    data,coface,old=map(read,(TAG,COFACE,MODEL))
    assert data['status']=='PASS_RAW357_COFACE_EXTENSION_EXACT_INTERVAL'
    for d in(data,coface,old):
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert data['optimizer_invocations']==0 and data['face_masses_free_in_model']is True
    assert data['deletion_violations']==[]and data['full_matrix_column_checks']==442226
    assert len(data['certificates'])==2
    prior=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    primal=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    dual=next(c for c in data['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert next(c for c in data['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')==primal
    assert data['candidate']==dual
    source={r['source_face']for r in coface['family']};assert len(source)==len(coface['family'])==47
    alpha=dict(prior['face_multipliers']);removed=[[f,alpha[f]]for f in sorted(source)if f in alpha]
    assert removed==data['zeroed_source_multipliers']and len(removed)==9
    assert dual['face_multipliers']==[[f,a]for f,a in prior['face_multipliers']if f not in source]
    assert dual['vertex_multipliers']==prior['vertex_multipliers']
    assert len(dict(dual['face_multipliers']))==len(dual['face_multipliers'])
    mass=sum((F(a)for _,a in removed),F());assert mass==F(83284683,10**10)
    assert F(data['deleted_required_mass'])==mass
    numerator=F(prior['objective_numerator'])-mass
    denominator=sum((F(a)for _,a in dual['vertex_multipliers']),F())
    assert numerator==F(dual['objective_numerator'])and denominator==F(prior['normalizer'])==F(dual['normalizer'])
    assert F(dual['bound'])==numerator/denominator and numerator>12*denominator
    assert dual['excludes_target12']is True
    p=data['prepared'];counts=p['counts'];oldcounts=old['prepared']['counts']
    assert p['inherited_model_sha256']==old['prepared']['model_sha256']
    assert p['inherited_counts']==oldcounts
    assert counts['inherited_variables']==442210 and counts['variables']==442226
    assert counts['faces']==oldcounts['faces']==82806 and counts['vertex_forests']==6558
    assert counts['zero_budget_signed_variables']==oldcounts['zero_budget_signed_variables']+16==96192
    assert counts['active_variables']==oldcounts['active_variables']+16
    assert counts['face_matrix_nonzeros']==oldcounts['face_matrix_nonzeros']+32
    assert counts['added_zero_budget_transfers']==16 and counts['absent_simplex_rows_uninstalled']==31
    assert counts['inherited_positive_budget_signed_variables']==52
    assert counts['preserved_multi_owner_columns']==len(p['inherited_multi_vertex_columns'])==14
    assert p['inherited_multi_vertex_columns']==old['prepared']['inherited_multi_vertex_columns']
    assert len(coface['existing_target_transfer_columns'])==len(p['new_columns'])==16
    assert len(coface['missing_simplex_rows'])==31
    a={f:F(x)for f,x in dual['face_multipliers']}
    for k,(c,r)in enumerate(zip(coface['existing_target_transfer_columns'],p['new_columns'])):
        assert r['source_column']==k and r['variable']==442210+k
        assert r['integer_vertex_capacities']==[]==c['integer_vertex_capacities']
        assert r['integer_face_coefficients']==c['integer_face_coefficients']
        assert sum((a.get(f,F())*z for f,z in r['integer_face_coefficients']),F())<=0
    lo,hi=F(dual['bound']),F(primal['maximum']);assert lo<=hi
    assert data['exact_interval']==dict(lower=str(lo),upper=str(hi),gap=str(hi-lo),lower_exceeds12=True)
    print('PASS_RAW357_COFACE_SAVED_NEW_LOADS_AND_INTERVAL',16,9,str(lo),str(hi),flush=True)


if __name__=='__main__':main()
