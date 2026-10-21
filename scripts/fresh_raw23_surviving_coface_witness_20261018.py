#!/usr/bin/env python3
"""Exact full-matrix raw23-only surviving Q witness; no optimization.

Every positive simplex and edge multiplier is retained. All Q multipliers
except the unambiguous raw23 source rows are deleted, including every merged
row. Feasibility is checked on the full frozen coface matrix, not inferred
from monotonicity of signed face coefficients.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_raw357_coface_full_model_20261018 as C
import fresh_all_marker_families_model_support_20261018 as S

ROOT=C.ROOT;P=C.P;M=C.M
TAG='fresh_raw23_surviving_coface_witness_20261018'

def run():
    paths={f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',f'results/{C.TAG}.md',
        f'scripts/{S.TAG}.py',f'results/{S.TAG}.json',f'results/{C.A.TAG}.md'}
    sources={p:M.digest(p)for p in sorted(paths)}
    current=M.read(C.TAG);support=M.read(S.TAG)
    assert current['status']=='PASS_RAW357_COFACE_EXTENSION_EXACT_INTERVAL'
    assert support['status']=='PASS_SAVED_ALL_MARKER_FAMILIES_MODEL_SUPPORT'
    for doc in(current,support):
        for p,h in doc['sources'].items():assert M.digest(p)==h,p
    # Reconstruct every inherited coefficient and all16 coface additions.
    d,prepared,_,_,_,_,_=C.build()
    assert prepared==current['prepared']
    prior=next(c for c in current['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    primal=next(c for c in current['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    P.verify(d,prior);P.verify(d,primal)
    oldalpha={f:F(x)for f,x in prior['face_multipliers']}
    rows={r['face']:r for r in support['positive_faces']}
    assert len(rows)==len(support['positive_faces']) and set(oldalpha)<=set(rows)
    kept=[];removed=[]
    for f,a in sorted(oldalpha.items()):
        r=rows[f];face=d['faces'][f]
        assert F(r['mass'])==a and r['required']==face['required']
        if r['kind']=='Q':
            assert face['kind']=='canonical_dotted_quad_aggregate'
            keep=r.get('raw_family')==23
            if keep:
                assert r['provenance_kind']in('explicit','unique_saved_context')
                assert r['canonical_core']==19 and r['required']==10
        else:
            assert r['kind']==face['kind'] and r['kind']in('maximal_simplex_aggregate','maximal_edge_incidence')
            keep=True
        (kept if keep else removed).append(r)
    keptq=[r for r in kept if r['kind']=='Q']
    assert len(keptq)==50 and all(r['raw_family']==23 for r in keptq)
    alpha={r['face']:oldalpha[r['face']]for r in kept}
    mu={v:F(x)for v,x in prior['vertex_multipliers']}
    loads=P.loads(d,alpha);costs=P.cost_rows(d)
    violations=[[j,str(loads[j]),str(sum((a*mu.get(v,F())for v,a in cc),F()))]
        for j,cc in enumerate(costs)if loads[j]>sum((a*mu.get(v,F())for v,a in cc),F())]
    numerator=sum((d['faces'][f]['required']*a for f,a in alpha.items()),F())
    denominator=sum(mu.values(),F())
    deleted=sum((F(r['mass'])*r['required']for r in removed),F())
    assert numerator==F(prior['objective_numerator'])-deleted
    cert=dict(kind='EXACT_PARTIAL_MODEL_DUAL',face_multipliers=[[f,str(x)]for f,x in sorted(alpha.items())],
        vertex_multipliers=prior['vertex_multipliers'],objective_numerator=str(numerator),
        normalizer=str(denominator),bound=str(numerator/denominator),excludes_target12=numerator>12*denominator,
        recovery='Retain all positive simplex/edge masses and only unambiguous raw23 Q masses. Check every full-matrix column, including signed zero-budget and positive-budget columns.')
    if not violations:P.verify(d,cert)
    counts=Counter(r['kind']for r in kept)
    neg=[[r['face'],sum(a<0 for j,a in d['faces'][r['face']]['coefficients'])]for r in removed]
    return dict(status='PASS_RAW23_ONLY_FULL_MATRIX_LOWER_WITNESS'if not violations else'PASS_RAW23_ONLY_INFEASIBLE_DELETION_DIAGNOSTIC',
        sources=sources,prepared=prepared,
        inherited_certificate_bound=prior['bound'],candidate=cert,
        certificates=[primal,cert]if not violations else[primal],full_matrix_column_checks=len(loads),
        deletion_violations=violations,kept_positive_counts=dict(counts),
        kept_Q_sources=keptq,deleted_Q_sources=removed,
        deleted_required_mass=str(deleted),deleted_rows_with_negative_coefficients=[r for r in neg if r[1]],
        retained_simplex_required_mass=str(sum((F(r['mass'])*r['required']for r in kept if r['kind']=='maximal_simplex_aggregate'),F())),
        retained_raw23_required_mass=str(sum((F(r['mass'])*r['required']for r in keptq),F())),
        retained_edge_required_mass='0',vertex_masses_unchanged=True,
        preserved_multi_owner_columns=14,optimizer_invocations=0,model_changed=False,
        scope='A feasible exact dual of the unchanged full finite coface matrix if all column checks pass. An above12 witness obstructs target12 in this finite matrix, not a geometric theorem or an optimum. No general new-coface or cross-simplex constraint is asserted here.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['kept_positive_counts'],out['candidate']['bound'],
        'violations',len(out['deletion_violations']),'negative_deleted_rows',len(out['deleted_rows_with_negative_coefficients']),flush=True)
