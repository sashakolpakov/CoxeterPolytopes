#!/usr/bin/env python3
"""Incremental frozen archive, exact incidence repairs, and restricted cusp bound.

Old dependency graphs are inherited only after every frozen byte digest is
rechecked. New JSON manifests and record references are inspected recursively.
No producer verifier, optimizer, or historical million-case scan is rerun.
"""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/forced_quadrilaterals_and_affine_carriers_20261017.md'
PREVIOUS='results/spectator_field_gates_and_boundary_replacements_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_spectator_field_gates_and_boundary_replacements_20261017.py'
OUT=ROOT/'results/forced_quadrilaterals_and_affine_carriers_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
    'fresh_sqrt5_core23_25_forced_quad_20261017',
    'fresh_sqrt5_core23_25_incidence_repair_20261017',
    'fresh_sqrt5_core23_25_spectator_incidence_20261017',
    'fresh_large_affine_factor_quad_payment_20261017',
    'fresh_cusp_link_rank_monotonicity_20261017',
    'fresh_private_affine_tree_exclusion_20261017',
    'fresh_affine_tree_carrier_exclusion_20261017',
    'fresh_large_affine_white_contact_counterexamples_20261017',
    'fresh_large_affine_white_continuation_scout_20261017',
    'fresh_cusp_minimum_factor_rank17_dimension29_20261017',
    'fresh_cusp_minimum_factor_rank17_dimension29_20261017_review',
    'fresh_nonarithmetic_release30_66_hinge_20261017',
    'fresh_nonarithmetic_alternative66_endcaps_20261017',
    'fresh_nonarithmetic_alternative66_release59_20261017',
    'fresh_nonarithmetic_alternative66_chain_patch_20261017',
)
EXTRA_REPORTS=(
    'results/fresh_sqrt5_core23_25_forced_quad_20261017_review.md',
    'results/fresh_private_affine_cycle_exclusion_20261017.md',
    'results/fresh_affine_cycle_carrier_exclusion_20261017.md',
)

@lru_cache(None)
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())

def masses(rows,limit):
    d={i:F(a)for i,a in rows}
    assert len(d)==len(rows)and all(0<=i<limit and a>=0 for i,a in d.items())
    return d

def small_checks():
    old=read('fresh_sqrt5_a3_d5_cap_joint_ports_20261017')
    gate=read('fresh_sqrt5_core23_25_forced_quad_20261017')
    first=read('fresh_sqrt5_core23_25_incidence_repair_20261017')
    final=read('fresh_sqrt5_core23_25_spectator_incidence_20261017')
    assert old['face_multipliers']==first['face_multipliers']==final['face_multipliers']
    alpha=masses(old['face_multipliers'],82806)
    mu0,mu1,mu2=[masses(d['vertex_multipliers'],6558)for d in(old,first,final)]
    def changed(a,b):return{i for i in set(a)|set(b)if a.get(i,F())!=b.get(i,F())}
    assert changed(mu0,mu1)=={5572}and changed(mu1,mu2)=={5413}
    assert mu1[5572]-mu0.get(5572,F())==F(3363761,3000000000)
    assert mu2[5413]-mu1.get(5413,F())==F(50071,5000000000)
    assert all(mu2.get(i,F())>=mu1.get(i,F())>=mu0.get(i,F())for i in set(mu0)|set(mu1)|set(mu2))
    assert first['vertex_updates']==[[5572,str(mu0.get(5572,F())),str(mu1[5572])]]
    assert final['vertex_updates']==[[5413,str(mu1.get(5413,F())),str(mu2[5413])]]
    def load(c):
        return sum((F(a)*alpha.get(i,F())for i,a in c['positive_faces']),F())-sum((F(a)*alpha.get(i,F())for i,a in c['negative_faces']),F())
    def capacity(c,mu):return sum((F(h)*mu.get(v,F())for v,h in c['vertex_capacity']),F())
    row=gate['inequality'];assert row==first['new_incidence_row']
    assert row['positive_faces']==[[26410,2],[26880,1]]and row['negative_faces']==[[580,2]]and row['vertex_capacity']==[[5572,3]]
    assert load(row)==F(3363761,1000000000)==F(row['face_lower_bound'])
    assert capacity(row,mu0)==F(row['vertex_upper_bound'])==0
    assert load(row)==capacity(row,mu1)<=capacity(row,mu2)
    cs=final['constraints'];assert len(cs)==52
    assert len({c['vertex_capacity'][0][0]for c in cs})==52
    violations=[];ks=set();raw_counts=Counter()
    for ix,c in enumerate(cs):
        assert c['index']==ix and c['simplex_dimension']==4
        k=c['A1_count'];ks.add(k);raw_counts[tuple(c['selected_raw_sources'])]+=1
        assert k==c['spectator_forest'].count('A1')>=2
        target=list(c['spectator_forest']);target.remove('A1');target.remove('A1');target+=['H3','H3'];target.sort()
        assert target==c['target_forest']
        h=target.count('H3');assert c['vertex_capacity'][0][1]==h>=2
        assert c['negative_faces']==[[c['simplex_face'],k]]
        profiles={q['raw_source']:q for q in c['actual_profiles']}
        assert sorted(profiles)==c['selected_raw_sources']
        expected=[[profiles[t]['face'],k if t==23 else 1]for t in(23,25)if t in profiles]
        assert c['positive_faces']==expected
        l=load(c);before=capacity(c,mu1)
        assert l==F(c['input_load'])and before==F(c['input_capacity'])and l-before==F(c['gap'])
        if l>before:violations.append(ix)
        assert l<=capacity(c,mu2)
    assert raw_counts=={(23,):39,(23,25):13}
    assert violations==final['violation_indices']and len(violations)==1
    bad=cs[violations[0]];assert bad['vertex_capacity']==[[5413,2]]and F(bad['gap'])==F(50071,2500000000)
    count_tables=[]
    for k in sorted(ks):
        for na in(0,1):
            for nb in range(k-na+1):
                assert na*nb>=k*na+nb-k
                count_tables.append([k,na,nb,na*nb,k*na+nb-k])
    assert count_tables==final['integer_count_tables']
    objective=F(old['objective_numerator']);assert objective==F(3842394855727,309375000000)
    ratios=[]
    for d,mu in((old,mu0),(first,mu1),(final,mu2)):
        assert F(d['objective_numerator'])==objective
        assert F(d['normalizer'])==sum(mu.values(),F())>0
        assert F(d['bound'])==objective/F(d['normalizer'])>12
        ratios.append(F(d['bound']))
    assert ratios[0]>ratios[1]>ratios[2]
    assert F(first['input_bound'])==ratios[0]and F(final['input_bound'])==ratios[1]
    for d in(first,final):assert d['optimizer_invocations']==d['repeated_large_scans']==0 and d['face_masses_unchanged']is True
    assert final['alpha_changed']is False
    theorem=read('fresh_cusp_minimum_factor_rank17_dimension29_20261017')
    review=read('fresh_cusp_minimum_factor_rank17_dimension29_20261017_review')
    assert theorem['status']=='PASS_ALL_RANK_MINIMUM_CUSP_FACTOR_RANK17_DIMENSION29_COMPOSITION'
    assert review['status']=='PASS_INDEPENDENT_RATIONAL_SCALARS_AND_MATHEMATICAL_COMPOSITION_REVIEW'
    assert theorem['dimension_bound']==29 and theorem['affine_order_exclusion_through']==17
    assert theorem['arithmetic_hypothesis']is False and theorem['factor_count_limit']is None
    assert theorem['enumeration_cutoff_used']is False and theorem['unrestricted163_tail_used']is False
    assert F(theorem['curvature_increment'])==F(7,76)
    # Small independent rational composition checks; the all-rank proof is in
    # the separately reviewed nonnegative polynomial identities.
    for n in range(30,81):
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>=3*n and 5*(n-1)>0
    assert 4563*2-8075==1051>0
    for stem in STEMS:
        d=read(stem)
        assert d['status'].startswith('PASS')or(stem.startswith('fresh_nonarithmetic_')and d['status'].startswith('BOUNDED_NUMERICAL'))
        if stem.startswith('fresh_nonarithmetic_'):assert d['improvements']==[]
    deletion=read('fresh_nonarithmetic_alternative66_release59_20261017')['exact_baseline_deletion59_cost']
    assert deletion['lost_finite_count']==7 and deletion['compact_faces_by_dimension']=={'2':6,'3':1}
    chain=read('fresh_nonarithmetic_alternative66_chain_patch_20261017')
    assert list(map(len,chain['groups']))==[16,27]and sum(chain['pair_relations'].values())==432
    return dict(scope='Fraction-only incidence and monotone-budget composition, theorem-status/scalar consistency and construction-scope counts; no rederived geometry or repeated large scans.',
                incidence_constraints=52,literal_profile_maps=65,changed_vertex_indices=[5572,5413],
                added_vertex_masses=[str(F(3363761,3000000000)),str(F(50071,5000000000))],
                alpha_unchanged=True,all_vertex_masses_monotone=True,
                exact_input_bound=str(ratios[0]),exact_first_repair_bound=str(ratios[1]),exact_final_bound=str(ratios[2]),
                target12_still_excluded_in_partial_model=True,optimizer_invocations=0,repeated_large_scans=0,
                unrestricted_dimension_bound_unchanged=163,minimum_cusp_factor_rank17_dimension_bound=29,
                factor_count_limit=None,arithmetic_hypothesis=False,new_finite_volume_example=None,
                construction_control=[63,199,10])

def run():
    prior=json.loads((ROOT/PREVIOUS).read_text());frozen={};preserved={}
    for section in('source_hashes','top_artifact_hashes','report_hashes'):
        for name,h in prior[section].items():
            assert digest(ROOT/name)==h,('frozen predecessor',name)
            assert name not in frozen or frozen[name]==h
            frozen[name]=h
        preserved[section]=len(prior[section])
    assert digest(ROOT/PREVIOUS_AUDITOR)==prior['auditor_sha256']
    frozen[PREVIOUS_AUDITOR]=prior['auditor_sha256'];frozen[PREVIOUS]=digest(ROOT/PREVIOUS)
    print('ARCHIVE_FROZEN_BYTES_PASS',len(frozen),flush=True)
    sources=dict(frozen);sources[HELPER]=digest(ROOT/HELPER)
    top=[ROOT/f'results/{s}.json'for s in STEMS]
    pending=list(top);seen=set();new_manifests=[];bibliography=list(prior['inherited_bibliographic_entries'])
    references=records_count=0;skipped=set();new_jsons=[]
    while pending:
        path=pending.pop();name=str(path.relative_to(ROOT))
        if name in seen:continue
        seen.add(name)
        if name in frozen:
            assert digest(path)==frozen[name];skipped.add(name);continue
        doc=json.loads(path.read_text());new_jsons.append(name)
        if not isinstance(doc,dict):continue
        bound=[]
        if 'database_sha256'in doc and name not in A.HISTORICAL_DATABASE_METADATA:
            assert isinstance(doc.get('database'),str);bound.append((doc['database'],doc['database_sha256']))
        mdocs=[doc]+([doc['prepared']]if isinstance(doc.get('prepared'),dict)else[])
        for md in mdocs:
            for field in('sources','source_hashes'):
                for source,value in md.get(field,{}).items():
                    if isinstance(value,dict):
                        assert field=='sources'and'title'in value
                        bibliography.append(dict(document=name,citation=source))
                        for key,h in value.items():
                            if key.endswith('_sha256'):
                                q=value.get(key[:-7])
                                if isinstance(q,str)and(ROOT/q).is_file():bound.append((q,h))
                    else:bound.append((source,value))
        if 'source_sha256'in doc:
            self_source='scripts/'+path.stem+'.py'
            assert(ROOT/self_source).is_file();bound.append((self_source,doc['source_sha256']))
        records=list(A.record_refs(doc));records_count+=len(records)
        for source,h in bound+records:
            assert re.fullmatch('[0-9a-f]{64}',h),(name,source)
            source=A.EXTERNAL_SOURCE_ALIASES.get(source,source);target=ROOT/source
            assert target.is_file()and digest(target)==h,(name,source)
            assert source not in sources or sources[source]==h
            sources[source]=h;references+=1
            if target.suffix=='.json':pending.append(target)
        if bound or records:new_manifests.append(dict(path=name,status=doc.get('status'),source_references=len(bound),record_references=len(records)))
    for path in top:sources[str(path.relative_to(ROOT))]=digest(path)
    print('ARCHIVE_NEW_GRAPH_PASS',len(new_manifests),'inherited_json_skips',len(skipped),flush=True)
    checks=small_checks();print('ARCHIVE_FRACTION_CHECKS_PASS',checks['incidence_constraints'],flush=True)
    reports=[ROOT/CHECKPOINT]+[ROOT/f'results/{s}.md'for s in STEMS]+[ROOT/p for p in EXTRA_REPORTS]
    links=[]
    for report in reports:
        assert report.is_file(),report
        for target in re.findall(r'\]\(([^)]+)\)',report.read_text()):
            if '://'in target or target.startswith('#'):continue
            resolved=(report.parent/target.split('#')[0]).resolve()
            assert resolved.is_file()or resolved==OUT,(report,target)
            links.append(dict(report=str(report.relative_to(ROOT)),target=target))
    syntax=[]
    for name in sorted(set(sources)-set(frozen)):
        if name.endswith('.py'):ast.parse((ROOT/name).read_text());syntax.append(name)
    ast.parse(Path(__file__).read_text());syntax.append(str(Path(__file__).relative_to(ROOT)))
    all_manifests={r['path']:r for r in prior['manifests']}
    for r in new_manifests:assert r['path']not in all_manifests;all_manifests[r['path']]=r
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_FORCED_QUADS_AND_AFFINE_CARRIERS',
                scope='Archive integrity, exact small rational repairs and reviewed theorem composition metadata; not formalization or a repeated historical mathematical verification.',
                checkpoint=CHECKPOINT,previous_archive=PREVIOUS,previous_archive_sha256=digest(ROOT/PREVIOUS),
                previous_frozen_sections_preserved=preserved,previous_auditor_preserved=True,
                frozen_json_graph_policy='Rehash all inherited files first; do not reparse their already verified dependency graph.',
                inherited_json_documents_not_reparsed=sorted(skipped),new_json_documents_parsed=sorted(new_jsons),
                historical_nonbinding_database_metadata=sorted(A.HISTORICAL_DATABASE_METADATA),
                recovered_historical_source_aliases=A.EXTERNAL_SOURCE_ALIASES,
                recovered_provenance_sha256=digest(ROOT/A.RECOVERED_PROVENANCE),
                inherited_bibliographic_entries=bibliography,small_certificate_checks=checks,
                report_count=len(reports),local_link_count=len(links),
                source_manifest_count=len(all_manifests),new_source_manifest_count=len(new_manifests),
                source_reference_count=prior['source_reference_count']+references,
                new_source_reference_count=references,
                nested_record_reference_count=prior['nested_record_reference_count']+records_count,
                unique_source_count=len(sources),manifests=sorted(all_manifests.values(),key=lambda r:r['path']),
                source_hashes=dict(sorted(sources.items())),
                source_hash_scope='Cumulative predecessor files and new recursive sources; entrypoint snapshots are historical.',
                top_artifact_hashes={str(p.relative_to(ROOT)):digest(p)for p in top},
                report_hashes={str(p.relative_to(ROOT)):digest(p)for p in reports},local_links=links,
                new_python_syntax_checks=syntax,
                entrypoint_hashes={p:digest(ROOT/p)for p in('README.md','RESEARCH_PLAN.md')},
                entrypoint_scope='Current entrypoints only; prior snapshots remain historical.',auditor_sha256=digest(Path(__file__)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--small-only',action='store_true');args=ap.parse_args()
    if args.small_only:
        print(json.dumps(small_checks(),indent=2));raise SystemExit(0)
    d=run()
    if args.verify:assert d==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],{k:d[k]for k in('report_count','local_link_count','source_manifest_count','new_source_manifest_count','unique_source_count')},flush=True)
