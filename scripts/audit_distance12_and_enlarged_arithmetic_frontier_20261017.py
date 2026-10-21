#!/usr/bin/env python3
"""Incremental archive of distance12 and the enlarged arithmetic frontier.

The recursive archive walker is inherited from the frozen previous auditor.
Every inherited file is rehashed before its earlier dependency graph is reused.
New small checks are independent Fraction compositions, not geometry proofs.
"""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import re
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/distance12_and_enlarged_arithmetic_frontier_20261017.md'
PREVIOUS='results/selective13_and_actual_corner_compatibility_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_selective13_and_actual_corner_compatibility_20261017.py'
OUT=ROOT/'results/distance12_and_enlarged_arithmetic_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_selective_distance12_dimension136_20261017',
 'fresh_selective_distance12_dimension136_20261017_review',
 'fresh_nonincident_carrier_incidence_20261017',
 'fresh_actual_corner_enlarged_model_20261017',
 'fresh_actual_corner_enlarged_model_20261017_support',
 'fresh_actual_corner_enlarged_model_20261017_local_cones',
 'fresh_actual_corner_enlarged_model_20261017_joint_ports',
 'fresh_new_cavity_construction_20261017',
 'fresh_new_cavity_construction_gate_20261017',
)
EXTRA_REPORTS=()

@lru_cache(None)
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())

def small_checks():
    p=read(STEMS[0]);assert p['dimension_bound']==136 and p['inherited_dimension_bound']==146
    review=read(STEMS[1]);assert review['status'].startswith('PASS')
    assert p['unrestricted_finite_volume'] and not p['arithmetic_hypothesis'] and p['cusp_restriction'] is None
    assert p['excluded_dimensions']==list(range(137,147))
    rows=p['reduced_profiles'];assert len(rows)==1175
    for r in rows:
        rs=r['ranks'];d=sum(rs);v=prod(t+1 for t in rs);h=[1]
        for t in rs:
            q=[0]*(len(h)+t)
            for i,x in enumerate(h):
                for j in range(t+1):q[i+j]+=x
            h=q
        G=sum(x*min(i,d-i)for i,x in enumerate(h))
        assert d==r['link_dimension'] and v==r['link_vertices'] and G==r['G']
        curv=F(5*G-d*v,2*v);cost=sum((min(F(12),F(t,2))for t in rs),F())
        assert curv==F(r['normalized_curvature']) and cost==F(r['normalized_cost'])
        assert curv-cost==F(r['margin'])>0
    assert min(F(r['margin'])for r in rows)==F(172096,1953125)
    matched={tuple(sorted(r['ranks'])):F(r['margin'])for r in review['budgets']['profiles']}
    assert len(matched)==1175 and matched=={tuple(sorted(r['ranks'])):F(r['margin'])for r in rows}
    assert F(p['boundary_at_n136']['margin'])==F(-59157,3906250)
    assert p['compact_facet_ambient_bound']==103
    assert p['compact_facet_parity_slack_times8']==[[12,208,4],[112,212,4]]
    c=p['census'];assert len(c['ports'])==75 and sum(r['allowed']for r in c['ports'])==18
    assert len(c['disjoint_joins'])==115 and c['overlap_generated']==208
    assert len(c['overlap_records'])==162 and c['overlap_counts']=={'inverse_rejection':120,'white_affine_rejection':6,'shared_vertex':27,'shared_edge':9}
    assert review['classification']['overlap_counts']==c['overlap_counts']
    def portmap(ports):return {(r['diagram'],tuple(r['black']),r['node']):r for r in ports}
    assert portmap(c['ports'])==portmap(review['classification']['all_distance_at_least5_ports'])
    assert Counter(r['type']for r in p['isolation']['cases'])=={'shared_vertex':27,'shared_edge':9}
    corners={(r['type'],tuple(r['left_black']),tuple(r['right_black'])):r for r in review['isolation']['preferred_corners']}
    assert len(corners)==36
    for r in p['isolation']['cases']:
        s=corners[(r['type'],tuple(r['left_black']),tuple(r['right_black']))]
        assert r['possible_contacts']==s['allowed_core_contacts']
        for key in ('preferred_corner','marked_component','marked_distance','shortcut_bound'):assert r[key]==s[key]
    for r in p['disjoint_obstructions']['simple_leaf_joins']:
        a=[[F(x)for x in row]for row in r['schur']]
        det=a[0][0]*a[1][1]-a[0][1]*a[1][0]
        assert a[0][0]<0 and det==F(r['determinant'])>0
    model=read(STEMS[3]);primal,dual=model['certificates'];assert model['solver_invocations']==1
    lower=F(dual['bound']);upper=F(primal['maximum'])
    assert lower==F(338427277204734900,27027022505771209)>12
    assert upper==F(3130455146279,250000000000)>=lower
    alpha=dual['face_multipliers'];mu={i:F(x)for i,x in dual['vertex_multipliers']}
    assert sum(mu.values(),F())==F(dual['normalizer'])
    assert F(dual['objective_numerator'])/sum(mu.values(),F())==lower
    ratios=[]
    for stem in STEMS[5:7]:
        d=read(stem);assert d['face_multipliers']==alpha
        initial={i:F(x)for i,x in d['initial_vertex_multipliers']};assert initial==mu
        mu={i:F(x)for i,x in d['vertex_multipliers']}
        assert all(mu.get(i,F())>=x for i,x in initial.items())
        norm=sum(mu.values(),F());assert norm==F(d['normalizer'])
        assert F(d['objective_numerator'])==F(dual['objective_numerator'])
        ratio=F(d['objective_numerator'])/norm;assert ratio==F(d['bound'])>12
        ratios.append(str(ratio))
    support=read(STEMS[4]);assert support['counts']['selected_actual_Q']==397
    assert support['counts']['unselected_merged_Q']==1897 and support['counts']['simplex_faces']==2638
    restricted=support['restricted_actual_Q_test']['repaired_dual'];assert F(restricted['bound'])<12
    local=read(STEMS[5]);joint=read(STEMS[6])
    assert local['complete_positive_simplex_forest_coverage']
    assert local['counts']['forests']==815 and local['counts']['block_orbits']==5011
    assert joint['counts']['total_inequalities']==3519737
    assert joint['counts']['unselected_Q_rows_with_zero_new_coefficients']==1897
    carrier=read(STEMS[2])['face_ledger'];assert carrier['counts']=={'quadrilaterals':1296,'triangles':2160}
    assert carrier['quad_weight_counts']=={'0':1048,'1':243,'2':5}
    assert carrier['compact_quad_weight_counts']=={'0':987,'1':240,'2':5}
    assert 1048>648 and 987>6*2*81
    scout=read(STEMS[7]);assert scout['R_bank_size']==94
    assert scout['counts']=={'low_norm_equations':242,'real_low_branches':484,'complete_high_label_shape_families':12,'continuous_high_samples':780}
    gate=read(STEMS[8]);assert gate['status']=='PASS_EXACT_ALL_LABEL_FOUR_CENTER_CHANGED_PENTAGON_ENDPOINT_GATE'
    assert gate['counts']=={'tested_boxes':2305,'certified_boxes':1295,'maximum_depth':21,'R_low_recipes':96,'R_high_families':10,'AR_cases':285}
    assert len(gate['records'])==106
    assert sum(len(r['complete_AR_cases'])for r in gate['records'])==285
    treecounts=Counter()
    def walk(node,depth=0):
        treecounts['nodes']+=1;treecounts['maximum_depth']=max(treecounts['maximum_depth'],depth)
        if 'children'in node:
            assert len(node['children'])==2
            for child in node['children']:walk(child,depth+1)
        else:assert 'box'in node and 'certificate'in node;treecounts['leaves']+=1
    for r in gate['records']:
        for case in r['complete_AR_cases']:walk(case['proof'])
    assert treecounts=={'nodes':2305,'leaves':1295,'maximum_depth':21}
    return dict(scalar_profile_checks=1175,independent_profile_matches=1175,independent_port_matches=75,
                independent_corner_matches=36,unrestricted_dimension_bound=136,compact_facet_ambient_bound=103,
                overlap_cases=162,isolated_marked_contexts=36,finite_model_exact_lower=str(lower),
                finite_model_exact_upper=str(upper),full_local_and_partial_joint_repaired_ratios=ratios,
                selected_actual_Q=397,unexpanded_merged_Q=1897,compact_unpaid_quad_counterexample=987,
                construction_all_label_cases=285,construction_interval_leaves=1295,
                scope='Exact new scalar composition and finite-model mass monotonicity; independent proof reviews and scoped construction certificates are source-bound. No repeated large LP or historical census.')

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
    checks=small_checks();print('ARCHIVE_FRACTION_CHECKS_PASS',checks['scalar_profile_checks'],flush=True)
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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_DISTANCE12_AND_ENLARGED_ARITHMETIC',
                scope='Archive integrity and exact composition of distance12 dimension136, enlarged finite arithmetic models and bounded construction gates; not formalization or repeated historical verification.',
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
