#!/usr/bin/env python3
"""Incremental archive of distance11, raw25 corner capacities and changed cap adjacency.

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
CHECKPOINT='results/distance11_and_raw25_corner_frontier_20261017.md'
PREVIOUS='results/distance12_and_enlarged_arithmetic_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_distance12_and_enlarged_arithmetic_frontier_20261017.py'
OUT=ROOT/'results/distance11_and_raw25_corner_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_selective_distance11_dimension126_20261017',
 'fresh_selective_distance11_dimension126_20261017_review',
 'fresh_distance11_extension_review_20261017',
 'fresh_distance11_extension_review_20261017_mixed',
 'fresh_raw25_actual_corner_compatibility_20261017',
 'fresh_broken_cap_adjacency_20261017',
 'fresh_broken_cap_adjacency_events_20261017',
 'fresh_broken_cap_adjacency_split_patch_20261017',
 'fresh_broken_cap_adjacency_gate_20261017',
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
    p=read(STEMS[0]);review=read(STEMS[1])
    assert p['dimension_bound']==126 and p['inherited_dimension_bound']==136
    assert p['unrestricted_finite_volume'] and not p['arithmetic_hypothesis'] and p['cusp_restriction'] is None
    assert p['excluded_dimensions']==list(range(127,137))
    assert review['status'].startswith('PASS')
    rows=p['reduced_profiles'];assert len(rows)==1095
    for r in rows:
        rs=r['ranks'];d=sum(rs);v=prod(t+1 for t in rs);h=[1]
        for t in rs:
            q=[0]*(len(h)+t)
            for i,x in enumerate(h):
                for j in range(t+1):q[i+j]+=x
            h=q
        G=sum(x*min(i,d-i)for i,x in enumerate(h))
        curvature=F(5*G-d*v,2*v);cost=sum((min(F(11),F(t,2))for t in rs),F())
        assert d==r['link_dimension'] and v==r['link_vertices'] and G==r['G']
        assert curvature==F(r['normalized_curvature']) and cost==F(r['normalized_cost'])
        assert curvature-cost==F(r['margin'])>0
    assert min(F(r['margin'])for r in rows)==F(40603,559682)
    assert F(p['boundary_at_n126']['margin'])==F(-238535,12872686)
    assert p['compact_facet_ambient_bound']==95
    assert p['compact_facet_parity_slack_times8']==[[12,192,4],[104,196,4]]
    c=p['census'];j=p['disjoint_obstructions'];iso=p['isolation']
    assert c['marked_contexts']==384 and c['allowed_ports']==55
    assert c['overlap_generated']==2137 and c['consistent_overlaps']==1010
    expected={'inverse_rejection':856,'edge_pd_rejection':40,'mixed_22_25':9,
              'shared_vertex':60,'shared_edge':27,'mixed_25_22':9,'shared_three_path':9}
    assert c['overlap_counts']==expected
    assert j['candidate_joins']==425 and j['white_rejections']==119
    assert j['remaining_unique_context_pairs']==306 and len(j['underlying_inertia_types'])==13
    assert all(r['signature'][1]>=2 for r in j['underlying_inertia_types'])
    assert len(iso['cases'])==114 and iso['component_order_cap']==17
    matched={tuple(sorted(r['ranks'])):F(r['margin'])for r in review['budgets']['profiles']}
    assert len(matched)==1095 and matched=={tuple(sorted(r['ranks'])):F(r['margin'])for r in rows}
    def portmap(ports):return {(r['diagram'],tuple(r['black']),r['node']):r for r in ports}
    assert len(c['ports'])==200 and portmap(c['ports'])==portmap(review['classification']['ports'])
    def overlap_key(r):return (r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['match'])))
    primary_overlap={overlap_key(r):r for r in c['overlap_records']}
    review_overlap={overlap_key(r):r for r in review['classification']['overlap_records']}
    assert len(primary_overlap)==len(review_overlap)==1010 and primary_overlap.keys()==review_overlap.keys()
    kinds={'edge_pd_rejection':'edge_rejection','shared_three_path':'shared_path',
           'mixed_22_25':'mixed','mixed_25_22':'mixed'}
    for key,r in primary_overlap.items():
        kind='inverse_rejection'if'inverse_rejection'in r else'edge_pd_rejection'if'edge_pd_rejection'in r else r['surviving_type']
        assert kinds.get(kind,kind)==review_overlap[key]['kind']
    pcorners={(r['type'],tuple(r['left_black']),tuple(r['right_black'])):r for r in iso['cases']}
    assert len(pcorners)==114 and len(review['isolation']['preferred_corners'])==114
    short=extra=0
    for r in review['isolation']['preferred_corners']:
        typ=r['type'];di=r['left'][0];dj=r['right'][0]
        pt='shared_three_path'if typ=='shared_path'else f'mixed_{di}_{dj}'if typ=='mixed'else typ
        s=pcorners[pt,tuple(r['left_black']),tuple(r['right_black'])]
        cross=s['extra_cross_bond_checks']
        assert cross['shortcut_rejections']==r['extra_cross_shortcut_rejections']
        assert len(cross['edge_pd_z_rejections'])==len(r['extra_cross_pd_rejections'])
        short+=cross['shortcut_rejections'];extra+=len(cross['edge_pd_z_rejections'])
        if typ in('shared_vertex','shared_edge'):
            for key in('preferred_corner','blocks','possible_contacts','marked_component','marked_distance'):
                assert r[key]==s[key]
        else:
            mapping=s['literal_isomorphism'];local=s['local_lemma_context']
            preferred=local['preferred']if typ=='mixed'else local['preferred']['preferred_corner']
            component=local['marked_component']if typ=='mixed'else local['preferred']['marked_component']
            available=local['allowed_contacts']if typ=='mixed'else local['preferred']['possible_contacts']
            assert set(mapping[i]for i in r['preferred_corner'])==set(preferred)
            assert set(mapping[i]for i in r['marked_component'])==set(component)
            assert set(mapping[i]for i in r['possible_contacts'])==set(available)
    assert short==8280 and extra==45
    assert review['saved_comparison']['independently_verified_saved_nonnegative_witnesses']==204
    triple=read(STEMS[2]);mixed=read(STEMS[3])
    assert triple['marked_component_order_bound']==mixed['marked_component_order_bound']==16
    assert len(triple['marked_contexts'])==len(mixed['contexts'])==9
    assert all(max(e['order']for e in r['one_extension_variants'])<=16 for r in triple['marked_contexts'])
    assert all(len(r['marked_component'])<=16 for r in mixed['contexts'])
    raw=read(STEMS[4]);counts=raw['counts'];repair=raw['fixed_witness_repair']
    baseline=read(repair['baseline'])
    assert raw['face_multipliers']==baseline['face_multipliers']
    initial={i:F(x)for i,x in raw['initial_vertex_multipliers']}
    mu={i:F(x)for i,x in raw['vertex_multipliers']}
    assert initial=={i:F(x)for i,x in baseline['vertex_multipliers']}
    assert all(mu.get(i,F())>=x for i,x in initial.items())
    alpha={i:F(x)for i,x in raw['face_multipliers']}
    violations=0
    for col in raw['columns']:
        load=sum((alpha[i]*F(k)for i,k in col['positive_face_coefficients']),F())
        vertex=col['vertex'];cost=F(col['vertex_cost'])
        assert load==F(col['load']) and F(col['baseline_vertex_mass'])==initial.get(vertex,F())
        assert cost>0 and load<=cost*mu.get(vertex,F())
        assert col['violated']==(load>cost*initial.get(vertex,F()))
        violations+=col['violated']
    norm=sum(mu.values(),F());numerator=F(repair['objective_numerator']);ratio=numerator/norm
    assert norm==F(repair['normalizer']) and numerator==F(baseline['objective_numerator'])
    assert ratio==F(repair['bound'])==F(67685455440946980000000000000000,5466425626097138727892456269673)>12
    assert F(repair['baseline_bound'])==F(baseline['bound'])>ratio
    assert repair['optimizer_invocations']==0 and repair['face_masses_unchanged']
    assert counts['actual_Q_classes']==55 and counts['new_unique_context_classes']==36
    assert counts['corner_occurrences']==440 and counts['literal_orientations']==6260
    assert counts['literal_pairs_tested']==149510 and counts['exact_literal_exclusions']==141850
    assert len(raw['columns'])==counts['columns']==500 and violations==counts['violated_columns']==28
    assert sum(mu.get(i,F())>initial.get(i,F())for i in set(initial)|set(mu))==25
    scout=read(STEMS[5]);events=read(STEMS[6]);split=read(STEMS[7]);gate=read(STEMS[8])
    assert scout['counts']['norm_equations']==50072 and not scout['accepted_bridges']
    assert events['counts']['R_event_norm_equations']==4476 and not events['coherent_ARX_triples']
    assert split['counts']['norm_equations']==5665 and not split['whole_patch_candidates']
    assert not any(d['unresolved_events']for d in(scout,events,split))
    assert gate['counts']=={'tested_boxes':1038,'certified_boxes':691,'maximum_depth':11,
        'X_recipe_shapes':148,'AX_cases':344,'maximum_high_parameters':2}
    assert len(gate['records'])==148 and sum(len(r['A_cases'])for r in gate['records'])==344
    treecounts=Counter()
    def walk(node,depth=0):
        treecounts['nodes']+=1;treecounts['maximum_depth']=max(treecounts['maximum_depth'],depth)
        if 'children'in node:
            assert len(node['children'])==2
            for child in node['children']:walk(child,depth+1)
        else:assert 'box'in node and 'certificate'in node;treecounts['leaves']+=1
    for r in gate['records']:
        for case in r['A_cases']:walk(case['proof'])
    assert treecounts=={'nodes':1038,'leaves':691,'maximum_depth':11}
    return dict(scalar_profile_checks=1095,independent_profile_matches=1095,
        independent_port_matches=200,independent_overlap_matches=1010,independent_corner_matches=114,
        extra_cross_shortcuts=8280,extra_cross_pd_obstructions=45,
        unrestricted_dimension_bound=126,compact_facet_ambient_bound=95,
        marked_contexts=384,allowed_ports=55,disjoint_joins=425,overlap_cases=1010,
        isolated_or_controlled_marked_contexts=114,raw25_actual_classes=55,raw25_columns=500,
        raw25_monotone_repair=str(ratio),construction_all_label_cases=344,construction_interval_leaves=691,
        scope='Exact scalar composition and finite-model mass monotonicity; source-bound independent mathematical reviews and scoped construction certificates. No repeated optimizer or historical census.')

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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_DISTANCE11_AND_RAW25_CORNERS',
                scope='Archive integrity and exact composition of distance11 dimension126, fixed-field raw25 corner capacities and changed-incidence construction searches/gates; not formalization or repeated historical verification.',
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
