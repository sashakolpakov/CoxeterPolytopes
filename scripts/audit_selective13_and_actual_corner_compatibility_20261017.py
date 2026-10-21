#!/usr/bin/env python3
"""Incremental archive of selective13 and full finite-vertex corner geometry.

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
CHECKPOINT='results/selective13_and_actual_corner_compatibility_20261017.md'
PREVIOUS='results/distance14_and_small_cusp_factor_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_distance14_and_small_cusp_factor_frontier_20261017.py'
OUT=ROOT/'results/selective13_and_actual_corner_compatibility_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_selective_distance13_dimension146_20261017',
 'fresh_selective_distance13_dimension146_20261017_review',
 'fresh_small_cusp_defects_20261017',
 'fresh_dominant_cap_compatibility_20261017',
 'fresh_dominant_cap_compatibility_family_20261017',
 'fresh_dominant_cap_compatibility_family_20261017_review',
 'fresh_three_cap_pentagon_20261017',
 'fresh_three_cap_pentagon_changed_seed_20261017',
 'fresh_three_cap_pentagon_changed_seed_gate_20261017',
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
    p=read(STEMS[0]);q=read(STEMS[1])
    assert p['dimension_bound']==146 and p['inherited_dimension_bound']==156
    assert p['unrestricted_finite_volume'] and not p['arithmetic_hypothesis'] and p['cusp_restriction']is None
    assert p['excluded_dimensions']==list(range(147,157))
    primary={tuple(sorted(r['ranks'])):r for r in p['reduced_profiles']}
    review={tuple(sorted(r['ranks'])):r for r in q['budgets']['profiles']}
    assert len(primary)==len(review)==1255 and primary.keys()==review.keys()
    for ranks,r in primary.items():
        d=sum(ranks);v=prod(t+1 for t in ranks)
        assert d==r['link_dimension'] and v==r['link_vertices']
        curv=F(5*r['G']-d*v,2*v)
        cost=sum((min(F(13),F(t,2))for t in ranks),F())
        assert curv==F(r['normalized_curvature']) and cost==F(r['normalized_cost'])
        assert curv-cost==F(r['margin'])==F(review[ranks]['margin'])>0
    assert min(F(r['margin'])for r in primary.values())==F(830407,9565938)
    assert F(p['boundary_at_n146']['margin'])==F(-515185,19131876)<0
    assert p['compact_facet_ambient_bound']==111
    assert p['compact_facet_parity_slack_times8']==[list(map(int,r))for r in q['budgets']['compact_facet_parity_slack_times8']]==[[12,224,4],[120,228,4]]
    ports=[]
    for r in p['radius_contexts']:
        for t in r['ports']:ports.append(dict(diagram=r['diagram'],black=r['black'],**t))
    assert ports==q['ports']['high_distance_ports'] and len(ports)==21
    assert q['ports']['marked_pairs']==384 and q['ports']['high_distance_contexts']==16
    for a,b in zip(p['supporting_white_isolation']['cases'],q['core']['preferred_corners'],strict=True):
        for key in ('left_black','right_black','preferred_corner','possible_core_contacts','type'):assert a[key]==b[key]
    assert p['known_sharp_distance14_pyramid']['vertices']=={'finite':80,'ideal':2}
    assert q['known_pyramid']['finite']==80 and q['known_pyramid']['ideal']==2
    corner=read('fresh_dominant_cap_compatibility_20261017')
    old=read('fresh_sqrt5_scaled_core_cap_edge_repair_20261017')
    alpha={i:F(a)for i,a in old['face_multipliers']};mu={i:F(a)for i,a in old['vertex_multipliers']}
    c=corner['inequality'];assert c['positive_face_coefficients']==[[26448,1],[26474,2]]
    load=alpha[26448]+2*alpha[26474];assert load==F(c['load'])==F(310256599,10**10)
    gap=load-mu[3137];assert gap==F(c['gap'])==F(55089993,5*10**9)>0
    fix=corner['fixed_witness_repair'];normalizer=F(old['normalizer'])+gap
    assert F(fix['normalizer'])==normalizer
    ratio=F(old['objective_numerator'])/normalizer;assert ratio==F(fix['bound'])<12
    assert fix['vertex_updates']==[[3137,str(mu[3137]),str(load)]] and fix['face_masses_unchanged']
    reasons=Counter(r['reason']for r in corner['incompatibilities'])
    assert reasons=={'distinct_full_stars_share_actual_quad_port_pair':7,'forced_distinct_wall_product_in_open_minus1_zero_gap':8}
    family=read('fresh_dominant_cap_compatibility_family_20261017')
    mu[3137]=load;baseline_mu=dict(mu);violated=[]
    assert sum(mu.values(),F())==normalizer
    for i,c in enumerate(family['columns']):
        assert c['index']==i and c['vertex_cost']==1
        v=c['vertex'];a=sum((m*alpha.get(f,F())for f,m in c['positive_face_coefficients']),F())
        assert all(isinstance(m,int)and m>0 for f,m in c['positive_face_coefficients'])
        assert a==F(c['load']) and F(c['baseline_vertex_mass'])==baseline_mu.get(v,F())
        assert c['violated']==(a>baseline_mu.get(v,F()))
        if c['violated']:violated.append(i)
        mu[v]=max(mu.get(v,F()),a)
    updates=[[v,str(baseline_mu.get(v,F())),str(mu[v])]for v in sorted(mu)if mu[v]!=baseline_mu.get(v,F())]
    repair=family['fixed_witness_repair']
    assert updates==repair['vertex_updates'] and len(updates)==len(violated)==120
    assert len(family['columns'])==294 and family['counts']['corner_occurrences']==570
    assert family['counts']['actual_Q_classes']==57 and family['counts']['forests']==299
    assert family['counts']['literal_orientations']==6010 and family['counts']['capacity_one_orbits']==560
    assert family['counts']['literal_pairs_tested']==134445
    assert family['counts']['exact_literal_exclusions']==134435
    assert family['counts']['distinct_full_stars_share_actual_quad_port_pair']==22397
    assert family['counts']['actual_wall_product_in_open_minus1_zero']==94464
    assert family['counts']['positive_actual_wall_product']==17574
    assert 22397+94464+17574==134435
    family_normalizer=sum(mu.values(),F())
    assert family_normalizer==F(repair['normalizer'])
    family_ratio=F(old['objective_numerator'])/family_normalizer
    assert family_ratio==F(repair['bound'])<ratio<12
    assert repair['baseline']=='fresh_dominant_cap_compatibility_20261017'
    assert F(repair['baseline_bound'])==ratio and repair['face_masses_unchanged']
    assert repair['optimizer_invocations']==0
    gate=read('fresh_three_cap_pentagon_20261017')
    assert gate['counts']=={'mandatory_frames':3,'low_recipes':9,'low_norm_branches':18,'complete_high_label_families':2,'whole_interval_branches':4,'rotations':5}
    assert gate['deleted_walls']==[12,47,52]
    scout=read('fresh_three_cap_pentagon_changed_seed_20261017')
    changed=read('fresh_three_cap_pentagon_changed_seed_gate_20261017')
    assert scout['deleted_walls']==changed['deleted_walls']==[12,47,52,64,65]
    assert scout['boundary_cycle']==changed['boundary_cycle']==[7,8,44,63,30]
    assert scout['exact_old_deletion_cost']==changed['exact_old_deletion_cost']
    ledger=scout['exact_old_deletion_cost']
    assert ledger['lost_finite_count']==len(ledger['lost_finite_supports'])==11
    assert ledger['compact_faces_by_dimension']=={'2':6}
    assert len(scout['accepted_low_caps'])==94
    assert all(c['frame']==[8,44,63]for c in scout['accepted_low_caps'])
    assert len(changed['all_five_four_contact_gates'])==5
    assert sum(len(c['low'])for c in changed['all_five_four_contact_gates'])==282
    assert sum(len(c['tails'])for c in changed['all_five_four_contact_gates'])==41
    assert len(changed['mandatory_three_contact_gates'])==3
    small=read('fresh_small_cusp_defects_20261017')['cusp_balance']
    assert len(small['general_profiles'])==128 and len(small['marked_small_profiles'])==988
    return dict(independent_profile_matches=1255,independent_port_matches=21,
                unrestricted_dimension_bound=146,compact_facet_ambient_bound=111,
                fixed_witness_new_cap_load=str(load),fixed_witness_allowance_raise=str(gap),
                fixed_witness_repaired_bound=str(ratio),partial_feasibility_or_new_arithmetic_theorem=False,
                actual_corner_family_columns=294,actual_corner_family_allowance_raises=120,
                actual_corner_family_repaired_bound=str(family_ratio),
                three_cap_low_branches=18,three_cap_high_interval_branches=4,
                changed_seed_four_contact_low_recipes=282,changed_seed_complete_high_families=41,
                changed_seed_low_caps_numerical_only=94,
                scope='Exact new scalar/mass composition and independent-certificate matching; new geometry is supplied by separately written reviews. Family and changed-seed packages are source-bound and retain their own completeness scopes.')

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
    checks=small_checks();print('ARCHIVE_FRACTION_CHECKS_PASS',checks['independent_profile_matches'],flush=True)
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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_SELECTIVE13_AND_ACTUAL_CORNER_COMPATIBILITY',
                scope='Archive integrity and exact composition of selective13 dimension146, actual finite-corner compatibility and bounded construction gates; not formalization or repeated historical verification.',
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
