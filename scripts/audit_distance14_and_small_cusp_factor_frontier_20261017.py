#!/usr/bin/env python3
"""Incremental frozen archive: distance14, small cusp factors, and edge-cap repair.

Old dependency graphs are inherited only after every frozen byte digest is
rechecked. New JSON manifests and record references are inspected recursively.
No producer verifier, optimizer, or historical million-case scan is rerun.
"""
import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, prod
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/distance14_and_small_cusp_factor_frontier_20261017.md'
PREVIOUS='results/forced_quadrilaterals_and_affine_carriers_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_forced_quadrilaterals_and_affine_carriers_20261017.py'
OUT=ROOT/'results/distance14_and_small_cusp_factor_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
    'fresh_cusp_rank7_fine_curvature_20261017',
    'fresh_affine_carriers_order10_20261017',
    'fresh_affine_order10_quad_payment_20261017',
    'fresh_cusp_minimum_factor_rank9_dimension29_20261017',
    'fresh_distance14_unrestricted_dimension156_20261017',
    'fresh_distance14_unrestricted_dimension156_20261017_review',
    'fresh_sqrt5_scaled_core_cap_incompatibilities_20261017',
    'fresh_sqrt5_scaled_core_cap_edge_repair_20261017',
    'fresh_nonarithmetic_release12_47_52_pentagon_20261017',
    'fresh_nonarithmetic_pentagon_two_center_gate_20261017',
)
EXTRA_REPORTS=(
    'results/fresh_vinberg_finite_quad_affine_cutoff9_20261017.md',
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
    base=read('fresh_sqrt5_core23_25_spectator_incidence_20261017')
    cover=read('fresh_sqrt5_a3_d5_cap_local_cover_update_20261017')
    caps=read('fresh_sqrt5_scaled_core_cap_incompatibilities_20261017')
    new=read('fresh_sqrt5_scaled_core_cap_edge_repair_20261017')
    alpha0=masses(base['face_multipliers'],82806);alpha=masses(new['face_multipliers'],82806)
    mu0=masses(base['vertex_multipliers'],6558);mu=masses(new['vertex_multipliers'],6558)
    def changed(a,b):return{i for i in set(a)|set(b)if a.get(i,F())!=b.get(i,F())}
    da={i:alpha[i]-alpha0[i]for i in changed(alpha0,alpha)}
    assert da=={80010:F(9,2500000000),82048:F(2991,200000000)}
    assert set(alpha)==set(alpha0) and all(mu.get(i,F())>=v for i,v in mu0.items())
    assert changed(mu0,mu)=={1614,1656,1875}
    assert new['simplex_mass_raises']==[[i,str(alpha0[i]),str(alpha[i])]for i in sorted(da)]
    assert new['simplex_mass_raises']==caps['required_simplex_masses']
    assert new['vertex_mass_raises']==[[v,str(mu0.get(v,F())),str(mu[v])]for v in sorted(changed(mu0,mu))]
    for key in('positive_simplex_features','block_orbit_columns','literal_block_orbits'):
        assert new[key]==cover[key]
    features=new['positive_simplex_features'];cols=new['block_orbit_columns']
    assert len(features)==2188 and len(cols)==4169
    byface={f['face']:f for f in features};assert len(byface)==len(features)
    assert set(byface)<=set(alpha)
    # A singleton's local intersection with any subset is0 or1. This proves
    # the unchanged joint loads; no historical joint enumeration is repeated.
    for i in da:
        f=byface[i];assert f['simplex_dimension']==1
        assert sum(m for j,m in f['block_orbits'])==2
        assert all(cols[j]['block_size']==1 for j,m in f['block_orbits'])
    assert new['joint_singleton_zero_checks']==[[k,comb(13,k)]for k in(3,4,5)]
    demand0=defaultdict(F);demand=defaultdict(F)
    for f in features:
        for j,m in f['block_orbits']:
            assert cols[j]['index']==j and isinstance(m,int)and m>0
            demand0[j]+=m*alpha0[f['face']];demand[j]+=m*alpha[f['face']]
    assert {j for j in demand if demand[j]!=demand0[j]}=={2643,3137,4043,4044}
    oldcovers={r['vertex']:r for r in cover['local_covers']}
    assert len(oldcovers)==len(new['local_covers'])==648
    expected_mu=dict(mu0);added=[];oldterms=0
    for r in new['local_covers']:
        v=r['vertex'];p=oldcovers[v]
        assert r['forest']==p['forest'] and r['literal_block_orbits']==p['literal_block_orbits']
        old=p['inherited_remapped_packing_mixture']+p['singleton_packing_increments']
        assert r['inherited_remapped_packing_mixture']==old
        oldterms+=len(old)
        # Inherited literal disjointness is already certified and byte-bound.
        # Recheck the four new terms and exact marginal arithmetic here.
        covered={j:F(a)for j,a in p['covered_marginals']}
        if not covered:
            covered={j:F()for j,a in p['required_marginals']}
        oldmass=sum((F(t['mass'])for t in old),F())
        assert oldmass==F(p['new_cover_mass'])==F(r['old_cover_mass'])
        assert {j:F(a)for j,a in p['required_marginals']}=={j:demand0[j]for j in covered}
        assert {j:F(a)for j,a in r['required_marginals']}=={j:demand[j]for j in covered}
        incmass=F()
        for t in r['singleton_packing_increments']:
            assert len(t['blocks'])==1
            j,mask=t['blocks'][0];a=F(t['mass'])
            assert a>0 and a==max(F(),demand[j]-covered[j])
            assert cols[j]['vertex']==v and cols[j]['block_size']==1
            assert mask.bit_count()==1 and mask in new['literal_block_orbits'][j]['literal_masks']
            covered[j]+=a;incmass+=a;added.append([v,j,mask,str(a)])
        assert all(covered[j]>=demand[j]for j in covered)
        assert {j:F(a)for j,a in r['covered_marginals']}==covered
        assert F(r['new_cover_mass'])==oldmass+incmass
        expected_mu[v]=max(expected_mu.get(v,F()),oldmass+incmass)
    assert oldterms==3712 and len(added)==4 and expected_mu==mu
    # The15 new columns use the common actual provenance map. Their geometry
    # is inherited from the separately reviewed producer, not rederived here.
    assert caps['counts']['columns']==15 and len(caps['columns'])==15
    violations=[]
    for ix,c in enumerate(caps['columns']):
        assert c['index']==ix and c['capacity']==1
        load0=sum((F(m)*alpha0.get(i,F())for i,m in c['positive_face_coefficients']),F())
        load=sum((F(m)*alpha.get(i,F())for i,m in c['positive_face_coefficients']),F())
        assert load==load0 # only negative edge-side masses changed
        if load0>alpha0.get(c['simplex_face'],F()):violations.append(ix)
        assert load<=alpha.get(c['simplex_face'],F())
    assert violations==caps['violating_columns']==[5,13]
    incidence=[]
    for c in base['constraints']:
        l=sum((F(m)*alpha.get(i,F())for i,m in c['positive_faces']),F())
        l-=sum((F(m)*alpha.get(i,F())for i,m in c['negative_faces']),F())
        r=sum((F(m)*mu.get(v,F())for v,m in c['vertex_capacity']),F())
        assert l<=r;incidence.append([c['index'],str(l),str(r)])
    assert len(incidence)==52 and incidence==new['forced_incidence_checks']
    objective=F(base['objective_numerator']);assert objective==F(new['objective_numerator'])==F(3842394855727,309375000000)
    assert F(new['normalizer'])==sum(mu.values(),F())
    bound=objective/F(new['normalizer']);assert bound==F(new['bound'])>12
    assert bound<F(base['bound'])==F(new['initial_bound'])
    assert new['optimizer_invocations']==new['repeated_large_joint_scans']==0
    assert new['counts']['cap_checks']==2071 and new['counts']['inherited_signed_columns']==94722
    assert new['counts']['joint_families_preserved']==17
    # Small exact scalar/status checks, without importing any producer.
    fine=read('fresh_cusp_rank7_fine_curvature_20261017')
    assert fine['minimum_factor_rank']==7 and fine['ambient_dimension_at_least']==30
    assert F(fine['curvature_increment'])==F(1,18) and F(fine['minimum_base_margin'])>0
    assert F(fine['rank6_counterexample']['margin_over_3s'])<0
    theorem=read('fresh_cusp_minimum_factor_rank9_dimension29_20261017')
    assert theorem['status']=='PASS_ALL_RANK_MINIMUM_CUSP_FACTOR_RANK9_DIMENSION29_COMPOSITION'
    assert theorem['dimension_bound']==29 and theorem['minimum_cusp_factor_rank']==9
    assert theorem['arithmetic_hypothesis']is False and theorem['factor_count_limit']is None
    assert theorem['enumeration_cutoff_used']is False and theorem['unrestricted163_tail_used']is False
    assert list(map(F,theorem['even_finite_slack_times8_coefficients']))==[F(0),F(60),F(4)]
    assert list(map(F,theorem['odd_finite_slack_times8_coefficients']))==[F(36),F(64),F(4)]
    universal=read('fresh_distance14_unrestricted_dimension156_20261017')
    assert universal['status']=='PASS_DISTANCE14_AND_CONCENTRATED_CUSP_BUDGET_DIMENSION156'
    assert universal['dimension_bound']==156 and universal['unrestricted_finite_volume']is True
    assert universal['arithmetic_hypothesis']is False and universal['simple_cusp_hypothesis']is False
    assert universal['factor_rank_or_count_restriction']is None
    assert universal['inherited_unrestricted_dimension_bound']==163
    assert universal['excluded_ambient_dimensions']==list(range(157,164))
    assert universal['triangle_distance_bound']==8 and universal['quadrilateral_distance_bound']==14
    ps=universal['reduced_profiles'];assert len(ps)==924
    byn=defaultdict(list)
    for p in ps:
        d=sum(p['ranks']);V=prod(r+1 for r in p['ranks'])
        assert d==p['link_dimension'] and V==p['link_vertices']
        curv=F(5*p['G'],2*V)-F(d,2)
        cost=sum((min(F(14),F(r,2))for r in p['ranks']),F())
        assert curv==F(p['normalized_curvature']) and cost==F(p['normalized_cost_bound'])
        gap=curv-cost;assert gap==F(p['margin'])>0;byn[d+1].append(p)
    assert sorted(byn)==list(range(157,164))
    for m in universal['dimension_minima']:
        ps= byn[m['ambient_dimension']]
        assert len(ps)==m['profiles'] and min(F(p['margin'])for p in ps)==F(m['minimum_margin'])
    assert min(F(p['margin'])for p in universal['reduced_profiles'])==F(universal['minimum_profile_margin'])
    assert F(universal['distance14_budget_obstruction_at_n156']['margin'])<0
    independent=read('fresh_distance14_unrestricted_dimension156_20261017_review')
    assert independent['status']=='PASS_INDEPENDENT_DISTANCE14_UNRESTRICTED_DIMENSION156_REVIEW'
    assert independent['dimension_bound']==156 and independent['unrestricted_finite_volume']is True
    assert independent['excluded_dimensions']==universal['excluded_ambient_dimensions']
    assert independent['inherited_upper_tail_bound']==163
    # The independent review enumerates rank partitions in a different order
    # and does not save the producer's intermediate G. Compare intrinsic data.
    primary_profiles={tuple(sorted(p['ranks'])):p for p in universal['reduced_profiles']}
    review_profiles={tuple(sorted(p['ranks'])):p for p in independent['reduced_profiles']}
    assert len(primary_profiles)==len(review_profiles)==924 and primary_profiles.keys()==review_profiles.keys()
    for ranks,p in primary_profiles.items():
        q=review_profiles[ranks]
        for key in('link_dimension','link_vertices','normalized_curvature','normalized_cost_bound','margin','mean_absolute_deviation'):
            assert p[key]==q[key],(ranks,key)
    for p,q in zip(universal['dimension_minima'],independent['dimension_minima'],strict=True):
        for key in('ambient_dimension','profiles','minimum_margin'):assert p[key]==q[key]
        assert sorted(tuple(sorted(r))for r in p['minimizing_ranks'])==sorted(tuple(sorted(r))for r in q['minimizing_ranks'])
    scout=read('fresh_nonarithmetic_release12_47_52_pentagon_20261017')
    gate=read('fresh_nonarithmetic_pentagon_two_center_gate_20261017')
    assert scout['improvements']==[] and scout['deleted']==gate['deleted_walls']==[12,47,52]
    assert scout['status']=='BOUNDED_NUMERICAL_TWO_CENTER_PENTAGON_REPLACEMENT'
    assert gate['status']=='PASS_EXACT_ALL_LABEL_PENTAGON_AT_MOST_TWO_CENTER_GATE'
    assert gate['counts']=={'four_contact_low_recipes':32,'complete_high_label_tails':6,'all_right_branches':2}
    return dict(scope='Independent Fraction arithmetic and frozen-data consistency only: no producer import, geometric rederivation, old matrix rebuild or historical large scan.',
        new_cap_columns=15,new_cap_violations=violations,changed_simplex_rows=sorted(da),
        singleton_increments=added,old_literal_terms_byte_preserved=oldterms,
        new_literal_terms=oldterms+len(added),changed_vertex_rows=sorted(changed(mu0,mu)),
        all_vertex_masses_monotone=True,forced_incidence_checks=52,
        all17_joint_families_preserved_by_singleton_zero_identity=True,
        exact_input_bound=base['bound'],exact_repaired_bound=str(bound),partial_model_bound_above12=True,
        optimizer_invocations=0,repeated_historical_scans=0,
        unrestricted_dimension_bound=156,minimum_cusp_factor_rank9_dimension_bound=29,
        reduced_profile_scalar_checks=924,construction_improvements=[],new_finite_volume_example=None)

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
    checks=small_checks();print('ARCHIVE_FRACTION_CHECKS_PASS',checks['forced_incidence_checks'],flush=True)
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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_DISTANCE14_AND_SMALL_CUSP_FACTOR_FRONTIER',
                scope='Archive integrity, exact edge-cap repairs and reviewed dimension156/rank9-cusp composition metadata; not formalization or a repeated historical mathematical verification.',
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
