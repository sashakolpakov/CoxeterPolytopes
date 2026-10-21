#!/usr/bin/env python3
"""Incremental archive of ideal5 payments, scalar comparison and corner geometry."""
import argparse,ast,json,re
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from itertools import combinations
from math import prod
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/ideal5_local_and_raw23_model_frontier_20261017.md'
PREVIOUS='results/dimension80_and_joint_corner_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_dimension80_and_joint_corner_frontier_20261017.py'
OUT=ROOT/'results/ideal5_local_and_raw23_model_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_finite8_leaf9_half_coupon_budget_20261017',
 'fresh_finite8_a17_d11_coupon_budget_20261017',
 'fresh_ideal5_finite8_triangle_payments_20261017',
 'fresh_finite8_cycle_replacement_20261017',
 'fresh_finite8_ideal5_local_replacement_20261017',
 'fresh_small_spectator_difference_gap_20261017',
 'fresh_cycle14_ideal5_transfer_20261017',
 'fresh_affine_d13_long_pair_transfer_20261017',
 'fresh_ideal5_small_factor_payments_20261017',
 'fresh_ideal5_combined_payments_20261017',
 'fresh_ideal5_piecewise_dimension77_20261017',
 'fresh_ideal5_piecewise_dimension77_20261017_review',
 'fresh_mixed_raw87_large_core_corners_20261017',
 'fresh_raw23_mixed_current_corners_20261017',
 'fresh_raw23_mixed_full_model_20261017',
 'fresh_raw23_mixed_full_model_20261017_support',
 'fresh_current220_actual_corner_union_20261017',
 'fresh_r_circle_admissibility_20261017',
 'fresh_joint_rz_circle_scout_20261017',
 'fresh_joint_rz_high_z8_gate_20261017',
 'fresh_joint_rz_low_z8_gate_20261017',
 'fresh_dotted_c10_domain_20261017',
 'fresh_dotted_c10_joint_scout_20261017',
 'fresh_dotted_c10_split_ac_scout_20261017',
)
REPORT_STEMS=STEMS
EXTRA_REPORTS=()
PENDING_INPUTS=()
COMPLETED_UNRESTRICTED_BOUND=77

@lru_cache(None)
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())


def payment_and_dimension_checks():
    census=read('fresh_ideal5_small_factor_payments_20261017')
    combined=read('fresh_ideal5_combined_payments_20261017')
    primary=read('fresh_ideal5_piecewise_dimension77_20261017')
    review=read('fresh_ideal5_piecewise_dimension77_20261017_review')
    assert len(census['patterns'])==571 and len(census['records'])==16
    assert census['order_bound_counts']=={'14':134,'15':177,'16':180,'17':80}
    expected_costs={'14':{'A':70,'B':87,'C':84,'D':82},
      '15':{'A':75,'B':100,'C':95,'D':96},'16':{'A':80,'B':107,'C':103,'D':111},
      '17':{'A':85,'B':111,'C':108,'D':127}}
    assert combined['small_order_costs']==expected_costs
    assert combined['small_order_pair_sets']==census['records']
    sets={}
    for r in census['records']:
        basic={tuple(p)for p in r['basic_pairs']};extra={tuple(p)for p in r['extra_pairs']}
        dropped={tuple(p)for p in r['dropped_pairs']};paid={tuple(p)for p in r['paid_pairs']}
        assert paid==(basic|extra)-dropped
        assert len(paid)==len(r['paid_pairs'])==r['cost']==expected_costs[str(r['order'])][r['kind']]
        assert F(r['normalized_cost'])==F(r['cost'],r['order'])
        assert all(0<=a<b<r['order']for a,b in paid)
        sets[r['order'],r['kind']]=paid
    for r in combined['triangle_pair_retention']:
        assert r['all_retained'] and {tuple(p)for p in r['elementary_triangle_pairs']}<=sets[r['order'],r['kind']]
    assert F(combined['triangle_minimum'])==2 and F(combined['quadrilateral_minimum'])==1
    assert combined['finite_vertex_budget']=='9n' and combined['finite_budget_dimension_maximum']==80
    assert combined['low_core_maximum_factor_order']==12 and combined['inherited_high_overlap_cases']==949
    for r in combined['finite9_increment']['rows']:
        delta=11 if r['type']=='A_or_B'else 12
        assert r['count']==len(r['added_pairs'])==r['rank']-delta<=r['rank']
    for kind,order in [('D13',12),('D14',13)]:
        assert len(combined['immediate_ideal_alternatives'][kind])==2
        assert all(r['order']==order for r in combined['immediate_ideal_alternatives'][kind])
    replacements=combined['old_order16_17_replacements']
    assert replacements['all_terminal_payments_finite']
    assert len(replacements['order16_cases'])==len(replacements['order17_cases'])==9
    assert primary['ideal_rank_envelope']==combined['ideal_rank_envelope']
    assert primary['dimension_bound']==77 and primary['inherited_dimension_bound']==80
    assert primary['excluded_dimensions']==[78,79,80] and primary['unrestricted_finite_volume']
    assert primary['arithmetic_hypothesis']is False and primary['cusp_restriction']is None
    assert primary['effective_compact_facet_ambient_bound']==77 and primary['unchanged_compact_facet_ambient_bound']==79
    assert review['dimension_bound']==77 and review['inherited_dimension_bound']==80
    scalars=review
    p={tuple(r['ranks']):r for r in primary['reduced_profiles']}
    q={tuple(r['ranks']):r for r in scalars['profiles']}
    assert len(p)==len(q)==3722 and set(p)==set(q)
    def cost(r):return F(r,2)if r<=12 else {13:F(87,14),14:F(20,3),15:F(111,16),16:F(127,17)}.get(r,F(163,23))
    groups={n:[]for n in(78,79,80)}
    for key,r in p.items():
        other=q[key]
        assert all(other[k]==v for k,v in r.items())
        ranks=r['ranks'];d=sum(ranks);h=r['h_polynomial'];V=prod(a+1 for a in ranks)
        assert d==r['link_dimension'] and V==r['link_vertices']==sum(h)
        assert len(h)==d+1 and h==h[::-1]
        g=sum(a*min(j,d-j)for j,a in enumerate(h));assert g==r['G']==other['independent_cumulative_G']
        curvature=F(5*g-d*V,2*V);budget=sum(map(cost,ranks),F())
        assert curvature==F(r['normalized_curvature']) and budget==F(r['normalized_cost'])
        assert F(r['margin'])==curvature-budget>0
        groups[d+1].append(r)
    for r,s in zip(primary['dimension_minima'],scalars['dimension_minima']):
        assert all(r[k]==v for k,v in s.items())
        n=r['ambient_dimension'];rows=groups[n];minimum=min(rows,key=lambda x:F(x['margin']))
        assert len(rows)==r['profiles'] and minimum['ranks']==r['minimizing_ranks']
        assert F(minimum['margin'])==F(r['minimum_margin'])
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert F(r['finite_slack'])==beta-9*n>=0 and r['strict_global_slack']==5*(n-1)>0
    assert [len(groups[n])for n in(78,79,80)]==[1173,1239,1310]
    assert min(F(r['margin'])for r in p.values())==F(409611,2171546)
    assert all(scalars['boundary'][k]==v for k,v in primary['boundary_at_n77'].items())
    assert F(primary['boundary_at_n77']['margin'])==F(-16664,1085773)
    assert F(primary['concentration']['double13_budget_increase'])==F(5,21)
    assert F(primary['concentration']['double15_budget_increase'])==F(107,408)
    branches=combined['large_order_linear_branches']
    independent=scalars['composition_review']['independent_linear_branches']
    feature_branches=[b for r in scalars['independent_affine_feature_maxima']for b in r['linear_branches']]
    fields=('kind','lower_order','upper_order','slope','intercept')
    assert {tuple(r[k]for k in fields)for r in branches}=={tuple(r[k]for k in fields)for r in independent}
    assert {tuple(r[k]for k in fields)for r in independent}=={tuple(r[k]for k in fields)for r in feature_branches}
    assert len(branches)==15
    for r in branches:
        slope=F(163,23)-r['slope'];lo,hi=r['lower_order'],r['upper_order']
        at=lo if slope>=0 else hi;assert at is not None
        assert at==r['minimum_at'] and slope*at-r['intercept']==F(r['minimum_slack'])>=0
    return dict(preferred_patterns=571,small_factor_pair_sets=16,matched_scalar_rows=3722,
      matched_all_rank_branches=15,dimension_profile_counts=[1173,1239,1310],
      minimum_ideal_margin='409611/2171546',finite_slacks=['0','21/2','20'],
      boundary77_margin='-16664/1085773',dimension77_review_bound=True)


def current220_checks():
    doc=read('fresh_current220_actual_corner_union_20261017');counts=doc['counts']
    assert counts['actual_Q_classes']==220 and counts['corner_occurrences']==1860
    assert counts['literal_orientations']==33526 and counts['forests']==510
    assert counts['columns']==len(doc['columns'])==443 and counts['violated_columns']==61
    assert counts['vertex_mass_increases']==46
    assert doc['restricted_witness_repair']['all_historical_capacity_columns_checked_at_initial_mu']==4266
    schema=doc['certificate_schema']['inherited_pair_domains']
    assert 'FIRST'in schema and 'unknown'in schema and 'Later maps cannot override'in schema
    # Check complete-domain overlap counts without recomputing wall pairs or
    # treating the new-only exclusion list as the full inherited relation.
    reused=excluded=total=0
    for row in doc['forest_records']:
        previous=[];owned_total=owned_excluded=0
        for dom in row['inherited_pair_domains']:
            mp=dict(dom['new_to_old_star_indices']);assert len(mp)==len(set(mp.values()))
            cur=set(mp);owned=len(cur)*(len(cur)-1)//2
            for k in range(1,len(previous)+1):
                for sub in combinations(previous,k):
                    m=len(cur.intersection(*sub));owned+=(-1)**k*m*(m-1)//2
            assert owned==dom['pairs_reused'] and 0<=dom['exclusions_reused']<=owned
            previous.append(cur);owned_total+=owned;owned_excluded+=dom['exclusions_reused']
        n=len(row['literal_stars']);pairs=n*(n-1)//2
        assert pairs==row['tested_literal_pairs']==owned_total+row['newly_tested_pairs']
        assert pairs-row['unknown_literal_pairs']==owned_excluded+len(row['incompatibility_certificates'])
        reused+=owned_total;excluded+=owned_excluded;total+=pairs
    assert total==counts['literal_pairs_tested']==4785958
    assert reused==counts['pairs_reused']==4321030 and excluded==counts['exclusions_reused']
    alpha={i:F(a)for i,a in doc['face_multipliers']}
    for name,initial_key,final_key in [('fixed_witness_repair','initial_vertex_multipliers','vertex_multipliers'),
                                     ('restricted_witness_repair','initial_vertex_multipliers','vertex_multipliers')]:
        cert=doc[name];holder=doc if name=='fixed_witness_repair'else cert
        a=alpha if name=='fixed_witness_repair'else{i:F(v)for i,v in cert['face_multipliers']}
        initial={i:F(v)for i,v in holder[initial_key]};final={i:F(v)for i,v in holder[final_key]}
        expected=dict(initial)
        for col in doc['columns']:
            load=sum((F(m)*a[i]for i,m in col['positive_face_coefficients']),F())
            assert load==F(col['load']) and col['vertex_cost']>0
            v=col['vertex'];expected[v]=max(expected.get(v,F()),load/col['vertex_cost'])
        assert final==expected
        assert sum(final.values())==F(cert['normalizer'])
        assert F(cert['objective_numerator'])/sum(final.values())==F(cert['bound'])
    full=F(doc['fixed_witness_repair']['bound']);restricted=F(doc['restricted_witness_repair']['bound'])
    assert full==F(893073974017141280,72281373619198223)>12
    assert restricted==F(2600114203061924640,216842261971938119)<12
    return dict(classes=220,corners=1860,new_columns=443,full_ratio=str(full),restricted_ratio=str(restricted),
      historical_restricted_capacity_checks=4266,pair_scope='Saved ordered-domain ownership counts only; no physical pair census or compatibility recomputation.')


def construction_successor_checks():
    joint=read('fresh_dotted_c10_joint_scout_20261017');split=read('fresh_dotted_c10_split_ac_scout_20261017')
    def cone(record):
        bd,ic=record['boundary'],record['PD_incidence']
        assert not bd['warnings'] and ic['status']=='numerically_separated_acute_incidence'
        assert bd['open_ends']==ic['open_ends'] and bd['counts']['finite']==ic['finite']
        assert {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
    assert len(joint['whole_control_C47_events'])==8 and not joint['unresolved']
    for r in joint['whole_control_C47_events']:cone(r)
    cone(joint['whole_control_C3'])
    assert [r['boundary']['open_ends']for r in joint['whole_control_C47_events']]==[10,10,10,13,13,13,13,13]
    assert joint['whole_control_C3']['boundary']['open_ends']==11
    assert len(joint['high_AC_finite_A53_events'])==13
    assert sum(r['X_trials']for r in joint['high_AC_finite_A53_events'])==234
    assert all(not r['native_X']for r in joint['high_AC_finite_A53_events'])
    assert not joint['coherent_six_prefixes'] and not joint['optimizations']
    assert len(split['four_prefix_events'])==205 and not split['seven_prefix_extensions'] and not split['unresolved']
    histogram={}
    for r in split['four_prefix_events']:
        cone(r);bd=r['boundary'];assert bd['counts']=={'finite':22,'exterior':8}and bd['open_ends']==12
        jr=[q for q in bd['rays']if 69 in q['active']]
        assert len(jr)==4 and all(q['kind']=='finite'for q in jr)
        key=tuple(r['profiles']);histogram[key]=histogram.get(key,0)+1
    assert histogram=={(0,0,0,0):138,(0,0,4,0):36,(0,0,1,0):29,(0,0,1,1):1,(0,0,0,1):1}
    assert split['counts']['Z_equations']==15652 and split['counts']['Z_branches']==31304
    for d in(joint,split):
        cost=d['exact_old_deletion_cost'];assert cost['deleted_walls']==[12,47,52,64,65]
        assert cost['lost_finite_count']==11 and cost['compact_faces_by_dimension']=={'2':6}
    return dict(whole_control_C47_events=8,high_AC_event_pairs=13,high_AC_X_recipes=234,
      split_AC_events=205,split_AC_scores=[22,8,12],whole_J_quadrilaterals=205,
      split_AC_low_Z_branches=31304,numerical_only=True,
      scope='Saved cone/PD support agreement and counts only; no numerical root/cone rerun, physical re-enclosure or exact chamber promotion.')

def small_checks():
    # Cross-package invariants only. These do not duplicate the underlying
    # geometric classifications, full model solve or physical enclosures.
    budget=read('fresh_finite8_a17_d11_coupon_budget_20261017')
    assert budget['finite_vertex_budget']=='8n' and budget['dimension_maximum']==80
    assert len(budget['dimension_endpoints'])==1487
    for r in budget['dimension_endpoints']:
        a,t,u=r['small_rank'],r['old_large_count'],r['A17_count']
        assert a+11*t+17*u<=r['dimension']<=80
        assert F(r['excess'])==a*(F(3*t,2)+u-6)-8*t<=0
    assert budget['new_A17_component']['eligible_path_positions']==[5,11]
    spectrum=read('fresh_small_spectator_difference_gap_20261017')
    assert len(spectrum['literal_patterns'])==7 and spectrum['excluded_cycle14_energy']=='4/7'
    cycle=read('fresh_cycle14_ideal5_transfer_20261017')
    assert cycle['spectrum']==spectrum['possible_energies']
    assert cycle['cycle_pair_rule']['cost']==70
    assert [r['required_D']for r in cycle['isolated_cases']]==['4/7','1/2']
    assert sum(F(r['coupon'])for r in cycle['equality']['corners'])==1
    d13=read('fresh_affine_d13_long_pair_transfer_20261017')
    assert d13['counts']['marked_pairs']==9 and d13['counts']['extreme_A3_cases']==72
    assert len(d13['nonextreme_zero_cases'])==12
    assert all(sum(F(c['unit'])+F(c['half_coupon'])for c in r['corners'])>=1
        for r in d13['extreme_equality']['rows'])
    arith=read('fresh_raw23_mixed_current_corners_20261017')
    prior=read('fresh_mixed_raw87_large_core_corners_20261017')
    assert arith['face_multipliers']==prior['face_multipliers']
    assert arith['initial_vertex_multipliers']==prior['vertex_multipliers']
    ratio=F(arith['fixed_witness_repair']['bound'])
    assert ratio==F(1342109395844304360,116030537828718707)<12
    assert len(arith['columns'])==739 and arith['counts']['violated_columns']==272
    assert arith['counts']['vertex_mass_increases']==148
    model=read('fresh_raw23_mixed_full_model_20261017')
    assert model['solver_invocations']==1 and model['solver']['success']
    prep=model['prepared'];counts=prep['counts']
    assert prep['dimension']==13 and prep['target']==12
    assert counts['inherited_variables']==436262 and counts['distinct_added_columns']==2004
    assert counts['variables']==438266 and counts['faces']==82806 and counts['vertex_forests']==6558
    assert counts['zero_budget_signed_variables']==96176 and counts['vertex_cost_signed_variables']==52
    certificates={c['kind']:c for c in model['certificates']}
    lower=F(certificates['EXACT_PARTIAL_MODEL_DUAL']['bound'])
    upper=F(certificates['EXACT_PARTIAL_MODEL_PRIMAL']['maximum'])
    assert lower==F(893073974017141280,72072030735616573)>12
    assert upper==F(6195706752053,500000000000) and lower<=upper<F(124,10)
    circle=read('fresh_r_circle_admissibility_20261017')
    high=read('fresh_joint_rz_high_z8_gate_20261017')
    low=read('fresh_joint_rz_low_z8_gate_20261017')
    assert circle['exact_circle']['necessary_labels']==list(range(2,16))
    assert high['status']=='PASS_ALL_42_HIGH_Z8_FAMILIES_EXCLUDED' and not high['unresolved']
    assert low['counts']['frame_recipes']==392 and low['counts']['native_rejected_branches']==784
    assert len(low['profile_domain'])==126 and len(low['records'])==14
    curve=read('fresh_dotted_c10_domain_20261017')
    assert curve['retained_control_walls']==63 and curve['full_control_tail']['deleted_walls']==[]
    assert len(curve['full_control_tail']['dotted_gates'])==55
    for r in curve['retained_dotted_gates']+curve['deleted_wall_diagnostic']:
        assert F(r['strict_derivative_upper'][1])<0
    for r in curve['full_control_tail']['dotted_gates']:
        assert F(r['value_at_u3'][1])<-2
    assert F(curve['common_center']['value_at_u2'][1])<0
    assert F(curve['common_center']['strict_derivative_upper'][1])<0
    assert F(curve['saved_control_escape']['C3_product'][0])>0
    return dict(finite_budget='8n for n<=80',budget_endpoint_consistency_checks=1487,
        cycle14_retained_pairs=70,cycle14_transferred_pairs=21,
        raw23_fixed_witness=str(ratio),finite_model_interval=[str(lower),str(upper)],
        finite_model_variables=438266,finite_model_new_columns=2004,
        construction_R_labels=list(range(2,16)),construction_low_branches=784,
        dotted_curve_whole_control_tail='u>=3',dotted_curve_old_walls_preserved=63,
        dotted_curve_saved_escape_cut=True,
        payment_and_dimension=payment_and_dimension_checks(),
        current220=current220_checks(),construction_continuations=construction_successor_checks(),
        completed_unrestricted_dimension_bound=COMPLETED_UNRESTRICTED_BOUND,
        effective_compact_facet_bound=77,older_auxiliary_compact_facet_bound=79,
        compact_arithmetic_bound_unchanged=21,
        scope='Cross-package rational and scope consistency only; underlying producer replays, written geometric reviews and exact matrix certificates retain their separately recorded scopes.')


def run(allow_pending=False):
    assert allow_pending or not PENDING_INPUTS, ('archive inputs not final', PENDING_INPUTS)
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
    checks=small_checks();print('ARCHIVE_FRACTION_CHECKS_PASS',len(STEMS),flush=True)
    reports=[ROOT/CHECKPOINT]+[ROOT/f'results/{s}.md'for s in REPORT_STEMS]+[ROOT/p for p in EXTRA_REPORTS]
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
    return dict(status='PASS_INCREMENTAL_IDEAL5_LOCAL_AND_RAW23_MODEL_ARCHIVE',
                scope='Preserves the completed unrestricted80 proof and binds the complete unrestricted77 proof, its independent review, the enlarged raw23 arithmetic model, current220 capacities and scoped exact/numerical construction results. Cross-package checks do not replace geometric proofs, contact arithmetic or matrix certificate replays. No new finite-volume example or arithmetic dimension bound.',
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
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--small-only',action='store_true');ap.add_argument('--inventory',action='store_true');args=ap.parse_args()
    if args.small_only:
        print(json.dumps(small_checks(),indent=2));raise SystemExit(0)
    d=run(allow_pending=args.inventory)
    if args.inventory:
        print('READ_ONLY_ARCHIVE_INVENTORY',{k:d[k]for k in('report_count','local_link_count','source_manifest_count','new_source_manifest_count','source_reference_count','unique_source_count')},flush=True)
        raise SystemExit(0)
    if args.verify:assert d==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],{k:d[k]for k in('report_count','local_link_count','source_manifest_count','new_source_manifest_count','unique_source_count')},flush=True)
