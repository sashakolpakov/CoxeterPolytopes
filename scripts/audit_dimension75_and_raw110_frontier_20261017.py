#!/usr/bin/env python3
"""Incremental dimension75 and raw110 archive; no underlying census is repeated."""
import argparse,ast,json,re
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from itertools import combinations
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A
ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/dimension75_and_raw110_frontier_20261017.md'
PREVIOUS='results/spectator_spectrum_and_c3_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_spectator_spectrum_and_c3_frontier_20261017.py'
OUT=ROOT/'results/dimension75_and_raw110_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_finite8_a4_d13_neighbor_coupon_budget_20261017',
 'fresh_d17_singleton_pair_coupon_20261017',
 'fresh_ideal5_piecewise_dimension75_20261017',
 'fresh_ideal5_piecewise_dimension75_20261017_review',
 'fresh_raw110_mixed_actual_corners_20261017',
)
REPORT_STEMS=STEMS
EXTRA_REPORTS=()
PENDING_INPUTS=()
COMPLETED_UNRESTRICTED_BOUND=75
@lru_cache(None)
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())
def positive_definite(g):
    g=[row[:]for row in g]
    for k in range(len(g)):
        pivot=g[k][k];assert pivot>0
        for i in range(k+1,len(g)):
            for j in range(k+1,len(g)):g[i][j]-=g[i][k]*g[k][j]/pivot
    return True

def coupon_checks():
    d=read('fresh_finite8_a4_d13_neighbor_coupon_budget_20261017')
    assert d['finite_vertex_budget']=='8n'and d['dimension_maximum']==80
    assert d['counts']==dict(distinct_endpoints_all_dimensions=3465,distinct_endpoints_at_80=138,
      endpoints_at_80_with_duplicate_zero_endpoints=140,positive_coefficient_cases_at_80=38)
    assert F(d['component_count']['D13_unit_cost'])==13*12//2==78
    assert F(d['component_count']['D13_base_slack'])==8*13-78==26
    assert F(d['component_count']['additional_D13_slack'])==18
    assert F(d['component_count']['cost_per_A4_D13_component_pair'])==2*F(1,2)==1
    assert d['analytic_dominance']['maximum_small_rank_if_b_positive']==80-13==67
    assert F(d['analytic_dominance']['per_D13_difference_upper'])==F(67,4)-18==F(-5,4)
    expected=set()
    for n in range(1,81):
        for t in range(n//11+1):
            for b in range(t+1):
                for u in range((n-11*t-2*b)//17+1):
                    for a in {0,n-11*t-2*b-17*u}:expected.add((n,a,t,b,u))
    found=set();counts=Counter()
    for r in d['dimension_endpoints']:
        n,a,t,b,u=(r[k]for k in('dimension','small_rank','old_large_count','D13_count','A17_count'))
        key=n,a,t,b,u;assert key not in found;found.add(key);counts[str(n)]+=1
        old=a*(F(3*t,2)+u-6)-8*t;delta=b*(F(a,4)-18)
        assert old==F(r['old_excess'])<=0 and delta==F(r['additional_excess'])
        assert old+delta==F(r['excess'])<=0
        if b:assert a<=67 and delta<=-F(5*b,4)
        else:assert delta==0
    assert found==expected and len(found)==3465 and counts==d['dimension_endpoint_counts']
    assert counts['80']==138 and len(d['positive_coefficient_cases_at_80'])==38
    return dict(endpoints=3465,endpoint_dimensions=80,dimension80_endpoints=138,finite_budget='8n',
      scope='Independent Fraction check of the saved complete integer endpoint domain and analytic slack identity; all-rank base costs remain the inherited theorem.')


def local_pair_checks():
    d=read('fresh_d17_singleton_pair_coupon_20261017')
    old=read('fresh_distance8_ideal5_small_factor_20261017')
    target=next(r for r in old['records']if r['kind']=='D'and r['order']==17)
    paid=set(map(tuple,d['new_paid_pairs']));before=set(map(tuple,target['retained_paid_pairs']))
    assert d['affine_type']=='D16'and d['diagram_order']==17
    assert len(before)==d['old_cost']==127 and len(paid)==d['new_cost']==126
    assert before-paid=={(1,11)}and paid<before and d['removed_pair']==[1,11]
    assert d['affected_overlap_cases']==[2518]
    assert len(d['automorphisms'])==8 and len(d['extra_pair_orbits'])==6
    targetG=target['target'];pairset=set();contextcounts=[]
    for a in d['automorphisms']:
        assert sorted(a)==list(range(17))
        assert all(targetG[i][j]==targetG[a[i]][a[j]]for i in range(17)for j in range(17))
        assert tuple(sorted((a[1],a[11])))==(1,11)
        assert {tuple(sorted((a[i],a[j])))for i,j in paid}==paid
    cases={r['index']:r for r in read('fresh_distance8_high_overlap_20261017')['isolation']['cases']}
    for orbit in d['extra_pair_orbits']:
        rep=orbit['representative'];pairs=set(map(tuple,orbit['pairs']))
        assert pairs=={tuple(sorted((a[rep[0]],a[rep[1]])))for a in d['automorphisms']}
        assert not(pairset&pairs);pairset|=pairs
        for r in orbit['contexts']:
            c=cases[r['index']]
            assert all(r[k]==c[k]for k in ('left','right','match','preferred_corner','component_order_bound','method'))
        contextcounts.append(len(orbit['contexts']))
    assert len(pairset)==17 and contextcounts==[2,2,1,4,4,16]
    oldrule=next(r for r in read('fresh_ideal5_small_factor_payments_20261017')['records']if r['kind']=='D'and r['order']==17)
    basic=set(map(tuple,oldrule['basic_pairs']))
    assert len(basic)==d['basic_pair_count']==110 and pairset==before-basic and basic<=paid
    singleton=next(r for r in d['extra_pair_orbits']if r['representative']==[1,11])
    assert singleton['distance']==10 and [r['index']for r in singleton['contexts']]==[2518]
    affected=[]
    for r in target['case_embeddings']:
        if [1,11]in r['marked_images']:affected.append(r['case'])
    assert affected==[2518]
    assert set(map(tuple,d['triangle_pairs_retained']))<=paid
    triangle=next(r for r in read('fresh_finite8_ideal5_combined_payments_20261017')['triangle_pair_retention']if r['kind']=='D'and r['order']==17)
    assert d['triangle_pairs_retained']==triangle['required_pairs']
    c=d['literal_core'];g=[list(map(F,row))for row in c['gram']]
    assert len(g)==19 and c['signature']==[17,1,1]
    radical=list(map(F,c['radical']));assert all(sum((x*y for x,y in zip(row,radical)),F())==0 for row in g)
    assert c['omitted_radical_coefficients']=={'3':-2,'13':2}and radical[3]==-2 and radical[13]==2
    ids=c['affine_nodes'];kernel=list(map(F,c['affine_positive_kernel']))
    assert len(ids)==len(kernel)==17 and all(x>0 for x in kernel)
    assert all(sum((g[i][j]*v for j,v in zip(ids,kernel)),F())==0 for i in ids)
    white=c['white'];assert len(white)==15 and positive_definite([[g[i][j]for j in white]for i in white])
    for mark in [3,4,13,14]:
        edge=white+[mark];positive_definite([[g[i][j]for j in edge]for i in edge])
    assert len(c['corners'])==4
    for r in c['corners']:
        if r['kind']=='finite':
            corner=white+r['marks'];positive_definite([[g[i][j]for j in corner]for i in corner])
            for block in r['components']:
                nodes=block['nodes'];m=block['mark'];adj={i:[j for j in nodes if i!=j and g[i][j]]for i in nodes}
                if block['type']=='A4':assert len(nodes)==4 and len(adj[m])==1 and block['role']=='endpoint'
                elif block['type']=='D13':
                    assert len(nodes)==13 and block['role']=='long_arm_neighbor_of_fork'
                    forks=[i for i in nodes if len(adj[i])==3];assert len(forks)==1
                    assert m in adj[forks[0]]and len(adj[m])==2
                else:raise AssertionError(block['type'])
    assert [F(r['half_coupon'])for r in c['corners']]==[F(0),F(1,2),F(1,2),F(0)]
    assert sum((F(r['half_coupon'])for r in c['corners']),F())==1
    assert d['new_finite_coupon']['weight']=='1/2'
    return dict(removed_pair=[1,11],old_cost=127,new_cost=126,automorphisms=8,extra_orbits=6,
      affected_context=2518,strict_edge_checks=4,finite_corner_checks=3,new_half_total='1',
      scope='Saved set/map interfaces, rational kernels, strict principal Grams and literal coupon roles. Universal extra-white sealing remains the reviewed local proof.')


def dimension_checks():
    p=read('fresh_ideal5_piecewise_dimension75_20261017');q=read('fresh_ideal5_piecewise_dimension75_20261017_review')
    assert p['dimension_bound']==q['dimension_bound']==75 and p['inherited_dimension_bound']==77
    assert p['excluded_dimensions']==[76,77]
    assert p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    assert p['finite_vertex_budget']=='8n'
    assert p['ideal_rank_envelope']['16']=='126/17'
    binder=read('fresh_finite8_ideal5_combined_payments_20261017')
    envelope=dict(binder['ideal_rank_envelope']);envelope['16']='126/17'
    costs={m:dict(v)for m,v in binder['small_order_costs'].items()};costs['17']['D']=126
    assert p['ideal_rank_envelope']==envelope and p['small_order_costs']==costs
    assert p['local_change']['removed_pair']==[1,11]and p['local_change']['affected_overlap_case']==2518
    assert p['local_change']['finite_coupon']==read('fresh_d17_singleton_pair_coupon_20261017')['new_finite_coupon']
    assert p['unchanged_compact_arithmetic_bound']==21 and p['unchanged_arbitrary_compact_bound']==29
    assert p['effective_compact_facet_ambient_bound']==75 and p['older_auxiliary_compact_facet_bound']==79
    pp={tuple(r['ranks']):r for r in p['reduced_profiles']};qq={tuple(r['ranks']):r for r in q['profiles']}
    assert len(pp)==len(qq)==2164 and set(pp)==set(qq)
    counts=Counter()
    def cost(r):return F(r,2)if r<=12 else{13:F(87,14),14:F(20,3),15:F(111,16),16:F(126,17)}.get(r,F(163,23))
    for ranks,r in pp.items():
        other=qq[ranks];assert all(other[k]==v for k,v in r.items())
        d=sum(ranks);h=r['h_polynomial'];V=sum(h);G=sum(v*min(i,d-i)for i,v in enumerate(h))
        assert len(h)==d+1 and h==h[::-1]and V==r['link_vertices']and G==r['G']==other['independent_cumulative_G']
        curvature=F(5*G-d*V,2*V);budget=sum(map(cost,ranks),F())
        assert curvature==F(r['normalized_curvature'])and budget==F(r['normalized_cost'])
        assert curvature-budget==F(r['margin'])>0;counts[d+1]+=1
    assert counts=={76:1053,77:1111}and p['dimension_minima']==q['dimension_minima']
    assert [r['minimum_margin']for r in p['dimension_minima']]==['11615/448188','41558/485537']
    for r in p['dimension_minima']:
        n=r['ambient_dimension'];rows=[v for k,v in pp.items()if sum(k)==n-1]
        minimum=min(rows,key=lambda x:F(x['margin']))
        assert len(rows)==r['profiles']and minimum['ranks']==r['minimizing_ranks']and minimum['margin']==r['minimum_margin']
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta-8*n==F(r['finite_slack'])>0 and r['strict_global_slack']==5*(n-1)
    assert all(q['boundary'][k]==v for k,v in p['boundary_at_n75'].items())
    assert F(p['boundary_at_n75']['margin'])==F(-16873,171366)
    assert F(p['concentration']['double13_budget_increase'])==F(5,21)>0
    assert F(p['concentration']['double15_budget_increase'])==F(83,408)>0
    assert {r['kind']:F(r['maximum_ratio'])for r in q['independent_affine_feature_maxima']}=={'B':F(160,23),'C':F(75,11),'D':F(163,23)}
    review=q['composition_review'];prior_review=read('fresh_ideal5_piecewise_dimension77_20261017_review')
    branches=[b for r in q['independent_affine_feature_maxima']for b in r['linear_branches']]
    assert len(branches)==15 and branches==review['independent_all_rank_branches']
    assert branches==prior_review['composition_review']['independent_linear_branches']
    assert review['small_factor_costs']==p['small_order_costs']and review['elementary_triangle_inclusions']==16
    assert review['removed_pair']==[1,11]and review['affected_saved_map_cases']==[2518]
    local=read('fresh_d17_singleton_pair_coupon_20261017')['literal_core']
    literal=review['independently_checked_literal_core']
    assert literal['core_signature']==local['signature']and literal['radical']==local['radical']
    assert literal['positive_kernel']==local['affine_positive_kernel']
    assert literal['white_signature']==[15,0,0]and literal['affine_signature']==[16,0,1]
    assert literal['strict_edge_checks']==4 and literal['independently_checked_coupon_corners']==[[3,14],[4,13]]
    assert F(literal['total_new_payment'])==1
    return dict(bound=75,excluded_dimensions=[76,77],matched_profiles=2164,profile_counts=[1053,1111],
      minimum_margins=['11615/448188','41558/485537'],finite_slacks=['57','68'],boundary75_margin='-16873/171366',
      unchanged_all_rank_branches=15,triangle_inclusions=16,independent_literal_core_interface=True,
      scope='Primary/independent-review equality for every saved scalar row plus Fraction formula/minimum checks; no repetition of inherited all-rank geometric classifications.')

def raw110_checks():
    d=read('fresh_raw110_mixed_actual_corners_20261017');counts=d['counts']
    assert(counts['actual_Q_classes'],counts['corner_occurrences'],counts['literal_orientations'],counts['forests'])==(189,1612,30109,491)
    assert counts['columns']==len(d['columns'])==263 and counts['violated_columns']==21 and counts['vertex_mass_increases']==14
    assert counts['mixed_columns']==45 and counts['violated_mixed_columns']==17 and counts['mixed_additional_vertex_raises']==11
    old=read('fresh_weighted_actual_corner_capacities_20261017')
    assert d['face_multipliers']==old['face_multipliers']and d['initial_vertex_multipliers']==old['repairs'][0]['vertex_multipliers']
    del old
    faces={r['face']:r for r in d['selected_faces']};origins={}
    for r in faces.values():
        for a in r['corners']:
            key=(a['face'],tuple(a['omitted']),a['vertex']);assert key not in origins;origins[key]=a
    assert len(origins)==1612
    total=reused=excluded=0;seen_cols=set()
    for row in d['forest_records']:
        previous=[];owned_total=owned_excluded=0
        for dom in row['inherited_pair_domains']:
            mp=dict(dom['new_to_old_star_indices']);assert len(mp)==len(set(mp.values()))
            cur=set(mp);owned=len(cur)*(len(cur)-1)//2
            for k in range(1,len(previous)+1):
                for sub in combinations(previous,k):
                    m=len(cur.intersection(*sub));owned+=(-1)**k*m*(m-1)//2
            assert owned==dom['pairs_reused']and 0<=dom['exclusions_reused']<=owned
            previous.append(cur);owned_total+=owned;owned_excluded+=dom['exclusions_reused']
        n=len(row['literal_stars']);pairs=n*(n-1)//2
        assert pairs==row['tested_literal_pairs']==owned_total+row['newly_tested_pairs']
        assert pairs-row['unknown_literal_pairs']==owned_excluded+len(row['incompatibility_certificates'])
        groups=row['corner_orbits'];assert sorted(i for g in groups for i in g['literal_stars'])==list(range(n))
        for g in groups:
            coeff=Counter()
            for a in g['origins']:
                assert origins[a['face'],tuple(a['omitted']),a['vertex']]==a;coeff[a['face']]+=1
            assert dict(coeff)==dict(g['positive_face_coefficients'])
        for ci in row['columns']:
            assert ci not in seen_cols;seen_cols.add(ci);c=d['columns'][ci]
            assert c['vertex']==row['vertex']
            indices=c['orbit_clique']if c['kind']=='capacity_one_orbit_clique'else c['orbit_indices']
            coeff=Counter()
            for j in indices:coeff.update(dict(groups[j]['positive_face_coefficients']))
            assert sorted(coeff.items())==list(map(tuple,c['positive_face_coefficients']))
            if c['kind']=='capacity_one_orbit_clique':assert c['vertex_cost']==1
            else:
                assert c['vertex_cost']==len(c['literal_coloring'])
                assert sorted(i for part in c['literal_coloring']for i in part)==sorted(i for j in indices for i in groups[j]['literal_stars'])
        total+=pairs;reused+=owned_total;excluded+=owned_excluded
    assert seen_cols==set(range(263))and total==counts['literal_pairs_tested']==4254990
    assert reused==counts['pairs_reused']==4182709 and excluded==counts['exclusions_reused']==4098437
    alpha={i:F(a)for i,a in d['face_multipliers']};initial={i:F(a)for i,a in d['initial_vertex_multipliers']}
    expected=dict(initial);single=dict(initial);violations=mixed=0
    for c in d['columns']:
        load=sum((F(m)*alpha[i]for i,m in c['positive_face_coefficients']),F());v=c['vertex'];cost=F(c['vertex_cost'])
        assert load==F(c['load']) and cost>0
        violation=load>cost*initial.get(v,F());assert violation==c['violated'];violations+=violation
        expected[v]=max(expected.get(v,F()),load/cost)
        if c['cross_cohort']:mixed+=violation
        else:single[v]=max(single.get(v,F()),load/cost)
    final={i:F(a)for i,a in d['vertex_multipliers']};cert=d['fixed_witness_repair']
    assert expected==final and violations==21 and mixed==17
    assert sum(expected.values())==F(cert['normalizer'])
    objective=F(cert['objective_numerator']);full=objective/sum(final.values());without=objective/sum(single.values())
    assert full==F(cert['bound'])<12<without==F(cert['without_cross_cohort_bound'])
    assert sum(final.get(v,F())>initial.get(v,F())for v in final)==14
    return dict(classes=189,corners=1612,columns=263,full_ratio=str(full),without_cross_ratio=str(without),
      saved_pairs=total,first_owned_pairs=reused,scope='Saved first-owner inclusion-exclusion, occurrence coefficients, coloring partitions and exact loads. No repeated wall-pair tests or independent graph-clique arithmetic.')


def small_checks():
    return dict(completed_unrestricted_bound=75,coupon=coupon_checks(),local_pair=local_pair_checks(),
      dimension75=dimension_checks(),raw110=raw110_checks(),finite_volume_example=False,
      scope='Completed75 composition and its independent scalar review plus selected raw110 fixed-witness capacities. No new construction, arithmetic dimension bound, optimizer or target12 primal feasibility.')
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
    return dict(status='PASS_INCREMENTAL_DIMENSION75_AND_RAW110_ARCHIVE',
                scope='Preserves the completed unrestricted77 and13-package spectrum/C3 checkpoint, binds the completed unrestricted75 proof and its independent review, and records raw110 fixed-witness capacities. No new finite-volume example or arithmetic dimension theorem. Saved checks do not replace the underlying geometric or pair censuses.',
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
