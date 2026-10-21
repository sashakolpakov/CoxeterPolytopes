#!/usr/bin/env python3
"""Incremental archive: affine distance8 budget and mixed actual-corner model.

Rehash all frozen predecessor bytes; only parse the new dependency graph.
Small checks independently compose saved rational certificates and counts.
They do not repeat geometry censuses, optimization or numerical construction.
The completed dimension97 composition is separate from exploratory core
lemmas and the finite arithmetic-model obstruction.
"""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import lcm,prod
from pathlib import Path
import re
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/distance8_feature_and_mixed_capacity_frontier_20261017.md'
PREVIOUS='results/distance9_and_raw39_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_distance9_and_raw39_frontier_20261017.py'
OUT=ROOT/'results/distance8_feature_and_mixed_capacity_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
MIXED='fresh_mixed_raw39_actual_corners_20261017'
MODEL='fresh_mixed_raw39_full_model_20261017'
LOW='fresh_distance8_low_core_20261017'
AFFINE='fresh_affine_distance8_feature_payment_20261017'
CONSTRUCTION=tuple('fresh_two_connector_'+s+'_20261017'for s in(
    'continuous_pair','angle_deformation','native53_deformation',
    'prefix_witnesses','coupled_disk','flipped_disk'))
REPLACE='fresh_distance8_ideal_corner_replacement_20261017'
HIGH='fresh_distance8_high_disjoint_20261017'
DIM97='fresh_distance8_ideal_dimension97_20261017'
STEMS=(MIXED,MODEL,LOW,AFFINE,REPLACE,REPLACE+'_review',HIGH,HIGH+'_review',DIM97,DIM97+'_review')+CONSTRUCTION
# Construction supplies six source/data packages and one combined report.
REPORT_STEMS=(MIXED,MODEL,LOW,AFFINE,REPLACE,REPLACE+'_review',HIGH,HIGH+'_review',DIM97,DIM97+'_review',CONSTRUCTION[0])
EXTRA_REPORTS=()
PENDING_INPUTS=()

@lru_cache(None)
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())
def masses(rows):
    d={i:F(a)for i,a in rows}
    assert len(d)==len(rows)and all(a>0 for a in d.values())
    return d
def load(alpha,coeff):return sum((alpha.get(i,F())*F(a)for i,a in coeff),F())
def updates(old,new):return [[i,str(old.get(i,F())),str(new[i])]for i in sorted(new)if new[i]>old.get(i,F())]

def mixed_checks(d):
    old=read('fresh_raw39_actual_corner_compatibility_20261017')
    assert d['face_multipliers']==old['face_multipliers']
    alpha=masses(d['face_multipliers']);assert len(alpha)==4661
    initial=masses(d['initial_vertex_multipliers']);mu=masses(d['vertex_multipliers'])
    assert initial==masses(old['vertex_multipliers'])
    rebuilt=dict(initial);without=dict(initial);violations=mixed_violations=mixed_columns=0
    for c in d['columns']:
        v=c['vertex'];cost=F(c['vertex_cost']);z=load(alpha,c['positive_face_coefficients'])
        assert cost>0 and z==F(c['load'])and F(c['baseline_vertex_mass'])==initial.get(v,F())
        assert c['violated']==(z>cost*initial.get(v,F()))
        violations+=c['violated'];rebuilt[v]=max(rebuilt.get(v,F()),z/cost)
        mixed=39 in c['raw_families']and len(c['raw_families'])>1
        mixed_columns+=mixed;mixed_violations+=mixed and c['violated']
        if not mixed:without[v]=max(without.get(v,F()),z/cost)
    assert rebuilt==mu and len(d['columns'])==1440
    assert (violations,mixed_columns,mixed_violations)==(100,168,74)
    rep=d['fixed_witness_repair'];assert rep['face_masses_unchanged']and rep['optimizer_invocations']==0
    assert updates(initial,mu)==rep['vertex_updates']and len(rep['vertex_updates'])==86
    assert updates(initial,without)==rep['without_new_raw39_mixed_updates']
    assert sum(mu.get(v,F())>without.get(v,F())for v in mu)==65
    numerator=F(rep['objective_numerator'])
    assert numerator==F(old['fixed_witness_repair']['objective_numerator'])
    assert numerator/sum(initial.values(),F())==F(rep['baseline_bound'])==F(old['fixed_witness_repair']['bound'])
    assert sum(mu.values(),F())==F(rep['normalizer'])
    assert numerator/sum(mu.values(),F())==F(rep['bound'])==F(247898037929594052000000000000000,21396665711915633070810190213789)<12
    assert numerator/sum(without.values(),F())==F(rep['without_new_raw39_mixed_bound'])>12
    assert rep['bound_below12']and not d['omitted_ambiguous_selected_rows']
    assert d['raw_sources']==[23,25,39]and d['canonical_cores']==[19,21,29]
    assert d['radical_bases']==[['7/2','1/2'],['41/2','17/2']]
    # Independent squareclass determinant/norm calculation in K.
    norm=lambda a,b:a*a-5*b*b
    assert norm(F(7,2),F(1,2))==11 and norm(F(41,2),F(17,2))==59
    assert norm(F(101,22),F(39,22))==F(59,11)
    counts=d['counts'];families={}
    for raw in(23,25,39):
        fs=[f for f in d['selected_faces']if f['raw_family']==raw]
        cs=[c for f in fs for c in f['corners']]
        families[str(raw)]=dict(classes=len(fs),corners=len(cs),forests=len({c['vertex']for c in cs}),orientations=sum(c['literal_orientations']for c in cs))
    assert families==d['family_counts']
    assert [(families[str(r)]['classes'],families[str(r)]['corners'])for r in(23,25,39)]==[(63,630),(54,432),(31,124)]
    records=d['forest_records'];assert len(records)==409
    assert sum(len(r['literal_stars'])for r in records)==14232
    assert sum(len(r['corner_orbits'])for r in records)==1186
    assert sum(r['tested_literal_pairs']for r in records)==627684
    assert sum(len(r['incompatibility_certificates'])for r in records)==608207
    assert sum(r['unknown_literal_pairs']for r in records)==19477
    assert sum(r['cross_family_pairs']for r in records)==301644
    assert sum(r['cross_family_exclusions']for r in records)==289390
    assert sum(r['raw39_cross_pairs']for r in records)==92676
    assert sum(r['raw39_cross_exclusions']for r in records)==90792
    assert len(d['selected_faces'])==148 and counts['columns']==1440
    return dict(fixed_faces=4661,mixed_actual_Q_classes=148,mixed_actual_corner_occurrences=1186,
        mixed_literal_orientations=14232,mixed_literal_pairs=627684,mixed_raw39_cross_pairs=92676,
        mixed_columns=1440,genuinely_mixed_raw39_columns=168,fixed_witness_violations=100,
        fixed_witness_vertex_raises=86,fixed_witness_bound=rep['bound'],
        without_new_raw39_mixed_bound=rep['without_new_raw39_mixed_bound'],
        scope='Saved actual-corner counts and independent exact load/max-mass composition; no duplicate pair/orbit geometry census.')

def model_checks(d):
    p=d['prepared'];base=read('fresh_mixed_actual_corner_full_model_20261017')
    assert p['inherited_model_sha256']==base['prepared']['model_sha256']
    c=p['counts'];assert c==dict(active_variables=167927,distinct_added_columns=1732,
        exact_duplicate_exports=152,exported_input_columns=1884,face_matrix_nonzeros=813374,
        faces=82806,inherited_variables=432021,variables=433753,vertex_cost_signed_variables=52,
        vertex_forests=6558,zero_budget_signed_variables=96176)
    assert p['provenance_to_face']==base['prepared']['provenance_to_face']
    assert len({a for a,b in p['provenance_to_face']})==len({b for a,b in p['provenance_to_face']})==1858
    inputs=(('fresh_mixed_actual_corner_full_model_20261017_local_cones','new_exact_charge_columns','positive_simplex_coefficients',132),
        ('fresh_mixed_actual_corner_full_model_20261017_joint_ports','critical_columns','positive_face_coefficients',183),
        ('fresh_raw39_actual_corner_compatibility_20261017','columns','positive_face_coefficients',129),
        (MIXED,'columns','positive_face_coefficients',1440))
    expected={};input_alpha=None
    for tag,field,cf,n in inputs:
        x=read(tag);assert len(x[field])==n
        if input_alpha is None:input_alpha=x['face_multipliers']
        else:assert input_alpha==x['face_multipliers']
        for j,r in enumerate(x[field]):
            cost=F(r.get('vertex_cost',1));coeff=tuple((i,F(a))for i,a in r[cf])
            assert cost>0 and all(a>0 and 0<=i<c['faces']for i,a in coeff)
            key=(r['vertex'],cost,coeff);expected.setdefault(key,[]).append([tag,j])
    assert len(expected)==1732
    assert p['input_column_counts']==[[a,n]for a,b,cc,n in inputs]
    seen=set()
    for j,r in enumerate(p['new_columns']):
        coeff=tuple((i,F(a))for i,a in r['positive_face_coefficients']);cost=F(r['vertex_cost'])
        key=(r['vertex'],cost,coeff);assert key not in seen;seen.add(key)
        assert r['variable']==432021+j and r['sources']==expected[key]
        den=lcm(cost.denominator,*(a.denominator for i,a in coeff))
        assert den==r['positive_denominator']and cost*den==r['scaled_vertex_cost']
    assert seen==set(expected)
    primal=next(x for x in d['certificates']if x['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    dual=next(x for x in d['certificates']if x['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha=masses(dual['face_multipliers']);mu=masses(dual['vertex_multipliers'])
    assert len(alpha)==5059 and len(mu)==1068 and all(i<c['faces']for i in alpha)and all(i<c['vertex_forests']for i in mu)
    for r in p['new_columns']:assert load(alpha,r['positive_face_coefficients'])<=F(r['vertex_cost'])*mu.get(r['vertex'],F())
    assert sum(mu.values(),F())==F(dual['normalizer'])
    lower=F(dual['objective_numerator'])/F(dual['normalizer'])
    assert lower==F(dual['bound'])==F(44991485646505450,3603606437217783)>12 and dual['excludes_target12']
    budgets=masses(primal['vertex_budgets']);upper=max(budgets.values())
    assert upper==F(primal['maximum'])==F(6242573637821,500000000000)>12
    assert primal['meets_target12']==(upper<=12)and lower<=upper
    rv=dict(primal['rounded_variable_numerators']);assert len(rv)==len(primal['rounded_variable_numerators'])and all(0<=i<c['variables']and a>=0 for i,a in rv.items())
    assert F(primal['repair_delta'])>=0 and primal['integer_denominator']>0
    assert all(0<=i<432021 for i in primal['repair_variable_indices'])
    assert len(set(primal['repair_variable_indices']))==len(primal['repair_variable_indices'])
    assert d['solver_invocations']==1 and d['solver']['success']and d['solver']['status']==0
    return dict(finite_model_counts=c,finite_model_lower=str(lower),finite_model_upper=str(upper),
        model_new_columns_independently_recounted=1732,model_positive_faces=5059,
        scope='Exact new-column deduplication/scaling/dual loads, saved scalar and budget consistency; inherited full matrix/certificate replay is hash-bound, not rerun. Changed optimized alpha has no automatic implicit local/joint coverage.')

def low_checks(d):
    assert d['status']=='PASS_COMPLETE_DISTANCE8_SMALL_CORE_DOMAIN'
    assert len(d['small_contexts'])==178
    t=d['triangle_all_label_gate'];assert len(t['high_ports'])==3 and F(t['high_inverse_diagonal'])==0
    assert all(F(a)==0 for row in t['schur_correction']for a in row)
    dis=d['disjoint'];assert len(dis['high_ports'])==18 and len(dis['rows'])==726 and len(dis['multiple'])==321
    counts=Counter(str((r['low'][0],r['cross_label'],r['kind']))for r in dis['rows']);assert dict(counts)==dis['counts']
    assert sum(r['kind']=='two_negative'for r in dis['rows'])==668
    for r in dis['types'].values():assert r['signature'][-2:]==[1,0] # Low Schur block plus the high negative direction.
    ov=d['overlaps'];assert ov['candidate_matches']==630 and ov['consistent_matches']==len(ov['records'])==523
    assert Counter(r['kind']for r in ov['records'])==dict(inverse_rejection=422,actual_edge_PD_rejection=35,survive=66)
    geom=d['extension_geometry'];cl=d['extension_closure'];assert len(geom['cases'])==len(cl['cases'])==66 and len(geom['cores'])==14
    key=lambda r:(tuple(r['low'][:2]),tuple(r['low'][2]),r['high'][0],tuple(r['high'][1]),tuple(map(tuple,r['match'])))
    assert {key(r)for r in geom['cases']}=={key(r)for r in cl['cases']}=={key(r)for r in ov['records']if r['kind']=='survive'}
    assert Counter(r['method']for r in cl['cases'])==cl['method_counts']==dict(radical_isolation=39,no_relevant_contact=24,unique_sealed_label3_root=3)
    assert max(r['component_order_bound']for r in cl['cases'])==cl['component_order_bound']==12
    return dict(small_contexts=178,small_disjoint_records=726,small_multiple_edge_rejections=321,
        small_overlap_candidates=630,small_consistent_overlaps=523,small_surviving_overlaps=66,
        small_extension_core_types=14,small_marked_component_bound=12,
        scope='Saved complete small-core classification/count/closure composition; no new all-high payment or dimension conclusion.')

def affine_checks(d):
    assert d['graph_cost_lemma_complete']and not d['face_payment_proved_here']and d['ambient_dimension_bound']is None
    assert d['counts']==dict(classical_bfs_checks=444,small_affine_checks=63,endpoint_interval_checks=450)
    assert {(r['kind'],r['order'])for r in d['classical_checks']}=={(k,m)for k in'ABCD'for m in range(18,129)}
    for r in d['classical_checks']:
        k=r['kind'];m=r['order']
        short=8*m-{'A':0,'B':35,'C':36,'D':34}[k]
        extra={'A':0,'B':2+min(m-11,31),'C':min(m-9,32),'D':4+min(m-13,30)}[k]
        assert r['distance_through8']==short and r['feature_near_distance9']==extra
        assert r['paid_pairs']==short+extra<=min(8*m,m*(m-1)//2)==r['budget']
        assert sum(n for dist,n in r['pair_distance_histogram'])==m*(m-1)//2
        assert sum(n for dist,n in r['pair_distance_histogram']if dist<=8)==short
    for r in d['small_affine_catalogue']:
        m=r['order'];assert m<=17 and r['paid_pairs']==m*(m-1)//2==r['budget']<=8*m
    windows={'C':(16,16),'B':(15,16),'D':(15,15)}
    for k,n,total in d['interval_union_checks']:
        a,b=windows[k];assert total==len(set(range(min(a,n)))|set(range(max(0,n-b),n)))==min(n,a+b)
    bad=d['radius15_failure'];assert bad['kind']=='D'and bad['order']==100 and bad['radius']==15 and bad['paid_pairs']==802>bad['budget']==800
    return dict(classical_affine_saved_formula_checks=444,small_affine_budget_checks=63,
        independent_interval_union_checks=450,affine_cost_only=True,
        scope='Independent integer formula/count checks; all-rank graph proof and inherited classification are bound, and no triangle/quad payment is inferred.')

def construction_checks(ds):
    assert len(ds)==6
    for d,n in zip(ds[:3],(3438,4104,4968)):
        c=d['counts'];assert c.get('A_evaluations',c.get('outer_evaluations'))==n
        assert c.get('complete_conic_rectangles',c.get('conic_rectangles'))==n
        assert c['R_seeds']==len(d['records'])==12 and not d['potential_domains']
    for d in ds:
        loss=d['exact_old_deletion_cost'];assert loss['lost_finite_count']==len(loss['lost_finite_supports'])==11
        assert len(loss['compact_faces_touched'])==6 and loss['compact_faces_by_dimension']=={'2':6}
        assert d['deleted_walls']==loss['deleted_walls']==ds[0]['deleted_walls']
    p=ds[3];assert p['count']==len(p['records'])==9
    for r in p['records']:
        b=r['boundary'];assert b['counts']==dict(finite=21,exterior=9)and b['open_ends']==13 and not b['warnings']
        v=r['PD_incidence'];assert v['finite']==len(v['finite_supports'])==21 and v['geodesic_edges']==len(v['edge_supports'])==38 and v['open_ends']==13
    a,b=ds[4:]
    assert a['counts']==dict(native_R=7,B_equations=700,B_branches=1400,B_wrong_center=319,native_B=422,A_events=3376)
    assert len(a['B_bank'])==422
    assert b['counts']==dict(Y_equations=1864,native_Y=615,X_equations=4353,native_XYB=888,A_equations=3552,X_R_not_dotted=1232)
    for d in(a,b):assert not d['prefixes']and not d['complete_disks']and not d['unresolved']
    assert (b['topology']['vertices'],b['topology']['edges'],b['topology']['triangles'])==(11,25,15)
    return dict(bounded_outer_counts=[3438,4104,4968],numerical_prefixes=9,prefix_finite=21,
        prefix_exterior=9,prefix_open_ends=13,coupled_equations=700,flipped_Y_equations=1864,
        complete_disks_found=0,scope='Bounded numerical diagnostics and saved incidence/count consistency only; no continuous all-domain exclusion or new chamber.')

def high_checks(p,r):
    assert p['status']=='PASS_DISTANCE8_ALL_LABEL_HIGH_DISJOINT_EXCLUSION'
    assert r['status']=='PASS_INDEPENDENT_DISTANCE8_HIGH_DISJOINT_CERTIFICATE'
    assert p['marked_contexts']==r['marked_contexts']==384
    assert len(p['ports'])==r['ports']==2404 and sum(x['allowed']for x in p['ports'])==p['allowed_ports']==r['allowed_ports']==754
    assert len(p['single_joins'])==r['single_joins']==12135
    assert len(p['multiple_edge_subsets'])==r['multiple_subsets']==10990
    c=p['counts'];q=r['checks']
    for a,b in(('single_label_index2','single_index2_cases'),('white_rejected','white_witnesses'),
        ('single_label_bad_corner','strict_corner_witnesses'),('multiple_monotone_z','monotone_multiple_witnesses'),
        ('multiple_needing_labels','residual_multiple_subsets'),('multiple_labeled_index2','residual_index2_inertias'),
        ('multiple_labeled_edge','residual_edge_witnesses')):assert c[a]==q[b]
    assert q==dict(zero_cofactor_schur_types=262,explicit_single_inertias=4,single_index2_cases=9256,
        white_witnesses=3037,strict_corner_witnesses=20,monotone_multiple_witnesses=10987,
        residual_multiple_subsets=3,residual_index2_inertias=3,residual_edge_witnesses=3)
    assert len(p['single_types'])==q['zero_cofactor_schur_types']+q['explicit_single_inertias']==266
    assert r['input_core_checks']==dict(lorentz_cores=26,finite_deletions=142,ideal_deletions=49)
    return dict(high_marked_contexts=384,high_ports=2404,high_allowed_ports=754,
        high_single_joins=12135,high_multiple_subsets=10990,independent_checks=q,
        scope='Matched primary and independent exact disjoint-core certificates; separate local lemma, not a prerequisite for dimension97.')

def replacement_checks(p,r):
    assert p['status']=='PASS_UNIVERSAL_IDEAL_DISTANCE9_CORNER_REPLACEMENT'
    assert 'No simple-cusp or global small-affine-exclusion hypothesis.'in p['local_theorem']
    assert r['status']=='PASS_INDEPENDENT_SYMPY_IDEAL_DISTANCE9_CORNER_REPLACEMENT_REVIEW'
    x=r['reconstruction']
    for k in('cycle_profiles','cycle_rectangle_checks','tree_profiles','tree_pair_checks','cycle_uniform_gap','tree_uniform_gap'):assert p[k]==x[k]
    assert x['formal_schur_terms']==p['formal_schur_determinant']['terms']
    assert len(x['cycle_rectangle_checks'])==256 and len(x['tree_pair_checks'])==64
    assert x['cycle_uniform_gap']==['5/81','0']and x['tree_uniform_gap']==['1/9','0']
    assert x['symbolic_path_inverse_identities']==4
    assert p['hypotheses']==dict(finite_sided_finite_volume=True,actual_quadrilateral=True,
        varying_affine_factor_order_greater_than=17,marked_distance=9)
    return dict(cycle_values=256,tree_pair_values=64,cycle_gap='5/81',tree_gap='1/9',
        scope='Exact primary/independent reconstruction match and explicit universal local hypothesis check.')

def dimension97_checks(p,r):
    assert p['status']=='PASS_FINITE9_IDEAL8_FEATURE_DIMENSION97'
    assert r['status']=='PASS_INDEPENDENT_FRACTION_FINITE9_IDEAL8_DIMENSION97_REVIEW'
    assert p['dimension_bound']==r['dimension_bound']==97 and p['inherited_dimension_bound']==108
    assert p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    assert p['excluded_dimensions']==list(range(98,109))and p['convex_order_cap']==16
    assert p['unchanged_compact_facet_ambient_bound']==r['unchanged_compact_facet_ambient_bound']==79
    assert p['finite_vertex_budget']=='9n'and p['ideal_normalized_budget']=='sum_i min(8,r_i/2)'
    expected=set()
    for d in range(97,108):
        for h in(0,*range(17,d+1)):
            for a in range(16):
                b=d-a-h
                if b>=0 and b%16==0:expected.add(tuple(sorted(([a]if a else[])+[16]*(b//16)+([h]if h else[]))))
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8);assert beta>9*n
    rows=p['reduced_profiles'];other={tuple(z['ranks']):z for z in r['reconstruction']['profiles']}
    assert len(rows)==len(other)==len(expected)==957 and {tuple(z['ranks'])for z in rows}==expected
    margins=[]
    for row in rows+[p['boundary_at_n97']]:
        rs=row['ranks'];d=sum(rs);h=row['h_polynomial'];V=prod(a+1 for a in rs)
        assert len(h)==d+1 and h==h[::-1]and sum(h)==V==row['link_vertices']
        g=sum(a*min(j,d-j)for j,a in enumerate(h));assert g==row['G']
        curvature=F(5*g-d*V,2*V);cost=sum((min(F(8),F(a,2))for a in rs),F())
        assert curvature==F(row['normalized_curvature'])and cost==F(row['normalized_cost'])
        assert curvature-cost==F(row['margin'])
        if d>=97:
            assert curvature-cost>0;margins.append(curvature-cost)
            rr=other[tuple(rs)];assert rr['independent_cumulative_G']==g
            assert {k:v for k,v in rr.items()if k!='independent_cumulative_G'}==row
    assert min(margins)==F(112101,1419857)
    assert F(p['boundary_at_n97']['margin'])==F(-185352,1419857)
    assert r['reconstruction']['dimension_minima']==p['dimension_minima']
    return dict(dimension_bound=97,unchanged_compact_facet_ambient_bound=79,
        scalar_profiles=957,minimum_margin='112101/1419857',boundary_margin='-185352/1419857',
        scope='Independent small scalar/count composition with frozen geometric proofs and exact review; no high-core census rerun or arithmetic dimension claim.')

def small_checks():
    return dict(mixed=mixed_checks(read(MIXED)),finite_model=model_checks(read(MODEL)),
        small_core=low_checks(read(LOW)),affine_budget=affine_checks(read(AFFINE)),
        high_disjoint=high_checks(read(HIGH),read(HIGH+'_review')),
        ideal_replacement=replacement_checks(read(REPLACE),read(REPLACE+'_review')),
        dimension97=dimension97_checks(read(DIM97),read(DIM97+'_review')),
        construction=construction_checks([read(s)for s in CONSTRUCTION]),
        pending_inputs=list(PENDING_INPUTS),scope='Independent saved-certificate compositions of the finalized97 theorem and separate arithmetic/local/numerical packages. No formalization, optimizer, new arithmetic dimension conclusion or repeated historical geometry census.')


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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_DISTANCE8_FEATURE_AND_MIXED_CAPACITY',
                scope='Archive integrity and exact saved-certificate composition of the finalized distance8 frontier, mixed actual-corner repair and enlarged finite model, and scoped numerical two-connector diagnostics; no formalization or new arithmetic dimension conclusion.',
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
