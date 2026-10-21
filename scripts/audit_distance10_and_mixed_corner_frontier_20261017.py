#!/usr/bin/env python3
"""Incremental archive of distance10, mixed actual corners and zero-arc cap gates.

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
CHECKPOINT='results/distance10_and_mixed_corner_frontier_20261017.md'
PREVIOUS='results/distance11_and_raw25_corner_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_distance11_and_raw25_corner_frontier_20261017.py'
OUT=ROOT/'results/distance10_and_mixed_corner_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_selective_distance10_dimension116_20261017',
 'fresh_selective_distance10_dimension116_20261017_review',
 'fresh_distance10_extension_20261017',
 'fresh_mixed_raw23_raw25_corners_20261017',
 'fresh_mixed_actual_corner_full_model_20261017',
 'fresh_zero_arc_cap_continuous_20261017',
 'fresh_zero_arc_cap_two_connectors_20261017',
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
    p,review,extension,mixed,model,continuous,scout=map(read,STEMS)
    assert p['dimension_bound']==116 and p['inherited_dimension_bound']==126
    assert p['unrestricted_finite_volume'] and not p['arithmetic_hypothesis'] and p['cusp_restriction']is None
    assert p['excluded_dimensions']==list(range(117,127))
    assert review['status']=='PASS_INDEPENDENT_FRACTION_SELECTIVE_DISTANCE10_DIMENSION116_REVIEW'
    assert review['dimension_bound']==116 and review['compact_facet_ambient_bound']==87
    rows=p['reduced_profiles'];assert len(rows)==1015
    expected_profiles=set()
    for d in range(116,126):
        for k in range(d//20+1):
            for low in range(20):
                high=d-20*k-low
                if high==0 or high>=21:
                    expected_profiles.add(tuple(sorted([20]*k+([low]if low else[])+([high]if high else[]))))
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>10*n
    assert len(expected_profiles)==1015
    assert {tuple(sorted(r['ranks']))for r in rows}==expected_profiles
    def scalar(r):
        rs=r['ranks'];d=sum(rs);v=prod(t+1 for t in rs);h=[1]
        for t in rs:
            q=[0]*(len(h)+t)
            for i,x in enumerate(h):
                for j in range(t+1):q[i+j]+=x
            h=q
        assert h==h[::-1]and sum(h)==v
        G=sum(x*min(i,d-i)for i,x in enumerate(h))
        curvature=F(5*G-d*v,2*v);cost=sum((min(F(10),F(t,2))for t in rs),F())
        assert d==r['link_dimension']and v==r['link_vertices']and G==r['G']
        assert curvature==F(r['normalized_curvature'])and cost==F(r['normalized_cost'])
        assert curvature-cost==F(r['margin'])
        return curvature-cost
    assert min(map(scalar,rows))==F(145234,4084101)>0
    assert scalar(p['boundary_at_n116'])==F(-673055,16336404)
    matched={tuple(sorted(r['ranks'])):(F(r['normalized_curvature']),F(r['normalized_cost']),F(r['margin']))for r in review['budgets']['profiles']}
    assert len(matched)==1015 and matched=={tuple(sorted(r['ranks'])):(F(r['normalized_curvature']),F(r['normalized_cost']),F(r['margin']))for r in rows}
    assert F(review['budgets']['boundary_margin'])==F(-673055,16336404)
    assert p['compact_facet_ambient_bound']==87
    facet=[]
    for start in(87,88):
        values=[]
        for t in range(4):
            d=start+2*t;beta=F(d*(d-6),8)if d%2==0 else F((d-1)*(d-5),8)
            values.append(8*(beta-10*(d+1)))
        c2=(values[2]-2*values[1]+values[0])/2;co=[values[0],values[1]-values[0]-c2,c2]
        assert all(x>0 for x in co)
        assert all(values[t]==sum(co[i]*t**i for i in range(3))for t in range(4))
        facet.append(list(map(int,co)))
    assert facet==p['compact_facet_parity_slack_times8']==[[12,176,4],[96,180,4]]
    assert facet==[list(map(int,r))for r in review['budgets']['compact_facet_slacks']]
    c,j,low,iso=p['census'],p['disjoint_obstructions'],p['low_order_disjoint'],p['isolation']
    rc=review['classification'];ri=review['isolation'];rl=review['low_order_disjoint']
    assert c['marked_contexts']==rc['marked_contexts']==384 and c['allowed_ports']==rc['allowed_ports']==161
    assert len(c['ports'])==len(rc['ports'])==527 and c['ports']==rc['ports']
    assert c['overlap_generated']==rc['candidate_matchings']==22881
    assert c['consistent_overlaps']==rc['consistent_matchings']==4862
    expected={'inverse_rejection':4296,'edge_pd_rejection':254,'L10_L25_vertex':6,'L18_L25_vertex':6,
        'L22_L25_short_vertex':18,'L22_L25_long_vertex':54,'L22_L25_edge':18,
        'L25_vertex':114,'L25_edge':60,'L25_three_path':27,'L25_four_path':9}
    assert c['overlap_counts']==expected
    assert rc['overlap_counts']=={('edge_rejection'if k=='edge_pd_rejection'else k):v for k,v in expected.items()}
    def okey(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['match'])))
    left={okey(r):r for r in c['overlap_records']};right={okey(r):r for r in rc['overlap_records']}
    assert len(left)==len(right)==4862 and left.keys()==right.keys()
    for key,r in left.items():
        kind='inverse_rejection'if'inverse_rejection'in r else'edge_rejection'if'edge_pd_rejection'in r else r['surviving_type']
        assert right[key]['kind']==kind
    assert j['candidate_joins']==len(j['rows'])==len(rc['disjoint_joins'])==1478
    assert j['white_rejections']==sum('white_rejection'in r for r in j['rows'])==363
    assert j['remaining_joins']==1115 and j['remaining_context_pair_multiplicities']=={'1':905,'2':96,'3':6}
    for a,b in zip(j['rows'],rc['disjoint_joins']):
        assert a['left']==b['left']and a['right']==b['right']
        assert('white_rejection'in a)==('white_obstruction'in b)
        if'white_rejection'not in a:assert a['white_three_node_tails']==b['white_tails']
    assert len(j['underlying_inertia_types'])==len(rc['surviving_disjoint_types'])==35
    assert all(r['signature'][1]==2 for r in j['underlying_inertia_types'])
    def mkey(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['edges'])))
    assert len(j['multiple_edge_obstructions'])==len(rc['multiple_edge_obstructions'])==120
    assert{mkey(r)for r in j['multiple_edge_obstructions']}=={mkey(r)for r in rc['multiple_edge_obstructions']}
    expected_low={'proper_indefinite':209,'not_indefinite':14,'black_not_PD':4,'port_rejected':12,'survive':4}
    assert low['counts']==rl['counts']==expected_low
    assert len(low['records'])==len(rl['records'])==243 and len(low['shapes'])==3
    lowkeys={(r['shape'],tuple(r['labels']))for r in low['records']}
    from itertools import product
    assert lowkeys=={(i,labels)for i in range(3)for labels in product((3,4,5),repeat=4)}
    def lowkind(r):
        return'proper_indefinite'if'proper_indefinite'in r else'not_indefinite'if'nonnegative_full_determinant'in r else'black_not_PD'if'black_not_PD'in r else'port_rejected'if'port_rejected'in r else'survive'
    assert Counter(map(lowkind,low['records']))==expected_low
    for a,b in zip(low['records'],rl['records']):assert(a['shape'],a['labels'],lowkind(a))==(b['shape'],b['labels'],b['kind'])
    assert len(low['surviving_low_cores'])==len(rl['survivors'])==len(low['exact_join_types'])==len(rl['joined_inertias'])==4
    assert low['ordered_marked_joins']==rl['ordered_marked_joins']==24
    for a,b in zip(low['exact_join_types'],rl['joined_inertias']):
        assert(a['low_shape'],a['low_labels'],a['signature'])==(b['shape'],b['labels'],b['signature'])
        assert a['signature']==[13,2,0]and a['high_contexts']==[[1,2],[1,3],[2,3]]and a['ordered_marked_multiplicity']==6
    assert len(iso['cases'])==len(ri['preferred_corners'])==312 and iso['component_order_cap']==17
    pc={okey(r):r for r in iso['cases']};rr={okey(r):r for r in ri['preferred_corners']}
    assert len(pc)==len(rr)==312 and pc.keys()==rr.keys()
    short=extra=0
    for key,a in pc.items():
        b=rr[key]
        for field in('type','preferred_corner','omitted_black','blocks','possible_contacts','marked_component','method','component_order_bound'):assert a[field]==b[field],(key,field)
        assert a['component_order_bound']<=17
        if a['method']=='radical':assert a['shortcut_bound']==b['shortcut_bound']<=10
        else:
            for field in('end_sets','single_contact_candidates','candidate_pair_tests'):assert a[field]==b[field]
        cross=a['extra_cross_bond_checks'];assert cross['shortcut_rejections']==b['extra_cross_shortcuts']
        assert len(cross['edge_pd_z_rejections'])==len(b['extra_cross_pd_rejections'])
        short+=cross['shortcut_rejections'];extra+=len(cross['edge_pd_z_rejections'])
    assert short==iso['extra_cross_shortcuts']==ri['extra_cross_shortcuts']==21393
    assert extra==iso['extra_cross_pd_rejections']==ri['extra_cross_pd_rejections']==168
    assert review['saved_comparison']['independently_verified_saved_nonnegative_witnesses']==363+120+254+168==905
    assert review['saved_comparison']['producer_high_and_low_classification_and_all_extensions_match']
    assert extension['contexts']==18 and extension['extension_variants']==78 and extension['marked_component_order_bound']==17
    assert len(extension['families'])==2 and all(len(f['contexts'])==9 for f in extension['families'])
    assert sum(len(r['allowed_extension_variants'])for f in extension['families']for r in f['contexts'])==78
    assert all(len(v['marked_component'])<=17 for f in extension['families']for r in f['contexts']for v in r['allowed_extension_variants'])
    baseline=read(mixed['fixed_witness_repair']['baseline']);repair=mixed['fixed_witness_repair'];counts=mixed['counts']
    assert mixed['face_multipliers']==baseline['face_multipliers']and len(mixed['face_multipliers'])==4932
    alpha={i:F(x)for i,x in mixed['face_multipliers']};initial={i:F(x)for i,x in mixed['initial_vertex_multipliers']};mu={i:F(x)for i,x in mixed['vertex_multipliers']}
    assert initial=={i:F(x)for i,x in baseline['vertex_multipliers']}
    assert all(mu.get(i,F())>=x for i,x in initial.items())
    rebuilt=dict(initial);single=dict(initial);violations=mixed_violations=0
    for col in mixed['columns']:
        load=sum((alpha[i]*F(k)for i,k in col['positive_face_coefficients']),F());v=col['vertex'];cost=F(col['vertex_cost'])
        assert load==F(col['load'])and F(col['baseline_vertex_mass'])==initial.get(v,F())and cost>0
        assert col['violated']==(load>cost*initial.get(v,F()))
        assert load<=cost*mu.get(v,F());rebuilt[v]=max(rebuilt.get(v,F()),load/cost)
        if len(col['raw_families'])==1:single[v]=max(single.get(v,F()),load/cost)
        violations+=col['violated'];mixed_violations+=col['violated']and len(col['raw_families'])==2
    assert {i:x for i,x in rebuilt.items()if x}=={i:x for i,x in mu.items()if x}
    assert all(single.get(i,F())==initial.get(i,F())for i in set(single)|set(initial))
    norm=sum(mu.values(),F());numerator=F(repair['objective_numerator']);ratio=numerator/norm
    assert norm==F(repair['normalizer'])and numerator==F(baseline['fixed_witness_repair']['objective_numerator'])
    assert ratio==F(repair['bound'])==F(67685455440946980000000000000000,5705370740772271626674000656033)<12
    assert F(repair['single_family_only_bound'])==F(repair['baseline_bound'])==F(baseline['fixed_witness_repair']['bound'])>12
    assert repair['single_family_only_updates']==[]and repair['optimizer_invocations']==0 and repair['face_masses_unchanged']
    assert counts['actual_Q_classes']==112 and counts['corner_occurrences']==1010 and counts['literal_orientations']==12682
    assert mixed['common_forests']==225 and counts['forests']==379
    assert counts['cross_family_pairs']==227064 and counts['cross_family_exclusions']==216446
    assert counts['literal_pairs_tested']==535065 and counts['exact_literal_exclusions']==516727
    assert len(mixed['columns'])==counts['columns']==1251 and violations==mixed_violations==counts['violated_columns']==67
    assert sum(mu.get(i,F())>initial.get(i,F())for i in set(mu)|set(initial))==counts['vertex_mass_increases']==65
    prepared=model['prepared'];mc=prepared['counts'];oldmodel=read('fresh_actual_corner_enlarged_model_20261017')
    assert prepared['inherited_model_sha256']==oldmodel['prepared']['model_sha256']
    assert prepared['provenance_to_face']==oldmodel['prepared']['provenance_to_face']
    assert mc==dict(active_variables=166195,distinct_added_columns=1355,exact_duplicate_exports=396,
        exported_input_columns=1751,face_matrix_nonzeros=809767,faces=82806,inherited_variables=430666,
        variables=432021,vertex_cost_signed_variables=52,vertex_forests=6558,zero_budget_signed_variables=96176)
    assert model['solver_invocations']==1 and model['solver']['success']and model['solver']['time_limit_seconds']==300
    primal=next(x for x in model['certificates']if x['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    dual=next(x for x in model['certificates']if x['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert len(model['certificates'])==2 and not primal['meets_target12']and dual['excludes_target12']
    upper=F(primal['maximum']);lower=F(dual['bound'])
    assert upper==F(6253797331083,500000000000)and lower==F(661061434478917472,52852828351105989)>12
    assert lower<=upper and upper-lower<F(1,100000)
    da={i:F(a)for i,a in dual['face_multipliers']};dm={v:F(a)for v,a in dual['vertex_multipliers']}
    assert len(da)==len(dual['face_multipliers'])==4661 and len(dm)==len(dual['vertex_multipliers'])==780
    assert all(0<=i<82806 and a>0 for i,a in da.items())and all(0<=i<6558 and a>0 for i,a in dm.items())
    assert sum(dm.values(),F())==F(dual['normalizer'])and lower==F(dual['objective_numerator'])/sum(dm.values(),F())
    assert max(F(a)for v,a in primal['vertex_budgets'])==upper
    assert all(n>=0 and isinstance(n,int)for j,n in primal['rounded_variable_numerators'])
    assert F(primal['repair_delta'])>=0
    def ckey(c):return(c['vertex'],F(c['vertex_cost']),tuple(sorted((i,F(m))for i,m in c['positive_face_coefficients'])))
    exported={}
    for tag,doc in((STEMS[3],mixed),(repair['baseline'],baseline)):
        for c in doc['columns']:exported.setdefault(ckey(c),[]).append((tag,c['index']))
    assert len(exported)==1355 and sum(map(len,exported.values()))==1751
    assert len(prepared['new_columns'])==1355
    for offset,c in enumerate(prepared['new_columns']):
        assert c['variable']==430666+offset and sorted(c['sources'])==sorted(map(list,exported.pop(ckey(c))))
        den=c['positive_denominator'];assert den>0 and F(c['scaled_vertex_cost'],den)==F(c['vertex_cost'])
        assert all((F(x)*den).denominator==1 for i,x in c['positive_face_coefficients'])
        load=sum((da.get(i,F())*F(x)for i,x in c['positive_face_coefficients']),F())
        assert load<=F(c['vertex_cost'])*dm.get(c['vertex'],F())
    assert not exported
    # Scope-level construction checks and exact saved-sign composition, not
    # a duplicate conic construction or numerical search.
    assert continuous['status']=='PASS_EXACT_CONTINUOUS_ZERO_ARC_LOW_R_BRIDGE_GATE'
    assert len(continuous['records'])==96
    cc=continuous['counts'];assert cc==dict(R_branches=94,X_cases=1246,boxes=1288,leaves=1266,max_depth=4,excluded_R_branches=98,angle_gap_leaves=1)
    tc=Counter();maxdepth=0
    def interval(pair):
        lo,hi=map(F,pair);assert lo<=hi;return lo,hi
    def tree(node,depth=0):
        nonlocal maxdepth
        tc['nodes']+=1;maxdepth=max(maxdepth,depth)
        if'children'in node:
            assert len(node['children'])==2 and'split'in node and'mid'in node
            tc['split']+=1
            for child in node['children']:tree(child,depth+1)
        elif'assigned_gap'in node:
            tc['gap']+=1;g=node['assigned_gap'];assert g['kind']=='strict_consecutive_angle_gap'and g['labels']==[6,7]
            lo,hi=interval(g['product']);assert interval(g['lower_threshold'])[1]<lo<=hi<interval(g['upper_threshold'])[0]
        else:
            tc['leaves']+=1;proof=node['proof'];assert len(proof['branches'])==2
            reasons={b['reason']:b for b in proof['branches']}
            assert len(reasons)==2 and set(reasons)=={'strict_positive_native_quadratic','AX_strictly_dotted_on_entire_parameter_quadrant'}
            a,c,b=map(interval,reasons['strict_positive_native_quadratic']['witness']['coefficients'])
            assert a[0]>0 and b[0]>=0 and c[0]>=0
            z=reasons['AX_strictly_dotted_on_entire_parameter_quadrant']
            for field in('F_shifted_coefficients','t_derivative_shifted_coefficients'):assert all(interval(x)[1]<0 for x in z[field])
            for field in('s_derivative_leading','s_derivative_at_boundary'):assert interval(z[field])[1]<0
            assert interval(z['maximum_AX'])[1]<-2
    rb=xc=excluded=0
    for rec in continuous['records']:
        assert len(rec['branches'])==2 and {b['sign']for b in rec['branches']}=={-1,1}
        for branch in rec['branches']:
            if'native_obstruction'in branch:excluded+=1
            else:
                rb+=1
                for case in branch['X_cases']:xc+=1;tree(case['certificate'])
    assert(rb,xc,excluded)==(94,1246,98)and tc==dict(nodes=1288,split=21,leaves=1266,gap=1)and maxdepth==4
    assert scout['status']=='BOUNDED_NUMERICAL_TWO_ZERO_ARC_CONNECTOR_DOMAIN'
    assert scout['counts']==dict(Y_equations=13888,Y_branches=27776,Y_wrong_center=15693,R_low_seeds=94,X_seeds=4)
    assert not scout['accepted_prefixes']and not scout['unresolved']
    blockers=dict((tuple(k),n)for k,n in scout['native_blockers'])
    assert sum(blockers.values())==27776-15693==12083 and blockers['30','63']==1210
    topology=scout['topology'];assert(topology['vertices'],topology['edges'],topology['triangles'])==(11,25,15)
    assert(topology['x'],topology['y'],topology['z'])==(5,6,4)
    assert scout['deleted_walls']==continuous['deleted_walls']and scout['exact_old_deletion_cost']==continuous['exact_old_deletion_cost']
    return dict(scalar_profile_checks=1015,independent_profile_matches=1015,unrestricted_dimension_bound=116,
        compact_facet_ambient_bound=87,independent_port_matches=527,independent_overlap_matches=4862,
        independent_corner_matches=312,low_order_recipes=243,low_order_marked_joins=24,
        high_disjoint_joins=1478,multiple_edge_obstructions=120,reviewed_nonnegative_witnesses=905,
        extra_cross_shortcuts=21393,extra_cross_pd_obstructions=168,local_extension_variants=78,
        mixed_columns=1251,mixed_monotone_repair=str(ratio),mixed_vertex_increases=65,
        finite_model_added_columns=1355,finite_model_bracket=[str(lower),str(upper)],
        construction_continuous_leaves=1266,construction_angle_gap_leaves=1,
        construction_six_cap_numerical_branches=27776,construction_six_cap_blocked_branches=12083,
        scope='Independent exact scalar and saved-mass/column composition, matched source-bound primary/review classifications, and scoped construction sign/count checks. No optimizer, historical census, full geometry replay or new coverage claim.')

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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_DISTANCE10_AND_MIXED_CORNERS',
                scope='Archive integrity and exact composition of distance10 dimension116, mixed fixed-field actual-corner constraints and the enlarged finite model, with scoped continuous and numerical zero-arc cap results; not formalization or repeated historical verification.',
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
