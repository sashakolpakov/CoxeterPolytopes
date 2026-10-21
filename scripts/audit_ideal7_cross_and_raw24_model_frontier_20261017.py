#!/usr/bin/env python3
"""Incremental frozen archive with independent small checks and dimension91 composition."""
import argparse
import ast
from collections import Counter,deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import re
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/ideal7_cross_and_raw24_model_frontier_20261017.md'
PREVIOUS='results/ideal7_and_raw24_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_ideal7_and_raw24_frontier_20261017.py'
OUT=ROOT/'results/ideal7_cross_and_raw24_model_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
LOCAL='fresh_ideal7_finite_cross_payments_20261017'
COUPON='fresh_h17_half_coupon_transfer_20261017_review'
DIM='fresh_ideal7_cross_dimension91_20261017'
HALL='fresh_h17_contextual_cusp_payments_20261017'
MODEL='fresh_mixed_raw24_full_model_20261017'
SCOUT='fresh_old30_connector_flip_20261017'
GATE='fresh_old30_connector_native_gate_20261017'
STEMS=(LOCAL,COUPON,DIM,DIM+'_review',HALL,MODEL,MODEL+'_support',SCOUT,GATE)
REPORT_STEMS=STEMS
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
    assert len(d)==len(rows)and all(a>=0 for a in d.values())
    return d

def local_checks():
    d=read(LOCAL);review=read(COUPON);pairs=set(combinations(range(17),2))
    assert len(d['order17']['patterns'])==80
    assert len({r['index']for r in d['order17']['patterns']})==80
    assert d['order17']['dropped_case_indices']==[383,393,400,647,662,675,813,828,841]
    costs={}
    for kind in 'ABCD':
        g=[set()for _ in range(17)]
        end=12 if kind=='D'else 14 if kind=='B'else 16
        edges=[(i,i+1)for i in range(end)]
        if kind=='A':edges.append((0,16))
        if kind=='B':edges.extend([(0,15),(0,16)])
        if kind=='D':edges.extend([(0,13),(0,14),(12,15),(12,16)])
        for a,b in edges:g[a].add(b);g[b].add(a)
        dist={}
        for start in range(17):
            row={start:0};todo=deque([start])
            while todo:
                v=todo.popleft()
                for w in g[v]:
                    if w not in row:row[w]=row[v]+1;todo.append(w)
            assert len(row)==17
            for w,l in row.items():dist[start,w]=l
        if kind=='A':paid={p for p in pairs if dist[p]<=7}
        elif kind=='B':paid={p for p in pairs if dist[p]<=9}|{(9,15),(9,16)}
        elif kind=='C':paid={p for p in pairs if dist[p]<=9}
        else:paid={p for p in pairs if dist[p]<13 and p!=(0,12)}
        target=d['order17']['targets'][kind]
        assert paid==set(map(tuple,target['paid_pairs']))and len(paid)==target['cost']
        assert len(paid)<=F(15,2)*17;costs[kind]=len(paid)
    assert costs==dict(A=119,B=111,C=108,D=127)
    independent={tuple(r['black']):r for r in review['cases']};needed=[]
    for r in d['order17']['replacement_cases']:
        q=independent[tuple(r['black'])];qr={tuple(c['marks']):c for c in q['corners']}
        total=F();base=F()
        for c in r['corners']:
            a=qr[tuple(c['marks'])]
            assert F(c['base'])==F(a['base_payment'])and F(c['cross_coupon'])==F(a['coupon'])
            assert F(c['weight'])==F(c['base'])+F(c['cross_coupon'])
            if c['finite']:total+=F(c['weight']);base+=F(c['base'])
        assert total==F(r['finite_total'])>=1
        if base<1:needed.append(r['case'])
    assert needed==[383,400,813,841]
    # Independent endpoint enumeration: for every b,c,d, the expression
    # is affine in a, so only its two rank-constrained endpoints are needed.
    budget_checks=0
    for n in range(92,109):
        worst=F()
        for c in range(n//14+1):
            for e in range((n-14*c)//16+1):
                remaining=n-14*c-16*e
                for b in range(remaining//2+1):
                    for a in{0,remaining-2*b}:
                        excess=F(3*a*e,2)+3*b*c-9*a-17*b-35*c-24*e
                        assert excess<=0;worst=max(worst,excess);budget_checks+=1
        saved=next(r for r in d['finite_cost']['rows']if r['dimension']==n)
        assert worst==F(saved['max_excess_over_9n'])==0
    for r in d['affine_cost']['classical_checks']:
        m=r['order'];kind=r['kind']
        expected={'A':7*m,'B':7*m-23+min(m-10,17)+min(m-11,17),
            'C':7*m-28+min(m-8,18)+min(m-9,18),
            'D':7*m-18+min(m-12,16)+min(m-13,16)}[kind]
        assert r['cost']==expected<=F(15*m,2)
    for r in d['affine_cost']['all_rank_linear_branches']:
        slope=F(15,2)-r['slope'];lo=r['lower_order'];hi=r['upper_order']
        assert hi is not None or slope>=0
        at=hi if slope<0 else lo
        assert at==r['minimum_at']and slope*at-r['intercept']==F(r['minimum_slack'])>=0
    return dict(order17_intrinsic_paid_sets=costs,exceptional_cases=9,coupon_needed_cases=needed,
        independent_finite_budget_endpoint_checks=budget_checks,large_rank_formula_comparisons=444,
        scope='Independent graph distance/rule reconstruction, exact primary/review coupon agreement, finite-budget endpoints and all-rank linear-branch comparisons. The80-pattern complete embeddings and physical sealing proof retain their explicit producer/reviewer scope.')

def dimension_checks():
    d=read(DIM);r=read(DIM+'_review');p=read('fresh_distance8_ideal_dimension97_20261017')
    assert d['dimension_bound']==r['dimension_bound']==91 and p['dimension_bound']==97
    assert d['excluded_dimensions']==list(range(92,98))
    assert d['unrestricted_finite_volume']and not d['arithmetic_hypothesis']and d['cusp_restriction']is None
    assert d['convex_order_cap']==15 and d['finite_vertex_budget']=='9n'
    assert d['ideal_normalized_budget']=='sum_i min(15/2,r_i/2)'
    assert len(d['reduced_profiles'])==len(r['reconstruction']['profiles'])==477
    for a,b in zip(d['reduced_profiles'],r['reconstruction']['profiles']):
        assert a=={k:v for k,v in b.items()if k!='independent_cumulative_G'}
        assert F(a['margin'])>0
    assert d['dimension_minima']==r['reconstruction']['dimension_minima']
    assert min(F(a['margin'])for a in d['reduced_profiles'])==F(10121,4194304)
    assert F(d['boundary_at_n91']['margin'])==F(-428955,2097152)
    return dict(dimension_bound=91,independently_reconstructed_profiles=477,
        minimum_margin='10121/4194304',boundary_margin='-428955/2097152',
        scope='Composition with the completed97 tail and comparison with the separate no-project-arithmetic reconstruction; this archive does not duplicate that reconstruction or inherited geometry.')

def hall_checks():
    d=read(HALL);faces=d['faces'];assert len(faces)==3456
    indices={r['face_index']:r for r in faces};assert len(indices)==3456
    for face in faces:
        fixed=sum(c['fixed_finite_payment']for c in face['corners']if c['cusp']is None)
        assert fixed==face['fixed_finite_payment']
        assert face['ideal_deficit']==max(0,5-face['polygon']-fixed)
        assert len(face['eligible_cusps'])==len(set(face['eligible_cusps']))
    cuts={}
    for r in d['hall_constraints']:
        U=set(r['cusps'])
        demand=sum(max(0,f['ideal_deficit']-len(set(f['eligible_cusps'])-U))for f in faces)
        assert demand==r['demand'];cuts[','.join(sorted(U))]=demand
    assert cuts=={'':0,'D16':109,'E8xE8':612,'D16,E8xE8':723}
    for c in d['certificates']:
        pay={r['face_index']:r['paid_cusps']for r in c['nonzero_payments']};counts=Counter()
        assert len(pay)==len(c['nonzero_payments'])
        for f in faces:
            paid=pay.get(f['face_index'],[])
            assert len(paid)==len(set(paid))and set(paid)<=set(f['eligible_cusps'])
            assert len(paid)>=f['ideal_deficit'];counts.update(paid)
        assert dict(counts)==c['cusp_costs']and sum(counts.values())==723
    assert max(F(109,17),F(612,162),F(723,179))==F(d['best_uniform_normalized_cost'])
    return dict(face_deficit_reconstructions=3456,hall_cut_demands=cuts,exact_allocations=3,
        scope='Exact saved-face deficit and allocation reconstruction for the fixed known chamber; no universal Hall-demand or dimension theorem.')

def arithmetic_checks():
    d=read(MODEL);s=read(MODEL+'_support');counts=d['prepared']['counts']
    assert counts['variables']==435793 and counts['distinct_added_columns']==2040
    assert counts['exported_input_columns']==2319 and counts['exact_duplicate_exports']==279
    assert counts['faces']==82806 and counts['vertex_forests']==6558
    assert counts['zero_budget_signed_variables']==96176 and counts['vertex_cost_signed_variables']==52
    dual=next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    primal=next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    alpha=masses(dual['face_multipliers']);mu=masses(dual['vertex_multipliers'])
    assert sum(mu.values(),F())==F(dual['normalizer'])
    assert F(dual['objective_numerator'])/F(dual['normalizer'])==F(dual['bound'])
    assert F(dual['bound'])==F(5755387047017860,462000155393947)>12
    assert F(primal['maximum'])==F(3114386878673,250000000000)>F(dual['bound'])
    assert s['face_multipliers']==dual['face_multipliers']and s['vertex_multipliers']==dual['vertex_multipliers']
    assert d['prepared']['model_sha256']==s['model_sha256']
    assert len(alpha)==4647 and len(mu)==807 and s['optimizer_invocations']==0 and d['solver_invocations']==1
    qs={r['face']:r for r in s['quad_rows']};ss=set(alpha)-set(qs)
    assert len(qs)==1900 and len(ss)==2747
    restrictions=[]
    for r in s['exact_restrictions']:
        a=r['certificate'];ra=masses(a['face_multipliers']);rm=masses(a['vertex_multipliers'])
        assert ss<=set(ra)and all(alpha[i]==v for i,v in ra.items())
        assert all(v<=mu.get(i,F())for i,v in rm.items())
        kept={i for i,q in qs.items()if q['in_mixed_export']}if r['name']=='already_mixed_exported_Q'else{i for i,q in qs.items()if q['provenance_kind']!='merged_unexpanded'}
        assert set(ra)==ss|kept and len(kept)==r['kept_Q']
        # Face requirements include aggregate simplex/edge multiplicities;
        # their reconstruction belongs to the full model verifier.
        numerator=F(a['objective_numerator'])
        assert sum(rm.values(),F())==F(a['normalizer'])
        assert numerator/sum(rm.values(),F())==F(a['bound'])
        restrictions.append(a['bound'])
    assert F(restrictions[0])<12<F(restrictions[1])
    loads=0
    for group in s['appended_export_loads']:
        assert len(group['loads'])==group['columns']
        for c in group['loads']:assert F(c['capacity'])-F(c['load'])==F(c['slack'])>=0;loads+=1
    assert loads==2319
    cohorts={r['raw_family']:r for r in s['uncovered_exact_provenance_cohorts']}
    assert [(cohorts[i]['classes'],cohorts[i]['actual_corner_occurrences'])for i in(178,296,297)]==[(18,144),(16,128),(13,104)]
    return dict(model_optimum_lower=dual['bound'],model_optimum_upper=primal['maximum'],
        saved_export_loads=loads,restricted_bounds=restrictions,next_cohort_classes=47,
        scope='Saved exact certificate composition and support restrictions. Full435793-column feasibility is checked by producer generation/fresh replay; no independent matrix rebuild or optimization here.')

def construction_checks():
    s=read(SCOUT);d=read(GATE)
    assert not s['complete']and not s['unresolved']
    assert s['counts']['direct_YB_tests']==1985 and s['counts'].get('direct_valid_finite_YB',0)==0
    for key in('direct_subdivision_topology','flipped_topology'):
        t=s[key];edges=Counter(tuple(sorted(e))for tri in t['triangles']for e in combinations(tri,2))
        boundary={tuple(sorted(e))for e in zip(t['boundary'],t['boundary'][1:]+t['boundary'][:1])}
        assert {e for e,c in edges.items()if c==1}==boundary and all(c in(1,2)for c in edges.values())
        vertices={v for tri in t['triangles']for v in tri};new=vertices-set(t['boundary'])
        edge_kinds=Counter(sum(v in new for v in e)for e in edges)
        assert edge_kinds=={0:5,1:12,2:14}
        assert len(edges)==31 and len(t['triangles'])==19 and len(vertices)==t['vertices']==13 and 13-31+19==1
        assert 2*(edge_kinds[1]+edge_kinds[2])-3*len(t['triangles'])==t['conditional_delta_ends']==-5
    assert s['predecessor_prose_correction']['actual_profile_counts']==d['predecessor_profile_counts']=={'0,0':14,'0,1':2,'0,4':2,'1,0':2}
    count=0;high=0;reasons=Counter()
    for r in d['records']:
        assert len(r['exact_prefix_recipes'])==5
        for case in r['complete_Z_cases']:
            high+=bool(case['high_positions']);proof=case['proof']
            assert 'children'not in proof
            c=proof['certificate'];assert c['reason']=='both_native_branches_excluded'
            assert {b['sign']for b in c['branches']}=={-1,1}
            for b in c['branches']:
                q=b['rejection'];lo,hi=map(F,q['product']);assert lo<=hi;reasons[q['kind']]+=1
                if q['kind']=='strict_obtuse':assert lo>0
                else:
                    assert q['kind']=='strict_consecutive_angle_gap'and q['labels'][1]==q['labels'][0]+1
                    assert F(q['lower_threshold'][1])<lo<=hi<F(q['upper_threshold'][0])
                count+=1
    assert len(d['records'])==20 and high==2 and count==594
    assert reasons==dict(strict_obtuse=390,strict_consecutive_angle_gap=204)
    return dict(exact_disk_triangles=19,literal_prefixes=20,exact_branch_rejections=594,high_label_families=2,
        scope='Independent abstract edge counts and saved rational rejection inequalities. No independent physical-interval reconstruction or finite-volume promotion; exact gate only covers the20 stated algebraic prefixes.')

def small_checks():
    return dict(local=local_checks(),dimension=dimension_checks(),contextual=hall_checks(),
        arithmetic=arithmetic_checks(),construction=construction_checks(),pending_inputs=list(PENDING_INPUTS),
        scope='Unrestricted dimension91 proof composition and scoped independent small certificate checks; arithmetic-compact21 and construction frontier unchanged.')

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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_IDEAL7_CROSS_DIMENSION91_AND_RAW24_MODEL',
                scope='Archive integrity, independent small local-payment checks, completed unrestricted dimension91 composition, exact arithmetic-model certificates and scoped construction diagnostics. No new finite-volume example.',
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
