#!/usr/bin/env python3
"""Frozen incremental archive and small exact composition checks for dimension86."""
import argparse,ast,json,re
from collections import Counter,deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/order16_and_large_core_frontier_20261017.md'
PREVIOUS='results/ideal7_cross_and_raw24_model_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_ideal7_cross_and_raw24_model_frontier_20261017.py'
OUT=ROOT/'results/order16_and_large_core_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
EMBED='fresh_order16_preferred_embeddings_20261017'
CYCLE='fresh_cycle16_half_coupon_transfer_20261017'
LOCAL='fresh_order16_finite_cross_payments_20261017'
DIM='fresh_order16_piecewise_dimension86_20261017'
ARITH='fresh_raw178_296_297_actual_corners_20261017'
DEFORM='fresh_dotted_a53_disk_20261017'
SCOUT='fresh_c10_endpoint_disk_20261017'
STEMS=(EMBED,CYCLE,LOCAL,DIM,DIM+'_review',ARITH,DEFORM,SCOUT)
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
def mass(rows):
    out={i:F(a)for i,a in rows};assert len(out)==len(rows)and all(a>=0 for a in out.values())
    return out

def distances(g,start,ids=None):
    ids=set(range(len(g)))if ids is None else set(ids)
    ds={start:0};q=deque([start])
    while q:
        i=q.popleft()
        for j in ids-set(ds):
            if i!=j and g[i][j]:ds[j]=ds[i]+1;q.append(j)
    return ds

def inverse(a):
    n=len(a);work=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(a)]
    for i in range(n):
        j=next(j for j in range(i,n)if work[j][i]);work[i],work[j]=work[j],work[i]
        p=work[i][i];work[i]=[v/p for v in work[i]]
        for j in range(n):
            if j==i:continue
            p=work[j][i];work[j]=[a-p*b for a,b in zip(work[j],work[i])]
    return [row[n:]for row in work]

def pd(g,ids):
    a=[[F(g[i][j])for j in ids]for i in ids]
    while a:
        p=a[0][0];assert p>0
        a=[[a[i][j]-a[i][0]*a[0][j]/p for j in range(1,len(a))]for i in range(1,len(a))]

def embedding_checks():
    d=read(EMBED);old=read('fresh_selective_distance9_dimension108_20261017');patterns={r['index']:r for r in d['patterns']}
    assert len(patterns)==260 and set(patterns)=={i for i,r in enumerate(old['isolation']['cases'])if r['component_order_bound']>=16}
    for i,p in patterns.items():
        r=old['isolation']['cases'][i]
        assert p['nodes']==r['marked_component']and [p['nodes'][j]for j in p['marks']]==r['preferred_corner']
        assert p['order_bound']==r['component_order_bound']
    totals={};maps=0
    for kind,t in d['targets'].items():
        g=t['target'];assert len(g)==16
        ds=[distances(g,i)for i in range(16)]
        features={i for i in range(16)if sum(bool(v)for v in g[i])>=3 or any(v>=4 for v in g[i])}
        near={i for i in range(16)if any(ds[i][j]<=7 for j in features)}
        basic={p for p in combinations(range(16),2)if ds[p[0]][p[1]]<=7 or(ds[p[0]][p[1]]in(8,9)and set(p)&near)}
        assert basic==set(map(tuple,t['basic_pairs']));extra=set()
        for row in t['case_embeddings']:
            p=patterns[row['case']];images=set()
            for mp in row['embeddings']:
                assert len(mp)==len(p['nodes'])and len(set(mp))==len(mp)and all(0<=v<16 for v in mp)
                assert all(p['labels'][i][j]==g[mp[i]][mp[j]]for i in range(len(mp))for j in range(len(mp)))
                images.add(tuple(sorted(mp[i]for i in p['marks'])));maps+=1
            assert images==set(map(tuple,row['marked_images']))and images-basic==set(map(tuple,row['extra_images']))
            extra|=images-basic
        assert extra==set(map(tuple,t['extra_pairs']))and basic|extra==set(map(tuple,t['paid_pairs']))
        assert len(basic|extra)==t['cost'];totals[kind]=t['cost']
    assert totals==dict(A=112,B=107,C=103,D=120)
    return dict(preferred_domains=260,independently_validated_saved_maps=maps,costs_before_replacement=totals,
        scope='Independent domain-filter, emitted-map validity and union checks. Completeness of injective backtracking is its written proof and producer replay, not a duplicate search here.')

def local_checks():
    c=read(CYCLE);d=read(LOCAL);cycle=c['equality_core'];g=[[F(v)for v in row]for row in cycle['gram']]
    assert len(g)==18 and len(cycle['ideal_short7_pairs'])==112 and len(cycle['dropped_antipodal_pairs'])==8
    assert all(sum(g[i][j]*cycle['radical'][j]for j in range(18))==0 for i in range(18))
    for r in cycle['corners']:
        if r['kind']!='finite':continue
        corner_ids=sorted(set(cycle['white'])|set(r['marks']));pd(g,corner_ids)
        center_parts=0;component_sizes=[]
        for p in r['marked_components']:
            ids=p['nodes'];mark=p['mark'];ds=distances(g,mark,ids)
            assert len(ds)==len(ids)and set(ids)==set(distances(g,mark,corner_ids))
            component_sizes.append(len(ids))
            if p['unique_middle']:
                leaves=[i for i in ids if sum(i!=j and bool(g[i][j])for j in ids)==1]
                assert len(ids)==15 and len(leaves)==2 and all(ds[i]==7 for i in leaves)
                center_parts+=1
        expected_coupon=F(1,2)if sorted(component_sizes)==[1,15]and center_parts==1 else F()
        assert F(r['coupon'])==expected_coupon
    assert sum(F(r['coupon'])for r in cycle['corners'])==1
    values=0
    for r in c['residual_profiles']['rows']:
        p=[(F(2),F(),F(1))if kind=='isolated'else(F(),F(1),F(1,2))for kind in(r['c'],r['d'])]
        actual=[]
        for row in r['rectangle_values']:
            u,v=[flag*z[2]for flag,z in zip(row['upper_reciprocals'],p)]
            x=1-(p[0][1]*u+p[1][1]*v)/2-(p[0][0]*u*u+p[1][0]*v*v)/4
            assert x==F(row['value']);actual.append(x);values+=1
        assert min(actual)==F(r['minimum'])>=0
    assert values==16
    e=d['exceptional_D15'];g=[[F(v)for v in row]for row in e['gram']];iv=inverse(g)
    assert [[iv[i][j]for j in(3,13)]for i in(3,13)]==[[4,-4],[-4,4]]
    z=list(map(F,e['affine_positive_kernel']));ids=e['affine_nodes']
    assert all(v>0 for v in z)and all(sum(g[i][j]*z[k]for k,j in enumerate(ids))==0 for i in ids)
    corners=0;needed=[]
    for row in e['cases']:
        white=set(range(18))-set(row['black']);base=F();total=F()
        for p in row['corners']:
            if not p['finite']:assert F(p['weight'])==0;continue
            a,b=p['marks'];ids=sorted(white|{a,b});pd(g,ids);corners+=1
            ca=set(distances(g,a,ids));cb=set(distances(g,b,ids));unit=F(b in ca);coupon=F()
            if ca!=cb:
                for small,large,x,y in((ca,cb,a,b),(cb,ca,b,a)):
                    if(len(small),len(large))not in((1,15),(2,13)):continue
                    adj={i:{j for j in large if i!=j and g[i][j]}for i in large}
                    forks=[i for i in large if len(adj[i])==3]
                    if len(forks)==1:
                        f=forks[0];short={i for i in adj[f]if len(adj[i])==1}
                        if len(short)==2 and y in short|{f}:coupon=F(1,2)
            assert unit==F(p['base'])and coupon==F(p['cross_coupon'])and unit+coupon==F(p['weight'])
            base+=unit;total+=unit+coupon
        assert base==F(row['base_total'])and total==F(row['total'])>=1
        if base<1:needed.append(row['case'])
    assert corners==27 and needed==[384,401,814,842]
    assert d['order16_costs']==dict(A=112,B=107,C=103,D=111)
    assert len(e['remaining_D_pairs'])==111
    for r in d['large_order_linear_branches']:
        slope=F(217,29)-r['slope'];at=r['upper_order']if slope<0 else r['lower_order']
        assert at is not None and at==r['minimum_at']and slope*at-r['intercept']==F(r['minimum_slack'])>=0
    assert {r['dimension']for r in d['finite_cost']['rows']}==set(range(87,109))
    assert all(F(r['maximum_excess'])==0 for r in d['finite_cost']['rows'])
    return dict(independent_cycle_rectangle_values=16,independent_finite_D15_corners=27,
        independently_rebuilt_inverse_contact_block=[[4,-4],[-4,4]],coupon_needed_cases=needed,
        finite_budget_saved_endpoints=sum(r['endpoint_checks']for r in d['finite_cost']['rows']),
        scope='Independent small rational kernels, inverse block, finite corner PD/types/payments, cycle rectangles and affine linear bounds. Universal support arguments, complete profile classification and finite-budget case split retain their written review scope.')

def dimension_checks():
    d=read(DIM);r=read(DIM+'_review');assert d['dimension_bound']==r['dimension_bound']==86
    assert read('fresh_ideal7_cross_dimension91_20261017')['dimension_bound']==d['inherited_dimension_bound']==91
    assert d['excluded_dimensions']==list(range(87,92))and d['unrestricted_finite_volume']and not d['arithmetic_hypothesis']and d['cusp_restriction']is None
    assert len(d['reduced_profiles'])==len(r['reconstruction']['profiles'])==3355
    for a,b in zip(d['reduced_profiles'],r['reconstruction']['profiles']):
        assert a=={k:v for k,v in b.items()if k!='independent_cumulative_G'}and F(a['margin'])>0
    assert d['dimension_minima']==r['reconstruction']['dimension_minima']
    assert min(F(a['margin'])for a in d['reduced_profiles'])==F(146417,911250)
    assert F(d['boundary_at_n86']['margin'])==F(-112409,1518750)
    return dict(dimension_bound=86,independently_reconstructed_profiles=3355,
        minimum_margin='146417/911250',boundary_margin='-112409/1518750',
        scope='Comparison with the separate no-project-arithmetic reconstruction and completed91 tail; no third reconstruction or duplicate inherited geometry here.')

def arithmetic_checks():
    d=read(ARITH);base=read('fresh_mixed_raw24_full_model_20261017')
    cert=next(c for c in base['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert d['face_multipliers']==cert['face_multipliers']and d['initial_vertex_multipliers']==cert['vertex_multipliers']
    alpha=mass(d['face_multipliers']);old=mass(d['initial_vertex_multipliers']);new=mass(d['vertex_multipliers']);rebuilt=dict(old)
    columns={r['index']:r for r in d['columns']};assert len(columns)==469
    paircerts=colors=cliques=0;visited=set()
    for forest in d['forest_records']:
        stars={r['index']:r for r in forest['literal_stars']};orbits={r['index']:r for r in forest['corner_orbits']}
        for orbit in orbits.values():
            counts=Counter(o['face']for o in orbit['origins']);saved=dict(orbit['positive_face_coefficients'])
            assert len(saved)==len(orbit['positive_face_coefficients'])and counts==saved
            assert all(o['vertex']==forest['vertex']for o in orbit['origins'])
        exclusions={tuple(r[:2])for r in forest['incompatibility_certificates']}
        assert len(exclusions)==len(forest['incompatibility_certificates']);paircerts+=len(exclusions)
        assert forest['tested_literal_pairs']==len(stars)*(len(stars)-1)//2
        for j in forest['columns']:
            assert j not in visited;visited.add(j)
            c=columns[j];assert c['vertex']==forest['vertex']
            coeff=Counter()
            if c['kind']=='capacity_one_orbit_clique':
                selected=c['orbit_clique'];nodes=[i for q in selected for i in orbits[q]['literal_stars']]
                assert all(tuple(sorted(p))in exclusions for p in combinations(nodes,2));assert c['vertex_cost']==1;cliques+=1
            else:
                assert c['kind']=='proper_compatibility_coloring';selected=c['orbit_indices']
                nodes=[i for q in selected for i in orbits[q]['literal_stars']];parts=c['literal_coloring']
                assert sorted(i for part in parts for i in part)==sorted(nodes)and len(parts)==c['vertex_cost']
                assert all(tuple(sorted(p))in exclusions for part in parts for p in combinations(part,2));colors+=1
            assert len(nodes)==len(set(nodes))
            for q in selected:
                for i,a in orbits[q]['positive_face_coefficients']:coeff[i]+=a
            assert sorted(coeff.items())==[tuple(r)for r in c['positive_face_coefficients']]
            load=sum(alpha[i]*a for i,a in coeff.items());cost=F(c['vertex_cost']);v=c['vertex']
            assert load==F(c['load'])and old[v]==F(c['baseline_vertex_mass'])and c['violated']==(load>cost*old[v])
            rebuilt[v]=max(rebuilt.get(v,F()),load/cost)
    assert visited==set(columns)and rebuilt==new and sum(new[v]>old.get(v,F())for v in new)==65
    assert paircerts==760379 and(cliques,colors)==(268,201)
    repair=d['fixed_witness_repair'];num=F(cert['objective_numerator'])
    assert num==F(repair['objective_numerator'])and sum(new.values(),F())==F(repair['normalizer'])
    assert num/sum(new.values(),F())==F(repair['bound'])==F(2877693523508930,247081190427751)<12
    assert repair['optimizer_invocations']==0
    assert len(d['selected_faces'])==47 and sum(len(r['corners'])for r in d['selected_faces'])==376
    assert {r['raw_core_order']for r in d['selected_faces']}=={9,10}
    assert all(r['raw_core_order']+r['spectator_rank']==15 for r in d['selected_faces'])
    return dict(independent_capacity_cliques=cliques,independent_capacity_colorings=colors,
        retained_pair_exclusion_inputs=paircerts,vertex_allowance_increases=65,bound=repair['bound'],
        scope='Independent capacity/multiplicity and exact fixed-witness composition conditional on the bound pair exclusions. Native pair projections/field tests and complete orientation census use producer generation/replay and scoped source review, not independent re-enclosure here.')

def construction_checks():
    d=read(DEFORM);s=read(SCOUT)
    assert len(d['parameters'])==19 and not d['accepted']and not d['unresolved']
    assert d['counts']['A_native']==19 and d['counts']['C_branches']==760 and d['counts']['C_native']==20
    assert len(s['AC_bank'])==27 and len(s['prefixes'])==254 and not s['Z_bank']and not s['unresolved']
    assert s['counts']['Z_equations']==2865 and s['counts']['Z_branches']==5730
    assert all(r['A53']<-2 for r in s['AC_bank'])
    pairs=Counter()
    for r in s['prefixes']:
        bd,ic=r['boundary'],r['PD_incidence'];assert not bd['warnings']
        assert bd['counts']==dict(finite=24,exterior=8)and bd['open_ends']==ic['open_ends']==12 and ic['finite']==24
        assert len(r['native_pairs'])==10;pairs[r['X_profile'],r['Y_profile']]+=1
    assert pairs=={(0,0):168,(0,1):42,(0,4):42,(0,2):1,(0,5):1}
    assert s['deleted_walls']==d['deleted_walls']==[12,47,52,64,65]
    return dict(parameter_samples=19,new_endpoint_pairs=27,numerical_prefixes=254,blocked_Z_branches=5730,
        scope='Saved numerical acceptance/count/boundary consistency only; no independent numerical geometry, all-label domain or exact promotion. Control63/199/10 unchanged.')

def small_checks():
    return dict(embeddings=embedding_checks(),local=local_checks(),dimension=dimension_checks(),
        arithmetic=arithmetic_checks(),construction=construction_checks(),pending_inputs=list(PENDING_INPUTS),
        scope='Completed unrestricted dimension86 composition, scoped independent small exact checks, new fixed arithmetic witness repair and changed-prefix construction diagnostics.')

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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_ORDER16_DIMENSION86_AND_LARGE_CORE',
                scope='Archive integrity, independent small order16 payment and capacity checks, completed unrestricted dimension86 composition and changed-prefix construction diagnostics. No new finite-volume example.',
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
