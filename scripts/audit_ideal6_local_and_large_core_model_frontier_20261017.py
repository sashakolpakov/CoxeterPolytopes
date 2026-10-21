#!/usr/bin/env python3
"""Frozen incremental archive for ideal6 local lemmas and the enlarged arithmetic model."""
import argparse,ast,json,re
from collections import Counter,deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/ideal6_local_and_large_core_model_frontier_20261017.md'
PREVIOUS='results/order16_and_large_core_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_order16_and_large_core_frontier_20261017.py'
OUT=ROOT/'results/ideal6_local_and_large_core_model_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
TRI='fresh_ideal6_triangle_payments_20261017'
EMBED='fresh_ideal6_small_factor_frontier_20261017'
LOCAL='fresh_ideal_distance6_mixed_replacement_20261017'
MODEL='fresh_raw178_296_297_full_model_20261017'
SUPPORT=MODEL+'_support'
FIRST='fresh_old8_triangle_flip_20261017'
SECOND='fresh_old8_triangle_double_flip_20261017'
STEMS=(TRI,EMBED,LOCAL,MODEL,SUPPORT,FIRST,SECOND)
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

def small_checks():
    tri=read(TRI);embedding=read(EMBED);local=read(LOCAL)
    ds={d['index']:d for d in read('fresh_high_rank_quasilanner_census_20261017')['diagrams']}
    checked=0;exceptions=[]
    for row in tri['rows']:
        d=ds[row['diagram']];n=d['order'];g=[[0]*n for _ in range(n)]
        for i,j,m in d['edges']:g[i][j]=g[j][i]=m
        paid=0
        for corner in row['corners']:
            a,b=corner['pair'];ids=set(range(n))-(set(row['marks'])-{a,b})
            dist=distances(g,a,ids)
            p=corner['path']
            if b not in dist:assert p is None and not corner['elementary_payment'];continue
            assert p[0]==a and p[-1]==b and set(p)<=ids
            assert all(g[i][j] for i,j in zip(p,p[1:]))
            assert len(p)-1==corner['length']==dist[b]
            features={i for i in dist if sum(bool(g[i][j]) for j in dist)>=3
                      or any(g[i][j]>=4 for j in dist)}
            near=any(min(dist[f],distances(g,b,ids)[f])<=6 for f in features)
            expected=dist[b]<=6 or(dist[b]<=9 and near)
            assert expected==corner['elementary_payment'];paid+=expected
        assert paid==row['elementary_payment_count']
        if paid<2:exceptions.append([row['diagram'],row['marks']])
        checked+=1
    assert checked==920 and exceptions==[[25,[1,2,9]]]
    exc=tri['exceptional_resolution'];d=ds[25]
    g=[[F(2 if i==j else 0)for j in range(10)]for i in range(10)]
    for i,j,m in d['edges']:assert m==3;g[i][j]=g[j][i]=F(-1)
    z=list(map(F,exc['affine_kernel']));assert all(v>0 for v in z)
    assert all(sum(g[i][j]*z[j]for j in range(9))==0 for i in range(9))
    pd(g,list(range(8)));pd(g,[i for i in range(10)if i!=2])
    patterns={r['index']:r for r in embedding['patterns']}
    old=read('fresh_selective_distance9_dimension108_20261017')
    assert set(patterns)=={i for i,r in enumerate(old['isolation']['cases'])if r['component_order_bound']>=15}
    assert len(patterns)==437
    maps=0;costs={};drop15=[]
    for row in embedding['records']:
        m=row['order'];g=row['target'];dd=[distances(g,i)for i in range(m)]
        fs={i for i in range(m)if sum(bool(v)for v in g[i])>=3 or any(v>=4 for v in g[i])}
        near={i for i in range(m)if any(dd[i][j]<=6 for j in fs)}
        basic={p for p in combinations(range(m),2)if dd[p[0]][p[1]]<=6 or(dd[p[0]][p[1]]<=9 and set(p)&near)}
        assert basic==set(map(tuple,row['basic_pairs']));extra=set();affected=[]
        proposed=set(map(tuple,row['proposed_unproved_order15_drop']))
        for r in row['case_embeddings']:
            p=patterns[r['case']];assert p['order_bound']>=m
            images=set()
            for mp in r['embeddings']:
                assert len(mp)==len(p['nodes'])and len(set(mp))==len(mp)and all(0<=i<m for i in mp)
                assert all(p['labels'][i][j]==g[mp[i]][mp[j]]for i in range(len(mp))for j in range(len(mp)))
                images.add(tuple(sorted(mp[i]for i in p['marks'])));maps+=1
            assert images==set(map(tuple,r['marked_images']))
            assert images-basic==set(map(tuple,r['extra_images']));extra|=images-basic
            if images&proposed:affected.append(r['case'])
        assert extra==set(map(tuple,row['extra_pairs']))
        dropped=set(map(tuple,row['inherited_dropped_pairs']))
        if dropped:
            saved=read('fresh_order16_finite_cross_payments_20261017')['exceptional_D15']['dropped_pairs'] if m==16 else read('fresh_ideal7_finite_cross_payments_20261017')['order17']['targets']['D']['dropped_pairs']
            assert dropped==set(map(tuple,saved))
        paid=(basic|extra)-dropped
        assert paid==set(map(tuple,row['paid_pairs']))and len(paid)==row['cost']
        assert F(len(paid),m)==F(row['normalized_cost'])
        costs[str(m)+row['kind']]=len(paid)
        assert affected==row['proposed_drop_affected_cases']
        if proposed:assert len(proposed)==9 and len(affected)==63;drop15=affected
    assert len(costs)==12 and costs['15D']==105
    for r in local['basic_affine_cost']['all_rank_branches']:
        slope=F(r['normalized_bound'])-r['slope']
        at=r['upper_order']if slope<0 else r['lower_order']
        assert at is not None and at==r['minimum_slack_at']
        assert slope*at-r['intercept']==F(r['minimum_slack'])>=0
    required={m:4-1/(F(1,7)+F(1,m-7))for m in(14,15,16)}
    assert required=={int(m):F(v)for m,v in local['small_cycle_checks']['required_unscaled_spectator_difference'].items()}
    core=local['cycle14_equality'];g=[[F(v)for v in row]for row in core['gram']]
    z=list(map(F,core['radical']));h=list(map(F,core['strict_center']))
    assert all(sum(g[i][j]*z[j]for j in range(17))==0 for i in range(17))
    assert sum(h[i]*g[i][j]*h[j]for i in range(17)for j in range(17))==-204
    assert all(sum(g[i][j]*h[j]for j in range(17))<0 for i in range(17))
    for r in core['corners']:
        if r['kind']=='finite':pd(g,sorted(set(core['white'])|set(r['marks'])))
    assert sum(F(r['coupon'])for r in core['corners'])==1
    model=read(MODEL);support=read(SUPPORT)
    primal=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha=mass(dual['face_multipliers']);mu=mass(dual['vertex_multipliers'])
    assert len(alpha)==4410 and len(mu)==984
    assert F(dual['objective_numerator'])/sum(mu.values(),F())==F(dual['bound'])
    assert sum(mu.values(),F())==F(dual['normalizer'])
    assert F(12)<F(dual['bound'])<=F(primal['maximum'])==max(mass(primal['vertex_budgets']).values())
    assert model['solver_invocations']==1 and support['optimizer_invocations']==0
    assert support['face_multipliers']==dual['face_multipliers']and support['vertex_multipliers']==dual['vertex_multipliers']
    columns=read('fresh_raw178_296_297_actual_corners_20261017')['columns']
    loads=support['appended_export_loads'][0]['loads'];assert len(columns)==len(loads)==469
    for c,r in zip(columns,loads):
        load=sum(F(a)*alpha.get(i,F())for i,a in c['positive_face_coefficients'])
        capacity=F(c['vertex_cost'])*mu.get(c['vertex'],F())
        assert r['index']==c['index']and r['vertex']==c['vertex']
        assert load==F(r['load'])<=capacity==F(r['capacity'])
        assert capacity-load==F(r['slack'])
    for r in support['exact_restrictions']:
        c=r['certificate'];aa=mass(c['face_multipliers']);mm=mass(c['vertex_multipliers'])
        assert all(alpha[i]==v for i,v in aa.items())and all(v<=mu.get(i,F())for i,v in mm.items())
        assert F(c['objective_numerator'])/sum(mm.values(),F())==F(c['bound'])
    construction=[]
    for stem in(FIRST,SECOND):
        d=read(stem);assert not d['complete']and not d['unresolved']
        if stem==FIRST:
            assert not d['prefixes']and not d['YZ_bank']and d['counts']['Z_branches']==14560
            assert d['Z_obstruction_types']==dict(obtuse=13712,acute_angle_gap=848)
        else:
            assert len(d['prefixes'])==6 and d['counts']['Z_branches']==23562
            assert d['W_stops']==dict(WZ_dotted=548,WZ_None=10)
            for r in d['prefixes']:
                bd,ic=r['boundary'],r['PD_incidence']
                assert not bd['warnings']and bd['counts']==dict(finite=26,exterior=8)
                assert bd['open_ends']==ic['open_ends']==12 and ic['finite']==26
                rays=[v for v in bd['rays']if 71 in v['active']]
                assert len(rays)==3 and all(v['kind']=='finite'for v in rays)
                assert sorted(sorted(v['active'])for v in rays)==r['Y_actual_finite_corners']
        construction.append(dict(package=stem,counts=d['counts'],scope='Saved-domain and empty-completion consistency only; geometry uses generation/replay and source review.'))
    return dict(triangle_rows=checked,small_factor_patterns=437,validated_labelled_maps=maps,
                small_factor_costs=costs,unresolved_order15_contexts=len(drop15),
                local_affine_branches=12,spectator_required_values={str(k):str(v)for k,v in required.items()},
                arithmetic_lower=dual['bound'],arithmetic_upper=primal['maximum'],recounted_new_capacity_loads=469,
                construction=construction,pending_inputs=list(PENDING_INPUTS),
                scope='Independent saved path/payment, labelled-map, exact local kernel/PD/branch and469-capacity-load checks. Does not repeat minimal-core classification, universal supporting-white geometry, full matrix recovery or numerical construction enclosures. No dimension ceiling below86 is established.')

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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_IDEAL6_LOCAL_AND_LARGE_CORE_MODEL',
                scope='Archive integrity and scoped independent exact checks of new local lemmas, small-factor embeddings, finite-model certificates and changed construction incidence. Completed dimension86 unchanged; no new finite-volume example.',
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
