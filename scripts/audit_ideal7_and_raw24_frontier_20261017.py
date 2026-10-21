#!/usr/bin/env python3
"""Incremental archive and exact local/certificate composition.

No new ambient dimension bound, global optimization or geometric census.
"""
import argparse
import ast
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import re
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/ideal7_and_raw24_frontier_20261017.md'
PREVIOUS='results/distance8_feature_and_mixed_capacity_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_distance8_feature_and_mixed_capacity_frontier_20261017.py'
OUT=ROOT/'results/ideal7_and_raw24_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
IDEAL='fresh_ideal_distance7_mixed_replacement_20261017'
BOTTLENECK='fresh_h17_d16_payment_bottleneck_20261017'
LOCAL='fresh_mixed_raw39_full_model_20261017_local_cones'
JOINT='fresh_mixed_raw39_full_model_20261017_joint_ports'
RAW24='fresh_raw24_actual_corner_compatibility_20261017'
MIXED='fresh_mixed_raw24_actual_corners_20261017'
IDENTITY='fresh_old9_terminal_identity_20261017'
SCOUT='fresh_interior_connector_disk_20261017'
STEMS=(IDEAL,BOTTLENECK,LOCAL,JOINT,RAW24,MIXED,IDENTITY,SCOUT)
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
def load(alpha,rows):return sum((alpha[i]*F(a)for i,a in rows),F())
def q(x):return tuple(map(F,x))
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def div(x,y):
    n=y[0]**2-2*y[1]**2
    assert n
    z=mul(x,(y[0],-y[1]));return z[0]/n,z[1]/n
def rational(x):
    assert x[1]==0
    return x[0]

def ideal_checks(d):
    c=d['cycle_mixed_replacement'];ps=c['profiles'];expected=[];values=[]
    assert len(ps)==8
    for p,r in product(ps,repeat=2):
        lo=[8,8]
        for s in(p,r):
            if s['name']not in('isolated','both_ends3'):
                for j,k in enumerate(s['position']):
                    if k:lo[j]=max(lo[j],9+k)
        boxes=[lo]if sum(lo)>=17 else[[8,9],[9,8]]
        for a,b in boxes:
            for u,v in product((0,1),repeat=2):
                U=q(p['umax'])if u else(F(),F());V=q(r['umax'])if v else(F(),F())
                moments=[rational(mul(q(p['moment'][i]),U))+rational(mul(q(r['moment'][i]),V))for i in(0,1)]
                norm=rational(mul(q(p['c0']),mul(U,U)))+rational(mul(q(r['c0']),mul(V,V)))
                L,M=moments
                z=1-2*(L/a+M/b)-(F(1,a)+F(1,b))*norm-F(1,a*b)*(L-M)**2
                assert z>=F(1,18);values.append(z)
                expected.append(dict(c=p['name'],d=r['name'],arc_lower_bounds=[a,b],upper_reciprocals=[u,v],normalized_determinant=[str(z),'0']))
    assert expected==c['rectangle_checks']and len(expected)==272 and min(values)==F(1,18)
    t=d['tree_mixed_replacement'];costs={}
    for p in t['profiles']:
        h=q(p['minimum_height'])
        costs[p['name']]=2*rational(div(q(p['moment']),h))+rational(div(q(p['c0']),mul(h,h)))
    assert max(costs.values())==4 and max(v for k,v in costs.items()if k!='inner_next3')==3
    seen=set();counts=Counter()
    for r in t['rows']:
        key=r['P'],r['c'],r['d'];assert key not in seen;seen.add(key);counts[r['kind']]+=1
        if r['P']==8 and 'inner_next3'in(r['c'],r['d']):assert r['kind']=='immediate_path7'
        else:
            z=1-F(costs[r['c']]+costs[r['d']],r['P'])
            assert r['kind']=='positive_gap'and q(r['gap'])==(z,F())
            assert z>=F(1,4)if r['P']==8 else z>=F(1,9)
    assert len(seen)==128 and counts==dict(immediate_path7=15,positive_gap=113)
    cat=read('fresh_high_rank_quasilanner_census_20261017')['diagrams'];graphs={}
    for x in cat:
        if x['order']!=10:continue
        g=[set()for _ in range(10)]
        for i,j,m in x['edges']:g[i].add(j);g[j].add(i)
        graphs[x['index']]=g
    triples=d['triangle_two_path7']['rows'];seen=set();worst=0
    for r in triples:
        key=r['diagram'],tuple(r['marks']);assert key not in seen;seen.add(key);g=graphs[r['diagram']];lengths=[]
        for path in r['paths']:
            start,end=path['pair'];forbidden=set(r['marks'])-{start,end};dist={start:0};queue=deque([start])
            while queue:
                i=queue.popleft()
                for j in g[i]-forbidden:
                    if j not in dist:dist[j]=dist[i]+1;queue.append(j)
            assert path['length']==dist.get(end)
            if end in dist:lengths.append(dist[end])
            if path['path']is not None:
                p=path['path'];assert p[0]==start and p[-1]==end and not(set(p)&forbidden)
                assert all(b in g[a]for a,b in zip(p,p[1:]))
        assert len(lengths)>=2;second=sorted(lengths)[1];assert second<=7;worst=max(worst,second)
    assert seen=={(i,marks)for i in graphs for marks in combinations(range(10),3)}and len(seen)==360 and worst==7
    o=d['order16_local_obstruction'];assert o['finite_volume']is False
    assert o['signature']==[16,1,1]and o['ray_counts']==dict(finite=35,ideal=7,exterior=24)
    assert len(o['extreme_rays'])==66 and Counter(r['kind']for r in o['face_rays'])==dict(finite=3,ideal=1)
    assert all(r['ideal7_whole15_payment']==0 for r in o['corners'])
    return dict(independent_rational_cycle_values=272,minimum_cycle_gap='1/18',independent_tree_cases=128,
        independent_graph_shortest_path_triples=360,order16_saved_ray_counts=o['ray_counts'],
        scope='Independent Q(sqrt2) rectangle arithmetic and tree costs, independent BFS distances on bound graphs; order16 cone census is producer replay plus saved-count consistency, not duplicated here. No ambient bound.')

def local_checks(d):
    base=read('fresh_mixed_raw39_full_model_20261017');cert=next(c for c in base['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert d['face_multipliers']==cert['face_multipliers']and d['initial_vertex_multipliers']==cert['vertex_multipliers']
    alpha=masses(d['face_multipliers']);old=masses(d['initial_vertex_multipliers']);new=masses(d['vertex_multipliers']);rebuilt=dict(old)
    weights=defaultdict(F);cols={c['index']:c for c in d['literal_block_orbits']}
    for r in d['positive_simplex_features']:
        for j,m in r['block_orbits']:weights[j]+=alpha[r['face']]*m
    assert set(cols)==set(weights)and len(cols)==5373
    terms=blocks=0;seen=set()
    for r in d['local_cones']:
        v=r['vertex'];assert v not in seen;seen.add(v)
        cs={c['index']:c for c in r['literal_block_orbits']};assert all(cols[j]==c and c['vertex']==v for j,c in cs.items())
        assert masses(r['required_marginals'])=={j:weights[j]for j in cs}
        cover=defaultdict(F);mass=F()
        for term in r['certificate']['packing_mixture']:
            a=F(term['mass']);assert a>0;mass+=a;used=0;terms+=1
            for j,mask in term['blocks']:
                assert mask in cs[j]['literal_masks']and 0<mask<(1<<13)and not used&mask
                used|=mask;cover[j]+=a;blocks+=1
        assert all(cover[j]>=weights[j]for j in cs)
        assert mass==F(r['certificate']['exact_cover_mass'])
        rebuilt[v]=max(old.get(v,F()),mass)
        assert rebuilt[v]==F(r['repaired_vertex_mass'])and old.get(v,F())==F(r['old_vertex_mass'])
    assert rebuilt==new and len(seen)==804 and(terms,blocks)==(4658,20074)
    assert not d['open_local_optima']and d['global_optimizer_invocations']==d['replay_optimizer_invocations']==0
    for r in d['new_exact_charge_columns']:
        z=load(alpha,r['positive_simplex_coefficients']);assert z==F(r['exact_load'])>old[r['vertex']]
        assert z<=new[r['vertex']]
    assert len(d['new_exact_charge_columns'])==78 and sum(new[i]>old.get(i,F())for i in new)==164
    assert F(d['objective_numerator'])==F(cert['objective_numerator'])
    assert sum(new.values(),F())==F(d['normalizer'])and F(d['objective_numerator'])/sum(new.values(),F())==F(d['bound'])>12
    return dict(independently_checked_mixtures=804,packing_terms=terms,literal_blocks=blocks,new_columns=78,bound=d['bound'],
        scope='Literal masks, disjointness, reconstructed orbit marginal coverage and exact mass/column composition; orbit completeness and all433753 model columns use bound producer replay.')

def joint_checks(d):
    base=read(LOCAL);assert d['face_multipliers']==base['face_multipliers']and d['initial_vertex_multipliers']==base['vertex_multipliers']
    alpha=masses(d['face_multipliers']);old=masses(d['initial_vertex_multipliers']);new=masses(d['vertex_multipliers']);rebuilt=dict(old)
    stages=d['complete_family_stages'];assert[len(s['vertex_required_masses'])for s in stages]==[821]*3
    for s in stages:
        assert s['literal_subsets']==821*comb(13,s['subset_size'])
        assert s['inequalities']==s['subset_orbits']*(15 if s['subset_size']==5 else 1)
        for v,a in s['vertex_required_masses']:rebuilt[v]=max(rebuilt[v],F(a))
    assert rebuilt==new and sum(s['inequalities']for s in stages)==3318122
    assert len(d['critical_columns'])==68 and sum(new[v]>old[v]for v in old)==21
    for r in d['critical_columns']:
        z=load(alpha,r['positive_face_coefficients']);cost=F(r['vertex_cost']);v=r['vertex']
        assert z==F(r['load'])and z/cost==F(r['required_vertex_mass'])>old[v]
        assert z<=cost*new[v]
    qs=d['actual_quad_features'];assert len(qs)==519 and sum(len(r['corners'])for r in qs)==3504
    assert sum(r['provenance_kind']=='existing_explicit_map'for r in qs)==461
    outside=d['zero_coefficient_unselected_Q_rows'];assert len(outside)==1654 and not(set(outside)&{r['face']for r in qs})
    assert F(d['objective_numerator'])==F(base['objective_numerator'])
    assert sum(new.values(),F())==F(d['normalizer'])and F(d['objective_numerator'])/sum(new.values(),F())==F(d['bound'])>12
    return dict(selected_Q_classes=519,saved_joint_inequalities=3318122,independently_checked_violated_loads=68,bound=d['bound'],
        scope='Saved family counts, all exported loads and simultaneous max-mass composition; no duplicate3318122-inequality geometry scan.')

def raw24_checks(d):
    base=read(JOINT);assert d['face_multipliers']==base['face_multipliers']and d['initial_vertex_multipliers']==base['vertex_multipliers']
    alpha=masses(d['face_multipliers']);old=masses(d['initial_vertex_multipliers']);new=masses(d['vertex_multipliers']);rebuilt=dict(old)
    assert d['raw_source']==24 and d['canonical_core']==20 and d['dotted_square']==['23/2','9/2']
    assert F(23,2)**2-5*F(9,2)**2==31 and not d['omitted_ambiguous_raw24_rows']
    records=d['forest_records'];assert len(records)==219
    assert len(d['selected_faces'])==45 and sum(len(f['corners'])for f in d['selected_faces'])==360
    assert sum(len(r['literal_stars'])for r in records)==4522
    assert sum(r['tested_literal_pairs']for r in records)==126921
    assert sum(len(r['incompatibility_certificates'])for r in records)==126903
    assert sum(r['unknown_literal_pairs']for r in records)==18
    for r in records:
        assert r['tested_literal_pairs']==comb(len(r['literal_stars']),2)
        assert r['tested_literal_pairs']==len(r['incompatibility_certificates'])+r['unknown_literal_pairs']
    violations=0
    for r in d['columns']:
        z=load(alpha,r['positive_face_coefficients']);v=r['vertex'];cost=F(r['vertex_cost']);assert cost>0
        assert z==F(r['load'])and F(r['baseline_vertex_mass'])==old.get(v,F())
        assert r['violated']==(z>cost*old.get(v,F()));violations+=r['violated']
        rebuilt[v]=max(rebuilt.get(v,F()),z/cost)
    assert len(d['columns'])==380 and violations==9 and rebuilt==new
    assert sum(new[v]>old.get(v,F())for v in new)==7
    r=d['fixed_witness_repair'];assert r['face_masses_unchanged']and r['optimizer_invocations']==0
    assert r['vertex_updates']==[[v,str(old.get(v,F())),str(new[v])]for v in sorted(new)if new[v]>old.get(v,F())]
    assert F(r['objective_numerator'])==F(base['objective_numerator'])
    assert sum(new.values(),F())==F(r['normalizer'])and F(r['objective_numerator'])/sum(new.values(),F())==F(r['bound'])>12
    return dict(saved_actual_Q_classes=45,saved_literal_pairs=126921,independent_column_loads=380,violations=9,bound=r['bound'],
        scope='Saved literal-domain counts, norm31 and all380 exact loads plus max-mass repair; no independent pair-field or orbit census.')

def mixed_checks(d):
    base=read(RAW24);assert d['face_multipliers']==base['face_multipliers']and d['initial_vertex_multipliers']==base['vertex_multipliers']
    alpha=masses(d['face_multipliers']);old=masses(d['initial_vertex_multipliers']);new=masses(d['vertex_multipliers']);rebuilt=dict(old);without=dict(old)
    violations=mixed_columns=mixed_violations=0
    for c in d['columns']:
        z=load(alpha,c['positive_face_coefficients']);v=c['vertex'];cost=F(c['vertex_cost']);assert cost>0
        assert z==F(c['load'])and F(c['baseline_vertex_mass'])==old.get(v,F())
        assert c['violated']==(z>cost*old.get(v,F()));violations+=c['violated']
        rebuilt[v]=max(rebuilt.get(v,F()),z/cost)
        mixed=24 in c['raw_families']and len(c['raw_families'])>1
        mixed_columns+=mixed;mixed_violations+=mixed and c['violated']
        if not mixed:without[v]=max(without.get(v,F()),z/cost)
    assert rebuilt==new and(len(d['columns']),violations,mixed_columns,mixed_violations)==(1793,148,357,139)
    assert sum(new[v]>old.get(v,F())for v in new)==120 and sum(new[v]>without.get(v,F())for v in new)==115
    rep=d['fixed_witness_repair'];assert rep['face_masses_unchanged']and rep['optimizer_invocations']==0
    update=lambda mu:[[v,str(old.get(v,F())),str(mu[v])]for v in sorted(mu)if mu[v]>old.get(v,F())]
    assert update(new)==rep['vertex_updates']and update(without)==rep['without_new_raw24_mixed_updates']
    num=F(rep['objective_numerator']);assert num==F(base['fixed_witness_repair']['objective_numerator'])
    assert sum(new.values(),F())==F(rep['normalizer'])and num/sum(new.values(),F())==F(rep['bound'])<12
    assert num/sum(without.values(),F())==F(rep['without_new_raw24_mixed_bound'])>12
    assert not d['omitted_ambiguous_selected_rows']and d['raw_sources']==[23,24,25,39]
    families={}
    for raw in(23,24,25,39):
        fs=[f for f in d['selected_faces']if f['raw_family']==raw];cs=[c for f in fs for c in f['corners']]
        families[str(raw)]=dict(classes=len(fs),corners=len(cs),forests=len({c['vertex']for c in cs}),orientations=sum(c['literal_orientations']for c in cs))
    assert families==d['family_counts']and[(families[str(r)]['classes'],families[str(r)]['corners'])for r in(23,24,25,39)]==[(71,710),(45,360),(54,432),(11,44)]
    rs=d['forest_records'];assert len(rs)==458 and sum(len(r['literal_stars'])for r in rs)==17888
    assert sum(r['tested_literal_pairs']for r in rs)==950774 and sum(len(r['incompatibility_certificates'])for r in rs)==933174
    assert sum(r['raw24_cross_pairs']for r in rs)==339222 and sum(r['raw24_cross_exclusions']for r in rs)==335420
    for raw,n in((23,11),(24,31),(39,59)):
        a,b=map(F,d['dotted_squares'][str(raw)]);assert a*a-5*b*b==n
    return dict(actual_Q_classes=181,saved_raw24_cross_pairs=339222,independent_column_loads=1793,
        mixed_columns=357,mixed_violations=139,bound=rep['bound'],without_new_mixed_bound=rep['without_new_raw24_mixed_bound'],
        scope='Independent exact all-column load/max-mass composition and saved domain counts; no duplicate pair-field geometry census and no optimized model claim.')

def bottleneck_checks(d):
    nodes=set(d['cusp']['nodes']);assert len(nodes)==17
    rows=d['faces'];assert len(rows)==136 and{tuple(r['marked_pair'])for r in rows}==set(combinations(sorted(nodes),2))
    mandatory=set();comparison=set();counts=Counter();changed=0
    for r in rows:
        fixed=short=0
        for c in r['corners']:
            if c['target_D16']:continue
            if c['ideal']:a=b=1;assert len(c['component'])==9
            else:
                b=int(c['distance']is not None and c['distance']<=9)
                a=int(c['distance']is not None and(c['distance']<=9 or len(c['component'])<=17))
            assert(a,b)==(c['fixed_payment'],c['distance9_only_payment']);fixed+=a;short+=b
        need=5-r['polygon'];deficit=max(0,need-fixed);other=max(0,need-short)
        assert need==r['required_face_payment']and deficit==r['D16_payment_deficit']and other==r['distance9_only_D16_deficit']
        pair=tuple(r['marked_pair']);counts[r['polygon'],deficit]+=1;changed+=deficit!=other
        if deficit:mandatory.add(pair)
        if other:comparison.add(pair)
    assert counts=={(3,1):60,(3,0):12,(4,1):49,(4,0):15}and len(mandatory)==109 and len(comparison)==122 and changed==13
    actions=[{int(k):v for k,v in g.items()}for g in d['automorphisms']];assert len(actions)==8
    # Every action is checked on the literal double-fork graph.
    core=[(0,1),(0,2),(0,4),(2,3),(4,5),(5,6),(6,7),(7,8),(8,9)]
    graph={tuple(sorted(e))for e in core+[(a+10 if a<9 else 9,b+10 if b<9 else 9)for a,b in core]if set(e)<=nodes}
    for g in actions:
        assert set(g)==set(g.values())==nodes
        assert{tuple(sorted((g[a],g[b])))for a,b in graph}==graph
    assert len(d['orbit_demands'])==136
    for r in d['orbit_demands']:
        g=actions[r['automorphism']];pair=tuple(sorted(g[i]for i in r['pair']))
        assert pair==tuple(r['mandatory_image'])and pair in mandatory
    assert d['invariant_D16_cost_lower_bound']==136 and F(d['normalized_cost'])==8
    return dict(independently_checked_face_deficits=136,full_baseline_mandatory_pairs=109,
        distance9_only_mandatory_pairs=122,automorphisms=8,orbit_demands=136,invariant_cost=136,
        scope='Independent fixed-payment recomputation and literal automorphism/closure checks; actual finite-volume chamber and face completeness are inherited.')

def construction_checks(exact,scout):
    assert exact['finite_type']=='H4+H4'and exact['exact_coordinate_identities']==10 and exact['native_squares']==[2,2]
    assert exact['support']==[0,1,2,3,5,6,7,8]and exact['second_endpoint_support']==[0,1,2,3,5,6,7,30]
    assert [r['lost_finite_count']for r in exact['deletion_debits']]==[47,114,0,11]
    assert exact['deletion_debits'][0]['compact_H8_certificates']==[dict(support=[11],finite_vertices=41)]
    t=scout['topology'];e=Counter(tuple(sorted(p))for tri in t['triangles']for p in combinations(tri,2))
    boundary={tuple(sorted(p))for p in zip(t['boundary'],t['boundary'][1:]+t['boundary'][:1])}
    assert {p for p,c in e.items()if c==1}==boundary and all(c in(1,2)for c in e.values())
    assert len(e)==28 and len(t['triangles'])==17 and t['vertices']==12 and 12-28+17==1
    assert not scout['unresolved']and not scout['complete']and len(scout['seeds'])==20
    assert scout['counts']['matching_B_tests']==scout['counts']['matching_B_C_dotted']==1985
    for r in scout['seeds']:
        assert r['boundary']['counts']==dict(finite=24,exterior=8)and not r['boundary']['warnings']
        assert r['boundary']['open_ends']==r['PD_incidence']['open_ends']==12 and r['PD_incidence']['finite']==24
        assert len(r['native_pair_checks'])==10 and all(p['relation']not in(None,'parallel_or_unresolved')for p in r['native_pair_checks'])
    return dict(saved_exact_coordinate_identities=10,deletion_certificate_debits=[47,114,0,11],
        exact_disk_triangles=17,numerical_prefixes=20,numerical_blocked_B_joins=1985,
        scope='Exact saved identity/debit consistency and abstract disk counts; numerical ledgers remain numerical and no new chamber is promoted.')

def small_checks():
    return dict(ideal7=ideal_checks(read(IDEAL)),bottleneck=bottleneck_checks(read(BOTTLENECK)),local=local_checks(read(LOCAL)),joint=joint_checks(read(JOINT)),
        raw24=raw24_checks(read(RAW24)),mixed=mixed_checks(read(MIXED)),construction=construction_checks(read(IDENTITY),read(SCOUT)),
        pending_inputs=list(PENDING_INPUTS),scope='Independent small exact checks and saved certificate composition. Unrestricted bound97 and compact-arithmetic21 unchanged; no finite-volume construction.')


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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_IDEAL7_AND_RAW24',
                scope='Archive integrity, independent small exact local-payment checks, fixed arithmetic witness repairs and scoped construction diagnostics; no new dimension bound or finite-volume example.',
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
