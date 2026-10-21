#!/usr/bin/env python3
"""Frozen incremental archive for dimension80 and joint actual-corner capacities."""
import argparse,ast,json,re
from collections import Counter,deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/dimension80_and_joint_corner_frontier_20261017.md'
PREVIOUS='results/ideal6_local_and_large_core_model_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_ideal6_local_and_large_core_model_frontier_20261017.py'
OUT=ROOT/'results/dimension80_and_joint_corner_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
D14='fresh_affine_d14_long_pair_transfer_20261017'
FINITE='fresh_broad_finite_half_coupons_20261017'
LOCAL='fresh_ideal6_combined_payments_20261017'
BOUND='fresh_ideal6_piecewise_dimension80_20261017'
REVIEW=BOUND+'_review'
ARITH='fresh_raw87_88_89_actual_corners_20261017'
CONSTRUCTION='fresh_double_flip_completion_gate_20261017'
STEMS=(D14,FINITE,LOCAL,BOUND,REVIEW,ARITH,CONSTRUCTION)
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

def components(g,ids):
    todo=set(ids);out=[]
    while todo:
        ds=distances(g,min(todo),ids);todo-=ds.keys();nodes=set(ds)
        degree={i:sum(bool(g[i][j])for j in nodes if i!=j)for i in nodes}
        branch={i for i in nodes if degree[i]==3}
        leaves={i for i in nodes if degree[i]==1 and any(g[i][j]for j in branch)}
        if max(degree.values())<=2:kind='A'+str(len(nodes))
        elif len(branch)==1 and len(leaves)==2:kind='D'+str(len(nodes))
        else:kind='affineD'+str(len(nodes)-1)
        out.append(dict(type=kind,nodes=sorted(nodes),branch=sorted(branch),short_leaves=sorted(leaves)))
    return out


def local_checks():
    d=read(D14);cases=d['equality']['rows'];seen=set();counts=Counter();finite_corners=0
    for r in cases:
        a,b=r['marks'];L,R=r['leaf_attachments'];side=r['spectator_attached_to']
        key=(a,b,L,R,side);assert key not in seen;seen.add(key)
        assert a in(0,11,12)and b in(10,13,14)and L in(11,12)and R in(13,14)and side in(15,16)
        g=[[F(2 if i==j else 0)for j in range(18)]for i in range(18)]
        edges=[(i,i+1)for i in range(10)]+[(0,11),(0,12),(10,13),(10,14),(15,L),(16,R),(17,side)]
        for i,j in edges:g[i][j]=g[j][i]=F(-1)
        assert g==[[F(v)for v in row]for row in r['gram']]
        z=list(map(F,r['radical']));assert any(z)
        assert all(sum(g[i][j]*z[j]for j in range(18))==0 for i in range(18))
        W=sorted(set(range(15))-{a,b})+[17];assert W==r['white'];pd(g,W)
        units=0;coupon=F()
        assert [c['marks']for c in r['corners']]==[[a,b],[b,15],[15,16],[16,a]]
        for j,c in enumerate(r['corners']):
            marks=c['marks'];ids=sorted(set(W)|set(marks));cc=components(g,ids)
            for part in cc:
                part['marks']=[i for i in marks if i in part['nodes']]
                part['fork_marks']=[i for i in part['marks']if i in part['branch']+part['short_leaves']]
            assert cc==c['components']
            if not j:
                assert cc[0]['type']=='affineD14'and len(cc[0]['nodes'])==15
                nu=[F(1)]*11+[F(1,2)]*4
                assert all(sum(g[i][k]*nu[k]for k in range(15))==0 for i in range(15))
                pd(g,list(range(14)))
                assert c['finite_small17_unit']==0 and F(c['half_coupon'])==0
                continue
            pd(g,ids);finite_corners+=1
            p=next(p for p in cc if marks[0]in p['nodes']);q=next(p for p in cc if marks[1]in p['nodes'])
            unit=int(p is q and len(p['nodes'])<=17)
            cp=F()
            if p is not q:
                small=any(len(t['nodes'])<=3 for t in(p,q))
                large=any(t['type']in('D12','D14')and t['fork_marks']for t in(p,q))
                if small and large:cp=F(1,2)
            assert c['finite_small17_unit']==unit and F(c['half_coupon'])==cp
            units+=unit;coupon+=cp
        assert units==r['finite_units']and coupon==F(r['coupon_total'])and units+coupon>=1
        counts[(units,str(coupon))]+=1
    assert len(seen)==72 and finite_corners==216
    assert counts==Counter({(0,'1'):32,(1,'1/2'):32,(3,'0'):8})
    budget=read(FINITE);checks=0
    for r in budget['small_components']+budget['large_components']:
        n=r['rank'];base=F(n*(n-1),2);assert base==F(r['base_cost'])
        assert 9*n-base==F(r['slack'])
        if n<=4:assert F(r['slack'])>=F(15*n,2)
        else:assert n>=12 and F(r['slack'])>=24 and r['paid_mark_count']<=3
    for r in budget['rows']:
        n=r['dimension'];expected=[]
        for t in range(n//12+1):
            for a in sorted({0,n-12*t}):
                val=F(a*(3*t-15),2)-24*t;assert val<=0
                expected.append([a,t,str(val)]);checks+=1
        assert expected==r['endpoint_checks']
    assert [r['dimension']for r in budget['rows']]==list(range(1,109))
    old=read('fresh_ideal6_small_factor_frontier_20261017');local=read(LOCAL)
    records={(r['order'],r['kind']):r for r in old['records']};costs={};removed=0
    for r in local['small_order_pair_sets']:
        key=r['order'],r['kind'];prior=records[key];pairs=set(map(tuple,prior['paid_pairs']));drop=set()
        if key==(15,'D'):
            g=prior['target'];branch=[i for i,row in enumerate(g)if sum(bool(v)for v in row)==3]
            assert len(branch)==2
            ends=[[i]+[j for j in range(15)if g[i][j]and sum(bool(v)for v in g[j])==1]for i in branch]
            assert list(map(len,ends))==[3,3]
            drop={tuple(sorted((a,b)))for a in ends[0]for b in ends[1]};assert len(drop)==9
            assert drop<=pairs;pairs-=drop;removed+=len(drop)
        assert pairs==set(map(tuple,r['paid_pairs']))and drop==set(map(tuple,r['newly_removed_pairs']))
        assert len(pairs)==r['cost'];costs.setdefault(str(key[0]),{})[key[1]]=len(pairs)
    assert removed==9 and costs==local['small_order_costs']=={'15':dict(A=90,B=100,C=95,D=96),'16':dict(A=96,B=107,C=103,D=111),'17':dict(A=102,B=111,C=108,D=127)}
    return dict(D14_cases=72,white_PD_checks=72,finite_corner_PD_checks=finite_corners,
        coupon_classes=[dict(units=u,half_total=q,cases=n)for(u,q),n in sorted(counts.items())],
        finite_budget_dimensions=108,finite_budget_endpoint_checks=checks,
        factor_pair_sets=12,new_D14_removed_pairs=9,factor_costs=costs)


def scalar_checks():
    p=read(BOUND);q=read(REVIEW);rs=q['reconstruction']['profiles']
    assert p['dimension_bound']==q['dimension_bound']==80
    assert p['unchanged_compact_facet_ambient_bound']==q['unchanged_compact_facet_ambient_bound']==79
    assert p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    assert p['excluded_dimensions']==list(range(81,87))and p['inherited_dimension_bound']==86
    assert len(p['reduced_profiles'])==len(rs)==5008
    pp={tuple(r['ranks']):r for r in p['reduced_profiles']};assert len(pp)==5008
    for r in rs:
        assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==pp[tuple(r['ranks'])]
        h=r['h_polynomial'];d=r['link_dimension'];V=r['link_vertices']
        assert len(h)==d+1 and h==h[::-1]and sum(h)==V
        G=sum(a*min(j,d-j)for j,a in enumerate(h));assert G==r['G']==r['independent_cumulative_G']
        curv=F(5*G-d*V,2*V);cost=sum((F(x,2)if x<=13 else F(20,3)if x==14 else F(111,16)if x==15 else F(127,17)if x==16 else F(95,13)for x in r['ranks']),F())
        assert curv==F(r['normalized_curvature'])and cost==F(r['normalized_cost'])
        assert curv-cost==F(r['margin'])>0
    assert p['dimension_minima']==q['reconstruction']['dimension_minima']
    assert min(F(r['margin'])for r in rs)==F(940675,11294304)
    assert F(p['boundary_at_n80']['margin'])==F(-1131729,7529536)
    assert F(13,2)+F(111,16)-2*F(20,3)==F(5,48)>0
    assert F(20,3)+F(127,17)-2*F(111,16)==F(107,408)>0
    return dict(dimension_bound=80,compact_facet_ambient_bound=79,matched_independent_profiles=5008,
        minimum_margin='940675/11294304',boundary_margin='-1131729/7529536')


def arithmetic_checks():
    d=read(ARITH);cols=d['columns'];alpha=mass(d['face_multipliers']);old=mass(d['initial_vertex_multipliers']);new=mass(d['vertex_multipliers'])
    model=read('fresh_raw178_296_297_full_model_20261017')
    baseline=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert d['face_multipliers']==baseline['face_multipliers']and d['initial_vertex_multipliers']==baseline['vertex_multipliers']
    assert len(alpha)==4410 and all(new.get(i,F())>=a for i,a in old.items())
    chosen={r['face']:r for r in d['selected_faces']};assert len(chosen)==89
    assert Counter(r['provenance_kind']for r in chosen.values())==dict(existing_explicit_map=51,unique_saved_context=38)
    assert all(r['spectator_rank']==7 and r['raw_core_order']==8 and len(r['corners'])==8 for r in chosen.values())
    for r in chosen.values():
        if r['provenance_kind']=='unique_saved_context':assert len(r['source_contexts'])==1
    assert len(d['forest_records'])==380
    seen_cols=set();reasons=Counter();star_count=orbit_count=0;pair_refs=0
    for r in d['forest_records']:
        stars=r['literal_stars'];groups=r['corner_orbits'];vi=r['vertex'];N=len(stars)
        assert [a['index']for a in stars]==list(range(N))
        assert [a['index']for a in groups]==list(range(len(groups)))
        assert sorted(i for g in groups for i in g['literal_stars'])==list(range(N))
        occupied=[set(map(tuple,a['occupied_pairs']))for a in stars]
        bad=set()
        for cert in r['incompatibility_certificates']:
            i,j,code=cert[:3];assert 0<=i<j<N and(i,j)not in bad;bad.add((i,j));reasons[code]+=1
            if code==0:assert len(cert)==3 and occupied[i]&occupied[j]
            else:
                assert len(cert)==6 and code in(1,2,3,4,5)and cert[3]in(0,1)and cert[4]in(0,1)
                e=r['radical_expressions'][cert[5]]
                if code==1:assert e['sign_product']>0
                elif code==2:assert e['sign_product']<0<e['sign_product_plus_one']
                elif code==3:assert e['arithmetic_gate']['two_cycle_integral']is False
                elif code==4:assert e['arithmetic_gate']['two_cycle_in_fixed_K']is False
        assert r['tested_literal_pairs']==N*(N-1)//2
        assert r['unknown_literal_pairs']==N*(N-1)//2-len(bad)
        def all_bad(nodes):
            nonlocal pair_refs
            for i,j in combinations(sorted(nodes),2):assert(i,j)in bad;pair_refs+=1
        for g in groups:
            coeff=Counter(a['face']for a in g['origins'])
            assert dict(g['positive_face_coefficients'])==dict(coeff)and len(coeff)==1
            if g['capacity_one_proved']:all_bad(g['literal_stars'])
        gedges=set(map(tuple,r['group_incompatibility_edges']))
        for i,j in gedges:
            assert i<j and groups[i]['capacity_one_proved']and groups[j]['capacity_one_proved']
            for a in groups[i]['literal_stars']:
                for b in groups[j]['literal_stars']:assert tuple(sorted((a,b)))in bad;pair_refs+=1
        for ci in r['columns']:
            c=cols[ci];assert ci not in seen_cols and c['index']==ci and c['vertex']==vi;seen_cols.add(ci)
            if c['kind']=='capacity_one_orbit_clique':
                ix=c['orbit_clique'];assert c['vertex_cost']==1
                assert all(groups[i]['capacity_one_proved']for i in ix)
                assert all(tuple(sorted(p))in gedges for p in combinations(ix,2))
            else:
                assert c['kind']=='proper_compatibility_coloring';ix=c['orbit_indices']
                parts=c['literal_coloring'];nodes=sorted(i for j in ix for i in groups[j]['literal_stars'])
                assert sorted(i for part in parts for i in part)==nodes
                assert len(parts)==c['vertex_cost']
                for part in parts:all_bad(part)
            coeff=Counter()
            for i in ix:coeff.update(dict(groups[i]['positive_face_coefficients']))
            assert sorted(coeff.items())==list(map(tuple,c['positive_face_coefficients']))
            assert sorted({chosen[i]['raw_family']for i in coeff})==c['raw_families']
        star_count+=N;orbit_count+=len(groups)
    assert seen_cols==set(range(1002))and star_count==19820 and orbit_count==712
    assert reasons==Counter({0:389179,1:741367,2:495111,3:574369,4:330544})
    assert sum(reasons.values())==2530570==d['counts']['exact_literal_exclusions']
    expected=dict(old);single=dict(old);violated=mixed=0
    for c in cols:
        load=sum((F(a)*alpha[i]for i,a in c['positive_face_coefficients']),F());cost=F(c['vertex_cost']);vi=c['vertex']
        assert cost>0 and load==F(c['load'])and F(c['baseline_vertex_mass'])==old.get(vi,F())
        v=load>cost*old.get(vi,F());assert v==c['violated'];violated+=v;mixed+=len(c['raw_families'])>1
        expected[vi]=max(expected.get(vi,F()),load/cost)
        if len(c['raw_families'])==1:single[vi]=max(single.get(vi,F()),load/cost)
        assert load<=cost*new.get(vi,F())
    assert expected==new
    rep=d['fixed_witness_repair'];updates=[[i,str(old.get(i,F())),str(a)]for i,a in sorted(new.items())if a>old.get(i,F())]
    assert updates==rep['vertex_updates']and len(updates)==75
    assert violated==113 and mixed==286
    assert sum(new.values(),F())==F(rep['normalizer'])
    assert F(rep['objective_numerator'])==F(baseline['objective_numerator'])
    ratio=F(rep['objective_numerator'])/sum(new.values(),F())
    assert ratio==F(rep['bound'])>12 and rep['optimizer_invocations']==0 and not rep['bound_below12']
    assert F(rep['without_mixed_family_bound'])==F(rep['objective_numerator'])/sum(single.values(),F())
    assert sum(new.get(i,F())>single.get(i,F())for i in new)==21
    return dict(columns=1002,forests=380,literal_stars=star_count,corner_orbits=orbit_count,
        saved_incompatibilities=sum(reasons.values()),combinatorial_pair_references=pair_refs,
        violated_columns=violated,vertex_raises=75,mixed_additional_raises=21,bound=str(ratio),
        scope='Independent coefficient, exclusion-reference, clique/coloring and exact monotone-repair checks; does not recompute biquadratic product signs, active embeddings or optimization.')


def interval(xs):
    lo,hi=map(F,xs);assert lo<=hi;return lo,hi


def theta_interval(label,xs):
    # Saved intervals enclose the NEGATIVE norm-2 Coxeter product.
    lo,hi=interval(xs);a,b=-hi,-lo
    if label==2:assert a==b==0
    elif label==3:assert a==b==1
    elif label in(4,6):assert 0<a and a*a<=(2 if label==4 else 3)<=b*b
    elif label==5:assert a>1 and(2*a-1)**2<=5<=(2*b-1)**2
    elif label==7:
        def f(x):return x*x*x-x*x-2*x+1
        assert a>F(3,2)and f(a)<=0<=f(b)
    else:raise AssertionError(label)
    return lo,hi


def construction_checks():
    d=read(CONSTRUCTION);top=d['topology'];tris=top['triangles'];boundary=top['boundary']
    vertices=set(x for t in tris for x in t);edges=Counter(tuple(sorted(e))for t in tris for e in combinations(t,2))
    be={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert len(vertices)==13 and len(edges)==31 and len(tris)==19
    assert set(e for e,n in edges.items()if n==1)==be and all(n==(1 if e in be else 2)for e,n in edges.items())
    for v in vertices:
        link={}
        for tri in tris:
            if v in tri:
                a,b=[x for x in tri if x!=v];link.setdefault(a,set()).add(b);link.setdefault(b,set()).add(a)
        todo=[next(iter(link))];seen=set()
        while todo:
            u=todo.pop()
            if u not in seen:seen.add(u);todo.extend(link[u]-seen)
        assert seen==set(link)
        assert sorted(map(len,link.values()))==([1,1]+[2]*(len(link)-2)if v in boundary else[2]*len(link))
    assert tuple(sorted(('X','B')))in edges and tuple(sorted(('W','Z')))not in edges
    counts=Counter();rejections=Counter()
    for rec in d['records']:
        for which in('W_cases','B_cases'):
            for row in rec[which]:
                high=6 in row['labels'];counts['B'if which=='B_cases'else'W_high'if high else'W_low']+=1
                proof=row['proof'];assert 'children'not in proof
                if high:
                    lo,hi=interval(proof['cosine_box']);assert lo>0 and lo*lo<=3 and hi>=2
                else:assert proof['cosine_box']is None
                cert=proof['certificate'];assert cert['reason']=='both_norm_branches_excluded'
                info=cert['frame_certificate'];assert interval(info['kernel_square'])[0]>0 and interval(info['remaining_square'])[0]>=0
                for r in info['inverse_pivots']:
                    lo,hi=interval(r['pivot']);assert lo>0 or hi<0
                assert sorted(b['sign']for b in cert['branches'])==[-1,1]
                for b in cert['branches']:
                    counts['branches']+=1
                    if b['reason']=='old_native_contact':
                        r=b['rejection'];lo,hi=interval(r['product']);rejections[r['kind']]+=1
                        assert r['wall']in d['targets']
                        if r['kind']=='strict_obtuse':assert lo>0
                        else:
                            assert r['kind']=='strict_consecutive_angle_gap';a,beta=r['labels'];assert beta==a+1
                            lower=theta_interval(beta,r['lower_threshold']);upper=theta_interval(a,r['upper_threshold'])
                            assert lower[1]<lo<=hi<upper[0]
                    else:
                        assert b['reason']=='required_finite_WZ_angle_absent';rejections['finite_angle_absent']+=1
                        lo,hi=interval(b['product']);assert b['allowed_labels']==row['allowed_WZ']
                        assert [r['label']for r in b['separations']]==row['allowed_WZ']
                        for r in b['separations']:
                            lower,upper=theta_interval(r['label'],r['threshold'])
                            assert (hi<lower if r['side']=='strictly_below'else lo>upper)
    assert counts==dict(W_low=368,W_high=22,B=60,branches=900)
    assert rejections==dict(strict_obtuse=502,strict_consecutive_angle_gap=70,finite_angle_absent=328)
    return dict(fixed_prefixes=6,disk_counts=dict(vertices=13,edges=31,triangles=19),recipe_counts=dict(counts),
        rational_saved_branch_inequalities=dict(rejections),
        scope='Independent disk links and rational saved inequality/threshold checks only; physical interval enclosure construction and all-label domain completeness retain producer generation/replay and proof-review scope.')


def small_checks():
    local=local_checks();print('SMALL_D14_AND_BUDGET_PASS',flush=True)
    scalar=scalar_checks();print('SMALL_5008_SCALAR_MATCH_PASS',flush=True)
    arith=arithmetic_checks();print('SMALL_1002_CAPACITY_CHECKS_PASS',flush=True)
    construction=construction_checks();print('SMALL_900_SAVED_BRANCH_INEQUALITIES_PASS',flush=True)
    return dict(local=local,dimension=scalar,arithmetic=arith,construction=construction,
        scope='Independent saved-certificate consistency and rational checks. No repeat of universal supporting-white geometry, complete preferred embeddings, biquadratic contact algebra, physical interval enclosures or old millions of constraints.')


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
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_DIMENSION80_AND_JOINT_CORNERS',
                scope='Archive integrity and scoped independent exact checks of the completed unrestricted dimension80 theorem, D14 coupon corners, broad finite budget, actual-corner capacities and fixed-prefix interval exclusions. No new finite-volume example; arithmetic obstruction remains above12.',
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
