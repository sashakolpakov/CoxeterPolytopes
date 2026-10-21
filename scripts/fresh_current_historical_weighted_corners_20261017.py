#!/usr/bin/env python3
"""Exact weighted capacities at the current dual on the complete saved ancestry.

No wall-pair computation or global optimizer. Local LPs propose exact charge /
cover certificates; --verify re-enumerates cliques and checks saved rationals.
Every pair uses its first complete saved domain, even when that outcome is
unknown. A pair with no saved owner is also allowed in the conservative graph.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json,time
import networkx as nx
import fresh_weighted_actual_corner_capacities_20261017 as W

ROOT=W.ROOT
TAG='fresh_current_historical_weighted_corners_20261017'
GEOMETRY='fresh_raw110_mixed_actual_corners_20261017'
MODEL='fresh_weighted_raw110_full_model_20261017'
SUPPORT=MODEL+'_support'
NATIVE_FIELDS=('contacts','occupied_pairs','squared_timelike_heights','native_omitted_product')

def frozen(x):
    if isinstance(x,list):return tuple(map(frozen,x))
    if isinstance(x,dict):return tuple((k,frozen(v))for k,v in sorted(x.items()))
    return x

def origin_key(o):
    return (o['face'],o['provenance'],tuple(o['omitted']),o['vertex'],o['raw_family'])

def source_order():
    """Preorder, preserving each producer's explicit first-owner order."""
    order=[]
    def visit(tag):
        if tag in order:return
        order.append(tag)
        for row in W.read(tag)['forest_records']:
            for dom in row.get('inherited_pair_domains',[]):visit(dom['source'])
    visit(GEOMETRY)
    return order

def prepare():
    model=W.read(MODEL);support=W.read(SUPPORT)
    cert=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert cert['face_multipliers']==support['face_multipliers']
    assert cert['vertex_multipliers']==support['vertex_multipliers']
    alpha={i:F(x)for i,x in cert['face_multipliers']}
    rr=support['exact_restrictions'][0]
    assert rr['name']=='latest_248_actual_corner_domain_Q'
    restricted=rr['certificate'];ralpha={i:F(x)for i,x in restricted['face_multipliers']}
    assert all(alpha[i]==x for i,x in ralpha.items())
    order=source_order();qt={r['face']:r for r in support['quad_rows']}
    named=set(support['explicit_geometry_domains']['combined_rows'])
    selected={};owners={};historical=set();orders={}
    for tag in order:
        for family,n in W.read(tag).get('raw_core_orders',{}).items():
            family=int(family)
            assert family not in orders or orders[family]==n
            orders[family]=n
    for tag in order:
        for raw_record in W.read(tag)['selected_faces']:
            r=dict(raw_record)
            fi=r['face'];historical.add(fi)
            if fi not in alpha:continue
            q=qt[fi]
            assert q['provenance_kind']in('explicit','unique_saved_context')
            assert r['source_context']==q['source_context']
            assert r['canonical_prism_core']==q['canonical_core']
            assert r['raw_family']==q['raw_family']
            assert r['provenance']==q.get('provenance')
            # The face obligation is binomial(d,2), not the 2d prism vertices.
            nc=len(r['corners']);assert nc%2==0 and 2<=nc//2<=5
            assert q['required']==(nc//2)*(nc//2-1)//2
            assert len({tuple(o['omitted'])for o in r['corners']})==len(r['corners'])
            assert all(o['face']==fi and o['provenance']==r['provenance'] and
                       o['raw_family']==r['raw_family']for o in r['corners'])
            n=orders[r['raw_family']]
            sr=sum(2 if name.startswith('I2(')else int(name[1:])for name in r['spectator_forest'])
            assert r.get('raw_core_order',n)==n and r.get('spectator_rank',sr)==sr
            assert n+sr==15
            r.update(raw_core_order=n,spectator_rank=sr)
            if fi in selected:assert r==selected[fi],(fi,tag)
            else:selected[fi]=r;owners[fi]=tag
    assert len(named)==248 and len(named&set(alpha))==209
    assert named&set(alpha)<=set(selected)
    assert len(selected)==248 and len(set(selected)-named)==39
    assert len(alpha)==5053
    sources={f'scripts/{TAG}.py':W.digest(f'scripts/{TAG}.py')}
    for tag in order+[MODEL,SUPPORT,W.TAG]:
        d=W.read(tag)
        for p,h in d.get('sources',d.get('prepared',{}).get('sources',{})).items():assert W.digest(p)==h,p
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{tag}.{ext}';sources[p]=W.digest(p)
    return cert,restricted,alpha,ralpha,order,selected,owners,named,historical,sources

def union_forest(v,order,selected):
    """Merge complete orbits, not independently chosen literal presentations."""
    rows=[(tag,W.forests(tag)[v])for tag in order if v in W.forests(tag)]
    forest=rows[0][1]['forest'];assert all(r['forest']==forest for _,r in rows)
    refs=[];stars=[];keys={};groups=[];orbit_keys={};star_group={};present_origins={}
    for tag,row in rows:
        for g in row['corner_orbits']:
            ff={fi for fi,m in g['positive_face_coefficients']}
            assert len(ff)==1
            fi=next(iter(ff))
            if fi not in selected:continue
            origins={origin_key(o):o for o in g['origins']}
            assert len(origins)==len(g['origins'])
            expected={origin_key(o):o for o in selected[fi]['corners']if o['vertex']==v}
            assert set(origins)<=set(expected)
            assert all(o==expected[k]for k,o in origins.items())
            assert g['positive_face_coefficients']==[[fi,len(origins)]]
            nodes=[]
            for oldi in g['literal_stars']:
                star=row['literal_stars'][oldi];assert star['index']==oldi
                key=frozen(star['contacts'])
                if key not in keys:
                    i=len(stars);keys[key]=i;stars.append(star);refs.append([tag,oldi])
                else:
                    i=keys[key]
                    assert all(star[k]==stars[i][k]for k in NATIVE_FIELDS),(v,tag,oldi)
                nodes.append(i)
            assert len(nodes)==len(set(nodes)) and nodes
            ok=frozenset(nodes)
            if ok in orbit_keys:
                gi=orbit_keys[ok];oldg=groups[gi]
                assert oldg['face']==fi and oldg['origins']==origins
                oldg['source_orbits'].append([tag,g['index']])
            else:
                assert not any(i in star_group for i in nodes),(v,tag,'partial orbit overlap')
                gi=len(groups);orbit_keys[ok]=gi
                for i in nodes:star_group[i]=gi
                groups.append(dict(index=gi,face=fi,literal_stars=sorted(nodes),origins=origins,
                    positive_face_coefficients=g['positive_face_coefficients'],source_orbits=[[tag,g['index']]]))
                for k,o in origins.items():
                    assert k not in present_origins
                    present_origins[k]=o
    expected={origin_key(o):o for r in selected.values()for o in r['corners']if o['vertex']==v}
    assert present_origins==expected,(v,'missing or duplicated actual corners')
    if not groups:return None
    assert set(star_group)==set(range(len(stars)))
    # A source owns a pair whenever it contains both full native star keys.
    # Membership, rather than exclusion, determines precedence.
    domains=[];bits=[0]*len(stars)
    for tag,row in rows:
        mp={}
        for s in row['literal_stars']:
            key=frozen(s['contacts'])
            if key in keys:
                i=keys[key];assert i not in mp
                assert all(s[k]==stars[i][k]for k in NATIVE_FIELDS)
                mp[i]=s['index']
        if not mp:continue
        t=len(domains)
        for i in mp:bits[i]|=1<<t
        domains.append(dict(source=tag,new_to_old_star_indices=sorted(mp.items())))
    bad=set();owned_counts=Counter();owned_bad=Counter()
    for t,dom in enumerate(domains):
        inv={j:i for i,j in dom['new_to_old_star_indices']}
        assert len(inv)==len(dom['new_to_old_star_indices'])
        for a,b in W.relation(dom['source'],v):
            if a not in inv or b not in inv:continue
            i,j=sorted((inv[a],inv[b]));common=bits[i]&bits[j]
            if common&-common==1<<t:
                assert(i,j)not in bad;bad.add((i,j));owned_bad[t]+=1
    graph=nx.Graph();graph.add_nodes_from(range(len(stars)));unowned=0
    for i,j in combinations(range(len(stars)),2):
        common=bits[i]&bits[j]
        if common:owned_counts[(common&-common).bit_length()-1]+=1
        else:unowned+=1
        if(i,j)not in bad:graph.add_edge(i,j)
    for t,d in enumerate(domains):
        d['first_owned_pairs']=owned_counts[t];d['first_owned_exclusions']=owned_bad[t]
    assert sum(owned_counts.values())+unowned==len(stars)*(len(stars)-1)//2
    assert sum(owned_bad.values())==len(bad)
    W.relation.cache_clear()
    for g in groups:g['origins']=[g['origins'][k]for k in sorted(g['origins'],key=str)]
    return dict(vertex=v,forest=forest,literal_references=refs,corner_orbits=groups,
        first_owner_domains=domains,unowned_allowed_pairs=unowned,
        excluded_pairs=len(bad),inherited_unknown_pairs=sum(owned_counts.values())-len(bad)),graph,star_group

def vectors_for(graph,owner,n):
    patterns={};count=0;start=time.monotonic()
    for clique in nx.find_cliques(graph):
        count+=1;clique=tuple(sorted(clique));c=Counter(owner[i]for i in clique)
        vec=tuple(c.get(j,0)for j in range(n))
        if vec not in patterns or clique<patterns[vec]:patterns[vec]=clique
        if count>W.MAX_CLIQUES or time.monotonic()-start>W.ENUM_SECONDS:
            return None,dict(reason='literal_clique_enumeration_guard',enumerated_cliques=count)
    vv=sorted(patterns)
    maximal=[a for a in vv if not any(a!=b and all(x<=y for x,y in zip(a,b))for b in vv)]
    return (maximal,[patterns[a]for a in maximal],count,len(patterns)),None

def exact_profile(vectors,demand):
    positive=[i for i,x in enumerate(demand)if x]
    if not positive:return dict(weights=['0']*len(demand),mixture=[],optimal_mass='0')
    c=W.exact_lp([[v[i]for i in positive]for v in vectors],[demand[i]for i in positive])
    w=['0']*len(demand)
    for i,x in zip(positive,c['weights']):w[i]=x
    c['weights']=w
    return c

def run(verify=False,resume=False,prepare_only=False):
    start=time.monotonic()
    cert,restricted,alpha,ralpha,order,selected,owners,named,historical,sources=prepare()
    print('DOMAIN',len(historical),'historical classes;',len(selected),'positive;',
          len(set(selected)-named),'additional; sources',len(order),flush=True)
    work=sorted({o['vertex']for r in selected.values()for o in r['corners']})
    if prepare_only:
        counts=Counter()
        for v in work:
            u,graph,owner=union_forest(v,order,selected)
            counts.update(forests=1,orbits=len(u['corner_orbits']),literal_stars=len(owner),
                unowned_pairs=u['unowned_allowed_pairs'],inherited_unknown_pairs=u['inherited_unknown_pairs'])
        print('PREPARED',dict(counts),flush=True);return
    baseline=[cert,restricted];alphas=[alpha,ralpha]
    mus=[{v:F(x)for v,x in c['vertex_multipliers']}for c in baseline]
    newmus=[dict(m)for m in mus];labels=['full','named_248_class_Q_restriction']
    saved=W.read(TAG)if verify or resume else None
    if saved:assert saved['sources']==sources and saved['face_multipliers']==cert['face_multipliers']
    savedrows={r['vertex']:r for r in saved['forest_records']}if saved else {}
    rows=[];columns=[];opened=[];counts=Counter();columnkeys=set()
    for step,v in enumerate(work,1):
        u,graph,owner=union_forest(v,order,selected);groups=u['corner_orbits']
        enum,err=vectors_for(graph,owner,len(groups))
        if err:opened.append(dict(vertex=v,**err));continue
        vectors,reps,ncliques,npatterns=enum
        demands=[[sum((a.get(fi,F())*F(m)for fi,m in g['positive_face_coefficients']),F())for g in groups]for a in alphas]
        assert all(x>0 for x in demands[0])
        ccs=[]
        try:
            for k,demand in enumerate(demands):
                if v in savedrows:
                    old=savedrows[v]
                    assert old['union']==json.loads(json.dumps(u))
                    assert old['maximal_count_vectors']==list(map(list,vectors))
                    assert old['literal_representatives']==list(map(list,reps))
                    assert old['demands'][k]==list(map(str,demand))
                    c=old['certificates'][k]
                else:
                    assert not verify
                    c=exact_profile(vectors,demand)
                W.check_certificate(vectors,demand,c);ccs.append(c)
        except Exception as e:
            if verify:raise
            opened.append(dict(vertex=v,reason='exact_local_certificate_recovery',detail=str(e),vectors=len(vectors)))
            print('OPEN',v,repr(e),flush=True);continue
        for vec,clique in zip(vectors,reps):
            assert all(graph.has_edge(i,j)for i,j in combinations(clique,2))
            assert tuple(Counter(owner[i]for i in clique).get(j,0)for j in range(len(groups)))==vec
        for k,c in enumerate(ccs):
            value=F(c['optimal_mass']);newmus[k][v]=max(newmus[k].get(v,F()),value)
            if value>mus[k].get(v,F()):counts[labels[k]+'_violations']+=1
            co=Counter()
            for w,g in zip(map(F,c['weights']),groups):
                for fi,m in g['positive_face_coefficients']:co[fi]+=w*F(m)
            terms=sorted((fi,x)for fi,x in co.items()if x)
            loads=[sum((a.get(fi,F())*x for fi,x in terms),F())for a in alphas]
            assert loads[k]==value
            assert all(loads[j]<=F(ccs[j]['optimal_mass'])for j in range(2))
            key=(v,tuple(terms))
            if key not in columnkeys and any(loads[j]>mus[j].get(v,F())for j in range(2)):
                columnkeys.add(key);columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',
                    positive_face_coefficients=[[fi,str(x)]for fi,x in terms],
                    orbit_weights=c['weights'],charge_domain=labels[k],loads=list(map(str,loads))))
        rows.append(dict(vertex=v,union=u,demands=[list(map(str,d))for d in demands],
            compatibility_edges=graph.number_of_edges(),maximal_literal_cliques=ncliques,
            distinct_count_vectors=npatterns,maximal_count_vectors=vectors,literal_representatives=reps,
            certificates=ccs,baseline_mu=[str(m.get(v,F()))for m in mus]))
        counts.update(forests=1,orbits=len(groups),literal_stars=len(owner),
            actual_corners=sum(len(g['origins'])for g in groups),maximal_literal_cliques=ncliques,
            distinct_count_vectors=npatterns,maximal_count_vectors=len(vectors),
            full_mixture_terms=len(ccs[0]['mixture']),restricted_mixture_terms=len(ccs[1]['mixture']),
            unowned_allowed_pairs=u['unowned_allowed_pairs'],inherited_unknown_pairs=u['inherited_unknown_pairs'],
            excluded_pairs=u['excluded_pairs'])
        if step%25==0 or step==len(work):print('FORESTS',step,'/',len(work),'complete',len(rows),
            'open',len(opened),'vectors',counts['maximal_count_vectors'],'raises',
            counts['full_violations'],counts['named_248_class_Q_restriction_violations'],flush=True)
    repairs=[]
    for k,(base,oldm,newm)in enumerate(zip(baseline,mus,newmus)):
        assert all(newm.get(v,F())>=x for v,x in oldm.items())
        objective=F(base['objective_numerator']);normalizer=sum(newm.values(),F())
        repairs.append(dict(domain=labels[k],baseline_bound=base['bound'],objective_numerator=str(objective),
            normalizer=str(normalizer),bound=str(objective/normalizer),
            vertex_updates=[[v,str(oldm.get(v,F())),str(x)]for v,x in sorted(newm.items())if x>oldm.get(v,F())],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(newm.items())if x]))
    result=dict(status='PASS_COMPLETE_CURRENT_HISTORICAL_WEIGHTED_CORNER_CONES'if not opened else'PARTIAL_OPEN_FORESTS',
        sources=sources,source_precedence=order,historical_classes=len(historical),named_geometry_classes=sorted(named),
        selected_faces=[dict(source=owners[i],**selected[i])for i in sorted(selected)],
        additional_positive_faces=sorted(set(selected)-named),face_multipliers=cert['face_multipliers'],
        restricted_face_multipliers=restricted['face_multipliers'],counts=dict(counts),forest_records=rows,
        open_forests=opened,columns=columns,repairs=repairs,
        theorem='Every actual simultaneous corner set is a clique of this conservative graph. Exact nonnegative orbit charges and equal-mass nonnegative covering mixtures certify the entire orbit-weight cone at each stated fixed demand. Complete automorphism orbits and group-averaging justify orbit-count marginals.',
        preservation='All current face masses stay fixed. Only vertex multipliers increase. Every439410 inherited finite-model column, including96176 zero-budget and52 positive-cost signed columns, is preserved. These variables are nonnegative; signed face coefficients do not impose a dual equality. No old-alpha implicit simplex/joint cover is reused.',
        scope='Complete only for the stated saved historical corner orbits and conservative graph, when no forests remain open. First-owner unknowns and never-tested cross-source pairs remain allowed. Covering cliques need not be physically realizable. No global optimizer, all-field coverage, target12 feasibility or dimension theorem.')
    result=json.loads(json.dumps(result,sort_keys=True));out=ROOT/f'results/{TAG}.json'
    if verify:assert result==saved
    else:out.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(result['status'],dict(counts),'bounds',[(r['domain'],float(F(r['bound'])))for r in repairs],
          'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--resume',action='store_true')
    p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    assert sum((a.verify,a.resume,a.prepare_only))<=1
    run(a.verify,a.resume,a.prepare_only)
