#!/usr/bin/env python3
"""Updated counterfactual unowned-pair gains at the frozen290-class witness.

All new raw125/raw205 owned outcomes, including unknowns, remain fixed.
Only96 changed forests are recomputed;475 unchanged hypothetical certificates
are inherited by hash. No physical exclusion or global optimizer is used.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_raw125_raw205_actual_corner_frontier_20261017 as N

O=N.O;W=N.W;H=N.H;ROOT=N.ROOT
TAG='fresh_raw125_raw205_unowned_pair_frontier_20261017'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.json')


def normal(x):return json.loads(json.dumps(x,sort_keys=True))
def datahash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def current(v,row,p):
    oldcount=row['old_literal_stars'];n=len(row['literal_stars']);groups=row['corner_orbits']
    if oldcount:
        u,g,oldowner,unowned,demand=O.current(v,p,p['expanded_domains'],p['baseline_rows'],p['later'])
        assert normal(u)==row['inherited_union']
        oldnodes=[normal(W.forests(s)[v]['literal_stars'][i])for s,i in u['literal_references']]
        assert oldnodes==row['literal_stars'][:oldcount]
        assert normal(u['corner_orbits'])==groups[:len(u['corner_orbits'])]
    else:
        assert v not in p['baseline_rows']and row['inherited_union']is None
        g=nx.Graph();unowned={};oldowner={};demand=[]
    owner={i:j for j,z in enumerate(groups)for i in z['literal_stars']}
    assert set(owner)==set(range(n))and all(owner[i]==j for i,j in oldowner.items())
    oldedges={tuple(sorted(e))for e in g.edges()};g.add_nodes_from(range(n))
    newpairs={(i,j):c for i,j,c in row['newly_owned_pair_domain']}
    assert len(newpairs)==len(row['newly_owned_pair_domain'])
    assert set(newpairs)=={(i,j)for i,j in combinations(range(n),2)if j>=oldcount}
    for(i,j),code in newpairs.items():
        if code is None:g.add_edge(i,j)
    assert {tuple(sorted(e))for e in g.edges()if max(e)<oldcount}==oldedges
    assert all(j<oldcount and(i,j)not in newpairs and g.has_edge(i,j)for pairs in unowned.values()for i,j in pairs)
    newdemand=[sum((p['alpha'].get(fi,F())*F(m)for fi,m in z['positive_face_coefficients']),F())for z in groups]
    assert all(x>0 for x in newdemand)and list(map(str,newdemand))==row['demands'][0]
    assert newdemand[:len(demand)]==demand
    W.check_certificate(row['maximal_count_vectors'],newdemand,row['certificates'][0])
    for vi,x in row['certificates'][0]['mixture']:
        cl=row['literal_representatives'][vi]
        assert all(g.has_edge(i,j)for i,j in combinations(cl,2))
        assert row['maximal_count_vectors'][vi]==[Counter(owner[i]for i in cl).get(j,0)for j in range(len(groups))]
    return g,owner,unowned,newdemand


def shortlist(ranks,need):
    # A certificate-backed lower gain is additive across distinct forests.
    # This DP minimizes pair count within this restricted selection rule,
    # not over arbitrary simultaneous deletions or physical compatibility.
    positive=[r for r in ranks if F(r['increase_lower'])>0]
    by={v:[r for r in positive if r['vertex']==v]for v in sorted({r['vertex']for r in positive})}
    order=sorted(positive,key=lambda r:(F(r['increase_lower'])/r['pairs'],F(r['increase_lower'])),reverse=True)
    used=set();seed=[];gain=F()
    for r in order:
        if r['vertex']in used:continue
        used.add(r['vertex']);seed.append(r);gain+=F(r['increase_lower'])
        if gain>=need:break
    if gain<need:return dict(status='no_certified_disjoint_forest_combination',required_increase=str(need))
    cap=sum(r['pairs']for r in seed);states={0:(F(),[])};tables=[]
    for v,choices in sorted(by.items()):
        nxt=dict(states)
        for cost,(old,items)in states.items():
            for r in choices:
                c=cost+r['pairs']
                if c>cap:continue
                z=old+F(r['increase_lower']);ids=items+[[r['vertex'],r['corner_orbits']]]
                if c not in nxt or z>nxt[c][0]:nxt[c]=(z,ids)
        # Delete dominated states; a cheaper state with at least as much
        # certified gain remains no worse under every future forest choice.
        best=F(-1);states={}
        for cost,val in sorted(nxt.items()):
            if val[0]>best:states[cost]=val;best=val[0]
        feasible=[c for c,(z,ids)in states.items()if z>=need]
        if feasible:cap=min(cap,min(feasible));states={c:z for c,z in states.items()if c<=cap}
        tables.append(dict(vertex=v,states=[[c,str(z),ids]for c,(z,ids)in sorted(states.items())]))
    cost=min(c for c,(z,ids)in states.items()if z>=need);value,ids=states[cost]
    chosen=[next(r for r in positive if r['vertex']==v and r['corner_orbits']==pair)for v,pair in ids]
    assert len({r['vertex']for r in chosen})==len(chosen)
    assert sum(r['pairs']for r in chosen)==cost and sum((F(r['increase_lower'])for r in chosen),F())==value
    return dict(status='PASS_EXACT_DISJOINT_FOREST_SHORTLIST',required_increase=str(need),
        complete_pair_count=cost,certified_increase_lower=str(value),chosen=chosen,
        dynamic_tables=tables,scope='Minimum pair count among saved positive lower-gain single-domain scenarios, at most one domain from each forest. No claim of globally smallest physical test or arbitrary-union optimum.')


def run(verify=False,resume=False):
    start=time.monotonic();p=N.prepare();base=W.read(N.TAG);prior=W.read(O.TAG)
    assert not base['open_forests']and len(base['forest_records'])==96
    assert base['face_multipliers']==prior['face_multipliers']==p['baseline']['face_multipliers']
    sources={}
    for tag in(N.TAG,O.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(N,O,H,W):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    binding=dict(sources=sources,baseline_repair=base['repairs'][0])
    saved=W.read(TAG)if verify else None;cache={}
    if saved is not None:
        assert saved['sources']==sources;cache={r['vertex']:r for r in saved['forest_records']}
    elif resume and PROGRESS.exists():
        x=json.loads(PROGRESS.read_text());assert x['binding']==normal(binding)
        cache={r['vertex']:r for r in x['forest_records']}
    mu={i:F(x)for i,x in base['repairs'][0]['vertex_multipliers']};oldmu={i:F(x)for i,x in prior['baseline_full_repair']['vertex_multipliers']}
    changed={r['vertex']:r for r in base['forest_records']};pr={r['vertex']:r for r in prior['forest_records']}
    inherited=[];allrecords={};records=[]
    for v,r in sorted(pr.items()):
        if v in changed:continue
        assert mu.get(v,F())==oldmu.get(v,F())
        inherited.append(dict(vertex=v,source=O.TAG,source_record_sha256=datahash(r)))
        allrecords[v]=r
    assert len(inherited)==475
    for step,(v,old)in enumerate(sorted(changed.items()),1):
        if v in cache and not verify:rec=cache[v]
        else:
            g,owner,unowned,demand=current(v,old,p);pairs=sorted(pair for pp in unowned.values()for pair in pp)
            previous=cache.get(v)
            grand=O.scenario(v,g,owner,old,demand,pairs,mu.get(v,F()),previous['all_unowned_deletion']if previous else None)
            candidates=[];no_gain=F(grand['required_mu_increase_upper_bound'])==0
            oldc={tuple(c['corner_orbits']):c for c in previous['candidate_domains']}if previous else {}
            for key,pp in sorted(unowned.items()):
                if no_gain:
                    candidates.append(dict(corner_orbits=list(key),complete_domain_pairs=len(pp),pair_domain_sha256=O.pairhash(pp),
                        method='all_unowned_cover_proves_zero_witness_gain',required_mu_increase_lower_bound='0',required_mu_increase_upper_bound='0'))
                else:
                    cert=oldc[key]['scenario']if previous else None
                    result=O.scenario(v,g,owner,old,demand,pp,mu.get(v,F()),cert,key[0],grand)
                    candidates.append(dict(corner_orbits=list(key),complete_domain_pairs=len(pp),
                        face_coefficients=[old['corner_orbits'][j]['positive_face_coefficients']for j in key],
                        raw_families=[sorted({x['raw_family']for x in old['corner_orbits'][j]['origins']})for j in key],
                        endpoint_demands=[str(demand[j])for j in key],scenario=result))
            rec=dict(vertex=v,geometry_source=N.TAG,forest=old['forest'],corner_orbits=old['corner_orbits'],
                new_owned_domain_sha256=datahash(old['newly_owned_pair_domain']),demands=list(map(str,demand)),
                complete_unowned_pairs=len(pairs),all_unowned_deletion=grand,candidate_domains=candidates)
            if previous:assert normal(rec)==previous
        records.append(rec);allrecords[v]=rec
        if not verify:PROGRESS.write_text(json.dumps(normal(dict(binding=binding,forest_records=records)),separators=(',',':'))+'\n')
        grand=rec['all_unowned_deletion']
        print('FOREST',step,'/96',v,'pairs',rec['complete_unowned_pairs'],'delta',grand['required_mu_increase_lower_bound'],grand['required_mu_increase_upper_bound'],flush=True)
    assert len(allrecords)==571
    counts=Counter();low=dict(mu);high=dict(mu);ranks=[]
    objective=F(base['repairs'][0]['objective_numerator']);den=F(base['repairs'][0]['normalizer']);assert den==sum(mu.values(),F())
    for v,r in sorted(allrecords.items()):
        q=r['all_unowned_deletion'];low[v]=max(mu.get(v,F()),F(q['optimal_mass_lower_bound']));high[v]=max(mu.get(v,F()),F(q['optimal_mass_upper_bound']))
        counts.update(forests=1,unowned_pairs=r['complete_unowned_pairs'],candidate_domains=len(r['candidate_domains']),
            positive_grand_forests=int(F(q['required_mu_increase_lower_bound'])>0),
            zero_gain_forests=int(F(q['required_mu_increase_upper_bound'])==0),
            exact_grand_cones=int(q['method']=='exact_charge_and_cover'))
        for c in r['candidate_domains']:
            s=c.get('scenario',c);lo=F(s['required_mu_increase_lower_bound']);hi=F(s['required_mu_increase_upper_bound'])
            if hi:ranks.append(dict(vertex=v,corner_orbits=c['corner_orbits'],pairs=c['complete_domain_pairs'],
                increase_lower=str(lo),increase_upper=str(hi),hypothetical_ratio_lower=str(objective/(den+hi)),
                hypothetical_ratio_upper=str(objective/(den+lo)),certified_exact=lo==hi,
                raw_families=c.get('raw_families'),face_coefficients=c.get('face_coefficients')))
    assert counts['unowned_pairs']==prior['counts']['unowned_pairs']==570524
    ranks.sort(key=lambda r:(-F(r['increase_lower']),-F(r['increase_upper']),r['pairs'],r['vertex'],r['corner_orbits']))
    need=objective/12-den;assert need>0
    chosen=shortlist(ranks,need)
    if chosen['status'].startswith('PASS'):
        chosen['hypothetical_ratio_upper']=str(objective/(den+F(chosen['certified_increase_lower'])))
        assert F(chosen['hypothetical_ratio_upper'])<=12
    out=normal(dict(status='PASS_UPDATED_290_CLASS_UNOWNED_PAIR_DIAGNOSTIC',sources=sources,baseline_source=N.TAG,
        baseline_full_repair=base['repairs'][0],face_multipliers=base['face_multipliers'],
        inherited_forests=inherited,forest_records=records,counts=dict(counts),ranked_candidate_domains=ranks,shortlist=chosen,
        optimistic_all_unowned=dict(required_normalizer_lower=str(sum(low.values(),F())),required_normalizer_upper=str(sum(high.values(),F())),
            hypothetical_ratio_lower=str(objective/sum(high.values(),F())),hypothetical_ratio_upper=str(objective/sum(low.values(),F()))),
        actual_new_incompatibilities=0,actual_exported_columns=[],actual_witness_unchanged=True,
        scope='Updated fixed290-class witness only. All571 forests covered,475 same-demand hypothetical certificates inherited and96 changed domains rebuilt. Every new raw125/raw205 tested outcome, including unknowns, remains owned. Only genuinely unowned old pairs may be deleted counterfactually; no physical exclusion, actual repaired bound, global optimum or dimension theorem. Saved replay never invokes an optimizer.'))
    path=ROOT/f'results/{TAG}.json'
    if verify:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'all_unowned',out['optimistic_all_unowned'],'shortlist',{k:v for k,v in chosen.items()if k not in('dynamic_tables','chosen')},'seconds',round(time.monotonic()-start,2),flush=True)
    print('CHOSEN',chosen.get('chosen'),'TOP',ranks[:8],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--resume',action='store_true');a=ap.parse_args()
    assert not(a.verify and a.resume);run(a.verify,a.resume)
