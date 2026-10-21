#!/usr/bin/env python3
"""Restricted-alpha optimistic genuinely-unowned deletions in571 current forests.

All full-source and explicit pair owners are immutable. Native geometry is
inherited, not recomputed. Only changed hypothetical local cones may need
new charge proposals; saved replay never invokes an optimizer.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_2599_raw110_raw24_pair_geometry_20261017 as S

T=S.T;N=S.N;U=S.U;O=N.O;H=N.H;W=N.W;ROOT=N.ROOT
TAG='fresh_restricted_290_unowned_pair_frontier_20261017'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.json')


def normal(x):return json.loads(json.dumps(x,sort_keys=True))


def prepare():
    p=N.prepare();n=W.read(N.TAG);t=W.read(T.TAG);s=W.read(S.TAG)
    assert not n['open_forests']
    for key in('face_multipliers','restricted_face_multipliers'):assert n[key]==t[key]==s[key]
    changed={r['vertex']:r for r in n['forest_records']};rows=dict(p['baseline_rows']);rows.update(changed)
    patches=defaultdict(list)
    for tag,d in((T.TAG,t),(S.TAG,s)):
        for r in d['forest_records']:
            patches[r['vertex']].append((tag,r));rows[r['vertex']]=r
    assert len(rows)==571 and len(patches)==9
    assert sum(len(r['newly_owned_pair_domain'])for pp in patches.values()for tag,r in pp)==126
    assert sum(x[2]is None for pp in patches.values()for tag,r in pp for x in r['newly_owned_pair_domain'])==2
    return p,changed,rows,patches,s


def current(v,p,changed,rows,patches):
    previous,g,owner,unowned,dd,nodes,groups,forest,refs,source=T.current(v,p,changed)
    unowned={key:list(pp)for key,pp in unowned.items()};remaining={pair:key for key,pp in unowned.items()for pair in pp}
    allpatch=[]
    for tag,r in patches.get(v,[]):
        assert r['literal_stars']==normal(nodes)and r['corner_orbits']==normal(groups)
        assert r['forest']==normal(forest)and r['demands']==[[str(x)for x in d]for d in dd]
        for i,j,code in r['newly_owned_pair_domain']:
            assert(i,j)in remaining and g.has_edge(i,j),'patch must acquire a previously unowned pair'
            key=remaining.pop((i,j));unowned[key].remove((i,j))
            if code is not None:g.remove_edge(i,j)
            allpatch.append([i,j,code,tag])
    unowned={k:pp for k,pp in unowned.items()if pp}
    assert all(g.has_edge(*pair)for pp in unowned.values()for pair in pp)
    old=rows[v];demand=dd[1]
    if 'demands'in old:assert old['demands'][1]==list(map(str,demand))
    W.check_certificate(old['maximal_count_vectors'],demand,old['certificates'][1])
    for i,x in old['certificates'][1]['mixture']:
        clique=old['literal_representatives'][i]
        assert all(g.has_edge(a,b)for a,b in combinations(clique,2))
        assert old['maximal_count_vectors'][i]==[Counter(owner[z]for z in clique).get(j,0)for j in range(len(groups))]
    return old,g,owner,unowned,demand,groups,forest,allpatch


def positive_cover(old,owner,demand):
    """Remove only zero-demand stars before the hypothetical splitting step.

The charge remains valid on the COMPLETE old vector domain; the projected
cover has the same objective. This avoids spurious splitting costs from
irrelevant zero-demand stars. These cover vectors need not be maximal.
"""
    c=old['certificates'][1];w=[F(x)if demand[j]else F()for j,x in enumerate(c['weights'])]
    assert all(sum(x*y for x,y in zip(v,w))<=1 for v in old['maximal_count_vectors'])
    cover=Counter()
    for vi,x in c['mixture']:
        clique=tuple(i for i in old['literal_representatives'][vi]if demand[owner[i]]>0)
        assert clique,'an optimal positive-mass cover cannot have an entirely zero-demand term'
        cover[clique]+=F(x)
    cliques=sorted(cover);vectors=[[Counter(owner[i]for i in cl).get(j,0)for j in range(len(demand))]for cl in cliques]
    cert=dict(weights=list(map(str,w)),mixture=[[i,str(cover[cl])]for i,cl in enumerate(cliques)],optimal_mass=c['optimal_mass'])
    W.check_certificate(vectors,demand,cert)
    return dict(certificates=[cert],maximal_count_vectors=vectors,literal_representatives=list(map(list,cliques))),dict(
        weights=cert['weights'],cliques=list(map(list,cliques)),count_vectors=vectors,mixture=cert['mixture'],optimal_mass=cert['optimal_mass'])


def bindings():
    sources={}
    for tag in(S.TAG,T.TAG,N.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(S,T,N,U,O,H,W):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    return sources


def run(verify=False,resume=False,prepare_only=False):
    start=time.monotonic();p,changed,rows,patches,base=prepare();sources=bindings()
    alpha=base['restricted_face_multipliers'];baseline=base['repairs'][1]
    mu={i:F(x)for i,x in baseline['vertex_multipliers']};binding=dict(sources=sources,alpha=alpha,baseline_repair=baseline)
    saved=W.read(TAG)if verify else None;cache={}
    if saved:
        assert saved['sources']==sources;cache={r['vertex']:r for r in saved['forest_records']}
    elif resume and PROGRESS.exists():
        c=json.loads(PROGRESS.read_text());assert c['binding']==normal(binding);cache={r['vertex']:r for r in c['forest_records']}
    records=[];counts=Counter();lo_mu=dict(mu);hi_mu=dict(mu)
    for step,v in enumerate(sorted(rows),1):
        if v in cache and not verify and not prepare_only:rec=cache[v]
        else:
            old,g,owner,unowned,demand,groups,forest,applied=current(v,p,changed,rows,patches)
            shadow,projected=positive_cover(old,owner,demand)
            pairs=sorted(pair for pp in unowned.values()for pair in pp);positive=[pair for pair in pairs if all(demand[owner[i]]>0 for i in pair)]
            c=dict(unowned_pairs=len(pairs),orbit_pair_domains=len(unowned),positive_demand_pairs=len(positive),
                positive_demand_domains=sum(demand[k[0]]>0 and demand[k[1]]>0 for k in unowned),
                applied_recent_tests=len(applied),applied_recent_unknowns=sum(x[2]is None for x in applied))
            rec=dict(vertex=v,forest=forest,corner_orbits=groups,demands=list(map(str,demand)),
                recent_patch_domains=applied,projected_baseline_cover=projected,counts=c)
            if not prepare_only:
                prior=cache.get(v)
                grand=O.scenario(v,g,owner,shadow,demand,pairs,mu.get(v,F()),prior['all_unowned_deletion']if prior else None)
                candidates=[];no_gain=F(grand['required_mu_increase_upper_bound'])==0
                savedc={tuple(x['corner_orbits']):x for x in prior['candidate_domains']}if prior else {}
                for key,pp in sorted(unowned.items()):
                    if no_gain or not(demand[key[0]]and demand[key[1]]):
                        candidates.append(dict(corner_orbits=list(key),complete_domain_pairs=len(pp),pair_domain_sha256=O.pairhash(pp),
                            method='verified_zero_gain',required_mu_increase_lower_bound='0',required_mu_increase_upper_bound='0',
                            reason='all_unowned_cover'if no_gain else 'one_endpoint_has_zero_demand'))
                    else:
                        sc=savedc[key]['scenario']if prior else None
                        result=O.scenario(v,g,owner,shadow,demand,pp,mu.get(v,F()),sc,key[0],grand)
                        candidates.append(dict(corner_orbits=list(key),complete_domain_pairs=len(pp),
                            face_coefficients=[groups[j]['positive_face_coefficients']for j in key],
                            raw_families=[sorted({x['raw_family']for x in groups[j]['origins']})for j in key],
                            endpoint_demands=[str(demand[j])for j in key],scenario=result))
                rec.update(all_unowned_deletion=grand,candidate_domains=candidates)
                if prior:assert normal(rec)==prior
        records.append(rec);counts.update(rec['counts']);counts['forests']+=1
        if not prepare_only:
            q=rec['all_unowned_deletion'];lo_mu[v]=max(mu.get(v,F()),F(q['optimal_mass_lower_bound']));hi_mu[v]=max(mu.get(v,F()),F(q['optimal_mass_upper_bound']))
            counts.update(exact_grand_cones=int(q['method']=='exact_charge_and_cover'),
                zero_gain_forests=int(F(q['required_mu_increase_upper_bound'])==0),
                positive_grand_forests=int(F(q['required_mu_increase_lower_bound'])>0))
            if not verify:PROGRESS.write_text(json.dumps(normal(dict(binding=binding,forest_records=records)),separators=(',',':'))+'\n')
        if step%40==0 or rec['counts']['positive_demand_pairs']:
            print('FOREST',step,'/571',v,rec['counts'],
                'grand_delta',rec.get('all_unowned_deletion',{}).get('required_mu_increase_upper_bound','not computed'),flush=True)
    assert counts['unowned_pairs']==570524-126==570398
    assert counts['forests']==571 and counts['applied_recent_tests']==126 and counts['applied_recent_unknowns']==2
    if prepare_only:
        path=Path('/private/tmp/restricted290_unowned_domain.json');path.write_text(json.dumps(normal(dict(binding=binding,counts=dict(counts),forest_records=records)),separators=(',',':'))+'\n')
        print('PASS_RESTRICTED_290_CURRENT_DOMAIN',dict(counts),'seconds',round(time.monotonic()-start,2),flush=True);return
    numerator=F(baseline['objective_numerator']);den=F(baseline['normalizer']);assert den==sum(mu.values(),F())
    ranks=[]
    for r in records:
        for c in r['candidate_domains']:
            q=c.get('scenario',c);lo=F(q['required_mu_increase_lower_bound']);hi=F(q['required_mu_increase_upper_bound'])
            if hi:ranks.append(dict(vertex=r['vertex'],corner_orbits=c['corner_orbits'],pairs=c['complete_domain_pairs'],
                increase_lower=str(lo),increase_upper=str(hi),hypothetical_ratio_lower=str(numerator/(den+hi)),hypothetical_ratio_upper=str(numerator/(den+lo)),
                certified_exact=lo==hi,raw_families=c.get('raw_families'),face_coefficients=c.get('face_coefficients')))
    ranks.sort(key=lambda r:(-F(r['increase_lower']),-F(r['increase_upper']),r['pairs'],r['vertex'],r['corner_orbits']))
    need=numerator/12-den;assert need>0
    choice=U.shortlist(ranks,need)
    if choice['status'].startswith('PASS'):
        choice['hypothetical_ratio_upper']=str(numerator/(den+F(choice['certified_increase_lower'])))
    out=normal(dict(status='PASS_RESTRICTED_290_CURRENT_UNOWNED_PAIR_FRONTIER',sources=sources,baseline_source=S.TAG,
        face_multipliers=alpha,baseline_restricted_repair=baseline,forest_records=records,counts=dict(counts),ranked_candidate_domains=ranks,shortlist=choice,
        optimistic_all_unowned=dict(required_normalizer_lower=str(sum(lo_mu.values(),F())),required_normalizer_upper=str(sum(hi_mu.values(),F())),
            hypothetical_ratio_lower=str(numerator/sum(hi_mu.values(),F())),hypothetical_ratio_upper=str(numerator/sum(lo_mu.values(),F()))),
        actual_witness_unchanged=True,actual_new_incompatibilities=0,actual_exported_columns=[],
        scope='Current571-forest conservative graph at restricted alpha. All historical owners, raw125/raw205 tests and126 recent explicit outcomes including2 unknowns remain fixed. Only genuinely unowned pairs may be deleted hypothetically. Positive-coordinate projection avoids spurious zero-demand splitting costs. No physical pair/covector reconstruction, global optimizer, new actual inequality, physical feasibility or dimension theorem. Saved replay uses rational certificates without an optimizer.'))
    path=ROOT/f'results/{TAG}.json'
    if saved:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'optimistic',out['optimistic_all_unowned'],'shortlist',{k:v for k,v in choice.items()if k not in('dynamic_tables','chosen')},'seconds',round(time.monotonic()-start,2),flush=True)
    print('TOP',ranks[:10],flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');a.add_argument('--resume',action='store_true');a.add_argument('--prepare',action='store_true');z=a.parse_args()
    assert sum((z.verify,z.resume,z.prepare))<=1;run(z.verify,z.resume,z.prepare)
