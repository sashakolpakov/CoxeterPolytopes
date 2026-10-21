#!/usr/bin/env python3
"""Counterfactual complete deletion of genuinely unowned compatibility pairs.

No physical incompatibility is asserted and no actual witness is repaired.
Generation may propose small local cone certificates; fresh replay uses saved
rationals. Existing owned outcomes, including unknowns, are immutable.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_historical22_corner_domain_20261017 as D
import fresh_historical22_corner_cones_20261017 as C

ROOT=D.ROOT;H=D.H;W=D.W
TAG='fresh_unowned_pair_optimistic_frontier_20261017'
BASE='fresh_historical22_raw88_raw25_pair_geometry_20261017'
EIGHT='fresh_historical22_eight_pair_geometry_20261017'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.json')


def normal(x):return json.loads(json.dumps(x,sort_keys=True))
def pairhash(pairs):return sha256(json.dumps(pairs,separators=(',',':')).encode()).hexdigest()


def current(v,p,domains,rows,later):
    if v in domains:
        u,g,owner=D.assemble_forest(v,p);assert normal(u)==domains[v]
        ds=u['source_domains']
        tr={j:i for i,j in u['explicit_patch']['new_to_patch_star_indices']}
    else:
        u,g,owner=H.union_forest(v,p['order'],p['selected'])
        assert normal(u)==p['original_rows'][v]['union']
        ds=u['first_owner_domains'];tr={i:i for i in range(len(owner))}
    bits=[0]*len(owner)
    for k,d in enumerate(ds):
        for i,j in d['new_to_old_star_indices']:bits[i]|=1<<k
    explicit={}
    if v in p['patch_rows']:
        for i,j,code in p['patch_rows'][v]['newly_owned_pair_domain']:
            a,b=sorted((tr[i],tr[j]));assert not(bits[a]&bits[b])
            explicit[a,b]=code
            if v not in domains and code is not None:g.remove_edge(a,b)
            assert g.has_edge(a,b)==(code is None)
    for patch in later.get(v,[]):
        assert patch['literal_references']==normal(u['literal_references'])
        for i,j,code in patch['newly_owned_pair_domain']:
            assert i<j and not(bits[i]&bits[j])and(i,j)not in explicit
            assert g.has_edge(i,j)
            explicit[i,j]=code
            if code is not None:g.remove_edge(i,j)
    unowned=defaultdict(list)
    for i,j in combinations(range(len(owner)),2):
        if not(bits[i]&bits[j])and(i,j)not in explicit:
            assert g.has_edge(i,j)
            key=tuple(sorted((owner[i],owner[j])));assert key[0]!=key[1]
            unowned[key].append((i,j))
    old=rows[v];groups=u['corner_orbits']
    demand=[sum((p['alpha'].get(fi,F())*F(m)for fi,m in z['positive_face_coefficients']),F())for z in groups]
    assert all(x>0 for x in demand)
    if 'demands'in old:assert old['demands'][0]==list(map(str,demand))
    assert len(old['certificates'][0]['weights'])==len(groups)
    assert W.check_certificate(old['maximal_count_vectors'],demand,old['certificates'][0])==F(old['certificates'][0]['optimal_mass'])
    for vi,amount in old['certificates'][0]['mixture']:
        clique=old['literal_representatives'][vi]
        assert all(g.has_edge(i,j)for i,j in combinations(clique,2))
        assert old['maximal_count_vectors'][vi]==[Counter(owner[i]for i in clique).get(j,0)for j in range(len(groups))]
    return u,g,owner,unowned,demand


def split_cover(graph,owner,old,demand,selected_orbit=None):
    """Explicit feasible cover after deletions; never assumes realizability."""
    cover=Counter();affected=F()
    for vi,amount in old['certificates'][0]['mixture']:
        mass=F(amount);clique=old['literal_representatives'][vi]
        if all(graph.has_edge(i,j)for i,j in combinations(clique,2)):parts=[clique]
        elif selected_orbit is not None:
            # Deleted edges run between exactly two distinct corner orbits.
            parts=[[i for i in clique if owner[i]==selected_orbit],
                   [i for i in clique if owner[i]!=selected_orbit]]
            assert all(parts);affected+=mass
        else:
            parts=[];affected+=mass
            for i in sorted(clique,key=lambda i:(-graph.degree(i),i)):
                target=next((part for part in parts if all(graph.has_edge(i,j)for j in part)),None)
                if target is None:parts.append([i])
                else:target.append(i)
        for part in parts:
            key=tuple(sorted(part));assert key and all(graph.has_edge(i,j)for i,j in combinations(key,2))
            cover[key]+=mass
    marginal=[F()for _ in demand]
    for clique,mass in cover.items():
        for i in clique:marginal[owner[i]]+=mass
    assert all(x>=y for x,y in zip(marginal,demand))
    terms=[[list(clique),str(x)]for clique,x in sorted(cover.items())if x]
    value=sum(cover.values(),F());assert value>=F(old['certificates'][0]['optimal_mass'])
    return dict(terms=terms,total_mass=str(value),baseline_cover_mass_using_deleted_edges=str(affected))


def scenario(v,graph,owner,old,demand,pairs,mu,saved=None,orbit=None,upper_cap=None):
    gg=graph.copy();assert all(gg.has_edge(i,j)for i,j in pairs);gg.remove_edges_from(pairs)
    lower=F(old['certificates'][0]['optimal_mass'])
    cover=split_cover(gg,owner,old,demand,orbit)
    upper=F(cover['total_mass']);assert upper>=lower
    out=dict(deleted_pairs=len(pairs),pair_domain_sha256=pairhash(pairs),
        baseline_optimal_mass=str(lower),baseline_mu=str(mu),split_cover=cover)
    if upper_cap is not None:
        # A cover for the graph with ALL unowned edges deleted is valid here.
        assert F(upper_cap['optimal_mass_upper_bound'])>=lower
        upper=min(upper,F(upper_cap['optimal_mass_upper_bound']))
    if upper<=mu or upper==lower:
        out.update(method='verified_cover_upper_bound',optimal_mass_lower_bound=str(lower),optimal_mass_upper_bound=str(upper))
    elif saved is not None and saved['method']=='verified_cover_upper_bound_after_guard':
        out.update(method=saved['method'],optimal_mass_lower_bound=str(lower),optimal_mass_upper_bound=str(upper),
                   limitation='Exact complete clique enumeration or certificate proposal was not completed; explicit cover remains a rigorous upper bound.')
    else:
        enum,err=C.vectors_for(v,gg,owner,len(demand))
        if err:
            assert saved is None
            out.update(method='verified_cover_upper_bound_after_guard',optimal_mass_lower_bound=str(lower),optimal_mass_upper_bound=str(upper),
                       limitation='Exact complete clique enumeration or certificate proposal was not completed; explicit cover remains a rigorous upper bound.')
        else:
            vectors,reps,ncliques,npatterns=enum
            try:
                cert=saved['certificate']if saved is not None else H.exact_profile(vectors,demand)
                value=W.check_certificate(vectors,demand,cert)
            except Exception:
                if saved is not None:raise
                out.update(method='verified_cover_upper_bound_after_guard',optimal_mass_lower_bound=str(lower),optimal_mass_upper_bound=str(upper),
                           limitation='Exact complete clique enumeration or certificate proposal was not completed; explicit cover remains a rigorous upper bound.')
            else:
                assert lower<=value<=upper
                for vec,clique in zip(vectors,reps):
                    assert all(gg.has_edge(i,j)for i,j in combinations(clique,2))
                    assert list(vec)==[Counter(owner[i]for i in clique).get(j,0)for j in range(len(demand))]
                out.update(method='exact_charge_and_cover',optimal_mass_lower_bound=str(value),optimal_mass_upper_bound=str(value),
                    maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,
                    maximal_count_vectors=list(map(list,vectors)),literal_representatives=list(map(list,reps)),certificate=cert)
    lo,hi=map(F,(out['optimal_mass_lower_bound'],out['optimal_mass_upper_bound']))
    out.update(required_mu_increase_lower_bound=str(max(F(),lo-mu)),required_mu_increase_upper_bound=str(max(F(),hi-mu)))
    if saved is not None:assert normal(out)==saved,(v,orbit,'saved hypothetical certificate')
    return out


def run(verify=False,resume=False):
    start=time.monotonic();p=D.prepare();domain=W.read(D.TAG);cones=W.read(C.TAG);eight=W.read(EIGHT);base=W.read(BASE)
    assert not cones['open_forests']
    for key in('face_multipliers','restricted_face_multipliers'):
        assert base[key]==eight[key]==cones[key]==domain[key]
    rows=dict(p['original_rows']);rows.update(p['patch_rows']);rows.update({r['vertex']:r for r in cones['forest_records']})
    rows[eight['vertex']]=eight;rows.update({r['vertex']:r for r in base['forest_records']});assert len(rows)==553
    domains={r['vertex']:r for r in domain['forest_records']}
    later=defaultdict(list);later[eight['vertex']].append(eight)
    for r in base['forest_records']:later[r['vertex']].append(r)
    mu={i:F(x)for i,x in base['repairs'][0]['vertex_multipliers']}
    sources={}
    for tag in(D.BASE,D.PATCH,D.TAG,C.TAG,EIGHT,BASE):
        d=W.read(tag)
        for path,h in d['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for pp in(Path(__file__),Path(H.__file__),Path(W.__file__)):
        path=str(pp.relative_to(ROOT));sources[path]=W.digest(path)
    binding=dict(sources=sources,face_multipliers=base['face_multipliers'])
    saved=W.read(TAG)if verify else None;cached={}
    if verify:
        assert saved['sources']==sources;cached={r['vertex']:r for r in saved['forest_records']}
    elif resume and PROGRESS.exists():
        progress=json.loads(PROGRESS.read_text());assert progress['binding']==binding
        cached={r['vertex']:r for r in progress['forest_records']}
    records=[];counts=Counter();grand_mu_low=dict(mu);grand_mu_high=dict(mu)
    for step,v in enumerate(sorted(rows),1):
        u,g,owner,unowned,demand=current(v,p,domains,rows,later)
        old=rows[v];assert F(old['certificates'][0]['optimal_mass'])<=mu.get(v,F())
        allpairs=sorted(pair for pairs in unowned.values()for pair in pairs)
        assert len(allpairs)==len(set(allpairs))
        prior=cached.get(v)
        grand=scenario(v,g,owner,old,demand,allpairs,mu.get(v,F()),prior['all_unowned_deletion']if prior else None)
        candidates=[];all_no_gain=F(grand['required_mu_increase_upper_bound'])==0
        if all_no_gain:
            for key,pairs in sorted(unowned.items()):
                candidates.append(dict(corner_orbits=list(key),complete_domain_pairs=len(pairs),pair_domain_sha256=pairhash(pairs),
                    method='all_unowned_cover_proves_zero_witness_gain',required_mu_increase_lower_bound='0',required_mu_increase_upper_bound='0'))
        else:
            saved_candidates={tuple(r['corner_orbits']):r for r in prior['candidate_domains']}if prior else {}
            for key,pairs in sorted(unowned.items()):
                oldscenario=saved_candidates[key]['scenario']if prior else None
                result=scenario(v,g,owner,old,demand,pairs,mu.get(v,F()),oldscenario,key[0],grand)
                candidates.append(dict(corner_orbits=list(key),complete_domain_pairs=len(pairs),
                    face_coefficients=[u['corner_orbits'][j]['positive_face_coefficients']for j in key],
                    raw_families=[sorted({o['raw_family']for o in u['corner_orbits'][j]['origins']})for j in key],
                    endpoint_demands=[str(demand[j])for j in key],scenario=result))
        rec=dict(vertex=v,forest=u['forest'],literal_references=u['literal_references'],corner_orbits=u['corner_orbits'],
            baseline_certificate_source=BASE if v in later and v!=eight['vertex']else EIGHT if v==eight['vertex']else C.TAG if v in domains else D.PATCH if v in p['patch_rows']else D.BASE,
            demands=list(map(str,demand)),complete_unowned_pairs=len(allpairs),all_unowned_deletion=grand,candidate_domains=candidates)
        if prior:assert normal(rec)==prior
        records.append(rec)
        grand_mu_low[v]=max(mu.get(v,F()),F(grand['optimal_mass_lower_bound']))
        grand_mu_high[v]=max(mu.get(v,F()),F(grand['optimal_mass_upper_bound']))
        counts.update(forests=1,unowned_pairs=len(allpairs),candidate_domains=len(candidates),
            zero_gain_forests=int(all_no_gain),potential_gain_forests=int(not all_no_gain),
            exact_grand_cones=int(grand['method']=='exact_charge_and_cover'),
            individual_scenarios=sum('scenario'in c for c in candidates),
            exact_individual_cones=sum(c.get('scenario',{}).get('method')=='exact_charge_and_cover'for c in candidates),
            positive_individual_gain=sum(F(c.get('scenario',c)['required_mu_increase_lower_bound'])>0 for c in candidates))
        if not verify:PROGRESS.write_text(json.dumps(dict(binding=binding,forest_records=normal(records)),sort_keys=True,separators=(',',':'))+'\n')
        if not all_no_gain or step%30==0:
            print('FOREST',step,'/553',v,'unowned',len(allpairs),'domains',len(candidates),'grand_delta',grand['required_mu_increase_lower_bound'],grand['required_mu_increase_upper_bound'],
                  'single_positive',sum(F(c.get('scenario',c)['required_mu_increase_lower_bound'])>0 for c in candidates),flush=True)
    objective=F(base['repairs'][0]['objective_numerator']);den=F(base['repairs'][0]['normalizer']);assert den==sum(mu.values(),F())
    ranks=[]
    for r in records:
        for c in r['candidate_domains']:
            q=c.get('scenario',c);lo=F(q['required_mu_increase_lower_bound']);hi=F(q['required_mu_increase_upper_bound'])
            if hi:
                ranks.append(dict(vertex=r['vertex'],corner_orbits=c['corner_orbits'],pairs=c['complete_domain_pairs'],
                    increase_lower=str(lo),increase_upper=str(hi),hypothetical_ratio_lower=str(objective/(den+hi)),
                    hypothetical_ratio_upper=str(objective/(den+lo)),certified_exact=lo==hi))
    ranks.sort(key=lambda x:(-F(x['increase_lower']),-F(x['increase_upper']),x['pairs'],x['vertex'],x['corner_orbits']))
    out=normal(dict(status='PASS_COUNTERFACTUAL_COMPLETE_UNOWNED_PAIR_DELETION_FRONTIER',sources=sources,
        baseline_source=BASE,counts=dict(counts),forest_records=records,ranked_candidate_domains=ranks,
        face_multipliers=base['face_multipliers'],baseline_full_repair=base['repairs'][0],
        optimistic_all_unowned=dict(required_normalizer_lower=str(sum(grand_mu_low.values(),F())),
            required_normalizer_upper=str(sum(grand_mu_high.values(),F())),
            hypothetical_ratio_lower=str(objective/sum(grand_mu_high.values(),F())),
            hypothetical_ratio_upper=str(objective/sum(grand_mu_low.values(),F()))),
        actual_new_incompatibilities=0,actual_exported_columns=[],actual_witness_unchanged=True,
        scope='Counterfactual deletion only. Every current553 forest and all complete genuinely unowned orbit-pair domains are covered. Earlier owned outcomes including explicit unknowns are immutable. Split covers give rigorous improvement ceilings; complete graph charges/covers give exact hypothetical local optima where computed. No physical exclusion, actual witness repair, global model solve or dimension theorem. Fresh replay never invokes an optimizer.'))
    path=ROOT/f'results/{TAG}.json'
    if verify:assert out==saved
    else:path.write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'optimistic',out['optimistic_all_unowned'],'seconds',round(time.monotonic()-start,2),flush=True)
    print('TOP',ranks[:12],flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');a.add_argument('--resume',action='store_true')
    z=a.parse_args();assert not(z.verify and z.resume);run(z.verify,z.resume)
