#!/usr/bin/env python3
"""Actual exact tests of the frozen118-pair shortlist, with local certificates.

Only genuinely unowned pairs may change. Unknowns remain compatible.
No global solve; fresh replay uses saved rational local covers.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_raw125_raw205_unowned_pair_frontier_20261017 as U
import fresh_historical22_raw88_raw25_pair_geometry_20261017 as B

N=U.N;O=N.O;H=N.H;W=N.W;R=N.R;ROOT=N.ROOT
TAG='fresh_290_class_shortlist_pair_geometry_20261017'
TARGETS=((788,1,22,72),(1505,15,18,16),(1683,0,11,8),(1763,3,8,4),
         (2622,5,14,4),(2948,3,4,4),(3324,6,7,2),(4621,2,3,8))


def normal(x):return json.loads(json.dumps(x,sort_keys=True))


def current(v,p,changed):
    if v in changed:
        old=changed[v];g,owner,unowned,demand=U.current(v,old,p)
        nodes=old['literal_stars'];groups=old['corner_orbits'];forest=old['forest'];refs=None
        source=N.TAG
    else:
        u,g,owner,unowned,demand=O.current(v,p,p['expanded_domains'],p['baseline_rows'],p['later'])
        old=p['baseline_rows'][v];nodes=[normal(W.forests(s)[v]['literal_stars'][i])for s,i in u['literal_references']]
        groups=u['corner_orbits'];forest=u['forest'];refs=u['literal_references'];source=N.BASE
    demands=[[sum((a.get(i,F())*F(x)for i,x in z['positive_face_coefficients']),F())for z in groups]for a in(p['alpha'],p['ralpha'])]
    assert demands[0]==demand
    for k,c in enumerate(old['certificates']):W.check_certificate(old['maximal_count_vectors'],demands[k],c)
    return old,g,owner,unowned,demands,nodes,groups,forest,refs,source


def analyze(v,oi,oj,size,p,changed,priority,saved):
    old,g,owner,unowned,demands,nodes,groups,forest,refs,source=current(v,p,changed)
    key=(oi,oj);pairs=sorted(unowned[key]);assert len(pairs)==size
    assert pairs==sorted(tuple(sorted((i,j)))for i in groups[oi]['literal_stars']for j in groups[oj]['literal_stars'])
    pr=next(r for r in priority['shortlist']['chosen']if r['vertex']==v)
    assert pr['corner_orbits']==list(key)and pr['pairs']==size
    prior=next(r for r in priority['forest_records']if r['vertex']==v)if v in changed else next(r for r in W.read(O.TAG)['forest_records']if r['vertex']==v)
    scenario=next(c for c in prior['candidate_domains']if c['corner_orbits']==list(key))['scenario']
    assert O.pairhash(pairs)==scenario['pair_domain_sha256']
    oldedges={tuple(sorted(e))for e in g.edges()};assert all(pair in oldedges for pair in pairs)
    pp=[N.native(s)for s in nodes];target=R.P.forest_graph(forest,dict(R.C.models()))
    J,active,components=R.projection_data(target);checked=sorted(set(i for pair in pairs for i in pair));hs={};native_checks=[]
    for i in checked:
        star=nodes[i];cov=pp[i];edge=tuple(star['native_omitted_product']);field=N.pair_field(cov,(edge,));heights=[]
        assert all(j in active for q in cov for j,e in q)
        for q in cov:
            z=R.projection(q,q,J,field);assert not any(z[1:]);height=z[0]-2;assert height.sign()>0;heights.append(height)
        assert [z.out()for z in heights]==star['squared_timelike_heights']
        cross=R.projection(cov[0],cov[1],J,field);native=field.negative_contact(edge)
        assert field.signs_difference(tuple(a-b for a,b in zip(cross,native)),heights[0]*heights[1])==0
        hs[i]=heights;native_checks.append(dict(index=i,projection=R.qout(cross),pair_field_radicals=[field.tau.out(),field.sigma.out()],
            squared_timelike_heights=[z.out()for z in heights],native_omitted_product=list(edge)))
    cache={};expressions=[];exprids={};tested=[];bad=[];reasons=Counter()
    for i,j in pairs:
        shared=set(map(tuple,nodes[i]['occupied_pairs']))&set(map(tuple,nodes[j]['occupied_pairs']))
        why=[i,j,0]if shared else None
        if why is None:
            for a,b in((0,0),(0,1),(1,0),(1,1)):
                pi,pj=pp[i][a],pp[j][b]
                if pi==pj:continue
                ck=tuple(sorted((pi,pj)))
                if ck not in cache:
                    field=N.pair_field((pi,pj));projection=R.projection(pi,pj,J,field);h=hs[i][a]*hs[j][b]
                    sg=field.signs_difference(projection,h);sg1=field.signs_difference(R.qaddk(projection,1),h)
                    code=1 if sg>0 else 2 if sg<0 and sg1>0 else None
                    expr=dict(pair_field_radicals=[field.tau.out(),field.sigma.out()],projection=R.qout(projection),height_product=h.out(),
                        sign_product=sg,sign_product_plus_one=sg1)
                    if code is None:
                        code,extra=field.arithmetic_wall_gate(tuple((z.a,z.b)for z in projection),h.a,h.b)
                        if code is not None:expr['arithmetic_gate']=extra
                    cache[ck]=(code,expr)
                code,expr=cache[ck]
                if code is not None:
                    ek=json.dumps(expr,sort_keys=True,separators=(',',':'))
                    if ek not in exprids:exprids[ek]=len(expressions);expressions.append(expr)
                    why=[i,j,code,a,b,exprids[ek]];break
        if why is None:tested.append([i,j,None])
        else:bad.append(why);tested.append([i,j,why[2]]);reasons[why[2]]+=1;g.remove_edge(i,j)
    assert oldedges-{tuple(sorted(e))for e in g.edges()}=={tuple(x[:2])for x in bad}
    if bad:
        enum,error=O.C.vectors_for(v,g,owner,len(groups));assert error is None,error
        vectors,reps,ncliques,npatterns=enum
    else:vectors,reps,ncliques,npatterns=(old[k]for k in('maximal_count_vectors','literal_representatives','maximal_literal_cliques','distinct_count_vectors'))
    vectors=list(map(list,vectors));reps=list(map(list,reps))
    for vec,cl in zip(vectors,reps):
        assert all(g.has_edge(i,j)for i,j in combinations(cl,2))
        assert vec==[Counter(owner[i]for i in cl).get(j,0)for j in range(len(groups))]
    if saved:
        assert vectors==saved['maximal_count_vectors']and reps==saved['literal_representatives']
        full=saved['certificates'][0]
    elif len(bad)==size and scenario['method']=='exact_charge_and_cover':
        assert vectors==scenario['maximal_count_vectors']and reps==scenario['literal_representatives']
        full=scenario['certificate']
    elif bad:full=H.exact_profile(vectors,demands[0])
    else:full=old['certificates'][0]
    try:
        restricted,transport=B.projected_cover(old,vectors,demands[1]);rmode='zero_demand_projected_cover'
    except AssertionError:
        restricted=saved['certificates'][1]if saved else H.exact_profile(vectors,demands[1]);transport=[];rmode='exact_new_restricted_cover'
    certs=[full,restricted]
    for k,c in enumerate(certs):
        value=W.check_certificate(vectors,demands[k],c);assert value>=F(old['certificates'][k]['optimal_mass'])
    out=normal(dict(vertex=v,forest=forest,baseline_geometry_source=source,literal_references=refs,literal_stars=nodes,corner_orbits=groups,
        selected_orbit_pairs=[list(key)],complete_domain_sha256=O.pairhash(pairs),priority_record=pr,
        active_components=components,checked_star_indices=checked,native_checks=native_checks,
        newly_owned_pair_domain=tested,incompatibility_certificates=bad,wall_expressions=expressions,reason_counts=dict(reasons),
        maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,maximal_count_vectors=vectors,literal_representatives=reps,
        demands=[[str(x)for x in d]for d in demands],certificates=certs,restricted_cover_transport=transport,restricted_cover_method=rmode,
        baseline_optimal_masses=[c['optimal_mass']for c in old['certificates']],
        counts=dict(new_pair_tests=size,new_exclusions=len(bad),new_unknown_pairs=size-len(bad),tested_stars=len(checked),
            distinct_wall_expressions=len(expressions),literal_stars=len(nodes),corner_orbits=len(groups))))
    if saved:assert out==saved
    print('FOREST',v,out['counts'],'masses',[c['optimal_mass']for c in certs],flush=True)
    return out


def run(verify=False):
    start=time.monotonic();p=N.prepare();base=W.read(N.TAG);priority=W.read(U.TAG)
    assert not base['open_forests']and priority['baseline_full_repair']==base['repairs'][0]
    assert priority['shortlist']['complete_pair_count']==118
    assert sorted((r['vertex'],*r['corner_orbits'],r['pairs'])for r in priority['shortlist']['chosen'])==list(TARGETS)
    sources={}
    for tag in(N.TAG,U.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(N,U,O,H,W,B,R):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    saved=W.read(TAG)if verify else None
    if saved:assert saved['sources']==sources
    changed={r['vertex']:r for r in base['forest_records']};rows=[]
    for v,oi,oj,n in TARGETS:
        sr=next(r for r in saved['forest_records']if r['vertex']==v)if saved else None
        rows.append(analyze(v,oi,oj,n,p,changed,priority,sr))
    mus=[{i:F(x)for i,x in r['vertex_multipliers']}for r in base['repairs']];newmus=[dict(m)for m in mus];alphas=[p['alpha'],p['ralpha']]
    columns=[];seen=set()
    for row in rows:
        v=row['vertex']
        for k,c in enumerate(row['certificates']):
            newmus[k][v]=max(newmus[k].get(v,F()),F(c['optimal_mass']));co=Counter()
            for w,g in zip(map(F,c['weights']),row['corner_orbits']):
                for i,x in g['positive_face_coefficients']:co[i]+=w*F(x)
            terms=sorted((i,x)for i,x in co.items()if x);loads=[sum((a.get(i,F())*x for i,x in terms),F())for a in alphas]
            assert loads[k]==F(c['optimal_mass'])and all(loads[j]<=F(row['certificates'][j]['optimal_mass'])for j in range(2))
            key=(v,tuple(terms))
            if key not in seen and any(loads[j]>mus[j].get(v,F())for j in range(2)):
                seen.add(key);columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',orbit_weights=c['weights'],
                    positive_face_coefficients=[[i,str(x)]for i,x in terms],loads=list(map(str,loads)),charge_domain=base['repairs'][k]['domain']))
    repairs=[]
    for k,old in enumerate(base['repairs']):
        assert all(newmus[k].get(i,F())>=x for i,x in mus[k].items())
        num=F(old['objective_numerator']);den=sum(newmus[k].values(),F())
        repairs.append(dict(domain=old['domain'],baseline_bound=old['bound'],objective_numerator=str(num),normalizer=str(den),bound=str(num/den),
            vertex_updates=[[i,str(mus[k].get(i,F())),str(x)]for i,x in sorted(newmus[k].items())if x>mus[k].get(i,F())],
            vertex_multipliers=[[i,str(x)]for i,x in sorted(newmus[k].items())if x]))
    counts=Counter()
    for r in rows:counts.update(r['counts'])
    out=normal(dict(status='PASS_290_CLASS_SHORTLIST_ACTUAL_PAIR_GEOMETRY',sources=sources,baseline_source=N.TAG,selection_source=U.TAG,
        forest_records=rows,columns=columns,repairs=repairs,face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        counts=dict(counts),changed_forests=[r['vertex']for r in rows if r['counts']['new_exclusions']],inherited_forests=571-sum(bool(r['counts']['new_exclusions'])for r in rows),
        preservation='All older complete-owner and explicit tested outcomes, including unknowns and every raw125/raw205 test, are preserved. Only selected genuinely unowned pairs receive new outcomes. Fixed alpha and monotone mu preserve439410 inherited model columns, signed costs and previous capacities.',
        scope='Eight complete genuinely unowned orbit-pair domains, exact common finite vertex wall arithmetic and actual occurrence capacities. Unknowns remain compatible; only eight changed graph cones recomputed. No hypothetical exclusion assumed. Fresh replay has no optimizer. No global optimum, full physical feasibility, missing implicit coverage or dimension theorem.'))
    path=ROOT/f'results/{TAG}.json'
    if saved:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'columns',len(columns),'bounds',[float(F(x['bound']))for x in repairs],'seconds',round(time.monotonic()-start,2),flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');z=a.parse_args();run(z.verify)
