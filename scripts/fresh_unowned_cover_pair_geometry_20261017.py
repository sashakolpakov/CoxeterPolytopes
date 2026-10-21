#!/usr/bin/env python3
"""New physical tests on complete unowned pair domains at five cover forests.

Frozen first-owner outcomes are never changed. New exact physical exclusions
are applied only where no historical source owned the pair. Local graph cones
are then certified at the unchanged face masses; no global optimization.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_current_historical_weighted_corners_20261017 as H
import fresh_raw110_mixed_actual_corners_20261017 as R

ROOT=H.ROOT
TAG='fresh_unowned_cover_pair_geometry_20261017'
BASE=H.TAG
TARGETS=(982,2524,2190,3324,3165)
K=R.K
REASONS={0:'distinct full stars share an actual quadrilateral port pair',
    1:'positive actual wall product',2:'actual wall product strictly between minus one and zero',
    3:'nonintegral two-cycle in K',4:'two-cycle outside fixed K',5:'forbidden finite Coxeter product'}

def native(star):
    return tuple(tuple((i,tuple(e))for i,e in p)for p in star['contacts'])

def pair_owner(maps,i,j):
    return next((k for k,m in enumerate(maps)if i in m and j in m),None)

def main(verify=False):
    start=time.monotonic();base=H.W.read(BASE)
    assert base['status']=='PASS_COMPLETE_CURRENT_HISTORICAL_WEIGHTED_CORNER_CONES'
    assert not base['open_forests']
    for path,h in base['sources'].items():assert H.W.digest(path)==h,path
    # Reuse the full frozen preparation to validate current provenance and ranks.
    cert,restricted,alpha,ralpha,order,selected,owners,named,historical,sources=H.prepare()
    assert base['face_multipliers']==cert['face_multipliers']
    assert base['restricted_face_multipliers']==restricted['face_multipliers']
    assert base['source_precedence']==order
    by={r['vertex']:r for r in base['forest_records']}
    assert set(TARGETS)<=set(by)
    inputs=[cert,restricted];alphas=[alpha,ralpha]
    mus=[{v:F(x)for v,x in r['vertex_multipliers']}for r in base['repairs']]
    newmus=[dict(m)for m in mus]
    models=dict(R.C.models());rows=[];columns=[];columnkeys=set();counts=Counter()
    previous=H.W.read(TAG)if verify else None
    oldrows={r['vertex']:r for r in previous['forest_records']}if previous else {}
    for v in TARGETS:
        original=by[v];u,graph,owner=H.union_forest(v,order,selected)
        assert json.loads(json.dumps(u))==original['union']
        maps=[set(i for i,j in d['new_to_old_star_indices'])for d in u['first_owner_domains']]
        nodes=[]
        for tag,i in u['literal_references']:
            nodes.append(H.W.forests(tag)[v]['literal_stars'][i])
        assert len({H.frozen(r['contacts'])for r in nodes})==len(nodes)
        target=R.P.forest_graph(tuple(u['forest']),models)
        assert len(target)==13
        J,active,components=R.projection_data(target)
        covectors=[native(star)for star in nodes]
        heights=[];native_checks=[]
        for i,(star,pp)in enumerate(zip(nodes,covectors)):
            assert all(j in active for p in pp for j,e in p)
            edge=tuple(star['native_omitted_product']);field=R.pair_field(pp,(edge,))
            hs=[]
            for p in pp:
                q=R.projection(p,p,J,field)
                assert not any(q[1:]);h=q[0]-2;assert h.sign()>0;hs.append(h)
            assert [h.out()for h in hs]==star['squared_timelike_heights']
            p=R.projection(pp[0],pp[1],J,field);native_product=field.negative_contact(edge)
            assert field.signs_difference(tuple(a-b for a,b in zip(p,native_product)),hs[0]*hs[1])==0
            heights.append(hs)
            native_checks.append(dict(index=i,pair_field_radicals=[field.tau.out(),field.sigma.out()],
                projection=R.qout(p),squared_timelike_heights=[h.out()for h in hs],native_omitted_product=list(edge)))
        candidate=[(i,j)for i,j in combinations(range(len(nodes)),2)if pair_owner(maps,i,j)is None]
        assert len(candidate)==u['unowned_allowed_pairs']
        assert all(graph.has_edge(i,j)for i,j in candidate)
        wall_cache={};expressions=[];expression_ids={};tested=[];excluded=[];reasons=Counter()
        for i,j in candidate:
            shared=sorted(set(map(tuple,nodes[i]['occupied_pairs']))&set(map(tuple,nodes[j]['occupied_pairs'])))
            reason=None
            if shared:reason=[i,j,0];reasons[0]+=1
            else:
                for a,b in((0,0),(0,1),(1,0),(1,1)):
                    pi,pj=covectors[i][a],covectors[j][b]
                    if pi==pj:continue
                    key=tuple(sorted((pi,pj)))
                    if key not in wall_cache:
                        field=R.pair_field((pi,pj));p=R.projection(pi,pj,J,field)
                        h=heights[i][a]*heights[j][b]
                        sg=field.signs_difference(p,h);sg1=field.signs_difference(R.qaddk(p,1),h)
                        code=1 if sg>0 else 2 if sg<0 and sg1>0 else None
                        expr=dict(pair_field_radicals=[field.tau.out(),field.sigma.out()],projection=R.qout(p),
                            height_product=h.out(),sign_product=sg,sign_product_plus_one=sg1)
                        if code is None:
                            code,extra=field.arithmetic_wall_gate(tuple((z.a,z.b)for z in p),h.a,h.b)
                            if code is not None:expr['arithmetic_gate']=extra
                        wall_cache[key]=(code,expr)
                    code,expr=wall_cache[key]
                    if code is not None:
                        ek=json.dumps(expr,sort_keys=True,separators=(',',':'))
                        if ek not in expression_ids:expression_ids[ek]=len(expressions);expressions.append(expr)
                        reason=[i,j,code,a,b,expression_ids[ek]];reasons[code]+=1;break
            if reason is None:tested.append([i,j,None])
            else:
                excluded.append(reason);tested.append([i,j,reason[2]])
                graph.remove_edge(i,j)
        enum,error=H.vectors_for(graph,owner,len(u['corner_orbits']))
        assert error is None,(v,error)
        vectors,reps,ncliques,npatterns=enum;ccs=[]
        for k in range(2):
            demand=list(map(F,original['demands'][k]))
            # All actual occurrence coefficients and both face masses are unchanged.
            expected=[sum((alphas[k].get(fi,F())*F(m)for fi,m in g['positive_face_coefficients']),F())for g in u['corner_orbits']]
            assert demand==expected
            if verify:
                old=oldrows[v]
                assert old['maximal_count_vectors']==list(map(list,vectors))
                assert old['literal_representatives']==list(map(list,reps))
                c=old['certificates'][k]
            else:c=H.exact_profile(vectors,demand)
            value=H.W.check_certificate(vectors,demand,c);ccs.append(c)
            assert value>=F(original['certificates'][k]['optimal_mass'])
            newmus[k][v]=max(newmus[k].get(v,F()),value)
        for vec,clique in zip(vectors,reps):
            assert all(graph.has_edge(i,j)for i,j in combinations(clique,2))
            assert tuple(Counter(owner[i]for i in clique).get(j,0)for j in range(len(u['corner_orbits'])))==vec
        for k,c in enumerate(ccs):
            co=Counter()
            for w,g in zip(map(F,c['weights']),u['corner_orbits']):
                for fi,m in g['positive_face_coefficients']:co[fi]+=w*F(m)
            terms=sorted((fi,x)for fi,x in co.items()if x)
            loads=[sum((a.get(fi,F())*x for fi,x in terms),F())for a in alphas]
            assert loads[k]==F(c['optimal_mass'])
            assert all(loads[j]<=F(ccs[j]['optimal_mass'])for j in range(2))
            key=(v,tuple(terms))
            if key not in columnkeys and any(loads[j]>mus[j].get(v,F())for j in range(2)):
                columnkeys.add(key);columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',
                    positive_face_coefficients=[[fi,str(x)]for fi,x in terms],orbit_weights=c['weights'],
                    charge_domain=base['repairs'][k]['domain'],loads=list(map(str,loads))))
        rows.append(dict(vertex=v,forest=u['forest'],literal_references=u['literal_references'],
            active_components=components,native_checks=native_checks,
            newly_owned_pair_domain=tested,incompatibility_certificates=excluded,wall_expressions=expressions,
            reason_counts=dict(reasons),maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,
            maximal_count_vectors=vectors,literal_representatives=reps,certificates=ccs,
            baseline_optimal_masses=[c['optimal_mass']for c in original['certificates']],
            baseline_mu=[str(m.get(v,F()))for m in mus]))
        counts.update(forests=1,literal_stars=len(nodes),new_pair_tests=len(candidate),new_exclusions=len(excluded),
            new_unknown_pairs=len(candidate)-len(excluded),distinct_wall_expressions=len(expressions),
            maximal_literal_cliques=ncliques,maximal_count_vectors=len(vectors),
            full_mu_increases=int(newmus[0][v]>mus[0].get(v,F())),restricted_mu_increases=int(newmus[1][v]>mus[1].get(v,F())))
        print('FOREST',v,'tested',len(candidate),'excluded',len(excluded),'codes',dict(reasons),
            'old/new full optimal',original['certificates'][0]['optimal_mass'],ccs[0]['optimal_mass'],flush=True)
    assert counts['new_pair_tests']==2286
    repairs=[]
    for k,(oldr,newm)in enumerate(zip(base['repairs'],newmus)):
        objective=F(oldr['objective_numerator']);normalizer=sum(newm.values(),F())
        assert all(newm.get(v,F())>=x for v,x in mus[k].items())
        repairs.append(dict(domain=oldr['domain'],baseline_bound=oldr['bound'],objective_numerator=str(objective),
            normalizer=str(normalizer),bound=str(objective/normalizer),
            vertex_updates=[[v,str(mus[k].get(v,F())),str(x)]for v,x in sorted(newm.items())if x>mus[k].get(v,F())],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(newm.items())if x]))
    paths=[Path(__file__),Path(H.__file__),Path(R.__file__),Path(R.T.__file__),Path(R.P.__file__),Path(R.M.__file__),Path(R.V.__file__)]
    sources={str(p.relative_to(ROOT)):H.W.digest(str(p.relative_to(ROOT)))for p in paths}
    for folder,ext in(('scripts','py'),('results','json'),('results','md')):
        p=f'{folder}/{BASE}.{ext}';sources[p]=H.W.digest(p)
    result=dict(status='PASS_EXACT_NEW_UNOWNED_COVER_PAIR_GEOMETRY',sources=sources,
        target_forests=list(TARGETS),counts=dict(counts),forest_records=rows,columns=columns,repairs=repairs,
        face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        pair_reason_codes=REASONS,fixed_ambient_ground_field='Q(sqrt(5))',
        preservation='Same full5053 face masses and same restricted masses. Only vertex multipliers increase; all439410 inherited finite-model columns and all78 prior weighted columns remain valid. The unchanged545 graph cones retain their same-mass exact covers; the five changed graphs receive new exact covers.',
        scope='Complete newly owned pair domains only in the five specified forests. Inherited first-owner outcomes, including unknowns, are immutable. All unexcluded and all other never-owned pairs stay compatible. New physical exclusions use common actual finite-vertex contact and future-height data. No global optimizer, complete physical realizability or dimension theorem.')
    result=json.loads(json.dumps(result,sort_keys=True));path=ROOT/f'results/{TAG}.json'
    if verify:assert result==previous
    else:path.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(result['status'],dict(counts),'columns',len(columns),'bounds',
        [float(F(r['bound']))for r in repairs],'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();main(a.verify)
