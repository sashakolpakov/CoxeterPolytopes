#!/usr/bin/env python3
"""Two complete unowned raw88/raw25 pair domains; exact local cone repair.

The baseline includes the frozen eight-pair3324 successor. Only changed
forests receive new clique enumerations; --verify uses saved rational covers.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_historical22_corner_domain_20261017 as D
import fresh_raw110_mixed_actual_corners_20261017 as R

ROOT=D.ROOT;H=D.H;W=D.W
TAG='fresh_historical22_raw88_raw25_pair_geometry_20261017'
BASE='fresh_historical22_eight_pair_geometry_20261017'
CONES='fresh_historical22_corner_cones_20261017'
PRIORITY='fresh_historical22_active_pair_diagnostic_20261017'
TARGETS=((944,2,20,72),(788,2,23,288))
REASONS={0:'distinct full stars share an actual quadrilateral port pair',
 1:'positive actual wall product',2:'actual wall product strictly between minus one and zero',
 3:'nonintegral two-cycle in K',4:'two-cycle outside fixed K',5:'forbidden finite Coxeter product'}


def normalize(x):return json.loads(json.dumps(x,sort_keys=True))
def native(star):return tuple(tuple((i,tuple(e))for i,e in p)for p in star['contacts'])


def projected_cover(old,vectors,demands):
    """Delete zero-demand coordinates; select dominating complete new vectors."""
    original=old['certificates'][1];mix=Counter();transport=[]
    for vi,amount in original['mixture']:
        reduced=[n if demands[j]>0 else 0 for j,n in enumerate(old['maximal_count_vectors'][vi])]
        matches=[i for i,v in enumerate(vectors)if all(a>=b for a,b in zip(v,reduced))]
        assert matches,'positive-demand projected old clique missing from new domain'
        chosen=min(matches);mix[chosen]+=F(amount)
        transport.append(dict(old_vector=vi,positive_projection=reduced,new_dominating_vector=chosen,mass=amount))
    cert=dict(weights=original['weights'],mixture=[[i,str(x)]for i,x in sorted(mix.items())if x],
              optimal_mass=original['optimal_mass'])
    W.check_certificate(vectors,demands,cert)
    return cert,transport


def analyze(v,oi,oj,npairs,p,domain,cones,priority,saved):
    u,graph,owner=D.assemble_forest(v,p)
    assert normalize(u)==next(r for r in domain['forest_records']if r['vertex']==v)
    old=next(r for r in cones['forest_records']if r['vertex']==v)
    assert u['demands']==old['demands']
    groups=u['corner_orbits'];refs=u['literal_references'];oldedges={tuple(sorted(e))for e in graph.edges()}
    gi,gj=groups[oi],groups[oj]
    assert (gi['face'],gj['face'])==(45649,26798)
    assert {c['raw_family']for c in gi['origins']}=={88}and {c['raw_family']for c in gj['origins']}=={25}
    assert all(c['provenance']==982 for c in gi['origins'])
    assert p['selected'][26798]['provenance_kind']=='unique_saved_context'
    candidate=sorted(tuple(sorted((i,j)))for i in gi['literal_stars']for j in gj['literal_stars'])
    assert len(candidate)==len(set(candidate))==npairs
    domains=[{i for i,j in d['new_to_old_star_indices']}for d in u['source_domains']]
    explicit={}
    if v in p['patch_rows']:
        tr={j:i for i,j in u['explicit_patch']['new_to_patch_star_indices']}
        explicit={tuple(sorted((tr[i],tr[j]))):code for i,j,code in p['patch_rows'][v]['newly_owned_pair_domain']}
    for pair in candidate:
        assert pair in oldedges and pair not in explicit
        assert not any(all(i in ss for i in pair)for ss in domains),'already-owned outcome'
    pr=next(c for c in priority['candidate_domains']if c['vertex']==v and c['kind']=='unowned'and c['corner_orbits']==[oi,oj])
    assert sha256(json.dumps(candidate,separators=(',',':')).encode()).hexdigest()==pr['pair_domain_sha256']
    assert pr['complete_domain_pairs']==npairs and pr['both_endpoint_demands_positive']==[True,False]
    nodes=[W.forests(tag)[v]['literal_stars'][i]for tag,i in refs]
    covectors=[native(s)for s in nodes]
    target=R.P.forest_graph(tuple(u['forest']),dict(R.C.models()));assert len(target)==13
    J,active,components=R.projection_data(target)
    checked=sorted(set(gi['literal_stars']+gj['literal_stars']));heights={};native_checks=[]
    for i in checked:
        star=nodes[i];pp=covectors[i]
        assert all(j in active for q in pp for j,e in q)
        edge=tuple(star['native_omitted_product']);field=R.pair_field(pp,(edge,));hs=[]
        for q in pp:
            z=R.projection(q,q,J,field);assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;hs.append(h)
        assert [h.out()for h in hs]==star['squared_timelike_heights']
        cross=R.projection(pp[0],pp[1],J,field);native_product=field.negative_contact(edge)
        assert field.signs_difference(tuple(a-b for a,b in zip(cross,native_product)),hs[0]*hs[1])==0
        heights[i]=hs
        native_checks.append(dict(index=i,pair_field_radicals=[field.tau.out(),field.sigma.out()],
            projection=R.qout(cross),squared_timelike_heights=[h.out()for h in hs],native_omitted_product=list(edge)))
    cache={};expressions=[];exprids={};tested=[];excluded=[];reasons=Counter()
    for i,j in candidate:
        shared=sorted(set(map(tuple,nodes[i]['occupied_pairs']))&set(map(tuple,nodes[j]['occupied_pairs'])))
        why=None
        if shared:why=[i,j,0]
        else:
            for a,b in((0,0),(0,1),(1,0),(1,1)):
                pi,pj=covectors[i][a],covectors[j][b]
                if pi==pj:continue
                key=tuple(sorted((pi,pj)))
                if key not in cache:
                    field=R.pair_field((pi,pj));projection=R.projection(pi,pj,J,field);h=heights[i][a]*heights[j][b]
                    sg=field.signs_difference(projection,h);sg1=field.signs_difference(R.qaddk(projection,1),h)
                    code=1 if sg>0 else 2 if sg<0 and sg1>0 else None
                    expr=dict(pair_field_radicals=[field.tau.out(),field.sigma.out()],projection=R.qout(projection),
                        height_product=h.out(),sign_product=sg,sign_product_plus_one=sg1)
                    if code is None:
                        code,extra=field.arithmetic_wall_gate(tuple((z.a,z.b)for z in projection),h.a,h.b)
                        if code is not None:expr['arithmetic_gate']=extra
                    cache[key]=(code,expr)
                code,expr=cache[key]
                if code is not None:
                    ek=json.dumps(expr,sort_keys=True,separators=(',',':'))
                    if ek not in exprids:exprids[ek]=len(expressions);expressions.append(expr)
                    why=[i,j,code,a,b,exprids[ek]];break
        if why is None:tested.append([i,j,None])
        else:excluded.append(why);tested.append([i,j,why[2]]);reasons[why[2]]+=1;graph.remove_edge(i,j)
    assert oldedges-{tuple(sorted(e))for e in graph.edges()}=={tuple(r[:2])for r in excluded}
    if excluded:
        enum,error=H.vectors_for(graph,owner,len(groups));assert error is None,error
        vectors,reps,ncliques,npatterns=enum
    else:
        vectors=old['maximal_count_vectors'];reps=old['literal_representatives']
        ncliques=old['maximal_literal_cliques'];npatterns=old['distinct_count_vectors']
    vectors=list(map(list,vectors));reps=list(map(list,reps))
    for vec,clique in zip(vectors,reps):
        assert all(graph.has_edge(i,j)for i,j in combinations(clique,2))
        assert vec==[Counter(owner[i]for i in clique).get(j,0)for j in range(len(groups))]
    if saved is not None:
        assert saved['maximal_count_vectors']==vectors and saved['literal_representatives']==reps
        full=saved['certificates'][0]
    elif excluded:full=H.exact_profile(vectors,list(map(F,u['demands'][0])))
    else:full=old['certificates'][0]
    restricted,transport=projected_cover(old,vectors,list(map(F,u['demands'][1])))
    certificates=[full,restricted]
    for k,c in enumerate(certificates):
        value=W.check_certificate(vectors,list(map(F,u['demands'][k])),c)
        assert value>=F(old['certificates'][k]['optimal_mass'])
    assert certificates[1]['optimal_mass']==old['certificates'][1]['optimal_mass']
    record=dict(vertex=v,forest=u['forest'],corner_orbits=groups,literal_references=refs,
        selected_orbit_pairs=[[oi,oj]],priority_record=pr['index'],complete_domain_sha256=pr['pair_domain_sha256'],
        active_components=components,checked_star_indices=checked,native_checks=native_checks,
        newly_owned_pair_domain=tested,incompatibility_certificates=excluded,wall_expressions=expressions,
        reason_counts=dict(reasons),maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,
        maximal_count_vectors=vectors,literal_representatives=reps,demands=u['demands'],certificates=certificates,
        restricted_cover_transport=transport,baseline_optimal_masses=[c['optimal_mass']for c in old['certificates']],
        prior_owner_counts=u['counts'],counts=dict(new_pair_tests=npairs,new_exclusions=len(excluded),
            new_unknown_pairs=npairs-len(excluded),distinct_wall_expressions=len(expressions),
            tested_stars=len(checked),literal_stars=len(nodes),corner_orbits=len(groups)))
    print('FOREST',v,record['counts'],'codes',dict(reasons),'masses',[c['optimal_mass']for c in certificates],flush=True)
    return record


def run(verify=False):
    start=time.monotonic();p=D.prepare();baseline=W.read(BASE);cones=W.read(CONES);priority=W.read(PRIORITY);domain=W.read(D.TAG)
    assert baseline['status']=='PASS_HISTORICAL22_EIGHT_UNOWNED_PHYSICAL_PAIR_TESTS'
    assert cones['status']=='PASS_COMPLETE_HISTORICAL22_CHANGED_CORNER_CONES'and not cones['open_forests']
    assert baseline['vertex']==3324 and baseline['vertex']not in[v for v,*_ in TARGETS]
    for key in('face_multipliers','restricted_face_multipliers'):
        assert baseline[key]==cones[key]==domain[key]
    assert priority['fixed_repairs']==cones['repairs']
    sources={}
    for tag in(BASE,CONES,PRIORITY,D.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for path in(Path(__file__),Path(D.__file__),Path(H.__file__),Path(W.__file__),Path(R.__file__),Path(R.T.__file__),Path(R.P.__file__),Path(R.M.__file__),Path(R.V.__file__)):
        rel=str(path.relative_to(ROOT));sources[rel]=W.digest(rel)
    saved=W.read(TAG)if verify else None;records=[]
    for v,oi,oj,npairs in TARGETS:
        old=next(r for r in saved['forest_records']if r['vertex']==v)if verify else None
        records.append(analyze(v,oi,oj,npairs,p,domain,cones,priority,old))
    alphas=[p['alpha'],p['ralpha']]
    mus=[{v:F(x)for v,x in r['vertex_multipliers']}for r in baseline['repairs']]
    newmus=[dict(m)for m in mus];columns=[];seen=set()
    for r in records:
        v=r['vertex'];groups=r['corner_orbits']
        for k,c in enumerate(r['certificates']):
            newmus[k][v]=max(newmus[k].get(v,F()),F(c['optimal_mass']))
            coef=Counter()
            for w,g in zip(map(F,c['weights']),groups):
                for fi,m in g['positive_face_coefficients']:coef[fi]+=w*F(m)
            terms=sorted((fi,x)for fi,x in coef.items()if x)
            loads=[sum((a.get(fi,F())*x for fi,x in terms),F())for a in alphas]
            assert loads[k]==F(c['optimal_mass'])and all(loads[j]<=F(r['certificates'][j]['optimal_mass'])for j in range(2))
            key=(v,tuple(terms))
            if key not in seen and any(loads[j]>mus[j].get(v,F())for j in range(2)):
                seen.add(key);columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',orbit_weights=c['weights'],
                    positive_face_coefficients=[[fi,str(x)]for fi,x in terms],loads=list(map(str,loads)),charge_domain=baseline['repairs'][k]['domain']))
    repairs=[]
    for k,(old,nm)in enumerate(zip(baseline['repairs'],newmus)):
        assert all(nm.get(v,F())>=x for v,x in mus[k].items())
        numerator=F(old['objective_numerator']);den=sum(nm.values(),F())
        repairs.append(dict(domain=old['domain'],baseline_bound=old['bound'],objective_numerator=str(numerator),normalizer=str(den),bound=str(numerator/den),
            vertex_updates=[[v,str(mus[k].get(v,F())),str(x)]for v,x in sorted(nm.items())if x>mus[k].get(v,F())],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(nm.items())if x]))
    assert repairs[1]['bound']==baseline['repairs'][1]['bound']and not repairs[1]['vertex_updates']
    totals=Counter()
    for r in records:totals.update(r['counts'])
    changed=[r['vertex']for r in records if r['counts']['new_exclusions']]
    out=normalize(dict(status='PASS_HISTORICAL22_RAW88_RAW25_COMPLETE_UNOWNED_PAIR_DOMAINS',sources=sources,
        baseline_source=BASE,cone_predecessor=CONES,forest_records=records,columns=columns,repairs=repairs,
        face_multipliers=baseline['face_multipliers'],restricted_face_multipliers=baseline['restricted_face_multipliers'],
        counts=dict(totals),changed_forests=changed,inherited_forests=553-len(changed),pair_reason_codes=REASONS,
        preserved_explicit_pair_source=BASE,
        preservation='Each new domain is outside all historical owners and all earlier explicit tested pairs. Owned unknowns remain allowed; the frozen3324 eight-pair successor is inherited. Fixed face masses and monotone vertex increases preserve all439410 model columns and prior weighted capacities.',
        scope='Two complete raw88/raw25 orbit-pair domains in full common finite rank13 frames. Only changed graph cones recomputed. Complete projected restricted cover transported without a local optimizer. Fresh replay checks saved rational covers, no optimizer. No full physical realization, unchanged-global reoptimization, implicit old-alpha simplex/joint coverage or dimension theorem.'))
    path=ROOT/f'results/{TAG}.json'
    if verify:assert out==saved
    else:path.write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'columns',len(columns),'bounds',[float(F(x['bound']))for x in repairs],
          'seconds',round(time.monotonic()-start,2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
