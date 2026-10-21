#!/usr/bin/env python3
"""Exact raw33/raw39 tests on two complete four-pair domains at forest3324.

Only genuinely unowned pairs are tested. Frozen unknown outcomes are retained.
The changed forest receives exact local covers; all552 others are inherited.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_historical22_corner_domain_20261017 as D
import fresh_raw110_mixed_actual_corners_20261017 as R

ROOT=D.ROOT;H=D.H;W=D.W
TAG='fresh_historical22_eight_pair_geometry_20261017'
BASE='fresh_historical22_corner_cones_20261017'
PRIORITY='fresh_historical22_active_pair_diagnostic_20261017'
VERTEX=3324;ORBIT_PAIRS=((5,7),(5,8))
REASONS={0:'distinct full stars share an actual quadrilateral port pair',
 1:'positive actual wall product',2:'actual wall product strictly between minus one and zero',
 3:'nonintegral two-cycle in K',4:'two-cycle outside fixed K',5:'forbidden finite Coxeter product'}

def native(star):return tuple(tuple((i,tuple(e))for i,e in p)for p in star['contacts'])
def normalize(x):return json.loads(json.dumps(x,sort_keys=True))

def run(verify=False):
    start=time.monotonic();p=D.prepare();base=W.read(BASE);priority=W.read(PRIORITY);domain=W.read(D.TAG)
    assert base['status']=='PASS_COMPLETE_HISTORICAL22_CHANGED_CORNER_CONES'and not base['open_forests']
    assert priority['status']=='PASS_HISTORICAL22_EXACT_COVER_PAIR_PRIORITY_DIAGNOSTIC'
    sources={}
    for tag in(BASE,PRIORITY,D.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for path in(Path(__file__),Path(D.__file__),Path(H.__file__),Path(W.__file__),Path(R.__file__),Path(R.T.__file__),Path(R.P.__file__),Path(R.M.__file__),Path(R.V.__file__)):
        rel=str(path.relative_to(ROOT));sources[rel]=W.digest(rel)
    assert base['face_multipliers']==domain['face_multipliers']
    assert base['restricted_face_multipliers']==domain['restricted_face_multipliers']
    assert priority['fixed_repairs']==base['repairs']
    saved=W.read(TAG)if verify else None
    u,graph,owner=D.assemble_forest(VERTEX,p)
    assert normalize(u)==next(r for r in domain['forest_records']if r['vertex']==VERTEX)
    groups=u['corner_orbits'];refs=u['literal_references'];oldedges={tuple(sorted(e))for e in graph.edges()}
    original=next(r for r in base['forest_records']if r['vertex']==VERTEX)
    maps=[set(i for i,j in d['new_to_old_star_indices'])for d in u['source_domains']]
    tr={j:i for i,j in u['explicit_patch']['new_to_patch_star_indices']}
    explicit={tuple(sorted((tr[i],tr[j])))for i,j,code in p['patch_rows'][VERTEX]['newly_owned_pair_domain']}
    candidate=[];subdomains=[]
    for oi,oj in ORBIT_PAIRS:
        gi,gj=groups[oi],groups[oj]
        assert (gi['face'],gj['face'])==(27388,28798)
        assert {c['raw_family']for c in gi['origins']}=={33}and {c['raw_family']for c in gj['origins']}=={39}
        pairs=sorted(tuple(sorted((i,j)))for i in gi['literal_stars']for j in gj['literal_stars'])
        assert len(pairs)==4 and len(set(pairs))==4
        for pair in pairs:
            assert pair in oldedges and pair not in explicit
            assert not any(all(i in ss for i in pair)for ss in maps),'already-owned pair'
        pr=next(c for c in priority['candidate_domains']if c['vertex']==VERTEX and c['kind']=='unowned'and c['corner_orbits']==[oi,oj])
        from hashlib import sha256
        assert sha256(json.dumps(pairs,separators=(',',':')).encode()).hexdigest()==pr['pair_domain_sha256']
        assert pr['complete_domain_pairs']==4 and pr['both_endpoint_demands_positive']==[True,False]
        candidate.extend(pairs);subdomains.append(dict(corner_orbits=[oi,oj],pairs=pairs,priority_record=pr['index']))
    candidate=sorted(candidate);assert len(candidate)==len(set(candidate))==8
    nodes=[W.forests(tag)[VERTEX]['literal_stars'][i]for tag,i in refs]
    target=R.P.forest_graph(tuple(u['forest']),dict(R.C.models()));assert len(target)==13
    J,active,components=R.projection_data(target);covectors=[native(s)for s in nodes]
    heights=[];native_checks=[]
    for i,(star,pp)in enumerate(zip(nodes,covectors)):
        assert all(j in active for q in pp for j,e in q)
        edge=tuple(star['native_omitted_product']);field=R.pair_field(pp,(edge,));hs=[]
        for q in pp:
            z=R.projection(q,q,J,field);assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;hs.append(h)
        assert [h.out()for h in hs]==star['squared_timelike_heights']
        cross=R.projection(pp[0],pp[1],J,field);native_product=field.negative_contact(edge)
        assert field.signs_difference(tuple(a-b for a,b in zip(cross,native_product)),hs[0]*hs[1])==0
        heights.append(hs)
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
                    key=json.dumps(expr,sort_keys=True,separators=(',',':'))
                    if key not in exprids:exprids[key]=len(expressions);expressions.append(expr)
                    why=[i,j,code,a,b,exprids[key]];break
        if why is None:tested.append([i,j,None])
        else:excluded.append(why);tested.append([i,j,why[2]]);reasons[why[2]]+=1;graph.remove_edge(i,j)
    assert oldedges-{tuple(sorted(e))for e in graph.edges()}=={tuple(r[:2])for r in excluded}
    enum,error=H.vectors_for(graph,owner,len(groups));assert error is None,error
    vectors,reps,ncliques,npatterns=enum
    alphas=[p['alpha'],p['ralpha']];ccs=[]
    for k in range(2):
        demand=list(map(F,u['demands'][k]));assert u['demands'][k]==original['demands'][k]
        if verify:
            assert saved['maximal_count_vectors']==list(map(list,vectors))and saved['literal_representatives']==list(map(list,reps))
            c=saved['certificates'][k]
        elif not excluded:c=original['certificates'][k]
        else:c=H.exact_profile(vectors,demand)
        value=W.check_certificate(vectors,demand,c)
        assert value>=F(original['certificates'][k]['optimal_mass']);ccs.append(c)
    # Both tested orbit pairs have a zero-demand restricted endpoint. Their
    # deletion cannot change the projection onto positive-demand coordinates.
    assert ccs[1]['optimal_mass']==original['certificates'][1]['optimal_mass']
    for vec,clique in zip(vectors,reps):
        assert all(graph.has_edge(i,j)for i,j in combinations(clique,2))
        assert tuple(Counter(owner[i]for i in clique).get(j,0)for j in range(len(groups)))==vec
    mus=[{v:F(x)for v,x in r['vertex_multipliers']}for r in base['repairs']];newmus=[dict(m)for m in mus]
    columns=[];seen=set()
    for k,c in enumerate(ccs):
        newmus[k][VERTEX]=max(newmus[k].get(VERTEX,F()),F(c['optimal_mass']))
        coef=Counter()
        for w,g in zip(map(F,c['weights']),groups):
            for fi,m in g['positive_face_coefficients']:coef[fi]+=w*F(m)
        terms=sorted((fi,x)for fi,x in coef.items()if x)
        loads=[sum((a.get(fi,F())*x for fi,x in terms),F())for a in alphas]
        assert loads[k]==F(c['optimal_mass'])and all(loads[j]<=F(ccs[j]['optimal_mass'])for j in range(2))
        key=tuple(terms)
        if key not in seen and any(loads[j]>mus[j].get(VERTEX,F())for j in range(2)):
            seen.add(key);columns.append(dict(index=len(columns),vertex=VERTEX,vertex_cost='1',orbit_weights=c['weights'],
                positive_face_coefficients=[[fi,str(x)]for fi,x in terms],loads=list(map(str,loads)),charge_domain=base['repairs'][k]['domain']))
    repairs=[]
    for k,(old,nm)in enumerate(zip(base['repairs'],newmus)):
        assert all(nm.get(v,F())>=x for v,x in mus[k].items())
        numerator=F(old['objective_numerator']);den=sum(nm.values(),F())
        repairs.append(dict(domain=old['domain'],baseline_bound=old['bound'],objective_numerator=str(numerator),normalizer=str(den),bound=str(numerator/den),
            vertex_updates=[[v,str(mus[k].get(v,F())),str(x)]for v,x in sorted(nm.items())if x>mus[k].get(v,F())],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(nm.items())if x]))
    out=dict(status='PASS_HISTORICAL22_EIGHT_UNOWNED_PHYSICAL_PAIR_TESTS',sources=sources,baseline_source=BASE,
        vertex=VERTEX,forest=u['forest'],literal_references=refs,corner_orbits=groups,subdomains=subdomains,
        active_components=components,native_checks=native_checks,newly_owned_pair_domain=tested,
        incompatibility_certificates=excluded,wall_expressions=expressions,reason_counts=dict(reasons),pair_reason_codes=REASONS,
        maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,maximal_count_vectors=vectors,literal_representatives=reps,
        demands=u['demands'],certificates=ccs,baseline_optimal_masses=[c['optimal_mass']for c in original['certificates']],
        columns=columns,repairs=repairs,face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        counts=dict(new_pair_tests=len(candidate),new_exclusions=len(excluded),new_unknown_pairs=len(candidate)-len(excluded),
            distinct_wall_expressions=len(expressions),literal_stars=len(nodes),corner_orbits=len(groups),inherited_forests=552),
        preservation='Only eight genuinely unowned pairs become explicitly owned. Every frozen owned exclusion/unknown is preserved. Both face-mass lists are unchanged and only vertex masses may increase; all439410 model columns and earlier weighted capacities stay valid. The other552 complete graph cones inherit same-demand certificates.',
        scope='Two complete raw33/raw39 orbit-pair domains at one actual finite-vertex forest. Per-wall-pair immutable squareclass fields and full native contact/future-height data. Unknowns remain allowed. No global optimization, full physical realizability, implicit old-alpha simplex/joint coverage or dimension theorem.')
    out=normalize(out);path=ROOT/f'results/{TAG}.json'
    if verify:assert out==saved
    else:path.write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'reasons',dict(reasons),'columns',len(columns),'bounds',[float(F(x['bound']))for x in repairs],
        'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');x=a.parse_args();run(x.verify)
