#!/usr/bin/env python3
"""Exact graph-cone certificates for the75 changed historical22 forests.

All478 unchanged forests inherit their same-demand frozen certificates.
Generation may solve bounded local covering LPs; --verify uses only saved
rational charges/covers and never invokes an optimizer. No global solve.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_historical22_corner_domain_20261017 as D

ROOT=D.ROOT;H=D.H;W=D.W
TAG='fresh_historical22_corner_cones_20261017'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.json')

def normalize(x):return json.loads(json.dumps(x,sort_keys=True))

def vectors_for(v,graph,owner,n):
    if v!=143:return H.vectors_for(graph,owner,n)
    # This one literal domain has many symmetry presentations but only a
    # small count-vector domain. Its explicit finite guard is larger; no
    # partial enumeration is ever accepted as a certificate.
    patterns={};count=0;start=time.monotonic()
    for clique in D.nx.find_cliques(graph):
        count+=1;clique=tuple(sorted(clique));c=Counter(owner[i]for i in clique)
        vec=tuple(c.get(j,0)for j in range(n))
        if vec not in patterns or clique<patterns[vec]:patterns[vec]=clique
        if count>3000000 or time.monotonic()-start>30:
            return None,dict(reason='literal_clique_enumeration_guard',enumerated_cliques=count)
    vv=sorted(patterns)
    maximal=[a for a in vv if not any(a!=b and all(x<=y for x,y in zip(a,b))for b in vv)]
    return(maximal,[patterns[a]for a in maximal],count,len(patterns)),None

def run(verify=False,resume=False):
    start=time.monotonic();p=D.prepare();domain=W.read(D.TAG)
    assert domain['status']=='PASS_HISTORICAL22_CURRENT_PROVENANCE_AND_NATIVE_DOMAINS'
    for path,h in domain['sources'].items():assert W.digest(path)==h,path
    assert domain['face_multipliers']==p['base']['face_multipliers']==p['patch']['face_multipliers']
    assert domain['restricted_face_multipliers']==p['base']['restricted_face_multipliers']==p['patch']['restricted_face_multipliers']
    assert domain['fixed_baseline_repairs']==p['patch']['repairs']
    assert p['patch']['status']=='PASS_EXACT_NEW_UNOWNED_COVER_PAIR_GEOMETRY'
    assert not p['base']['open_forests']
    dr={r['vertex']:r for r in domain['forest_records']}
    assert sorted(dr)==p['affected'] and len(dr)==75
    previous=W.read(TAG)if verify else None
    sources=dict(p['sources'])
    for tag in(D.TAG,D.BASE,D.PATCH):
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for path in(Path(__file__),Path(H.__file__),Path(W.__file__)):
        rel=str(path.relative_to(ROOT));sources[rel]=W.digest(rel)
    binding=dict(sources=sources,face_multipliers=domain['face_multipliers'],
        restricted_face_multipliers=domain['restricted_face_multipliers'])
    cached={}
    if verify:
        assert previous['sources']==sources
        cached={r['vertex']:r for r in previous['forest_records']}
    elif resume and PROGRESS.exists():
        progress=json.loads(PROGRESS.read_text());assert progress['binding']==binding
        cached={r['vertex']:r for r in progress['forest_records']}
    alphas=[p['alpha'],p['ralpha']]
    mus=[{v:F(x)for v,x in r['vertex_multipliers']}for r in p['patch']['repairs']]
    newmus=[dict(m)for m in mus]
    labels=[r['domain']for r in p['patch']['repairs']]
    rows=[];opened=[];columns=[];columnkeys=set();counts=Counter()
    inherited=[]
    for v,old in sorted(p['original_rows'].items()):
        if v in dr:continue
        source=D.PATCH if v in p['patch_rows']else D.BASE
        ccs=(p['patch_rows'][v]if source==D.PATCH else old)['certificates']
        assert all(mus[k].get(v,F())>=F(c['optimal_mass'])for k,c in enumerate(ccs))
        inherited.append(dict(vertex=v,source=source,optimal_masses=[c['optimal_mass']for c in ccs]))
    assert len(inherited)==478
    assert {r['vertex']for r in inherited if r['source']==D.PATCH}=={2190}
    print('DOMAIN75 changed;478 inherited;270 positive represented classes',flush=True)
    for step,v in enumerate(p['affected'],1):
        u,graph,owner=D.assemble_forest(v,p)
        assert normalize(u)==dr[v],(v,'frozen domain reconstruction')
        groups=u['corner_orbits']
        enum,error=vectors_for(v,graph,owner,len(groups))
        if error:
            opened.append(dict(vertex=v,**error));print('OPEN',v,error,flush=True);continue
        vectors,reps,ncliques,npatterns=enum
        demands=[[sum((a.get(fi,F())*F(m)for fi,m in g['positive_face_coefficients']),F())for g in groups]for a in alphas]
        assert all(x>0 for x in demands[0])
        assert [list(map(str,d))for d in demands]==u['demands']
        old=cached.get(v)
        if old:
            assert old['maximal_count_vectors']==list(map(list,vectors))
            assert old['literal_representatives']==list(map(list,reps))
            assert old['demands']==u['demands']
        ccs=[]
        try:
            for k,demand in enumerate(demands):
                if old:c=old['certificates'][k]
                else:
                    assert not verify,(v,'missing saved exact certificate')
                    c=H.exact_profile(vectors,demand)
                value=W.check_certificate(vectors,demand,c);assert value==F(c['optimal_mass'])
                ccs.append(c)
        except Exception as e:
            if verify:raise
            opened.append(dict(vertex=v,reason='exact_local_certificate_recovery',detail=repr(e),vectors=len(vectors)))
            print('OPEN',v,repr(e),flush=True);continue
        # Every saved mixture column is an actual literal clique in the
        # conservative graph; counts use complete physical corner orbits.
        for vec,clique in zip(vectors,reps):
            assert len(clique)==len(set(clique))
            assert all(graph.has_edge(i,j)for i,j in combinations(clique,2))
            assert tuple(Counter(owner[i]for i in clique).get(j,0)for j in range(len(groups)))==vec
        baseline_optimal=[]
        if v in p['original_rows']:
            oldrow=p['patch_rows'].get(v,p['original_rows'][v])
            baseline_optimal=[c['optimal_mass']for c in oldrow['certificates']]
            assert all(F(ccs[k]['optimal_mass'])>=F(x)for k,x in enumerate(baseline_optimal))
        for k,c in enumerate(ccs):
            value=F(c['optimal_mass']);newmus[k][v]=max(newmus[k].get(v,F()),value)
            co=Counter()
            for w,g in zip(map(F,c['weights']),groups):
                assert w>=0
                for fi,m in g['positive_face_coefficients']:co[fi]+=w*F(m)
            terms=sorted((fi,x)for fi,x in co.items()if x)
            loads=[sum((a.get(fi,F())*x for fi,x in terms),F())for a in alphas]
            assert loads[k]==value
            assert all(loads[j]<=F(ccs[j]['optimal_mass'])for j in range(2))
            key=(v,tuple(terms))
            if key not in columnkeys and any(loads[j]>mus[j].get(v,F())for j in range(2)):
                columnkeys.add(key);columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',
                    positive_face_coefficients=[[fi,str(x)]for fi,x in terms],orbit_weights=c['weights'],
                    charge_domain=labels[k],loads=list(map(str,loads))))
        row=dict(vertex=v,forest=u['forest'],domain_record_index=p['affected'].index(v),
            compatibility_edges=graph.number_of_edges(),maximal_literal_cliques=ncliques,
            distinct_count_vectors=npatterns,maximal_count_vectors=vectors,literal_representatives=reps,
            demands=u['demands'],certificates=ccs,baseline_optimal_masses=baseline_optimal,
            baseline_mu=[str(m.get(v,F()))for m in mus])
        rows.append(row)
        counts.update(forests=1,literal_stars=len(owner),corner_orbits=len(groups),
            actual_corners=sum(len(g['origins'])for g in groups),maximal_literal_cliques=ncliques,
            distinct_count_vectors=npatterns,maximal_count_vectors=len(vectors),
            full_mixture_terms=len(ccs[0]['mixture']),restricted_mixture_terms=len(ccs[1]['mixture']),
            full_mu_increases=int(newmus[0][v]>mus[0].get(v,F())),
            restricted_mu_increases=int(newmus[1][v]>mus[1].get(v,F())))
        if not verify:PROGRESS.write_text(json.dumps(dict(binding=binding,forest_records=normalize(rows)),sort_keys=True,separators=(',',':'))+'\n')
        print('FOREST',step,'/75',v,'stars',len(owner),'cliques',ncliques,'vectors',len(vectors),
            'raises',newmus[0][v]>mus[0].get(v,F()),newmus[1][v]>mus[1].get(v,F()),flush=True)
    if verify:assert opened==previous['open_forests']
    assert set(r['vertex']for r in rows)|{r['vertex']for r in opened}==set(p['affected'])
    assert not(set(r['vertex']for r in rows)&{r['vertex']for r in opened})
    repairs=[]
    for k,(base,newm)in enumerate(zip(p['patch']['repairs'],newmus)):
        assert all(newm.get(v,F())>=x for v,x in mus[k].items())
        objective=F(base['objective_numerator']);normalizer=sum(newm.values(),F())
        for col in columns:
            load=sum((alphas[k].get(fi,F())*F(x)for fi,x in col['positive_face_coefficients']),F())
            assert load==F(col['loads'][k]) and load<=newm.get(col['vertex'],F())
        repairs.append(dict(domain=labels[k],baseline_bound=base['bound'],objective_numerator=str(objective),
            normalizer=str(normalizer),bound=str(objective/normalizer),
            vertex_updates=[[v,str(mus[k].get(v,F())),str(x)]for v,x in sorted(newm.items())if x>mus[k].get(v,F())],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(newm.items())if x]))
    out=dict(status='PASS_COMPLETE_HISTORICAL22_CHANGED_CORNER_CONES'if not opened else'PARTIAL_HISTORICAL22_OPEN_FORESTS',
        sources=sources,domain_source=D.TAG,affected_forests=p['affected'],inherited_forests=inherited,
        complete_represented_forests=len(inherited)+len(rows),combined_positive_classes=270,
        counts=dict(counts),forest_records=rows,open_forests=opened,columns=columns,repairs=repairs,
        enumeration_guards=dict(default_maximal_cliques=W.MAX_CLIQUES,forest143_maximal_cliques=3000000,seconds_per_forest=30),
        face_multipliers=domain['face_multipliers'],restricted_face_multipliers=domain['restricted_face_multipliers'],
        theorem='Complete maximal literal cliques yield all maximal orbit-count vectors. Nonnegative charge vectors and equal-value nonnegative covering mixtures certify the full orbit-weight cone at each fixed demand. Actual simultaneous stars form a graph clique; covering cliques may be physically unrealizable.',
        preservation='Same5053 full face masses and same named restricted masses; vertex multipliers only increase from the frozen unowned-pair successor. All439410 inherited finite-model columns, including96176 zero-budget and52 positive-cost signed columns, and every prior weighted capacity stay valid. All478 unchanged conservative graph cones inherit their same-demand frozen certificates.',
        scope='Only75 changed forests recomputed. Historical first-owner outcomes, the explicit tested-pair patch including unknowns, and then added historical domains retain the frozen diagnostic precedence. No new wall-pair arithmetic, global optimization, complete physical realizability, old-alpha simplex/joint-cover claim, or dimension theorem.')
    out=normalize(out);path=ROOT/f'results/{TAG}.json'
    if verify:assert out==previous
    else:path.write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'columns',len(columns),'bounds',[float(F(r['bound']))for r in repairs],
        'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');a.add_argument('--resume',action='store_true')
    x=a.parse_args();assert not(x.verify and x.resume);run(x.verify,x.resume)
