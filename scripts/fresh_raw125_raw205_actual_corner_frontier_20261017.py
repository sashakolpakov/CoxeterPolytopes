#!/usr/bin/env python3
"""New actual raw125/raw205 corners at the unchanged5053-face witness.

Every contact is reconstructed in one complete labelled finite-vertex frame.
Historical and explicit owner outcomes remain fixed, including unknowns.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_unowned_pair_optimistic_frontier_20261017 as O
import fresh_raw110_mixed_actual_corners_20261017 as R

ROOT=O.ROOT;D=O.D;H=O.H;W=O.W;P=R.P
M,V,K,C,X,B=R.M,R.V,R.K,R.C,R.X,R.B
TAG='fresh_raw125_raw205_actual_corner_frontier_20261017'
BASE=O.BASE
PROGRESS=Path('/private/tmp')/(TAG+'_progress.json')
RIGHT=R.RIGHT
FAMILIES={125:dict(core=105,order=8,dimension=3,tau=K(13,5)),
          205:dict(core=176,order=9,dimension=3,tau=K(6,2))}


def normal(x):return json.loads(json.dumps(x,sort_keys=True))
def native(s):return tuple(tuple((i,tuple(e))for i,e in p)for p in s['contacts'])


@lru_cache(None)
def pair_field(covectors,extra_contacts=()):
    # Native squared magnitudes choose the squareclasses. No vector encoded
    # in a previous pair's field is ever transported to this field.
    ts=sorted({(t.a,t.b)for e in [e for p in covectors for i,e in p]+list(extra_contacts)
               if (t:=V.squared(e)) and M.ksqrt(t)is None})
    selected=[]
    for a,b in ts:
        t=K(a,b)
        if all(M.ksqrt(t/u)is None for u in selected):selected.append(t)
    assert len(selected)<=2,'more than two native squareclasses in one wall pair'
    for t in (R.TAUS[23],R.TAUS[24],R.TAUS[39],K(2)):
        if len(selected)==2:break
        assert M.ksqrt(t)is None
        if all(M.ksqrt(t/u)is None for u in selected):selected.append(t)
    assert len(selected)==2
    selected.sort(key=lambda t:(t.a,t.b))
    field=R.canonical_field(tuple((t.a,t.b)for t in selected))
    for e in [e for p in covectors for i,e in p]+list(extra_contacts):field.negative_contact(e)
    return field


def prepare():
    p=D.prepare();base=W.read(BASE);cones=W.read(O.C.TAG);domain=W.read(D.TAG);eight=W.read(O.EIGHT)
    assert base['face_multipliers']==domain['face_multipliers']==cones['face_multipliers']
    assert base['restricted_face_multipliers']==domain['restricted_face_multipliers']
    rows=dict(p['original_rows']);rows.update(p['patch_rows']);rows.update({r['vertex']:r for r in cones['forest_records']})
    rows[eight['vertex']]=eight;rows.update({r['vertex']:r for r in base['forest_records']})
    assert len(rows)==553
    later=defaultdict(list);later[eight['vertex']].append(eight)
    for r in base['forest_records']:later[r['vertex']].append(r)
    p.update(baseline=base,baseline_rows=rows,expanded_domains={r['vertex']:r for r in domain['forest_records']},later=later)
    return p


def reconstruct(p):
    support=W.read(H.SUPPORT);model=W.read(H.MODEL)
    actual=P.read(M.TAGS['actual']);quads=P.read(M.TAGS['quad'])
    raw=P.read(M.TAGS['raw']);spectators=P.read(B.Q.SPECTATORS);old=P.read(P.FORESTS)
    models=dict(C.models());finite={C.code(g):n for n,g in models.items()}
    lookup={tuple(v['forest']):v['index']for v in old['vertices']}
    p['forest_names']={v['index']:v['forest']for v in old['vertices']}
    inverse={fi:j for j,fi in model['prepared']['provenance_to_face']}
    assert len(inverse)==len(model['prepared']['provenance_to_face'])
    wanted={r['face']:r for r in support['quad_rows']if r.get('raw_family')in FAMILIES and r['provenance_kind']!='merged_unexpanded'}
    assert Counter(r['raw_family']for r in wanted.values())=={125:9,205:11}
    assert not(set(wanted)&set(p['selected']))
    assert FAMILIES[125]['tau']==R.TAUS[23]*R.TAUS[295]
    assert M.ksqrt(FAMILIES[125]['tau'])is None
    assert M.ksqrt(FAMILIES[205]['tau'])==K(1,1)
    selected=[];occurrences=defaultdict(list);targets={}
    for fi,q in sorted(wanted.items()):
        assert fi in p['alpha']and F(q['mass'])==p['alpha'][fi]
        family=q['raw_family'];meta=FAMILIES[family]
        if q['provenance_kind']=='explicit':
            j=inverse[fi];assert j==q['provenance']
            pr=actual['provenance_masses'][j]
            assert pr['source_context']==q['source_context']and pr['canonical_prism_core']==q['canonical_core']
            assert pr['prism_dimension']==3
            provenance_kind='existing_explicit_map';contexts=None
        else:
            assert fi not in inverse and q['source_context_count']==1
            source=old['faces'][fi]
            assert source['kind']=='canonical_dotted_quad_aggregate'
            assert source['source_contexts']==[[q['source_context'],q['canonical_core']]]
            j=None;contexts=source['source_contexts'];provenance_kind='unique_saved_context'
        cc=quads['canonical_cores'][q['canonical_core']]
        assert cc['index']==meta['core']and cc['source_core']==family and cc['nodes']==meta['order']
        assert q['required']==cc['required']==3
        ctx=spectators['contexts'][q['source_context']]
        assert quads['source_core_maps'][ctx['core']]['canonical_core']==cc['index']
        sp=ctx['spectator_forest'];sr=sum(len(models[name])for name in sp)
        assert sr+meta['order']==15
        rg=V.raw_graph(raw['records'][family]);assert len(rg)==meta['order']
        assert all(e[0]in(0,2,3,5)for row in rg for e in row)
        assert all(V.squared(e)==meta['tau']for row in rg for e in row if e[0]==0)
        corners=X.quad_features(dict(quadrilateral_black_partitions=cc['quad_partitions']))
        assert len(corners)==len({tuple(c['vertex_omitted'])for c in corners})==6
        records=[]
        for c in corners:
            omitted=tuple(c['vertex_omitted']);ids=[i for i in range(len(rg))if i not in omitted]
            h={i:{j:rg[i][j][0]for j in ids if j!=i and rg[i][j]!=RIGHT}for i in ids}
            assert all(m in(3,5)for row in h.values()for m in row.values())
            vf=tuple(sorted(sp+[finite[C.code(a)]for a in C.pieces(h)]));vi=lookup[vf]
            if vi not in targets:targets[vi]=P.forest_graph(vf,models)
            assert len(targets[vi])==13
            maps=R.active_embeddings(h,targets[vi]);oris={}
            for mp in maps:
                es=[]
                for x in omitted:
                    es.append(tuple(rg[x][next(i for i,z in mp.items()if z==j)]if j in mp.values()else RIGHT for j in range(13)))
                key=tuple(sorted(es));assert key[0]!=key[1]
                occupied=tuple(sorted(tuple(sorted(mp[z]for z in pair))for pair in c['occupied_pairs']))
                assert len(occupied)==len(set(occupied))==2
                if key in oris:assert oris[key][0]==occupied
                else:oris[key]=(occupied,sorted(mp.items()))
            assert oris
            rec=dict(face=fi,provenance=j,raw_family=family,omitted=omitted,vertex=vi,
                     active_embeddings=len(maps),literal_orientations=len(oris))
            records.append(rec);occurrences[vi].append((rec,oris,rg[omitted[0]][omitted[1]]))
        selected.append(dict(face=fi,provenance=j,raw_family=family,canonical_prism_core=cc['index'],
            source_context=q['source_context'],provenance_kind=provenance_kind,source_contexts=contexts,
            raw_core_order=meta['order'],spectator_rank=sr,spectator_forest=sp,corners=records))
    assert len(selected)==20 and sum(len(r['corners'])for r in selected)==120
    return selected,occurrences,targets


def merge(v,orecs,target,p):
    if v in p['baseline_rows']:
        u,graph,owner,unowned,demand=O.current(v,p,p['expanded_domains'],p['baseline_rows'],p['later'])
        nodes=[normal(W.forests(tag)[v]['literal_stars'][i])for tag,i in u['literal_references']]
        groups=normal(u['corner_orbits']);oldcount=len(nodes)
    else:
        u=None;graph=nx.Graph();owner={};nodes=[];groups=[];oldcount=0
    keys={H.frozen(s['contacts']):i for i,s in enumerate(nodes)}
    assert len(keys)==len(nodes)
    J,active,components=R.projection_data(target)
    added_orbits={};new_checks=[]
    for rec,oris,native_product in orecs:
        indices=[]
        for full,(occupied,mp)in sorted(oris.items()):
            covectors=tuple(R.sparse_covector(q)for q in full)
            contacts=normal([[[i,list(e)]for i,e in q]for q in covectors]);key=H.frozen(contacts)
            assert all(i in active for q in covectors for i,e in q)
            field=pair_field(covectors,(native_product,));heights=[]
            for q in covectors:
                z=R.projection(q,q,J,field);assert not any(z[1:]),(v,rec,'height outside K')
                height=z[0]-2;assert height.sign()>0;heights.append(height)
            projection=R.projection(covectors[0],covectors[1],J,field)
            prescribed=field.negative_contact(native_product)
            assert field.signs_difference(tuple(a-b for a,b in zip(projection,prescribed)),heights[0]*heights[1])==0
            star=dict(index=len(nodes),contacts=contacts,occupied_pairs=normal(occupied),
                active_embedding=mp,origin=normal(rec),squared_timelike_heights=[z.out()for z in heights],
                native_omitted_product=list(native_product))
            if key in keys:
                i=keys[key]
                assert i>=oldcount,'new actual class duplicates an old complete native star'
                assert all(star[k]==nodes[i][k]for k in H.NATIVE_FIELDS)
            else:
                i=len(nodes);keys[key]=i;star['index']=i;nodes.append(star);graph.add_node(i)
                new_checks.append(dict(index=i,pair_field_radicals=[field.tau.out(),field.sigma.out()],
                    projection=R.qout(projection),squared_timelike_heights=[z.out()for z in heights]))
            indices.append(i)
        orbit=tuple(sorted(indices));assert len(orbit)==len(set(orbit))==rec['literal_orientations']
        if orbit not in added_orbits:
            assert not any(i in owner for i in orbit),'partially overlapping complete native orbits'
            gi=len(groups)
            for i in orbit:owner[i]=gi
            groups.append(dict(index=gi,face=rec['face'],literal_stars=list(orbit),origins=[],positive_face_coefficients=[]))
            added_orbits[orbit]=gi
        gi=added_orbits[orbit];g=groups[gi]
        assert g['face']==rec['face']and H.origin_key(normal(rec))not in{H.origin_key(c)for c in g['origins']}
        g['origins'].append(normal(rec));g['positive_face_coefficients']=[[rec['face'],len(g['origins'])]]
    assert set(owner)==set(range(len(nodes)))
    candidates=[(i,j)for i,j in combinations(range(len(nodes)),2)if j>=oldcount]
    # A new native star cannot belong to an older full domain, even if that
    # old producer omitted its class from the current positive selection.
    for tag in p['order']+list(D.EXTRA):
        if v in W.forests(tag):
            oldkeys={H.frozen(s['contacts'])for s in W.forests(tag)[v]['literal_stars']}
            assert not(oldkeys&{H.frozen(s['contacts'])for s in nodes[oldcount:]})
    demand=[[str(sum((a.get(fi,F())*F(m)for fi,m in g['positive_face_coefficients']),F()))for g in groups]
            for a in(p['alpha'],p['ralpha'])]
    return dict(vertex=v,forest=p['forest_names'][v],
        inherited_union=u,literal_stars=nodes,corner_orbits=groups,native_checks=new_checks,demands=demand,
        old_literal_stars=oldcount,new_literal_stars=len(nodes)-oldcount,new_pair_domain=candidates,
        counts=dict(literal_stars=len(nodes),new_literal_stars=len(nodes)-oldcount,new_pair_tests=len(candidates),
                    corner_orbits=len(groups),new_actual_corners=len(orecs))),graph,owner,J


def physical_pairs(row,graph,J):
    nodes=row['literal_stars'];pp=[native(s)for s in nodes]
    heights={};hs=[]
    for i,star in enumerate(nodes):
        hh=[]
        for q in pp[i]:
            if q not in heights:
                field=pair_field((q,));z=R.projection(q,q,J,field)
                assert not any(z[1:]);h=z[0]-2;assert h.sign()>0
                heights[q]=h
            hh.append(heights[q])
        assert [h.out()for h in hh]==star['squared_timelike_heights']
        hs.append(hh)
    cache={};exprids={};expressions=[];tested=[];excluded=[];reasons=Counter()
    oldedges=set(map(lambda e:tuple(sorted(e)),graph.edges()))
    for i,j in row['new_pair_domain']:
        assert not graph.has_edge(i,j)
        shared=set(map(tuple,nodes[i]['occupied_pairs']))&set(map(tuple,nodes[j]['occupied_pairs']))
        why=[i,j,0]if shared else None
        if why is None:
            for a,b in((0,0),(0,1),(1,0),(1,1)):
                pi,pj=pp[i][a],pp[j][b]
                if pi==pj:continue
                key=tuple(sorted((pi,pj)))
                if key not in cache:
                    field=pair_field((pi,pj));projection=R.projection(pi,pj,J,field);h=hs[i][a]*hs[j][b]
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
        if why is None:graph.add_edge(i,j);tested.append([i,j,None])
        else:excluded.append(why);tested.append([i,j,why[2]]);reasons[why[2]]+=1
    assert {e for e in map(lambda e:tuple(sorted(e)),graph.edges())if e[1]<row['old_literal_stars']}==oldedges
    row.update(newly_owned_pair_domain=tested,incompatibility_certificates=excluded,
        wall_expressions=expressions,reason_counts=dict(reasons))
    row['counts'].update(new_exclusions=len(excluded),new_unknown_pairs=len(tested)-len(excluded),
                         distinct_wall_expressions=len(expressions),recomputed_native_walls=len(heights))
    del row['new_pair_domain']


def restricted_certificate(row,vectors,p):
    v=row['vertex'];groups=row['corner_orbits'];demand=list(map(F,row['demands'][1]))
    if v not in p['baseline_rows']:
        assert not any(demand)
        return dict(weights=['0']*len(groups),mixture=[],optimal_mass='0'),[]
    old=p['baseline_rows'][v];c=old['certificates'][1]
    n=len(c['weights']);assert all(x==0 for x in demand[n:])
    mix=Counter();transport=[]
    for i,amount in c['mixture']:
        vec=old['maximal_count_vectors'][i]+[0]*(len(groups)-n)
        choices=[j for j,w in enumerate(vectors)if all(a>=b for a,b in zip(w,vec))]
        assert choices;chosen=min(choices);mix[chosen]+=F(amount)
        transport.append([i,chosen,amount])
    out=dict(weights=c['weights']+['0']*(len(groups)-n),
        mixture=[[i,str(x)]for i,x in sorted(mix.items())],optimal_mass=c['optimal_mass'])
    W.check_certificate(vectors,demand,out)
    return out,transport


def complete_cone(row,graph,owner,p,saved):
    if saved is not None and saved['cone_status']=='open_guard':
        row.update(cone_status='open_guard',cone_limitation=saved['cone_limitation']);return
    enum,error=O.C.vectors_for(row['vertex'],graph,owner,len(row['corner_orbits']))
    if error:
        assert saved is None
        row.update(cone_status='open_guard',cone_limitation='Complete clique enumeration exceeded its guard; no partial vector domain or capacity is promoted.');return
    vectors,reps,ncliques,npatterns=enum;vectors=list(map(list,vectors));reps=list(map(list,reps))
    for vec,clique in zip(vectors,reps):
        assert all(graph.has_edge(i,j)for i,j in combinations(clique,2))
        assert vec==[Counter(owner[i]for i in clique).get(j,0)for j in range(len(row['corner_orbits']))]
    if saved is not None:
        assert saved['maximal_count_vectors']==vectors and saved['literal_representatives']==reps
        full=saved['certificates'][0]
    else:full=H.exact_profile(vectors,list(map(F,row['demands'][0])))
    restricted,transport=restricted_certificate(row,vectors,p)
    certificates=[full,restricted]
    for k,c in enumerate(certificates):W.check_certificate(vectors,list(map(F,row['demands'][k])),c)
    row.update(cone_status='complete',maximal_literal_cliques=ncliques,distinct_count_vectors=npatterns,
        maximal_count_vectors=vectors,literal_representatives=reps,certificates=certificates,
        restricted_cover_transport=transport)


def bindings():
    tags=(H.MODEL,H.SUPPORT,D.BASE,D.PATCH,D.TAG,O.C.TAG,O.EIGHT,BASE)
    sources={}
    for tag in tags:
        doc=W.read(tag)
        for path,h in doc.get('sources',doc.get('prepared',{}).get('sources',{})).items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for tag in(M.TAGS['actual'],M.TAGS['quad'],M.TAGS['raw'],B.Q.SPECTATORS,P.FORESTS):
        path=f'results/{tag}.json';sources[path]=W.digest(path)
    for mod in(D,H,W,O,R,R.T,P,M,V,C,X,B):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    return sources


def main(verify=False,resume=False,prepare_only=False):
    start=time.monotonic();p=prepare();selected,occurrences,targets=reconstruct(p);sources=bindings()
    binding=dict(sources=sources,face_multipliers=p['baseline']['face_multipliers'])
    saved=W.read(TAG)if verify else None;cache={}
    if saved is not None:
        assert saved['sources']==sources;cache={r['vertex']:r for r in saved['forest_records']}
    elif resume and PROGRESS.exists():
        progress=json.loads(PROGRESS.read_text());assert progress['binding']==normal(binding)
        cache={r['vertex']:r for r in progress['forest_records']}
    print('NEW_DOMAIN',len(selected),'classes',sum(map(len,occurrences.values())),'corners',len(occurrences),'forests',flush=True)
    rows=[];counts=Counter()
    for v,orecs in sorted(occurrences.items()):
        if v in cache and not verify and not prepare_only:
            row=cache[v]
        else:
            row,g,owner,J=merge(v,orecs,targets[v],p)
            if not prepare_only:
                physical_pairs(row,g,J);complete_cone(row,g,owner,p,cache.get(v))
                if verify:assert normal(row)==cache[v],v
        rows.append(row);counts.update(row['counts'])
        print('FOREST',v,row['counts'],row.get('cone_status','domain only'),flush=True)
        if not verify and not prepare_only:
            PROGRESS.write_text(json.dumps(normal(dict(binding=binding,forest_records=rows)),separators=(',',':'))+'\n')
    if prepare_only:
        out=normal(dict(binding=binding,selected_faces=selected,forest_records=rows,counts=dict(counts)))
        Path('/private/tmp/raw125_raw205_domain_probe.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
        print('PASS_RAW125_RAW205_DOMAIN_PROBE',dict(counts),'seconds',round(time.monotonic()-start,2),flush=True);return
    oldmus=[{i:F(x)for i,x in r['vertex_multipliers']}for r in p['baseline']['repairs']]
    newmus=[dict(m)for m in oldmus];columns=[];seen=set();newfaces={r['face']for r in selected}
    alphas=[p['alpha'],p['ralpha']]
    for row in rows:
        if row['cone_status']!='complete':continue
        v=row['vertex']
        for k,c in enumerate(row['certificates']):newmus[k][v]=max(newmus[k].get(v,F()),F(c['optimal_mass']))
        c=row['certificates'][0];coef=Counter()
        for weight,g in zip(map(F,c['weights']),row['corner_orbits']):
            for fi,m in g['positive_face_coefficients']:coef[fi]+=weight*F(m)
        terms=sorted((fi,a)for fi,a in coef.items()if a);key=(v,tuple(terms))
        loads=[sum((a.get(fi,F())*x for fi,x in terms),F())for a in alphas]
        assert loads[0]==F(c['optimal_mass'])
        assert all(loads[k]<=F(row['certificates'][k]['optimal_mass'])for k in range(2))
        if any(fi in newfaces for fi,x in terms)and key not in seen:
            seen.add(key);columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',
                orbit_weights=c['weights'],positive_face_coefficients=[[fi,str(x)]for fi,x in terms],loads=list(map(str,loads))))
    repairs=[]
    for k,old in enumerate(p['baseline']['repairs']):
        assert all(newmus[k].get(v,F())>=x for v,x in oldmus[k].items())
        numerator=F(old['objective_numerator']);den=sum(newmus[k].values(),F())
        repairs.append(dict(domain=old['domain'],baseline_bound=old['bound'],objective_numerator=str(numerator),normalizer=str(den),
            bound=str(numerator/den),vertex_updates=[[v,str(oldmus[k].get(v,F())),str(x)]for v,x in sorted(newmus[k].items())if x>oldmus[k].get(v,F())],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(newmus[k].items())if x]))
    assert not repairs[1]['vertex_updates']and repairs[1]['bound']==p['baseline']['repairs'][1]['bound']
    oldvertices=set(p['baseline_rows']);changed=set(occurrences)
    out=normal(dict(status='PASS_RAW125_RAW205_NEW_ACTUAL_CORNERS_AND_EXACT_CAPACITIES',sources=sources,
        baseline_source=BASE,selected_faces=selected,added_classes=20,combined_positive_classes=290,
        affected_forests=sorted(changed),new_forests=sorted(changed-oldvertices),inherited_forests=len(oldvertices-changed),
        total_represented_forests=len(oldvertices|changed),forest_records=rows,columns=columns,repairs=repairs,
        counts=dict(counts),open_forests=[r['vertex']for r in rows if r['cone_status']!='complete'],
        face_multipliers=p['baseline']['face_multipliers'],restricted_face_multipliers=p['baseline']['restricted_face_multipliers'],
        preservation='Every old complete owner and explicit tested pair, including unknowns, is preserved. New tests involve at least one genuinely new native star absent from all older source domains. Fixed face masses and monotone vertex increases preserve439410 old model columns and earlier capacities.',
        scope='Twenty explicit/unique-provenance actual classes, complete labelled native corner orientations, immutable pairwise fields and necessary physical pair gates. Only affected conservative graph cones recomputed; open guards, if any, are explicit and promote no partial cone. Saved replay uses rational certificates without an optimizer. No full physical realization, global reoptimization, stale-alpha implicit coverage, finite-model optimum or dimension theorem.'))
    path=ROOT/f'results/{TAG}.json'
    if verify:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'columns',len(columns),'bounds',[float(F(r['bound']))for r in repairs],
        'open',out['open_forests'],'seconds',round(time.monotonic()-start,2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--resume',action='store_true');ap.add_argument('--prepare',action='store_true')
    args=ap.parse_args();main(args.verify,args.resume,args.prepare)
