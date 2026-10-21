#!/usr/bin/env python3
"""New-mass support inventory and a bounded exact raw76 corner preflight.

No compatibility pair, clique, local LP or global model is recomputed.
The singleton-cover ceiling is a possible-benefit bound, not an actual repair.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_current324_new5073_local_covers_20261018 as L
import fresh_raw125_raw205_actual_corner_frontier_20261017 as A

ROOT=L.ROOT; M=L.M; R=A.R; C=A.C; P=A.P
TAG='fresh_new5073_raw76_support_frontier_20261018'
FACES=(43028,43066)


def normal(x):return json.loads(json.dumps(x,sort_keys=True))


class Inputs:
    def __init__(self):self.paths=set()
    def read(self,tag):
        p=ROOT/f'results/{tag}.json';d=json.loads(p.read_text())
        assert d['status'].startswith('PASS'),(tag,d['status'])
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):
            q=f'{folder}/{tag}.{ext}'
            if (ROOT/q).exists():self.paths.add(q)
        return d


def support(inp,local):
    model=inp.read(L.MODEL);old=inp.read(M.MODEL)
    cert=next(c for c in model['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    oldcert=next(c for c in old['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={i:F(x) for i,x in cert['face_multipliers']}
    oa={i:F(x) for i,x in oldcert['face_multipliers']}
    represented=set(local['represented_classes'])
    assert len(alpha)==5073 and len(oa)==5121 and len(represented)==324
    assert cert['face_multipliers']==local['face_multipliers']
    fresh=set(alpha)-set(oa);missing=fresh-represented
    assert len(fresh)==1002 and len(missing)==982
    faces=list(inp.read(M.FORESTS)['faces']);forced=inp.read(M.FORCED)
    assert len(faces)==forced['inherited_counts']['faces']==78335
    for r in forced['new_faces']:
        assert r['index']==len(faces);faces.append(r)
    assert len(faces)==82530
    provenance=inp.read(M.ACTUAL)['provenance_masses']
    quads=inp.read(M.QUAD);cores=quads['canonical_cores']
    pmap={f:j for j,f in model['prepared']['provenance_to_face']}
    assert len(pmap)==len(model['prepared']['provenance_to_face'])
    historical_order=list(dict.fromkeys(inp.read(M.H)['source_precedence']+list(M.EXTRA)+[M.RAW125]))
    historical=defaultdict(list)
    # Read and discard each historical document; do not retain all old graphs.
    for tag in historical_order:
        d=inp.read(tag)
        for r in d['selected_faces']:
            if r['face'] in missing:historical[r['face']].append(dict(source=tag,record=r))
        del d
    rows=[]
    for fi in sorted(missing):
        mass=alpha[fi]
        if fi in pmap:
            pi=pmap[fi];p=provenance[pi];ci=p['canonical_prism_core'];ctx=p['source_context']
            r=dict(face=fi,kind='unambiguous_Q',provenance_kind='explicit',provenance=pi,
                   source_context=ctx,canonical_core=ci,raw_family=cores[ci]['source_core'],required=cores[ci]['required'])
        else:
            assert fi<len(faces);ff=faces[fi]
            if ff['kind']=='canonical_dotted_quad_aggregate':
                contexts=ff['source_contexts']
                if len(contexts)==1:
                    ctx,ci=contexts[0]
                    r=dict(face=fi,kind='unambiguous_Q',provenance_kind='unique_saved_context',provenance=None,
                           source_context=ctx,canonical_core=ci,raw_family=cores[ci]['source_core'],required=ff['required'])
                else:r=dict(face=fi,kind='merged_Q',source_contexts=contexts,required=ff['required'])
            else:
                assert 'simplex' in ff['kind'] or 'edge' in ff['kind']
                r=dict(face=fi,kind=ff['kind'],required=ff['required'])
        r.update(mass=mass,demand_mass=mass*r['required'],saved_historical_sources=[])
        for h in historical.get(fi,[]):
            assert r['kind']=='unambiguous_Q';z=h['record']
            assert (z['source_context'],z['canonical_prism_core'],z['raw_family'],z['provenance'])==(
                r['source_context'],r['canonical_core'],r['raw_family'],r['provenance'])
            r['saved_historical_sources'].append(h['source'])
        rows.append(r)
    counts=Counter(r['kind'] for r in rows)
    assert counts==dict(unambiguous_Q=425,merged_Q=123,maximal_simplex_aggregate=245,maximal_edge_incidence=189)
    ranked=sorted((r for r in rows if r['kind']=='unambiguous_Q'),key=lambda r:(-r['demand_mass'],r['face']))
    assert tuple(r['face'] for r in ranked[:2])==FACES
    selected=ranked[:2]
    assert all(r['raw_family']==76 and r['canonical_core']==65 and not r['saved_historical_sources'] for r in selected)
    return dict(counts=dict(newly_positive=1002,newly_positive_outside324=982,
                    outside_kinds=dict(counts),positive_represented=len(set(alpha)&represented),
                    newly_positive_saved_historical=sum(bool(r['saved_historical_sources']) for r in rows)),
                historical_scan_order=historical_order,rows=rows,unambiguous_ranking=[r['face'] for r in ranked],
                selected_faces=selected),quads,alpha


def geometry(inp,selected,quads,alpha):
    raw=inp.read(A.M.TAGS['raw']);spectators=inp.read(A.B.Q.SPECTATORS);forests=inp.read(P.FORESTS)
    models=dict(C.models());finite={C.code(g):n for n,g in models.items()}
    lookup={tuple(v['forest']):v['index'] for v in forests['vertices']}
    names={v['index']:v['forest'] for v in forests['vertices']}
    core=quads['canonical_cores'][65];graph=A.V.raw_graph(raw['records'][76])
    assert (core['source_core'],core['nodes'],core['required'])==(76,7,1)
    assert len(graph)==7 and raw['records'][76]['spectator_rank']==8
    corners=A.X.quad_features(dict(quadrilateral_black_partitions=core['quad_partitions']))
    assert len(corners)==len({tuple(c['vertex_omitted']) for c in corners})==4
    targets={};inverse={};stars=defaultdict(dict);occurrences=defaultdict(list);records=[]
    nmap=norientation=0
    for q in selected:
        fi=q['face'];face=forests['faces'][fi]
        assert face['kind']=='canonical_dotted_quad_aggregate'
        assert q['provenance'] is None and face['source_contexts']==[[q['source_context'],65]]
        context=spectators['contexts'][q['source_context']]
        assert quads['source_core_maps'][context['core']]['canonical_core']==65
        spectator=context['spectator_forest'];assert sum(len(models[n]) for n in spectator)==8
        cs=[]
        for c in corners:
            omitted=tuple(c['vertex_omitted']);ids=[i for i in range(7) if i not in omitted]
            h={i:{j:graph[i][j][0] for j in ids if i!=j and graph[i][j]!=A.RIGHT} for i in ids}
            assert all(m in (3,5) for row in h.values() for m in row.values())
            vf=tuple(sorted(spectator+[finite[C.code(a)] for a in C.pieces(h)]));v=lookup[vf]
            if v not in targets:
                targets[v]=P.forest_graph(vf,models);assert len(targets[v])==13
                inverse[v]=R.projection_data(targets[v])
            maps=R.active_embeddings(h,targets[v]);orientations={}
            for mp in maps:
                reverse={j:i for i,j in mp.items()}
                full=tuple(sorted(tuple(graph[x][reverse[j]] if j in reverse else A.RIGHT for j in range(13)) for x in omitted))
                assert full[0]!=full[1]
                occupied=tuple(sorted(tuple(sorted(mp[z] for z in pair)) for pair in c['occupied_pairs']))
                assert len(occupied)==1
                if full in orientations:assert orientations[full][0]==occupied
                else:orientations[full]=(occupied,sorted(mp.items()))
            assert orientations
            o=dict(face=fi,provenance=None,raw_family=76,omitted=omitted,vertex=v,
                   active_embeddings=len(maps),literal_orientations=len(orientations))
            nmap+=len(maps);norientation+=len(orientations);indices=[]
            J,active,blocks=inverse[v]
            for full,(occupied,mp) in sorted(orientations.items()):
                cov=tuple(R.sparse_covector(es) for es in full)
                assert all(i in active for p in cov for i,e in p)
                native=graph[omitted[0]][omitted[1]];field=A.pair_field(cov,(native,));heights=[]
                for p in cov:
                    norm=R.projection(p,p,J,field);assert not any(norm[1:])
                    height=norm[0]-2;assert height.sign()>0;heights.append(height)
                projection=R.projection(cov[0],cov[1],J,field);prescribed=field.negative_contact(native)
                assert field.signs_difference(tuple(a-b for a,b in zip(projection,prescribed)),heights[0]*heights[1])==0
                key=A.H.frozen(full)
                data=normal(dict(contacts=[[[i,list(e)] for i,e in p] for p in cov],occupied_pairs=occupied,
                    squared_timelike_heights=[x.out() for x in heights],native_omitted_product=list(native)))
                if key in stars[v]:assert stars[v][key]['native']==data
                else:stars[v][key]=dict(index=len(stars[v]),native=data,active_embedding=mp,
                    pair_field_radicals=[field.tau.out(),field.sigma.out()],projection=R.qout(projection))
                indices.append(stars[v][key]['index'])
            assert len(indices)==len(set(indices))==o['literal_orientations']
            cs.append(o);occurrences[v].append(dict(origin=o,literal_star_indices=sorted(indices)))
        records.append(dict(**q,spectator_forest=spectator,corners=cs))
    rows=[dict(vertex=v,forest=names[v],target=targets[v],active_projection_blocks=inverse[v][2],
               occurrences=occurrences[v],literal_stars=list(stars[v].values())) for v in sorted(targets)]
    return dict(selected_faces=records,forest_records=rows,counts=dict(faces=2,actual_corner_occurrences=8,
        vertex_types=len(rows),active_embeddings=nmap,literal_orientation_occurrences=norientation,
        distinct_native_stars=sum(len(z) for z in stars.values())))


def inventory(domain,local,alpha,inp):
    mu={v:F(x) for v,x in local['vertex_multipliers']};old={r['vertex']:r for r in local['forest_records']}
    num,den=F(local['objective_numerator']),F(local['normalizer']);assert sum(mu.values(),F())==den
    rows=[];increase=F();total=F();geometry_docs={}
    def read_geometry(tag):
        if tag not in geometry_docs:geometry_docs[tag]=inp.read(tag)
        return geometry_docs[tag]
    for raw in domain['forest_records']:
        v=raw['vertex'];r=old.get(v);cover=F(r['exact_local_optimum']) if r else F()
        oldmu=mu.get(v,F());assert cover<=oldmu
        extra=sum((alpha[o['origin']['face']] for o in raw['occurrences']),F())
        # Native merging with old nodes is deliberately not assumed. This is
        # an upper bound, since duplicate native stars only reduce workload.
        oldn=0
        if r:
            d=read_geometry(r['geometry_source'])
            if 'forest_records' in d:g=next(z for z in d['forest_records'] if z['vertex']==v)
            else:g=d;assert g['vertex']==v
            groups=(g.get('domain') or {}).get('corner_orbits') or g.get('corner_orbits') or (g.get('union') or {}).get('corner_orbits')
            if groups is None:
                gd=read_geometry(r['group_source']);gr=next(z for z in gd['forest_records'] if z['vertex']==v)
                groups=gr['corner_orbits']
            oldn=sum(len(z['literal_stars']) for z in groups)
            assert len(groups)==len(r['demands'])
        n=len(raw['literal_stars']);upper=max(oldmu,cover+extra);gain=upper-oldmu
        rows.append(dict(vertex=v,forest=raw['forest'],existing=v in old,old_literal_stars=oldn,
            candidate_native_stars=n,new_pair_positions_upper=oldn*n+n*(n-1)//2,
            added_actual_corners=len(raw['occurrences']),extra_singleton_mass=extra,
            old_local_cover=cover,old_vertex_mass=oldmu,maximum_mu_increase=gain,
            hypothetical_vertex_mass_upper=upper))
        increase+=gain;total+=extra
    assert total==4*sum((alpha[f] for f in FACES),F())
    return dict(forests=rows,total_extra_singleton_mass=total,maximum_mu_increase=increase,
        candidate_new_pair_positions_upper=sum(r['new_pair_positions_upper'] for r in rows),
        numerator=num,baseline_normalizer=den,baseline_ratio=num/den,
        singleton_cover_ratio_lower_bound=num/(den+increase),required_increment_to_reach12=num/12-den,
        claim='Old attaining covers plus one singleton mass per added actual corner form a cover for every extension that preserves all old-old outcomes and only removes pairs involving new stars. Even overlapping identical native orbits are safely overcovered. This is a ceiling on possible monotone multiplier increases, not a new graph or actual repaired witness.')


def run():
    inp=Inputs();local=inp.read(L.TAG);sup,quads,alpha=support(inp,local)
    domain=geometry(inp,sup['selected_faces'],quads,alpha);inv=inventory(domain,local,alpha,inp)
    paths=set(inp.paths)
    for mod in (L,M,A,R,C,P,A.M,A.V,A.X,A.B,R.T):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    return normal(M.enc(dict(status='PASS_NEW5073_SUPPORT_AND_RAW76_TWO_CLASS_PREFLIGHT',
        sources={p:M.digest(p) for p in sorted(paths)},support=sup,domain=domain,inventory=inv,
        current_represented_classes=324,current_represented_forests=600,
        pair_tests=0,clique_enumerations=0,local_LPs=0,global_optimizer_invocations=0,
        current_graphs_and_owners_changed=False,new_dimension_bound=None,
        scope='Complete saved newly-positive support classification, bounded historical selection scan, and two selected raw76 full corner domains with exact native projections. Saved old local optima supply a singleton-cover benefit ceiling. No native merge with old graphs, compatibility exclusions, new covers at enlarged demands, changed ratio, physical realization or dimension theorem.')))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
    p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['domain']['counts'],'pair_upper',out['inventory']['candidate_new_pair_positions_upper'],
          'ratio_lower',float(F(out['inventory']['singleton_cover_ratio_lower_bound'])),flush=True)
