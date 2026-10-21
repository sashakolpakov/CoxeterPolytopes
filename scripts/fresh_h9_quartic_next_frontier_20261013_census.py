#!/usr/bin/env python3
"""Complete near-corner segment census in shifted integral wall-pairing cosets."""
import argparse
import hashlib
import json
import time
from collections import Counter
from itertools import product

from fresh_h9_quartic_edge_walk_20261013_next import (
    ROOT, K, Z, O, T, ALPHA, NORMS, Host, enc, menc, pair, mv, dot,
    inverse, alphabet, complement, square_root, finite_vertex)
from fresh_h9_quartic_mixed_20261013 import balanced_alphabet
from fresh_h9_glued_host_20261013 import block_vectors
from fresh_h9_quartic_edge_walk_20261013_next import nullspace
from fractions import Fraction

BASE='results/fresh_h9_quartic_first_endpoint_20261013_completion.json'
PROPOSAL='results/fresh_h9_quartic_next_frontier_20261013.json'
OUT=ROOT/'results/fresh_h9_quartic_next_frontier_20261013_census.json'
STRUCTURE='results/fresh_h9_glued_host_20261013_structure.json'


def line_key(r):
    first=next(z for z in r if z)
    return tuple(z/first for z in r)


def reflection(r,e,g):
    t=2*pair(r,g,e)/pair(e,g,e)
    return tuple(x-t*y for x,y in zip(r,e))


def finite_orbit(roots,generators,g):
    found={line_key(r):r for r in roots};queue=list(found.values())
    for r in queue:
        for e in generators:
            s=reflection(r,e,g);key=line_key(s)
            if key not in found:
                found[key]=s;queue.append(s)
                assert len(found)<1000,'finite reflection orbit guard'
    return found


def root_test(r,q,g):
    return (all(z.is_integral for z in r) and
            not all((z/T).is_integral for z in r) and
            all((2*z/q).is_integral for z in mv(g,r)))


def run(edge_index=2,probe_fraction=None):
    start=time.monotonic();g=Host(1).g
    base_path,proposal_path=BASE,PROPOSAL
    output=OUT.with_name(OUT.stem+f'_edge{edge_index}.json')
    if probe_fraction is not None:output=output.with_name(output.stem+'_probe.json')
    base=json.loads((ROOT/base_path).read_text());proposal=json.loads((ROOT/proposal_path).read_text())
    er=proposal['edges'][edge_index];edge_ids=er['edge_indices'];start_wall=er['released_wall']
    walls=[tuple(map(K.from_json,r))for r in base['forced_H5_face']['all_walls']]
    inside=tuple(map(K.from_json,base['forced_H5_face']['strict_interior']))
    edge=[walls[i]for i in edge_ids];w=walls[start_wall]
    initial_vertex=finite_vertex(walls,g,inside,edge_ids+[start_wall])
    k0=tuple(map(K.from_json,initial_vertex['vector']));k0=tuple(z/k0[0]for z in k0)
    candidates=[(i,row)for i,row in enumerate(proposal['compatible_dictionary'])if row['edge_index']==edge_index]
    proposal_index,proposed=min(candidates,key=lambda t:(K.from_json(t[1]['edge_parameter']),K.from_json(t[1]['height'])))
    k1=tuple(map(K.from_json,proposed['point']))
    if probe_fraction is not None:
        lam=K(Fraction(probe_fraction));assert 0<lam<1
        k1=tuple(x+lam*(y-x)for x,y in zip(k0,k1))
    assert k0[0]==k1[0]==1 and all(pair(k,g,e)==0 for k in(k0,k1)for e in edge)
    assert pair(k1,g,w)<0
    h0=ALPHA/(-pair(k0,g,k0))-1;h1=ALPHA/(-pair(k1,g,k1))-1;H=max(h0,h1)
    eg=tuple(tuple(pair(r,g,s)for s in edge)for r in edge);ei=inverse(eg)
    ww=mv(ei,tuple(pair(w,g,e)for e in edge))
    residual_axis=tuple(w[j]-sum((c*e[j]for c,e in zip(ww,edge)),Z)for j in range(6))
    B=pair(residual_axis,g,residual_axis);N=pair(k0,g,k0)
    C=pair(k0,g,k1)**2/(N*pair(k1,g,k1))
    assert C>=1 and all(B.sign(e)>0 for e in range(4))
    assert pair(residual_axis,g,k0)==0 and all(pair(residual_axis,g,e)==0 for e in edge)
    census=[];profiles=[];folded={};span={};stats=Counter();bounds=[]
    for q in NORMS:
        alph=[alphabet((q*pair(e,g,e),)*4)for e in edge]
        qrows=[]
        for fixed in product(*alph):
            coeff=mv(ei,tuple(-z for z in fixed))
            r0=tuple(sum((c*e[j]for c,e in zip(coeff,edge)),Z)for j in range(6))
            residual=q-pair(r0,g,r0)
            if residual==0:
                if root_test(r0,q,g):span[line_key(r0)]=r0
                continue
            if any(residual.sign(e)<=0 for e in range(4)):continue
            profile_index=len(profiles)
            record={'norm':q.to_json(),'fixed_magnitudes':enc(fixed),
                    'fixed_component':enc(r0),'residual_norm':residual.to_json(),
                    'residual_square':square_root(residual) is not None}
            offset=pair(r0,g,w);caps=(residual*B*C,*(residual*B for _ in range(3)))
            rounded=K(*(round(Fraction(z))for z in offset.to_json()));shift=rounded-offset
            loose=tuple(2*cap+2*shift*shift for cap in caps)
            absolute=balanced_alphabet(loose)
            ps=sorted({rounded+z for u in absolute for z in(u,-u)
                       if all((caps[e]-(rounded+z-offset)**2).sign(e)>=0 for e in range(4))})
            record.update(pairing_offset=offset.to_json(),pairing_caps=enc(caps),
                          integral_center=rounded.to_json(),loose_caps=enc(loose),pairing_alphabet=enc(ps))
            profiles.append(record);qrows.append(profile_index)
            for value in ps:
                anchor=tuple(z+(value-offset)*bb/B for z,bb in zip(r0,residual_axis))
                assert pair(anchor,g,k0)==0 and pair(anchor,g,w)==value
                square=(q-pair(anchor,g,anchor))/N;stats['quadratics']+=1
                radical=square_root(square)
                if radical is None:continue
                stats['square_discriminants']+=1
                for t in {radical,-radical}:
                    r=tuple(a+t*d for a,d in zip(anchor,k0))
                    assert pair(r,g,w)==value and pair(r,g,r)==q
                    assert tuple(-pair(r,g,e)for e in edge)==fixed
                    if not root_test(r,q,g):continue
                    stats['primitive_integral_solutions']+=1
                    a,b=pair(k0,g,r),pair(k1,g,r)
                    if a*b>0:continue
                    assert a or b
                    lam=-a/(b-a)
                    assert 0<=lam<=1
                    key=line_key(r)
                    if key in folded:continue
                    folded[key]=r
                    census.append({'profile_index':profile_index,'norm':q.to_json(),
                                   'vector':enc(r),'start_pairing':a.to_json(),
                                   'end_pairing':b.to_json(),'segment_parameter':lam.to_json()})
        print({'norm':q.to_json(),'profiles':len(qrows),'shifted_pairings':sum(len(profiles[i]['pairing_alphabet'])for i in qrows),
               'folded_segment_mirrors':len(census),'seconds':round(time.monotonic()-start,2)},flush=True)
        record={'norm':q.to_json(),'pairing_alphabets':[enc(a)for a in alph],
                'profile_indices':qrows}
        # Each completed norm packet is preserved if a later computation is stopped.
        partial={'status':'RUNNING','completed_norms':len(bounds)+1,'last_norm':q.to_json(),
                 'profiles':profiles,'folded_segment_roots':census}
        output.with_name(output.stem+'_partial.json').write_text(json.dumps(partial,indent=2)+'\n')
        bounds.append(record)
    edge_mirrors=finite_orbit(edge,edge,g)
    complete_span=finite_orbit(list(span.values()),edge,g)
    assert set(complete_span)==set(edge_mirrors)
    expanded=finite_orbit(list(folded.values()),edge,g)
    initial={key:r for key,r in expanded.items()if pair(k0,g,r)==0}
    initial.update(edge_mirrors)
    initial_expected=finite_orbit(edge+[walls[start_wall]],edge+[walls[start_wall]],g)
    initial_complete=set(initial)==set(initial_expected)
    positive=[row for row in census if K.from_json(row['segment_parameter'])>0]
    if not positive:
        assert probe_fraction is not None and initial_complete
        lifting=[]
        for nq in NORMS:
            vectors,_=block_vectors(nq)
            assert all(en==0 or en==nq for en,v in vectors)
            lifting.append({'norm':nq.to_json(),'bounded_dominant_F4_vectors':[
                {'norm':en.to_json(),'vector':enc(v)}for en,v in vectors]})
        out={'status':'COMPLETE_PROBE_SEGMENT_NO_NEW_MIRROR','edge_index':edge_index,
             'probe_fraction':probe_fraction,'proposal_dictionary_index':proposal_index,
             'base_source':base_path,'base_source_case':None,'proposal_source':proposal_path,
             'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in(base_path,proposal_path,STRUCTURE)},
             'gram':menc(g),'edge_indices':edge_ids,'edge_roots':menc(edge),'edge_gram':menc(eg),
             'start_corner_wall':start_wall,'start_corner':enc(k0),'proposed_endpoint':enc(k1),
             'start_height_cap':h0.to_json(),'end_height_cap':h1.to_json(),'height_cap':H.to_json(),
             'released_wall_residual':enc(residual_axis),'released_wall_residual_norm':B.to_json(),
             'endpoint_cosh_squared':C.to_json(),'norm_bounds':bounds,'profiles':profiles,
             'span_folded_roots':menc(list(span.values())),'folded_segment_roots':census,
             'expanded_segment_mirror_count':len(expanded),'edge_stabilizer_mirrors':len(edge_mirrors),
             'complete_start_stabilizer':True,'start_stabilizer_mirrors':len(initial),
             'H9_lifting_block_census':lifting,'H9_start_stabilizer_mirrors':len(initial)+24,
             'stats':dict(stats),'seconds':round(time.monotonic()-start,3),
             'scope':'Complete all-norm census on the explicitly stated strict initial portion of a proposed edge segment. No new mirror meets that portion. The unsearched remaining edge and other residual norm orbits remain open. Independent F4 energy exhaustion gives the same no-crossing conclusion in full H9.'}
        output.write_text(json.dumps(out,indent=2)+'\n')
        output.with_name(output.stem+'_partial.json').write_text(json.dumps({'status':'SUPERSEDED_BY_COMPLETE_CERTIFICATE','complete_certificate':str(output.relative_to(ROOT))},indent=2)+'\n')
        print({'status':out['status'],'quadratics':stats['quadratics'],'seconds':out['seconds']},flush=True)
        return out
    nearest_parameter=min(K.from_json(row['segment_parameter'])for row in positive)
    nearest_rows=[row for row in positive if K.from_json(row['segment_parameter'])==nearest_parameter]
    oriented=[]
    for row in nearest_rows:
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm']);a=pair(k0,g,r)
        if a>0 and all(pair(r,g,e)==0 for e in edge):
            r=tuple(-z for z in r);a=-a
        if a<0:oriented.append((a*a/q,q,r,row))
    assert oriented
    _,q,nearest,nearest_row=min(oriented,key=lambda row:(row[0],row[1]))
    point=tuple((1-nearest_parameter)*a+nearest_parameter*b for a,b in zip(k0,k1))
    new_incident=edge+[nearest];ng=tuple(tuple(pair(r,g,s)for s in new_incident)for r in new_incident)
    nearest_group=finite_orbit(new_incident,new_incident,g)
    nearest_full={key:r for key,r in expanded.items()if pair(point,g,r)==0};nearest_full.update(edge_mirrors)
    assert set(nearest_group)==set(nearest_full)
    proposed_root=tuple(map(K.from_json,proposed['vector']))
    proposed_root=(*proposed_root[:5],proposed_root[9])
    lifting=[]
    for nq in NORMS:
        vectors,_=block_vectors(nq)
        assert all(en==0 or en==nq for en,v in vectors)
        assert all(not any(v)for en,v in vectors if all((nq-en).sign(e)>0 for e in range(4)))
        lifting.append({'norm':nq.to_json(),'bounded_dominant_F4_vectors':[
            {'norm':en.to_json(),'vector':enc(v)}for en,v in vectors]})
    g9=Host(2).g
    lift=lambda r:(r[0],*r[1:5],*(Z,)*4,r[5])
    other=[tuple(-O if i==j else Z for i in range(10))for j in range(5,9)]
    f4mirrors=finite_orbit(other,other,g9);assert len(f4mirrors)==24
    for generators,expected in((edge,edge_mirrors),(edge+[walls[start_wall]],initial),(new_incident,nearest_full)):
        full=finite_orbit([*map(lift,generators),*other],[*map(lift,generators),*other],g9)
        assert set(full)=={line_key(lift(r))for r in expected.values()}|set(f4mirrors)
    out={'status':'COMPLETE_SHIFTED_PAIRING_EDGE_CENSUS','edge_index':edge_index,'proposal_dictionary_index':proposal_index,
         'base_source':base_path,'base_source_case':None,'proposal_source':proposal_path,
         'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in(base_path,proposal_path,STRUCTURE)},
         'gram':menc(g),'edge_indices':edge_ids,'edge_roots':menc(edge),'edge_gram':menc(eg),
         'start_corner_wall':start_wall,'start_corner':enc(k0),'proposed_endpoint':enc(k1),
         'start_height_cap':h0.to_json(),'end_height_cap':h1.to_json(),'height_cap':H.to_json(),
         'released_wall_residual':enc(residual_axis),'released_wall_residual_norm':B.to_json(),
         'endpoint_cosh_squared':C.to_json(),'norm_bounds':bounds,'profiles':profiles,
         'span_folded_roots':menc(list(span.values())),'folded_segment_roots':census,
         'expanded_segment_mirror_count':len(expanded),'edge_stabilizer_mirrors':len(edge_mirrors),
         'complete_start_stabilizer':initial_complete,'start_stabilizer_mirrors':len(initial),
         'expected_start_stabilizer_mirrors':len(initial_expected),
         'nearest_parameter':nearest_parameter.to_json(),'nearest_root':enc(nearest),
         'nearest_norm':q.to_json(),'nearest_point':enc(point),'nearest_point_norm':pair(point,g,point).to_json(),
         'nearest_incident_gram':menc(ng),'nearest_full_stabilizer_mirrors':len(nearest_full),
         'nearest_is_proposed_mirror':line_key(nearest)==line_key(proposed_root),
         'nearest_old_wall_pairings':enc(tuple(pair(nearest,g,r)for r in walls)),
         'H9_lifting_block_census':lifting,'other_F4_mirrors':24,
         'H9_edge_stabilizer_mirrors':len(edge_mirrors)+24,
         'H9_start_stabilizer_mirrors':len(initial)+24,
         'H9_end_stabilizer_mirrors':len(nearest_full)+24,
         'stats':dict(stats),'seconds':round(time.monotonic()-start,3),
         'scope':'Complete census of mirrors meeting the closed specified H5 edge segment, modulo its finite edge reflection group, using all necessary primitive norms, all bounded integral edge pairings, all embeddings, and every shifted integral pairing with the released starting wall in a rigorous finite alphabet. Original time bounds are diagnostic only. No unit-orbit restriction or residual-square restriction. The additional bounded F4 census excludes mixed-block crossing mirrors and lifts the conclusion to H9. Full starting and nearest-corner stabilizers are checked by finite root-orbit comparisons. One edge step, not finite volume or completion of H5/H9.'}
    output.write_text(json.dumps(out,indent=2)+'\n')
    output.with_name(output.stem+'_partial.json').write_text(json.dumps({
        'status':'SUPERSEDED_BY_COMPLETE_CERTIFICATE','completed_norms':len(bounds),
        'complete_certificate':str(output.relative_to(ROOT))},indent=2)+'\n')
    print({k:out[k]for k in('status','expanded_segment_mirror_count','edge_stabilizer_mirrors','complete_start_stabilizer','nearest_parameter','nearest_is_proposed_mirror','stats','seconds')},flush=True)
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--edge',type=int,choices=(0,1,2,3),default=2)
    parser.add_argument('--probe-fraction')
    args=parser.parse_args();run(args.edge,args.probe_fraction)
