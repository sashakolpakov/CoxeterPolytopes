#!/usr/bin/env python3
"""Complete bounded-segment root census, up to the finite edge Weyl group.

All integral pairing profiles are included, rather than just finite Coxeter
labels. All norm classes and binary residual square classes are included.
"""
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

BASE='results/fresh_h9_quartic_mixed_20261013_endpoint.json'
PROPOSAL='results/fresh_h9_quartic_edge_walk_20261013_profiles.json'
OUT=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013.json'
STRUCTURE='results/fresh_h9_glued_host_20261013_structure.json'
STAGES={1:(BASE,None,PROPOSAL,[1,2,4,12],6,''),
        2:(PROPOSAL,0,'results/fresh_h9_quartic_edge_walk_20261013_next.json',[1,2,4,70],12,'_second'),
        3:('results/fresh_h9_quartic_edge_walk_20261013_next.json',0,
           'results/fresh_h9_quartic_edge_walk_20261013_f4.json',[1,2,70,72],4,'_third')}


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


def run(stage=1):
    start=time.monotonic();g=Host(1).g
    base_path,base_case,proposal_path,edge_ids,start_wall,suffix=STAGES[stage]
    output=OUT.with_name(OUT.stem+suffix+'.json')
    base=json.loads((ROOT/base_path).read_text());proposal=json.loads((ROOT/proposal_path).read_text())
    if base_case is not None:base=base['cases'][base_case]
    walls=[tuple(map(K.from_json,r))for r in base['forced_H5_face']['all_walls']]
    inside=tuple(map(K.from_json,base['forced_H5_face']['strict_interior']))
    edge=[walls[i]for i in edge_ids]
    initial_vertex=finite_vertex(walls,g,inside,edge_ids+[start_wall])
    final_vertex=proposal['cases'][0]['forced_H5_face']['finite_vertex']
    k0=tuple(map(K.from_json,initial_vertex['vector']));k0=tuple(z/k0[0]for z in k0)
    k1=tuple(map(K.from_json,final_vertex['vector']));k1=tuple(z/k1[0]for z in k1)
    assert k0[0]==k1[0]==1 and all(pair(k,g,e)==0 for k in(k0,k1)for e in edge)
    h0=ALPHA/(-pair(k0,g,k0))-1;h1=ALPHA/(-pair(k1,g,k1))-1;H=max(h0,h1)
    assert pair(k1,g,walls[start_wall])<0
    eg=tuple(tuple(pair(r,g,s)for s in edge)for r in edge);ei=inverse(eg)
    basis,_,_,_,_=complement(edge,g)
    first,second=basis
    if not first[0]:first,second=second,first
    axis=tuple(z/first[0]for z in first)
    direction=tuple(z-second[0]*y for z,y in zip(second,axis))
    assert axis[0]==1 and direction[0]==0
    qa=pair(direction,g,direction);assert all(qa.sign(e)>0 for e in range(4))
    census=[];profiles=[];folded={};span={};stats=Counter();bounds=[]
    for q in NORMS:
        alph=[alphabet((q*pair(e,g,e),)*4)for e in edge]
        caps=(H*q/ALPHA,*(-q/ALPHA for _ in range(3)))
        positive=balanced_alphabet(caps)
        times=sorted({x for z in positive for x in (z,-z)})
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
            profiles.append(record);qrows.append(profile_index)
            for x in times:
                anchor=tuple(z+(x-r0[0])*a for z,a in zip(r0,axis))
                qb=2*pair(anchor,g,direction);qc=pair(anchor,g,anchor)-q
                discriminant=qb*qb-4*qa*qc;stats['quadratics']+=1
                radical=square_root(discriminant)
                if radical is None:continue
                stats['square_discriminants']+=1
                for t in {(-qb+radical)/(2*qa),(-qb-radical)/(2*qa)}:
                    r=tuple(a+t*d for a,d in zip(anchor,direction))
                    assert r[0]==x and pair(r,g,r)==q
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
        print({'norm':q.to_json(),'profiles':len(qrows),'time_coordinates':len(times),
               'folded_segment_mirrors':len(census),'seconds':round(time.monotonic()-start,2)},flush=True)
        record={'norm':q.to_json(),'pairing_alphabets':[enc(a)for a in alph],
                'time_caps':enc(caps),'time_alphabet':enc(times),'profile_indices':qrows}
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
    proposed_root=tuple(map(K.from_json,proposal['compatible_dictionary'][proposal['cases'][0]['dictionary_index']]['vector']))
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
    out={'status':'COMPLETE_EDGE_SEGMENT_ROOT_CENSUS','stage':stage,
         'base_source':base_path,'base_source_case':base_case,'proposal_source':proposal_path,
         'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in(base_path,proposal_path,STRUCTURE)},
         'gram':menc(g),'edge_indices':edge_ids,'edge_roots':menc(edge),'edge_gram':menc(eg),
         'start_corner_wall':start_wall,'start_corner':enc(k0),'proposed_endpoint':enc(k1),
         'start_height_cap':h0.to_json(),'end_height_cap':h1.to_json(),'height_cap':H.to_json(),
         'time_one_plane_vector':enc(axis),'time_zero_plane_vector':enc(direction),
         'time_zero_norm':qa.to_json(),'norm_bounds':bounds,'profiles':profiles,
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
         'scope':'Complete census of mirrors meeting the closed specified H5 edge segment, modulo its finite edge reflection group, using all necessary primitive norms, all bounded integral pairings, all embeddings, and every time coordinate in a rigorous finite alphabet. No unit-orbit restriction or residual-square restriction. The additional bounded F4 census excludes mixed-block crossing mirrors and lifts the conclusion to H9. Full starting and nearest-corner stabilizers are checked by finite root-orbit comparisons. One edge step, not finite volume or completion of H5/H9.'}
    output.write_text(json.dumps(out,indent=2)+'\n')
    output.with_name(output.stem+'_partial.json').write_text(json.dumps({
        'status':'SUPERSEDED_BY_COMPLETE_CERTIFICATE','completed_norms':len(bounds),
        'complete_certificate':str(output.relative_to(ROOT))},indent=2)+'\n')
    print({k:out[k]for k in('status','expanded_segment_mirror_count','edge_stabilizer_mirrors','complete_start_stabilizer','nearest_parameter','nearest_is_proposed_mirror','stats','seconds')},flush=True)
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',type=int,choices=(1,2,3),default=1)
    args=parser.parse_args();run(args.stage)
