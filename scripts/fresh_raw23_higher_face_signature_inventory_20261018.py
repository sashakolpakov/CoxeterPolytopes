#!/usr/bin/env python3
"""Exact complete inventories of original raw23 faces of dimensions5..9.

This classifies spans only. It does not exclude any new normal direction
or assume that an all-spacelike ray set has positive-definite span.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
import fresh_raw23_four_face_contact_preflight_20261018 as F

ROOT=F.ROOT;P=F.P;H=F.H;Q=F.Q
TAG='fresh_raw23_higher_face_signature_inventory_20261018'

def one(old):
    rows=H.mat(old['initial_cone']['normal_rows']);metric=H.mat(old['metric'])
    rays=old['initial_cone']['rays'];xyz={r['index']:F.vec(r['coordinates'])for r in rays}
    zeros={r['index']:set(r['zero_facets'])for r in rays}
    assert len(rows)==12 and len(metric)==10 and len(rays)==len(xyz)==50
    cache={}
    def basis(ids):
        if ids not in cache:
            bs=[]
            for i in ids:
                if F.rank([xyz[j]for j in bs+[i]])>len(bs):bs.append(i)
            cache[ids]=bs
        return cache[ids]
    inventories=[]
    for dim in range(5,10):
        codim=10-dim;decisions=[];faces={}
        for active in combinations(range(12),codim):
            arank=F.rank([rows[i]for i in active])
            full=tuple(i for i in xyz if set(active)<=zeros[i]);bs=basis(full)
            found=arank==codim and len(bs)==dim
            r=dict(active=active,active_rank=arank,all_extreme_rays=full,ray_span_rank=len(bs),target_face=found)
            decisions.append(r)
            if not found:continue
            common=tuple(sorted(set.intersection(*(zeros[i]for i in full))))
            assert F.rank([rows[i]for i in common])==codim;r['closed_active_walls']=common
            if common in faces:
                assert faces[common]['all_extreme_rays']==full
                faces[common]['active_bases'].append(active);continue
            gram=[[P.inner(xyz[i],metric,xyz[j])for j in bs]for i in bs]
            sig=H.inertia(gram);signs=[next(r['norm_sign']for r in rays if r['index']==i)for i in full]
            faces[common]=dict(active_walls=common,active_bases=[active],all_extreme_rays=full,
                basis_rays=bs,plane_gram=gram,plane_inertia=sig,
                all_rays_spacelike=all(s>0 for s in signs),positive_definite_plane=sig==[dim,0,0])
        ordered=[faces[k]for k in sorted(faces)]
        types=Counter((x['all_rays_spacelike'],tuple(x['plane_inertia']),len(x['all_extreme_rays']))for x in ordered)
        print('FACE_SIGNATURES',old['kinds'],dim,len(ordered),dict(types),flush=True)
        inventories.append(dict(dimension=dim,decisions=decisions,faces=ordered,
            counts=dict(active_subsets=len(decisions),distinct_faces=len(ordered),
                entirely_spacelike=sum(x['all_rays_spacelike']for x in ordered),
                positive_definite=sum(x['positive_definite_plane']for x in ordered)),
            signature_classes=[dict(all_rays_spacelike=s,inertia=i,rays=n,count=c)for(s,i,n),c in sorted(types.items())]))
    return dict(kinds=old['kinds'],inventories=inventories)

def run():
    old=H.read(P.TAG);four=H.read(F.TAG);H.bind(old);H.bind(four)
    cases=[one(o)for o in old['cases']]
    paths={f'scripts/{TAG}.py',f'scripts/{F.TAG}.py',f'results/{F.TAG}.json',f'results/{F.TAG}.md',
        f'scripts/{P.TAG}.py',f'results/{P.TAG}.json',f'results/{P.TAG}.md'}
    return P.enc(dict(status='PASS_RAW23_HIGHER_FACE_SIGNATURE_INVENTORY',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},cases=cases,
        dimensions=[5,6,7,8,9],all_original_faces_in_these_dimensions_classified=True,
        new_normal_direction_exclusion=False,arbitrary_completion_exclusion=False,
        model_changed=False,new_dimension_bound=None,
        scope='Exact complete active-wall face inventory and Gram signatures for dimensions5..9 of the original frozen raw23 cones. Frozen ray completeness remains input; no contact bound or new nonexistence theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();out=run()
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],flush=True)
