#!/usr/bin/env python3
"""Saved higher-face inventories, independent arithmetic/rank implementation.

Uses the frozen standalone four-face checker's Q(sqrt5) arithmetic, reverse
column ranks and congruence inertia, without importing either producer.
"""
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import check_raw23_four_face_saved_20261018 as H

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_higher_face_signature_inventory_20261018'
OLD='fresh_raw23_union_polar_completion_preflight_20261018'

def run():
    read=lambda tag:json.loads((ROOT/f'results/{tag}.json').read_text())
    data=read(TAG);old=read(OLD)
    for d in(data,old):
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert data['status']=='PASS_RAW23_HIGHER_FACE_SIGNATURE_INVENTORY'
    assert data['dimensions']==list(range(5,10))and len(data['cases'])==len(old['cases'])==3
    for case,o in zip(data['cases'],old['cases'],strict=True):
        assert case['kinds']==o['kinds'];rays=o['initial_cone']['rays'];rows=H.mat(o['initial_cone']['normal_rows'])
        xyz={r['index']:H.vec(r['coordinates'])for r in rays};zeros={r['index']:set(r['zero_facets'])for r in rays}
        signs={r['index']:r['norm_sign']for r in rays};metric=H.mat(o['metric']);cache={}
        def rank(ids):
            if ids not in cache:cache[ids]=H.rank([xyz[i]for i in ids])
            return cache[ids]
        assert len(case['inventories'])==5
        for dim,inv in zip(range(5,10),case['inventories'],strict=True):
            assert inv['dimension']==dim;codim=10-dim;expected=list(combinations(range(12),codim))
            assert len(inv['decisions'])==len(expected);found={};fulls={}
            for a,rr in zip(expected,inv['decisions'],strict=True):
                full=tuple(i for i in xyz if set(a)<=zeros[i]);arank=H.rank([rows[i]for i in a]);r=rank(full)
                isface=arank==codim and r==dim
                assert rr['active']==list(a)and rr['active_rank']==arank
                assert rr['all_extreme_rays']==list(full)and rr['ray_span_rank']==r and rr['target_face']==isface
                if not isface:continue
                common=tuple(sorted(set.intersection(*(zeros[i]for i in full))))
                assert H.rank([rows[i]for i in common])==codim and rr['closed_active_walls']==list(common)
                found.setdefault(common,[]).append(list(a));fulls[common]=full
            assert [f['active_walls']for f in inv['faces']]==[list(k)for k in sorted(found)]
            counts=Counter();space=positive=0
            for f in inv['faces']:
                k=tuple(f['active_walls']);full=fulls[k];bs=tuple(f['basis_rays'])
                assert f['all_extreme_rays']==list(full)and f['active_bases']==found[k]
                assert len(bs)==len(set(bs))==dim and set(bs)<=set(full)and rank(bs)==dim
                gram=[[H.bil(xyz[i],metric,xyz[j])for j in bs]for i in bs];sig=H.inertia(gram)
                assert H.mat(f['plane_gram'])==gram and f['plane_inertia']==sig
                sp=all(signs[i]>0 for i in full);pd=sig==[dim,0,0]
                assert f['all_rays_spacelike']==sp and f['positive_definite_plane']==pd
                counts[sp,tuple(sig),len(full)]+=1;space+=sp;positive+=pd
            assert inv['counts']==dict(active_subsets=len(expected),distinct_faces=len(found),entirely_spacelike=space,positive_definite=positive)
            classes=[dict(all_rays_spacelike=s,inertia=list(i),rays=n,count=c)for(s,i,n),c in sorted(counts.items())]
            assert inv['signature_classes']==classes
            print('CHECKED_HIGHER_FACES',case['kinds'],dim,inv['counts'],flush=True)
    assert data['all_original_faces_in_these_dimensions_classified']
    assert all(data[k]is False for k in('new_normal_direction_exclusion','arbitrary_completion_exclusion','model_changed'))
    assert data['new_dimension_bound']is None
    print('PASS_RAW23_HIGHER_FACE_SAVED',flush=True)

if __name__=='__main__':run()
