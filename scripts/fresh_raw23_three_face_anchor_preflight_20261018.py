#!/usr/bin/env python3
"""Classify original three-dimensional faces and finite-anchor coverage.

No native125-choice candidate enumeration. Every face found from three
independent spacelike rays is expanded to its entire frozen ray set.
Every PD face is covered by all independent three-ray positive cones.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
import fresh_raw23_union_polar_completion_preflight_20261018 as P
Q=P.Q;H=P.H;ROOT=P.ROOT
TAG='fresh_raw23_three_face_anchor_preflight_20261018'
OLD=P.TAG
PINS={f'scripts/{OLD}.py':'801198f128ee0e8d3ff9587fbd65d656edaedc0d413e4b4665c094c290010ee8',
 f'results/{OLD}.json':'eee8c3836a16c6342b2882c3b65d84999bb7fbc7363413fbc59626c739ff8361',
 f'results/{OLD}.md':'42bd1f4baeecd9c8188a1590adf1e689c5f9582ae5df378a75d94dc5116da973',
 'scripts/check_raw23_same_family_saved_20261018.py':'5b1ab9873b2e069f18966a7052856d8a9df8a8c1b0eef4794522a812a54616a7'}

def rank(rows):return len(P.rref(rows)[1])if rows else 0
def getq(x):return[Q(*a)for a in x]
def cone_bound(h,p,s):
    stationary=[]
    for size in(1,2,3):
        for ids in combinations(range(3),size):
            hh=H.principal(h,ids);j=H.inverse(hh);pp=[p[k]for k in ids]
            z=[sum((j[a][b]*pp[b]for b in range(size)),Q())for a in range(size)]
            admissible=all(q.sign()>0 for q in z)
            val=2*sum((a*b for a,b in zip(pp,z)),Q())/s if admissible else None
            stationary.append(dict(support=ids,coefficients=z,strictly_positive=admissible,bound=val))
    values=[r['bound']for r in stationary if r['strictly_positive']];assert values
    maximum=values[0]
    for v in values[1:]:
        if(v-maximum).sign()>0:maximum=v
    maximizers=[r['support']for r in stationary if r['strictly_positive']and r['bound']==maximum]
    interior=(0,1,2)in maximizers
    return dict(supports=stationary,maximum_square=maximum,maximizing_supports=maximizers,
        attained_in_relative_interior=interior,
        strictly_below_four_in_relative_interior=(maximum-4).sign()<0 or maximum==4 and not interior)

def full_cone_bounds(full,xyz,metric,byid,sc):
    """Conic Caratheodory cover; never infer full-face strictness from a subcone boundary."""
    triples=[]
    for ids in combinations(full,3):
        if rank([xyz[i]for i in ids])!=3:continue
        h=[[P.inner(xyz[i],metric,xyz[j])for j in ids]for i in ids]
        assert H.positive(h)
        triples.append((ids,h))
    assert triples
    bounds=[]
    for wall,s in enumerate(sc):
        p=[-Q(*byid[r]['products'][wall])for r in full]
        assert all(x.sign()>=0 for x in p)
        if all(x==0 for x in p):continue
        subs=[]
        for ids,h in triples:
            pp=[-Q(*byid[r]['products'][wall])for r in ids]
            if all(x==0 for x in pp):
                subs.append(dict(rays=ids,zero_contact=True,maximum_square=Q()))
            else:subs.append(dict(rays=ids,zero_contact=False,**cone_bound(h,pp,s)))
        maximum=subs[0]['maximum_square']
        for sub in subs[1:]:
            if(sub['maximum_square']-maximum).sign()>0:maximum=sub['maximum_square']
        bounds.append(dict(wall=wall,contact_on_all_rays=p,subcones=subs,
            maximum_square=maximum,strictly_below_one=(maximum-1).sign()<0))
    excluded=[r['wall']for r in bounds if r['strictly_below_one']]
    return dict(independent_triple_cover=[dict(rays=ids,gram=h)for ids,h in triples],
        full_cone_contact_bounds=bounds,native_excluding_walls=excluded,
        relative_interior_native_excluded=bool(excluded))

def run():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=H.read(OLD);H.bind(old);cases=[];summary=Counter()
    for ci,o in enumerate(old['cases']):
        c=o['initial_cone'];rays=c['rays'];byid={r['index']:r for r in rays}
        assert len(byid)==len(rays)==50
        rows=H.mat(c['normal_rows']);metric=H.mat(o['metric']);sc=getq(o['original_squared_scales'])
        space=[r for r in rays if r['norm_sign']>0];assert len(space)==32
        xyz={r['index']:getq(r['coordinates'])for r in rays}
        decisions=[];faces={}
        for triple in combinations(space,3):
            ids=tuple(r['index']for r in triple)
            common=sorted(set.intersection(*(set(r['zero_facets'])for r in triple)))
            nr=rank([rows[i]for i in common]);independent=rank([xyz[i]for i in ids])==3
            found=nr==7 and independent
            decisions.append(dict(rays=ids,common_active_walls=common,common_rank=nr,
                linearly_independent=independent,three_dimensional_face=found))
            if not found:continue
            full=tuple(r['index']for r in rays if set(common)<=set(r['zero_facets']))
            actual_common=sorted(set.intersection(*(set(byid[i]['zero_facets'])for i in full)))
            assert actual_common==common and rank([xyz[i]for i in full])==3
            key=tuple(actual_common)
            if key in faces:
                assert faces[key]['all_extreme_rays']==full
                faces[key]['spacelike_generating_triples'].append(ids);continue
            h=[[P.inner(xyz[i],metric,xyz[j])for j in ids]for i in ids]
            inertia=H.inertia(h);pd=inertia==[3,0,0]
            allspace=all(byid[i]['norm_sign']>0 for i in full);simplicial=len(full)==3
            record=dict(active_walls=actual_common,dimension=3,all_extreme_rays=full,
                extreme_ray_signs=[byid[i]['norm_sign']for i in full],
                spacelike_generating_triples=[ids],basis_rays=ids,plane_gram=h,
                plane_inertia=inertia,all_rays_spacelike=allspace,simplicial=simplicial,
                positive_definite_plane=pd,anchor_status=None,anchor_bounds=[],anchors=None)
            if not allspace:record['anchor_status']='contains_nonspacelike_ray'
            elif not pd:record['anchor_status']='non_positive_definite_span'
            else:
                record.update(full_cone_bounds(full,xyz,metric,byid,sc))
                if not simplicial:
                    record['anchor_status']='nonsimplicial_strict_native_exclusion'if record['relative_interior_native_excluded']else'nonsimplicial_no_strict_native_exclusion'
                    faces[key]=record;continue
                contacts=[[-Q(*byid[r]['products'][i])for r in ids]for i in range(len(rows))]
                assert all(q.sign()>=0 for p in contacts for q in p)
                record['contact_magnitudes']=contacts
                eligible=[]
                for i,p in enumerate(contacts):
                    if all(q==0 for q in p):continue
                    rec=dict(wall=i,contact=p,**cone_bound(h,p,sc[i]));record['anchor_bounds'].append(rec)
                    if rec['strictly_below_four_in_relative_interior']:eligible.append(i)
                record['eligible_anchor_walls']=eligible
                erank=rank([contacts[i]for i in eligible]);record['eligible_anchor_rank']=erank
                if erank==3:
                    record['anchors']=next(ids for ids in combinations(eligible,3)if rank([contacts[i]for i in ids])==3)
                    record['anchor_status']='three_independent_finite_anchors'
                else:record['anchor_status']='insufficient_finite_anchor_rank'
            faces[key]=record
        records=[faces[k]for k in sorted(faces)]
        counts=Counter(r['anchor_status']for r in records)
        counts['triple_domain']=len(decisions)
        counts['found_triples']=sum(d['three_dimensional_face']for d in decisions)
        counts['distinct_faces']=len(records)
        counts['all_spacelike_faces']=sum(r['all_rays_spacelike']for r in records)
        counts['simplicial_faces']=sum(r['simplicial']for r in records)
        counts['PD_faces']=sum(r['positive_definite_plane']for r in records)
        counts['native_excluded_faces']=sum(r.get('relative_interior_native_excluded',False)for r in records)
        counts['independent_subcone_count']=sum(len(r.get('independent_triple_cover',[]))for r in records)
        assert len(decisions)==4960
        summary.update(counts)
        cases.append(dict(kinds=o['kinds'],triple_decisions=decisions,faces=records,counts=dict(counts)))
        print('CASE',o['kinds'],dict(counts),flush=True)
    paths=set(PINS)|{f'scripts/{TAG}.py'}
    return P.enc(dict(status='PASS_RAW23_THREE_FACE_ANCHOR_PREFLIGHT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        cases=cases,summary=dict(summary),native_candidate_enumeration_performed=False,
        source_domain='All independent triples of original spacelike rays; every found face expanded to all original extreme rays.',
        all_entirely_spacelike_three_faces_covered=True,
        all_three_faces_regardless_of_ray_signs_claimed=False,
        nonsimplicial_faces_treated_as_unresolved=False,all_PD_faces_use_complete_triple_cover=True,
        strict_less_than_one_only_for_full_cone_exclusion=True,
        arbitrary_passive_lifts_relaxed=True,unknown_passive_contacts_assumed_zero=False,
        geometric_realization_claimed=False,new_normal_claimed=False,pair_exclusion_claimed=False,
        exported_columns=[],optimizer_invocations=0,
        scope='Classification and exact native-contact exclusion for original three-cone3D faces reached by independent spacelike triples. Every entirely spacelike3D face has such a triple and is covered. Complete independent-triple covers handle nonsimplicial PD faces; only a strict full-cone contact bound below1 excludes its relative interior. Non-PD and faces not reached by spacelike triples remain outside the exclusion; no125-choice native enumeration or completion/capacity claim.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--replace-draft',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or args.replace_draft or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
