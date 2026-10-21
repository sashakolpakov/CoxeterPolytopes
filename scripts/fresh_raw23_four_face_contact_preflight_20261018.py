#!/usr/bin/env python3
"""Complete original4-face inventory and strict contact bounds on PD spans.

All six-wall active sets are tested, including faces with timelike rays.
Bounds cover all real coefficients and arbitrary positive passive lifts;
no native625 enumeration, arithmetic-field restriction on candidates,
capacity, completion, model change or optimization is asserted.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
import fresh_raw23_union_polar_completion_preflight_20261018 as P
Q=P.Q;H=P.H;ROOT=P.ROOT
TAG='fresh_raw23_four_face_contact_preflight_20261018';OLD=P.TAG
THREE='fresh_raw23_three_face_anchor_preflight_20261018'
PINS={f'scripts/{OLD}.py':'801198f128ee0e8d3ff9587fbd65d656edaedc0d413e4b4665c094c290010ee8',
 f'results/{OLD}.json':'eee8c3836a16c6342b2882c3b65d84999bb7fbc7363413fbc59626c739ff8361',
 f'results/{OLD}.md':'42bd1f4baeecd9c8188a1590adf1e689c5f9582ae5df378a75d94dc5116da973',
 f'scripts/{THREE}.py':'30ff587a47bf671dc0ce47ce4c5c4924c6eb8eeab713efaf569b08dfb046d0d6',
 f'results/{THREE}.json':'748f41703fc8b68347544e53bd2efcbc804b001aa7cc61f275c08badbde382e3',
 f'results/{THREE}.md':'2262069090075514fd06011e8ae32e0755a51334699ed8b7947e46ca31ec814a',
 'scripts/check_raw23_same_family_saved_20261018.py':'5b1ab9873b2e069f18966a7052856d8a9df8a8c1b0eef4794522a812a54616a7'}
def rank(rows):return len(P.rref(rows)[1])if rows else 0
def vec(xs):return[Q(*x)for x in xs]

def one_case(old):
    cone=old['initial_cone'];rays=cone['rays'];byid={r['index']:r for r in rays}
    assert len(rays)==len(byid)==50
    rows=H.mat(cone['normal_rows']);metric=H.mat(old['metric']);sc=vec(old['original_squared_scales'])
    assert len(rows)==len(sc)==12 and len(metric)==10
    xyz={i:vec(r['coordinates'])for i,r in byid.items()}
    contacts={i:[-q for q in vec(r['products'])]for i,r in byid.items()}
    assert all(q.sign()>=0 for p in contacts.values()for q in p)
    rcache={};gcache={}
    def rayrank(ids):
        ids=tuple(ids)
        if ids not in rcache:rcache[ids]=rank([xyz[i]for i in ids])
        return rcache[ids]
    def gram(ids):
        out=[]
        for i in ids:
            rr=[]
            for j in ids:
                key=tuple(sorted((i,j)))
                if key not in gcache:gcache[key]=P.inner(xyz[i],metric,xyz[j])
                rr.append(gcache[key])
            out.append(rr)
        return out
    decisions=[];faces={}
    for active in combinations(range(12),6):
        arank=rank([rows[i]for i in active])
        full=tuple(i for i in byid if set(active)<=set(byid[i]['zero_facets']))
        dimension=rayrank(full);found=arank==6 and dimension==4
        rec=dict(active=active,active_rank=arank,all_extreme_rays=full,ray_span_rank=dimension,four_face=found)
        decisions.append(rec)
        if not found:continue
        common=tuple(sorted(set.intersection(*(set(byid[i]['zero_facets'])for i in full))))
        assert rank([rows[i]for i in common])==6
        rec['closed_active_walls']=common
        if common in faces:
            assert faces[common]['all_extreme_rays']==full
            faces[common]['active_bases'].append(active);continue
        basis=next(ids for ids in combinations(full,4)if rayrank(ids)==4)
        h=gram(basis);inertia=H.inertia(h)
        signs=[byid[i]['norm_sign']for i in full]
        faces[common]=dict(active_walls=common,active_bases=[active],dimension=4,
            all_extreme_rays=full,ray_signs=signs,all_rays_spacelike=all(s>0 for s in signs),
            basis_rays=basis,plane_gram=h,plane_inertia=inertia,
            positive_definite_plane=inertia==[4,0,0],simplicial=len(full)==4)
    face_records=[faces[k]for k in sorted(faces)]
    preliminary=Counter((r['all_rays_spacelike'],tuple(r['plane_inertia']),len(r['all_extreme_rays']))for r in face_records)
    print('INVENTORY',old['kinds'],'active_sets',len(decisions),'faces',len(face_records),dict(preliminary),flush=True)
    # Independent <=4-ray supports cover every positive cone by conic
    # Caratheodory. Cache their equations once, sharing across all faces.
    supports={}
    for face in face_records:
        if not(face['all_rays_spacelike']and face['positive_definite_plane']):
            face.update(status='outside_PD_spacelike_bound',support_ids=[],contact_bounds=[],excluding_walls=[],relative_interior_excluded=False)
            continue
        idslist=[]
        for size in range(1,5):
            for ids in combinations(face['all_extreme_rays'],size):
                if rayrank(ids)!=size:continue
                idslist.append(ids)
                if ids in supports:continue
                h=gram(ids);assert H.positive(h);j=H.inverse(h);stationary=[]
                for wall,s in enumerate(sc):
                    p=[contacts[i][wall]for i in ids]
                    z=[sum((j[a][b]*p[b]for b in range(size)),Q())for a in range(size)]
                    positive=all(q.sign()>0 for q in z)
                    maximum=2*sum((x*y for x,y in zip(p,z)),Q())/s if positive else None
                    stationary.append(dict(wall=wall,contact=p,coefficients=z,strictly_positive=positive,candidate_square=maximum))
                supports[ids]=dict(rays=ids,gram=h,stationary=stationary)
        face['support_ids']=idslist;bounds=[];excluded=[]
        for wall,s in enumerate(sc):
            fullp=[contacts[i][wall]for i in face['all_extreme_rays']]
            if all(p==0 for p in fullp):continue
            candidates=[(ids,supports[ids]['stationary'][wall]['candidate_square'])for ids in idslist
                if supports[ids]['stationary'][wall]['strictly_positive']]
            assert candidates
            maximum=candidates[0][1]
            for ids,v in candidates[1:]:
                if(v-maximum).sign()>0:maximum=v
            witnesses=[ids for ids,v in candidates if v==maximum]
            below=(maximum-1).sign()<0
            if below:excluded.append(wall)
            bounds.append(dict(wall=wall,contact_on_all_rays=fullp,maximum_square=maximum,
                maximizing_supports=witnesses,strictly_below_one=below))
        face.update(contact_bounds=bounds,excluding_walls=excluded,relative_interior_excluded=bool(excluded),
            status='strict_contact_exclusion'if excluded else'PD_no_strict_contact_exclusion')
    counts=Counter(r['status']for r in face_records)
    counts.update(active_six_sets=len(decisions),found_active_bases=sum(r['four_face']for r in decisions),
        distinct_faces=len(face_records),all_spacelike_faces=sum(r['all_rays_spacelike']for r in face_records),
        PD_faces=sum(r['positive_definite_plane']for r in face_records),
        simplicial_faces=sum(r['simplicial']for r in face_records),
        excluded_faces=sum(r['relative_interior_excluded']for r in face_records),
        unique_stationary_supports=len(supports),support_occurrences=sum(len(r['support_ids'])for r in face_records))
    assert len(decisions)==924
    print('CASE',old['kinds'],dict(counts),flush=True)
    return dict(kinds=old['kinds'],active_decisions=decisions,faces=face_records,
        stationary_supports=[supports[k]for k in sorted(supports,key=lambda k:(len(k),k))],counts=dict(counts))

def run():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=H.read(OLD);H.bind(old);H.bind(H.read(THREE))
    cases=[one_case(o)for o in old['cases']];summary=Counter()
    for c in cases:summary.update(c['counts'])
    paths=set(PINS)|{f'scripts/{TAG}.py'}
    return P.enc(dict(status='PASS_RAW23_FOUR_FACE_CONTACT_PREFLIGHT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},cases=cases,summary=dict(summary),
        all_original_four_faces_classified=True,all_PD_entirely_spacelike_faces_bounded=True,
        arbitrary_real_candidate_coefficients=True,candidate_ground_field_restriction=False,
        arbitrary_passive_lifts_relaxed=True,unknown_passive_contacts_assumed_zero=False,
        native_tuple_enumeration_performed=False,geometric_realization_claimed=False,new_normal_claimed=False,
        pair_exclusion_claimed=False,exported_columns=[],optimizer_invocations=0,
        scope='All original4faces classified from complete six-active-wall domain and full ray membership. Entirely-spacelike PD faces have exact full-cone bounds from all independent ray supports of size at most4. Strict nonzero contact<1 excludes relative-interior Coxeter normals over any field, including arbitrary passive lifts. Other faces and arbitrary completions remain open; no capacity/export/model claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();out=run()
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
