#!/usr/bin/env python3
"""A field-independent contact-sum obstruction on the full spacelike face.

All1024 active subsets of the ten nontrivial walls are classified. Every
face interior of dimension2..8 and30 of32 ray interiors are excluded by
one exact Cauchy certificate. Timelike-generator directions stay open.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
import fresh_raw23_four_face_contact_preflight_20261018 as F

ROOT=F.ROOT;P=F.P;H=F.H;Q=F.Q
TAG='fresh_raw23_spacelike_contact_sum_20261018'

def one(old):
    rays=old['initial_cone']['rays'];by={r['index']:r for r in rays}
    xyz={i:F.vec(r['coordinates'])for i,r in by.items()};zeros={i:set(r['zero_facets'])for i,r in by.items()}
    rows=H.mat(old['initial_cone']['normal_rows']);metric=H.mat(old['metric']);sc=F.vec(old['original_squared_scales'])
    contacts={i:[-p for p in F.vec(r['products'])]for i,r in by.items()}
    assert all(p.sign()>=0 for ps in contacts.values()for p in ps)
    space=tuple(i for i,r in by.items()if r['norm_sign']>0);common=tuple(sorted(set.intersection(*(zeros[i]for i in space))))
    assert len(space)==32 and common==(10,11)
    assert tuple(i for i in by if set(common)<=zeros[i])==space
    remaining=tuple(i for i in range(12)if i not in common)
    lower=[2 if(s-Q(4)).sign()>=0 else 1 for s in sc]
    assert all((s-Q(l*l)).sign()>=0 for s,l in zip(sc,lower))
    cache={}
    def basis(ids):
        if ids not in cache:
            bs=[]
            for i in ids:
                if F.rank([xyz[j]for j in bs+[i]])>len(bs):bs.append(i)
            cache[ids]=bs
        return cache[ids]
    fullbasis=basis(space);assert len(fullbasis)==8
    fullgram=[[P.inner(xyz[i],metric,xyz[j])for j in fullbasis]for i in fullbasis];assert H.positive(fullgram)
    decisions=[];faces={}
    for size in range(11):
        for subset in combinations(remaining,size):
            active=common+subset;full=tuple(i for i in space if set(active)<=zeros[i]);bs=basis(full)
            arank=F.rank([rows[i]for i in active])
            rec=dict(additional_active_walls=subset,active_rank=arank,all_extreme_rays=full,dimension=len(bs));decisions.append(rec)
            if not full:
                rec['origin_only']=True;continue
            closed=tuple(sorted(set.intersection(*(zeros[i]for i in full))))
            assert F.rank([rows[i]for i in closed])+len(bs)==10
            rec.update(origin_only=False,closed_active_walls=closed)
            if closed in faces:
                assert faces[closed]['all_extreme_rays']==full
                faces[closed]['generating_subsets'].append(subset);continue
            gram=[[P.inner(xyz[i],metric,xyz[j])for j in bs]for i in bs];assert H.positive(gram)
            nonzero=[j for j in range(12)if any(contacts[i][j]!=0 for i in full)]
            assert nonzero==sorted(set(range(12))-set(closed))
            p=[sum((contacts[i][j]for j in nonzero),Q())for i in bs];inv=H.inverse(gram)
            z=[sum((inv[a][b]*p[b]for b in range(len(bs))),Q())for a in range(len(bs))]
            q=sum((a*b for a,b in zip(p,z)),Q());L=sum(lower[j]for j in nonzero);margin=Q(L*L)-2*q
            assert q.sign()>0
            faces[closed]=dict(active_walls=closed,generating_subsets=[subset],dimension=len(bs),
                all_extreme_rays=full,basis_rays=bs,plane_gram=gram,
                nonzero_wall_forms=nonzero,individual_contact_lower_bounds=[lower[j]for j in nonzero],
                summed_contact_on_basis=p,projection_coefficients=z,summed_form_squared_norm=q,
                required_contact_sum=L,Cauchy_margin=margin,relative_interior_excluded=margin.sign()>0)
    ordered=[faces[k]for k in sorted(faces)];dims=Counter(f['dimension']for f in ordered)
    assert len(decisions)==1024 and dict(dims)=={1:32,2:112,3:184,4:180,5:112,6:44,7:10,8:1}
    assert all(f['relative_interior_excluded']for f in ordered if f['dimension']>=2)
    survivors=[f for f in ordered if not f['relative_interior_excluded']]
    assert len(survivors)==2 and all(f['dimension']==1 for f in survivors)
    assert sorted(f['all_extreme_rays']for f in survivors)==[(6,),(31,)]
    print('CONTACT_SUM',old['kinds'],'faces',len(ordered),'excluded',sum(f['relative_interior_excluded']for f in ordered),'surviving_rays',[6,31],flush=True)
    return dict(kinds=old['kinds'],spacelike_rays=space,common_active_walls=common,
        full_face_basis=fullbasis,full_face_gram=fullgram,contact_lower_bounds=lower,
        active_decisions=decisions,faces=ordered,counts=dict(active_subsets=1024,nonzero_faces=len(ordered),
            excluded_faces=sum(f['relative_interior_excluded']for f in ordered),excluded_ray_interiors=30,
            excluded_higher_face_interiors=643,unexcluded_ray_interiors=2),
        remaining_projected_ray_directions=[6,31])

def run():
    old=H.read(P.TAG);four=H.read(F.TAG);H.bind(old);H.bind(four)
    cases=[one(o)for o in old['cases']]
    paths={f'scripts/{TAG}.py',f'scripts/{F.TAG}.py',f'results/{F.TAG}.json',f'results/{F.TAG}.md',
        f'scripts/{P.TAG}.py',f'results/{P.TAG}.json',f'results/{P.TAG}.md'}
    return P.enc(dict(status='PASS_RAW23_SPACELIKE_FACE_CONTACT_SUM',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},cases=cases,
        active_subsets_checked=3072,classified_nonzero_faces=2025,
        excluded_higher_face_interiors=1929,excluded_ray_interiors=90,
        remaining_ray_directions=6,full_spacelike_face_dimension=8,
        extension_ground_field_restricted=False,arbitrary_passive_components_allowed=True,
        all_nonzero_spacelike_ray_combinations_reduced_to_six_ray_directions=True,
        zero_projection_excluded=False,
        directions_with_timelike_generators_excluded=False,arbitrary_completion_excluded=False,
        capacity_or_model_changed=False,new_dimension_bound=None,
        scope='Complete nonzero face lattice of each original maximal spacelike8face from all1024 active subsets, with exact summed-contact Cauchy certificates. Nonzero relative interiors except two ray directions per cone are excluded over arbitrary extension fields and passive lifts. Zero projection, remaining ray directions and all directions using timelike rays stay open; no arbitrary completion, pair capacity, model or dimension theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();out=run()
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],flush=True)
