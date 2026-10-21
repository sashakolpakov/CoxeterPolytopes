#!/usr/bin/env python3
"""Exact extreme-ray polar truncation preflight, not arbitrary completion.

All rays of the three twelve-wall cones are enumerated. Candidate Lorentz
projections are restricted to spacelike extreme-ray directions; arbitrary
passive lifts are tested by five exact finite-label scales per direction.
Known-wall native gates and the conjugate known-span Schur condition are
necessary only. No optimizer, column, or chamber claim.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import check_raw23_same_family_saved_20261018 as H
ROOT=H.ROOT;Q=H.Q
TAG='fresh_raw23_union_polar_completion_preflight_20261018'
UNION='fresh_raw23_same_family_union_support_pilot_20261018'
REL='fresh_raw23_same_family_relative_pilot_20261018'
PINS={
 f'scripts/{UNION}.py':'3c1beadb17d87913699f3900dd6332465c4dbf9e0beae9845a31f526b09f3d30',
 f'results/{UNION}.json':'f7af708a8f54851c1d01904cc5c9c4aa41a889ed6b7f6a688493c92deb80d83b',
 f'results/{REL}.md':'bb69bd6a80a1f918f6e42bfad6c4d1f3b3c114992ae3491cfcb84e878db7ed45',
 'scripts/check_raw23_same_family_saved_20261018.py':'5b1ab9873b2e069f18966a7052856d8a9df8a8c1b0eef4794522a812a54616a7',
}

def enc(x):
    if isinstance(x,Q):return[str(x.a),str(x.b)]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return[enc(v)for v in x]
    return x

def inner(x,g,y):return sum((x[i]*g[i][j]*y[j]for i in range(len(x))for j in range(len(y))),Q())

def rref(rows):
    a=[r[:]for r in rows];piv=[];k=0
    for j in range(len(a[0])):
        i=next((i for i in range(k,len(a))if a[i][j]!=0),None)
        if i is None:continue
        a[k],a[i]=a[i],a[k];z=a[k][j];a[k]=[x/z for x in a[k]]
        for i in range(len(a)):
            if i!=k:
                z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[k])]
        piv.append(j);k+=1
        if k==len(a):break
    return a,piv

def cone(normals,metric):
    """Complete pointed ten-dimensional cone, one inequality per normal."""
    n=len(normals);rays=[];seen={};decisions=[]
    assert len(metric)==10 and all(len(r)==10 for r in normals)
    for active in combinations(range(n),9):
        a,piv=rref([normals[i]for i in active]);rec=dict(active=active,rank=len(piv),ray=None)
        if len(piv)!=9:decisions.append(rec);continue
        free=next(j for j in range(10)if j not in piv)
        x=[Q(int(j==free))for j in range(10)]
        for i,j in enumerate(piv):x[j]=-a[i][free]
        prod=[sum((z*v for z,v in zip(row,x)),Q())for row in normals]
        signs=[z.sign()for z in prod]
        if max(signs)<=0:pass
        elif min(signs)>=0:x=[-z for z in x];prod=[-z for z in prod]
        else:decisions.append(rec);continue
        fac=-next(z for z in prod if z!=0)
        x=[z/fac for z in x];prod=[z/fac for z in prod]
        key=tuple((z.a,z.b)for z in prod)
        if key not in seen:
            index=len(rays);seen[key]=index;norm=inner(x,metric,x)
            rays.append(dict(index=index,coordinates=x,products=prod,norm=norm,norm_sign=norm.sign(),
                zero_facets=[i for i,z in enumerate(prod)if z==0]))
        rec['ray']=seen[key];decisions.append(rec)
    finite=[r for r in rays if r['norm_sign']<0]
    assert finite and all(inner(a['coordinates'],metric,b['coordinates']).sign()<0 for a in finite for b in finite)
    center=[sum((r['coordinates'][i]for r in finite),Q())for i in range(10)]
    centerproducts=[sum((a*b for a,b in zip(row,center)),Q())for row in normals]
    assert inner(center,metric,center).sign()<0 and all(z.sign()<0 for z in centerproducts)
    return dict(wall_count=n,dimension=10,normal_rows=normals,strict_timelike_center=center,
        center_products=centerproducts,complete_active_sets=decisions,rays=rays,
        ray_counts={str(s):sum(r['norm_sign']==s for r in rays)for s in(-1,0,1)})

FINITE_SQUARES=(Q(),Q(1),Q(2),Q('3/2','1/2'),Q(3),Q('5/2','1/2'))
def square_gate(s):
    assert s.sign()>=0
    if (s-4).sign()<0 and s not in FINITE_SQUARES:return'non_Coxeter_finite_angle'
    if (2*s.a).denominator!=1 or(s.a*s.a-5*s.b*s.b).denominator!=1:return'nonintegral_two_cycle'
    conjugate=Q(s.a,-s.b)
    if conjugate.sign()<0 or(conjugate-4).sign()>0:return'conjugate_pair_minor'
    return None

def poles(c,metric,squared_scales):
    finite=[r for r in c['rays']if r['norm_sign']<0]
    outside=[r for r in c['rays']if r['norm_sign']>0];out=[]
    for r in outside:
        x,n=r['coordinates'],r['norm']
        fdots=[inner(x,metric,f['coordinates'])for f in finite]
        preserve=all(z.sign()<0 for z in fdots)
        contacts=[2*z*z/(n*s)for z,s in zip(r['products'],squared_scales)]
        gates=[square_gate(s)for s in contacts]
        passed=preserve and all(g is None for g in gates)
        cuts=[z['index']for z in outside if inner(x,metric,z['coordinates']).sign()>0]
        anchor=next(i for i,s in enumerate(contacts)if s.sign()>0 and(s-4).sign()<0)
        lifts=[]
        for native_square in FINITE_SQUARES[1:]:
            t=native_square/contacts[anchor];within=(t-1).sign()<=0
            conjugate_t=Q(t.a,-t.b)
            conjugate_bound=conjugate_t.sign()>0 and(conjugate_t-1).sign()<=0
            squares=[t*s for s in contacts]if within else None
            failures=[square_gate(s)for s in squares]if within else None
            lifts.append(dict(anchor_native_square=native_square,projection_norm_fraction=t,
                within_physical_norm_bound=within,passive_squared_norm=2*(Q(1)-t),
                conjugate_projection_norm_fraction=conjugate_t,
                within_conjugate_known_span_bound=conjugate_bound,
                scaled_known_contact_squares=squares,gate_failures=failures,
                passes_necessary_lift_gates=within and conjugate_bound and preserve and not any(failures)))
        surviving=[z['projection_norm_fraction']for z in lifts if z['passes_necessary_lift_gates']]
        assert surviving==([Q(1)]if passed else[])
        out.append(dict(ray=r['index'],preserves_all_finite_vertices=preserve,
            finite_vertex_products=fdots,physical_squared_contacts=contacts,gate_failures=gates,
            passes_necessary_polar_gates=passed,cut_spacelike_rays=cuts,
            passive_lift_anchor_wall=anchor,passive_lift_anchor_square=contacts[anchor],
            passive_lift_tests=lifts,passing_projection_norm_fractions=surviving))
    return dict(candidates=out,passing_rays=[r['ray']for r in out if r['passes_necessary_polar_gates']],
        all_spacelike_extreme_rays_tested=len(out),arbitrary_spacelike_combinations_tested=False,
        higher_cycle_arithmetic_sufficiency_claimed=False,arbitrary_passive_lift_relaxation_checked=True,
        unknown_passive_contacts_assumed_zero=False)

def run():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    source=H.read(UNION);H.bind(source);records=[]
    for old in source['cases']:
        gram=H.mat(old['scaled_gram']);scales=[Q(*x)for x in old['squared_scales']]
        basis=next(ids for ids in combinations(range(12),10)if H.inertia(H.principal(gram,ids))==[9,1,0])
        metric=H.principal(gram,basis);assert H.positive([[Q(x.a,-x.b)for x in row]for row in metric])
        normals=[[gram[i][j]for j in basis]for i in range(12)]
        c=cone(normals,metric);assert c['ray_counts']=={'-1':18,'0':0,'1':32}
        assert {tuple(r['zero_facets'])for r in c['rays']if r['norm_sign']<0}==set(map(tuple,old['all_positive_definite_nine_supports']))
        p=poles(c,metric,scales);rec=dict(kinds=old['kinds'],relative_index=old['relative_index'],
            original_basis=basis,metric=metric,original_scaled_gram=gram,original_squared_scales=scales,
            initial_cone=c,initial_polar_tests=p,after_unique_polar=None)
        if old['kinds']==[1,6]:assert p['passing_rays']==[]
        else:
            assert len(p['passing_rays'])==1
            idx=p['passing_rays'][0];r=c['rays'][idx];x,n=r['coordinates'],r['norm']
            assert next(t for t in p['candidates']if t['ray']==idx)['cut_spacelike_rays']==[s['index']for s in c['rays']if s['norm_sign']>0]
            row=[sum((x[i]*metric[i][j]for i in range(10)),Q())for j in range(10)]
            augmented=cone(normals+[row],metric);assert augmented['ray_counts']=={'-1':52,'0':0,'1':30}
            again=poles(augmented,metric,scales+[n/2]);assert again['passing_rays']==[]
            rec['after_unique_polar']=dict(added_initial_ray=idx,added_normal_coordinates=x,
                added_squared_scale=n/2,cone=augmented,polar_tests=again)
        records.append(rec)
        print('CASE',old['kinds'],c['ray_counts'],'initial_passing',p['passing_rays'],
            None if rec['after_unique_polar']is None else rec['after_unique_polar']['cone']['ray_counts'],flush=True)
    paths=set(PINS)|{f'scripts/{TAG}.py',f'results/{REL}.json'}
    return enc(dict(status='PASS_RAW23_EXTREME_RAY_POLAR_METHOD_PREFLIGHT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},cases=records,
        initial_support_sets=660,additional_support_sets=1430,initial_spacelike_poles=96,
        subsequent_spacelike_poles=60,method_passes_required_vertex_retention=True,
        finite_passive_lift_scale_candidates=780,passing_nonzero_passive_lifts=0,
        sequential_extreme_ray_polar_completion_found=False,arbitrary_acute_completion_excluded=False,
        native_pair_tests_are_necessary_only=True,geometric_realization_claimed=False,
        exported_columns=[],old_matrix_changed=False,optimizer_invocations=0,
        scope='Complete exact cone and necessary native tests for sequential truncation whose Lorentz projection is a current spacelike extreme-ray direction. Five exact finite-label scales per direction also exclude every nonzero passive lift without setting passive contacts to zero. The unique same-kind first pole cuts all old outside rays but creates30 new outside rays; none of their projections passes. Mixed kinds have no first pole. Non-extreme projections and arbitrary acute completion remain open; no occurrence charge, capacity, model improvement or chamber exclusion.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['initial_support_sets']+out['additional_support_sets'],flush=True)
