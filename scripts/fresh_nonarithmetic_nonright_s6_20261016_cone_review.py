#!/usr/bin/env python3
"""Independent complete55 cone from two exhaustive clips of reviewed53."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
import fresh_nonarithmetic_nonright_s6_20261016_caps_review as C

B,ROOT,I=C.B,C.ROOT,C.I
TAG='fresh_nonarithmetic_nonright_s6_20261016_cone_review'
PRIMARY='results/fresh_nonarithmetic_nonright_s6_20261016_cone.json'
PRIOR='results/fresh_nonarithmetic_escape52_20261016_cone_review.json'
LOCAL=f'results/{C.TAG}.json'


def run():
    target,prior,local=map(B.read,(PRIMARY,PRIOR,LOCAL))
    assert target['status']=='PASS_COMPLETE_EXACT_55_WALL_CONE_AFTER_TWO_NONRIGHT_S6_CUTTERS'
    assert prior['status']=='PASS_INDEPENDENT_COMPLETE_53_WALL_CONE_AND_RESIDUAL_H3_FACE'
    assert local['status']=='PASS_INDEPENDENT_TWO_NONRIGHT_S6_CUTTERS_AND_ZERO_NEW_FINITE_VERTICES'
    assert local['sources'][C.PRIMARY]==B.digest(C.PRIMARY)
    g=C.geometry();cov,Q,h=g['covI'],g['QI'],g['hI']
    rays=[dict(active=row['active'],mask=sum(1<<j for j in row['active']),kind=row['kind'],
               v=[I(*map(F,x)) for x in row['normalized_coordinate_bounds']]) for row in prior['final_rays']]
    assert len(rays)==2511
    oldcut=[set() for _ in range(2)]
    for idx,row in enumerate(rays):
        for k in range(2):
            if B.idot(row['v'],cov[53+k]).strict_sign()>0:oldcut[k].add(idx)
    assert len(oldcut[0]|oldcut[1])==target['old53_cut_union_count']==257
    stages=[]
    for wall,expectedcounts in ((53,(2388,123,217)),(54,(2462,143,226))):
        signs=[B.idot(row['v'],cov[wall]).strict_sign() for row in rays]
        negative=[i for i,s in enumerate(signs) if s<0]
        positive=[i for i,s in enumerate(signs) if s>0]
        assert len(negative)+len(positive)==len(rays)
        assert (len(negative),len(positive))==expectedcounts[:2]
        assert all(rays[i]['kind']=='hyperideal_rays' for i in positive)
        out={rays[i]['mask']:dict(rays[i]) for i in negative};crossings=[]
        for i in negative:
            for j in positive:
                a,b=rays[i],rays[j];common=a['mask']&b['mask'];support=B.ids(common)
                if len(support)<8:continue
                minor=B.R.rank_minor([cov[k] for k in support]);assert minor['rank_lower_bound']==8,support
                ga,gb=B.idot(a['v'],cov[wall]),B.idot(b['v'],cov[wall]);assert ga.hi<0<gb.lo
                ca,cb=gb/(gb-ga),-ga/(gb-ga);assert ca.lo>0 and cb.lo>0
                v=[ca*x+cb*y for x,y in zip(a['v'],b['v'])]
                q=B.idot(v,B.imv(Q,v));assert q.strict_sign()>0
                mask=common|(1<<wall);assert mask not in out
                out[mask]=dict(active=B.ids(mask),mask=mask,kind='hyperideal_rays',v=v)
                crossings.append(dict(negative_support=a['active'],positive_support=b['active'],
                    new_support=B.ids(mask),kind='hyperideal_rays',rank8_minor=minor,norm=q.saved(),
                    positive_convex_coefficients=[ca.saved(),cb.saved()]))
        assert len(crossings)==expectedcounts[2]
        expected=next(x for x in target['clipping_stages'] if x['wall']==wall)
        keys=('negative_support','positive_support','new_support','kind')
        brief=lambda row:{k:row[k] for k in keys}
        order=lambda row:(row['negative_support'],row['positive_support'])
        assert sorted(map(brief,crossings),key=order)==sorted(expected['crossings'],key=order)
        rays=sorted(out.values(),key=lambda x:x['active'])
        stages.append(dict(wall=wall,negative=len(negative),positive=len(positive),crossings=crossings,
            counts=dict(Counter(x['kind'] for x in rays))))
        print('independent55 clip',wall,len(rays),'rays',flush=True)
    final=[]
    for count,row in enumerate(rays,1):
        v=row['v'];q=B.idot(v,B.imv(Q,v));time=B.idot(v,h)
        kind='finite_vertices' if q.strict_sign()<0 else 'hyperideal_rays';assert kind==row['kind']
        if kind=='finite_vertices':assert time.hi<0
        for i,c in enumerate(cov):
            z=B.idot(v,c)
            if row['mask']>>i&1:assert z.lo<=0<=z.hi
            else:assert z.hi<0,(row['active'],i)
        final.append(dict(active=row['active'],kind=kind,norm_bounds=q.saved(),center_pair_bounds=time.saved(),
            normalized_coordinate_bounds=[x.saved() for x in v]))
        if count%500==0:print('independent55 all final inequalities',count,flush=True)
    census={kind:sorted(row['active'] for row in final if row['kind']==kind)
        for kind in ('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    assert census==target['complete_census']
    counts={k:len(v) for k,v in census.items()}
    assert counts==dict(finite_vertices=179,cusps=0,hyperideal_rays=2509,opposite_sheet_rays=0)
    assert census['finite_vertices']==prior['complete_census']['finite_vertices']
    E=list(map(B.dec,target['strict_escape']['point']))
    qE=B.pair(E,E,g['Q']);hE=B.dot(E,g['hcov']);assert B.sg(qE)>0 and B.sg(hE)<0
    EI=[B.iv(x) for x in E];ep=[B.idot(EI,c) for c in cov];assert all(x.hi<0 for x in ep)
    assert target['strict_escape']['active'] in census['hyperideal_rays']
    paths=[f'scripts/{TAG}.py',PRIMARY,PRIOR,LOCAL,f'scripts/{C.TAG}.py',C.PRIMARY]
    return dict(status='PASS_INDEPENDENT_COMPLETE_55_WALL_CONE_AND_RENEWED_ESCAPE',
        sources={p:B.digest(p) for p in paths},counts=counts,complete_census=census,
        old53_cut_union_count=len(oldcut[0]|oldcut[1]),clipping_stages=stages,final_rays=final,
        full_final_inequality_checks=55*len(final),
        strict_escape=dict(norm_bounds=B.iv(qE).saved(),center_pair_bounds=B.iv(hE).saved(),all55_pairs=[x.saved() for x in ep]),
        proof='Independent reviewed53cone and two exhaustive clips, with443rank8-certified crossings. All147840finalinequalities checked. Complete cone has179finite,0ideal,2509exterior and0opposite rays; alloldfinitevertices persist. Bothcaps jointly remove257oldexteriorrays, but create443crossingrays across the two stages. Renewed exactstrictescape proves infinitevolume.',
        scope='Complete55wall cone and renewed infinite-volume certificate; removingoldrays is not a finite-volume claim. No new construction dimension.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
