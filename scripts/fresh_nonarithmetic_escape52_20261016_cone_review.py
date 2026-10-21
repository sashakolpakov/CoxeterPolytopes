#!/usr/bin/env python3
"""Independent complete53-wall cone by one rank-certified clip of reviewed52."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
import fresh_nonarithmetic_escape52_20261016_cap_review as C

B,ROOT,I=C.B,C.ROOT,C.I
TAG='fresh_nonarithmetic_escape52_20261016_cone_review'
PRIMARY='results/fresh_nonarithmetic_escape52_20261016_cone.json'
PRIOR='results/fresh_nonarithmetic_affine_filtered_20261016_cone_review.json'
LOCAL=f'results/{C.TAG}.json'


def run():
    target,prior,local=map(B.read,(PRIMARY,PRIOR,LOCAL))
    assert prior['status']=='PASS_INDEPENDENT_COMPLETE_FILTERED52_CONE_AND_ESCAPE'
    assert local['status']=='PASS_INDEPENDENT_53RD_CAP_AND_COMPLETE_NEW_FINITE_VERTEX_CENSUS'
    assert local['sources'][C.PRIMARY]==B.digest(C.PRIMARY)
    g=C.geometry();cov,Q,h=g['covI'],g['QI'],g['hI']
    rays=[dict(active=row['active'],mask=sum(1<<j for j in row['active']),
               kind=row['kind'],v=[I(*map(F,x)) for x in row['normalized_coordinate_bounds']])
          for row in prior['final_rays']]
    assert len(rays)==2457
    signs=[B.idot(row['v'],cov[52]).strict_sign() for row in rays]
    negative=[i for i,s in enumerate(signs) if s<0]
    positive=[i for i,s in enumerate(signs) if s>0]
    assert len(negative)==2307 and len(positive)==150
    assert all(rays[i]['kind']=='hyperideal_rays' for i in positive)
    out={rays[i]['mask']:dict(rays[i]) for i in negative}
    crossings=[]
    for i in negative:
        for j in positive:
            a,b=rays[i],rays[j]
            common=a['mask']&b['mask'];support=B.ids(common)
            if len(support)<8:continue
            minor=B.R.rank_minor([cov[k] for k in support])
            assert minor['rank_lower_bound']==8,support
            ga,gb=B.idot(a['v'],cov[52]),B.idot(b['v'],cov[52])
            assert ga.hi<0<gb.lo
            ca,cb=gb/(gb-ga),-ga/(gb-ga)
            assert ca.lo>0 and cb.lo>0
            v=[ca*x+cb*y for x,y in zip(a['v'],b['v'])]
            q=B.idot(v,B.imv(Q,v))
            kind='finite_vertices' if q.strict_sign()<0 else 'hyperideal_rays'
            if kind=='finite_vertices':assert B.idot(v,h).hi<0
            mask=common|(1<<52);assert mask not in out
            out[mask]=dict(active=B.ids(mask),mask=mask,kind=kind,v=v)
            crossings.append(dict(negative_support=a['active'],positive_support=b['active'],
                                  new_support=B.ids(mask),kind=kind,rank8_minor=minor,
                                  norm=q.saved(),positive_convex_coefficients=[ca.saved(),cb.saved()]))
    assert len(crossings)==204
    expected=target['clipping_stages'][-1]
    keys=('negative_support','positive_support','new_support','kind')
    brief=lambda row:{k:row[k] for k in keys}
    order=lambda row:(row['negative_support'],row['positive_support'])
    assert sorted(map(brief,crossings),key=order)==sorted(expected['crossings'],key=order)
    final=[]
    for count,row in enumerate(sorted(out.values(),key=lambda x:x['active']),1):
        v=row['v'];q=B.idot(v,B.imv(Q,v));time=B.idot(v,h)
        kind='finite_vertices' if q.strict_sign()<0 else 'hyperideal_rays'
        assert kind==row['kind']
        if kind=='finite_vertices':assert time.hi<0
        for i,c in enumerate(cov):
            z=B.idot(v,c)
            if row['mask']>>i&1:assert z.lo<=0<=z.hi
            else:assert z.hi<0,(row['active'],i)
        final.append(dict(active=row['active'],kind=kind,norm_bounds=q.saved(),
                          center_pair_bounds=time.saved(),normalized_coordinate_bounds=[x.saved() for x in v]))
        if count%500==0:print('independent53 final inequalities',count,flush=True)
    census={kind:sorted(row['active'] for row in final if row['kind']==kind)
            for kind in ('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    assert census==target['complete_census']
    counts={k:len(v) for k,v in census.items()}
    assert counts==dict(finite_vertices=179,cusps=0,hyperideal_rays=2332,opposite_sheet_rays=0)
    newfinite=[s for s in census['finite_vertices'] if 52 in s]
    assert newfinite==[C.SUPPORT+[52]]
    assert [s for s in census['finite_vertices'] if 52 not in s]==prior['complete_census']['finite_vertices']
    E=list(map(B.dec,target['strict_escape']['point']))
    qE=B.pair(E,E,g['Q']);hE=B.dot(E,g['hcov'])
    assert B.sg(qE)>0 and B.sg(hE)<0
    EI=[B.iv(x) for x in E];ep=[B.idot(EI,c) for c in cov]
    assert all(x.hi<0 for x in ep)
    assert target['strict_escape']['active'] in census['hyperideal_rays']
    # The actual codimension6 face is read from the complete ambient cone;
    # this is not a claim that its induced dihedral angles are Coxeter.
    S6={0,1,2,3,5,6}
    face={kind:[s for s in sets if S6<=set(s)] for kind,sets in census.items()}
    assert {k:len(v) for k,v in face.items()}==dict(finite_vertices=13,cusps=0,hyperideal_rays=7,opposite_sheet_rays=0)
    paths=[f'scripts/{TAG}.py',PRIMARY,PRIOR,LOCAL,f'scripts/{C.TAG}.py',C.PRIMARY]
    return dict(status='PASS_INDEPENDENT_COMPLETE_53_WALL_CONE_AND_RESIDUAL_H3_FACE',
                sources={p:B.digest(p) for p in paths},counts=counts,complete_census=census,
                negative_rays=len(negative),positive_rays=len(positive),crossings=crossings,
                final_rays=final,full_final_inequality_checks=53*len(final),
                strict_escape=dict(norm_bounds=B.iv(qE).saved(),center_pair_bounds=B.iv(hE).saved(),
                                   all53_pairs=[x.saved() for x in ep]),
                actual_H3_face=dict(support=sorted(S6),complete_census=face),
                proof='Independent reviewed52cone plus an exhaustive53rd clip. All204crossings have certifiedrank8common minors andpositiveconvex coefficients. Complete finalcheck has179finite/0ideal/2332exterior rays; old178finitevertices persist andone newvertex agreeswithindependentlocalcensus. Exact strictescape provesinfinitevolume. The H4+A2 supportface has13finite and7exterior rays.',
                scope='Complete53-wall cone and inducedH3face census; projectedfacetangles ontheH3face arenotassumedCoxeter. Nofinitevolume/newdimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
