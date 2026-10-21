#!/usr/bin/env python3
"""Independent three-step rank clipping from the reviewed complete49 cone."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import fresh_nonarithmetic_affine_filtered_20261016_cluster_review as C

B=C.B
ROOT,I,T=B.ROOT,B.I,B.T
TAG='fresh_nonarithmetic_affine_filtered_20261016_cone_review'
PRIMARY='results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'
PRIOR='results/fresh_nonarithmetic_common_escape_20261016_cone_review.json'
LOCAL=f'results/{C.TAG}.json'


def run():
    target,prior,local=map(B.read,(PRIMARY,PRIOR,LOCAL))
    assert prior['status']=='PASS_INDEPENDENT_COMPLETE_49_WALL_CONE_AND_RENEWED_ESCAPE'
    assert local['status']=='PASS_INDEPENDENT_FILTERED_THREE_CAPS_AND_52_FACETS'
    assert local['sources'][C.PRIMARY]==B.digest(C.PRIMARY)
    g=C.geometry();covI,QI,hI=[g[k] for k in ('covI','QI','hI')]
    current=[dict(active=row['active'],mask=sum(1<<i for i in row['active']),
                  kind=row['kind'],v=[I(*map(F,x)) for x in row['normalized_coordinate_bounds']])
             for row in prior['final_rays']]
    assert len(current)==2377
    stages=[]
    for wall in range(49,52):
        signs=[B.idot(row['v'],covI[wall]).strict_sign() for row in current]
        negative=[i for i,s in enumerate(signs) if s<0]
        positive=[i for i,s in enumerate(signs) if s>0]
        out={current[i]['mask']:dict(current[i]) for i in negative}
        crossings=[];nonadjacent=[];candidates=0
        for i in negative:
            for j in positive:
                left,right=current[i],current[j]
                common=left['mask']&right['mask'];support=B.ids(common)
                if len(support)<8:continue
                candidates+=1
                minor=B.R.rank_minor([covI[k] for k in support])
                if minor['rank_lower_bound']<8:
                    face=[k for k,row in enumerate(current) if row['mask']&common==common]
                    found=None
                    for triple in combinations(face,3):
                        witness=B.R.rank_minor([current[k]['v'] for k in triple])
                        if witness['rank_lower_bound']>=3:
                            found=dict(rays=[current[k]['active'] for k in triple],rank3_minor=witness)
                            break
                    assert found is not None,(wall,support,'unresolved adjacency')
                    nonadjacent.append(dict(common_support=support,**found))
                    continue
                assert minor['rank_lower_bound']==8
                a,b=B.idot(left['v'],covI[wall]),B.idot(right['v'],covI[wall])
                assert a.hi<0<b.lo
                ca,cb=b/(b-a),-a/(b-a)
                assert ca.lo>0 and cb.lo>0
                v=[ca*x+cb*y for x,y in zip(left['v'],right['v'])]
                norm=B.idot(v,B.imv(QI,v));assert norm.lo>0
                mask=common|(1<<wall);assert mask not in out
                out[mask]=dict(mask=mask,active=B.ids(mask),v=v,kind='hyperideal_rays')
                crossings.append(dict(negative_support=left['active'],positive_support=right['active'],
                                      new_support=B.ids(mask),kind='hyperideal_rays',
                                      common_rank8_minor=minor,norm_bounds=norm.saved(),
                                      positive_convex_coefficients=[ca.saved(),cb.saved()]))
        current=sorted(out.values(),key=lambda row:row['active'])
        stage=dict(wall=wall,negative=len(negative),positive=len(positive),zero=0,
                   adjacency_candidates=candidates,crossings=crossings,
                   nonadjacent_certificates=nonadjacent,counts=dict(Counter(r['kind'] for r in current)))
        expected=next(s for s in target['clipping_stages'] if s['wall']==wall)
        brief=lambda row:{k:row[k] for k in ('negative_support','positive_support','new_support','kind')}
        order=lambda row:(row['negative_support'],row['positive_support'])
        assert sorted(map(brief,crossings),key=order)==sorted(expected['crossings'],key=order)
        assert stage['counts']==expected['counts']
        assert len(negative)==expected['negative'] and len(positive)==expected['positive'] and expected['zero']==0
        stages.append(stage)
        print('independent filtered clip',wall,len(crossings),'nonadjacent',len(nonadjacent),flush=True)
    final=[]
    for count,row in enumerate(current,1):
        v=row['v'];norm=B.idot(v,B.imv(QI,v));time=B.idot(v,hI)
        if norm.strict_sign()>0:kind='hyperideal_rays'
        else:
            assert time.strict_sign()<0
            kind='finite_vertices'
        assert kind==row['kind']
        for i,c in enumerate(covI):
            value=B.idot(v,c)
            if row['mask']>>i&1:assert value.lo<=0<=value.hi
            else:assert value.strict_sign()<0,(count,i)
        final.append(dict(active=row['active'],kind=kind,norm_bounds=norm.saved(),
                          center_pair_bounds=time.saved(),normalized_coordinate_bounds=[x.saved() for x in v]))
        if count%500==0:print('independent filtered52 final inequalities',count,flush=True)
    census={kind:sorted(r['active'] for r in final if r['kind']==kind)
            for kind in ('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    assert census==target['complete_census']
    assert census['finite_vertices']==prior['complete_census']['finite_vertices']
    assert all(max(s)<49 for s in census['finite_vertices'])
    counts={k:len(v) for k,v in census.items()}
    assert counts==dict(finite_vertices=178,cusps=0,hyperideal_rays=2279,opposite_sheet_rays=0)
    escape=target['strict_escape'];E=list(map(B.dec,escape['point']))
    Eiv=[B.iv(x) for x in E]
    norm=B.pair(E,E,g['Q']);time=B.dot(E,g['hcov'])
    assert B.sg(norm)>0 and B.sg(time)<0
    ep=[B.idot(Eiv,c) for c in covI];assert all(x.hi<0 for x in ep)
    assert escape['active'] in census['hyperideal_rays']
    paths=[f'scripts/{TAG}.py',PRIMARY,PRIOR,LOCAL,f'scripts/{C.TAG}.py',C.PRIMARY]
    return dict(status='PASS_INDEPENDENT_COMPLETE_FILTERED52_CONE_AND_ESCAPE',
                sources={p:B.digest(p) for p in paths},inherited_complete49_cone=PRIOR,
                clipping_stages=stages,counts=counts,complete_census=census,final_rays=final,
                full_final_inequality_checks=52*len(final),
                strict_escape=dict(norm_bounds=B.iv(norm).saved(),center_pair_bounds=B.iv(time).saved(),
                                   all52_pairs=[x.saved() for x in ep]),
                proof='Three exhaustive independent clips from the reviewed49cone. Every crossing has a certified rank8 common-covector minor and positive convex coefficients; all new norm and incidence signs are strict. All final inequalities are checked. Active zeros follow by exact inheritance and crossing identities. The178 oldfinitevertices are unchanged.',
                scope='Complete52wall cone with178finite/0ideal/2279exterior/0opposite extreme rays and a strict spacelike escape, hence infinite volume. No finite-volume completion claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
