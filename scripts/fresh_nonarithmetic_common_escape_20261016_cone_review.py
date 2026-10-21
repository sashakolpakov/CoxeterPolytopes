#!/usr/bin/env python3
"""Independent one-step rank clipping from the complete reviewed48-wall cone."""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json
import fresh_nonarithmetic_common_escape_20261016_cap_review as C

B=C.B
ROOT,I,T=B.ROOT,B.I,B.T
TAG='fresh_nonarithmetic_common_escape_20261016_cone_review'
PRIMARY='results/fresh_nonarithmetic_common_escape_20261016_cone.json'
PRIOR='results/fresh_nonarithmetic_joint_regions_20261016_cone_review.json'
CAP_REVIEW=f'results/{C.TAG}.json'
INCIDENCE='results/fresh_nonarithmetic_common_escape_20261016_incidence.json'


def run():
    target,prior,capreview=map(B.read,(PRIMARY,PRIOR,CAP_REVIEW))
    assert prior['status']=='PASS_INDEPENDENT_COMPLETE_48_WALL_CONE_BY_RANK_CLIPPING'
    assert capreview['status']=='PASS_INDEPENDENT_49TH_NONPOLAR_CAP_AND_ALL_FACETS'
    assert capreview['sources'][C.PRIMARY]==B.digest(C.PRIMARY)
    G=C.geometry();Q,cov,covI,QI,hI=[G[k] for k in ('Q','cov','covI','QI','hI')]
    initial={tuple(row['active']):row for row in prior['final_rays']}
    assert len(initial)==2281
    routes={}
    for stage in prior['clipping_stages']:
        wall=stage['wall']
        for row in stage['crossings']:
            key=tuple(row['new_support']);assert key not in routes
            routes[key]=('cross',wall,tuple(row['negative_support']),tuple(row['positive_support']))
        for row in stage['exact_pairing_certificates']:
            if row['sign']==0:
                key=tuple(sorted(row['active']+[wall]));assert key not in routes
                routes[key]=('zero',tuple(row['active']))

    @lru_cache(None)
    def exact48(support):
        if support in G['old_rays']:return G['exact_old'](support)
        route=routes[support]
        if route[0]=='zero':return exact48(route[1])
        _,wall,sa,sb=route
        a,b=exact48(sa),exact48(sb)
        ga,gb=B.dot(a,cov[wall]),B.dot(b,cov[wall])
        assert B.sg(ga)<0<B.sg(gb)
        return tuple(gb*x-ga*y for x,y in zip(a,b))

    nodes=[];fallback=[]
    for support,row in sorted(initial.items()):
        v=[I(*map(F,x)) for x in row['normalized_coordinate_bounds']]
        value=B.idot(v,covI[48])
        if value.lo>0:sign=1
        elif value.hi<0:sign=-1
        else:
            exact=B.dot(exact48(support),cov[48]);sign=B.sg(exact)
            fallback.append(dict(active=list(support),sign=sign,formal_zero=(exact==0)))
        nodes.append(dict(active=list(support),mask=sum(1<<j for j in support),v=v,kind=row['kind'],sign=sign))
    negative=[i for i,row in enumerate(nodes) if row['sign']<0]
    positive=[i for i,row in enumerate(nodes) if row['sign']>0]
    zero=[i for i,row in enumerate(nodes) if row['sign']==0]
    assert [len(negative),len(positive),len(zero)]==[2180,94,7]
    out={nodes[i]['mask']:dict(nodes[i]) for i in negative}
    for i in zero:
        row=dict(nodes[i]);row['mask']|=1<<48;row['active']+= [48]
        assert row['mask'] not in out;out[row['mask']]=row
    crossings=[];candidate_count=0
    for i in negative:
        for j in positive:
            left,right=nodes[i],nodes[j]
            common=left['mask']&right['mask'];support=B.ids(common)
            if len(support)<8:continue
            candidate_count+=1
            minor=B.R.rank_minor([covI[k] for k in support])
            assert minor['rank_lower_bound']==8
            a,b=B.idot(left['v'],covI[48]),B.idot(right['v'],covI[48])
            assert a.hi<0<b.lo
            ca,cb=b/(b-a),-a/(b-a)
            assert ca.lo>0 and cb.lo>0
            v=[ca*x+cb*y for x,y in zip(left['v'],right['v'])]
            norm=B.idot(v,B.imv(QI,v));time=B.idot(v,hI)
            kind='hyperideal_rays' if norm.strict_sign()>0 else 'finite_vertices' if time.strict_sign()<0 else 'opposite_sheet_rays'
            mask=common|(1<<48);assert mask not in out
            out[mask]=dict(mask=mask,active=B.ids(mask),v=v,kind=kind)
            crossings.append(dict(negative_support=left['active'],positive_support=right['active'],
                                  new_support=B.ids(mask),kind=kind,common_rank8_minor=minor,
                                  norm_bounds=norm.saved(),convex_coefficients=[ca.saved(),cb.saved()]))
    assert candidate_count==len(crossings)==190
    saved_stage=target['clipping_stages'][-1]
    brief=lambda row:{k:row[k] for k in ('negative_support','positive_support','new_support','kind')}
    order=lambda row:(row['negative_support'],row['positive_support'])
    assert sorted(map(brief,crossings),key=order)==sorted(saved_stage['crossings'],key=order)
    print('PASS independent49 rank clip,190 crossing minors and7 exact zeros',flush=True)
    final=[]
    for count,row in enumerate(sorted(out.values(),key=lambda r:r['active']),1):
        v=row['v'];norm=B.idot(v,B.imv(QI,v));time=B.idot(v,hI)
        kind='hyperideal_rays' if norm.strict_sign()>0 else 'finite_vertices' if time.strict_sign()<0 else 'opposite_sheet_rays'
        assert kind==row['kind']
        for i,c in enumerate(covI):
            value=B.idot(v,c)
            if row['mask']>>i&1:assert value.lo<=0<=value.hi
            else:assert value.strict_sign()<0,(count,i)
        final.append(dict(active=row['active'],kind=kind,norm_bounds=norm.saved(),
                          center_pair_bounds=time.saved(),normalized_coordinate_bounds=[x.saved() for x in v]))
        if count%500==0:print('independent49 full inequalities',count,flush=True)
    census={kind:sorted(row['active'] for row in final if row['kind']==kind)
            for kind in ('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    assert census==target['complete_census']
    counts={k:len(v) for k,v in census.items()}
    assert counts==dict(finite_vertices=178,cusps=0,hyperideal_rays=2199,opposite_sheet_rays=0)
    finite=[r for r in final if r['kind']=='finite_vertices'];edges=[]
    for i,a in enumerate(finite):
        assert len(a['active'])==9
        for j in range(i+1,len(finite)):
            common=sorted(set(a['active'])&set(finite[j]['active']))
            if len(common)!=8:continue
            minor=B.R.rank_minor([covI[k] for k in common])
            assert minor['rank_lower_bound']==8
            edges.append(dict(vertices=[i,j],support=common))
    assert len(edges)==409 and 9*len(finite)-2*len(edges)==784
    # Exact strict escape is checked in the independent nested field. Its
    # normalization need not equal the review's interval ray normalization.
    esc=target['strict_escape'];E=list(map(B.dec,esc['point']))
    ray=list(map(B.dec,esc['ray']));eps=T(K.convert(F(esc['epsilon'])))
    assert E==[v+eps*h for v,h in zip(ray,G['H'])]
    assert esc['active'] in census['hyperideal_rays']
    rp=[B.dot(ray,c) for c in cov]
    assert [i for i,v in enumerate(rp) if B.sg(v)==0]==esc['active']
    assert all(B.sg(v)<=0 for v in rp)
    norm=B.pair(E,E,Q);time=B.dot(E,G['hcov'])
    assert B.sg(norm)>0 and B.sg(time)<0
    ep=[B.dot(E,c) for c in cov];assert all(B.sg(v)<0 for v in ep)
    assert norm==B.dec(esc['norm']) and ep==list(map(B.dec,esc['pairings']))
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'scripts/{B.TAG}.py',PRIMARY,PRIOR,CAP_REVIEW,C.PRIMARY,INCIDENCE]
    return dict(status='PASS_INDEPENDENT_COMPLETE_49_WALL_CONE_AND_RENEWED_ESCAPE',
                sources={p:B.digest(p) for p in paths},recursive_integrity=B.R.integrity(paths),
                inherited_complete_48_cone=PRIOR,initial_rays=2281,
                clip_counts=dict(negative=2180,positive=94,zero=7,crossings=190),
                exact_pairing_fallbacks=fallback,crossings=crossings,
                counts=counts,complete_census=census,final_rays=final,
                full_final_inequality_checks=49*len(final),finite_edges=edges,
                finite_edge_count=len(edges),open_finite_germs=784,
                strict_escape=dict(active=esc['active'],epsilon=esc['epsilon'],norm_bounds=B.iv(norm).saved(),
                                   center_pair_bounds=B.iv(time).saved(),all_pair_bounds=[B.iv(v).saved() for v in ep]),
                proof='Single exhaustive clip of the independently complete48cone. Every190 crossing accepted by a nonzero8x8 common-covector minor. Fewer-than8 filters are rank-safe; seven zero contacts reconstructed exactly through independent prior clipping provenance. Positive convex combinations preserve zero sets and rank9 unique masks. All116473 final inequalities and409 finite edges checked.',
                scope='Complete49wall cone and finite graph, with178finite/0ideal/2199hyperideal/0opposite rays. Renewed exact strictescape proves infinitevolume; no new construction dimension.')


K=B.K
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==B.read(str(out.relative_to(ROOT)))
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
