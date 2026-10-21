#!/usr/bin/env python3
"""Independent rank-based exhaustive clipping of the twelve polar walls.

Starts from the frozen independent complete36cone. Strict signs use its
independent interval engine; new zero incidences use separate nested exact
radical expressions. No new producer arithmetic is imported.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json
import fresh_nonarithmetic_joint_regions_20261016_batch_review as B

ROOT, I, T = B.ROOT, B.I, B.T
TAG='fresh_nonarithmetic_joint_regions_20261016_cone_review'
PRIMARY='results/fresh_nonarithmetic_joint_regions_20261016_cone.json'
BATCH_REVIEW='results/fresh_nonarithmetic_joint_regions_20261016_batch_review.json'
INCIDENCE='results/fresh_nonarithmetic_joint_regions_20261016_incidence.json'


def run():
    target=B.read(PRIMARY);batchreview=B.read(BATCH_REVIEW)
    assert batchreview['status']=='PASS_INDEPENDENT_SIXTEEN_POLARS_AND_SIMULTANEOUS_48_FACETS'
    assert batchreview['sources'][B.BATCH]==B.digest(B.BATCH)
    G=B.geometry();cov,covI,QI,hI=[G[k] for k in ('cov','covI','QI','hI')]
    nodes=[];stats=Counter()
    for active,row in G['old_rays'].items():
        v=[I(*map(F,bounds)) for bounds in row['normalized_coordinate_bounds']]
        nodes.append(dict(mask=sum(1<<j for j in active),v=v,kind=row['kind'],initial=active,parents=None))
    assert len(nodes)==1168
    current=list(range(len(nodes)))

    @lru_cache(None)
    def exact_value(i,wall):
        row=nodes[i]
        if row['mask']>>wall&1:return T()
        if row['initial'] is not None:
            return B.dot(G['exact_old'](row['initial']),cov[wall])
        a,b=row['parents'];k=row['wall']
        return exact_value(b,k)*exact_value(a,wall)-exact_value(a,k)*exact_value(b,wall)

    def sign_value(i,wall):
        row=nodes[i]
        if row['mask']>>wall&1:return 0,None
        val=B.idot(row['v'],covI[wall])
        if val.lo>0:return 1,None
        if val.hi<0:return -1,None
        ex=exact_value(i,wall);sign=B.sg(ex)
        stats['exact_pairing_fallbacks']+=1
        if sign==0:stats['exact_zero_pairings']+=1
        else:stats['exact_nonzero_pairings']+=1
        return sign,dict(active=B.ids(row['mask']),wall=wall,sign=sign,
                         formal_nested_zero=(ex==0))

    stages=[]
    for wall in range(36,48):
        signs=[];fallbacks=[]
        for i in current:
            sign,cert=sign_value(i,wall);signs.append(sign)
            if cert:fallbacks.append(cert)
        neg=[j for j,s in enumerate(signs) if s<0]
        pos=[j for j,s in enumerate(signs) if s>0]
        zero=[j for j,s in enumerate(signs) if s==0]
        out={nodes[current[j]]['mask']:current[j] for j in neg}
        for j in zero:
            old=nodes[current[j]];row=dict(old);row['mask']|=1<<wall
            assert row['mask'] not in out
            out[row['mask']]=len(nodes);nodes.append(row)
        crossings=[];rejects=[];card=Counter()
        for ni in neg:
            a=current[ni]
            for pi in pos:
                b=current[pi];common=nodes[a]['mask']&nodes[b]['mask'];support=B.ids(common)
                card[len(support)]+=1
                if len(support)<8:continue
                minor=B.R.rank_minor([covI[i] for i in support])
                if minor['rank_lower_bound']<8:
                    # A third linearly independent common-null ray proves
                    # common rank<=7. This is a geometric rank certificate,
                    # not the producer's exact-two-face-rays acceptance test.
                    third=None
                    for c in current:
                        if c in (a,b) or nodes[c]['mask']&common!=common:continue
                        certificate=B.R.rank_minor([nodes[j]['v'] for j in (a,b,c)])
                        if certificate['rank_lower_bound']==3:
                            third=c;break
                    assert third is not None,('Unresolved rank',wall,support)
                    rejects.append(dict(negative_support=B.ids(nodes[a]['mask']),
                        positive_support=B.ids(nodes[b]['mask']),common_support=support,
                        third_common_ray=B.ids(nodes[third]['mask']),independent_three_ray_minor=certificate))
                    continue
                va,vb=B.idot(nodes[a]['v'],covI[wall]),B.idot(nodes[b]['v'],covI[wall])
                assert va.hi<0<vb.lo
                alpha,beta=vb/(vb-va),-va/(vb-va)
                assert alpha.lo>0 and beta.lo>0
                v=[alpha*x+beta*y for x,y in zip(nodes[a]['v'],nodes[b]['v'])]
                norm=B.idot(v,B.imv(QI,v));time=B.idot(v,hI)
                kind='hyperideal_rays' if norm.strict_sign()>0 else 'finite_vertices' if time.strict_sign()<0 else 'opposite_sheet_rays'
                mask=common|(1<<wall)
                assert mask not in out
                out[mask]=len(nodes)
                nodes.append(dict(mask=mask,v=v,kind=kind,initial=None,parents=(a,b),wall=wall))
                crossings.append(dict(negative_support=B.ids(nodes[a]['mask']),
                    positive_support=B.ids(nodes[b]['mask']),new_support=B.ids(mask),kind=kind,
                    common_rank8_minor=minor,norm_bounds=norm.saved(),
                    convex_coefficient_bounds=[alpha.saved(),beta.saved()]))
        current=sorted(out.values(),key=lambda i:B.ids(nodes[i]['mask']))
        counts=dict(Counter(nodes[i]['kind'] for i in current))
        stage=dict(wall=wall,negative=len(neg),positive=len(pos),zero=len(zero),counts=counts,
                   total_rays=len(current),crossings=crossings,rank_nonadjacency=rejects,
                   exact_pairing_certificates=fallbacks,
                   common_cardinality_counts={str(k):v for k,v in sorted(card.items())})
        prod=target['clipping_stages'][wall-36]
        assert all(stage[k]==prod[k] for k in ('wall','negative','positive','zero','counts','total_rays'))
        assert len(crossings)==prod['crossing_edges']
        assert len(crossings)+len(rejects)==prod['adjacency_candidates']
        simple=[{k:r[k] for k in ('negative_support','positive_support','new_support','kind')} for r in crossings]
        # The input ray order differs at the first clip; compare sets.
        key=lambda r:(r['negative_support'],r['positive_support'])
        assert sorted(simple,key=key)==sorted(prod['crossings'],key=key)
        stages.append(stage)
        print('independent48 rank clip',wall,counts,'crossings',len(crossings),
              'rank rejected',len(rejects),'zeros',len(zero),flush=True)
    final=[]
    for count,i in enumerate(current,1):
        row=nodes[i];v=row['v'];norm=B.idot(v,B.imv(QI,v));time=B.idot(v,hI)
        kind='hyperideal_rays' if norm.strict_sign()>0 else 'finite_vertices' if time.strict_sign()<0 else 'opposite_sheet_rays'
        assert kind==row['kind']
        for wall,c in enumerate(covI):
            value=B.idot(v,c)
            if row['mask']>>wall&1:assert value.lo<=0<=value.hi
            else:assert value.strict_sign()<0,(count,wall)
        final.append(dict(active=B.ids(row['mask']),kind=kind,norm_bounds=norm.saved(),
                          center_pair_bounds=time.saved(),normalized_coordinate_bounds=[x.saved() for x in v]))
        if count%400==0:print('independent48 full inequalities',count,flush=True)
    census={kind:sorted(r['active'] for r in final if r['kind']==kind)
            for kind in ('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    assert census==target['complete_census']
    counts={k:len(v) for k,v in census.items()}
    assert counts==dict(finite_vertices=172,cusps=0,hyperideal_rays=2109,opposite_sheet_rays=0)
    # Independently recount the finite graph by rank8 common faces. Every
    # finite ray here is simple, so exactly8 common rows have rank8; verify.
    finite=[r for r in final if r['kind']=='finite_vertices'];edges=[]
    for a in range(len(finite)):
        assert len(finite[a]['active'])==9
        for b in range(a+1,len(finite)):
            common=sorted(set(finite[a]['active'])&set(finite[b]['active']))
            if len(common)!=8:continue
            minor=B.R.rank_minor([covI[j] for j in common])
            assert minor['rank_lower_bound']==8
            edges.append(dict(vertices=[a,b],support=common))
    assert len(edges)==395 and 9*len(finite)-2*len(edges)==758
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',PRIMARY,BATCH_REVIEW,
           B.BATCH,B.CONE_REVIEW,B.OLD,B.OLD_REVIEW,INCIDENCE]
    return dict(status='PASS_INDEPENDENT_COMPLETE_48_WALL_CONE_BY_RANK_CLIPPING',
                sources={p:B.digest(p) for p in paths},recursive_integrity=B.R.integrity(paths),
                initial_ray_count=1168,initial_cone_review=B.CONE_REVIEW,
                clipping_stages=stages,arithmetic_counts=dict(stats),counts=counts,
                complete_census=census,final_rays=final,full_final_inequality_checks=48*len(final),
                finite_edges=edges,finite_edge_count=len(edges),open_finite_germs=758,
                proof='Fullrank10 pointedness and strict center inherited and batch-checked. Every crossing accepted by a nonzero8x8 common-covector minor; rejected >=8-wall candidates have three independent common-null rays. Exact zero-ray contacts enlarge masks. Positive convex combinations give exact inherited zero sets and new active rank9; all masks unique. These facts propagate exhaustive halfspace clipping from independently complete36cone through all12walls.',
                arithmetic='Independent512bit outward rational intervals, with initial192bit certified coordinate bounds; nested exact radical arithmetic reconstructs only ambiguous pairing signs. Projective exact expression DAG may use different positive scales from the interval convex DAG; every crossing is scale invariant.',
                scope='Complete cone and finite graph of the fixed48wall batch.172finite,0ideal,2109hyperideal,0opposite-sheet rays;395finiteedges and758open germs. Common strict escape inherited from independently checked batch proves infinite volume. All-label completeness of the16polar bank is separate.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==B.read(str(p.relative_to(ROOT)))
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],flush=True)
