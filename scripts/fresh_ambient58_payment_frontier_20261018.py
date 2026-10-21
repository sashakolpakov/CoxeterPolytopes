#!/usr/bin/env python3
"""Complete changed-dimension diagnostic under the frozen dimension58 rule.

Fresh rank58 additive tables and rank57 ideal convolutions. Both contextual
coupons and every ordinary charge are unchanged. No dimension57 theorem.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement,product
from math import prod
from pathlib import Path
import fresh_a15_d15_ambient59_forest_closure_20261018 as C
import fresh_ambient60_payment_frontier_20261018 as Z
ROOT=C.ROOT;N=58;D=57;BETA=F(377);TAG='fresh_ambient58_payment_frontier_20261018'
COMPLETE='fresh_a15_d15_dimension58_20261018'
FINAL_MDS={C.TAG:'021b4e474b3a15b805788d8736ff732bd102c0f93d8f5921c148a45720a2607c',
 COMPLETE:'c153463b2b9b043eb3f076db1e6504548e0fc001072dd1455c7bd06f2014f845'}
ITEMS=tuple(r for r in C.ITEMS if r<=N)
actual_base=lru_cache(None)(C.actual_base);upper_base=C.upper_base
actual_cross=lru_cache(None)(C.actual_cross);upper_cross=C.upper_cross;context=C.context;value=C.value;item=C.item;enc=C.enc


def fixed_domain():
    out=[]
    def visit(r,left,cs):
        if r==18:
            c=dict(zip(range(9,18),cs))
            key=(2*c[10]+3*sum(c[s]for s in range(11,17))+2*c[17],
                c[11]+c[12],c[13],c[11],c[14],c[9],c[10],c[16],c[15])
            out.append((tuple(cs),N-left,key));return
        for k in range(left//r+1):visit(r+1,left-r*k,cs+[k])
    visit(9,N,[])
    return sorted(out)


def table(key,k5,exclude2):
    ws={r:item(r,key,k5)for r in ITEMS if not exclude2 or r!=2};v=[0]+[-10**9]*N
    for r,w in ws.items():
        for n in range(r,N+1):v[n]=max(v[n],v[n-r]+w)
    last=[0]+[next(r for r,w in ws.items()if r<=n and v[n]==v[n-r]+w)for n in range(1,N+1)]
    return dict(key=list(key),fixed_A5_count=k5,exclude2=exclude2,values=v,last=last)


def check_table(t,key,k5,exclude2):
    assert t['key']==list(key)and t['fixed_A5_count']==k5 and t['exclude2']==exclude2
    v,last=t['values'],t['last'];assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    ws={r:item(r,key,k5)for r in ITEMS if not exclude2 or r!=2};checks=0
    for n in range(1,N+1):
        for r,w in ws.items():
            if r<=n:assert v[n]>=v[n-r]+w;checks+=1
        r=last[n];assert r in ws and r<=n and v[n]==v[n-r]+ws[r]
    return checks


def finite(saved=None):
    rep=C.C.C.P.Z.B.B.representative;types=[t for t in C.C.C.P.Z.B.B.types()if t[1]<=N]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    dom=fixed_domain();assert len(dom)==len({key for cs,u,key in dom})==885
    kd=sorted({(key,k5)for cs,u,key in dom for k5 in range((N-u)//5+1)})
    keys={(key,k5,e)for key,k5 in kd for e in(False,True)}
    tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in saved['tables']}if saved else{k:table(*k)for k in sorted(keys)}
    assert set(tabs)==keys and(saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tabs[k],*k)for k in sorted(keys))
    residual=[rep(('A'if r<=8 else'D',r))for r in ITEMS]
    assert 5 not in ITEMS and all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];expanded=[];identities=choices=0
    for cs,used,key in dom:
        fixedlarge=[('B'if r==9 else'A'if r==17 else'D',r)for r,c in zip(range(9,18),cs)for _ in range(c)]
        assert sum(n for k,n in fixedlarge)==used and key[2]==cs[4]and key[8]==cs[6]
        options5=[]
        for k5 in range((N-used)//5+1):
            fixed=fixedlarge+[('A',5)]*k5;rem=N-used-5*k5
            for r,t in zip(ITEMS,residual):
                assert item(r,key,k5)==upper_base(t)+sum(upper_cross(t,b)for b in fixed);identities+=1
            normal=tabs[key,k5,False];no2=tabs[key,k5,True]
            options=[dict(value=normal['values'][rem],fixed_A2_count=0,reserved_A1_count=0,exclude2=False,old_context=0,new_context=0)]
            for reserve in(1,2):
                oldguard=reserve==2 and bool(cs[1]);newguard=bool(cs[0]and k5)
                if rem<reserve or not(oldguard or newguard)or(reserve==2 and not oldguard):continue
                for j in range((rem-reserve)//2+1):
                    old=4*j*(j-1)if oldguard else 0;new=16*j*k5 if newguard else 0
                    v=reserve*item(1,key,k5)+j*item(2,key,k5)+no2['values'][rem-reserve-2*j]+old+new
                    options.append(dict(value=v,fixed_A2_count=j,reserved_A1_count=reserve,exclude2=True,old_context=old,new_context=new))
            choices+=len(options);win=max(options,key=lambda z:(z['value'],z['reserved_A1_count'],z['fixed_A2_count']))
            j,reserve=win['fixed_A2_count'],win['reserved_A1_count'];tab=tabs[key,k5,win['exclude2']]
            fs=fixed+[('A',2)]*j+[('A',1)]*reserve;rest=rem-2*j-reserve
            while rest:
                r=tab['last'][rest];fs.append(rep(('A'if r<=8 else'D',r)));rest-=r
            fixedcost=sum(upper_base(t)for t in fixed)+sum(upper_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            upper=fixedcost+win['value']
            # This equality verifies contextual payments for the chosen full
            # representative forest, including any relaxed extra A1 roots.
            assert sum(n for k,n in fs)==N and value(fs,True)==upper and value(fs)<=upper
            row=dict(counts=list(cs),key=list(key),fixed_A5_count=k5,forest=sorted(k+str(n)for k,n in fs),
                sixteenths_upper=upper,sixteenths_actual=value(fs),winning_option=win,context_options=options)
            expanded.append(row);options5.append(row)
        rows.append(max(options5,key=lambda z:(z['sixteenths_upper'],z['fixed_A5_count'])))
    maximum=max(r['sixteenths_upper']for r in rows)
    maximizers=[r for r in expanded if r['sixteenths_upper']==maximum]
    above=sorted((r for r in expanded if r['sixteenths_upper']>16*BETA),key=lambda r:(-r['sixteenths_upper'],r['forest']))
    assert len(expanded)==len(kd)==2349 and len(tabs)==4698
    assert len(above)==35 and maximum==6152
    return dict(certified_upper=F(maximum,16),curvature=str(BETA),upper_excess=F(maximum,16)-BETA,
        actual_maximum_claimed=False,actual_attainment_of_global_upper_claimed=False,
        representative_maximizers=maximizers,fixed_large_vectors=len(dom),
        fixed_large_and_A5_vectors=len(expanded),interaction_keys=len(kd),
        tables=[tabs[k]for k in sorted(keys)],simultaneous_type_pairs=len(types)**2,
        residual_item_checks=identities,Bellman_upper_inequalities=checks,attaining_predecessors=len(keys)*N,
        context_options=choices,alternatives=rows,expanded_alternatives=expanded,above_curvature=above,
        representative_additive_attainment_only=True,rank15_actual_costs=dict(A15='100',B15='97',D15='100'),
        rank15_representative_upper='100',fresh_residual_tables=True,
        contextual_witness_multiplicity_not_counted=True,residual_rank17_D17_fixed_rank17_A17=True,
        actual_component_costs={name:F(actual_cross(*pair),16)for name,pair in{
            'A1_A12':(('A',1),('A',12)),'A1_A13':(('A',1),('A',13)),
            'A4_A13':(('A',4),('A',13)),'A4_A11':(('A',4),('A',11))}.items()},
        representative_component_costs={name:F(upper_cross(*pair),16)for name,pair in{
            'A1_D12':(('A',1),('D',12)),'A1_D13':(('A',1),('D',13)),
            'A4_D13':(('A',4),('D',13)),'A4_D11':(('A',4),('D',11))}.items()})


def actual_typed_refinement(ff,saved=None):
    """Exact actual-type maxima in every expanded upper class above377.

    Rank5 and large-component multiplicities are fixed. Actual A1 and A2
    multiplicities are explicit, so both contextual guards are exact.
    The remaining actual types interact only with the fixed components.
    """
    engine=C.C.C.P.Z.B.B;rep=engine.representative
    types=sorted(t for t in engine.types()if t[1]<=N)
    fibers={}
    for t in types:fibers.setdefault(rep(t),[]).append(t)
    residual_reps={rep(('A'if r<=8 else'D',r))for r in ITEMS}
    residual=[t for t in types if rep(t)in residual_reps and t not in [('A',1),('A',2)]]
    assert all(actual_cross(a,b)==0 for a in residual+[('A',1),('A',2)]for b in residual+[('A',1),('A',2)])
    old={(tuple(r['counts']),r['fixed_A5_count']):r for r in saved['classes']}if saved else {}
    expected={(tuple(r['counts']),r['fixed_A5_count'])for r in ff['above_curvature']}
    if saved:assert set(old)==expected and len(old)==len(saved['classes'])
    classes=[];checks=options_count=attainers=0
    for upperrow in ff['above_curvature']:
        cs=tuple(upperrow['counts']);k5=upperrow['fixed_A5_count'];key=(cs,k5)
        groups=[(('B'if r==9 else'A'if r==17 else'D',r),n)for r,n in zip(range(9,18),cs)if n]
        if k5:groups.append((('A',5),k5))
        rem=N-sum(t[1]*n for t,n in groups)
        fixedsets=sorted(tuple(sorted(t for group in choice for t in group))for choice in
            product(*(combinations_with_replacement(fibers[t],n)for t,n in groups)))
        assert len(fixedsets)==len(set(fixedsets))
        savedrows={tuple(map(tuple,r['fixed_actual'])):r for r in old[key]['rows']}if saved else {}
        if saved:assert set(savedrows)==set(fixedsets)and len(savedrows)==len(old[key]['rows'])
        rows=[]
        for fixed in fixedsets:
            fixedcost=sum(actual_base(t)for t in fixed)+sum(actual_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            weights=[actual_base(t)+sum(actual_cross(t,b)for b in fixed)for t in residual]
            if saved:
                prior=savedrows[fixed];assert prior['residual_types']==[list(t)for t in residual]and prior['weights']==weights
                vals,last=prior['values'],prior['last']
            else:
                vals=[0]+[None]*rem
                for t,w in zip(residual,weights):
                    for n in range(t[1],rem+1):
                        if vals[n-t[1]]is not None:
                            z=vals[n-t[1]]+w
                            if vals[n]is None or z>vals[n]:vals[n]=z
                last=[None]+[next((i for i,(t,w)in enumerate(zip(residual,weights))if t[1]<=n and vals[n-t[1]]is not None and vals[n]==vals[n-t[1]]+w),None)for n in range(1,rem+1)]
            assert len(vals)==len(last)==rem+1 and vals[0]==0 and last[0]is None
            for n in range(1,rem+1):
                possible=[]
                for i,(t,w)in enumerate(zip(residual,weights)):
                    if t[1]<=n and vals[n-t[1]]is not None:
                        assert vals[n]is not None and vals[n]>=vals[n-t[1]]+w
                        checks+=1;possible.append(i)
                if possible:
                    i=last[n];assert i in possible and vals[n]==vals[n-residual[i][1]]+weights[i];attainers+=1
                else:assert vals[n]is None and last[n]is None
            w1=actual_base(('A',1))+sum(actual_cross(('A',1),t)for t in fixed)
            w2=actual_base(('A',2))+sum(actual_cross(('A',2),t)for t in fixed)
            co=Counter(fixed);choices=[]
            for a1 in range(rem+1):
                for a2 in range((rem-a1)//2+1):
                    left=rem-a1-2*a2
                    if vals[left]is None:continue
                    oldbonus=4*a2*(a2-1)if co['A',10]and a1>=2 else 0
                    newbonus=16*a2*co['D',5]if co['A',9]and a1>=1 else 0
                    choices.append((fixedcost+a1*w1+a2*w2+vals[left]+oldbonus+newbonus,a1,a2,left))
            options_count+=len(choices);best,a1,a2,left=max(choices)
            forest=list(fixed)+[('A',1)]*a1+[('A',2)]*a2
            while left:
                t=residual[last[left]];forest.append(t);left-=t[1]
            assert sum(n for k,n in forest)==N and value(forest)==best<=upperrow['sixteenths_upper']
            row=dict(fixed_actual=fixed,residual_types=residual,weights=weights,values=vals,last=last,
                sixteenths_actual_maximum=best,attaining_forest=sorted(forest),actual_A1_count=a1,actual_A2_count=a2)
            if saved:assert enc(row)==savedrows[fixed]
            rows.append(row)
        maximum=max(r['sixteenths_actual_maximum']for r in rows)
        classes.append(dict(counts=cs,fixed_A5_count=k5,residual_rank=rem,
            representative_upper=F(upperrow['sixteenths_upper'],16),actual_maximum=F(maximum,16),
            actual_fixed_type_assignments=len(rows),rows=rows,
            actual_obstruction_retained=maximum>16*BETA,
            maximizers=[r['attaining_forest']for r in rows if r['sixteenths_actual_maximum']==maximum]))
    maximum=max(F(r['actual_maximum'])for r in classes)
    witness=[('A',16),('A',16),('D',26)]
    assert F(value(witness),16)==maximum==F(ff['certified_upper'])==F(769,2)
    assert len(classes)==35 and sum(r['actual_fixed_type_assignments']for r in classes)==2330
    assert sum(r['actual_obstruction_retained']for r in classes)==31
    return dict(classes=classes,complete_expanded_upper_classes=len(classes),
        exact_actual_maximum=maximum,actual_attaining_forest=witness,
        actual_fixed_type_assignments=sum(r['actual_fixed_type_assignments']for r in classes),
        actual_residual_tables=sum(r['actual_fixed_type_assignments']for r in classes),
        Bellman_inequalities=checks,attaining_predecessors=attainers,explicit_A1_A2_options=options_count,
        classes_still_above=sum(r['actual_obstruction_retained']for r in classes),
        all_other_classes_at_most_curvature=True,both_context_guards_exact=True,
        actual_type_family_I2star_is_inherited_equivalence=True,
        every_obstructing_actual_forest_enumerated=False,
        scope='Complete exact actual-type maximum within every above-curvature expanded class. All actual large/rank5 types and residual types are covered; explicit A1/A2 counts evaluate both contextual charges exactly. Lists of maximizing witnesses do not enumerate every lower-cost obstructing forest or assert geometric realization.')


def scalar(rs,env):
    h=[1]
    for r in rs:
        out=[0]*(len(h)+r)
        for i,x in enumerate(h):
            for j in range(r+1):out[i+j]+=x
        h=out
    v=prod(r+1 for r in rs)
    assert sum(rs)==D and len(h)==D+1 and h==h[::-1]and sum(h)==v
    g=sum(min(j,D-j)*x for j,x in enumerate(h))
    assert 2*g==D*v-sum(abs(2*j-D)*x for j,x in enumerate(h))
    cumulative=0;deviation=F()
    for j in range(D):
        cumulative+=h[j];deviation+=F(cumulative if 2*j<D else v-cumulative,v)
    assert F(g,v)==F(D,2)-deviation
    curvature=F(5*g-D*v,2);cost=sum(env[r]for r in rs);margin=curvature/v-cost
    return dict(ranks=list(rs),h_polynomial=h,link_vertices=v,G=g,curvature=str(curvature),
        normalized_cost=str(cost),normalized_margin=str(margin),absolute_margin=str(v*margin))


def ideal(previous):
    env={int(k):F(v)for k,v in previous['current_envelopes'].items()}
    assert all(sum(env[r]for r in row['ranks'])==F(row['normalized_cost'])for row in previous['rows'])
    assert all(env[r]==max(F(Z.family_cost(k,r+1),r+1)for k in 'ABCD')for r in range(17,D+1))
    domain=sorted((low,)+p if low else p for low in range(9)for p in Z.high_parts(D-low))
    assert len(domain)==len(set(domain))==1687
    rows=[scalar(rs,env)for rs in domain];summary=Z.summarize(rows)
    active=sorted({n for row in summary['negative_rows']for n in row['ranks']if n>=9})
    thresholds=[Z.threshold(rows,r,{})for r in active]
    scenarios=[Z.conditional(rows,{r:delta for r in active if r>=10})for delta in(F(1,4),F(1,2),F(1))]
    for s in scenarios:
        s['low9_concentration_unchanged']=True;s['new_bound_claimed']=False
    for delta in(F(1,2),F(5,8),F(3,4)):
        scenario=Z.conditional(rows,{9:delta})
        scenario.update(low9_concentration_unchanged=False,low9_concentration_revalidation_required=True,
            terminal_domain_only=True,new_bound_claimed=False)
        scenarios.append(scenario)
    assert summary['negative_count']==31 and summary['zero_count']==0
    assert summary['normalized_minimum']=='-2261/10000'and summary['absolute_minimum']=='-363175'
    assert scenarios[-2]['summary']['negative_count']==scenarios[-2]['summary']['zero_count']==0
    assert scenarios[-2]['summary']['normalized_minimum']=='87/10240'
    for t in thresholds:t['low9_revalidation_required']=t['varied_affine_rank']<=9
    return dict(complete_profiles=len(rows),rows=rows,summary=summary,current_envelopes=env,
        exact_family_tail=[dict(order=r+1,costs={k:Z.family_cost(k,r+1)for k in 'ABCD'},envelope=env[r])for r in range(17,D+1)],
        complete_low9_domain=True,high_factor_concentration_used=False,fresh_convolutions=True,
        single_rank_thresholds=thresholds,conditional_scenarios=scenarios,
        negative_absolute_margin_not_extended_to_arbitrary_profiles=True,
        high_rank_changes_alone_leave_low_only_negative_profile=[3,9,9,9,9,9,9],
        scope='Complete terminal scalar diagnostic with unchanged all-family envelopes. Any hypothetical saving needs compatible geometry and finite accounting; a rank9 change also needs low9 concentration revalidation.')


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in FINAL_MDS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
        assert sha256((ROOT/f'results/{s}.md').read_bytes()).hexdigest()==FINAL_MDS[s]
    prior=docs[C.TAG];assert docs[COMPLETE]['dimension_bound']==58
    dom=fixed_domain();old={(tuple(r['counts']),tuple(r['key']))for r in prior['finite']['alternatives']if sum(n*k for n,k in zip(range(9,18),r['counts']))<=N}
    assert {(cs,key)for cs,u,key in dom}==old
    g0=sum(min(j,D-j)for j in range(D+1));beta=F(5*g0-D*(D+1),2)
    correction=sum(2*j-D for j in range(D+1)if 2*j>D)
    assert(g0,beta,correction)==(812,BETA,841)
    ff=finite(saved['finite']if saved else None);ii=ideal(prior['ideal'])
    print('COMPLETE_UPPER_AND_IDEAL',str(ff['certified_upper']),len(ff['above_curvature']),ii['summary']['negative_count'],flush=True)
    refined=actual_typed_refinement(ff,saved['actual_typed_refinement']if saved else None)
    ff['actual_maximum_claimed']=ff['actual_attainment_of_global_upper_claimed']=True
    ff['exact_actual_maximum']=refined['exact_actual_maximum']
    ff['actual_attaining_forest']=refined['actual_attaining_forest']
    # This changed dimension leaves every item weight unchanged. The fresh
    # tables must equal the matching prefix of the frozen rank59 certificates.
    oldtabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in prior['finite']['tables']}
    for t in ff['tables']:
        k=(tuple(t['key']),t['fixed_A5_count'],t['exclude2'])
        assert t['values']==oldtabs[k]['values'][:N+1]
    ff['frozen_rank59_prefix_comparisons']=len(ff['tables'])
    witnesses=[]
    for fs in [[('D',14)]*2+[('A',15)]*2,[('D',13)]+[('A',15)]*3]:
        assert sum(n for k,n in fs)==N
        v=F(value(fs),16);assert v>BETA
        witnesses.append(dict(forest=fs,cost=v,excess=v-BETA))
    assert [x['cost']for x in witnesses]==[F(381),F(378)]
    ff['initial_actual_obstructions']=witnesses
    paths={Path(__file__),Path(C.__file__),Path(Z.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}'for s in FINAL_MDS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return enc(dict(status='PASS_AMBIENT58_CURRENT_PAYMENT_FRONTIER',sources=sources,
        ambient_dimension=N,link_dimension=D,completed_unrestricted_bound=58,new_dimension_bound=None,
        all_local_weights_unchanged=True,finite_curvature=beta,strict_global_constant=5*correction,
        finite=ff,ideal=ii,actual_typed_refinement=refined,constant_counted_once=True,actual_forest_realizations_claimed=False,
        scope='Complete changed-dimension finite upper and ideal terminal diagnostic under the frozen dimension58 rule. No new payment theorem or dimension57 bound; no extrapolation of negative terminal surplus.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    if not args.verify:assert not path.exists()
    out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite upper',out['finite']['certified_upper'],'above',len(out['finite']['above_curvature']),
        'ideal negatives',out['ideal']['summary']['negative_count'],'ideal min',out['ideal']['summary']['normalized_minimum'],flush=True)
