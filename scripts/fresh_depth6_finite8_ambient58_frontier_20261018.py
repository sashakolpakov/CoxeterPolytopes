#!/usr/bin/env python3
"""Fresh complete rank58 count with depth6 finite8 zeros.

The conditional mode does not bind or assume an unfinished local theorem.
All representative and actual residual tables are new certificates.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement,product
from pathlib import Path
import fresh_ambient58_payment_frontier_20261018 as B
C=B.C;ROOT=B.ROOT;N=58;D=57;BETA=F(377)
TAG='fresh_depth6_finite8_ambient58_frontier_20261018'
LOCAL='fresh_depth6_finite8_zero_transfer_20261018'
BASE_PINS={'scripts/'+B.TAG+'.py':'f7c3d96ca1999c71ca20169cd22b9d15c8809edf611d021346231e7a8ef9eb64',
 'results/'+B.TAG+'.json':'c636d85ba7433f26d3473f29aec8e9e4c122c350a2828775578bcf8c49c41e5e',
 'results/'+B.TAG+'.md':'8067e8ae37e33e33b7362e6083add10dee6e611309bd3d7b29b9d39e86270586'}
LOCAL_PINS={'scripts/'+LOCAL+'.py':'73fbc982ac867d965508d36ea0f30b9560b2cdc8d9783ac9432f73254e2a7e42',
 'results/'+LOCAL+'.json':'a67803a9bea4a0b468a87fcd2ee87044060f5ffd01a3bb9776338aa7c429ec8c',
 'results/'+LOCAL+'.md':'294d550819cd749dac05194079d61f09f324802c76e5438578cfdbe332e19bb3'}
ITEMS=B.ITEMS;enc=B.enc
actual_cross=B.actual_cross;upper_cross=B.upper_cross;context=B.context


def deleted(t):
    k,n=t
    return min(4,max(0,n-{'A':20,'B':21,'D':22}[k]))if k in 'ABD'else 0


@lru_cache(None)
def actual_base(t):return B.actual_base(t)-16*deleted(t)
def upper_base(t):return B.upper_base(t)-16*deleted(t)
def item(r,key,k5):return B.item(r,key,k5)-16*deleted(('D',r))
def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)


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
    assert above and maximum>16*BETA
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
    witness=[('A',1)]*7+[('D',5)]*3+[('D',12)]*3
    assert sum(n for k,n in witness)==N
    assert F(value(witness),16)==maximum==F(ff['certified_upper'])==F(3063,8)
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


def bind(conditional):
    pins=dict(BASE_PINS)
    if not conditional:
        assert len(LOCAL_PINS)==3,'Local theorem is not frozen/bound yet'
        pins.update(LOCAL_PINS)
    for p,h in pins.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    prior=json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    assert prior['status']=='PASS_AMBIENT58_CURRENT_PAYMENT_FRONTIER'
    assert prior['ambient_dimension']==N and prior['completed_unrestricted_bound']==58
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    if not conditional:
        local=json.loads((ROOT/f'results/{LOCAL}.json').read_text())
        assert local['status']=='PASS_DEPTH6_FINITE8_ZERO_LOCAL_TRANSFER'
        assert (local['feature_depth'],local['finite_distance'],local['previous_feature_depth'],local['new_weight'])==(6,8,8,'0')
        assert local['new_coupons'] is False and local['all_ideal_weights_unchanged'] is True and local['other_finite_weights_unchanged'] is True
        assert local['triangle_minimum']==2 and local['quadrilateral_minimum']==1
        rows={(r['kind'],r['rank']):r for r in local['counts']['rows']}
        assert len(rows)==len(local['counts']['rows'])
        for k in 'ABD':
            for n in range(1,N+1):
                if n<{'A':1,'B':2,'D':4}[k]:
                    assert deleted((k,n))==0
                    continue
                row=rows[k,n]
                a={tuple(p)for p in row['depth6_pairs']};b={tuple(p)for p in row['depth8_pairs']}
                assert b<=a and len(a-b)==row['hypothetical_additional_deletions']==deleted((k,n))
        for p,h in local['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return prior,pins


def run(saved=None,conditional=False):
    prior,pins=bind(conditional)
    dom=fixed_domain()
    assert {(cs,key)for cs,u,key in dom}=={(tuple(r['counts']),tuple(r['key']))for r in prior['finite']['alternatives']}
    assert all(deleted((k,n))==0 for k in 'ABD'for n in range(1,18))
    g0=sum(min(j,D-j)for j in range(D+1));beta=F(5*g0-D*(D+1),2)
    correction=sum(2*j-D for j in range(D+1)if 2*j>D)
    assert (g0,beta,correction)==(812,BETA,841)
    ff=finite(saved['finite']if saved else None)
    print('REPRESENTATIVE',ff['certified_upper'],'above',len(ff['above_curvature']),flush=True)
    refined=actual_typed_refinement(ff,saved['actual_typed_refinement']if saved else None)
    ff['actual_maximum_claimed']=ff['actual_attainment_of_global_upper_claimed']=True
    ff['exact_actual_maximum']=refined['exact_actual_maximum'];ff['actual_attaining_forest']=refined['actual_attaining_forest']
    # These are changed tables; old tables are used only to check monotone
    # cost decrease and the identical complete domain, never as new optima.
    oldtabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in prior['finite']['tables']}
    changed=0
    for t in ff['tables']:
        k=(tuple(t['key']),t['fixed_A5_count'],t['exclude2']);o=oldtabs[k]
        assert all(x<=y for x,y in zip(t['values'],o['values']))
        changed+=t['values']!=o['values']
    oldclasses={(tuple(r['counts']),r['fixed_A5_count']):r for r in prior['actual_typed_refinement']['classes']}
    assert all((tuple(r['counts']),r['fixed_A5_count'])in oldclasses for r in refined['classes'])
    ff['changed_representative_tables']=changed
    ff['old_attainer_new_cost']=F(value([('A',16),('A',16),('D',26)]),16)
    assert ff['old_attainer_new_cost']==F(761,2)
    ii=prior['ideal'];assert ii['complete_profiles']==1687 and ii['summary']['negative_count']==31
    paths={Path(__file__),Path(C.__file__),Path(B.__file__)}|{ROOT/p for p in pins}
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return enc(dict(status='PASS_CONDITIONAL_DEPTH6_FINITE8_AMBIENT58_COUNT'if conditional else'PASS_DEPTH6_FINITE8_AMBIENT58_FRONTIER',
        sources=sources,geometry_inputs_final=not conditional,ambient_dimension=N,link_dimension=D,
        completed_unrestricted_bound=58,new_dimension_bound=None,finite_curvature=beta,strict_global_constant=5*correction,
        changed_weights='Only whole finite A/B/D distance8 units newly deep at threshold6 are deleted.',
        component_deletions=[dict(type=k,rank=n,deleted=deleted((k,n)))for k in 'ABD'for n in range(1,N+1)],
        ordinary_coupons_unchanged=True,contextual_coupons_unchanged=True,all_ideal_weights_unchanged=True,
        finite=ff,actual_typed_refinement=refined,ideal=ii,ideal_arithmetic_inherited_without_reconvolution=True,
        constant_counted_once=True,actual_forest_realizations_claimed=False,
        scope='Fresh complete representative and actual-type rank58 counting certificates. The optional conditional mode assumes the stated deletion only, with no geometric conclusion. Final mode binds the reviewed local theorem. No dimension57 bound; unchanged negative ideal margins are not extended to arbitrary profiles.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--conditional',action='store_true');args=ap.parse_args()
    path=Path('/private/tmp/depth6_finite8_ambient58_conditional_20261018.json')if args.conditional else ROOT/f'results/{TAG}.json'
    saved=json.loads(path.read_text())if args.verify else None
    if not args.verify:assert not path.exists(),path
    out=run(saved,args.conditional)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['certified_upper'],'actual',out['actual_typed_refinement']['exact_actual_maximum'],
        'upper classes',len(out['finite']['above_curvature']),'actual classes',out['actual_typed_refinement']['classes_still_above'],flush=True)
