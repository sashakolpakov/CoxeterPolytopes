#!/usr/bin/env python3
"""Complete changed ambient58 count for rank16 transfers and the witness guard.

Fresh representative tables and complete actual fibers in every upper class
above curvature. The A1/D5 witness guard is exact in actual refinement.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from functools import lru_cache
from itertools import combinations_with_replacement,product
from pathlib import Path
import fresh_quarter_center_rank13_ambient58_frontier_20261018 as OLD
import fresh_a16_d16_a4_a12_fractional_transfer_20261018 as L
C=OLD.C;W=L.W;ROOT=OLD.ROOT;N=58;SCALE=64;BETA=F(377)
ITEMS=OLD.ITEMS;TYPES=OLD.TYPES;ENGINE=OLD.ENGINE;enc=OLD.enc
TAG='fresh_a16_d16_witness_context_ambient58_frontier_20261018'
STEMS=(OLD.TAG,L.TAG,W.TAG)
FINAL_REPORTS={L.TAG:'2701281f8d5d3660a7ff792c5829c59c8e30b5f168ddea1aa33d509def66271f',W.TAG:'c8acf530f290cd3ea1ebae1bcf88af06b82f17423952ca0cb0d3a0819ca0e004'}
representative=OLD.representative;fixed_domain=OLD.fixed_domain;fixed_components=OLD.fixed_components
@lru_cache(None)
def actual_base(t):return L.actual_base(t)
@lru_cache(None)
def upper_base(t):return L.upper_base(t)
@lru_cache(None)
def actual_cross(a,b):return L.actual_cross(a,b)
@lru_cache(None)
def upper_cross(a,b):return L.upper_cross(a,b)
def context(fs,upper=False):return OLD.context(fs,True)if upper else W.context(fs)
@lru_cache(None)
def item(r,key,k5):return OLD.item(r,key,k5)+(160*key[9]if r==4 else 0)
def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)

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
    rep=representative;types=TYPES
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    dom=fixed_domain();assert len(dom)==len({key for cs,u,key in dom})
    assert {(cs,u,key[:9])for cs,u,key in dom}==set(C.C.L.fixed_domain())
    kd=sorted({(key,k5)for cs,u,key in dom for k5 in range((N-u)//5+1)})
    keys={(key,k5,e)for key,k5 in kd for e in(False,True)}
    tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in saved['tables']}if saved else{k:table(*k)for k in sorted(keys)}
    assert set(tabs)==keys and(saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tabs[k],*k)for k in sorted(keys))
    residual=[rep(('A'if r<=8 else'D',r))for r in ITEMS]
    assert 5 not in ITEMS and all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];expanded=[];identities=choices=0
    for cs,used,key in dom:
        fixedlarge=fixed_components(cs,key)
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
                    old=36*j*(j-1)if oldguard else 0;new=64*j*k5 if newguard else 0
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
            row=dict(counts=list(cs),key=list(key),fixed_A12_count=key[9],fixed_A5_count=k5,forest=sorted(k+str(n)for k,n in fs),
                sixtyfourths_upper=upper,sixtyfourths_actual=value(fs),winning_option=win,context_options=options)
            expanded.append(row);options5.append(row)
        rows.append(max(options5,key=lambda z:(z['sixtyfourths_upper'],z['fixed_A5_count'])))
    maximum=max(r['sixtyfourths_upper']for r in rows)
    maximizers=[r for r in expanded if r['sixtyfourths_upper']==maximum]
    above=sorted((r for r in expanded if r['sixtyfourths_upper']>64*BETA),key=lambda r:(-r['sixtyfourths_upper'],r['forest']))
    assert len(expanded)==len(kd) and len(tabs)==2*len(kd)
    assert above and maximum>64*BETA
    return dict(certified_upper=F(maximum,64),curvature=str(BETA),upper_excess=F(maximum,64)-BETA,
        actual_maximum_claimed=False,actual_attainment_of_global_upper_claimed=False,
        representative_maximizers=maximizers,fixed_large_vectors=len(dom),
        fixed_large_and_A5_vectors=len(expanded),interaction_keys=len(kd),
        tables=[tabs[k]for k in sorted(keys)],simultaneous_type_pairs=len(types)**2,
        residual_item_checks=identities,Bellman_upper_inequalities=checks,attaining_predecessors=len(keys)*N,
        context_options=choices,alternatives=rows,expanded_alternatives=expanded,above_curvature=above,
        representative_additive_attainment_only=True,rank15_actual_costs=dict(A15='201/2',B15='97',D15='201/2'),
        rank15_representative_upper='201/2',fresh_residual_tables=True,
        contextual_witness_multiplicity_not_counted=True,residual_rank17_D17_fixed_rank17_A17=True,
        rank12_actual_costs={k+'12':F(actual_base((k,12)),64)for k in'ABD'},
        rank12_representative_costs={k+'12':F(upper_base((k,12)),64)for k in'AD'},
        actual_component_costs={name:F(actual_cross(*pair),64)for name,pair in{
            'A1_A11':(('A',1),('A',11)),'A2_A11':(('A',2),('A',11)),
            'A1_A12':(('A',1),('A',12)),'A1_A13':(('A',1),('A',13)),
            'A4_A13':(('A',4),('A',13)),'A4_A11':(('A',4),('A',11))}.items()},
        representative_component_costs={name:F(upper_cross(*pair),64)for name,pair in{
            'A1_D11':(('A',1),('D',11)),'A2_D11':(('A',2),('D',11)),
            'A5_A12':(('A',5),('A',12)),'A5_D12':(('A',5),('D',12)),
            'A1_D12':(('A',1),('D',12)),'A1_D13':(('A',1),('D',13)),
            'A4_D13':(('A',4),('D',13)),'A4_D11':(('A',4),('D',11))}.items()})


def actual_typed_refinement(ff,saved=None):
    """Exact actual-type maxima in every expanded upper class above377.

    Rank5 and large-component multiplicities are fixed. Actual A1 and A2
    multiplicities are explicit, so all three contextual guards are exact.
    The remaining actual types interact only with the fixed components.
    """
    engine=ENGINE;rep=representative
    types=sorted(t for t in engine.types()if t[1]<=N)
    fibers={}
    for t in types:fibers.setdefault(rep(t),[]).append(t)
    residual_reps={rep(('A'if r<=8 else'D',r))for r in ITEMS}
    residual=[t for t in types if rep(t)in residual_reps and t not in [('A',1),('A',2)]]
    assert all(actual_cross(a,b)==0 for a in residual+[('A',1),('A',2)]for b in residual+[('A',1),('A',2)])
    old={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count']):r for r in saved['classes']}if saved else {}
    expected={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count'])for r in ff['above_curvature']}
    if saved:assert set(old)==expected and len(old)==len(saved['classes'])
    classes=[];checks=options_count=attainers=0
    for upperrow in ff['above_curvature']:
        cs=tuple(upperrow['counts']);k5=upperrow['fixed_A5_count'];a12=upperrow['fixed_A12_count'];key=(cs,k5,a12)
        groups=sorted(Counter(fixed_components(cs,upperrow['key'])).items())
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
                    oldbonus=36*a2*(a2-1)if co['A',10]and a1>=2 else 0
                    newbonus=64*a2*co['D',5]if co['A',9]and a1>=1 else 0
                    witnessbonus=24*a1*co['D',5]if co['A',11]or(co['A',9]and a2>0)else 0
                    choices.append((fixedcost+a1*w1+a2*w2+vals[left]+oldbonus+newbonus+witnessbonus,a1,a2,left))
            options_count+=len(choices);best,a1,a2,left=max(choices)
            forest=list(fixed)+[('A',1)]*a1+[('A',2)]*a2
            while left:
                t=residual[last[left]];forest.append(t);left-=t[1]
            assert sum(n for k,n in forest)==N and value(forest)==best<=upperrow['sixtyfourths_upper']
            row=dict(fixed_actual=fixed,residual_types=residual,weights=weights,values=vals,last=last,
                sixtyfourths_actual_maximum=best,attaining_forest=sorted(forest),actual_A1_count=a1,actual_A2_count=a2)
            if saved:assert enc(row)==savedrows[fixed]
            rows.append(row)
        maximum=max(r['sixtyfourths_actual_maximum']for r in rows)
        classes.append(dict(counts=cs,fixed_A5_count=k5,fixed_A12_count=a12,residual_rank=rem,
            representative_upper=F(upperrow['sixtyfourths_upper'],64),actual_maximum=F(maximum,64),
            actual_fixed_type_assignments=len(rows),rows=rows,
            actual_obstruction_retained=maximum>64*BETA,
            maximizers=[r['attaining_forest']for r in rows if r['sixtyfourths_actual_maximum']==maximum]))
    maximum=max(F(r['actual_maximum'])for r in classes)
    winning_class=next(r for r in classes if F(r['actual_maximum'])==maximum)
    witness=winning_class['maximizers'][0]
    assert sum(n for k,n in witness)==N
    assert F(value(witness),64)==maximum<=F(ff['certified_upper'])
    assert maximum>377
    return dict(classes=classes,complete_expanded_upper_classes=len(classes),
        exact_actual_maximum=maximum,actual_attaining_forest=witness,
        actual_fixed_type_assignments=sum(r['actual_fixed_type_assignments']for r in classes),
        actual_residual_tables=sum(r['actual_fixed_type_assignments']for r in classes),
        Bellman_inequalities=checks,attaining_predecessors=attainers,explicit_A1_A2_options=options_count,
        classes_still_above=sum(r['actual_obstruction_retained']for r in classes),
        all_other_classes_at_most_curvature=True,all_three_context_guards_exact=True,
        actual_type_family_I2star_is_inherited_equivalence=True,
        every_obstructing_actual_forest_enumerated=False,
        scope='Complete exact actual-type maximum within every above-curvature expanded class. All actual large/rank5 types and residual types are covered; explicit A1/A2 counts evaluate all three contextual charges exactly. Lists of maximizing witnesses do not enumerate every lower-cost obstructing forest or assert geometric realization.')


def run(saved=None):
    assert all(h!='PENDING'for h in FINAL_REPORTS.values()),'Frozen local reports required'
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert sha256((ROOT/f'results/{s}.md').read_bytes()).hexdigest()==h,s
    assert docs[L.TAG]['new_same_component_costs']==dict(A16='219/2',D16='219/2')
    assert actual_base(('A',15))==actual_base(('D',15))==64*F(201,2)
    for a in TYPES:
        assert actual_base(a)==W.actual_base(a)-({('A',16):48,('D',16):16}.get(a,0))
        for b in TYPES:
            assert actual_cross(a,b)==W.actual_cross(a,b)+(160 if L.pair(a,b)else 0)
            # The representative ordinary charge includes the old unconditional
            # A1/D5 credit. It dominates the new guarded total for every forest.
            reinstatement=24 if W.ispair(a,b)else 0
            assert actual_cross(a,b)+reinstatement<=upper_cross(representative(a),representative(b))
    ff=finite(saved['finite']if saved else None)
    print('UPPER',ff['certified_upper'],'classes',len(ff['above_curvature']),flush=True)
    ref=actual_typed_refinement(ff,saved['actual_typed_refinement']if saved else None)
    ff['actual_maximum_claimed']=True
    ff['actual_attainment_of_global_upper_claimed']=ref['exact_actual_maximum']==F(ff['certified_upper'])
    ff['exact_actual_maximum']=ref['exact_actual_maximum'];ff['actual_attaining_forest']=ref['actual_attaining_forest']
    idealpath=ROOT/f'results/{C.C.L.TAG}.json'
    ii=C.C.unchanged_ideal(json.loads(idealpath.read_text()))
    paths={Path(__file__),idealpath}|{ROOT/f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in
        (('scripts','py'),('results','json'),('results','md'))}
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    excess=F(ref['exact_actual_maximum'])-BETA
    return enc(dict(status='PASS_CONDITIONAL_A16_D16_WITNESS_CONTEXT_AMBIENT58_FRONTIER',sources=sources,
        ambient_dimension=58,completed_unrestricted_bound=58,new_dimension_bound=None,scale=64,
        conditional_on_joint_local_whole_face_composition=True,whole_face_composition_proved_here=False,
        finite_curvature=BETA,strict_global_constant=4205,finite=ff,actual_typed_refinement=ref,ideal=ii,
        actual_costs={k+str(n):F(actual_base((k,n)),64)for k,n in(('A',15),('D',15),('A',16),('D',16))},
        actual_A4_A12_cost=F(actual_cross(('A',4),('A',12)),64),rank16_upper=F(upper_base(('D',16)),64),
        rank15_upper=F(upper_base(('D',15)),64),ordinary_A1_D5_cost='0',guarded_A1_D5_weight='3/8',
        witness_context='whole A11 OR(whole A9 AND whole A2); once per marked pair',
        representative_A1_rank5_allowance_retained=True,
        residual_item_delta='160*a12 at rank4 only; fixedrank16 changes evaluated directly; rank15 quarters unchanged',
        three_contexts_exact_in_actual_fibers=True,
        necessary_ambient58_inequality_if_composed=f'({excess})*v_f >= 4205 + 36*v_i',
        scope='Complete fresh changed-cost tables and actual-type refinements of every above-curvature class, including all three contextual rules. Positive low8 ideal certificate unchanged. Conditional on final joint proof composition; no dimension57 theorem or realization.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert a.verify or not p.exists()
    saved=json.loads(p.read_text())if a.verify else None;out=run(saved)
    if a.verify:assert out==saved
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_typed_refinement']['exact_actual_maximum'],out['actual_typed_refinement']['classes_still_above'],flush=True)
