#!/usr/bin/env python3
"""Conditional complete ambient58 forest count with a separate D12 type.

This certifies arithmetic for the stated proposed payments. The separate
whole-face composition theorem is not assumed or claimed here.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement,product
from pathlib import Path
import fresh_order10_low8_ambient58_frontier_20261018 as L
ROOT=L.ROOT;N=58;SCALE=32;BETA=F(377);ITEMS=L.ITEMS;enc=L.enc
TAG='fresh_d12_type_split_ambient58_forest_frontier_20261018'
STEMS=(L.TAG,'fresh_d12_radical_saturated_half_replacements_20261018',
 'fresh_d12_remaining_selected_half_replacements_20261018',
 'fresh_d12_low_protected_adaptive_interfaces_20261018',
 'fresh_d12_universal_ideal_output_repairs_20261018')
ENGINE=L.C.C.C.P.Z.B.B
TYPES=tuple(sorted(t for t in ENGINE.types()if t[1]<=N))

def representative(t):
    return ('A',12)if t in(('A',12),('B',12))else ENGINE.representative(t)
FIBERS={}
for _t in TYPES:FIBERS.setdefault(representative(_t),[]).append(_t)

@lru_cache(None)
def actual_base(t):return L.actual_base(t)-(32 if t==('D',12)else 0)
@lru_cache(None)
def actual_cross(a,b):
    pair=frozenset((a,b))
    extra=32 if pair==frozenset((('A',1),('A',11)))else 64 if pair==frozenset((('A',2),('A',11)))else 0
    return L.actual_cross(a,b)+extra

def upper_base(t):
    if t==('D',12):return 32*65
    if t==('A',12):return max(actual_base(a)for a in FIBERS[t])
    return L.upper_base(t)

@lru_cache(None)
def upper_cross(a,b):
    # The new A12 representative covers both actual A12/B12. Its partners
    # must cover their full fibers too: in particular representative A5
    # includes actual D5, whose A12 interaction is nonzero.
    if ('A',12)in(a,b):return max(actual_cross(x,y)for x in FIBERS[a]for y in FIBERS[b])
    return L.upper_cross(a,b)

def context(fs,upper=False):return L.context(fs,upper)
def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)

def fixed_domain():
    # Split every old rank12 count into A/B12 and D12 counts, including
    # every mixture. The other fixed ranks and their type fibers are unchanged.
    return [(cs,used,key+(a12,))for cs,used,key in L.fixed_domain()for a12 in range(cs[3]+1)]

def fixed_components(cs,key):
    out=[]
    for r,c in zip(range(9,18),cs):
        if r==12:out += [('A',12)]*key[9]+[('D',12)]*(c-key[9])
        else:out += [('B'if r==9 else'A'if r==17 else'D',r)]*c
    return out

def item(r,key,k5):
    t=representative(('A'if r<=8 else'D',r))
    return L.item(r,key[:9],k5)+key[9]*(upper_cross(t,('A',12))-upper_cross(t,('D',12)))

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
    assert {(cs,u,key[:9])for cs,u,key in dom}==set(L.fixed_domain())
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
                    old=8*j*(j-1)if oldguard else 0;new=32*j*k5 if newguard else 0
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
                thirtyseconds_upper=upper,thirtyseconds_actual=value(fs),winning_option=win,context_options=options)
            expanded.append(row);options5.append(row)
        rows.append(max(options5,key=lambda z:(z['thirtyseconds_upper'],z['fixed_A5_count'])))
    maximum=max(r['thirtyseconds_upper']for r in rows)
    maximizers=[r for r in expanded if r['thirtyseconds_upper']==maximum]
    above=sorted((r for r in expanded if r['thirtyseconds_upper']>32*BETA),key=lambda r:(-r['thirtyseconds_upper'],r['forest']))
    assert len(expanded)==len(kd) and len(tabs)==2*len(kd)
    assert above and maximum>32*BETA
    return dict(certified_upper=F(maximum,32),curvature=str(BETA),upper_excess=F(maximum,32)-BETA,
        actual_maximum_claimed=False,actual_attainment_of_global_upper_claimed=False,
        representative_maximizers=maximizers,fixed_large_vectors=len(dom),
        fixed_large_and_A5_vectors=len(expanded),interaction_keys=len(kd),
        tables=[tabs[k]for k in sorted(keys)],simultaneous_type_pairs=len(types)**2,
        residual_item_checks=identities,Bellman_upper_inequalities=checks,attaining_predecessors=len(keys)*N,
        context_options=choices,alternatives=rows,expanded_alternatives=expanded,above_curvature=above,
        representative_additive_attainment_only=True,rank15_actual_costs=dict(A15='100',B15='97',D15='100'),
        rank15_representative_upper='100',fresh_residual_tables=True,
        contextual_witness_multiplicity_not_counted=True,residual_rank17_D17_fixed_rank17_A17=True,
        rank12_actual_costs={k+'12':F(actual_base((k,12)),32)for k in'ABD'},
        rank12_representative_costs={k+'12':F(upper_base((k,12)),32)for k in'AD'},
        actual_component_costs={name:F(actual_cross(*pair),32)for name,pair in{
            'A1_A11':(('A',1),('A',11)),'A2_A11':(('A',2),('A',11)),
            'A1_A12':(('A',1),('A',12)),'A1_A13':(('A',1),('A',13)),
            'A4_A13':(('A',4),('A',13)),'A4_A11':(('A',4),('A',11))}.items()},
        representative_component_costs={name:F(upper_cross(*pair),32)for name,pair in{
            'A1_D11':(('A',1),('D',11)),'A2_D11':(('A',2),('D',11)),
            'A5_A12':(('A',5),('A',12)),'A5_D12':(('A',5),('D',12)),
            'A1_D12':(('A',1),('D',12)),'A1_D13':(('A',1),('D',13)),
            'A4_D13':(('A',4),('D',13)),'A4_D11':(('A',4),('D',11))}.items()})


def actual_typed_refinement(ff,saved=None):
    """Exact actual-type maxima in every expanded upper class above377.

    Rank5 and large-component multiplicities are fixed. Actual A1 and A2
    multiplicities are explicit, so both contextual guards are exact.
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
                    oldbonus=8*a2*(a2-1)if co['A',10]and a1>=2 else 0
                    newbonus=32*a2*co['D',5]if co['A',9]and a1>=1 else 0
                    choices.append((fixedcost+a1*w1+a2*w2+vals[left]+oldbonus+newbonus,a1,a2,left))
            options_count+=len(choices);best,a1,a2,left=max(choices)
            forest=list(fixed)+[('A',1)]*a1+[('A',2)]*a2
            while left:
                t=residual[last[left]];forest.append(t);left-=t[1]
            assert sum(n for k,n in forest)==N and value(forest)==best<=upperrow['thirtyseconds_upper']
            row=dict(fixed_actual=fixed,residual_types=residual,weights=weights,values=vals,last=last,
                thirtyseconds_actual_maximum=best,attaining_forest=sorted(forest),actual_A1_count=a1,actual_A2_count=a2)
            if saved:assert enc(row)==savedrows[fixed]
            rows.append(row)
        maximum=max(r['thirtyseconds_actual_maximum']for r in rows)
        classes.append(dict(counts=cs,fixed_A5_count=k5,fixed_A12_count=a12,residual_rank=rem,
            representative_upper=F(upperrow['thirtyseconds_upper'],32),actual_maximum=F(maximum,32),
            actual_fixed_type_assignments=len(rows),rows=rows,
            actual_obstruction_retained=maximum>32*BETA,
            maximizers=[r['attaining_forest']for r in rows if r['thirtyseconds_actual_maximum']==maximum]))
    maximum=max(F(r['actual_maximum'])for r in classes)
    witness=[('A',16)]+[('D',14)]*3
    assert sum(n for k,n in witness)==N
    assert F(value(witness),32)==maximum==F(ff['certified_upper'])==F(1527,4)
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


def unchanged_ideal(prior):
    ii=prior['ideal'];rows=ii['rows'];D=57
    domain=sorted((lo,)+p if lo else p for lo in range(8)for p in L.high8_parts(D-lo))
    assert len(domain)==len(set(domain))==len(rows)==2552
    # Independently verify the saved polynomial, cost and odd-dimensional
    # curvature identities, without reconvolving an unchanged ideal rule.
    env={int(k):F(v)for k,v in ii['current_envelopes'].items()}
    for rs,row in zip(domain,rows):
        assert list(rs)==row['ranks'] and sum(rs)==D
        coeff=row['h_polynomial'];V=1
        for r in rs:V*=r+1
        assert len(coeff)==D+1 and coeff==coeff[::-1]and sum(coeff)==V
        assert sum(i*v for i,v in enumerate(coeff))==F(D*V,2)
        cost=sum(env[r]for r in rs)
        absolute=sum(abs(2*i-D)*v for i,v in enumerate(coeff))
        curvature=F(3*D,4)-F(5*absolute,4*V)
        assert F(row['normalized_cost'])==cost
        assert row['link_vertices']==V and F(row['curvature'])==V*curvature
        assert row['G']==sum(min(i,D-i)*v for i,v in enumerate(coeff))
        assert F(row['normalized_margin'])==curvature-cost>0
        assert F(row['absolute_margin'])==V*(curvature-cost)
    assert min(F(r['normalized_margin'])for r in rows)==F(87,10240)
    assert min(F(r['absolute_margin'])for r in rows)==36
    return dict(complete_profiles=2552,all_positive=True,normalized_minimum='87/10240',absolute_surplus='36',
        weights_and_envelopes_unchanged=True,saved_polynomial_identities_rechecked=True,
        inherited_positive_low8_concentration_and_factor_fiber=True,no_convolution_repeated=True)

def run(saved=None):
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    prior=docs[L.TAG]
    assert prior['finite']['certified_upper']=='3063/8'and prior['strict_global_constant']==4205
    assert actual_base(('A',12))==actual_base(('B',12))==32*66 and actual_base(('D',12))==32*65
    assert actual_cross(('A',1),('A',11))==32*2 and actual_cross(('A',2),('A',11))==32*3
    ff=finite(saved['finite']if saved else None)
    print('REPRESENTATIVE',ff['certified_upper'],'above',len(ff['above_curvature']),flush=True)
    refined=actual_typed_refinement(ff,saved['actual_typed_refinement']if saved else None)
    ff['actual_maximum_claimed']=ff['actual_attainment_of_global_upper_claimed']=True
    ff['exact_actual_maximum']=refined['exact_actual_maximum'];ff['actual_attaining_forest']=refined['actual_attaining_forest']
    ii=unchanged_ideal(prior)
    paths={Path(__file__),Path(L.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return enc(dict(status='PASS_CONDITIONAL_D12_TYPE_SPLIT_AMBIENT58_FOREST_FRONTIER',sources=sources,
        ambient_dimension=N,completed_unrestricted_bound=58,new_dimension_bound=None,
        conditional_on_complete_D12_whole_face_composition=True,universal_D12_theorem_claimed=False,
        finite_curvature=BETA,strict_global_constant=4205,scale=SCALE,finite=ff,
        actual_typed_refinement=refined,ideal=ii,
        ordinary_coupon_increments={'whole_A1_A11_endpoints':'1','whole_A2_A11_endpoints':'2'},
        ordinary_representative_allowances_unchanged=True,
        contextual_guards_unchanged=['whole A10 and at least two A1, A2/A2','whole A9 and at least one A1, A2/D5'],
        rank12_type_split='A12/B12 -> A12 with base66; D12 -> D12 with base65',
        necessary_ambient58_inequality_if_composed='(19/4)*v_f >= 4205 + 36*v_i',
        scope='Complete conditional changed-rule finite count, with actual-type maxima in every above-curvature upper class. Frozen ideal rule and positive surplus unchanged. No universal D12 theorem, dimension57 bound or geometric realization follows.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    if not args.verify:assert not path.exists()
    out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['certified_upper'],out['actual_typed_refinement']['classes_still_above'],flush=True)
