#!/usr/bin/env python3
"""Simultaneous order10 local payments, fresh rank58 finite tables and low8 cusps.

The changed rank9 cost is never passed through the old low9 concentration.
This is a complete new positive-cusp certificate, not dimension57 closure.
"""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import prod
from pathlib import Path
import fresh_order10_d_opposite_fork_fractional_transfer_20261018 as L
import fresh_depth6_finite8_ambient58_frontier_20261018 as BASE
import fresh_ambient58_payment_frontier_20261018 as P
ROOT=BASE.ROOT;TAG='fresh_order10_low8_ambient58_frontier_20261018'
N=58;D=57;SCALE=32;BETA=F(377);C=BASE.C;ITEMS=BASE.ITEMS;enc=BASE.enc
FINAL_LOCAL_MD='3f2154ec10e3782918c36f62157986e1831983dbbfb36efeec3a0ff6069564cb'
STEMS=(L.TAG,L.A10.TAG,L.B10.TAG,L.B10.CZERO,BASE.TAG,'fresh_a15_d15_dimension58_20261018')
actual_base=L.actual_base;upper_base=L.upper_base
actual_cross=lru_cache(None)(L.actual_cross);upper_cross=L.upper_cross
fixed_domain=BASE.fixed_domain

def context(fs,upper=False):return 2*BASE.context(fs,upper)
def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)
def item(r,key,k5):return 2*BASE.item(r,key,k5)+(16 if r==1 else 32 if r==3 else 20 if r in(6,7,8)else 0)*key[5]


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
            row=dict(counts=list(cs),key=list(key),fixed_A5_count=k5,forest=sorted(k+str(n)for k,n in fs),
                thirtyseconds_upper=upper,thirtyseconds_actual=value(fs),winning_option=win,context_options=options)
            expanded.append(row);options5.append(row)
        rows.append(max(options5,key=lambda z:(z['thirtyseconds_upper'],z['fixed_A5_count'])))
    maximum=max(r['thirtyseconds_upper']for r in rows)
    maximizers=[r for r in expanded if r['thirtyseconds_upper']==maximum]
    above=sorted((r for r in expanded if r['thirtyseconds_upper']>32*BETA),key=lambda r:(-r['thirtyseconds_upper'],r['forest']))
    assert len(expanded)==len(kd)==2349 and len(tabs)==4698
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
        actual_component_costs={name:F(actual_cross(*pair),32)for name,pair in{
            'A1_A12':(('A',1),('A',12)),'A1_A13':(('A',1),('A',13)),
            'A4_A13':(('A',4),('A',13)),'A4_A11':(('A',4),('A',11))}.items()},
        representative_component_costs={name:F(upper_cross(*pair),32)for name,pair in{
            'A1_D12':(('A',1),('D',12)),'A1_D13':(('A',1),('D',13)),
            'A4_D13':(('A',4),('D',13)),'A4_D11':(('A',4),('D',11))}.items()})


@lru_cache(None)
def high8_parts(n,lo=8):
    if n==0:return((),)
    return tuple((r,)+p for r in range(lo,n+1)for p in high8_parts(n-r,r))


def concentration():
    rows=[]
    for a in range(1,8):
        for b in range(a,8):
            def probabilities(x,y):
                return [F(sum(i+j==k for i in range(x+1)for j in range(y+1)),(x+1)*(y+1))for k in range(x+y+1)]
            old=probabilities(a,b);new=probabilities(a-1,b+1)
            assert sum(old)==sum(new)==1 and sum(i*(y-x)for i,(x,y)in enumerate(zip(old,new)))==0
            gaps=[sum(F(max(i-t,0))*(y-x)for i,(x,y)in enumerate(zip(old,new)))for t in range(a+b+1)]
            assert all(z>=0 for z in gaps)
            assert F(a,2)+F(b,2)==F(a-1,2)+F(b+1,2)
            assert(a+1)*(b+1)-a*(b+2)==b-a+1>0
            rows.append(dict(a=a,b=b,stop_loss_gaps=list(map(str,gaps)),old_product=(a+1)*(b+1),new_product=a*(b+2)))
    assert len(rows)==28
    return dict(complete_low_rank_transfers=rows,cutoff=8,changed_rank9_never_transferred=True,
        terminal_form='arbitrary ranks>=8 and at most one rank1..7',cost_preserved=True,
        normalized_curvature_nondecreasing_in_reverse=True,vertex_product_nondecreasing_in_reverse=True,
        absolute_surplus_extension_requires_positive_terminal_margin=True)


def ideal(previous,local):
    env={int(k):F(v)for k,v in previous['current_envelopes'].items()}
    assert all(env[r]==F(r,2)for r in range(1,10))
    families={k:F(v)for k,v in local['current_order10_family_costs'].items()}
    assert families==dict(A=F(85,2),B=F(44),C=F(44),D=F(355,8))
    env[9]=max(families.values())/10;assert env[9]==F(71,16)
    domain=sorted((lo,)+p if lo else p for lo in range(8)for p in high8_parts(D-lo))
    other=[]
    def counts(r,left,out):
        if r>D:
            if left<8:other.append(tuple(([left]if left else[])+out))
            return
        for k in range(left//r+1):counts(r+1,left-r*k,out+[r]*k)
    counts(8,D,[])
    assert domain==sorted(other)and len(domain)==len(set(domain))==2552
    rows=[P.scalar(rs,env)for rs in domain]
    assert all(F(row['normalized_margin'])>0 for row in rows)
    summary=P.Z.summarize(rows)
    assert summary['negative_count']==summary['zero_count']==0 and summary['normalized_minimum']=='87/10240'and summary['absolute_minimum']=='36'
    assert all(env[r]==max(F(P.Z.family_cost(k,r+1),r+1)for k in 'ABCD')for r in range(17,D+1))
    return dict(current_envelopes=env,current_order10_family_costs=families,complete_profiles=len(rows),rows=rows,summary=summary,
        complete_low8_domain=True,independent_domain_reconstructions=2,fresh_convolutions=True,
        old_low9_domain_not_reused_as_complete=True,concentration=concentration(),absolute_surplus='36',
        actual_charge_bounded_above_by_envelope=True,absolute_surplus_valid_for_every_actual_cusp=True)


def run(saved=None):
    assert FINAL_LOCAL_MD!='PENDING'
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,doc in docs.items():
        assert doc['status'].startswith('PASS')
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert sha256((ROOT/f'results/{L.TAG}.md').read_bytes()).hexdigest()==FINAL_LOCAL_MD
    assert docs[BASE.TAG]['finite']['certified_upper']=='3063/8'and docs[BASE.TAG]['finite_curvature']=='377'and docs[BASE.TAG]['strict_global_constant']==4205
    assert docs['fresh_a15_d15_dimension58_20261018']['dimension_bound']==58
    ff=finite(saved['finite']if saved else None)
    assert ff['certified_upper']==F(3063,8)and len(ff['above_curvature'])==32
    witness=[('A',1)]*7+[('D',5)]*3+[('D',12)]*3
    assert sum(n for k,n in witness)==N and F(value(witness),32)==ff['certified_upper']
    ff.update(actual_maximum_claimed=True,actual_attainment_of_global_upper_claimed=True,exact_actual_maximum=F(3063,8),actual_attaining_forest=witness,
        all_obstructing_actual_forests_enumerated=False,prior28_actual_class_count_not_reasserted=True)
    added_pairs=[(('A',r),('A',9))for r in range(1,5)]+[(('A',1),('B',9)),(('A',3),('B',9))]+[(('A',r),('D',9))for r in range(1,9)]
    increments=[dict(types=[a,b],old=F(BASE.actual_cross(a,b),16),new=F(actual_cross(a,b),32),increment=F(actual_cross(a,b),32)-F(BASE.actual_cross(a,b),16))for a,b in added_pairs]
    assert len(increments)==14 and all(row['increment']>0 for row in increments)
    ii=ideal(docs[BASE.TAG]['ideal'],docs[L.TAG])
    g0=sum(min(j,D-j)for j in range(D+1));beta=F(5*g0-D*(D+1),2);correction=sum(2*j-D for j in range(D+1)if 2*j>D)
    assert(g0,beta,correction)==(812,BETA,841)and 5*correction==4205
    assert ff['certified_upper']-beta==F(47,8)
    paths={Path(__file__),Path(L.__file__),Path(BASE.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/f'results/{tag}.md'for tag in('fresh_distance14_unrestricted_dimension156_20261017','fresh_cusp_factor_payment_framework_20261017','fresh_cusp_link_shape_curvature_20261017'))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return enc(dict(status='PASS_ORDER10_LOW8_AMBIENT58_FRONTIER',sources=sources,ambient_dimension=N,link_dimension=D,
        completed_unrestricted_bound=58,new_dimension_bound=None,all_ideal_consequence_ceiling=57,novelty_or_best_known_claim=False,
        finite_curvature=beta,strict_global_constant=5*correction,scale=32,finite=ff,ideal=ii,
        ordinary_coupon_changes=increments,A17_endpoint_restoration='5/32',
        both_prior_contextual_coupons_unchanged=True,finite_weight_changes_nondecreasing=True,
        current_face_priority='Independent cycle10, C10, B10, then D10 actual-corner proofs; if none occurs inherit the finite8/depth6 baseline pointwise. All new finite changes increase weights.',
        necessary_ambient58_inequality='(47/8)*v_f >= 4205 + 36*v_i',constant_counted_once=True,
        scope='Complete simultaneous finite budget and positive low8 cusp envelope. Finite obstruction remains; no unrestricted dimension57 theorem or new example.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    saved=json.loads(path.read_text())if args.verify else None
    if not args.verify:assert not path.exists()
    out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['certified_upper'],out['ideal']['complete_profiles'],out['ideal']['summary'],flush=True)
