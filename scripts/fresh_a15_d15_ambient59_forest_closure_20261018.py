#!/usr/bin/env python3
"""Complete ambient59 finite upper with simultaneous A15/D15 zero halves.

Fresh additive upper/attainment certificates and exhaustive contextual
classes. A representative upper need not be an actual formal maximum.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_a14_cd16_ambient59_frontier_20261018 as C
import fresh_finite_a15_three_coupon_zero_transfer_20261018 as A15
import fresh_finite_d15_a4_a11_half_zero_transfer_20261018 as D15

ROOT=C.ROOT;N=59;ITEMS=C.ITEMS;enc=C.enc
TAG='fresh_a15_d15_ambient59_forest_closure_20261018'
FINAL_MDS={C.TAG:'257219a7c1e324f44d4c7ed1c39d846a0994da807d2328f07b2387bad9ddb6cf',
 A15.TAG:'90b45e20ba0835656a32481f8d42eb47bf2e96d82c38da1c986e1f0fd30efaa2',
 D15.TAG:'5fc4e6622636fbce7dcc3f349028b657d2f868d2593ebfadaf51488c97adfe61'}
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,*A15.STEMS,A15.TAG,*D15.STEMS,D15.TAG)))


def actual_base(t):return C.actual_base(t)-(16 if t==('A',15)else 8 if t==('D',15)else 0)
def upper_base(t):return C.upper_base(t)-(16 if t==('D',15)else 0)
def actual_cross(a,b):return A15.actual_cross(a,b)+(32 if {a,b}=={('A',4),('A',11)}else 0)
def upper_cross(a,b):return A15.upper_cross(a,b)
def context(fs,upper=False):return C.context(fs,upper)
def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)
def item(r,key,k5):return C.item(r,key,k5)+(key[2]if r==1 else 16*key[2]if r==4 else 0)


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
    rep=C.C.P.Z.B.B.representative;types=[t for t in C.C.P.Z.B.B.types()if t[1]<=N]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    dom=C.fixed_domain();assert len(dom)==len({key for cs,u,key in dom})==961
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
    maximum=max(r['sixteenths_upper']for r in rows);assert maximum==6264
    maximizers=[r for r in expanded if r['sixteenths_upper']==maximum]
    above=[r for r in expanded if r['sixteenths_upper']>6264];assert not above
    assert len(expanded)==len(kd)==2574 and len(tabs)==5148
    return dict(certified_upper=F(maximum,16),curvature='783/2',upper_slack='0',
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


def literal(docs):
    key=C.C.R.key;old={key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    changes={s:{key(r['source'],r['case']):r for r in docs[s]['current_literal_interface']['rows']}for s in(A15.TAG,D15.TAG)}
    assert len(old)==216 and all(set(d)==set(old)and len(d)==216 for d in changes.values())
    rows=[];debits=credits=0
    for k,r in old.items():
        aa,dd=changes[A15.TAG][k],changes[D15.TAG][k];corners=[]
        assert len(r['corners'])==len(aa['corners'])==len(dd['corners'])
        for oldc,a,d in zip(r['corners'],aa['corners'],dd['corners']):
            assert oldc['marks']==a['marks']==d['marks'];before=F(oldc['retained_payment_lower_bound'])
            ad,ac,ddelta,dc=map(F,(a['debit'],a['increment'],d['debit'],d['increment']))
            assert before-ad+ac==F(a['retained_payment_lower_bound'])
            assert before-ddelta+dc==F(d['retained_payment_lower_bound'])
            after=before-ad-ddelta+ac+dc;assert after>=0
            debits+=bool(ad or ddelta);credits+=bool(ac or dc)
            corners.append(dict(marks=oldc['marks'],previous_lower=before,A15_debit=ad,D15_debit=ddelta,
                A15_coupon_increment=ac,D15_coupon_increment=dc,retained_payment_lower_bound=after))
        total=sum(c['retained_payment_lower_bound']for c in corners);assert total>=1
        rows.append(dict(source=r['source'],case=r['case'],corners=corners,total_lower_bound=total))
    return dict(configurations=216,changed_debit_corners=debits,changed_credit_corners=credits,
        rows=rows,minimum_total=min(r['total_lower_bound']for r in rows),native_Grams_reconstructed=False,
        exact_local_native_audits_compared=True)


def run(saved=None):
    assert all(h!='PENDING'for h in FINAL_MDS.values()),'Wait for frozen local proofs.'
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    for s,h in FINAL_MDS.items():assert sha256((ROOT/f'results/{s}.md').read_bytes()).hexdigest()==h,s
    aa,dd=docs[A15.TAG],docs[D15.TAG]
    assert aa['new_A15_cost']==dd['new_cost']=='100'
    assert aa['retained_A15_unit']==[1,13]and aa['retained_A16_halves']==[[1,13],[2,14]]
    assert dd['retained_D16_half']==[1,12]and dd['new_coupon']['component_pair_increment']=='2'
    assert aa['new_coupons']==A15.COUPONS
    assert aa['all_ideal_weights_unchanged']and dd['ideal_weights_unchanged']
    ff=finite(saved['finite']if saved else None)
    # Recheck all frozen ideal polynomials and positive margins; no new
    # convolution is needed because every ideal payment is unchanged.
    ii=C.ideal(docs);assert enc(ii)==docs[C.TAG]['ideal']
    assert ii['summary']['negative_count']==ii['summary']['zero_count']==0
    assert ii['summary']['absolute_minimum']==F(91,2)
    ll=literal(docs)
    assert F(5*((N-1)**2//4)-N*(N-1),2)==F(783,2)
    assert 5*(N*N//4)==int(docs[C.TAG]['strict_global_constant'])==4350
    paths={Path(__file__),Path(C.__file__),Path(A15.__file__),Path(D15.__file__),Path(C.C.__file__),Path(C.C.U.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_A15_D15_AMBIENT59_FINITE_UPPER_783_OVER2',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        finite=ff,ideal=ii,combined_literal_interface=ll,ambient_dimension=59,
        finite_curvature='783/2',strict_global_constant='4350',constant_counted_once=True,
        triangle_minimum=2,quadrilateral_minimum=1,ideal_actual_surplus='91/2',
        necessary_ambient59_inequality='0 >= 4350 + (91/2)*v_i',
        new_dimension_bound=None,separate_global_composition_required=True,
        scope='Complete simultaneous finite upper and unchanged positive cusp envelope; local theorems and full current priorities are bound inputs. No actual attainment of representative global upper is claimed; the dimension statement is supplied separately.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    if not args.verify:assert not path.exists()
    out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['certified_upper'],len(out['finite']['above_curvature']),out['ideal']['summary']['negative_count'],flush=True)
