#!/usr/bin/env python3
"""Fresh ambient59 tables with deep finite6 and both contextual coupons.

Final geometric bindings are required before authoritative generation. The
optional preflight performs only the new complete arithmetic calculation.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_order15_d16_d17_ambient59_frontier_20261018 as C
import fresh_depth6_finite6_zero_transfer_20261018 as S
import fresh_order13_bc_inner_terminal_transfer_20261018 as H
ROOT=C.ROOT;N=59;SCALE=16;ITEMS=C.ITEMS
TAG='fresh_depth6_finite6_order13_order17_ambient59_frontier_20261018'
CENTER='fresh_order17_center_contextual_zero_transfer_20261018'
FINAL_MDS={S.TAG:'cc2deec8ac4391dd9889b9d1591b25eb2074ce1ac904d0a99456d6232743f8a6',
 H.TAG:'4354857b1e4ce2b5b2a67e58d31ccd9cbe1cecf944e1f94502f00007483b6ea0',CENTER:'7d4e96b2d7be94925208016f00f16c3d8f00d3b2f2c918ea168cf32f5a2bf021'}
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,S.TAG,H.TAG,CENTER)))
enc=C.enc


def deleted(t):return max(0,t[1]-dict(A=18,B=19,D=20)[t[0]])if t[0]in('A','B','D')else 0

def actual_base(t):return C.actual_base(t)-16*deleted(t)
def upper_base(t):return C.upper_base(t)-16*deleted(t)

def actual_cross(a,b):
    p={a,b}
    return C.actual_cross(a,b)+(32 if p=={('A',1),('A',12)}else 8 if p=={('A',1),('B',12)}else 0)

def upper_cross(a,b):return C.upper_cross(a,b)

def context(fs,upper=False):
    co=Counter(fs);old=C.R.context_cost(fs,upper)
    witness=any(n==9 for k,n in fs)if upper else bool(co['A',9])
    new=16*co['A',2]*co['A'if upper else'D',5]if co['A',1]and witness else 0
    return old+new

def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)

def item(r,key,k5):return C.item(r,key,k5)-16*deleted(('A'if r<=8 else'D',r))


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
    rep=C.P.Z.B.B.representative;types=[t for t in C.P.Z.B.B.types()if t[1]<=N]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0 and deleted(a)==S.count(a)
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    dom=C.P.fixed_domain();assert len(dom)==961 and len({key for cs,u,key in dom})==960
    kd=sorted({(key,k5)for cs,used,key in dom for k5 in range((N-used)//5+1)})
    keys={(key,k5,e)for key,k5 in kd for e in(False,True)}
    tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in saved['tables']}if saved else{k:table(*k)for k in sorted(keys)}
    assert set(tabs)==keys and(saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tabs[k],*k)for k in sorted(keys));residual=[rep(('A'if r<=8 else'D',r))for r in ITEMS]
    assert 5 not in ITEMS and all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];expanded=[];identity_checks=choices=0
    for cs,used,key in dom:
        large=[('B'if r==9 else'A'if r==17 else'D',r)for r,c in zip(range(9,18),cs)for _ in range(c)]
        assert sum(n for k,n in large)==used and key[2]==cs[4]and key[4]==cs[5]and key[6]==cs[1]
        options5=[]
        for k5 in range((N-used)//5+1):
            fixed=large+[('A',5)]*k5;rem=N-used-5*k5
            for r,t in zip(ITEMS,residual):
                assert item(r,key,k5)==upper_base(t)+sum(upper_cross(t,b)for b in fixed);identity_checks+=1
            normal=tabs[key,k5,False];no2=tabs[key,k5,True]
            options=[dict(value=normal['values'][rem],fixed_A2_count=0,reserved_A1_count=0,exclude2=False,old_context=0,new_context=0)]
            # Every forest lies in the no-context class, old-context class
            # (>=2 A1 and rank10), or new-only class (>=1 A1 and rank9).
            # Reserving1 allows extra A1 in the table, a safe relaxation.
            # When both contexts hold, reserve2 includes both exact costs.
            for reserve in(1,2):
                oldguard=reserve==2 and bool(cs[1]);newguard=bool(cs[0]and k5)
                if rem<reserve or not(oldguard or newguard)or(reserve==2 and not oldguard):continue
                for j in range((rem-reserve)//2+1):
                    old=4*j*(j-1)if oldguard else 0;new=16*j*k5 if newguard else 0
                    cost=reserve*item(1,key,k5)+j*item(2,key,k5)+no2['values'][rem-reserve-2*j]+old+new
                    options.append(dict(value=cost,fixed_A2_count=j,reserved_A1_count=reserve,exclude2=True,old_context=old,new_context=new))
            choices+=len(options);winner=max(options,key=lambda z:(z['value'],z['reserved_A1_count'],z['fixed_A2_count']))
            j=winner['fixed_A2_count'];reserve=winner['reserved_A1_count'];tab=tabs[key,k5,winner['exclude2']]
            fs=fixed+[('A',2)]*j+[('A',1)]*reserve;rest=rem-2*j-reserve
            while rest:
                r=tab['last'][rest];fs.append(rep(('A'if r<=8 else'D',r)));rest-=r
            fixedcost=sum(upper_base(t)for t in fixed)+sum(upper_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            upper=fixedcost+winner['value']
            # The relaxed one-A1 branch cannot maximize with an unpaid old
            # contextual bonus: the two-A1 branch would have greater value.
            assert sum(n for k,n in fs)==N and value(fs,True)==upper and value(fs)<=upper
            row=dict(counts=list(cs),key=list(key),fixed_A5_count=k5,forest=sorted(k+str(n)for k,n in fs),
                sixteenths_upper=upper,sixteenths_actual=value(fs),winning_option=winner,context_options=options)
            expanded.append(row);options5.append(row)
        rows.append(max(options5,key=lambda z:(z['sixteenths_upper'],z['fixed_A5_count'])))
    maximum=max(r['sixteenths_upper']for r in rows);witness=[('A',14)]+[('A',15)]*3
    assert sum(n for k,n in witness)==N and value(witness)==maximum
    above=[r for r in rows if r['sixteenths_upper']>6264]
    return dict(exact_maximum=F(maximum,16),curvature='783/2',actual_attainer=['A14','A15','A15','A15'],
        fixed_large_vectors=len(dom),fixed_large_and_A5_vectors=len(expanded),interaction_keys=len(kd),
        tables=[tabs[k]for k in sorted(keys)],simultaneous_type_pairs=len(types)**2,residual_item_checks=identity_checks,
        Bellman_upper_inequalities=checks,attaining_predecessors=len(keys)*N,context_options=choices,
        alternatives=rows,expanded_alternatives=expanded,above_curvature=above,
        component_costs={name:F(actual_cross(*pair),16)for name,pair in{'A1_A12':(('A',1),('A',12)),'A1_B12':(('A',1),('B',12)),
            'A1_D12':(('A',1),('D',12)),'D5_D12':(('D',5),('D',12)),'D5_A12':(('D',5),('A',12))}.items()},
        contextual_old='A2_count*(A2_count-1)/4 if whole A10 and at least two whole A1; upper allows any rank10',
        contextual_new='A2_count*D5_count if whole A9 and at least one whole A1; upper allows any rank9 and substitutes A5 for D5',
        context_witness_multiplicities_not_counted=True,fixed_A5_counts_complete=True,residual_rank5_excluded=True,
        geometric_attainment_claimed=False,new_residual_tables=True)


def ideal(docs):
    old=docs[C.TAG]['ideal'];h={tuple(r['ranks']):r for r in docs[H.TAG]['ideal_frontier_update']['rows']}
    center=docs[CENTER]['ideal'];z={tuple(r['ranks']):r for r in center['rows']};rows=[]
    assert len(old['rows'])==len(h)==len(z)==1890
    env={int(k):F(x)for k,x in old['current_envelopes'].items()};assert env[12]==F(145,26)and env[16]==F(229,34)
    env[12]=F(143,26);env[16]=F(227,34)
    for r in old['rows']:
        rs=tuple(r['ranks']);saving13=F(rs.count(12),13);saving17=F(rs.count(16),17)
        assert F(r['normalized_margin'])+saving13==F(h[rs]['normalized_margin'])
        margin=F(r['normalized_margin'])+saving13+saving17;cost=sum(env[n]for n in rs)
        assert margin==F(r['curvature'])/r['link_vertices']-cost==F(z[rs]['normalized_margin'])
        assert r['link_vertices']==h[rs]['link_vertices']==z[rs]['link_vertices']
        rows.append(dict(ranks=rs,link_vertices=r['link_vertices'],G=r['G'],curvature=r['curvature'],normalized_cost=cost,
            normalized_margin=margin,absolute_margin=margin*r['link_vertices'],order13_saving=saving13,order17_saving=saving17))
    summary=C.U.summary(rows);assert summary['negative_count']==1 and summary['zero_count']==0
    assert summary['negative_rows'][0]['ranks']==(9,9,9,15,16)and summary['normalized_minimum']==-F(69,5440)
    return dict(complete_profiles=1890,current_envelopes=env,rows=rows,summary=summary,
        no_new_convolution=True,negative_absolute_margins_not_extended=True,
        conditional_order16_quarter=C.U.conditional(rows,{15:F(1,4)}),next_order16_threshold=C.U.threshold(rows,15))


def literal(docs):
    old={C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']};records={}
    for tag in(S.TAG,H.TAG,CENTER):
        rr=docs[tag]['current_literal_interface']['rows'];records[tag]={C.R.key(r['source'],r['case']):r for r in rr}
        assert len(rr)==len(records[tag])==216 and set(records[tag])==set(old)
    rows=[];matches=0
    for key,before in old.items():
        later=[records[tag][key]for tag in(S.TAG,H.TAG,CENTER)];corners=[]
        assert all(len(r['corners'])==len(before['corners'])for r in later)
        for i,c in enumerate(before['corners']):
            now=[r['corners'][i]for r in later];assert all(x['marks']==c['marks']for x in now)
            assert all(F(x['debit'])==F(x['increment'])==0 for x in now[:2])
            base=F(c['retained_payment_lower_bound']);assert all(F(x['retained_payment_lower_bound'])==base for x in now[:2])
            pay=base-F(now[2]['debit'])+F(now[2]['increment']);assert pay==F(now[2]['retained_payment_lower_bound'])>=0
            corners.append(dict(marks=c['marks'],previous_lower=base,finite6_delta='0',order13_delta='0',
                order17_debit=now[2]['debit'],order17_increment=now[2]['increment'],retained_payment_lower_bound=pay));matches+=3
        total=sum(c['retained_payment_lower_bound']for c in corners);assert total>=1
        rows.append(dict(source=before['source'],case=before['case'],corners=corners,total_lower_bound=total))
    return dict(configurations=216,local_corner_comparisons=matches,rows=rows,minimum_total=min(r['total_lower_bound']for r in rows),
        native_Grams_reconstructed=False,scope='Exact aligned saved literal interfaces, not a new historical corner census')


def run(saved=None):
    assert all(h!='PENDING'for h in FINAL_MDS.values()),'Wait for final geometric inputs.'
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    for tag,h in FINAL_MDS.items():assert sha256((ROOT/f'results/{tag}.md').read_bytes()).hexdigest()==h,tag
    assert docs[S.TAG]['feature_depth']==docs[S.TAG]['finite_distance']==6 and docs[S.TAG]['new_weight']=='0'
    assert docs[S.TAG]['new_coupon']['component_pair_increment']=='2'and docs[S.TAG]['all_ideal_weights_unchanged']
    assert docs[H.TAG]['new_order13_maximum']=='143/2'and docs[H.TAG]['new_coupon']['component_pair_increment']=='1/2'
    cp=docs[CENTER]['new_coupon'];assert cp['required_whole_unmarked_components']==['A1','A9']and cp['component_pair_increment']=='1'and cp['witnesses_not_multiplied']
    assert docs[CENTER]['local_proof']['new_cost']=='227/2'and docs[CENTER]['unchanged_D5_D12_center_weight']=='1/2'
    assert docs[S.TAG]['current_literal_interface']['changed']==docs[H.TAG]['current_literal_interface']['changed']==[]
    paths={Path(__file__),Path(C.__file__),Path(S.__file__),Path(H.__file__),Path(C.P.__file__),Path(C.R.__file__),Path(C.U.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_DEEP_FINITE6_ORDER13_ORDER17_AMBIENT59_FRONTIER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        finite=finite(saved['finite']if saved else None),ideal=ideal(docs),combined_literal_interface=literal(docs),
        ambient_dimension=59,completed_unrestricted_bound=59,new_dimension_bound=None,finite_curvature='783/2',strict_global_constant='4350',constant_counted_once=True,
        triangle_minimum=2,quadrilateral_minimum=1,geometric_inputs_final=True,
        scope='Complete fresh rank59 finite certificates and saved-domain ideal update under the simultaneous proved rule. No new ambient bound or arithmetic optimization.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--preflight',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if args.preflight:
        out=enc(finite());Path('/private/tmp/'+TAG+'_finite_preflight.json').write_text(json.dumps(out,indent=2)+'\n')
        print('CONDITIONAL_ARITHMETIC_PREFLIGHT',out['exact_maximum'],len(out['above_curvature']),flush=True)
    else:
        saved=json.loads(path.read_text())if args.verify else None;out=run(saved)
        if args.verify:assert out==saved
        else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
        print(out['status'],out['finite']['exact_maximum'],len(out['finite']['above_curvature']),out['ideal']['summary']['negative_count'],flush=True)
