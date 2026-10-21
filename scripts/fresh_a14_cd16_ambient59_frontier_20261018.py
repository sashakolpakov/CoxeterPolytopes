#!/usr/bin/env python3
"""Complete ambient59 certificate with A14 zeros and the C/D16 sixteenth.

Fresh c15-resolved additive certificates and full contextual coverage;
actual/representative bases, local priorities and positive ideal extension.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path
import fresh_depth6_finite6_order13_order17_ambient59_frontier_20261018 as M
import fresh_order16_cd_sparse_sixteenth_transfer_20261018 as I
import fresh_finite_a14_near_endpoint_zero_transfer_20261018 as J
ROOT=M.ROOT;C=M.C;S=M.S;N=59;SCALE=16;ITEMS=M.ITEMS;enc=M.enc
TAG='fresh_a14_cd16_ambient59_frontier_20261018'
FINAL_MDS={M.TAG:'0dab0a9fbbbc9cf9eed55e9006a444bd723cd626d719d60a8f41eb90122451ea',
 I.TAG:'2b690d6a7a4393af171276168f469d39771d5f13b6b4a8fa152869e331fc6675',
 J.TAG:'9a3315bd2a960ae6202d8f2ad9b7c7b4e9d0d47ee94c971a25350d134cf85c58'}
STEMS=tuple(dict.fromkeys((*M.STEMS,M.TAG,*I.STEMS,I.TAG,*J.STEMS,J.TAG)))

def deleted(t):return M.deleted(t)
def actual_base(t):return M.actual_base(t)-(32 if t==('A',14) else 0)
def upper_base(t):return M.upper_base(t)-(8 if t==('D',14) else 0)
def actual_cross(a,b):return M.actual_cross(a,b)+int({a,b}=={('A',1),('D',15)})
def upper_cross(a,b):return M.upper_cross(a,b)+int({a,b}=={('A',1),('D',15)})
def context(fs,upper=False):return M.context(fs,upper)
def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)
def fixed_domain():
    dom=[(cs,u,key+(cs[6],))for cs,u,key in C.P.fixed_domain()]
    assert len(dom)==len({cs for cs,u,key in dom})==961
    assert all(key[8]==cs[6] for cs,u,key in dom)
    return dom

def item(r,key,k5):return M.item(r,key[:8],k5)+(key[8]if r==1 else 0)

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
    dom=fixed_domain();assert len(dom)==961 and len({key for cs,u,key in dom})==961
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
    maximum=max(r['sixteenths_upper']for r in rows);witness=[('D',14)]+[('A',15)]*3
    assert sum(n for k,n in witness)==N and value(witness)==maximum==6296
    above=[r for r in rows if r['sixteenths_upper']>6264]
    assert sorted(r['sixteenths_upper'] for r in above)==[6272,6276,6296]
    return dict(exact_maximum=F(maximum,16),curvature='783/2',actual_attainer=['D14','A15','A15','A15'],
        fixed_large_vectors=len(dom),fixed_large_and_A5_vectors=len(expanded),interaction_keys=len(kd),
        tables=[tabs[k]for k in sorted(keys)],simultaneous_type_pairs=len(types)**2,residual_item_checks=identity_checks,
        Bellman_upper_inequalities=checks,attaining_predecessors=len(keys)*N,context_options=choices,
        alternatives=rows,expanded_alternatives=expanded,above_curvature=above,
        component_costs={name:F(actual_cross(*pair),16)for name,pair in{'A1_A12':(('A',1),('A',12)),'A1_B12':(('A',1),('B',12)),
            'A1_D12':(('A',1),('D',12)),'D5_D12':(('D',5),('D',12)),'D5_A12':(('D',5),('A',12)),'A1_D15':(('A',1),('D',15))}.items()},
        contextual_old='A2_count*(A2_count-1)/4 if whole A10 and at least two whole A1; upper allows any rank10',
        contextual_new='A2_count*D5_count if whole A9 and at least one whole A1; upper allows any rank9 and substitutes A5 for D5',
        context_witness_multiplicities_not_counted=True,fixed_A5_counts_complete=True,residual_rank5_excluded=True,
        geometric_attainment_claimed=False,new_residual_tables=True)


def ideal(docs):
    old=docs[M.TAG]['ideal'];prior={tuple(r['ranks']):r for r in old['rows']}
    local={tuple(r['ranks']):r for r in docs[I.TAG]['ideal_frontier_update']['rows']}
    original={tuple(r['ranks']):r for r in docs[C.U.FRONT]['ideal']['rows']}
    domain=sorted(((low,)+p if low else p)for low in range(9)for p in C.U.high_parts(58-low))
    assert len(domain)==len(set(domain))==len(prior)==len(local)==len(original)==1890
    assert set(domain)==set(prior)==set(local)==set(original)
    env={int(k):F(v)for k,v in old['current_envelopes'].items()};assert env[15]==F(103,16)
    env[15]=F(411,64);rows=[]
    for rs in domain:
        r=prior[rs];h=original[rs]['h_polynomial'];v=prod(x+1 for x in rs)
        assert sum(rs)==58 and len(h)==59 and h==h[::-1]and all(isinstance(x,int)and x>0 for x in h)
        assert sum(h)==v==r['link_vertices']==local[rs]['link_vertices']
        g=sum(min(j,58-j)*x for j,x in enumerate(h));curv=F(5*g-58*v,2)
        assert g==r['G']and curv==F(r['curvature'])
        assert 2*g==58*v-sum(abs(2*j-58)*x for j,x in enumerate(h))
        saving=F(rs.count(15),64);cost=sum(env[n]for n in rs);margin=curv/v-cost
        assert margin==F(r['normalized_margin'])+saving==F(local[rs]['normalized_margin'])
        rows.append(dict(ranks=rs,link_vertices=v,h_polynomial=h,G=g,curvature=curv,normalized_cost=cost,
            normalized_saving=saving,normalized_margin=margin,absolute_margin=v*margin))
    summary=C.U.summary(rows)
    assert summary['negative_count']==summary['zero_count']==0
    assert summary['normalized_minimum']==F(1,340)and summary['absolute_minimum']==F(91,2)
    return dict(complete_profiles=1890,current_envelopes=env,rows=rows,summary=summary,
        saved_polynomial_moments_checked=True,fresh_convolution=False,complete_low9_domain=True,
        positive_terminal_extension_to_actual_cusps=True,actual_charge_bounded_by_envelope=True,
        uniform_actual_cusp_surplus='91/2')


def literal(docs):
    old={C.R.key(r['source'],r['case']):r for r in docs[M.TAG]['combined_literal_interface']['rows']}
    locals={tag:{C.R.key(r['source'],r['case']):r for r in docs[tag]['current_literal_interface']['rows']}for tag in (I.TAG,J.TAG)}
    assert len(old)==216 and all(len(x)==216 and set(x)==set(old)for x in locals.values())
    assert docs[J.TAG]['current_literal_interface']['changed_corners']==0
    rows=[];debits=credits=0
    for key,b in old.items():
        i,j=locals[I.TAG][key],locals[J.TAG][key];corners=[]
        assert len(b['corners'])==len(i['corners'])==len(j['corners'])
        for p,c in enumerate(b['corners']):
            x,y=i['corners'][p],j['corners'][p];assert c['marks']==x['marks']==y['marks']
            before=F(c['retained_payment_lower_bound']);assert F(y['debit'])==0 and F(y['retained_payment_lower_bound'])==before
            debit,credit=F(x['debit']),F(x['increment']);now=before-debit+credit
            assert now==F(x['retained_payment_lower_bound'])>=0
            debits+=bool(debit);credits+=bool(credit)
            corners.append(dict(marks=c['marks'],previous_lower=before,A14_debit='0',ideal16_debit=debit,
                new_A1_D15_increment=credit,retained_payment_lower_bound=now))
        total=sum(x['retained_payment_lower_bound']for x in corners);assert total>=1
        rows.append(dict(source=b['source'],case=b['case'],corners=corners,total_lower_bound=total))
    assert (debits,credits)==(4,2)
    return dict(configurations=216,ideal_debits=debits,finite_credits=credits,A14_debits=0,rows=rows,
        minimum_total=min(r['total_lower_bound']for r in rows),native_Grams_reconstructed=False)


def run(saved=None):
    assert all(h!='PENDING'for h in FINAL_MDS.values()),'Wait for frozen local theorem.'
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    for s,h in FINAL_MDS.items():assert sha256((ROOT/f'results/{s}.md').read_bytes()).hexdigest()==h,s
    assert docs[J.TAG]['new_costs']==dict(A14='89',B14='87',D14='181/2')
    assert docs[J.TAG]['rank14_representative_upper']=='181/2'and not docs[J.TAG]['new_cross_component_coupons']
    assert docs[J.TAG]['retained_A15_pair']==[1,13]and docs[J.TAG]['retained_A15_weight']=='1'
    assert docs[I.TAG]['new_family_costs']==dict(A='80',B='102',C='101',D='411/4')
    cp=docs[I.TAG]['new_coupon'];assert cp['large_marks']==[11]and cp['component_pair_increment']=='1/16'
    assert cp['new_actual_component_cost']=='25/16'and docs[I.TAG]['finite_accounting']['representative_A1_D15_allowance']=='25/16'
    assert F((5*((N-1)**2//4)-N*(N-1)),2)==F(docs[M.TAG]['finite_curvature'])==F(783,2)
    assert 5*(N*N//4)==int(docs[M.TAG]['strict_global_constant'])==4350
    ff=finite(saved['finite']if saved else None);ii=ideal(docs);ll=literal(docs)
    paths={Path(__file__),Path(M.__file__),Path(I.__file__),Path(J.__file__),Path(C.__file__),Path(C.U.__file__),Path(C.P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_A14_CD16_AMBIENT59_POSITIVE_IDEAL_FINITE787_OVER2',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        finite=ff,ideal=ii,combined_literal_interface=ll,ambient_dimension=59,completed_unrestricted_bound=59,
        new_dimension_bound=None,all_ideal_dimension_upper_bound=58,all_ideal_novelty_or_best_known_claim=False,
        finite_curvature='783/2',strict_global_constant='4350',
        constant_counted_once=True,triangle_minimum=2,quadrilateral_minimum=1,
        necessary_ambient59_inequality='2*v_f >= 4350 + (91/2)*v_i',minimum_finite_vertices=2175,
        minimum_finite_vertices_with_an_ideal_vertex=2198,
        scope='Complete simultaneous local rule, fresh changed-key finite certificates and positive ideal terminal/product comparison; no unrestricted dimension58 claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--preflight',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if args.preflight:
        out=enc(finite());Path('/private/tmp/'+TAG+'_preflight.json').write_text(json.dumps(out,indent=2)+'\n')
        print('CONDITIONAL_PREFLIGHT',out['exact_maximum'],len(out['above_curvature']),flush=True)
    else:
        saved=json.loads(path.read_text())if args.verify else None
        if not args.verify:assert not path.exists()
        out=run(saved)
        if args.verify:assert out==saved
        else:path.write_text(json.dumps(out,indent=2)+'\n')
        print(out['status'],out['finite']['exact_maximum'],len(out['finite']['above_curvature']),out['ideal']['summary']['negative_count'],flush=True)
