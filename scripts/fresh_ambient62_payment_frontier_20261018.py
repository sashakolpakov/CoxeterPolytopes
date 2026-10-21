#!/usr/bin/env python3
"""Ambient62 current payment deficits and explicitly conditional targets.

No project arithmetic imports. Current finite certificates are independently
regenerated rank-first and compared with the frozen rank63 tables. Exact
rank-specific cusp tails use a larger complete low9-only domain.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ambient62_payment_frontier_20261018'
COMPLETED='fresh_sparse_coupon_order11_13_dimension62_20261018'
FINITE='fresh_sparse_coupon_rank63_forest_budget_20261018'
IDEAL='fresh_ambient63_ideal_margin_closure_20261018'
TAIL='fresh_ideal5_small_factor_payments_20261017'
STEMS=(COMPLETED,FINITE,IDEAL,TAIL)
N,D=62,61
BETA=F(434)

def original_base(t):
    kind,r=t
    if r==16:return 208 if kind=='B' else 228
    if r==17:return {'A':228,'B':210,'D':250}[kind]
    if r<=15:return r*(r-1)
    old=16*r-(66 if kind=='D' and r>=20 else 64 if kind=='D' else 68)
    deep=max(0,r-(22+{'A':8,'B':9,'D':10}[kind]))
    return old-2*deep


def representative(t):
    kind,r=t
    if t==('A',17):return t
    if r==6:return ('E',6)
    if r==9:return ('B',9)
    return ('A' if r<=8 else 'D',r)


def old_cross(a,b):
    out=0
    for small,large in ((a,b),(b,a)):
        if small[1]<=5:
            if large==('D',10):out+=2*small[1]
            elif large[0]=='D' and 11<=large[1]<=16:out+=3*small[1]
            elif large in (('A',13),('A',15)):out+=small[1]
            elif large==('A',17):out+=2*small[1]
        if (small,large) in ((('A',1),('D',11)),(('A',1),('D',12)),
                (('A',2),('B',10)),(('A',4),('D',13)),
                (('A',1),('A',11)),(('A',2),('A',11))):out+=2
        if (small,large) in ((('A',2),('D',11)),(('E',6),('D',11))):out+=4
        if (small,large)==(('A',1),('B',10)):out+=1
        if small in (('A',1),('A',2)) and large==('D',14):out+=small[1]
    return out


# These are doubled total costs, not the weight of an individual marked pair.
NEW_CROSS={
    (('A',4),('A',12)):8,
    (('A',4),('A',13)):8,
    (('D',5),('A',12)):4,
    (('A',2),('A',13)):4,
    (('A',2),('A',14)):4,
    (('E',6),('A',11)):8,
    (('A',2),('B',9)):2,
    (('A',2),('D',10)):2,
    (('A',3),('D',10)):1,
}


def actual_cross(a,b):
    return old_cross(a,b)+NEW_CROSS.get((a,b),0)+NEW_CROSS.get((b,a),0)


def upper_cross(a,b):
    # This additional allowance is ONLY in the representative upper bound.
    # The actual E6/D11 coupon remains total2; E6/A11's new total4 is absorbed.
    extra=4 if {a,b}=={('E',6),('D',11)} else 0
    return actual_cross(a,b)+extra


def domain():
    out=[]
    def visit(rank,left,counts):
        if rank==18:
            c=dict(zip(range(9,18),counts))
            key=(2*c[10]+3*sum(c[r] for r in range(11,17))+2*c[17],
                 c[11]+c[12],c[13],c[11],c[14],c[9],c[10])
            out.append(dict(counts=counts,used_rank=N-left,key=list(key)))
            return
        for count in range(left//rank+1):visit(rank+1,left-rank*count,counts+[count])
    visit(9,N,[])
    return sorted(out,key=lambda r:r['counts'])


def finite_scenario(frozen,changes,saved=None):
    assert set(changes)<={'rank15','rank16','D17','large_D33','D13_zero_coupons'}
    def base(t):
        k,r=t;drop=F()
        if r==15:drop+=F(changes.get('rank15',0))
        if r==16 and k!='B':drop+=F(changes.get('rank16',0))
        if t==('D',17):drop+=F(changes.get('D17',0))
        if k=='D' and r>=33:drop+=F(changes.get('large_D33',0))
        ans=F(original_base(t))-2*drop
        assert ans.denominator==1 and ans>=0
        return int(ans)
    def cross(a,b,upper=False):
        value=(upper_cross if upper else actual_cross)(a,b)
        if changes.get('D13_zero_coupons'):
            for small,large in ((a,b),(b,a)):
                if large==('D',10):value+=2 if small==('A',2) else 1 if small==('A',3) else 0
        return value
    def cost(forest,upper=False):
        return sum(map(base,forest))+sum(cross(a,b,upper) for i,a in enumerate(forest) for b in forest[i+1:])
    types=([('A',r) for r in range(1,N+1)]+[('B',r) for r in range(2,N+1)]+
           [('D',r) for r in range(4,N+1)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    for t in types:
        rep=representative(t);assert base(rep)>=base(t)
        assert all(cross(rep,representative(u),True)>=cross(t,u) for u in types)
    dom=domain();assert dom==[r for r in frozen['large_domain'] if r['used_rank']<=N]
    keys=sorted({tuple(r['key']) for r in dom});tables={};items={};checks=Counter()
    if saved is not None:assert [t['key'] for t in saved['tables']]==list(map(list,keys))
    for index,key in enumerate(keys):
        broad,extra,a4,d11,d14,b9,d10=key;weights={}
        for r in tuple(range(1,9))+tuple(range(17,N+1)):
            t=representative(('A' if r<=8 else 'D',r))
            value=base(t)+(r*broad if r<=5 else 0)
            value+=(2*extra+d14 if r==1 else 4*d11+2*d14+2*b9+2*d10 if r==2
                    else d10 if r==3 else 2*a4 if r==4 else 8*d11 if r==6 else 0)
            if changes.get('D13_zero_coupons'):value+=2*d10 if r==2 else d10 if r==3 else 0
            weights[r]=value
        items[key]=weights
        if saved is None:
            # Every-rank generation is separate from the frozen item-first DP.
            values,last=[0],[0]
            for n in range(1,N+1):
                value,r=max((values[n-r]+w,r) for r,w in weights.items() if r<=n)
                values.append(value);last.append(r)
        else:values,last=saved['tables'][index]['values'],saved['tables'][index]['last']
        assert len(values)==len(last)==N+1 and values[0]==last[0]==0
        assert all(type(x) is int for x in values+last)
        for n in range(1,N+1):
            for r,w in weights.items():
                if r<=n:
                    assert values[n]>=values[n-r]+w;checks['bellman_upper']+=1
            r=last[n];assert r in weights and r<=n and values[n]==values[n-r]+weights[r]
            checks['attaining_predecessors']+=1
        tables[key]=dict(key=list(key),values=values,last=last)
    if not changes:
        old={tuple(t['key']):t for t in frozen['tables']}
        assert all(t['values']==old[k]['values'][:N+1] for k,t in tables.items())
    alternatives=[]
    for row in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r)
               for r,c in zip(range(9,18),row['counts']) for _ in range(c)]
        assert all(cross(a,b,True)==0 for i,a in enumerate(large) for b in large[i+1:])
        key=tuple(row['key']);table=tables[key]
        for r,w in items[key].items():
            t=representative(('A' if r<=8 else 'D',r))
            assert w==base(t)+sum(cross(t,u,True) for u in large);checks['literal_interactions']+=1
        rem=N-row['used_rank'];value=sum(map(base,large))+table['values'][rem];forest=list(large)
        while rem:
            r=table['last'][rem];forest.append(representative(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N and cost(forest,True)==value and cost(forest)<=value
        alternatives.append(dict(counts=row['counts'],doubled_upper=value,actual_doubled=cost(forest),
                                 forest=sorted(k+str(r) for k,r in forest)))
        checks['global_alternatives']+=1
    maximum=max(r['doubled_upper'] for r in alternatives)
    actual_attainers=sorted({tuple(r['forest']) for r in alternatives if r['doubled_upper']==r['actual_doubled']==maximum})
    assert actual_attainers # Equality makes this an exact actual-rule maximum, not just an envelope.
    return dict(conditional_changes=changes,uses_current_rule=not changes,geometric_change_proved_here=False,
                exact_maximum=str(F(maximum,2)),curvature=str(BETA),finite_slack=str(BETA-F(maximum,2)),
                actual_attaining_forests=actual_attainers,complete_large_vectors=len(dom),complete_keys=len(keys),
                simultaneous_type_pairs=len(types)**2,checks=dict(checks),
                tables=[tables[k] for k in keys],alternatives=alternatives,
                largest_alternatives=sorted(alternatives,key=lambda r:(-r['doubled_upper'],r['forest']))[:24])


def profiles(left,minimum=1,prefix=(),low=False,high=False,cap_high=False):
    if not left:
        yield prefix;return
    for r in range(minimum,left+1):
        if (r<9 and low) or (cap_high and r>=18 and high):continue
        yield from profiles(left-r,r,prefix+(r,),low or r<9,high or r>=18,cap_high)


def family_cost(kind,m):
    if kind=='A':return 5*m
    offset,cap,first={'C':(-15,14,6),'B':(-6,13,8),'D':(3,12,10)}[kind]
    return 5*m+offset+sum(min(m-k,cap) for k in range(first,first+4))


def envelope(r,exact_tail):
    if r<=9:return F(r,2)
    small={10:F(54,11),11:F(21,4),12:F(147,26),13:F(6),14:F(94,15),15:F(107,16),16:F(118,17)}
    if r in small:return small[r]
    return max(F(family_cost(k,r+1),r+1) for k in 'ABCD') if exact_tail else F(163,23)


def scalar(rs,exact_tail):
    h=[1]
    for r in rs:
        out,window=[],0
        for j in range(len(h)+r):
            if j<len(h):window+=h[j]
            if 0<=j-r-1<len(h):window-=h[j-r-1]
            out.append(window)
        h=out
    V=prod(r+1 for r in rs);assert len(h)==D+1 and h==h[::-1] and sum(h)==V
    G=sum(v*min(j,D-j) for j,v in enumerate(h))
    # Both parity cases use the doubled absolute-deviation identity.
    assert 2*G==D*V-sum(v*abs(2*j-D) for j,v in enumerate(h))
    cumulative=0;distance=F()
    for j in range(D):
        cumulative+=h[j]
        if j<D//2:distance+=F(cumulative,V)
        elif j>D//2:distance+=F(V-cumulative,V)
        else:distance+=F(1,2) if D%2 else F(V-cumulative,V)
    assert F(G,V)==F(D,2)-distance
    curvature=F(5*G-D*V,2);cost=sum((envelope(r,exact_tail) for r in rs),F());margin=curvature/V-cost
    return dict(ranks=list(rs),h_polynomial=h,link_vertices=V,G=G,curvature=str(curvature),
                normalized_cost=str(cost),normalized_margin=str(margin),absolute_margin=str(V*margin))


def summarize(rows):
    order=sorted(rows,key=lambda r:(F(r['normalized_margin']),r['ranks']))
    return dict(count=len(rows),negative=sum(F(r['normalized_margin'])<0 for r in rows),
                zero=sum(F(r['normalized_margin'])==0 for r in rows),minimum=order[0],
                minimum_absolute=min(rows,key=lambda r:(F(r['absolute_margin']),r['ranks'])),
                negative_rows=[r for r in order if F(r['normalized_margin'])<0])


def ideal_scenario(rows,drops):
    assert all(10<=r<=16 for r in drops) # The unchanged low9 reduction remains valid.
    out=[]
    for row in rows:
        gain=sum(F(row['ranks'].count(r),r+1)*drop for r,drop in drops.items())
        margin=F(row['normalized_margin'])+gain
        out.append(dict(ranks=row['ranks'],normalized_gain=str(gain),normalized_margin=str(margin),
                        absolute_margin=str(row['link_vertices']*margin)))
    return dict(conditional_total_cost_drops_by_affine_rank={str(r):str(v) for r,v in drops.items()},
                summary=summarize(out),rows=out,geometric_coverage_proved=False,
                additional_finite_coupons_included=False)


def threshold(rows,varied,fixed):
    bounds=[];unfixable=[];unfixable_zeros=[]
    for row in rows:
        rs=row['ranks'];margin=F(row['normalized_margin'])+sum(F(rs.count(r),r+1)*v for r,v in fixed.items())
        if margin>0:continue
        if margin==0:
            if varied not in rs:unfixable_zeros.append(rs)
            continue
        if varied not in rs:unfixable.append(dict(ranks=rs,margin=str(margin)))
        else:bounds.append((-(varied+1)*margin/rs.count(varied),rs))
    bound,rs=max(bounds) if bounds else (F(),[])
    return dict(varied_affine_rank=varied,fixed_drops={str(r):str(v) for r,v in fixed.items()},
                threshold_for_nonnegative=str(bound),active_profile=rs,unfixable_profiles=unfixable,
                unfixable_zero_profiles=unfixable_zeros,
                strict_note='A strictly larger drop gives strict positivity only if no unfixable negative or zero row remains.')


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[COMPLETED]['dimension_bound']==62
    assert docs[FINITE]['finite_maximum']=='449'
    assert docs[FINITE]['representative_allowance']['actual_total']=='2'
    assert docs[FINITE]['representative_allowance']['upper_total']=='4'
    assert [(tuple(map(tuple,r['types'])),F(r['added_component_pair_cost'])*2)
            for r in docs[FINITE]['actual_cross_coupon_lists']]==list(NEW_CROSS.items())
    assert all(F(r['cost'])==envelope(r['rank'],False) for r in docs[IDEAL]['envelopes'])
    old_tail=docs[TAIL]['basic_affine_cost']
    for row in old_tail['classical_checks']:
        assert row['cost']==family_cost(row['kind'],row['order'])
    assert BETA==F(5*sum(min(j,D-j) for j in range(N))-N*D,2)
    current=[scalar(rs,False) for rs in profiles(D,cap_high=True)]
    exact=[scalar(rs,True) for rs in profiles(D)]
    assert len(current)==len({tuple(r['ranks']) for r in current})==2004
    assert len(exact)==len({tuple(r['ranks']) for r in exact})==2641
    s0,s1=summarize(current),summarize(exact)
    assert (s0['negative'],s1['negative'])==(64,62) and s0['zero']==s1['zero']==0
    assert s0['minimum']['normalized_margin']==s1['minimum']['normalized_margin']=='-3379/17408'
    assert s1['minimum']['ranks']==[15,15,15,16]
    assert all(max(r['ranks'])<=16 for r in s1['negative_rows'])
    singleton=next(r for r in exact if r['ranks']==[61]);assert singleton['normalized_margin']=='73/62'
    changes=[{},dict(D13_zero_coupons=True),dict(D17='4'),
             dict(rank15='1',rank16='2',D17='4',D13_zero_coupons=True),
             dict(rank15='1',rank16='2',D17='5',D13_zero_coupons=True),
             dict(rank15='1',rank16='2',D17='4',large_D33='1',D13_zero_coupons=True)]
    finite=[finite_scenario(docs[FINITE],c,None if saved is None else saved['finite_scenarios'][i])
            for i,c in enumerate(changes)]
    assert finite[0]['exact_maximum']=='441' and finite[0]['actual_attaining_forests']==[('D17','D17','D28')]
    assert finite[1]['exact_maximum']=='441'
    drops=[{12:F(1,2)},
           {12:F(1,2),15:F(1),16:F(1,2)},
           {12:F(1,2),13:F(1,2),15:F(1),16:F(1,2)},
           {10:F(1,2),12:F(1,2),13:F(1,2),15:F(1),16:F(1,2)},
           {10:F(1,2),12:F(1,2),13:F(1,2),15:F(1),16:F(3,4)},
           {10:F(1,2),11:F(1,2),12:F(1,2),13:F(1,2),15:F(1),16:F(1,2)}]
    ideals=[ideal_scenario(exact,c) for c in drops]
    assert [r['summary']['negative'] for r in ideals]==[45,8,3,1,0,0]
    assert ideals[4]['summary']['minimum']['normalized_margin']=='663/71680'
    tail_rows=[dict(order=m,by_family={k:family_cost(k,m) for k in 'ABCD'},
                    normalized_maximum=str(envelope(m-1,True))) for m in range(18,N+1)]
    thresholds=[threshold(exact,16,{10:F(1,2),12:F(1,2),13:F(1,2),15:F(1)}),
                threshold(exact,15,{12:F(1,2)})]
    assert thresholds[0]['threshold_for_nonnegative']=='87/170'
    assert not thresholds[0]['unfixable_profiles'] and not thresholds[0]['unfixable_zero_profiles']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
          for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/'results/fresh_distance14_unrestricted_dimension156_20261017.md',
              ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md']
    return dict(status='PASS_AMBIENT62_CURRENT_AND_CONDITIONAL_PAYMENT_FRONTIER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                ambient_dimension=N,inherited_completed_ceiling=62,finite_curvature=str(BETA),
                finite_scenarios=finite,constant_envelope_rows=current,constant_envelope_summary=s0,
                exact_tail_rows=exact,exact_tail_summary=s1,exact_tail_factor_costs=tail_rows,
                ideal_conditional_scenarios=ideals,
                ideal_thresholds=thresholds,
                strict_global_constant=str(5*(N*N//4)),new_dimension_bound=None,
                scope='Complete current rank62 finite count and rank61 ideal deficits. Constant-envelope2004 domain uses9/17 concentration; exact inherited family tails use the larger complete2641 low9-only domain, with no high concentration assumed. Conditional cost reductions and D13 coupon upgrades are counting scenarios only, not new geometric rules or a dimension61 theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    out=json.loads(json.dumps(run(saved)))
    if args.verify:assert out==saved
    else:
        assert not path.exists(),'Existing output; use --verify.'
        path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'finite',[r['exact_maximum'] for r in out['finite_scenarios']],
          'ideal',out['exact_tail_summary']['count'],out['exact_tail_summary']['negative'],flush=True)


