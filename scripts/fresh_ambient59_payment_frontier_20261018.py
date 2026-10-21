#!/usr/bin/env python3
"""Complete ambient59 finite/ideal frontier under the completed59 rule.

Fresh item-outer finite tables; replay checks all last-item inequalities and
attainers. Fresh ideal convolutions use the complete low9 terminal domain.
All suggested changes are explicitly conditional, without new geometry.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path
import fresh_depth6_sparse_ambient60_forest_frontier_20261018 as C

ROOT=C.ROOT;Z=C.Z;N,D=59,58
TAG='fresh_ambient59_payment_frontier_20261018'
COMPLETE='fresh_depth6_dimension59_20261018'
LOCAL='fresh_depth6_finite7_adaptive_transfer_20261018'
FAMILY14='fresh_order14_c_distance10_terminal_projection_20261018'
STEMS=(COMPLETE,C.TAG,C.S.TAG,LOCAL,FAMILY14)
ITEMS=tuple(range(1,9))+tuple(range(17,N+1))
SMALL={10:F(53,11),11:F(31,6),12:F(145,26),13:F(41,7),
       14:F(92,15),15:F(103,16),16:F(229,34)}


def envelope(r):return SMALL.get(r,Z.envelope(r))


def fixed_domain():
    out=[]
    def visit(r,left,cs):
        if r==18:
            c=dict(zip(range(9,18),cs))
            key=(2*c[10]+3*sum(c[s] for s in range(11,17))+2*c[17],
                 c[11]+c[12],c[13],c[11],c[14],c[9],c[10],c[16])
            out.append((tuple(cs),N-left,key));return
        for n in range(left//r+1):visit(r+1,left-r*n,cs+[n])
    visit(9,N,[])
    return sorted(out)


def table(key):
    ws={r:C.item(r,key) for r in ITEMS};v=[0]+[-10**9]*N
    for r,w in ws.items():
        for n in range(r,N+1):v[n]=max(v[n],v[n-r]+w)
    last=[0]+[next(r for r,w in ws.items() if r<=n and v[n]==v[n-r]+w) for n in range(1,N+1)]
    return dict(key=list(key),values=v,last=last)


def check_table(t,key):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    ws={r:C.item(r,key) for r in ITEMS};checks=0
    for n in range(1,N+1):
        for r,w in ws.items():
            if r<=n:assert v[n]>=v[n-r]+w;checks+=1
        r=last[n];assert r in ws and r<=n and v[n]==v[n-r]+ws[r]
    return checks


def finite(previous,beta,saved=None):
    rep=Z.B.B.representative;types=[t for t in Z.B.B.types() if t[1]<=N]
    for a in types:
        assert C.upper_base(rep(a))>=C.actual_base(a)>=0
        for b in types:assert C.upper_cross(rep(a),rep(b))>=C.actual_cross(a,b)==C.actual_cross(b,a)>=0
    dom=fixed_domain();assert len(dom)==961
    prior={tuple(r['counts']) for r in previous['alternatives'] if sum(i*c for i,c in zip(range(9,18),r['counts']))<=N}
    assert {cs for cs,used,key in dom}==prior
    keys=sorted({key for cs,used,key in dom})
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {k:table(k) for k in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tables[k],k) for k in keys)
    oldtabs={tuple(t['key']):t for t in previous['tables']}
    assert all(tables[k]['values']==oldtabs[k]['values'][:N+1] for k in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in ITEMS]
    assert all(C.upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];ic=0
    for cs,used,key in dom:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in fixed)==used
        assert all(C.upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(ITEMS,residual):
            assert C.item(r,key)==C.upper_base(t)+sum(C.upper_cross(t,b) for b in fixed);ic+=1
        rem=N-used;tab=tables[key];upper=sum(C.upper_base(t) for t in fixed)+tab['values'][rem]
        fs=fixed[:]
        while rem:
            r=tab['last'][rem];fs.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in fs)==N and C.value(fs,True)==upper and C.value(fs)<=upper
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in fs),
            octupled_upper=upper,octupled_actual=C.value(fs)))
    maximum=max(r['octupled_upper'] for r in rows);assert maximum==3180
    witness_lists=[['D16','D17','D26'],['A17','D16','D26'],['D16','D16','D27'],
        ['D17','D17','D25'],['A17','D17','D25'],['A17','A17','D25'],
        ['A1']*3+['D5']*4+['D12']*3,['A2']*10+['D13']*3,
        ['A15']*3+['D14']]
    witnesses=[]
    for names in witness_lists:
        ts=[(s[0],int(s[1:])) for s in names];assert sum(r for k,r in ts)==N
        v=F(C.value(ts),8);assert v>beta
        witnesses.append(dict(forest=sorted(names),cost=str(v),deficit=str(v-beta)))
    assert all(F(r['cost'])==F(maximum,8) for r in witnesses[:3])
    above=sorted((r for r in rows if F(r['octupled_upper'],8)>beta),key=lambda r:(-r['octupled_upper'],r['forest']))
    assert len(above)==15
    # Further deep7 reductions alone cannot affect these actual small forests.
    unchanged=[r for r in witnesses if r['forest'].count('D12')==3 or r['forest'].count('D13')==3]
    assert [r['cost'] for r in unchanged]==['791/2','394']
    assert all(all(int(s[1:])<=13 for s in r['forest']) for r in unchanged)
    assert C.actual_cross(('D',5),('D',12))==90 and C.actual_cross(('A',2),('D',13))==40
    return dict(exact_maximum=str(F(maximum,8)),curvature=str(beta),finite_deficit=str(F(maximum,8)-beta),
        complete_fixed_vectors=len(dom),interaction_keys=len(keys),simultaneous_type_pairs=len(types)**2,
        residual_ranks=list(ITEMS),fixed_rank17_type='A17',residual_rank17_type='D17',
        residual_item_checks=ic,Bellman_upper_inequalities=checks,attaining_predecessors=N*len(keys),
        tables=[tables[k] for k in keys],alternatives=rows,above_curvature=above,
        representative_alternatives_above_curvature=len(above),actual_obstructing_forests=witnesses,
        actual_singleton_D59=str(F(C.value([('D',59)]),8)),
        actual_rank15_costs=dict(A15='101',B15='97',D15='201/2'),rank15_upper='101',
        all_current_cross_component_coupons_counted=True,D5_D12_cost='45/4',A2_D13_cost='5',
        deep7_only_unchanged_actual_obstructions=unchanged,
        new_Bellman_generation=True,replay_optimizes=False,old_dimension_prefix_comparisons=len(keys),
        scope='Complete current-rule rank59 upper and actual formal attainment. Listed actual obstructing forests are witnesses, not an enumeration of every actual type assignment. No geometric realization asserted.')


def scalar(rs):
    h=[1]
    for r in rs:
        out=[0]*(len(h)+r)
        for i,x in enumerate(h):
            for j in range(r+1):out[i+j]+=x
        h=out
    v=prod(r+1 for r in rs);assert sum(rs)==D and len(h)==D+1 and h==h[::-1] and sum(h)==v
    g=sum(min(j,D-j)*x for j,x in enumerate(h))
    assert 2*g==D*v-sum(abs(2*j-D)*x for j,x in enumerate(h))
    cumulative=0;deviation=F()
    for j in range(D):
        cumulative+=h[j];deviation+=F(cumulative if 2*j<D else v-cumulative,v)
    assert F(g,v)==F(D,2)-deviation
    curvature=F(5*g-D*v,2);cost=sum(map(envelope,rs),F());margin=curvature/v-cost
    return dict(ranks=list(rs),h_polynomial=h,link_vertices=v,G=g,curvature=str(curvature),
        normalized_cost=str(cost),normalized_margin=str(margin),absolute_margin=str(v*margin))


def ideal(previous):
    # Bound every changed envelope to the entire frozen rank59 domain.
    assert all(sum(map(envelope,r['ranks']),F())==F(r['normalized_cost']) for r in previous['rows'])
    profiles=sorted((low,)+p if low else p for low in range(9) for p in Z.high_parts(D-low))
    assert len(profiles)==len(set(profiles))==1890
    rows=[scalar(p) for p in profiles];summary=Z.summarize(rows)
    assert summary['negative_count']==10 and summary['zero_count']==0
    assert summary['normalized_minimum']=='-1537/14000' and summary['absolute_minimum']=='-153700'
    active=sorted({r for row in rows if F(row['normalized_margin'])<0 for r in row['ranks'] if r>=10})
    coordinated={12:F(1),13:F(7,4),14:F(1,4),15:F(1,2),16:F(3,4)}
    scenarios=[Z.conditional(rows,{r:delta for r in active}) for delta in (F(1,4),F(1,2),F(1))]
    scenarios.append(Z.conditional(rows,coordinated))
    assert scenarios[-1]['summary']['negative_count']==scenarios[-1]['summary']['zero_count']==0
    thresholds=[Z.threshold(rows,r,{}) for r in active]
    coordinated_thresholds=[Z.threshold(rows,r,{k:v for k,v in coordinated.items() if k!=r}) for r in coordinated]
    return dict(complete_profiles=1890,rows=rows,summary=summary,
        current_small_rank_envelopes={str(r):str(envelope(r)) for r in range(1,17)},
        exact_family_tail=[dict(order=r+1,costs={k:Z.family_cost(k,r+1) for k in 'ABCD'},envelope=str(envelope(r))) for r in range(17,D+1)],
        complete_low9_domain=True,high_factor_concentration_used=False,fresh_convolutions=True,
        active_negative_ranks=active,conditional_scenarios=scenarios,
        single_rank_thresholds=thresholds,coordinated_thresholds=coordinated_thresholds,
        negative_absolute_margin_not_extended_to_arbitrary_profiles=True,
        scope='Complete terminal scalar frontier under unchanged all-family envelopes. Conditional drops require new compatible geometry and every added finite coupon must be budgeted. No new local rule or dimension bound.')


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[COMPLETE]['dimension_bound']==59
    assert docs[C.TAG]['finite']['exact_maximum']=='404'
    assert docs[C.S.TAG]['finite']['D5_D12_actual_cost']=='45/4'
    assert docs[FAMILY14]['new_order14_family_costs']==dict(A=70,B=82,C=81,D=82)
    assert docs[LOCAL]['counts']['formulas']==dict(A='max(0,n-19)',B='max(0,n-20)',D='max(0,n-21)')
    g0=sum(min(j,D-j) for j in range(D+1));beta=F(5*g0-D*(D+1),2)
    extremum=sum(2*j-D for j in range(D+1) if 2*j>D)
    assert (g0,beta,extremum)==(841,F(783,2),870)
    paths={Path(__file__),Path(C.__file__),Path(Z.__file__)}
    paths.update(ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return dict(status='PASS_AMBIENT59_PAYMENT_FRONTIER',sources=sources,
        ambient_dimension=N,link_dimension=D,completed_unrestricted_dimension_bound=59,
        new_dimension_bound=None,all_local_weights_unchanged=True,
        current_order14_family_costs=docs[FAMILY14]['new_order14_family_costs'],
        finite=finite(docs[C.TAG]['finite'],beta,saved['finite'] if saved else None),ideal=ideal(docs[C.TAG]['ideal']),
        finite_curvature=str(beta),strict_global_constant=str(5*extremum),constant_counted_once=True,
        current_rule_does_not_close_ambient59=True,
        scope='Exact changed-dimension diagnostic only: complete finite upper/attainment certificates and freshly convolved terminal ideal profiles. No new geometry, arbitrary-profile negative-surplus extension or dimension58 theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    if not args.verify:assert not path.exists(),'Use --verify for an existing certificate.'
    out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'keys',out['finite']['interaction_keys'],
        'ideal negatives',out['ideal']['summary']['negative_count'],flush=True)
