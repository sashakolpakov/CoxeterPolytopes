#!/usr/bin/env python3
"""Complete actual frontier for the coupled D17/A16 transfer and A13 guard.

Inherited global representative tables locate every possibly costly class.
Fresh actual residual tables certify all39 classes, with two separately
bounded increasing families. The eight B13/D13+three-rank15 forests remain.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_d17_leaf_a16_context_coupled_transfer_20261018 as P

ROOT=P.ROOT;C=P.C;N=P.N;S=P.S
TAG='fresh_d17_a16_coupled_ambient58_frontier_20261018'
STEMS=tuple(dict.fromkeys((P.TAG,*P.STEMS)))
enc=P.enc;actual_base=P.actual_base;actual_cross=P.actual_cross;context=P.context;value=P.value
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def quarter_family(fs):return all(Counter(fs)[t]for t in [('A',1),('A',2),('A',4),('A',10)])

def representative_cover(docs):
    base_guard={('A',13):64,('A',14):32,('D',14):64,('A',15):32,('D',15):32,('A',16):80,('D',16):64}
    gaps=[]
    for t in C.TYPES:
        upper=C.upper_base(C.representative(t))+(22 if t==('A',17)else 0)
        actual=actual_base(t)+base_guard.get(t,0);assert actual<=upper,(t,actual,upper)
        gaps.append(dict(type=t,maximum_guarded_base=actual,representative_upper=upper))
        for u in C.TYPES:
            pair=frozenset((t,u));extra=24 if pair==frozenset((('A',1),('D',5)))else 32 if pair==frozenset((('A',4),('A',12)))else 0
            assert actual_cross(t,u)+extra<=C.upper_cross(C.representative(t),C.representative(u))
    old=docs[C.TAG]['finite'];assert len(old['expanded_alternatives'])==3298
    costly={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count'])for r in old['expanded_alternatives']if r['sixtyfourths_upper']+22*r['counts'][-1]>377*64}
    assert costly=={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count'])for r in docs[N.TAG]['finite']['classes']}
    assert len(costly)==39
    # N's extra40*A1*D5 is the sole earlier context not dominated by C.
    # Its complete conditional family is separately bounded, outside the
    # new quarter family; all other bases/crosses are already dominated.
    overlap=[dict(counts=r['counts'],upper=F(r['upper']+22*r['counts'][-1],64))for r in docs[N.TAG]['increasing_family_upper']['rows']]
    assert len(overlap)==759 and max(r['upper']for r in overlap)==F(2837,8)<377
    return dict(actual_type_base_domination=gaps,actual_type_pair_domination=len(C.TYPES)**2,
        inherited_expanded_classes=3298,complete_costly_classes=39,D17_context_family=overlap,
        D17_context_family_upper='2837/8')

def increasing_family(docs):
    old=docs[C.TAG]['finite'];tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in old['tables']}
    assert len(tabs)==6596
    rows=[];options=0
    for cs,used,key in C.fixed_domain():
        if not cs[1]:continue
        large=C.fixed_components(cs,key)
        for k5 in range((58-used)//5+1):
            fixed=large+[('A',5)]*k5;rem=58-used-5*k5
            if rem<7:continue
            fc=sum(C.upper_base(t)for t in fixed)+sum(C.upper_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            no2=tabs[key,k5,True];best=None
            for a1 in range(1,rem-5):
                for a4 in range(1,(rem-a1-2)//4+1):
                    for a2 in range(1,(rem-a1-4*a4)//2+1):
                        left=rem-a1-4*a4-2*a2
                        z=fc+a1*C.item(1,key,k5)+a2*C.item(2,key,k5)+a4*C.item(4,key,k5)+no2['values'][left]
                        if a1>=2:z+=36*a2*(a2-1)
                        if cs[0]:z+=64*a2*k5+40*a1*k5
                        z+=22*cs[-1]+64*a2*a4;options+=1
                        q=(z,a1,a2,a4,left)
                        if best is None or q>best:best=q
            assert best is not None
            rows.append(dict(counts=cs,key=key,k5=k5,residual_rank=rem,upper=F(best[0],64),explicit_counts=best[1:]))
    assert len(rows)==539 and options==37985
    high=[r for r in rows if r['upper']>377];assert len(high)==14
    classes={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count']):r for r in docs[N.TAG]['finite']['classes']}
    refined=[];fibers=options2=0
    for u in high:
        cl=classes[tuple(u['counts']),u['k5'],u['key'][9]];rr=[];rem=cl['residual_rank']
        for r in cl['rows']:
            fixed=list(map(tuple,r['fixed_actual']))
            if ('A',10)not in fixed:continue
            residual=list(map(tuple,r['residual_types']));fibers+=1
            assert not set(residual)&{('A',1),('A',2),('A',9),('A',10),('A',11),('A',12),('D',5),('A',13),('A',14),('D',14),('A',15),('D',15),('A',16),('D',16),('A',17)}
            for t,w in zip(residual,r['weights'],strict=True):
                assert actual_base(t)+sum(actual_cross(t,b)for b in fixed)<=w
                assert all(actual_cross(t,a)==0 for a in [('A',1),('A',2),('A',4)])
            best=None
            for a1 in range(1,rem-5):
                for a4 in range(1,(rem-a1-2)//4+1):
                    for a2 in range(1,(rem-a1-4*a4)//2+1):
                        left=rem-a1-4*a4-2*a2
                        if r['values'][left]is None:continue
                        fs=fixed+[('A',1)]*a1+[('A',2)]*a2+[('A',4)]*a4
                        assert quarter_family(fs)
                        bonus=context(fs)
                        assert all(context(fs+[t])==bonus for t in residual if t!=('A',4))
                        z=value(fs)+r['values'][left];options2+=1
                        q=(z,a1,a2,a4,left)
                        if best is None or q>best:best=q
            assert best is not None
            rr.append(dict(fixed_actual=fixed,upper=F(best[0],64),explicit_counts=best[1:]))
        assert rr
        refined.append(dict(counts=u['counts'],k5=u['k5'],a12=u['key'][9],rows=rr,upper=max(r['upper']for r in rr)))
    assert fibers==77 and options2==23376 and max(r['upper']for r in refined)==F(2861,8)<377
    return dict(complete_fixed_classes=539,complete_explicit_options=options,coarse_rows=rows,
        coarse_classes_above_curvature=14,actual_fixed_fibers=77,actual_explicit_options=options2,
        actual_refinements=refined,family_upper='2861/8',attainment_claimed=False,
        remaining_table_scope='Inherited actual residual upper relaxation; all actual A1/A2/A4 fixed explicitly, remaining actual types exclude them.')

def complement(docs,saved=None):
    classes=[];fibers=tables=options=inequalities=predecessors=0
    for ci,cl in enumerate(docs[N.TAG]['finite']['classes']):
        rows=[];rem=cl['residual_rank']
        for ri,old in enumerate(cl['rows']):
            prior=None if saved is None else saved['classes'][ci]['rows'][ri]
            fixed=list(map(tuple,old['fixed_actual']));co=Counter(fixed);types=list(map(tuple,old['residual_types']))
            assert ('A',1)not in types and ('A',2)not in types
            assert all(actual_cross(a,b)==0 for a in types for b in types)
            fc=sum(actual_base(t)for t in fixed)+sum(actual_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            w1=actual_base(('A',1))+sum(actual_cross(('A',1),t)for t in fixed)
            w2=actual_base(('A',2))+sum(actual_cross(('A',2),t)for t in fixed)
            cache={};certs=[]
            def table(has1,no4,no3):
                nonlocal tables,inequalities,predecessors
                key=(has1,no4,no3)
                if key in cache:return cache[key]
                ts=[t for t in types if not(no4 and t==('A',4))and not(no3 and t==('A',3))]
                ws=[actual_base(t)+sum(actual_cross(t,b)for b in fixed)+(32*co['A',12]if has1 and t==('A',4)else 0)for t in ts]
                if prior is None:
                    v=[0]+[None]*rem
                    for t,w in zip(ts,ws):
                        for rank in range(t[1],rem+1):
                            if v[rank-t[1]]is not None:
                                z=v[rank-t[1]]+w
                                if v[rank]is None or z>v[rank]:v[rank]=z
                    last=[None]+[next((i for i,(t,w)in enumerate(zip(ts,ws))if t[1]<=rank and v[rank-t[1]]is not None and v[rank]==v[rank-t[1]]+w),None)for rank in range(1,rem+1)]
                else:
                    p=next(p for p in prior['tables']if p['key']==list(key));v=p['values'];last=p['last']
                    assert p['types']==enc(ts)and p['weights']==ws
                assert len(v)==len(last)==rem+1 and v[0]==0 and last[0]is None
                for rank in range(1,rem+1):
                    possibilities=[]
                    for i,(t,w)in enumerate(zip(ts,ws)):
                        if t[1]<=rank and v[rank-t[1]]is not None:
                            assert v[rank]is not None and v[rank]>=v[rank-t[1]]+w
                            possibilities.append(i);inequalities+=1
                    if possibilities:
                        assert last[rank]in possibilities
                        i=last[rank];assert v[rank]==v[rank-ts[i][1]]+ws[i];predecessors+=1
                    else:assert v[rank]is None and last[rank]is None
                tables+=1;cache[key]=(ts,ws,v,last)
                certs.append(dict(key=key,types=ts,weights=ws,values=v,last=last))
                return cache[key]
            best=None
            for a1 in range(rem+1):
                for a2 in range((rem-a1)//2+1):
                    left=rem-a1-2*a2;fs=fixed+[('A',1)]*a1+[('A',2)]*a2
                    no4=bool(a1 and a2 and co['A',10])
                    modes=[(False,0)]if a1 or a2 or not co['D',14]else[(True,0),(False,1)]
                    for no3,res3 in modes:
                        if left<3*res3:continue
                        ts,ws,v,last=table(bool(a1),no4,no3);rest=left-3*res3
                        if v[rest]is None:continue
                        extra=0
                        if res3:
                            assert ('A',3)in ts;extra=ws[ts.index(('A',3))]+64*co['D',14]
                        z=fc+a1*w1+a2*w2+context(fs)+extra+v[rest];options+=1
                        if best is None or z>best[0]:best=(z,a1,a2,no4,no3,res3,rest)
            assert best is not None
            z,a1,a2,no4,no3,res3,left=best;fs=fixed+[('A',1)]*a1+[('A',2)]*a2+[('A',3)]*res3
            ts,ws,v,last=table(bool(a1),no4,no3)
            while left:
                t=ts[last[left]];fs.append(t);left-=t[1]
            assert sum(t[1]for t in fs)==58 and not quarter_family(fs)and value(fs)==z
            if prior is not None:assert len(certs)==len(prior['tables'])
            rows.append(dict(fixed_actual=fixed,tables=certs,maximum=F(z,64),attaining_forest=fs,winning_option=best[1:]));fibers+=1
        classes.append(dict(counts=cl['counts'],k5=cl['fixed_A5_count'],a12=cl['fixed_A12_count'],residual_rank=rem,
            rows=rows,maximum=max(r['maximum']for r in rows)))
    assert len(classes)==39 and fibers==748 and tables==1549 and options==96528
    above=[r for r in classes if r['maximum']>377]
    assert len(above)==1 and above[0]['counts']==[0,0,0,0,1,0,3,0,0]and above[0]['residual_rank']==0
    witnesses=[r['attaining_forest']for r in above[0]['rows']if r['maximum']>377]
    assert len(witnesses)==8 and all(value(fs)==378*64 for fs in witnesses)
    for fs in witnesses:
        co=Counter(fs);assert co['B',13]+co['D',13]==1 and co['A',15]+co['D',15]==3
        assert context(fs)==0 and all(actual_cross(a,b)==0 for i,a in enumerate(fs)for b in fs[i+1:])
    return dict(classes=classes,actual_fixed_fibers=fibers,fresh_actual_residual_tables=tables,explicit_options=options,
        Bellman_inequalities=inequalities,attaining_predecessors=predecessors,exact_maximum='378',
        classes_above_curvature=above,exact_obstructing_forests=witnesses,
        all_other_38_classes_at_most_curvature=True)

def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    cover=representative_cover(docs);inc=increasing_family(docs);finite=complement(docs,None if saved is None else saved['finite'])
    assert docs[P.TAG]['ideal_weights_changed']is False
    assert docs[S.TAG]['positive_ideal_profiles']==2552 and docs[S.TAG]['actual_cusp_surplus']=='36'
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_D17_A16_COUPLED_AMBIENT58_FRONTIER',sources={p:digest(p)for p in sorted(paths)},
        local_theorem=P.TAG,scale=64,representative_cover=cover,new_context_family=inc,finite=finite,
        exact_actual_finite_maximum='378',positive_ideal_profiles=2552,actual_cusp_surplus='36',
        finite_curvature='377',strict_global_constant='4205',necessary_ambient58_inequality='v_f >= 4205 + 36*v_i',
        minimum_finite_vertices=4205,completed_dimension_bound=58,new_dimension_bound=None,
        scope='Complete actual frontier under the five sealed exchanges and A13 guard. Fresh actual residual certificates, two complete increasing-family bounds, no new global representative DP or ideal convolution. Eight B13/D13 plus three rank15 forests remain.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();saved=json.loads(p.read_text())if args.verify else None
    out=run(saved)
    if args.verify:assert out==saved
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['exact_actual_finite_maximum'],len(out['finite']['exact_obstructing_forests']),flush=True)
