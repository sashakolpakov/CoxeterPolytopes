#!/usr/bin/env python3
"""Complete ambient57 diagnostic for the frozen dimension57 payment rule.

All final dimension58-only guards are vacuous at rank57. Fresh representative
and actual residual certificates retain every other coupled contextual cost.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement,product
from math import prod
import fresh_d17_a16_coupled_ambient58_frontier_20261018 as OLD
import fresh_b13_endpoint_retention_20261018 as Z
import fresh_d13_context_dimension57_20261018 as END

ROOT=OLD.ROOT;C=OLD.C;P=OLD.P;N=57;D=56;SCALE=64
TAG='fresh_ambient57_payment_frontier_20261018';BETA=F(364)
STEMS=(END.TAG,Z.TAG,OLD.TAG,END.IDEAL,C.TAG,P.TAG)
FINAL_ENDPOINT_MD='3b4c0ee9b8ba8163eddd996c18872b92e037f025f36f93e6df2b8619ecb0a144'
enc=OLD.enc;representative=C.representative;fixed_components=C.fixed_components
TYPES=sorted(t for t in C.ENGINE.types()if t[1]<=N)
ITEMS=tuple(r for r in C.ITEMS if r<=N)
@lru_cache(None)
def actual_base(t):return P.actual_base(t)-(64 if t==('B',13)else 0)
actual_cross=lru_cache(None)(P.actual_cross);context=P.context

def value(fs):
    assert sum(n for k,n in fs)<=N
    return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def fixed_domain():
    out=[]
    def visit(r,left,cs):
        if r==18:
            c=dict(zip(range(9,18),cs));key=(2*c[10]+3*sum(c[s]for s in range(11,17))+2*c[17],c[11]+c[12],c[13],c[11],c[14],c[9],c[10],c[16],c[15])
            for a12 in range(c[12]+1):out.append((tuple(cs),N-left,key+(a12,)))
            return
        for k in range(left//r+1):visit(r+1,left-r*k,cs+[k])
    visit(9,N,[]);out.sort()
    # A second domain construction by bounded Cartesian counts.
    other=[]
    for cs in product(*(range(N//r+1)for r in range(9,18))):
        used=sum(r*k for r,k in zip(range(9,18),cs))
        if used<=N:
            for a12 in range(cs[3]+1):other.append((cs,used,a12))
    assert [(cs,u,key[9])for cs,u,key in out]==sorted(other)
    assert len(out)==len({key for cs,u,key in out})
    return out

def item(r,key,k5):return C.item(r,key,k5)+(40*k5 if r==1 and key[5]else 0)
def upper_value(fs):
    co=Counter(fs)
    return C.value(fs,True)+22*co['A',17]+(40*co['A',1]*co['A',5]if co['B',9]else 0)
def additive(ts,ws,cap,saved=None):
    if saved is None:
        v=[0]+[None]*cap
        for t,w in zip(ts,ws):
            for n in range(t[1],cap+1):
                if v[n-t[1]]is not None:
                    z=v[n-t[1]]+w
                    if v[n]is None or z>v[n]:v[n]=z
        last=[None]+[next((i for i,(t,w)in enumerate(zip(ts,ws))if t[1]<=n and v[n-t[1]]is not None and v[n]==v[n-t[1]]+w),None)for n in range(1,cap+1)]
    else:
        assert saved['types']==enc(ts)and saved['weights']==ws
        v,last=saved['values'],saved['last']
    assert len(v)==len(last)==cap+1 and v[0]==0 and last[0]is None
    checks=attain=0
    for n in range(1,cap+1):
        possible=[]
        for i,(t,w)in enumerate(zip(ts,ws)):
            if t[1]<=n and v[n-t[1]]is not None:
                assert v[n]is not None and v[n]>=v[n-t[1]]+w;possible.append(i);checks+=1
        if possible:
            i=last[n];assert i in possible and v[n]==v[n-ts[i][1]]+ws[i];attain+=1
        else:assert v[n]is None and last[n]is None
    return dict(types=ts,weights=ws,values=v,last=last),checks,attain

def coarse(saved=None):
    guards={('A',13):64,('A',14):32,('D',14):64,('A',15):32,('D',15):32,('A',16):80,('D',16):64}
    for a in TYPES:
        assert actual_base(a)+guards.get(a,0)<=C.upper_base(representative(a))+(22 if a==('A',17)else 0)
        for b in TYPES:
            extra=24 if frozenset((a,b))==frozenset((('A',1),('D',5)))else 32 if frozenset((a,b))==frozenset((('A',4),('A',12)))else 0
            assert actual_cross(a,b)+extra<=C.upper_cross(representative(a),representative(b))
    dom=fixed_domain();tabs=[];rows=[];quarter=[];checks=attain=opts=qopts=0
    old={(tuple(t['key']),t['k5'],t['exclude2']):t for t in saved['tables']}if saved else{}
    expected=set()
    residual=[representative(('A'if r<=8 else'D',r))for r in ITEMS]
    assert all(C.upper_cross(a,b)==0 for a in residual for b in residual)
    for cs,used,key in dom:
        large=fixed_components(cs,key)
        for k5 in range((N-used)//5+1):
            fixed=large+[('A',5)]*k5;rem=N-used-5*k5;tables={}
            for r,t in zip(ITEMS,residual):
                assert item(r,key,k5)==C.upper_base(t)+sum(C.upper_cross(t,b)for b in fixed)+(40*k5 if r==1 and cs[0]else 0)
            for no2 in(False,True):
                kk=(key,k5,no2);expected.add(kk);ts=[t for r,t in zip(ITEMS,residual)if not no2 or r!=2]
                ws=[item(t[1],key,k5)for t in ts];tab,a,b=additive(ts,ws,N,old.get(kk))
                checks+=a;attain+=b;tab.update(key=key,k5=k5,exclude2=no2);tabs.append(tab);tables[no2]=tab
            normal,no2=tables[False],tables[True];fc=upper_value(fixed)
            options=[(normal['values'][rem],0,0,rem)]
            for reserve in(1,2):
                oldguard=reserve==2 and bool(cs[1]);newguard=bool(cs[0]and k5)
                if rem<reserve or not(oldguard or newguard)or(reserve==2 and not oldguard):continue
                for a2 in range((rem-reserve)//2+1):
                    left=rem-reserve-2*a2
                    bonus=(36*a2*(a2-1)if oldguard else 0)+(64*a2*k5 if newguard else 0)
                    options.append((reserve*item(1,key,k5)+a2*item(2,key,k5)+no2['values'][left]+bonus,reserve,a2,left))
            opts+=len(options);z,a1,a2,left=max(options);fs=fixed+[('A',1)]*a1+[('A',2)]*a2;tab=no2 if a1 else normal
            while left:
                t=tab['types'][tab['last'][left]];fs.append(t);left-=t[1]
            assert sum(t[1]for t in fs)==N and upper_value(fs)==fc+z
            row=dict(counts=cs,key=key,k5=k5,a12=key[9],residual_rank=rem,upper=F(fc+z,64),upper_forest=fs)
            rows.append(row)
            if cs[1]and rem>=7:
                best=None
                for a1 in range(1,rem-5):
                    for a4 in range(1,(rem-a1-2)//4+1):
                        for a2 in range(1,(rem-a1-4*a4)//2+1):
                            left=rem-a1-4*a4-2*a2
                            bonus=(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k5 if cs[0]else 0)+64*a2*a4
                            z=fc+a1*item(1,key,k5)+a2*item(2,key,k5)+a4*item(4,key,k5)+no2['values'][left]+bonus;qopts+=1
                            if best is None or z>best[0]:best=(z,a1,a2,a4,left)
                assert best is not None;quarter.append(dict(counts=cs,key=key,k5=k5,a12=key[9],residual_rank=rem,upper=F(best[0],64),explicit_counts=best[1:]))
    if saved:assert set(old)==expected and len(old)==len(saved['tables'])
    costly={(tuple(r['counts']),r['k5'],r['a12'])for r in rows+quarter if F(r['upper'])>BETA}
    return dict(fixed_large_vectors=len(dom),expanded_classes=len(rows),tables=tabs,rows=rows,quarter_rows=quarter,
        actual_type_pairs=len(TYPES)**2,Bellman_inequalities=checks,attaining_predecessors=attain,
        ordinary_context_options=opts,quarter_context_options=qopts,complete_costly_keys=sorted(costly),
        maximum_coarse_upper=max(r['upper']for r in rows+quarter),scope='Complete representative complement and explicit A1/A2/A4 increasing family; upper relaxations need not be attained.')

def actual_refinement(coarse,saved=None):
    fibers={}
    for t in TYPES:fibers.setdefault(representative(t),[]).append(t)
    reps={representative(('A'if r<=8 else'D',r))for r in ITEMS}
    residual=[t for t in TYPES if representative(t)in reps and t not in[('A',1),('A',2)]]
    assert all(actual_cross(a,b)==0 for a in residual+[('A',1),('A',2)]for b in residual+[('A',1),('A',2)])
    keyed={(tuple(r['counts']),r['k5'],r['a12']):r for r in coarse['rows']}
    old={(tuple(r['counts']),r['k5'],r['a12']):r for r in saved['classes']}if saved else{}
    classes=[];nrows=ntab=nopt=checks=attain=0
    for cs,k5,a12 in coarse['complete_costly_keys']:
        key=(tuple(cs),k5,a12);base=keyed[key];groups=sorted(Counter(fixed_components(cs,base['key'])+[('A',5)]*k5).items());rem=base['residual_rank']
        fixedsets=sorted(tuple(sorted(t for group in choice for t in group))for choice in product(*(combinations_with_replacement(fibers[t],n)for t,n in groups)))
        assert len(fixedsets)==len(set(fixedsets))
        priors={tuple(map(tuple,r['fixed_actual'])):r for r in old[key]['rows']}if saved else{}
        if saved:assert set(priors)==set(fixedsets)
        rows=[]
        for fixed in fixedsets:
            co=Counter(fixed);fc=value(fixed);w1=actual_base(('A',1))+sum(actual_cross(('A',1),t)for t in fixed);w2=actual_base(('A',2))+sum(actual_cross(('A',2),t)for t in fixed)
            cache={};certs=[];context_checked=set();prior=priors.get(fixed)
            def table(has1,a2,no3):
                nonlocal ntab,checks,attain
                a4bonus=(32*co['A',12]+64*a2*bool(co['A',10]))if has1 else 0
                kk=(a4bonus,no3)
                if kk not in cache:
                    ts=[t for t in residual if not(no3 and t==('A',3))]
                    ws=[actual_base(t)+sum(actual_cross(t,b)for b in fixed)+(a4bonus if t==('A',4)else 0)for t in ts]
                    p=next((p for p in prior['tables']if p['key']==list(kk)),None)if prior else None
                    if prior:assert p is not None
                    tab,a,b=additive(ts,ws,rem,p);checks+=a;attain+=b;ntab+=1;tab['key']=kk;certs.append(tab);cache[kk]=tab
                return cache[kk]
            best=None
            for a1 in range(rem+1):
                for a2 in range((rem-a1)//2+1):
                    fs=list(fixed)+[('A',1)]*a1+[('A',2)]*a2;left=rem-a1-2*a2
                    modes=[(False,0)]if a1 or a2 or not co['D',14]else[(True,0),(False,1)]
                    for no3,res3 in modes:
                        rest=left-3*res3
                        if rest<0:continue
                        tab=table(bool(a1),a2,no3);v=tab['values']
                        if v[rest]is None:continue
                        chosen=fs+[('A',3)]*res3
                        # All remaining contexts are linear in A4. A3 presence
                        # has been fixed by the reserved/no-A3 partition.
                        bonus=context(chosen)
                        ck=(bool(a1),a2,no3,res3)
                        if ck not in context_checked:
                            for t,w in zip(tab['types'],tab['weights']):
                                if t[1]>rest:continue
                                expected=(32*co['A',12]+64*a2*bool(co['A',10]))if a1 and t==('A',4)else 0
                                assert context(chosen+[t])-bonus==expected,(fixed,a1,a2,no3,res3,t)
                            context_checked.add(ck)
                        extra=0 if not res3 else actual_base(('A',3))+sum(actual_cross(('A',3),t)for t in fixed)
                        z=fc+a1*w1+a2*w2+extra+bonus-context(fixed)+v[rest];nopt+=1
                        if best is None or z>best[0]:best=(z,a1,a2,no3,res3,rest)
            assert best is not None
            z,a1,a2,no3,res3,left=best;fs=list(fixed)+[('A',1)]*a1+[('A',2)]*a2+[('A',3)]*res3;tab=table(bool(a1),a2,no3)
            while left:
                t=tab['types'][tab['last'][left]];fs.append(t);left-=t[1]
            assert sum(t[1]for t in fs)==N and value(fs)==z
            if prior:assert len(certs)==len(prior['tables'])
            rows.append(dict(fixed_actual=fixed,tables=certs,maximum=F(z,64),attaining_forest=sorted(fs),winning_option=best[1:]));nrows+=1
        classes.append(dict(counts=cs,k5=k5,a12=a12,residual_rank=rem,rows=rows,maximum=max(r['maximum']for r in rows)))
    if saved:assert set(old)=={(tuple(r['counts']),r['k5'],r['a12'])for r in classes}
    maximum=max(r['maximum']for r in classes);above=[r for r in classes if r['maximum']>BETA]
    assert maximum>BETA
    return dict(classes=classes,complete_costly_classes=len(classes),actual_fixed_fibers=nrows,fresh_actual_tables=ntab,
        Bellman_inequalities=checks,attaining_predecessors=attain,explicit_options=nopt,exact_actual_maximum=maximum,
        above_curvature_classes=[dict(counts=r['counts'],k5=r['k5'],a12=r['a12'],maximum=r['maximum'],
            maximizers=[q['attaining_forest']for q in r['rows']if q['maximum']==r['maximum']])for r in above],
        all_other_classes_at_most_curvature=True,all_obstructing_forests_enumerated=False)

@lru_cache(None)
def high_parts(n,least=8):
    if n==0:return((),)
    return tuple((r,)+p for r in range(least,n+1)for p in high_parts(n-r,r))
def ideal(previous):
    env={int(k):F(v)for k,v in previous['current_envelopes'].items()};domain=sorted((low,)+p if low else p for low in range(8)for p in high_parts(D-low))
    second=[]
    def counts(r,left,out):
        if r>D:
            if left<8:second.append(tuple(([left]if left else[])+out))
            return
        for k in range(left//r+1):counts(r+1,left-r*k,out+[r]*k)
    counts(8,D,[]);assert domain==sorted(second)and len(domain)==len(set(domain))
    rows=[]
    for rs in domain:
        h=END.polynomial(rs);v=prod(r+1 for r in rs);assert sum(rs)==D and len(h)==D+1 and h==h[::-1]and sum(h)==v
        g=sum(min(j,D-j)*a for j,a in enumerate(h));assert 2*g==D*v-sum(abs(2*j-D)*a for j,a in enumerate(h))
        curv=F(5*g-D*v,2);cost=sum(env[r]for r in rs);margin=curv/v-cost
        rows.append(dict(ranks=rs,h_polynomial=h,link_vertices=v,G=g,curvature=curv,normalized_cost=cost,normalized_margin=margin,absolute_margin=v*margin))
    negative=[r for r in rows if r['normalized_margin']<0];zero=[r for r in rows if r['normalized_margin']==0]
    return dict(current_envelopes=env,complete_profiles=len(rows),rows=rows,negative_profiles=negative,zero_profiles=zero,
        normalized_minimum=min(r['normalized_margin']for r in rows),absolute_minimum=min(r['absolute_margin']for r in rows),
        low8_terminal_domain=True,two_domain_reconstructions=True,fresh_convolutions=True,
        no_negative_surplus_extension=True,scope='Complete terminal rank56 scalar diagnostic; negative rows are not multiplied through concentration or asserted realized by cusps.')

def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    assert digest(f'results/{END.TAG}.md')==FINAL_ENDPOINT_MD
    for s,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    assert docs[END.TAG]['dimension_bound']==57
    assert 13+3*15>N and 28+2*15>N and 12+16+2*15>N
    ff=coarse(None if saved is None else saved['finite_upper']);print('COARSE',ff['expanded_classes'],len(ff['complete_costly_keys']),ff['maximum_coarse_upper'],flush=True)
    aa=actual_refinement(ff,None if saved is None else saved['actual_finite']);print('ACTUAL',aa['exact_actual_maximum'],len(aa['above_curvature_classes']),flush=True)
    ii=ideal(docs[END.IDEAL]['ideal'])
    g0=sum(min(j,D-j)for j in range(D+1));beta=F(5*g0-D*(D+1),2);ext=sum(2*j-D for j in range(D+1)if 2*j>D)
    assert(g0,beta,ext,5*ext)==(784,BETA,812,4060)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_COMPLETE_AMBIENT57_PAYMENT_FRONTIER_DIAGNOSTIC',sources={p:digest(p)for p in sorted(paths)},
        ambient_dimension=N,link_dimension=D,completed_dimension_bound=57,new_dimension_bound=None,scale=64,
        finite_curvature=beta,strict_global_constant=5*ext,dimension58_only_D13_guards_vacuous=True,
        finite_upper=ff,actual_finite=aa,ideal=ii,
        formalization_simplification='After fixing all large/rank5 types and A1/A2 counts, contexts are linear in residual A4; only D14/A3 needs a presence split. This is a cost normal form, not a replacement for coupled geometric output proofs.',
        scope='Fresh complete changed-dimension finite and terminal-ideal diagnostic; no new dimension56 theorem or realization claimed.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    saved=json.loads(path.read_text())if args.verify else None;out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_finite']['exact_actual_maximum'],len(out['actual_finite']['above_curvature_classes']),out['ideal']['complete_profiles'],len(out['ideal']['negative_profiles']),flush=True)
