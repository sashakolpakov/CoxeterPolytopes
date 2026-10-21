#!/usr/bin/env python3
"""Sparse proposal: all inherited bonus coefficients are free, not retained.

Complete marked finite-forest budgets and inherited zero-corner quadrilateral
obligations; selected triangles only. New full affected-face coverage is open.
"""
import argparse,json,math
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
import fresh_sqrt5_marked_compensation_proposal_20261016 as P
import fresh_sqrt5_marked_h4_remainder_free_bonus_20261016 as B

ROOT=P.ROOT;TAG='fresh_sqrt5_marked_a4_d4_cross_free_bonus_20261016'
H,R=P.H,P.R;MODELS=P.MODELS;HUBS=H.HUBS;RESERVED=P.RESERVED+('A4','D4')
RANKS=(1,2,3,4,4,4)
LEAF_COUNTS={'A4':2,'D4':3}

@lru_cache(None)
def contexts(rank):
    out=[]
    for base in P.contexts(rank):
        rem=rank-sum(m*c for m,c in zip(RANKS,base))
        for a4 in range(rem//4+1):
            for d4 in range((rem-4*a4)//4+1):
                if rem-4*a4-4*d4!=1:out.append(base+(a4,d4))
    return tuple(out)

COMP_ORDER=['long_distance2_discount','long_adjacent3_discount'];REMAP={(0,None):0,(1,None):1}
H4_INDICES=set()
for j,channel in enumerate(P.CHANNELS,2):
    if channel['type']=='H4':
        H4_INDICES.add(j)
        for count in range(1,5):
            for a1presence in range(2):
                for marked_h4 in([None]if channel['facet_internal']else[False,True]):
                    for d_remainder in ([False,True]if count==4 else[None]):
                        REMAP[j,count,a1presence,marked_h4,d_remainder]=len(COMP_ORDER)
                        COMP_ORDER.append(dict(**channel,H4_count=count,A1_presence=a1presence,
                            facet_component_is_H4=True if channel['facet_internal']else marked_h4,
                            non_H4_remainder_contains_D4_or_D5=d_remainder))
    else:REMAP[j,None]=len(COMP_ORDER);COMP_ORDER.append(channel)
SMALL=len(COMP_ORDER);OLDCOUNT=len(P.RULE.x)-1
SAME={(a,c):SMALL+OLDCOUNT+i for i,(a,c)in enumerate((a,c)for a in range(18)for c in range(2)if a+4+3*c<=21)}
TWIN_KEYS=sorted({(name,ctx[0],min(ctx[1],2),min(ctx[2],1),ctx[RESERVED.index(name)])
    for ctx in contexts(21)for name in LEAF_COUNTS if ctx[RESERVED.index(name)]>=2})
TWIN={key:SMALL+OLDCOUNT+len(SAME)+i for i,key in enumerate(TWIN_KEYS)}
assert SMALL+OLDCOUNT+len(SAME)==B.N
N=SMALL+OLDCOUNT+len(SAME)+len(TWIN);TARGET=F(35)


def plus(v,w,m=1):
    out=v.copy()
    for j,z in w.items():
        out[j]=out.get(j,F(0))+m*z
        if not out[j]:del out[j]
    return out
def value(v,x):return v.get(-1,0)+sum(z*x[j]for j,z in v.items()if j>=0)
def key(v):return tuple(sorted(v.items()))
def small(v,count,a1presence,marked_h4,d_remainder):
    out={}
    for i,z in enumerate(v):
        if not z:continue
        j=i-1
        if j==-1:out[j]=z
        elif j in H4_INDICES:
            if count in range(1,5):
                channel=P.CHANNELS[j-2]
                out[REMAP[j,count,a1presence,None if channel['facet_internal']else marked_h4,bool(d_remainder)if count==4 else None]]=z
        else:out[REMAP[j,None]]=z
    return out


def bonus(kind,k,ctx):
    a,b,c=ctx[:3];j=P.RULE.index.get((kind,k,a+1,min(b,2),min(c,1)))
    return {}if j is None else{SMALL+j:F(1)}
def same(ctx):
    a,b,c=ctx[:3]
    return {SAME[a,min(c,1)]:F(1)}if b>=2 else{}


def twin(name,ctx):
    h=ctx[RESERVED.index(name)]
    if h<2:return {}
    return {TWIN[name,ctx[0],min(ctx[1],2),min(ctx[2],1),h]:F(1)}


@lru_cache(None)
def static(name,count,a1presence,marked_h4,d_remainder):
    total,mvs,connected,markkeys,rootkeys=P.static_component(name,count==5)
    return small(total,count,a1presence,marked_h4,d_remainder),{e:small(v,count,a1presence,marked_h4,d_remainder)for e,v in mvs.items()},connected,markkeys,rootkeys


@lru_cache(None)
def components(ctx,mh):
    a,b,c,k,a4,d4,d_remainder=ctx
    # All pre-existing coefficients retain their literal indices. Share the
    # frozen base profiles across the newly counted A4/D4 spectator states.
    oldmh=mh if mh in B.RESERVED else None
    basevals,basemarks=B.components((a,b,c,k,d_remainder),oldmh)
    vals=basevals.copy()
    marks={(name,e):v for(name,e),v in basemarks.items()
        if(name not in RESERVED if mh is None else name==mh)}
    for name,m in LEAF_COUNTS.items():
        h=ctx[RESERVED.index(name)]
        if not h:continue
        t=twin(name,ctx)
        vals[name]=plus(vals[name],t,F(1,h))
        for(typ,e),v in list(marks.items()):
            if typ==name:
                # Among m^2*C(h,2) eligible leaf pairs, a marked leaf
                # removes m*(h-1), i.e. fraction 2/(m*h) of the total.
                factor=F(1,h)*(1-F(2,m)if len(MODELS[name][e])==1 else 1)
                marks[typ,e]=plus(v,t,factor)
    return vals,marks


@lru_cache(None)
def numerical(ctx,mh):
    vals,marks=components(ctx,mh)
    def pack(v):return float(v.get(-1,0)),[(j,float(z))for j,z in v.items()if j>=0]
    return {n:pack(v)for n,v in vals.items()},{n:pack(v)for n,v in marks.items()}


def oracle(x,exact=False):
    fronts=[]
    states=[base+(flag,)for base in contexts(21)
        for flag in([0,1]if base==(0,0,0,4,0,0)else[int(base[3]==4 and base[5]>0)])]
    for ctx in states:
        a,b,c,k,a4,d4,d_remainder=ctx;remaining=21-sum(m*t for m,t in zip(RANKS,ctx[:6]))
        # D4 is now counted explicitly. With four H4s, an unreserved D
        # component can only be D5; the separate D4 case leaves no rank.
        def allowed(name):
            return k!=4 or ((name=='D5')==bool(d_remainder))
        for mh in [None]+[hub for hub,count in zip(RESERVED,ctx[:6])if count]:
            vals,marks=components(ctx,mh)
            if exact:costs={n:value(v,x)for n,v in vals.items()};mc={n:value(v,x)for n,v in marks.items()}
            else:
                nv,nm=numerical(ctx,mh)
                evaluate=lambda p:p[0]+sum(z*x[j]for j,z in p[1])
                costs={n:evaluate(v)for n,v in nv.items()};mc={n:evaluate(v)for n,v in nm.items()}
            dp={0:(F(0)if exact else 0.0,())}
            for rank in range(1,remaining+1):
                choices=[]
                for name,g in MODELS.items():
                    if name in RESERVED or not allowed(name)or rank-len(g)not in dp:continue
                    old,forest=dp[rank-len(g)];choices.append((old+costs[name],forest+(name,)))
                if choices:dp[rank]=max(choices,key=lambda z:z[0])
            const=sum(count*costs[name]for name,count in zip(RESERVED,ctx[:6]));best=None;winner=None
            for(name,e),mv in marks.items():
                if name not in RESERVED and not allowed(name):continue
                if name in RESERVED:
                    if not ctx[RESERVED.index(name)]:continue
                    rank=remaining;score=const-costs[name]+mc[name,e]
                else:rank=remaining-len(MODELS[name]);score=const+mc[name,e]
                if rank not in dp:continue
                score+=dp[rank][0]
                if best is None or score>best:
                    best=score;winner=dict(context=ctx,marked_hub=mh,marked_component=name,facet_root=e,residual_components=dp[rank][1])
            if best is None:continue
            p=marks[winner['marked_component'],winner['facet_root']]
            counts=Counter(dict(zip(RESERVED,ctx[:6])))
            if winner['marked_component']in RESERVED:counts[winner['marked_component']]-=1
            counts.update(winner['residual_components'])
            for name,count in counts.items():p=plus(p,vals[name],count)
            fronts.append(dict(**winner,maximum=best,profile=p))
    return fronts


def compile_corner(vertices,edges,pair,marks):
    a,b=pair;pcs=R.components(R.graph(vertices,edges));names=[R.type_of(g)for g in pcs]
    ctx=tuple(names.count(name)for name in RESERVED)
    lc=next(g for g in pcs if a in g);rc=next(g for g in pcs if b in g)
    kind,ky=None,None
    if b in lc:
        kind='connected';ky=R.tree_key(lc,pair);name=R.type_of(lc)
        vs={e:P.pair_vector(name,lc,a,b,e,0)for e in marks}
        vs5={e:P.pair_vector(name,lc,a,b,e,5)for e in marks}
    else:
        vs={e:(F(0),)*(P.N+1)for e in marks};vs5=vs;lt,rt=R.type_of(lc),R.type_of(rc)
        if lt==rt=='I2(5)':kind='same'
        if lt==rt and lt in LEAF_COUNTS and len(lc[a])==len(rc[b])==1:kind='twin_'+lt
        for kd,hub,excluded in [('a1cross','A1',('A1',)),('hubcross','I2(5)',('A1','I2(5)')),('h3cross','H3',HUBS)]:
            if lt==hub and rt not in excluded:kind,ky=kd,R.tree_key(rc,[b])
            elif rt==hub and lt not in excluded:kind,ky=kd,R.tree_key(lc,[a])
    mark_h4={e:e==-2 or any(e in g and R.type_of(g)=='H4'for g in pcs)for e in marks}
    return ctx,kind,ky,vs,vs5,mark_h4,any(n in('D4','D5')for n in names)


def add_bonus(compiled,ext):
    out={};tables=[]
    for local,kind,ky,vs,vs5,mark_h4,local_d in compiled:
        base=tuple(a+b for a,b in zip(local,ext[:6]))
        ctx=base+(int(local_d or ext[5] or ext[6])if base[3]==4 else 0,)
        if kind=='same':out=plus(out,same(ctx),F(1,2*ctx[1]*(ctx[1]-1)))
        elif kind and kind.startswith('twin_'):
            name=kind[5:];h=ctx[RESERVED.index(name)];m=LEAF_COUNTS[name]
            out=plus(out,twin(name,ctx),F(2,m*m*h*(h-1)))
        elif kind:
            mult={'connected':1,'a1cross':ctx[0],'hubcross':2*ctx[1],'h3cross':3*ctx[2]}[kind]
            assert mult>0;out=plus(out,bonus(kind,ky,ctx),F(1,mult))
        tables.append({e:small(v,ctx[3],int(ctx[0]>0),mark_h4[e],ctx[6])for e,v in(vs5 if ctx[3]==5 else vs).items()})
    return out,tables


def faces():
    obligations={};counts={};tri=H.read('fresh_sqrt5_triangle_weights_20261015');literal=[]
    for i,row in enumerate(tri['markings']):
        literal.append((f'core{i}',row['order'],row['white'],row['black'],[(a,b,m)for(a,b),m in zip(combinations(range(row['order']),2),row['labels'])]))
    for i in(46,47,51,55):
        row=tri['markings'][i]
        for tail in range(1,18):
            es=[(a,b,m)for(a,b),m in zip(combinations(range(5),2),row['labels'])]+[(0,5,3)]+[(u,u+1,3)for u in range(5,4+tail)]
            literal.append((f'core{i}_tail{tail}',5+tail,row['white']+list(range(5,5+tail)),row['black'],es))
    for tail in range(1,20):
        es=[(0,2,3),(1,2,10),(0,3,3)]+[(u,u+1,3)for u in range(3,2+tail)]
        literal.append((f'I210_tail{tail}',3+tail,list(range(3,3+tail)),[0,1,2],es))
    cases=0
    for label,n,white,black,es in literal:
        remaining=22-n;marks=white+([-1,-2]if remaining else[])
        compiled=[compile_corner(white+[b for b in black if b!=o],es,[b for b in black if b!=o],marks)for o in black]
        # A4/D4 counts are exact common spectators. Only a residual
        # D5 needs an additional flag, and only residual rank five can
        # affect a corner whose total H4 count is four.
        external=[]
        for base in contexts(remaining):
            residue=remaining-sum(m*c for m,c in zip(RANKS,base))
            external.extend(base+(flag,)for flag in([0,1]if residue==5 else[0]))
        for ext in external:
            common,tables=add_bonus(compiled,ext)
            valid_marks=white+([-1]if remaining-4*ext[3]>0 else[])+([-2]if ext[3]>0 else[])
            for mark in valid_marks:
                p=common
                for table in tables:p=plus(p,table[mark])
                obligations.setdefault((key(p),F(2)),dict(kind='triangle',family=label,external_counts=ext,facet_root=mark));cases+=1
    counts['triangle_families']=len(literal);counts['triangle_marked_contexts']=cases
    old=H.read('fresh_sqrt5_parameter_weight_stress_20261015');queue=H.read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    rows=[('prior',i,r['row'])for i,r in enumerate(old['faces'])if r['row']['nodes']<=23]
    rows += [('finite_queue',r['index'],r['row'])for r in queue['records']if r['result']['singular']and r['row']['nodes']<=23]
    cases=0
    for family,index,row in rows:
        es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        compiled=[compile_corner(row['white']+[a,b],es,[a,b],[None])for a,b in product((0,2),(1,3))]
        assert all(not any(v[3][None])for v in compiled)
        for ext in contexts(23-row['nodes']):
            p,_=add_bonus(compiled,ext+(0,))
            obligations.setdefault((key(p),F(1)),dict(kind='quad',family=family,index=index,external_counts=ext));cases+=1
    counts['quad_rows']=len(rows);counts['quad_common_contexts']=cases
    return [(dict(p),target,json.loads(json.dumps(w)))for(p,target),w in obligations.items()],counts


def literal_profile(front):
    """Recount actual root pairs; do not use marked-hub debit factors."""
    ctx=tuple(front['context']);counts=Counter(dict(zip(RESERVED,ctx[:6])))
    name=front['marked_component'];mark=front['facet_root']
    if name in RESERVED:counts[name]-=1
    assert all(n not in RESERVED for n in front['residual_components'])
    counts.update(front['residual_components'])
    forest=[(name,mark)]+[(n,None)for n,count in counts.items()for _ in range(count)]
    assert sum(len(MODELS[n])for n,_ in forest)==21
    assert tuple(sum(n==typ for n,_ in forest)for typ in RESERVED)==ctx[:6]
    assert ctx[6]==(int(any(n in('D4','D5')for n,_ in forest))if ctx[3]==4 else 0)
    answer={};roots=[]
    for name,e in forest:
        total,mvs,con,mc,rk=static(name,ctx[3],int(ctx[0]>0),front['marked_component']=='H4',ctx[6])
        answer=plus(answer,total if e is None else mvs[e])
        for ky,count in(con if e is None else mc[e]).items():answer=plus(answer,bonus('connected',ky,ctx),count)
        roots.append((name,[(rk[u],len(MODELS[name][u])==1)for u in MODELS[name]if u!=e]))
    for(lt,ls),(rt,rs)in combinations(roots,2):
        for lk,ll in ls:
            for rk,rl in rs:
                if lt==rt and lt in LEAF_COUNTS and ll and rl:
                    h=ctx[RESERVED.index(lt)];m=LEAF_COUNTS[lt]
                    answer=plus(answer,twin(lt,ctx),F(2,m*m*h*(h-1)))
                if lt==rt=='I2(5)':answer=plus(answer,same(ctx),F(1,2*ctx[1]*(ctx[1]-1)))
                for kind,hub,mult,excluded in [('a1cross','A1',ctx[0],('A1',)),('hubcross','I2(5)',2*ctx[1],('A1','I2(5)')),('h3cross','H3',3*ctx[2],HUBS)]:
                    if lt==hub and rt not in excluded:answer=plus(answer,bonus(kind,rk,ctx),F(1,mult))
                    elif rt==hub and lt not in excluded:answer=plus(answer,bonus(kind,lk,ctx),F(1,mult))
    return answer


def dual_certificate(rows,rhs,meta,lp):
    for denominator in(10**6,10**9,10**12):
        mu=[F(round(float(-v)*denominator),denominator)if v<-(1e-12 if m['kind']=='vertex'else 1e-5) else F(0)for v,m in zip(lp.ineqlin.marginals,meta)]
        mass=sum(v for v,m in zip(mu,meta)if m['kind']=='vertex')
        assert mass>0
        mu=[v/mass if m['kind']=='vertex'else v for v,m in zip(mu,meta)]
        capacity=Counter()
        for row,v,m in zip(rows,mu,meta):
            if v and m['kind']=='vertex':
                for j,z in row.items():
                    if 2<=j<N:capacity[j]+=v*z
        for i,(row,m)in enumerate(zip(rows,meta)):
            if m['kind']!='vertex'and any(j>=2 and z<0 and not capacity[j]for j,z in row.items()):mu[i]=F(0)
        demand=Counter()
        for row,v,m in zip(rows,mu,meta):
            if v and m['kind']!='vertex':
                for j,z in row.items():
                    if 2<=j<N:demand[j]-=v*z
        factor=min([F(1)]+[capacity[j]/z for j,z in demand.items()if z>0])
        mu=[v*factor if m['kind']!='vertex'else v for v,m in zip(mu,meta)]
        columns=Counter()
        for row,v in zip(rows,mu):
            if v:
                for j,z in row.items():columns[j]+=v*z
        upper=[max(F(0),-columns[j])for j in range(2)]
        columns[0]+=upper[0];columns[1]+=upper[1];columns[N]+=1
        assert min(columns.values(),default=F(0))>=0
        lower=-sum(v*b for v,b in zip(mu,rhs))-sum(upper)
        print('EXACT_DUAL',denominator,float(lower),flush=True)
        if lower>TARGET:
            return dict(kind='exact_restricted_model_obstruction',lower_bound=str(lower),discount_upper_multipliers=list(map(str,upper)),
                terms=[dict(multiplier=str(v),row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(b),metadata=m)for row,b,m,v in zip(rows,rhs,meta,mu)if v])
    raise AssertionError('Exact dual recovery unresolved')


def solve(obligations):
    rows=[{j:-z for j,z in p.items()if j>=0}for p,_,_ in obligations]
    rhs=[p.get(-1,F(0))-target for p,target,_ in obligations]
    meta=[w for _,_,w in obligations];cuts=set();steps=[]
    def cut(f):
        k=key(f['profile'])
        if k in cuts:return
        cuts.add(k);r={j:z for j,z in f['profile'].items()if j>=0};r[N]=F(-1)
        rows.append(r);rhs.append(-f['profile'].get(-1,F(0)))
        meta.append(dict(kind='vertex',**{k:v for k,v in f.items()if k not in('profile','maximum')}))
    for f in oracle([0.0]*N):cut(f)
    for iteration in range(150):
        rr=[];cc=[];vv=[]
        for i,row in enumerate(rows):
            for j,z in row.items():rr.append(i);cc.append(j);vv.append(float(z))
        A=coo_matrix((vv,(rr,cc)),shape=(len(rows),N+1)).tocsr();objective=np.zeros(N+1);objective[N]=1
        lp=linprog(objective,A_ub=A,b_ub=np.array(list(map(float,rhs))),bounds=[(0,1),(0,1)]+[(0,None)]*(N-1),method='highs')
        assert lp.success,lp.message
        print('LP',iteration,'lower',lp.fun,'cuts',len(cuts),'discounts',lp.x[:2],'upper_duals',lp.upper.marginals[:2],flush=True)
        steps.append(dict(iteration=iteration,lower=float(lp.fun),cuts=len(cuts),discounts=list(map(float,lp.x[:2])),discount_upper_duals=list(map(float,lp.upper.marginals[:2]))))
        if lp.fun>35+1e-7:return dual_certificate(rows,rhs,meta,lp),steps
        fronts=oracle(lp.x[:N]);maximum=max(f['maximum']for f in fronts)
        print('ORACLE',maximum,flush=True)
        if maximum<=35+1e-7:
            for den in(10**6,10**8,10**10):
                x=[F(max(0,math.floor(float(v)*den-1e-7)),den)if j<2 else F(math.ceil(float(v)*den+1e-7),den)if v>1e-10 else F(0)for j,v in enumerate(lp.x[:N])]
                factor=F(1)
                for p,target,_ in obligations:
                    required=target-p.get(-1,F(0))-sum(p.get(j,0)*x[j]for j in range(2))
                    paid=sum(z*x[j]for j,z in p.items()if j>=2)
                    if required>paid:
                        assert paid>0;factor=max(factor,required/paid)
                x=[v if j<2 else v*factor for j,v in enumerate(x)]
                assert all(value(p,x)>=target for p,target,_ in obligations)
                exact=oracle(x,True);upper=max(f['maximum']for f in exact)
                print('EXACT_UPPER',den,float(upper),flush=True)
                if upper<=35:return dict(kind='exact_feasible_restricted_rule',primal=list(map(str,x)),maximum=str(upper)),steps
        old=len(cuts)
        for f in fronts:
            if f['maximum']>lp.fun+1e-7:cut(f)
        assert len(cuts)>old,'Numerical convergence but exact repair unresolved'
    raise AssertionError('Iteration cap')


def verify(obligations,cert):
    if cert['kind']=='exact_feasible_restricted_rule':
        x=list(map(F,cert['primal']));assert len(x)==N and min(x)>=0 and max(x[:2])<=1
        assert all(value(p,x)>=t for p,t,_ in obligations)
        fronts=oracle(x,True);maximum=max(f['maximum']for f in fronts);assert maximum==F(cert['maximum'])<=35
        assert all(literal_profile(f)==f['profile']for f in fronts)
        return dict(maximum=str(maximum),marked_contexts=len(fronts),fronts=[dict(**{k:v for k,v in f.items()if k not in('profile','maximum')},maximum=str(f['maximum']))for f in fronts])
    lookup={(key({j:-z for j,z in p.items()if j>=0}),p.get(-1,F(0))-t):w for p,t,w in obligations}
    columns=Counter();lower=F(0)
    for term in cert['terms']:
        row={j:F(z)for j,z in term['row']};rhs=F(term['rhs']);mu=F(term['multiplier']);m=term['metadata'];assert mu>=0
        if m['kind']=='vertex':
            ctx=tuple(m['context']);vals,marks=components(ctx,m['marked_hub']);p=marks[m['marked_component'],m['facet_root']]
            counts=Counter(dict(zip(RESERVED,ctx[:6])))
            if m['marked_component']in RESERVED:counts[m['marked_component']]-=1
            counts.update(m['residual_components'])
            for name,count in counts.items():p=plus(p,vals[name],count)
            assert p==literal_profile(m)
            wanted={j:z for j,z in p.items()if j>=0};wanted[N]=F(-1)
            assert row==wanted and rhs==-p.get(-1,F(0))
        else:assert lookup[key(row),rhs]==m
        for j,z in row.items():columns[j]+=mu*z
        lower-=mu*rhs
    upper=list(map(F,cert['discount_upper_multipliers']));assert min(upper)>=0
    for j,z in enumerate(upper):columns[j]+=z;lower-=z
    columns[N]+=1;assert min(columns.values(),default=F(0))>=0 and lower==F(cert['lower_bound'])>35
    return dict(lower_bound=str(lower),terms=len(cert['terms']),exact_weak_duality=True)


def normalization_audit():
    records=[]
    for name,m in LEAF_COUNTS.items():
        g=MODELS[name];leaves=[u for u in g if len(g[u])==1]
        assert len(leaves)==m
        for h in range(2,6):
            pairs=[((i,u),(j,v))for i,j in combinations(range(h),2)for u in leaves for v in leaves]
            denominator=m*m*h*(h-1)//2
            assert len(pairs)==denominator
            for e in g:
                retained=sum((0,e)not in p for p in pairs)
                expected=1-F(2,m*h)if e in leaves else F(1)
                assert F(retained,denominator)==expected
                allocated=F(h-1,h)+F(1,h)*(1-F(2,m)if e in leaves else 1)
                assert allocated==expected
                records.append(dict(component=name,copies=h,facet_root=e,eligible_facet_root=e in leaves,
                    eligible_pairs=denominator,retained_pairs=retained,remaining_gamma_coefficient=str(expected)))
    return records


def diagnostic_witnesses():
    # An explicit local release of the previous two-term barrier. This is
    # deliberately not a feasible assignment for the entire new LP.
    av={0:F(0),1:F(0),2:F(1),3:F(2),4:F(4),5:F(5)}
    dv={0:F(0),1:F(0),2:F(9,8),3:F(9,4),4:F(9,2),5:F(5)}
    face_values=[]
    for a in range(4):
        for d in range(4-a):
            h,k=a+2,d+2
            paid=F(3,4)+F(2)*av[h]/(4*h*(h-1))+F(2)*dv[k]/(9*k*(k-1))
            assert paid>=1
            face_values.append(dict(external_A4=a,external_D4=d,payment=str(paid)))
    vertices=[dict(A4_count=a,D4_count=5-a,payment=str(30+av[a]+dv[5-a]))for a in range(6)]
    assert all(F(v['payment'])<=35 for v in vertices)and 32+8*F(3,8)==35
    old=H.read('fresh_sqrt5_parameter_weight_stress_20261015');faces=[]
    for index in(172,1260,1246):
        row=old['faces'][index]['row']
        es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        corners=[]
        for a,b in product((0,2),(1,3)):
            pcs=R.components(R.graph(row['white']+[a,b],es));left=next(g for g in pcs if a in g);right=next(g for g in pcs if b in g)
            corners.append(dict(pair=[a,b],components=sorted(R.type_of(g)for g in pcs),same_component=b in left,
                left_component=R.type_of(left),right_component=R.type_of(right),
                left_root_degree=len(left[a]),right_root_degree=len(right[b]),
                left_root_key=R.tree_key(left,[a]),right_root_key=R.tree_key(right,[b])))
        faces.append(dict(prior_index=index,literal_row=row,corners=corners))
    return dict(prior64_local_release=dict(connected_D5_weight='3/8',A4_gamma={h:str(v)for h,v in av.items()},
        D4_gamma={h:str(v)for h,v in dv.items()},common_context='Actual A1 count1, I25/H3 count0; all feasible common A4/D4 spectator counts.',
        face_values=face_values,marked_A1_four_D5_payment='35',marked_A1_five_A4_D4_vertices=vertices,
        scope='Only prior64 and the marked-A1 rank21 forests built from A4,D4,D5. This is not a primal certificate for the full model.'),
        dominant_literal_faces=faces)


def run(saved=None):
    obligations,counts=faces();print('MODEL',N,'variables',len(obligations),'faceprofiles',counts,flush=True)
    cert,steps=solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=verify(obligations,cert)
    paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',f'results/{P.TAG}.json',f'scripts/{B.TAG}.py',f'results/{B.TAG}.json',f'scripts/{H.TAG}.py',f'results/{H.TAG}.json',
        f'results/{H.PRIMARY}_cap2.json','results/fresh_sqrt5_triangle_weights_20261015.json',
        'results/fresh_sqrt5_parameter_weight_stress_20261015.json','results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json']
    variable_order=COMP_ORDER[:]
    inherited_kinds={j:k for(k,orbit,a,b,c),j in P.RULE.index.items()}
    variable_order += [dict(kind='inherited_free_bonus',bonus_kind=inherited_kinds[j],**v,actual_A1_count=v['A1_count']-1)for j,v in enumerate(P.PRIMARY['variable_order'])]
    variable_order += [dict(kind='same_i25',actual_A1_count=a,I25_count_at_least=2,H3_presence=c)for a,c in SAME]
    variable_order += [dict(kind='same_component_leaf_cross',component=name,actual_A1_count=a,I25_context=b,H3_context=c,actual_component_count=h,eligible_roots_per_component=LEAF_COUNTS[name],pair_normalization=str(F(2,LEAF_COUNTS[name]**2*h*(h-1))))for name,a,b,c,h in TWIN]
    assert len(variable_order)==N
    return dict(status='PASS_EXACT_A4_D4_CROSS_FREE_BONUS_MARKED_PROPOSAL',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        dimension=21,facet_dimension=20,target='35',variables=N,variable_order=variable_order,face_counts=counts,
        unique_face_profiles=len(obligations),certificate=cert,verification=proof,steps=steps,
        new_cross_literal_normalization=normalization_audit(),
        diagnostic_witnesses=diagnostic_witnesses(),
        H4_compensation_context='Separate coefficients for actual H4 count1,2,3,4, A1 absence/presence, and facet component H4 context. At H4 count4 also distinguish whether the non-H4 remainder contains D4 or D5; the only positive states are D5 and A1+D4. Compensation stays zero at count5. Common external facet marks distinguish H4/non-H4 components; external D4/D5 presence is common to all corners.',
        vertex_context_order=['A1_count','I25_count','H3_count','H4_count','A4_count','D4_count','D4_or_D5_present_if_H4_count4'],
        triangle_external_context_order=['A1_count','I25_count','H3_count','H4_count','A4_count','D4_count','D5_presence_if_unreserved_residual_rank5'],
        quadrilateral_external_context_order=['A1_count','I25_count','H3_count','H4_count','A4_count','D4_count'],
        new_cross_normalization='Only A4 endpoint pairs across distinct A4 components, and D4 leaf pairs across distinct D4 components. For m eligible roots and h copies, each pair pays2*gamma/(m*m*h*(h-1)); the unmarked forest total is gamma. A marked eligible root removes2*gamma/(m*h). Coefficients distinguish actual h, exact A1 count, capped I25 count0/1/2+, and H3 absence/presence.',
        scope='All inherited bonus coefficients are independent nonnegative variables; no old primal weights are fixed. Complete marked finite-forest DP and all inherited baseline-zero quadrilateral obligations are retained. Only selected triangle families are covered. Newly affected triangles and paid-corner quadrilaterals still require complete coverage. No dimension-bound promotion.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    path=ROOT/f'results/{TAG}.json';old=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(old),default=str))
    if args.verify:assert result==old
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
