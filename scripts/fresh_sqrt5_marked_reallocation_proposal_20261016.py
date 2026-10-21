#!/usr/bin/env python3
"""Counterexample-guided proposal only; affected-face coverage remains open.

The old n22 rule is indexed with one extra common A1, but cross denominators
are actual rank21 multiplicities. This preserves the inherited zero-corner
quadrilateral obligations. Discounts invalidate other old face pruning.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import fresh_sqrt5_dimension22_h3_cross_20261016_review as H
import fresh_sqrt5_marked_discount_census_20261016 as D

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_marked_reallocation_proposal_20261016'
R=H.R
MODELS={n:g for n,g in R.MODELS.items() if len(g)<=21}
RESERVED=('A1','I2(5)','H3','H4')
PRIMARY=H.read(H.PRIMARY+'_cap2')
RULE=H.Rule(PRIMARY)
CHANNELS=[];CHANNEL_INDEX={}
for name in ('H4','H3','I2(5)','I2(10)'):
    g=MODELS[name]
    seen=set()
    for a,b in combinations(g,2):
        key=R.tree_key(g,[a,b])
        if key in seen:continue
        seen.add(key)
        for internal in ([False,True] if len(g)>2 else [False]):
            CHANNEL_INDEX[name,key,internal]=len(CHANNELS)+2
            CHANNELS.append(dict(type=name,pair=[a,b],facet_internal=internal))
COMPENSATION_N=2+len(CHANNELS)
SCALES=[dict(kind='hubcross',type='A5',context=(7,2,0)),
        dict(kind='a1cross',type='I2(5)',context=(7,2,0)),
        dict(kind='connected',type='A5',context=(7,2,0))]
SCALE_KEYS={}
for index,rec in enumerate(SCALES):
    g=MODELS[rec['type']]
    marks=combinations(g,2)if rec['kind']=='connected'else([u]for u in g)
    for mark in marks:SCALE_KEYS[rec['kind'],R.tree_key(g,mark),rec['context']]=COMPENSATION_N+index
ADDITIONS=[];ADDITION_KEYS={}
OLD=H.read('fresh_sqrt5_parameter_weight_stress_20261015')
for rowindex,cornerindex,ctx in [(10,0,(7,1,0)),(10,2,(5,1,0)),(10,3,(5,2,0)),
          (80,0,(7,1,0)),(80,1,(6,2,0)),(94,0,(8,1,1)),(94,2,(7,2,1)),(94,3,(6,3,0))]:
    corner=H.corners(OLD['faces'][rowindex]['row'])[cornerindex]
    kind=next(k for k in H.KINDS if corner[k] in RULE.orbits[k]);key=corner[kind]
    index=COMPENSATION_N+len(SCALES)+len(ADDITIONS)
    assert(kind,key,ctx)not in ADDITION_KEYS
    ADDITION_KEYS[kind,key,ctx]=index
    ADDITIONS.append(dict(kind=kind,orbit=RULE.orbits[kind][key]['label'],context=ctx,
                          prior_row=rowindex,corner_index=cornerindex))
SAME_INDEX=COMPENSATION_N+len(SCALES)+len(ADDITIONS)
ADDITIONS.append(dict(kind='same_i25',context=(6,3,0),prior_row=80,corner_index=3,
    denominator='2*b*(b-1)',marked_I25_factor='(b-1)/b'))
N=SAME_INDEX+1


def contexts(rank):
    return [(a,b,c,k) for a in range(rank+1) for b in range((rank-a)//2+1)
            for c in range((rank-a-2*b)//3+1) for k in range((rank-a-2*b-3*c)//4+1)
            if rank-a-2*b-3*c-4*k!=1]


def add(v,w,scale=1):
    if scale==0:return v
    return tuple(a if not b else a+b if scale==1 else a+scale*b for a,b in zip(v,w))
ZERO=(F(0),)* (N+1)


def pair_vector(name,g,a,b,mark,h4count):
    vec=list(ZERO);vec[0]=R.baseline(g,a,b)
    vec[1]=-int(D.eligible('long_local33_distance2',name,g,a,b))
    vec[2]=-int(D.eligible('long_adjacent3',name,g,a,b))
    if name in ('H4','H3','I2(5)','I2(10)') and not(name=='H4' and h4count==5):
        vec[1+CHANNEL_INDEX[name,R.tree_key(g,[a,b]),mark in g]]=1
    return tuple(vec)


def oldbonus(kind,key,ctx):
    a,b,c,_=ctx
    return RULE.bonus(kind,key,a+1,b,c)


def bonusvector(kind,key,ctx):
    vec=list(ZERO);value=oldbonus(kind,key,ctx);index=SCALE_KEYS.get((kind,key,ctx[:3]))
    if index is None:vec[0]=value
    else:vec[1+index]=value
    index=ADDITION_KEYS.get((kind,key,ctx[:3]))
    if index is not None:vec[1+index]+=1
    return tuple(vec)


def samevector(ctx):
    vec=list(ZERO)
    if ctx[:3]==(6,3,0):vec[1+SAME_INDEX]=1
    return tuple(vec)


@lru_cache(None)
def static_component(name,five):
    g=MODELS[name];roots={}
    for e in g:roots.setdefault(R.tree_key(g,[e]),e)
    pv={(u,v):pair_vector(name,g,u,v,None,5 if five else 0)for u,v in combinations(g,2)}
    keys={(u,v):R.tree_key(g,[u,v])for u,v in combinations(g,2)}
    connected=Counter(key for key in keys.values()if key in RULE.orbits['connected'])
    total=ZERO
    for p in pv.values():total=add(total,p)
    marked={};markkeys={}
    for e in roots.values():
        mv=ZERO
        for u,v in pv:
            if e not in(u,v):mv=add(mv,pair_vector(name,g,u,v,e,5 if five else 0))
        marked[e]=mv
        markkeys[e]=Counter(key for(u,v),key in keys.items()if e not in(u,v)and key in RULE.orbits['connected'])
    return total,marked,connected,markkeys,{u:R.tree_key(g,[u])for u in g}


@lru_cache(None)
def components(ctx,marked_hub):
    a,b,c,k=ctx
    values={};marked={}
    for name,g in MODELS.items():
        total,mvs,connected,markkeys,rootkeys=static_component(name,k==5)
        cb={key:bonusvector('connected',key,ctx)for key in connected}
        rootvals={u:ZERO for u in g}
        for kind,hub,mult,excluded in [('a1cross','A1',a,('A1',)),
                ('hubcross','I2(5)',2*b,('A1','I2(5)')),
                ('h3cross','H3',3*c,H.HUBS)]:
            if mult and name not in excluded:
                factor=F(mult-1,mult) if marked_hub==hub else F(1)
                for u in g:rootvals[u]=add(rootvals[u],bonusvector(kind,rootkeys[u],ctx),factor)
        for value in rootvals.values():total=add(total,value)
        for key,count in connected.items():total=add(total,cb[key],count)
        if name=='I2(5)'and b:total=add(total,samevector(ctx),F(1,b))
        values[name]=total
        for e,mv in mvs.items():
            if marked_hub is None and name in H.HUBS:continue
            if marked_hub is not None and name!=marked_hub:continue
            for u in g:
                if u!=e:mv=add(mv,rootvals[u])
            for key,count in markkeys[e].items():mv=add(mv,cb[key],count)
            # A marked I25 root removes2(b-1) of2b(b-1) same-hub pairs.
            # The b unmarked copies contribute gamma/b each; the marked
            # copy contributes zero after its exact gamma/b debit.
            marked[name,e]=tuple(mv)
    return values,marked


@lru_cache(None)
def floating_components(ctx,mh):
    vals,marks=components(ctx,mh)
    return np.array([[float(z)for z in p]for p in vals.values()]),np.array([[float(z)for z in p]for p in marks.values()])


def oracle(x,exact=False):
    best=None;winner=None;fronts=[]
    evaluate=lambda v: v[0]+sum(a*b for a,b in zip(v[1:],x))
    for ctx in contexts(21):
        a,b,c,k=ctx;remaining=21-a-2*b-3*c-4*k
        for mh in [None]+[hub for hub,count in zip(H.HUBS,(a,b,c)) if count]:
            vals,marks=components(ctx,mh)
            if exact:
                costs={n:evaluate(v) for n,v in vals.items()}
                markedcosts={n:evaluate(v)for n,v in marks.items()}
            else:
                va,ma=floating_components(ctx,mh);xx=np.array([1.0]+list(x))
                costs=dict(zip(vals,va@xx));markedcosts=dict(zip(marks,ma@xx))
            dp={0:(F(0) if exact else 0.0,())}
            for rank in range(1,remaining+1):
                choices=[]
                for name,g in MODELS.items():
                    if name in RESERVED or rank-len(g) not in dp:continue
                    score,forest=dp[rank-len(g)]
                    choices.append((score+costs[name],forest+(name,)))
                if choices:dp[rank]=max(choices,key=lambda z:z[0])
            constant=sum(count*costs[name]for name,count in zip(RESERVED,ctx))
            local=None;record=None
            for (name,e),mv in marks.items():
                if name in RESERVED:
                    if not ctx[RESERVED.index(name)]:continue
                    rank=remaining;score=constant-costs[name]+markedcosts[name,e]
                else:
                    rank=remaining-len(MODELS[name]);score=constant+markedcosts[name,e]
                if rank not in dp:continue
                score+=dp[rank][0]
                if local is None or score>local:
                    local=score;record=dict(context=ctx,marked_hub=mh,marked_component=name,
                        facet_root=e,residual_components=dp[rank][1])
            if local is None:continue
            p=marks[record['marked_component'],record['facet_root']]
            counts=Counter(dict(zip(RESERVED,ctx)))
            if record['marked_component'] in RESERVED:counts[record['marked_component']]-=1
            counts.update(record['residual_components'])
            for name,count in counts.items():p=add(p,vals[name],count)
            record['profile']=p
            if exact:fronts.append(dict(**record,maximum=str(local)))
            if best is None or local>best:best,winner=local,record
    return best,winner,fronts


def face_corner(vertices,edges,pair,mark,external):
    pcs=R.components(R.graph(vertices,edges));names=[R.type_of(g) for g in pcs]
    ctx=tuple(names.count(n)+e for n,e in zip(RESERVED,external))
    a,b=pair;lc=next(g for g in pcs if a in g);rc=next(g for g in pcs if b in g)
    if b in lc:
        p=list(pair_vector(R.type_of(lc),lc,a,b,mark,ctx[3]))
        p[0]+=oldbonus('connected',R.tree_key(lc,[a,b]),ctx)
        return tuple(p)
    p=list(ZERO)
    for kind,hub,mult,excluded in [('a1cross','A1',ctx[0],('A1',)),
         ('hubcross','I2(5)',2*ctx[1],('A1','I2(5)')),
         ('h3cross','H3',3*ctx[2],H.HUBS)]:
        lt,rt=R.type_of(lc),R.type_of(rc)
        if lt==hub and rt not in excluded:p[0]+=oldbonus(kind,R.tree_key(rc,[b]),ctx)/mult
        elif rt==hub and lt not in excluded:p[0]+=oldbonus(kind,R.tree_key(lc,[a]),ctx)/mult
    return tuple(p)


def triangle_constraints():
    tri=H.read('fresh_sqrt5_triangle_weights_20261015');literal=[]
    for index,row in enumerate(tri['markings']):
        literal.append((f'core{index}',row['order'],row['white'],row['black'],
            [(a,b,m)for(a,b),m in zip(combinations(range(row['order']),2),row['labels'])]))
    for index in (46,47,51,55):
        row=tri['markings'][index]
        for tail in range(1,18):
            edges=[(a,b,m)for(a,b),m in zip(combinations(range(5),2),row['labels'])]
            edges += [(0,5,3)]+[(i,i+1,3)for i in range(5,4+tail)]
            literal.append((f'core{index}_tail{tail}',5+tail,row['white']+list(range(5,5+tail)),row['black'],edges))
    for tail in range(1,20):
        edges=[(0,2,3),(1,2,10),(0,3,3)]+[(i,i+1,3)for i in range(3,2+tail)]
        literal.append((f'I210_tail{tail}',3+tail,list(range(3,3+tail)),[0,1,2],edges))
    constraints={};cases=0
    for label,n,white,black,edges in literal:
        remaining=22-n
        compiled=[]
        marks=white+([-1]if remaining else[])
        for omitted in black:
            pair=[b for b in black if b!=omitted];a,b=pair
            pcs=R.components(R.graph(white+pair,edges));names=[R.type_of(g)for g in pcs]
            local=tuple(names.count(name)for name in RESERVED)
            lc=next(g for g in pcs if a in g);rc=next(g for g in pcs if b in g)
            kind,key=None,None
            if b in lc:
                kind='connected';key=R.tree_key(lc,pair)
                vectors={mark:pair_vector(R.type_of(lc),lc,a,b,mark,0)for mark in marks}
                vectors5={mark:pair_vector(R.type_of(lc),lc,a,b,mark,5)for mark in marks}
            else:
                vectors={mark:ZERO for mark in marks};vectors5=vectors
                lt,rt=R.type_of(lc),R.type_of(rc)
                if lt==rt=='I2(5)':kind='same_i25'
                for kd,hub,excluded in [('a1cross','A1',('A1',)),('hubcross','I2(5)',('A1','I2(5)')),('h3cross','H3',H.HUBS)]:
                    if lt==hub and rt not in excluded:kind,key=kd,R.tree_key(rc,[b])
                    elif rt==hub and lt not in excluded:kind,key=kd,R.tree_key(lc,[a])
            compiled.append((local,kind,key,vectors,vectors5))
        for ext in contexts(remaining):
            constant=ZERO;tables=[]
            for local,kind,key,vectors,vectors5 in compiled:
                ctx=tuple(a+b for a,b in zip(local,ext))
                if kind=='same_i25':constant=add(constant,samevector(ctx),F(1,2*ctx[1]*(ctx[1]-1)))
                elif kind:
                    mult={'connected':1,'a1cross':ctx[0],'hubcross':2*ctx[1],'h3cross':3*ctx[2]}[kind]
                    assert mult>0
                    constant=add(constant,bonusvector(kind,key,ctx),F(1,mult))
                tables.append(vectors5 if ctx[3]==5 else vectors)
            for mark in marks:
                p=ZERO
                for table in tables:p=add(p,table[mark])
                p=add(p,constant)
                cases+=1
                if p not in constraints:constraints[p]=dict(family=label,external_counts=ext,facet_root=mark)
    return constraints,dict(literal_families=len(literal),marked_context_cases=cases,distinct_profiles=len(constraints))


def quad_constraints():
    old=H.read('fresh_sqrt5_parameter_weight_stress_20261015');q=H.read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    rows=[r['row']for r in old['faces']if r['row']['nodes']<=23]
    rows += [r['row']for r in q['records']if r['result']['singular'] and r['row']['nodes']<=23]
    assert len(rows)==1979
    count=0;profiles={}
    for rowindex,row in enumerate(rows):
        corners=H.corners(row)
        assert all(c['baseline']==0 for c in corners)
        matrix=R.label_matrix(row);edges=[(u,v,matrix[u][v])for u,v in combinations(range(row['nodes']),2)]
        same=[]
        for pair in product((0,2),(1,3)):
            pcs=R.components(R.graph(row['white']+list(pair),edges))
            lc=next(g for g in pcs if pair[0]in g);rc=next(g for g in pcs if pair[1]in g)
            same.append(pair[1]not in lc and R.type_of(lc)==R.type_of(rc)=='I2(5)')
        for a,b,c in H.contexts(23-row['nodes']):
            total=ZERO
            for corner,is_same in zip(corners,same):
                aa,bb,cc=corner['a']+a,corner['b']+b,corner['c']+c
                ctx=(aa,bb,cc,0)
                total=add(total,bonusvector('connected',corner['connected'],ctx))
                for kind,mult in [('a1cross',aa),('hubcross',2*bb),('h3cross',3*cc)]:
                    if corner[kind] in RULE.orbits[kind]:
                        assert mult>0
                        total=add(total,bonusvector(kind,corner[kind],ctx),F(1,mult))
                if is_same:total=add(total,samevector(ctx),F(1,2*bb*(bb-1)))
            count+=1;profiles.setdefault(total,dict(row_index=rowindex,external_counts=(a,b,c)))
    return profiles,dict(rows=len(rows),common_three_count_contexts=count,distinct_profiles=len(profiles),
        rationale='Every inherited zero-baseline quadrilateral is imposed as an exact affine constraint. The old coefficients may decrease in the three targeted families/context, and alternative corners may receive new payment. Discounts do not touch these rows. Other newly affected paid-corner quadrilaterals remain unproved.')


def saved(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {k:saved(w)for k,w in v.items()}
    if isinstance(v,(tuple,list)):return [saved(w)for w in v]
    return v


def literal_cross_check(witness):
    """Count individual allowed root pairs, independently of hub debit factors."""
    ctx=tuple(witness['context']);counts=Counter(dict(zip(RESERVED,ctx)))
    name=witness['marked_component'];e=witness['facet_root']
    if name in RESERVED:counts[name]-=1
    counts.update(witness['residual_components'])
    forest=[(name,e)]+[(n,None)for n,count in counts.items()for _ in range(count)]
    assert sum(len(MODELS[n])for n,_ in forest)==21
    answer=ZERO;roots=[]
    for name,mark in forest:
        g=MODELS[name];total,mvs,connected,markkeys,rootkeys=static_component(name,ctx[3]==5)
        v=total if mark is None else mvs[mark]
        keys=connected if mark is None else markkeys[mark]
        for key,count in keys.items():v=add(v,bonusvector('connected',key,ctx),count)
        answer=add(answer,v)
        roots.append((name,[(u,rootkeys[u])for u in g if u!=mark]))
    cross=ZERO
    for (lt,ls),(rt,rs) in combinations(roots,2):
        for _,lk in ls:
            for _,rk in rs:
                if lt==rt=='I2(5)':cross=add(cross,samevector(ctx),F(1,2*ctx[1]*(ctx[1]-1)))
                for kind,hub,mult,excluded in [('a1cross','A1',ctx[0],('A1',)),
                        ('hubcross','I2(5)',2*ctx[1],('A1','I2(5)')),
                        ('h3cross','H3',3*ctx[2],H.HUBS)]:
                    if lt==hub and rt not in excluded:cross=add(cross,bonusvector(kind,rk,ctx),F(1,mult))
                    elif rt==hub and lt not in excluded:cross=add(cross,bonusvector(kind,lk,ctx),F(1,mult))
    answer=add(answer,cross)
    assert answer==tuple(witness['profile'])
    return sum(len(rs)for _,rs in roots)


def run():
    triangles,tc=triangle_constraints();print('TRIANGLES',tc,flush=True)
    quads,quad=quad_constraints();print('QUADS',quad,flush=True)
    obligations=[(p,F(2),'triangle',w)for p,w in triangles.items()]+[(p,F(1),'quad',w)for p,w in quads.items()]
    rows=[[-float(v)for v in p[1:]]+[0.0]for p,_,_,_ in obligations]
    rhs=[float(p[0]-target)for p,target,_,_ in obligations]
    witnesses=[];seen=set();result=None
    for iteration in range(200):
        result=linprog([0.0]*N+[1.0],A_ub=rows,b_ub=rhs,
            bounds=[(0,1),(0,1)]+[(0,None)]*(N-2)+[(0,None)],method='highs')
        if not result.success:raise AssertionError(result.message)
        maximum,witness,_=oracle(result.x[:N]);key=tuple(witness['profile'])
        print('ITER',iteration,'LP',result.fun,'ORACLE',maximum,'DISCOUNTS',result.x[:2],
              'WITNESS',witness['context'],witness['marked_component'],witness['residual_components'],flush=True)
        if maximum-result.fun<1e-7:break
        assert key not in seen;seen.add(key);witnesses.append(witness)
        rows.append([float(v)for v in key[1:]]+[-1.0]);rhs.append(-float(key[0]))
    else:raise AssertionError('Separation did not converge')
    # Rational LP vertex recovery from tight independent rows, with exact
    # feasibility checked against every triangle and the complete marked DP.
    from sympy import Matrix, Rational
    active=[];arhs=[]
    for p,target,_,_ in obligations:
        if abs(float(p[0])+sum(float(v)*x for v,x in zip(p[1:],result.x[:N]))-float(target))<1e-7:
            active.append(list(p[1:])+[F(0)]);arhs.append(target-p[0])
    for w in witnesses:
        p=w['profile']
        if abs(float(p[0])+sum(float(v)*x for v,x in zip(p[1:],result.x[:N]))-result.fun)<1e-7:
            active.append(list(p[1:])+[F(-1)]);arhs.append(-p[0])
    for j,x in enumerate(result.x):
        if abs(x)<1e-7 or(j<2 and abs(x-1)<1e-7):
            v=[F(0)]*(N+1);v[j]=1;active.append(v);arhs.append(F(round(x)))
    A=Matrix([[Rational(v.numerator,v.denominator) if isinstance(v,F)else v for v in r]for r in active])
    pivots=list(A.T.rref()[1]);assert len(pivots)==N+1,(len(pivots),N+1)
    sol=A[pivots,:].inv()*Matrix([Rational(arhs[i].numerator,arhs[i].denominator)for i in pivots])
    x=[F(int(v.p),int(v.q))for v in sol]
    assert min(x)>=0 and max(x[:2])<=1
    assert all(p[0]+sum(v*y for v,y in zip(p[1:],x))>=target for p,target,_,_ in obligations)
    maximum,witness,fronts=oracle(x[:N],exact=True)
    assert maximum==x[-1]
    assert all(literal_cross_check(w)==20 for w in fronts)
    # Dual of min z with A*x<=b and x>=0. Nonnegative multipliers obey
    # c+A^T*mu>=0 and certify z>=-b*mu. Include discount upper bounds.
    exactrows=[([-v for v in p[1:]]+[F(0)],p[0]-target,
                 dict(kind=kind,witness=w,profile=p,target=target))for p,target,kind,w in obligations]
    exactrows += [(list(w['profile'][1:])+[F(-1)],-w['profile'][0],dict(kind='vertex',witness=w))for w in witnesses]
    support=[(exactrows[i],-v)for i,v in enumerate(result.ineqlin.marginals)if v<-1e-8]
    for j in range(2):
        if result.upper.marginals[j]<-1e-8:
            row=[F(0)]*(N+1);row[j]=1
            support.append(((row,F(1),dict(kind='discount_upper',variable=j)),-result.upper.marginals[j]))
    positive=[j for j,y in enumerate(x)if y]
    rat=lambda v:Rational(v.numerator,v.denominator)if isinstance(v,F)else Rational(v)
    da=Matrix([[rat(row[0][j])for row,_ in support]for j in positive])
    db=Matrix([-int(j==N)for j in positive]);ds=da.gauss_jordan_solve(db);assert not ds[1]
    dual=[F(int(v.p),int(v.q))for v in ds[0]];assert min(dual)>=0
    assert all(F(j==N)+sum(mu*row[0][j]for(row,_),mu in zip(support,dual))>=0 for j in range(N+1))
    assert -sum(mu*row[1]for(row,_),mu in zip(support,dual))==maximum
    paths=[f'scripts/{TAG}.py',f'scripts/{D.TAG}.py',f'results/{D.TAG}.json',
        f'scripts/{H.TAG}.py',f'results/{H.TAG}.json',f'results/{H.PRIMARY}_cap2.json',
        'results/fresh_sqrt5_triangle_weights_20261015.json',
        'results/fresh_sqrt5_parameter_weight_stress_20261015.json',
        'results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json',
        'results/fresh_sqrt5_marked_discount_triangles_20261016.json',
        'results/fresh_sqrt5_marked_old_bonus_context_scaling_20261016.json']
    return saved(dict(status='PASS_EXACT_RESTRICTED_REALLOCATION_PROPOSAL',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=21,
        facet_dimension=20,target='35',variable_order=['long_local33_distance2_discount','long_adjacent3_discount']+CHANNELS+SCALES+ADDITIONS,
        exact_primal=x,exact_maximum=maximum,target_met=maximum<=35,
        retained_old_rule='n22 H3 cap2, lookup A1_count+1, actual denominators; three targeted source/context coefficients may decrease and eight alternative orbits plus same-I25 pairs may increase',
        H4_exception='All H4 compensation channels are disabled only when the actual forest contains five H4 components.',
        triangle_diagnostics=tc,triangle_profiles=[dict(profile=p,witness=w)for p,w in triangles.items()],
        inherited_quad_obligations=quad,quad_profiles=[dict(profile=p,witness=w)for p,w in quads.items()],separation_iterations=iteration+1,
        generated_vertex_constraints=witnesses,exact_maximizer=witness,all_context_marked_maxima=fronts,
        literal_root_pair_recounts=len(fronts),exact_model_lower_bound=maximum,
        exact_dual=[dict(multiplier=mu,**row[2])for(row,_),mu in zip(support,dual)],
        scope='Counterexample-guided proposal with complete finite marked rank21 vertex DP and inherited baseline-zero quadrilaterals. Only the 76 minimal triangle cores and selected connected-white extensions are covered. All other affected triangles and all newly affected paid-corner quadrilaterals require fresh completion; no dimension bound follows.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(p.read_text())
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['exact_maximum'],'target_met',result['target_met'],flush=True)
