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

ROOT=P.ROOT;TAG='fresh_sqrt5_marked_h4_count_free_bonus_20261016'
H,R=P.H,P.R;MODELS=P.MODELS;HUBS=H.HUBS;RESERVED=P.RESERVED
COMP_ORDER=['long_distance2_discount','long_adjacent3_discount'];REMAP={(0,None):0,(1,None):1}
H4_INDICES=set()
for j,channel in enumerate(P.CHANNELS,2):
    if channel['type']=='H4':
        H4_INDICES.add(j)
        for count in range(1,5):
            REMAP[j,count]=len(COMP_ORDER);COMP_ORDER.append(dict(**channel,H4_count=count))
    else:REMAP[j,None]=len(COMP_ORDER);COMP_ORDER.append(channel)
SMALL=len(COMP_ORDER);OLDCOUNT=len(P.RULE.x)-1
SAME={(a,c):SMALL+OLDCOUNT+i for i,(a,c)in enumerate((a,c)for a in range(18)for c in range(2)if a+4+3*c<=21)}
N=SMALL+OLDCOUNT+len(SAME);TARGET=F(35)


def plus(v,w,m=1):
    out=v.copy()
    for j,z in w.items():
        out[j]=out.get(j,F(0))+m*z
        if not out[j]:del out[j]
    return out
def value(v,x):return v.get(-1,0)+sum(z*x[j]for j,z in v.items()if j>=0)
def key(v):return tuple(sorted(v.items()))
def small(v,count):
    out={}
    for i,z in enumerate(v):
        if not z:continue
        j=i-1
        if j==-1:out[j]=z
        elif j in H4_INDICES:
            if count in range(1,5):out[REMAP[j,count]]=z
        else:out[REMAP[j,None]]=z
    return out


def bonus(kind,k,ctx):
    a,b,c=ctx[:3];j=P.RULE.index.get((kind,k,a+1,min(b,2),min(c,1)))
    return {}if j is None else{SMALL+j:F(1)}
def same(ctx):
    a,b,c=ctx[:3]
    return {SAME[a,min(c,1)]:F(1)}if b>=2 else{}


@lru_cache(None)
def static(name,count):
    total,mvs,connected,markkeys,rootkeys=P.static_component(name,count==5)
    return small(total,count),{e:small(v,count)for e,v in mvs.items()},connected,markkeys,rootkeys


@lru_cache(None)
def components(ctx,mh):
    a,b,c,k=ctx;vals={};marks={}
    for name,g in MODELS.items():
        total,mvs,connected,markkeys,rootkeys=static(name,k)
        cb={ky:bonus('connected',ky,ctx)for ky in connected}
        rootvals={u:{}for u in g}
        for kind,hub,mult,excluded in [('a1cross','A1',a,('A1',)),
                ('hubcross','I2(5)',2*b,('A1','I2(5)')),('h3cross','H3',3*c,HUBS)]:
            if mult and name not in excluded:
                factor=F(mult-1,mult)if mh==hub else F(1)
                for u in g:rootvals[u]=plus(rootvals[u],bonus(kind,rootkeys[u],ctx),factor)
        for rv in rootvals.values():total=plus(total,rv)
        for ky,count in connected.items():total=plus(total,cb[ky],count)
        if name=='I2(5)'and b:total=plus(total,same(ctx),F(1,b))
        vals[name]=total
        for e,mv in mvs.items():
            if mh is None and name in HUBS:continue
            if mh is not None and name!=mh:continue
            for u in g:
                if u!=e:mv=plus(mv,rootvals[u])
            for ky,count in markkeys[e].items():mv=plus(mv,cb[ky],count)
            marks[name,e]=mv
    return vals,marks


@lru_cache(None)
def numerical(ctx,mh):
    vals,marks=components(ctx,mh)
    def pack(v):return float(v.get(-1,0)),[(j,float(z))for j,z in v.items()if j>=0]
    return {n:pack(v)for n,v in vals.items()},{n:pack(v)for n,v in marks.items()}


def oracle(x,exact=False):
    fronts=[]
    for ctx in P.contexts(21):
        a,b,c,k=ctx;remaining=21-a-2*b-3*c-4*k
        for mh in [None]+[hub for hub,count in zip(HUBS,(a,b,c))if count]:
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
                    if name in RESERVED or rank-len(g)not in dp:continue
                    old,forest=dp[rank-len(g)];choices.append((old+costs[name],forest+(name,)))
                if choices:dp[rank]=max(choices,key=lambda z:z[0])
            const=sum(count*costs[name]for name,count in zip(RESERVED,ctx));best=None;winner=None
            for(name,e),mv in marks.items():
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
            counts=Counter(dict(zip(RESERVED,ctx)))
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
        for kd,hub,excluded in [('a1cross','A1',('A1',)),('hubcross','I2(5)',('A1','I2(5)')),('h3cross','H3',HUBS)]:
            if lt==hub and rt not in excluded:kind,ky=kd,R.tree_key(rc,[b])
            elif rt==hub and lt not in excluded:kind,ky=kd,R.tree_key(lc,[a])
    return ctx,kind,ky,vs,vs5


def add_bonus(compiled,ext):
    out={};tables=[]
    for local,kind,ky,vs,vs5 in compiled:
        ctx=tuple(a+b for a,b in zip(local,ext))
        if kind=='same':out=plus(out,same(ctx),F(1,2*ctx[1]*(ctx[1]-1)))
        elif kind:
            mult={'connected':1,'a1cross':ctx[0],'hubcross':2*ctx[1],'h3cross':3*ctx[2]}[kind]
            assert mult>0;out=plus(out,bonus(kind,ky,ctx),F(1,mult))
        tables.append({e:small(v,ctx[3])for e,v in(vs5 if ctx[3]==5 else vs).items()})
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
        remaining=22-n;marks=white+([-1]if remaining else[])
        compiled=[compile_corner(white+[b for b in black if b!=o],es,[b for b in black if b!=o],marks)for o in black]
        for ext in P.contexts(remaining):
            common,tables=add_bonus(compiled,ext)
            for mark in marks:
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
        for ea,eb,ec in H.contexts(23-row['nodes']):
            p,_=add_bonus(compiled,(ea,eb,ec,0))
            obligations.setdefault((key(p),F(1)),dict(kind='quad',family=family,index=index,external_counts=(ea,eb,ec)));cases+=1
    counts['quad_rows']=len(rows);counts['quad_common_contexts']=cases
    return [(dict(p),target,json.loads(json.dumps(w)))for(p,target),w in obligations.items()],counts


def literal_profile(front):
    """Recount actual root pairs; do not use marked-hub debit factors."""
    ctx=tuple(front['context']);counts=Counter(dict(zip(RESERVED,ctx)))
    name=front['marked_component'];mark=front['facet_root']
    if name in RESERVED:counts[name]-=1
    assert all(n not in RESERVED for n in front['residual_components'])
    counts.update(front['residual_components'])
    forest=[(name,mark)]+[(n,None)for n,count in counts.items()for _ in range(count)]
    assert sum(len(MODELS[n])for n,_ in forest)==21
    assert tuple(sum(n==typ for n,_ in forest)for typ in RESERVED)==ctx
    answer={};roots=[]
    for name,e in forest:
        total,mvs,con,mc,rk=static(name,ctx[3])
        answer=plus(answer,total if e is None else mvs[e])
        for ky,count in(con if e is None else mc[e]).items():answer=plus(answer,bonus('connected',ky,ctx),count)
        roots.append((name,[rk[u]for u in MODELS[name]if u!=e]))
    for(lt,ls),(rt,rs)in combinations(roots,2):
        for lk in ls:
            for rk in rs:
                if lt==rt=='I2(5)':answer=plus(answer,same(ctx),F(1,2*ctx[1]*(ctx[1]-1)))
                for kind,hub,mult,excluded in [('a1cross','A1',ctx[0],('A1',)),('hubcross','I2(5)',2*ctx[1],('A1','I2(5)')),('h3cross','H3',3*ctx[2],HUBS)]:
                    if lt==hub and rt not in excluded:answer=plus(answer,bonus(kind,rk,ctx),F(1,mult))
                    elif rt==hub and lt not in excluded:answer=plus(answer,bonus(kind,lk,ctx),F(1,mult))
    return answer


def dual_certificate(rows,rhs,meta,lp):
    for denominator in(10**6,10**9,10**12):
        mu=[F(round(float(-v)*denominator),denominator)if v<-1e-10 else F(0)for v in lp.ineqlin.marginals]
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
            counts=Counter(dict(zip(RESERVED,ctx)))
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


def run(saved=None):
    obligations,counts=faces();print('MODEL',N,'variables',len(obligations),'faceprofiles',counts,flush=True)
    cert,steps=solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=verify(obligations,cert)
    paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',f'results/{P.TAG}.json',f'scripts/{H.TAG}.py',f'results/{H.TAG}.json',
        f'results/{H.PRIMARY}_cap2.json','results/fresh_sqrt5_triangle_weights_20261015.json',
        'results/fresh_sqrt5_parameter_weight_stress_20261015.json','results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json']
    variable_order=COMP_ORDER[:]
    variable_order += [dict(kind='inherited_free_bonus',**v,actual_A1_count=v['A1_count']-1)for v in P.PRIMARY['variable_order']]
    variable_order += [dict(kind='same_i25',actual_A1_count=a,I25_count_at_least=2,H3_presence=c)for a,c in SAME]
    assert len(variable_order)==N
    return dict(status='PASS_EXACT_H4_COUNT_FREE_BONUS_MARKED_PROPOSAL',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        dimension=21,facet_dimension=20,target='35',variables=N,variable_order=variable_order,face_counts=counts,
        unique_face_profiles=len(obligations),certificate=cert,verification=proof,steps=steps,
        H4_compensation_context='Separate coefficients for each actual H4 count1,2,3,4, additionally split by whether the facet root lies in that H4. Compensation stays zero at count5.',
        scope='All inherited bonus coefficients are independent nonnegative variables; no old primal weights are fixed. Complete marked finite-forest DP and all inherited baseline-zero quadrilateral obligations are retained. Only selected triangle families are covered. Newly affected triangles and paid-corner quadrilaterals still require complete coverage. No dimension-bound promotion.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    path=ROOT/f'results/{TAG}.json';old=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(old),default=str))
    if args.verify:assert result==old
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
