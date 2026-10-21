#!/usr/bin/env python3
"""Global dimension-21 payment proposal: complete vertices, target40."""
import argparse,json,os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import highspy,flint
import fresh_sqrt5_h4_e8_context_split_20261016 as X
M=X.M;S=X.S;ROOT=M.ROOT;N=X.N;TARGET=40
TAG='fresh_sqrt5_global_discount_20261016'
INITIAL=f'data/enumerations/{TAG}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json');RESUME=True
matrix=S.matrix;append_rows=S.append_rows

def literal_profile(front):
    """Enumerate all unordered root pairs in the full rank21 forest."""
    ctx=tuple(front['context']);counts=Counter(dict(zip(M.RESERVED,ctx[:7])))
    assert all(n not in M.RESERVED for n in front['residual_components'])
    counts.update(front['residual_components']);forest=[n for n,c in counts.items()for _ in range(c)]
    assert sum(len(M.MODELS[n])for n in forest)==21
    assert ctx[8]==int(ctx[3]in(1,2,3)and'E8'in forest)
    answer={};roots=[]
    for name in forest:
        total,mvs,con,mc,rk=M.static(name,ctx[3],int(ctx[0]>0),False,ctx[7])
        answer=M.plus(answer,total)
        for key,count in con.items():answer=M.plus(answer,M.bonus('connected',key,ctx),count)
        roots.append((name,[(rk[u],len(M.MODELS[name][u])==1)for u in M.MODELS[name]]))
    for(lt,ls),(rt,rs)in combinations(roots,2):
        for lk,ll in ls:
            for rk,rl in rs:
                if lt=='H3'and not ll:answer=M.plus(answer,M.center_bonus(rt,rk,ctx),F(1,ctx[2]))
                if rt=='H3'and not rl:answer=M.plus(answer,M.center_bonus(lt,lk,ctx),F(1,ctx[2]))
                mo=M.mixed_orbit(lt,lk,rt,rk)
                if mo:
                    name,key,m=M.MIX_ORBITS[mo];h=ctx[4];j=ctx[M.RESERVED.index(name)]
                    answer=M.plus(answer,M.mixed(mo,ctx),F(1,2*m*h*j))
                if lt==rt and lt in M.LEAF_COUNTS and ll and rl:
                    h=ctx[M.RESERVED.index(lt)];m=M.LEAF_COUNTS[lt]
                    answer=M.plus(answer,M.twin(lt,ctx),F(2,m*m*h*(h-1)))
                if lt==rt=='I2(5)':answer=M.plus(answer,M.same(ctx),F(1,2*ctx[1]*(ctx[1]-1)))
                for kind,hub,mult,excluded in [('a1cross','A1',ctx[0],('A1',)),('hubcross','I2(5)',2*ctx[1],('A1','I2(5)')),('h3cross','H3',3*ctx[2],M.HUBS)]:
                    if lt==hub and rt not in excluded:answer=M.plus(answer,M.bonus(kind,rk,ctx),F(1,mult))
                    elif rt==hub and lt not in excluded:answer=M.plus(answer,M.bonus(kind,lk,ctx),F(1,mult))
    return X.remap(answer,ctx[8])

def assembled_profile(front):
    vals,_=X.components(tuple(front['context']),None);p={}
    counts=Counter(dict(zip(M.RESERVED,front['context'][:7])));counts.update(front['residual_components'])
    for name,count in counts.items():p=M.plus(p,vals[name],count)
    return p

def oracle(x,exact=False):
    fronts=[]
    for base in M.contexts(21):
        oldctx=base+(int(base[3]==4 and(base[5]or base[6])),)
        remaining=21-sum(r*c for r,c in zip(M.RANKS,base));tracked=base[3]in(1,2,3)
        for flag in(range(2)if tracked and remaining>=8 else(0,)):
            ctx=oldctx+(flag,);vals,_=X.components(ctx,None)
            if exact:costs={n:M.value(p,x)for n,p in vals.items()}
            else:
                nv,_=M.numerical(ctx,None)
                costs={n:p[0]+sum(z*x[j]for j,z in p[1])for n,p in nv.items()}
            dp={(0,0):(F(0)if exact else 0.0,())}
            for rank in range(1,remaining+1):
                for name,g in M.MODELS.items():
                    if name in M.RESERVED:continue
                    for oldflag in(0,1):
                        entry=dp.get((rank-len(g),oldflag))
                        if entry is None:continue
                        newflag=int(bool(oldflag or(tracked and name=='E8')))
                        score=entry[0]+costs[name];state=(rank,newflag)
                        if state not in dp or score>dp[state][0]:dp[state]=(score,entry[1]+(name,))
            if(remaining,flag)not in dp:continue
            const=sum(c*costs[name]for name,c in zip(M.RESERVED,base))
            front=dict(context=ctx,residual_components=dp[remaining,flag][1])
            fronts.append(dict(**front,maximum=const+dp[remaining,flag][0],profile=assembled_profile(front)))
    return fronts

def install():
    X.install();M.N=N;M.TARGET=TARGET;M.oracle=oracle;M.literal_profile=literal_profile

def prepare_initial():
    path=ROOT/INITIAL
    if path.exists():return
    old=json.loads((ROOT/X.INITIAL).read_text());seen=set();vertices=[]
    for rec in old['vertex_cuts']:
        meta=rec['metadata'];remaining=list(meta['residual_components'])
        if meta['marked_component']not in M.RESERVED:remaining.append(meta['marked_component'])
        context=tuple(meta['context']);remaining.sort();key=(context,tuple(remaining))
        if key in seen:continue
        seen.add(key);front=dict(context=list(context),residual_components=remaining)
        p=literal_profile(front);assert p==assembled_profile(front)
        row={j:z for j,z in p.items()if j>=0};row[N]=F(-1)
        vertices.append(dict(row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(-p.get(-1,F(0))),metadata=front))
        if len(vertices)%1000==0:print('GLOBAL_LITERAL_FORESTS',len(vertices),flush=True)
    paths=[f'scripts/{TAG}.py',X.INITIAL]
    path.write_text(json.dumps(dict(status='GLOBAL_LITERAL_FOREST_HANDOFF_NOT_A_CERTIFICATE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},vertex_cuts=vertices,
        steps=[],latest_numeric_weights=None),separators=(',',':'))+'\n')

def faces():
    obligations={};old=M.H.read('fresh_sqrt5_parameter_weight_stress_20261015')
    queue=M.H.read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    rows=[('prior',i,r['row'])for i,r in enumerate(old['faces'])if r['row']['nodes']<=23]
    rows += [('finite_queue',r['index'],r['row'])for r in queue['records']if r['result']['singular']and r['row']['nodes']<=23]
    oldcases=0
    for family,index,row in rows:
        es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        compiled=[X.OLD_COMPILE(row['white']+[a,b],es,[a,b],[None])for a,b in product((0,2),(1,3))]
        assert all(not any(v[3][None])for v in compiled)
        for ext in M.contexts(23-row['nodes']):
            p,tables=X.OLD_ADD(compiled,ext)
            assert not any(j in X.SPLIT for j in p)and all(not t[None]for t in tables)
            obligations.setdefault((M.key(p),F(1)),dict(kind='quad',family=family,index=index,external_counts=list(ext)))
            oldcases+=1
    tri=json.loads((ROOT/f'results/{X.A.TRI}_rank22.json').read_text())
    quad=json.loads((ROOT/f'results/{X.A.QUAD}.json').read_text())
    assert len(tri['records'])==996 and tri['complete_through_bound']and not tri['pending']
    todo=[('triangle',r['index'],r['row'],r['status'])for r in tri['records']]
    todo += [('quad',r['class_index'],r['row'],'SEALED_SINGULAR')for r in quad['immediate_sealed_weight_obligations']]
    todo += [('quad',r['class_index'],r['row'],'SINGULAR_EXTENSIBLE_STOP')for r in quad['records']if r['class_index']in(4446,4448)]
    counts=dict(inherited_baseline_zero_quad_rows=len(rows),inherited_quad_common_contexts=oldcases,
                triangle_rows=996,paid_quad_stopping_rows=37,triangle_common_contexts=0,paid_quad_common_contexts=0)
    for pos,(kind,index,row,status)in enumerate(todo):
        n=row['nodes'];white=row['white'];remaining=(22 if kind=='triangle'else 23)-n
        es=[(a,b,m)for(a,b),m in zip(combinations(range(n),2),row['labels'])]
        pairs=[[b for b in row['black']if b!=o]for o in row['black']]if kind=='triangle'else list(product((0,2),(1,3)))
        compiled=[X.compile_corner(white+list(pair),es,pair,[None])for pair in pairs]
        for ext in X.external_contexts(remaining):
            p,tables=X.add_bonus(compiled,ext)
            for table in tables:p=M.plus(p,table[None])
            obligations.setdefault((M.key(p),F(2 if kind=='triangle'else 1)),dict(kind=kind,
                family='complete_triangles'if kind=='triangle'else'paid_quad_stopping',index=index,
                producer_status=status,external_counts=list(ext)))
            counts['triangle_common_contexts'if kind=='triangle'else'paid_quad_common_contexts']+=1
        if pos%100==0:print('GLOBAL_FACE_REBUILD',pos,'profiles',len(obligations),flush=True)
    obligations[M.key({0:F(-1)}),F(0)]=dict(kind='model_bound',name='distance2_discount_zero')
    obligations[M.key({-1:F(2,3),1:F(-1)}),F(0)]=dict(kind='model_bound',name='adjacent_discount_at_most_two_thirds')
    return [(dict(p),target,meta)for(p,target),meta in obligations.items()],counts

def checkpoint(vertices,steps,weights):
    paths=[f'scripts/{TAG}.py',f'scripts/{X.TAG}.py',f'results/{X.TAG}.json']
    doc=dict(status='PARTIAL_GLOBAL_CHECKPOINT_NOT_A_CERTIFICATE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},vertex_cuts=vertices,steps=steps,
        latest_numeric_weights=list(map(float,weights))if weights is not None else None)
    temporary=CHECKPOINT.with_suffix('.tmp');temporary.write_text(json.dumps(doc,separators=(',',':'))+'\n');os.replace(temporary,CHECKPOINT)

def verify(obligations,cert):
    if cert['kind']=='exact_feasible_restricted_rule':
        x=list(map(F,cert['primal']));assert len(x)==N and min(x)>=0 and max(x[:2])<=1
        assert all(M.value(p,x)>=t for p,t,_ in obligations)
        fronts=oracle(x,True);maximum=max(f['maximum']for f in fronts)
        assert maximum==F(cert['maximum'])<=TARGET
        assert all(literal_profile(f)==f['profile']for f in fronts)
        return dict(maximum=str(maximum),global_contexts=len(fronts),fronts=[dict(
            **{k:v for k,v in f.items()if k not in('profile','maximum')},maximum=str(f['maximum']))for f in fronts])
    lookup={(M.key({j:-z for j,z in p.items()if j>=0}),p.get(-1,F(0))-t):w for p,t,w in obligations}
    columns=Counter();lower=F(0)
    for term in cert['terms']:
        row={j:F(z)for j,z in term['row']};rhs=F(term['rhs']);mu=F(term['multiplier']);m=term['metadata'];assert mu>=0
        if m['kind']=='vertex':
            p=literal_profile(m);assert p==assembled_profile(m)
            wanted={j:z for j,z in p.items()if j>=0};wanted[N]=F(-1)
            assert row==wanted and rhs==-p.get(-1,F(0))
        else:assert lookup[M.key(row),rhs]==m
        for j,z in row.items():columns[j]+=mu*z
        lower-=mu*rhs
    upper=list(map(F,cert['discount_upper_multipliers']));assert min(upper)>=0
    for j,z in enumerate(upper):columns[j]+=z;lower-=z
    columns[N]+=1;assert min(columns.values(),default=F(0))>=0 and lower==F(cert['lower_bound'])>TARGET
    return dict(lower_bound=str(lower),terms=len(cert['terms']),exact_weak_duality=True)

# TARGET-parametrized copies of the frozen sparse solver follow. They keep
# fixed-budget minimum-total-bonus optimization and exact rational recovery.
def highs_model(rows,rhs,n):
    A=matrix(rows,n);lp=highspy.HighsLp();lp.num_col_=n;lp.num_row_=len(rows)
    lp.col_cost_=np.array([0,0]+[1]*(n-3)+[0],dtype=float)
    lp.col_lower_=np.array([0]*(n-1)+[TARGET],dtype=float)
    lp.col_upper_=np.array([1,1]+[highspy.kHighsInf]*(n-3)+[TARGET],dtype=float)
    lp.row_lower_=np.full(len(rows),-highspy.kHighsInf);lp.row_upper_=np.array(list(map(float,rhs)))
    lp.a_matrix_.format_=highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_=A.indptr;lp.a_matrix_.index_=A.indices;lp.a_matrix_.value_=A.data
    h=highspy.Highs();h.setOptionValue('output_flag',False);h.setOptionValue('solver','simplex')
    h.setOptionValue('threads',1);h.setOptionValue('random_seed',0)
    h.setOptionValue('primal_feasibility_tolerance',1e-9);h.setOptionValue('dual_feasibility_tolerance',1e-9)
    assert h.passModel(lp)==highspy.HighsStatus.kOk
    return h

def exact_basis(h,rows,rhs,n):
    import heapq
    basis=h.getBasis();B=highspy.HighsBasisStatus
    basic=[j for j,s in enumerate(basis.col_status)if s==B.kBasic]
    tight=[i for i,s in enumerate(basis.row_status)if s!=B.kBasic]
    assert len(basic)==len(tight),(len(basic),len(tight))
    x=[F(0)]*n
    for j,s in enumerate(basis.col_status):
        if s==B.kBasic:continue
        if j==n-1:x[j]=F(TARGET)
        elif s==B.kUpper:assert j<2;x[j]=F(1)
        else:assert s in(B.kLower,B.kZero,B.kNonbasic)
    active=set(basic);eqs={};right={};colrows={j:set()for j in basic};heap=[];version={}
    q=lambda z:flint.fmpq(z.numerator,z.denominator)
    for k,i in enumerate(tight):
        assert basis.row_status[i]in(B.kUpper,B.kNonbasic)
        eqs[k]={j:q(z)for j,z in rows[i].items()if j in active}
        right[k]=q(rhs[i]-sum(z*x[j]for j,z in rows[i].items()if j not in active))
        assert eqs[k]
        for j in eqs[k]:colrows[j].add(k)
        version[k]=0;heapq.heappush(heap,(len(eqs[k]),0,k))
    pivots=[];initial=len(active)
    print('EXACT_BASIS_SPARSE_START',initial,'nonzeros',sum(map(len,eqs.values())),flush=True)
    while eqs:
        while heap:
            size,v,k=heapq.heappop(heap)
            if k in eqs and version[k]==v:break
        else:raise AssertionError('Lost sparse pivot rows')
        assert size==len(eqs[k])and size>0
        # Sparse rows make exact substitution cheap. Leave a small dense
        # remainder to FLINT rather than create large fill in Python.
        if len(active)<=400 or size>16:break
        j=min(eqs[k],key=lambda j:(len(colrows[j]),j));a=eqs[k][j]
        normalized={v:z/a for v,z in eqs[k].items()if v!=j};b=right[k]/a
        neighbors=sorted(colrows[j]-{k})
        for v in eqs[k]:colrows[v].remove(k)
        del eqs[k];del right[k];active.remove(j);pivots.append((j,b,normalized))
        for i in neighbors:
            coefficient=eqs[i].pop(j);colrows[j].remove(i);right[i]-=coefficient*b
            for v,z in normalized.items():
                old=eqs[i].get(v,flint.fmpq(0));new=old-coefficient*z
                if new:
                    eqs[i][v]=new
                    if not old:colrows[v].add(i)
                elif old:del eqs[i][v];colrows[v].remove(i)
            assert eqs[i],('Dependent basis row',i,right[i])
            version[i]+=1;heapq.heappush(heap,(len(eqs[i]),version[i],i))
        assert not colrows[j]
        if len(pivots)%1000==0:print('EXACT_SPARSE_ELIMINATED',len(pivots),'remaining',len(active),flush=True)
    active=sorted(active);position={j:k for k,j in enumerate(active)};ids=sorted(eqs)
    assert len(active)==len(ids)
    print('EXACT_DENSE_REMAINDER',len(active),'nonzeros',sum(map(len,eqs.values())),flush=True)
    rational={}
    if active:
        A=flint.fmpq_mat(len(active),len(active));b=flint.fmpq_mat(len(active),1)
        for k,i in enumerate(ids):
            b[k,0]=right[i]
            for j,z in eqs[i].items():A[k,position[j]]=z
        solution=A.solve(b,algorithm='dixon');assert A*solution==b
        rational.update({j:solution[k,0]for j,k in position.items()})
    for j,b,row in reversed(pivots):rational[j]=b-sum((z*rational[v]for v,z in row.items()),flint.fmpq(0))
    for j,v in rational.items():x[j]=F(str(v))
    assert x[-1]==TARGET and min(x)>=0 and max(x[:2])<=1
    assert all(sum(z*x[j]for j,z in row.items())<=bound for row,bound in zip(rows,rhs))
    return x,dict(basic_variables=initial,exact_fraction_equations=len(tight),sparse_pivots=len(pivots),dense_remainder=len(active))

def solve(obligations):
    n=M.N+1;rows=[{j:-z for j,z in p.items()if j>=0}for p,t,w in obligations]
    rhs=[p.get(-1,F(0))-t for p,t,w in obligations];meta=[w for p,t,w in obligations]
    cuts=set();vertices=[];steps=[]
    def cut(front):
        profile=front['profile'];k=M.key(profile)
        if k in cuts:return False
        cuts.add(k);row={j:z for j,z in profile.items()if j>=0};row[M.N]=F(-1)
        bound=-profile.get(-1,F(0));m=dict(kind='vertex',**{k:v for k,v in front.items()if k not in('profile','maximum','kind')})
        rows.append(row);rhs.append(bound);meta.append(m)
        vertices.append(dict(row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(bound),metadata=m))
        return True
    if RESUME and CHECKPOINT.is_file():
        prior=json.loads(CHECKPOINT.read_text())
        for p,wanted in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==wanted
        for rec in prior['vertex_cuts']:
            profile={int(j):F(z)for j,z in rec['row']if int(j)!=M.N};
            if F(rec['rhs']):profile[-1]=-F(rec['rhs'])
            cut(dict(**rec['metadata'],profile=profile))
        steps=prior['steps'];print('RESUMED',len(cuts),'literal cuts',flush=True)
    else:
        for front in M.oracle([0.0]*M.N):cut(front)
    checkpoint(vertices,steps,None)
    h=highs_model(rows,rhs,n)
    for iteration in range(len(steps),200):
        assert h.run()==highspy.HighsStatus.kOk;status=h.getModelStatus()
        if status==highspy.HighsModelStatus.kInfeasible:
            A=matrix(rows,n);objective=np.zeros(n);objective[M.N]=1
            lp=M.linprog(objective,A_ub=A,b_ub=np.array(list(map(float,rhs))),
                bounds=[(0,1),(0,1)]+[(0,None)]*(M.N-1),method='highs')
            assert lp.success and lp.fun>TARGET+1e-7,(lp.message,lp.fun)
            cert=M.dual_certificate(rows,rhs,meta,lp);checkpoint(vertices,steps,lp.x)
            return cert,steps
        assert status==highspy.HighsModelStatus.kOptimal,h.modelStatusToString(status)
        weights=np.array(h.getSolution().col_value);fronts=M.oracle(weights[:M.N]);maximum=max(f['maximum']for f in fronts)
        steps.append(dict(iteration=iteration,fixed_target=TARGET,secondary_objective=h.getObjectiveValue(),
            cuts=len(cuts),maximum=float(maximum),discounts=list(map(float,weights[:2]))))
        print('SECONDARY',iteration,'sum',h.getObjectiveValue(),'cuts',len(cuts),'maximum',maximum,'discounts',weights[:2],flush=True)
        checkpoint(vertices,steps,weights)
        if maximum<=TARGET+1e-7:
            x,basis_record=exact_basis(h,rows,rhs,n);fronts=M.oracle(x[:M.N],True)
            upper=max(f['maximum']for f in fronts);print('EXACT_BASIS',basis_record,'maximum',float(upper),flush=True)
            if upper<=TARGET:
                assert all(M.value(p,x)>=t for p,t,w in obligations)
                return dict(kind='exact_feasible_restricted_rule',primal=list(map(str,x[:M.N])),maximum=str(upper),
                    exact_basis_recovery=basis_record),steps
        before=len(rows)
        for front in fronts:
            if front['maximum']>TARGET+(0 if isinstance(front['maximum'],F)else 1e-7):cut(front)
        assert len(rows)>before,'No new cut despite an unresolved exact feasibility check'
        append_rows(h,rows[before:],rhs[before:],n);checkpoint(vertices,steps,weights)
    raise AssertionError('Secondary iteration cap; literal cuts saved')

def run(saved=None):
    install();obligations,counts=faces()
    print('GLOBAL_MODEL',N,'variables',len(obligations),'profiles',counts,flush=True)
    cert,steps=solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=verify(obligations,cert)
    paths=[f'scripts/{TAG}.py',INITIAL,f'scripts/{X.TAG}.py',f'results/{X.TAG}.json',
           'scripts/fresh_sqrt5_marked_h3_center_sparse_recovery_20261016.py',
           f'results/{X.A.TRI}_rank22.json',f'results/{X.A.QUAD}.json']
    prior=json.loads((ROOT/f'results/{X.TAG}.json').read_text())
    return dict(status='PASS_EXACT_GLOBAL_DISCOUNT_PROPOSAL',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=21,
        target=str(TARGET),weight_scope='All unordered root pairs at each full rank21 vertex; no distinguished facet root.',
        variables=N,variable_order=prior['variable_order'],face_counts=counts,unique_face_profiles=len(obligations),
        certificate=cert,verification=proof,steps=steps,vertex_context_order=prior['vertex_context_order'],
        triangle_external_context_order=prior['triangle_external_context_order'],
        quadrilateral_external_context_order=prior['quadrilateral_external_context_order'],
        unused_marked_channels='Inherited coefficients requiring a distinguished facet root are identically unused; their nonnegative secondary-objective minimum is zero.',
        scope='Global rank21 vertex budget40 with the same adjacent discount,E6two-edge exemption,actualH3cross counts,H4/E8 context and nonnegative compensation family. Distance2discount0,adjacent<=2/3. All996triangle states,37paidquadstops and1979inheritedzeroquads are evaluated without facet marks. Complete finite-forest oracle retains actual E8 presence. Remaining newly affected paidquadextension coverage is open; no dimension claim.',
        optimizer_scope='T fixed40, minimize the sum of all nonnegative bonus coefficients. Discounts are bounded0..1 with additional explicit restrictions. Objective is bounded below0; every bonus coefficient is bounded on an objective sublevel. Exact rational recovery and full global forest DP are required for feasibility.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();install()
    if not args.verify:
        prepare_initial()
        if not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
