#!/usr/bin/env python3
"""Exact ordered residual-budget envelope; fresh continuation of frozen164.

Uses the frozen independent Fraction grounded-matrix helper read-only.
All new convex witnesses, support lists and packing decisions are recomputed.
"""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from math import isqrt
from pathlib import Path
import fresh_nikulin_height_charge_20261015_review as H

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_two_budget_frontier_20261015'
B=F(21,20)
CAP=F(3,2)

def sqrt_lower(q):
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    if a*a==q.numerator and b*b==q.denominator:return F(a,b)
    D=2**64;v=isqrt(q.numerator*D*D//q.denominator)
    assert F(v,D)**2<=q<F(v+1,D)**2
    return F(v,D)

def height(ds,U,V):
    s,t=[sum((sqrt_lower(ds[k])for k in ids),F())for ids in(U,V)]
    return max(min(s,t)**2,max(s,t)**2/B**2)

def solve(A,b):
    try:return H.mv(H.inverse(A),b)
    except(StopIteration,ZeroDivisionError):return None

@lru_cache(None)
def closest_matrix(K,nleft):
    """Exact KKT certificate for distance of two finite convex hulls."""
    n=len(K)
    for L in H.subsets(range(nleft)):
        if not L:continue
        for Q in H.subsets(range(nleft,n)):
            if not Q:continue
            active=L+Q
            A=[[F((j<nleft)==first)for j in active]for first in(True,False)]
            M=[[K[i][j]for j in active]+[A[z][k]for z in range(2)]for k,i in enumerate(active)]
            M += [row+[F(),F()]for row in A]
            sol=solve(M,[F()]*len(active)+[F(1),F(1)])
            if sol is None or min(sol[:-2])<0:continue
            x=[F()]*n
            for i,v in zip(active,sol[:-2]):x[i]=v
            lag=sol[-2:]
            slack=[v+lag[int(i>=nleft)]for i,v in enumerate(H.mv(K,x))]
            if min(slack)<0:continue
            value=H.dot(x,H.mv(K,x))
            assert value>=0 and value==-sum(lag)
            assert sum(x[:nleft])==sum(x[nleft:])==1
            assert all(a*b==0 for a,b in zip(x,slack))
            return value,tuple(x),tuple(lag),tuple(slack)
    raise AssertionError(('no convex certificate',K,nleft))

def closest(U,V,R):
    ids=U+V;sg=[1]*len(U)+[-1]*len(V)
    K=tuple(tuple(sg[i]*sg[j]*(R(a,0)+R(b,0)-R(a,b))/2 for j,b in enumerate(ids))for i,a in enumerate(ids))
    e,x,lag,slack=closest_matrix(K,len(U))
    return e,dict(weights=list(map(str,x)),multipliers=list(map(str,lag)),slacks=list(map(str,slack)))

def classical_model(typ,r):
    def potential(k):
        if typ=='C':return[F(int(i<k))for i in range(r)]
        if k==0:return[F()]*r
        if k==1:return[F(1)]+[F()]*(r-1)
        if typ=='D'and k==r-1:return[F(1,2)]*(r-1)+[F(-1,2)]
        return[F(1,2)*int(i<k)for i in range(r)]
    def R(i,j):
        return sum((a-b)**2 for a,b in zip(potential(i),potential(j)))/(2 if typ=='C'else 1)
    ds=([F(1)]+[F(2)]*(r-1)+[F(1)])if typ=='C'else([F(1),F(1)]+[F(4)]*(r-2)+[F(2)])if typ=='B'else([F(1),F(1)]+[F(4)]*(r-3)+[F(1),F(1)])
    return ds,R,potential

def models():
    out=[]
    def emit(name,r,modified,om,ds,R,**extra):
        out.append(dict(name=name,rank=r,modified=modified,om=om,ds=ds,R=R,**extra))
    # A_{m-1}, m2..17. Larger cycles use the stronger new residual-box tail.
    for r in range(1,17):
        m=r+1
        def R(i,j,m=m):return F(((i-j)%m)*(m-(i-j)%m),m)
        emit('A',r,False,(0,m//2),[F(1)]*m,R)
    # B/C tails start at8, D at18; these are proved by literal fixed coordinates.
    for typ,maxrank in [('B',7),('C',7),('D',17)]:
        for r in range(4 if typ=='D'else 2,maxrank+1):
            ds,R,potential=classical_model(typ,r)
            emit(typ,r,False,(0,r),ds,R)
            if typ=='B'and r>=3:emit(typ,r,True,(r-1,0),ds,R)
            if typ=='C':emit(typ,r,True,(1,r),ds,R)
    for name,arms in [('E6',[2,2,2]),('E7',[1,3,3]),('E8',[1,2,5])]:
        A=H.star(arms);ds,di,R,c=H.model(A,[1]*len(A));r=len(A)-1
        for f in[None]+list(range(r+1)):
            emit(name,r,f is not None,H.choices(name,r,f),ds,R,forbidden=f)
    for name,A,scales in [
        ('F4',[[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-2,0],[0,0,-2,4,-2],[0,0,0,-2,4]],[1,1,1,2,2]),
        ('G2',[[2,-1,0],[-1,2,-3],[0,-3,6]],[1,1,3]),
        ('B2modified',[[2,0,-2],[0,2,-2],[-2,-2,4]],[1,1,2])]:
        ds,di,R,c=H.model(A,scales);r=len(A)-1
        for f in([0]if name=='B2modified'else[None]+list(range(r+1))):
            emit(name,r,f is not None,(0,1)if name=='B2modified'else H.choices(name,r,f),ds,R,forbidden=f)
    return out

def support_cost(S,o,ds,R):
    cost=sum((ds[k]*R(o,k)for k in S),F())
    for k,l in combinations(S,2):
        v=R(o,k)+R(o,l)-R(k,l);assert v>=0
        if not v:continue
        # Irrational cross terms are only needed to rule a support out. All
        # surviving domains below3/2 have exact rational cost.
        root=sqrt_lower(ds[k]*ds[l]);cost+=root*v
    return cost

def finite_envelope():
    records=[];positive=[];support_count=0;case_count=0
    for model in models():
        name,r,modified,om,ds,R=[model[k]for k in('name','rank','modified','om','ds','R')]
        T=F(2*(r-int(modified)),9);domains=[]
        for o in om:
            nodes=[k for k in range(r+1)if k!=o and ds[k]*R(o,k)<=CAP]
            domain=[]
            for S in H.subsets(nodes):
                c=support_cost(S,o,ds,R)
                if c>CAP:continue
                # Independently verify surviving costs are exact, not merely
                # lower bounds obtained from a radical interval.
                exact=sum((ds[k]*R(o,k)for k in S),F())
                for k,l in combinations(S,2):
                    v=R(o,k)+R(o,l)-R(k,l)
                    if v:exact+=H.sqrtq(ds[k]*ds[l])*v
                assert c==exact
                domain.append((S,c))
            domains.append(domain)
        cases=[]
        for(L,cp),(Q,cq)in product(*domains):
            if not L and not Q:continue # reviewed zero-cost bound W>=T
            support_count+=1
            for a,b in product((False,True),repeat=2):
                U=L+((om[0],)if a else());V=Q+((om[1],)if b else())
                if not U or not V:continue
                e,dual=closest(U,V,R);h=height(ds,U,V);loss=T-h*e
                case=dict(left=U,right=V,costs=list(map(str,(cp,cq))),minimum_energy=str(e),height_lower=str(h),loss_upper=str(loss),dual=dual)
                cases.append(case);case_count+=1
                if loss>0:
                    positive.append(dict(model=len(records),case=len(cases)-1,name=name,rank=r,modified=modified,
                        left=U,right=V,p=cp,q=cq,loss=loss))
        records.append(dict(name=name,rank=r,modified=modified,forbidden=model.get('forbidden'),omissions=om,
            kernel_squares=list(map(str,ds)),target=str(T),domains=[[[list(s),str(c)]for s,c in d]for d in domains],cases=cases))
    # Exact Pareto table, respecting the single modified factor. Each retained
    # record bounds a full continuous contact domain, not merely a point row.
    unique={}
    for z in positive:unique.setdefault((z['p'],z['q'],z['loss'],z['modified']),z)
    candidates=list(unique.values());table=[]
    for z in candidates:
        if any(y is not z and y['modified']==z['modified']and y['p']<=z['p']and y['q']<=z['q']and y['loss']>=z['loss']and
               (y['p']<z['p']or y['q']<z['q']or y['loss']>z['loss'])for y in candidates):continue
        table.append(z)
    for z in positive:
        j=next(j for j,y in enumerate(table)if y['modified']==z['modified']and y['p']<=z['p']and y['q']<=z['q']and y['loss']>=z['loss'])
        records[z['model']]['cases'][z['case']]['dominator']=j
    assert all(z['left']==z['right']and len(z['left'])==1 for z in table)
    assert(len(records),support_count,case_count,len(positive),len(table))==(91,1116,3630,1395,28)
    return records,table,dict(models=len(records),retained_support_pairs=support_count,continuous_domains=case_count,
        positive_loss_domains=len(positive),pareto_records=len(table))

def packing(table):
    rows=[]
    for p,q,target in[(F(10,9),F(4,3),F(20,9)),(F(1),F(3,2),F(22,9)),(F(3,2),F(1),F(22,9))]:
        @lru_cache(None)
        def best(start,p,q,modified):
            value=F();witness=()
            for i in range(start,len(table)):
                z=table[i]
                if z['p']>p or z['q']>q or int(z['modified'])>modified:continue
                assert z['p']+z['q']>=F(1,2)
                v,w=best(i,p-z['p'],q-z['q'],modified-int(z['modified']))
                v+=z['loss']
                if v>value:value,witness=v,(i,)+w
            return value,witness
        value,w=best(0,p,q,1);assert value==target
        rows.append(dict(endpoint_budgets=list(map(str,(p,q))),maximum_loss=str(value),witness=w,states=best.cache_info().currsize))
    return rows

def all_rank_tails():
    rows=[]
    # Each row records rank r=start+t, fixed-coordinate energy and minimum
    # height-square. The target subtraction has nonnegative constant and
    # positive slope (the D18 constant is exactly zero).
    for name,start,shift,den,height2,modified,om,retained in[
        ('B',8,2,4,2/B**2,False,['0','r'],[['1'],['r-1']]),
        ('B',8,3,4,2/B**2,True,['r-1','0'],[['r-2','r'],['1']]),
        ('C',8,2,2,F(1),False,['0','r'],[['1'],['r-1']]),
        ('C',8,3,2,F(1),True,['1','r'],[['0','2'],['r-1']]),
        ('D',18,2,4,F(1),False,['0','r'],[['1'],['r-1']])]:
        margin=[height2*F(start-shift,den)-F(2*(start-int(modified)),9),height2/F(den)-F(2,9)]
        assert min(margin)>=0 and margin[1]>0
        guards=[]
        for r in range(start,start+3):
            def node(s):return r if s=='r'else r-1 if s=='r-1'else r-2 if s=='r-2'else int(s)
            ds,R,potential=classical_model(name,r);oi=list(map(node,om));allowed=[]
            for o,ids in zip(oi,retained):
                actual=[k for k in range(r+1)if k!=o and ds[k]*R(o,k)<=CAP]
                assert actual==sorted(map(node,ids));allowed.append(actual)
            full=[ids+[o]for ids,o in zip(allowed,oi)]
            begin=2 if name=='C'and modified else 1
            end=r-2 if name=='B'and modified else r-1
            expected=F(1)if name=='C'else F(1,2)
            for i,j in product(*full):
                assert all(abs(potential(i)[k]-potential(j)[k])==expected for k in range(begin,end))
            assert F(end-begin,2 if name=='C'else 4)==F(r-shift,den)
            guards.append(dict(rank=r,complete_diagonal_domains=allowed,fixed_coordinates=list(range(begin,end))))
        rows.append(dict(type=name,modified=modified,rank_start=start,omissions=om,possible_retained_nodes=retained,
            fixed_coordinate_energy=dict(numerator=['-'+str(shift),'r'],denominator=den),height_square=str(height2),
            W_minus_target_at_start_and_slope=list(map(str,margin)),literal_boundary_guards=guards))
    A=[]
    for start,eps in((18,0),(19,1)):
        m=[F(start),F(2)];square=H.pmul(m,m)
        # At cost<=3/2, only one neighbor of each omission may be retained.
        assert 2*F(start-2,start)>CAP and 2*F(start-1,start)>CAP
        # Singleton full supports: minimum cyclic distance h-2. With any
        # two-node full support, minheight>=2/B and the linear loss is at most
        # 2(m+eps)/m. These are separate cases, not the old W+C estimate.
        point=H.padd(H.padd(square,H.pscale(m,8)),[-9*(4+eps)**2])
        mixed=H.padd(H.padd(H.pscale(square,9-2*B**2),H.pscale(m,2*B**2-72)),[-81*eps])
        assert min(point)>=0 and point[0]>0 and min(mixed)>=0 and mixed[0]>0
        A.append(dict(cycle_size=[str(x)for x in m],parity=eps,possible_retained_offsets=[-1,1],
            maximum_retained_per_row=1,point_margin_numerator=list(map(str,point)),point_denominator='36*m',
            mixed_margin_numerator=list(map(str,mixed)),mixed_denominator='9*B^2*m'))
    return dict(BCD=rows,A=A,scope='For B/C/D the listed support reductions follow directly from delta_k^2 R(i,k)<=3/2. Fixed physical coordinates of every convex-hull point give the energy bounds. The new A tail separates singleton point rows from a row with positive omitted and retained entries. Rotations/fork swaps cover avoid-A/D and ordinary choices avoiding any other node, with a weaker ordinary target.')

def eliminate_same_row_E7():
    A=H.star([1,3,3]);ds,di,R,c=H.model(A,[1]*8);out=[]
    for side in(0,1):
        om=(3,6);T=F(14,9)
        for own_omitted in(False,True):
            z=(4,7)+((om[side],)if own_omitted else())
            U,V=(z,(om[1],))if side==0 else((om[0],),z)
            e,dual=closest(U,V,R);h=height(ds,U,V);lower=h*e+F(3,2)
            assert lower>T
            out.append(dict(side=side,own_omitted=own_omitted,left=U,right=V,energy=str(e),height=str(h),
                unit_margin=str(lower-T),dual=dual))
    return out

def failure_floor_census():
    old=json.loads((ROOT/'results/fresh_nikulin_height_charge_20261015_review.json').read_text())
    out=[];excluded=0
    for row in old['possible_defects']:
        name=row['type'];f=row['forbidden'];A=H.star([1,3,3]if name=='E7'else[1,2,5])
        ds,di,R,cross=H.model(A,[1]*len(A));om=H.choices(name,len(A)-1,f)
        nodes=row['nodes']
        if len({side for side,k in nodes})==1:
            assert name=='E7'and f is None and sorted(k for side,k in nodes)==[4,7]
            excluded+=1;continue
        costs=[sum((di[om[side]][k]for s,k in nodes if s==side),F())for side in(0,1)]
        assert sum(costs)==F(row['cost'])
        target=F(2*(len(A)-1-int(f is not None)),9)
        if name=='E8'and f is None:
            assert costs==[F(8,9),F(2,3)]and target==F(16,9)
        else:assert costs in([F(1,2),F(1)],[F(1),F(1,2)])and target==F(14,9)
        out.append(dict(type=name,forbidden=f,nodes=nodes,cost_floors=list(map(str,costs)),target=str(target)))
    assert excluded==2 and len(out)==13
    return out

def sharp_local_loss():
    ds,di,R,c=H.model(H.star([1,3,3]),[1]*8)
    assert ds[7]==1 and di[3][7]==1 and di[6][7]==F(1,2)
    bds,bR,p=classical_model('B',2)
    assert bds[0]==1 and bds[0]*bR(2,0)==F(1,2)
    loss=2*F(14,9)+2*F(4,9)
    assert loss==4 and 2*F(1)==2 and 2*F(1,2)+2*F(1,2)==2
    eds,edi,eR,ec=H.model(H.star([2,2,2]),[1]*7)
    assert eds[2]==eds[4]==1 and eR(2,4)==F(4,3)
    assert 2*7+2*2+1+24*6==163
    return dict(two_E7=dict(omissions=[3,6],both_rows_at=7,costs=['1','1/2'],W='0',target='14/9'),
        two_B2=dict(omissions=[0,2],both_rows_at=0,costs=['0','1/2'],W='0',target='4/9'),
        M='1',endpoint_cost_totals=['2','2'],loss='4',
        rank_padding='One avoid-A1 with zero target/loss and24 ordinary E6 with omitted point rows and W=T gives total rank163.',
        scope='Sharpness of this local numerical charge domain only, not an ambient Coxeter polytope or a claim of optimal dimension bound.')

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    records,table,counts=finite_envelope();pack=packing(table);tails=all_rank_tails();same=eliminate_same_row_E7();floors=failure_floor_census();sharp=sharp_local_loss()
    assert F(16,9)+F(20,9)==F(14,9)+F(22,9)==4
    assert(F(2*163,9)-F(2,9)-4)/2-2==14
    assert B*F(163-43,9)==14
    paths=[f'scripts/{TAG}.py']
    for stem in('fresh_nikulin_defect_frontier_20261015','fresh_nikulin_defect_frontier_20261015_review',
                 'fresh_nikulin_height_charge_20261015','fresh_nikulin_height_charge_20261015_review'):
        paths +=[f'scripts/{stem}.py',f'results/{stem}.json']
    for p in paths:
        if p.endswith('.json'):
            for z,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert sha(z)==h,(p,z)
    return dict(status='PASS_EXACT_ORDERED_RESIDUAL_BUDGET_CANDIDATE163_PENDING_INDEPENDENT_REVIEW',
        sources={p:sha(p)for p in paths},height_ratio_limit=str(B),residual_endpoint_cap=str(CAP),
        counts=counts,finite_models=records,pareto_table=[{k:str(v)if isinstance(v,F)else v for k,v in z.items()}for z in table],
        ordered_budget_packing=pack,all_rank_tails=tails,same_row_E7_exclusions=same,remaining_failure_floors=floors,sharp_local_loss=sharp,
        aggregate=dict(no_bad_factor_loss='4',ordinary_E8=dict(cost_floor=['8/9','2/3'],target='16/9',remaining=['10/9','4/3'],residual_loss='20/9'),
            ordinary_E7_or_avoid_E8=dict(cost_floors=[['1/2','1'],['1','1/2']],target='14/9',residual_loss='22/9'),
            total_loss='4',near_pairing='M*(N-37)/9',far_pairing='M*(N-43)/9',near_at_N163='14',far_at_N163_at_split='14',
            ideal_branch_dimension_bound=163,compact_branch_inherited_bound=143,combined_candidate_bound=163),
        scope='All-rank continuous residual support envelope with ordered endpoint costs. Finite cases carry exact rational convex KKT certificates and outward radical height bounds; no finite-label cutoff. Published geometric narrow-facet/endpoint/common-height and four-thirds inputs inherited from frozen independently reviewed164 chain. Internal research candidate, pending independent review; no publication-priority or formal-verification claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status']);print(result['counts'])
