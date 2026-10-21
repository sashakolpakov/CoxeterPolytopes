#!/usr/bin/env python3
"""Independent physical-Gram/Fraction audit of the ordered-budget163 proof.

No producer or project arithmetic imports. Rebuilds the finite support domains,
checks every convex dual, recomputes the ordered packing by a distinct dynamic
program, and proves the all-rank coordinate/energy inequalities symbolically.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from math import isqrt
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_nikulin_two_budget_frontier_20261015';TAG=BASE+'_review'
B=F(21,20);CAP=F(3,2)

def dot(x,y):return sum((a*b for a,b in zip(x,y)),F())
def mv(A,x):return[dot(row,x)for row in A]
def inverse(A):
    n=len(A);M=[list(row)+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
    for k in range(n):
        i=next(i for i in range(k,n)if M[i][k]);M[k],M[i]=M[i],M[k]
        z=M[k][k];M[k]=[x/z for x in M[k]]
        for i in range(n):
            if i!=k:
                z=M[i][k];M[i]=[x-z*y for x,y in zip(M[i],M[k])]
    return[row[n:]for row in M]
def pd(A):
    A=[list(row)for row in A]
    while A:
        d=A[0][0];assert d>0
        A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,len(A))]for i in range(1,len(A))]
def subsets(ids):
    for k in range(len(ids)+1):yield from combinations(ids,k)
def frac(x):return F(str(x))
def star(arms):
    G=2*s.eye(sum(arms)+1);start=1
    for length in arms:
        path=[0]+list(range(start,start+length));start+=length
        for a,b in zip(path,path[1:]):G[a,b]=G[b,a]=-1
    return G

@lru_cache(None)
def model(name,r):
    G=2*s.eye(r+1)
    def edge(i,j,w=1):G[i,j]=G[j,i]=-w
    if name=='A':
        if r==1:edge(0,1,2)
        else:
            for i in range(r+1):edge(i,(i+1)%(r+1))
    elif name=='B' or name=='B2modified':
        if r==2:edge(0,2,s.sqrt(2));edge(1,2,s.sqrt(2))
        else:
            edge(0,2);edge(1,2)
            for k in range(2,r):edge(k,k+1,s.sqrt(2)if k==r-1 else 1)
    elif name=='C':
        for k in range(r):edge(k,k+1,s.sqrt(2)if k in(0,r-1)else 1)
    elif name=='D':
        edge(0,2);edge(1,2)
        for k in range(2,r-2):edge(k,k+1)
        edge(r-2,r-1);edge(r-2,r)
    elif name.startswith('E'):G=star({'E6':[2,2,2],'E7':[1,3,3],'E8':[1,2,5]}[name])
    elif name=='F4':
        for k in range(4):edge(k,k+1,s.sqrt(2)if k==2 else 1)
    elif name=='G2':edge(0,1);edge(1,2,s.sqrt(3))
    else:raise AssertionError(name)
    kernel=G.nullspace();assert len(kernel)==1
    d=kernel[0]
    if d[0]<0:d=-d
    d=s.simplify(d/min(d));assert min(d)==1 and G*d==s.zeros(r+1,1)
    # Rational weighted Laplacian; its grounded inverse is an independent
    # Euclidean Gram for normalized point potentials, based at node0.
    L=[[frac(s.simplify(d[i]*G[i,j]*d[j]))for j in range(r+1)]for i in range(r+1)]
    assert all(sum(row)==0 for row in L)
    ground=[row[1:]for row in L[1:]];pd(ground);J=inverse(ground)
    K=[[F()]*(r+1)]+[[F()]+row for row in J]
    R=lambda i,j:K[i][i]+K[j][j]-2*K[i][j]
    ds=tuple(frac(s.simplify(x*x))for x in d)
    return dict(G=G,kernel=d,squares=ds,K=K,R=R)

def omissions(name,r,f):
    if name=='A':return(0,(r+1)//2)
    if name=='B':return(r-1,0)if f else(0,r)
    if name=='C':return(1,r)if f else(0,r)
    if name=='D':return(0,r)
    if name=='B2modified':return(0,1)
    if f is None:return{'E6':(2,4),'E7':(3,6),'E8':(1,6),'F4':(0,4),'G2':(0,2)}[name]
    if name=='E6':return tuple(k for k in(2,4,6)if k!=f)[:2]
    if name=='E7':return(2,5)if f in(4,7)else(4,7)
    if name=='E8':return(1,6)if f in(3,7)else(3,7)
    if name=='F4':return(0,4)if f in(1,3)else(1,3)
    return tuple(k for k in range(3)if k!=f)

def descriptors():
    out=[]
    for r in range(1,17):out.append(('A',r,False,None))
    for typ,end in(('B',7),('C',7),('D',17)):
        for r in range(4 if typ=='D'else 2,end+1):
            out.append((typ,r,False,None))
            if typ=='C'or typ=='B'and r>=3:out.append((typ,r,True,None))
    for typ,r in(('E6',6),('E7',7),('E8',8),('F4',4),('G2',2)):
        for f in[None]+list(range(r+1)):out.append((typ,r,f is not None,f))
    out.append(('B2modified',2,True,0));assert len(out)==91
    return out

def radical_floor(q):
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    if a*a==q.numerator and b*b==q.denominator:return F(a,b)
    D=2**64;t=isqrt(q.numerator*D*D//q.denominator)
    assert F(t,D)**2<=q<F(t+1,D)**2
    return F(t,D)
def height(ds,U,V):
    a,b=[sum((radical_floor(ds[k])for k in ids),F())for ids in(U,V)]
    return max(min(a,b)**2,max(a,b)**2/B**2)
def cost(m,o,ids):
    ds,R=m['squares'],m['R'];z=sum(ds[k]*R(o,k)for k in ids)
    for k,l in combinations(ids,2):
        cross=R(o,k)+R(o,l)-R(k,l);assert cross>=0
        z+=s.sqrt(s.Rational(ds[k]*ds[l]))*s.Rational(cross)
    return s.simplify(z)
def dual(m,U,V,raw):
    ids=U+V;sg=[1]*len(U)+[-1]*len(V);K=m['K']
    A=[[sg[i]*sg[j]*K[k][l]for j,l in enumerate(ids)]for i,k in enumerate(ids)]
    x=list(map(F,raw['weights']));lm=list(map(F,raw['multipliers']));assert len(x)==len(ids)and len(lm)==2
    assert min(x)>=0 and sum(x[:len(U)])==sum(x[len(U):])==1
    gx=mv(A,x);sl=[z+lm[int(i>=len(U))]for i,z in enumerate(gx)]
    assert sl==list(map(F,raw['slacks']))and min(sl)>=0 and all(a*b==0 for a,b in zip(x,sl))
    e=dot(x,gx);assert e>=0 and e==-sum(lm)
    return e

def finite_census(data):
    allrows=[];positive=[];counts=Counter();saved=data['finite_models'];assert len(saved)==91
    for index,(raw,desc)in enumerate(zip(saved,descriptors())):
        name,r,modified,f=desc;assert(name,r,modified,f)==(raw['name'],raw['rank'],raw['modified'],raw['forbidden'])
        m=model(name,r);ds,R=m['squares'],m['R'];om=omissions(name,r,modified if name in('B','C')else f)
        assert list(om)==raw['omissions']and list(map(str,ds))==raw['kernel_squares']
        T=F(2*(r-int(modified)),9);assert T==F(raw['target'])
        assert min(ds[o]for o in om)*R(*om)>=T
        domains=[]
        for o in om:
            ids=[k for k in range(r+1)if k!=o];G=m['G'].extract(ids,ids);J=s.simplify(G.inv())
            assert all(z>=0 for z in J)and all(J[k,k]>=s.Rational(1,2)for k in range(r))
            for i,k in enumerate(ids):
                for j,l in enumerate(ids):
                    v=s.sqrt(s.Rational(ds[k]*ds[l]))*s.Rational((R(o,k)+R(o,l)-R(k,l))/2)
                    assert s.simplify(J[i,j]-v)==0
            cheap=[k for k in ids if ds[k]*R(o,k)<=CAP];good=[]
            for group in subsets(cheap):
                z=cost(m,o,group)
                if z<=s.Rational(CAP):
                    assert z.is_Rational;good.append((group,frac(z)))
            domains.append(good)
        assert [[[list(ids),str(c)]for ids,c in row]for row in domains]==raw['domains']
        expected={}
        for(L,cp),(Q,cq)in product(*domains):
            if not L and not Q:continue
            counts['support_pairs']+=1
            for a,b in product((False,True),repeat=2):
                U=L+((om[0],)if a else());V=Q+((om[1],)if b else())
                if U and V:
                    assert(U,V)not in expected;expected[U,V]=(cp,cq)
        check=[]
        for ci,z in enumerate(raw['cases']):
            U,V=tuple(z['left']),tuple(z['right']);cp,cq=expected.pop((U,V));assert[cp,cq]==list(map(F,z['costs']))
            e=dual(m,U,V,z['dual']);h=height(ds,U,V);loss=T-h*e
            assert[e,h,loss]==[F(z[k])for k in('minimum_energy','height_lower','loss_upper')]
            if loss>0:
                positive.append((index,ci,modified,cp,cq,loss));counts['positive']+=1
            check.append(dict(left=U,right=V,energy=str(e),height=str(h),loss=str(loss)))
        assert not expected;counts['domains']+=len(check);counts['models']+=1
        allrows.append(dict(name=name,rank=r,modified=modified,forbidden=f,omissions=om,support_domain_sizes=list(map(len,domains)),verified_duals=len(check)))
    assert dict(counts)==dict(support_pairs=1116,positive=1395,domains=3630,models=91)
    # Derive the exact partial-order frontier without relying on saved order.
    unique={(mod,p,q,l)for _,_,mod,p,q,l in positive}
    frontier={z for z in unique if not any(y!=z and y[0]==z[0]and y[1]<=z[1]and y[2]<=z[2]and y[3]>=z[3]for y in unique)}
    table=data['pareto_table'];keys=[(z['modified'],F(z['p']),F(z['q']),F(z['loss']))for z in table]
    assert len(keys)==len(set(keys))==28 and set(keys)==frontier
    for mi,ci,mod,p,q,l in positive:
        j=saved[mi]['cases'][ci]['dominator'];a,b,c,d=keys[j]
        assert a==mod and b<=p and c<=q and d>=l
    for z in table:
        r=saved[z['model']];c=r['cases'][z['case']]
        assert c['left']==c['right']==z['left']==z['right']and len(c['left'])==1
        assert F(c['minimum_energy'])==0 and c['loss_upper']==z['loss']and c['costs']==[z['p'],z['q']]
    return dict(counts=dict(counts),models=allrows,pareto_records=len(keys)),keys

def ordered_pack(keys,data):
    rows=[]
    # Forward exact reachability of ordered cost pairs, allowing sequences.
    # This is different from the producer's nondecreasing multiset recursion.
    for P,Q,target in((F(10,9),F(4,3),F(20,9)),(F(1),F(3,2),F(22,9)),(F(3,2),F(1),F(22,9))):
        states={(F(),F(),0):F()};front=list(states)
        while front:
            a,b,k=front.pop();value=states[a,b,k]
            for modified,p,q,l in keys:
                assert p+q>=F(1,2)
                new=(a+p,b+q,k+int(modified))
                if new[0]>P or new[1]>Q or new[2]>1:continue
                if states.get(new,F(-1))<value+l:states[new]=value+l;front.append(new)
        maximum=max(states.values());assert maximum==target
        raw=next(z for z in data['ordered_budget_packing']if list(map(F,z['endpoint_budgets']))==[P,Q]);assert F(raw['maximum_loss'])==maximum
        witness=[keys[i]for i in raw['witness']]
        assert sum(z[0]for z in witness)<=1 and sum(z[1]for z in witness)<=P and sum(z[2]for z in witness)<=Q and sum(z[3]for z in witness)==maximum
        rows.append(dict(budgets=[str(P),str(Q)],maximum=str(maximum),reachable_cost_states=len(states)))
    return rows

def nonnegative_poly(expr,variables):
    p=s.Poly(s.expand(expr),*variables);assert min(p.coeffs())>=0
    return str(p.as_expr())

def tails(data):
    t,u=s.symbols('t u',nonnegative=True);out=[]
    # Exhaustive symbolic node classes from the literal prefix potentials.
    # Each expression is its inverse diagonal. Unlisted interior nodes are
    # represented by their distance d>=2, so their cost is at least2.
    table=[('B',False,8,2,4,2/B**2),('B',True,8,3,4,2/B**2),('C',False,8,2,2,F(1)),('C',True,8,3,2,F(1)),('D',False,18,2,4,F(1))]
    for typ,mod,start,shift,den,h in table:
        r=start+t;T=2*(r-int(mod))/9;E=(r-shift)/den;margin=s.expand(s.Rational(h)*E-T);nonnegative_poly(margin,[t])
        # Every omitted/retained convex-hull pair differs by the same fixed
        # coordinates, proving the E bound with arbitrary contact values.
        if typ=='B':
            excluded=[(r-1)/4]if mod else[r/4,r/2]
            allowed=[[F(1),F(1,2)],[F(1)]]if mod else[[F(1)],[F(1)]]
        elif typ=='C':excluded=[(r-1)/2]if mod else[r/2];allowed=[[F(1,2),F(1)],[F(1)]]if mod else[[F(1)],[F(1)]]
        else:excluded=[r/4];allowed=[[F(1)],[F(1)]]
        assert all(min(row)<=CAP and max(row)<=CAP for row in allowed)
        for z in excluded:assert s.expand(z-s.Rational(CAP)).subs(t,0)>0;nonnegative_poly(z-s.Rational(CAP),[t])
        assert 2+u-s.Rational(CAP)>0
        # Reconstruct physical boundary Grams and compare every inverse
        # diagonal and every coordinate distance, without producer helpers.
        boundaries=[]
        for rr in range(start,start+3):
            m=model(typ,rr);ds,R=m['squares'],m['R'];om=omissions(typ,rr,mod)
            actual=[[k for k in range(rr+1)if k!=o and ds[k]*R(o,k)<=CAP]for o in om]
            expected=([[rr-2,rr],[1]]if mod else[[1],[rr-1]])if typ=='B'else([[0,2],[rr-1]]if mod else[[1],[rr-1]])if typ=='C'else[[1],[rr-1]]
            assert actual==expected
            def v(k):
                if typ=='C':return[F(int(j<k))for j in range(rr)]
                if k==0:return[F()]*rr
                if k==1:return[F(1)]+[F()]*(rr-1)
                if typ=='D'and k==rr-1:return[F(1,2)]*(rr-1)+[F(-1,2)]
                return[F(int(j<k),2)for j in range(rr)]
            for i,j in product(range(rr+1),repeat=2):assert R(i,j)==sum((a-b)**2 for a,b in zip(v(i),v(j)))/(2 if typ=='C'else 1)
            full=[ids+[o]for ids,o in zip(actual,om)];begin=2 if typ=='C'and mod else 1;end=rr-2 if typ=='B'and mod else rr-1
            for i,j in product(*full):assert all(v(j)[k]-v(i)[k]==(F(-1,2)if typ=='B'and mod else F(1)if typ=='C'else F(1,2))for k in range(begin,end))
            assert F(end-begin,2 if typ=='C'else 4)==F(rr-shift,den)
            # For B one side's every allowed full-support node has delta>=sqrt2.
            if typ=='B':assert max(min(ds[i]for i in row)for row in full)>=2
            boundaries.append(dict(rank=rr,retained=actual))
        raw=next(z for z in data['all_rank_tails']['BCD']if z['type']==typ and z['modified']==mod)
        assert list(map(F,raw['W_minus_target_at_start_and_slope']))==[frac(margin.subs(t,0)),frac(s.diff(margin,t))]
        out.append(dict(type=typ,modified=mod,margin=str(margin),excluded_node_costs=list(map(str,excluded)),interior_cost='2+u',boundaries=boundaries))
    cyc=[]
    for start,eps in((18,0),(19,1)):
        m=start+2*t;h=(m-eps)/2;R=lambda z:z*(m-z)/m
        # Concavity reduces every cyclic distance>=2 to distance2.
        z=s.symbols('d',real=True)
        assert s.expand(z*(m-z)-2*(m-2)-(z-2)*(m-z-2))==0
        assert(2*R(1)-s.Rational(CAP)).subs(t,0)>0 and (R(2)-s.Rational(CAP)).subs(t,0)>0
        T=2*(m-1)/9;point=s.factor(R(h-2)-T);mixed=s.factor(4/s.Rational(B*B)*(R(h)-2*(m+eps)/m)-T)
        pointnum=s.expand(m*m+8*m-9*(4+eps)**2)
        mixednum=s.expand((9-2*s.Rational(B*B))*m*m+(2*s.Rational(B*B)-72)*m-81*eps)
        assert s.factor(point-pointnum/(36*m))==0 and s.factor(mixed-mixednum/(9*s.Rational(B*B)*m))==0
        nonnegative_poly(pointnum,[t]);nonnegative_poly(mixednum,[t]);assert pointnum.subs(t,0)>0 and mixednum.subs(t,0)>0
        # The mean may move at most one edge from each omission. Convex
        # expansion loses at most max |R(h+-1)-R(h)-R(1)| per endpoint.
        ell=(m+eps)/m
        for sign in(-1,1):assert s.factor(ell-(R(h)+R(1)-R(h+sign)))>=0
        cyc.append(dict(start=start,parity=eps,point_numerator=str(pointnum),mixed_numerator=str(mixednum),linear_loss_per_row=str(ell)))
    return dict(BCD=out,A=cyc,all_rank_scope='Literal prefix coordinate classes and cycle concavity/one-edge mean movement; finite boundary checks are illustrations, not the tail proof.')

def bad_factors(data):
    old=json.loads((ROOT/'results/fresh_nikulin_defect_frontier_20261015_review.json').read_text())
    oldrows=old['physical_exceptional']['defect_occurrences']
    m=model('E7',7);ds=m['squares'];same=[]
    for raw in data['same_row_E7_exclusions']:
        U,V=tuple(raw['left']),tuple(raw['right']);e=dual(m,U,V,raw['dual']);h=height(ds,U,V)
        assert e==F(raw['energy'])and h==F(raw['height'])and e*h+CAP-F(14,9)==F(raw['unit_margin'])>0
        same.append(dict(side=raw['side'],own_omitted=raw['own_omitted'],energy=str(e),height=str(h),margin=raw['unit_margin']))
    assert len(same)==4;floors=[];excluded=0
    for raw in oldrows:
        name,f,nodes=raw['type'],raw['forbidden'],raw['nodes'];m=model(name,int(name[1:]));om=omissions(name,int(name[1:]),f)
        if len(set(side for side,k in nodes))==1:
            assert name=='E7'and f is None and sorted(k for side,k in nodes)==[4,7];excluded+=1;continue
        c=[sum((m['squares'][k]*m['R'](om[side],k)for side,k in nodes if side==a),F())for a in(0,1)]
        z=dict(type=name,forbidden=f,nodes=nodes,cost_floors=list(map(str,c)),target=raw['target']);floors.append(z)
    assert excluded==2 and floors==data['remaining_failure_floors']and len(floors)==13
    assert F(16,9)+F(20,9)==F(14,9)+F(22,9)==4
    # Rebuild the sharp local witness, including a zero-target avoid-A1.
    E=model('E7',7);BB=model('B',2);EE=model('E6',6)
    costs=[2*E['squares'][7]*E['R'](o,7)+2*BB['squares'][0]*BB['R'](p,0)for o,p in zip((3,6),(0,2))]
    assert costs==[2,2]and 2*F(14,9)+2*F(4,9)==4
    assert EE['squares'][2]==EE['squares'][4]==1 and EE['R'](2,4)==F(4,3)
    assert 2*7+2*2+1+24*6==163
    return dict(same_row_exclusions=same,remaining_floors=floors,sharp_local_loss='4',sharp_rank_padding=163)

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    data=json.loads((ROOT/f'results/{BASE}.json').read_text())
    finite,keys=finite_census(data);pack=ordered_pack(keys,data);tail=tails(data);fail=bad_factors(data)
    t=s.symbols('t',nonnegative=True);near=s.expand((163+t-37)/9);far=s.expand(s.Rational(B)*(163+t-43)/9)
    assert near==14+t/9 and far==14+7*t/60
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json']
    for p in('fresh_nikulin_defect_frontier_20261015','fresh_nikulin_defect_frontier_20261015_review','fresh_nikulin_height_charge_20261015','fresh_nikulin_height_charge_20261015_review','fresh_nikulin_unit_charge_20261015'):
        paths +=[f'scripts/{p}.py',f'results/{p}.json']
    visited=set()
    def check(p):
        if p in visited:return
        visited.add(p)
        if p.endswith('.json'):
            for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert sha(q)==h,(p,q);check(q)
    for p in paths:check(p)
    return dict(status='PASS_INDEPENDENT_PHYSICAL_GRAMS_CONTINUOUS_ENVELOPE_ORDERED_PACKING_AND_163_SCOPE',sources={p:sha(p)for p in paths},recursive_manifest_files=len(visited),
        finite_continuous_envelope=finite,ordered_packing=pack,all_rank_tails=tail,actual_failing_factor_split=fail,
        global_threshold=dict(N='163+t',near=str(near),far=str(far),ideal_dimension=163,compact_inherited_dimension=143),
        scope='New arithmetic reconstructed without producer/project imports: physical affine Grams,91 complete support models,3630 convex KKT certificates,28-record dominance and ordered packing,all-rank tails and failing-factor split. Prior15-defect completeness and global geometric/narrow-facet/compact-branch inputs inherited from frozen independently reviewed164 chain. Research candidate; no novelty/current-best claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text());print('PASS saved independent163 review')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status']);print(data['finite_continuous_envelope']['counts'])
