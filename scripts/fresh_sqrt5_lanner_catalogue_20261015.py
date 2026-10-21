#!/usr/bin/env python3
"""Complete arithmetic Qsqrt5 Lanner inputs, including all dotted parameters.

Higher-order completeness uses the proved tree-or-spanning-cycle reduction;
all labelled trees/cycles are checked, with no geometric-diagram whitelist.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,combinations_with_replacement,permutations,product
from math import isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_lanner_catalogue_20261015'
VINBERG='https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt'
NIKULIN='https://arxiv.org/pdf/0708.3991'
# Qphi and Qphi(sqrt2), in the bases(1,phi) and(1,phi,sqrt2,phi*sqrt2).
def qa(x,y):return x[0]+y[0],x[1]+y[1]
def qn(x):return -x[0],-x[1]
def qm(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def qs(x,k):return x[0]*k,x[1]*k
def qc(x):return x[0]+x[1],-x[1]
def qsign(x):
    a,b=x;u=2*a+b
    if not b:return (u>0)-(u<0)
    if not u:return (b>0)-(b<0)
    if (u>0)==(b>0):return (u>0)-(u<0)
    z=u*u-5*b*b;assert z
    return ((u>0)-(u<0))*((z>0)-(z<0))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def mul(x,y):
    a,b=x[:2],x[2:];c,d=y[:2],y[2:]
    return qa(qm(a,c),qs(qm(b,d),2))+qa(qm(a,d),qm(b,c))
def conj(x):return qc(x[:2])+qc(x[2:])
def sign(x):
    a,b=x[:2],x[2:];sa,sb=qsign(a),qsign(b)
    if not sa:return sb
    if not sb or sa==sb:return sa
    z=qa(qm(a,a),qn(qs(qm(b,b),2)));assert z!=(0,0)
    return sa*qsign(z)
Z=(0,0,0,0);ONE=(1,0,0,0)
MAG={2:Z,3:ONE,4:(0,0,1,0),5:(0,1,0,0)}
def determinant(M):
    n=len(M)
    @lru_cache(None)
    def rec(row,cols):
        if not cols:return ONE
        out=Z
        for k,j in enumerate(cols):
            if M[row][j]==Z:continue
            z=mul(M[row][j],rec(row+1,cols[:k]+cols[k+1:]))
            out=add(out,neg(z)if k%2 else z)
        return out
    return rec(0,tuple(range(n)))
def gram(n,labels):
    M=[[(2,0,0,0)if i==j else Z for j in range(n)]for i in range(n)]
    for(i,j),m in zip(combinations(range(n),2),labels):M[i][j]=M[j][i]=neg(MAG[m])
    return M
def induced(n,labels,ids):
    e=dict(zip(combinations(range(n),2),labels))
    return tuple(e[tuple(sorted((i,j)))]for i,j in combinations(ids,2))
def canonical(n,labels):return min(induced(n,labels,p)for p in permutations(range(n)))
def components(adj):
    left=set(adj);out=[]
    while left:
        todo=[min(left)];seen=set()
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo.extend(adj[i])
        left-=seen;out.append(seen)
    return out
@lru_cache(None)
def finite(n,labels):
    """Complete finite diagrams of rank<=4 on labels2,3,4,5."""
    assert n<=4
    adj={i:{}for i in range(n)}
    for(i,j),m in zip(combinations(range(n),2),labels):
        if m!=2:adj[i][j]=m;adj[j][i]=m
    for ids in components(adj):
        r=len(ids)
        if r<=2:continue
        if sum(len(adj[i])for i in ids)!=2*(r-1):return False
        if max(len(adj[i])for i in ids)>3:return False
        branching=[i for i in ids if len(adj[i])==3]
        if branching:
            if r!=4 or any(m!=3 for i in ids for m in adj[i].values()):return False
            continue
        leaf=min(i for i in ids if len(adj[i])==1);previous=None;word=[]
        while True:
            nxt=[j for j in adj[leaf]if j!=previous]
            if not nxt:break
            assert len(nxt)==1
            j=nxt[0];word.append(adj[leaf][j]);previous,leaf=leaf,j
        if all(m==3 for m in word):continue
        if word.count(4)==1 and all(m in(3,4)for m in word):
            if word.index(4)in(0,len(word)-1)or word==[3,4,3]:continue
        if word.count(5)==1 and all(m in(3,5)for m in word)and word.index(5)in(0,len(word)-1):continue
        return False
    return True
def prufer_tree(n,word):
    degree=[1]*n
    for x in word:degree[x]+=1
    edges=[]
    for x in word:
        leaf=next(i for i in range(n)if degree[i]==1)
        edges.append(tuple(sorted((leaf,x))));degree[leaf]-=1;degree[x]-=1
    ends=[i for i in range(n)if degree[i]==1];edges.append(tuple(ends))
    return tuple(sorted(edges))
def topologies(n):
    trees={prufer_tree(n,w)for w in product(range(n),repeat=n-2)}
    assert len(trees)==n**(n-2)
    cycles=set()
    for p in permutations(range(1,n)):
        order=(0,)+p;cycles.add(tuple(sorted(tuple(sorted((order[i],order[(i+1)%n])))for i in range(n))))
    assert len(cycles)==(3 if n==4 else 12)
    return [('tree',e)for e in sorted(trees)]+[('cycle',e)for e in sorted(cycles)]
def higher(n):
    counts=Counter();accepted={};geometric={};rejected={}
    for kind,edges in topologies(n):
        for word in product((3,4,5),repeat=len(edges)):
            counts['tested']+=1;counts[kind+'_tested']+=1
            d=dict(zip(edges,word));labels=tuple(d.get(e,2)for e in combinations(range(n),2))
            if not all(finite(n-1,induced(n,labels,tuple(i for i in range(n)if i!=j)))for j in range(n)):
                counts['proper_not_finite']+=1;continue
            M=gram(n,labels);det=determinant(M)
            if sign(det)>=0:counts['not_indefinite']+=1;continue
            key=canonical(n,labels);geometric[key]=det
            # The only cycle is the spanning one; its product is inQphi iff
            # it contains an even number ofsqrt2 bonds. Trees have no cycles.
            if kind=='cycle'and word.count(4)%2:
                counts['cycle_field_failure']+=1;rejected[key]='cycle outside Qsqrt5';continue
            assert det[2:]==(0,0)
            cp=[[(conj(x))for x in r]for r in M]
            # Check every principal minor, independently of finite-type list.
            proper=[];ok=True
            for size in range(1,n):
                for ids in combinations(range(n),size):
                    a=determinant([[M[i][j]for j in ids]for i in ids])
                    b=determinant([[cp[i][j]for j in ids]for i in ids])
                    assert sign(a)>0
                    proper.append((ids,a,b))
                    if sign(b)<=0:ok=False
            cdet=determinant(cp)
            if sign(cdet)<=0 or not ok:
                counts['conjugate_failure']+=1;rejected[key]='conjugate not positive definite';continue
            counts['accepted_labelled']+=1
            accepted[key]=dict(labels=key,edges=[[i,j,m]for(i,j),m in zip(combinations(range(n),2),key)if m!=2],
                determinant=det[:2],conjugate_determinant=cdet[:2],proper_principal_minor_tests=len(proper),topology=kind)
    return dict(order=n,topology_counts=dict(Counter(k for k,e in topologies(n))),counts=dict(counts),
        geometric_classes=[dict(labels=k,determinant=v)for k,v in sorted(geometric.items())],
        rejected_classes=[dict(labels=k,reason=v)for k,v in sorted(rejected.items())],
        accepted_classes=[v for k,v in sorted(accepted.items())])
def triangles():
    squares={2:(0,0),3:(1,0),4:(2,0),5:(1,1),6:(3,0),10:(2,1)}
    rows=[];accepted=[]
    for labels in combinations_with_replacement(squares,3):
        if sum(F(1,m)for m in labels)>=1:
            rows.append(dict(labels=labels,status='not hyperbolic'));continue
        if 2 in labels:cycle=(0,0)
        else:
            if any(labels.count(m)%2 for m in(4,6,10)):
                rows.append(dict(labels=labels,status='cycle outside Qsqrt5'));continue
            cycle=(1,0)
            for _ in range(labels.count(5)):cycle=qm(cycle,(0,1))
            for m in(4,6,10):
                for _ in range(labels.count(m)//2):cycle=qm(cycle,squares[m])
        total=(0,0)
        for m in labels:total=qa(total,squares[m])
        det=qa((8,0),qn(qs(qa(total,cycle),2)))
        assert qsign(det)<0
        assert all(qsign(qa((4,0),qn(qc(squares[m]))))>0 for m in labels)
        row=dict(labels=labels,cycle_product=cycle,determinant=det,conjugate_determinant=qc(det),
            status='accepted'if qsign(qc(det))>0 else'conjugate not positive definite')
        rows.append(row)
        if row['status']=='accepted':accepted.append(row)
    assert len(rows)==56
    return dict(records=rows,accepted=accepted,counts=dict(Counter(r['status']for r in rows)))
def run():
    # Independently verify finite-type recognition over every rank<=4 label
    # matrix: every classified-finite matrix has all positive principal minors.
    finite_guards=[]
    for n in(3,4):
        checked=0
        tested=0
        for labels in product((2,3,4,5),repeat=n*(n-1)//2):
            tested+=1;M=gram(n,labels)
            exact_pd=all(sign(determinant([[M[i][j]for j in range(k)]for i in range(k)]))>0 for k in range(1,n+1))
            assert finite(n,labels)==exact_pd
            checked+=exact_pd
        finite_guards.append(dict(rank=n,tested_patterns=tested,positive_patterns=checked))
    tri=triangles();large=[higher(n)for n in(4,5)]
    assert [len(r['accepted_classes'])for r in large]==[6,4]
    dotted=(3,1);assert qsign(qa(dotted,(-4,0)))>0 and qsign(qc(dotted))>0 and qsign(qa((4,0),qn(qc(dotted))))>0
    assert qm((2,1),qc((2,1)))==(5,0)
    # squareclasses2,3,6 are not squares inQsqrt5: ab=0 would force
    # rational a²=d or5b²=d; their rational square tests fail.
    squareclass_checks=[]
    for d in(2,3,6):
        for z in(F(d),F(d,5)):
            assert isqrt(z.numerator)**2!=z.numerator or isqrt(z.denominator)**2!=z.denominator
        squareclass_checks.append(dict(radicand=d,rational_obstruction=[d,f'{d}/5']))
    for a,b in product(range(2),repeat=2):
        norm=5*(2**a*3**b)**2
        assert isqrt(norm)**2!=norm
    return json.loads(json.dumps(dict(status='PASS_COMPLETE_QSQRT5_ARITHMETIC_LANNER_CATALOGUE',
        sources={f'scripts/{TAG}.py':sha256((ROOT/f'scripts/{TAG}.py').read_bytes()).hexdigest()},
        field_basis='1,phi, with phi²=phi+1; higher determinant extension adds sqrt2',
        published_inputs=dict(Vinberg=VINBERG,Nikulin_arithmetic_criterion=NIKULIN,Lanner_maximum_order=5),
        finite_alphabet=[2,3,4,5,6,10],order_two=dict(parameter='t inZ[phi], t>4, 0<t_conjugate<4',example=dotted,no_dotted_cutoff=True),
        triangle_squareclass_checks=squareclass_checks,order_three=tri,higher_orders=large,finite_type_guards=finite_guards,
        scope='Complete arithmetic-compatible Lanner input overQsqrt5: parametric order2, all56 triangle triples, all labelled trees/spanningcycles oforders4/5. No ambient quadrilateral classification or dimension improvement.' )))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved Qsqrt5 catalogue replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
    print('triangles',len(d['order_three']['accepted']),'higher',[len(r['accepted_classes'])for r in d['higher_orders']])
