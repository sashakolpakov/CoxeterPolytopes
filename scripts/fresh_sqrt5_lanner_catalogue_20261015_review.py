#!/usr/bin/env python3
"""Independent Qsqrt5 Lanner replay: support masks and integer biquadratics.

No producer/project imports. No finite Dynkin classifier or Pruefer generator.
Gram matrices are doubled to have integer coefficients in1,sqrt2,sqrt5,sqrt10.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,combinations_with_replacement,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt5_lanner_catalogue_20261015';TAG=BASE+'_review'
Z=(0,0,0,0);ONE=(1,0,0,0)
MAG={2:Z,3:(2,0,0,0),4:(0,2,0,0),5:(1,0,1,0)}

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def mul(a,b):
    c=[0]*4
    for i,x in enumerate(a):
        if not x:continue
        for j,y in enumerate(b):
            if y:c[i^j]+=x*y*(2 if i&j&1 else 1)*(5 if i&j&2 else 1)
    return tuple(c)
def qsign(a,b):
    if b==0:return (a>0)-(a<0)
    if a==0 or (a>0)==(b>0):return (b>0)-(b<0)
    return ((a>0)-(a<0))*((a*a-2*b*b>0)-(a*a-2*b*b<0))
def sign(x):
    a,b,c,d=x;su,sv=qsign(a,b),qsign(c,d)
    if not su:return sv
    if not sv or su==sv:return su
    # (a+b sqrt2)^2 -5(c+d sqrt2)^2 is inQsqrt2.
    return su*qsign(a*a+2*b*b-5*c*c-10*d*d,2*a*b-10*c*d)
def conjugate(x):return x[0],x[1],-x[2],-x[3]

def gram(n,labels):
    A=[[(4,0,0,0) if i==j else Z for j in range(n)] for i in range(n)]
    for (i,j),lab in zip(combinations(range(n),2),labels):A[i][j]=A[j][i]=neg(MAG[lab])
    return A

def determinant(A):
    n=len(A);out=Z
    for p in permutations(range(n)):
        if any(A[i][p[i]]==Z for i in range(n)):continue
        z=ONE
        for i in range(n):z=mul(z,A[i][p[i]])
        parity=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        out=add(out,neg(z) if parity%2 else z)
    return out

def induced(n,labels,ids):
    d=dict(zip(combinations(range(n),2),labels))
    return tuple(d[tuple(sorted((i,j)))] for i,j in combinations(ids,2))

@lru_cache(None)
def det(n,labels):return determinant(gram(n,labels))

@lru_cache(None)
def positive(n,labels,other=False):
    if n<=1:return True
    prefix=induced(n,labels,tuple(range(n-1)))
    if not positive(n-1,prefix,other):return False
    value=det(n,labels)
    return sign(conjugate(value) if other else value)>0

def phi_encoding(value,n):
    a,b,c,d=value;scale=2**n
    return [F(a-c,scale),F(2*c,scale),F(b-d,scale),F(2*d,scale)]

def canonical(n,labels):return min(induced(n,labels,p) for p in permutations(range(n)))

def supports(n):
    pairs=list(combinations(range(n),2));out=[]
    for mask in range(1<<len(pairs)):
        if mask.bit_count() not in (n-1,n):continue
        es=[e for i,e in enumerate(pairs) if mask>>i&1];adj={i:set() for i in range(n)}
        for a,b in es:adj[a].add(b);adj[b].add(a)
        seen={0};queue=[0]
        for a in queue:
            for b in adj[a]:
                if b not in seen:seen.add(b);queue.append(b)
        if len(seen)<n:continue
        if len(es)==n-1:kind='tree'
        elif all(len(adj[i])==2 for i in range(n)):kind='cycle'
        else:continue
        out.append((kind,es))
    assert Counter(k for k,e in out)==Counter(tree=n**(n-2),cycle=3 if n==4 else 12)
    return out

def higher(n,producer):
    counts=Counter();geom={};accepted={};rejected={}
    for kind,edges in supports(n):
        for word in product((3,4,5),repeat=len(edges)):
            counts['tested']+=1;counts[kind+'_tested']+=1
            e=dict(zip(edges,word));labels=tuple(e.get(p,2) for p in combinations(range(n),2))
            deleted=[induced(n,labels,tuple(i for i in range(n) if i!=j)) for j in range(n)]
            if not all(positive(n-1,x) for x in deleted):counts['proper_not_finite']+=1;continue
            value=det(n,labels)
            if sign(value)>=0:counts['not_indefinite']+=1;continue
            key=canonical(n,labels);geom[key]=phi_encoding(value,n)
            if kind=='cycle' and word.count(4)%2:
                counts['cycle_field_failure']+=1;rejected[key]='cycle outside Qsqrt5';continue
            assert value[1]==value[3]==0
            cv=conjugate(value)
            if sign(cv)<=0 or not all(positive(n-1,x,True) for x in deleted):
                counts['conjugate_failure']+=1;rejected[key]='conjugate not positive definite';continue
            counts['accepted_labelled']+=1
            accepted[key]=(phi_encoding(value,n)[:2],phi_encoding(cv,n)[:2])
    assert dict(counts)==producer['counts'],(n,counts,producer['counts'])
    assert geom=={tuple(r['labels']):list(map(F,r['determinant'])) for r in producer['geometric_classes']}
    assert rejected=={tuple(r['labels']):r['reason'] for r in producer['rejected_classes']}
    assert accepted=={tuple(r['labels']):(list(map(F,r['determinant'])),list(map(F,r['conjugate_determinant']))) for r in producer['accepted_classes']}
    return dict(order=n,counts=dict(counts),support_masks=len(supports(n)),
                accepted_classes=[dict(labels=list(k),determinant=list(map(str,v[0])),conjugate_determinant=list(map(str,v[1]))) for k,v in sorted(accepted.items())])

def qa(x,y):return x[0]+y[0],x[1]+y[1]
def qm(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def qs(x):
    a,b=x;return sign((2*a+b,0,b,0))
def qc(x):return x[0]+x[1],-x[1]

def triangles(producer):
    square={2:(0,0),3:(1,0),4:(2,0),5:(1,1),6:(3,0),10:(2,1)}
    rows=[]
    for labels in combinations_with_replacement(square,3):
        rec=dict(labels=list(labels))
        if sum(F(1,m) for m in labels)>=1:rec['status']='not hyperbolic'
        elif 2 not in labels and any(labels.count(m)%2 for m in (4,6,10)):rec['status']='cycle outside Qsqrt5'
        else:
            cycle=(0,0) if 2 in labels else (1,0)
            if 2 not in labels:
                for _ in range(labels.count(5)):cycle=qm(cycle,(0,1))
                for m in (4,6,10):
                    for _ in range(labels.count(m)//2):cycle=qm(cycle,square[m])
            total=cycle
            for m in labels:total=qa(total,square[m])
            value=(8-2*total[0],-2*total[1]);cv=qc(value)
            assert qs(value)<0
            assert all(qs(qa((4,0),tuple(-x for x in qc(square[m]))))>0 for m in labels)
            rec.update(cycle_product=list(cycle),determinant=list(value),conjugate_determinant=list(cv),
                       status='accepted' if qs(cv)>0 else 'conjugate not positive definite')
        rows.append(rec)
    assert rows==producer['records']
    assert [r for r in rows if r['status']=='accepted']==producer['accepted']
    return dict(counts=dict(Counter(r['status'] for r in rows)),accepted=[r for r in rows if r['status']=='accepted'],tested=len(rows))

def run():
    data=json.loads((ROOT/'results'/f'{BASE}.json').read_text())
    assert data['finite_alphabet']==[2,3,4,5,6,10]
    # Check proper positivity independently for ALL rank3/4 label matrices,
    # without the producer's finite-classification function.
    pd_counts=[]
    for n in (3,4):
        number=sum(positive(n,word) for word in product((2,3,4,5),repeat=n*(n-1)//2))
        expected=next(r for r in data['finite_type_guards'] if r['rank']==n)
        assert number==expected['positive_patterns']
        pd_counts.append(dict(rank=n,tested=4**(n*(n-1)//2),positive=number))
    tr=triangles(data['order_three']);large=[higher(r['order'],r) for r in data['higher_orders']]
    assert [len(r['accepted_classes']) for r in large]==[6,4] and len(tr['accepted'])==12
    files=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json']
    return dict(status='PASS_INDEPENDENT_COMPLETE_QSQRT5_LANNER_CATALOGUE_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                independent_arithmetic='Integer basis1,sqrt2,sqrt5,sqrt10; doubled Gram; permutation determinants; Sylvester minors; no producer imports',
                independent_supports='All connected edge masks with n-1 edges, or n edges and degree2 everywhere; no Pruefer generator',
                finite_PD_checks=pd_counts,triangles=tr,higher_orders=large,
                dotted_scope='Symbolic complete parameter t inZphi, t>4,0<t_conjugate<4; not finitely enumerated',
                scope='Catalogue only. Full ambient attachments, quadrilateral classification and any new Qsqrt5 dimension bound are outside this review.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['finite_PD_checks']);print([(r['order'],r['counts']) for r in out['higher_orders']])
