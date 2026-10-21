#!/usr/bin/env python3
"""Exact principal-union field and signature filters. External whites stay open."""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt2_quad_exhaustion_20261015'
TAG=BASE+'_cores'
K=s.QQ.algebraic_field(s.sqrt(2+s.sqrt(2)),s.sqrt(3))
MINPOLY=[F(str(x)) for x in K.mod.to_list()]

def peval(cs,x):
    v=F(0)
    for c in cs:v=v*x+c
    return v

def root_interval(bits=160):
    # theta=sqrt(3)+sqrt(2+sqrt2) is the largest of its eight real conjugates.
    aa,bb=K.ext.minpoly.intervals(eps=s.Rational(1,10000))[-1][0]
    a,b=F(aa),F(bb)
    assert F(357,100)<a<b<F(359,100) and peval(MINPOLY,a)*peval(MINPOLY,b)<0
    for _ in range(bits):
        m=(a+b)/2
        if peval(MINPOLY,a)*peval(MINPOLY,m)<0:b=m
        else:a=m
    return a,b

INTERVAL=root_interval()

def enc(v):return list(map(str,v.to_list()))

def sign(v):
    if not v:return 0
    a,b=INTERVAL;lo=hi=F(0)
    for c in v.to_list():
        vals=(lo*a,lo*b,hi*a,hi*b);lo=min(vals)+F(str(c));hi=max(vals)+F(str(c))
    assert lo>0 or hi<0,('unresolved sign',enc(v),lo,hi)
    return 1 if lo>0 else -1

def matrix(row,conjugate=False):
    values={2:0,3:s.Rational(1,2),4:(-1 if conjugate else 1)*s.sqrt(2)/2,
            6:s.sqrt(3)/2,8:s.sqrt(2+(-1 if conjugate else 1)*s.sqrt(2))/2}
    n=row['nodes'];M=[[K.one if i==j else K.zero for j in range(n)] for i in range(n)]
    for (i,j),lab in zip(combinations(range(n),2),row['labels']):
        assert lab!=0
        M[i][j]=M[j][i]=K.from_sympy(s.sympify(-values[lab]))
    return DomainMatrix(M,(n,n),K)

def cyclic_field(row):
    n=row['nodes'];adj={i:[] for i in range(n)}
    for (a,b),lab in zip(combinations(range(n),2),row['labels']):
        if lab==2:continue
        bit=(lab==8)+2*(lab==6)
        adj[a].append((b,bit));adj[b].append((a,bit))
    color={}
    for start in range(n):
        if start in color:continue
        color[start]=0;queue=deque([start])
        while queue:
            a=queue.popleft()
            for b,bit in adj[a]:
                target=color[a]^bit
                if b not in color:color[b]=target;queue.append(b)
                elif color[b]!=target:return dict(passes=False,conflict_edge=[a,b],vertex_parities={str(k):v for k,v in color.items()})
    return dict(passes=True,vertex_parities={str(k):v for k,v in color.items()})

def exact_conditions(row):
    W=row['white'];corner=W+[0,1];extras=[2,3];out={}
    for conjugate in (False,True):
        G=matrix(row,conjugate);C=G.extract(corner,corner).det();assert sign(C)>0
        determinant=G.det()
        diagonals=[G.extract(corner+[i],corner+[i]).det()/C for i in extras]
        schur_det=determinant/C
        signs=[sign(x) for x in diagonals]+[sign(schur_det)]
        if not conjugate:assert signs[0]==signs[1]==-1
        passes=(all(v>=0 for v in signs) if conjugate else signs[-1]<=0)
        out['conjugate' if conjugate else 'physical']=dict(passes=passes,determinant=enc(determinant),
                corner_determinant=enc(C),schur_diagonals=list(map(enc,diagonals)),
                schur_determinant=enc(schur_det),signs=signs)
    return out

def run():
    data=json.loads((ROOT/'results'/f'{BASE}.json').read_text());rows=[];counts=Counter()
    for j,cl in enumerate(data['marked_isomorphism_classes']):
        index=cl['representative'];row=data['labelled_union_supports'][index]
        out=dict(class_index=j,union_index=index)
        if 0 in row['labels']:
            out['status']='DOTTED_PARAMETER_UNFILTERED';counts[out['status']]+=1;rows.append(out);continue
        cyc=cyclic_field(row);out['cyclic_field']=cyc
        if not cyc['passes']:
            out['status']='REJECT_CYCLE_FIELD';counts[out['status']]+=1;rows.append(out);continue
        out.update(exact_conditions(row))
        if not out['physical']['passes']:out['status']='REJECT_TWO_NEGATIVE_DIRECTIONS'
        elif not out['conjugate']['passes']:out['status']='REJECT_CONJUGATE_NOT_PSD'
        else:out['status']='PASS_PRINCIPAL_UNION_NECESSARY_TESTS'
        counts[out['status']]+=1;rows.append(out)
        print(j,out['status'],flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json']
    return dict(status='PASS_EXACT_PRINCIPAL_UNION_ARITHMETIC_FILTERS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        field=dict(primitive='sqrt3+sqrt(2+sqrt2)',minpoly=list(map(str,MINPOLY)),isolating_interval=list(map(str,INTERVAL))),
        counts=dict(counts),records=rows,
        scope='These necessary tests reject whole union-support classes independently of all outside vertices. Passing classes are NOT proved quadrilateral stars. Dotted variables and arbitrary external white attachments remain open.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'])
