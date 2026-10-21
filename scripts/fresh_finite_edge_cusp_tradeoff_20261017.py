#!/usr/bin/env python3
"""Exact rank-ten star checks and finite-edge triangle coefficient14.

Standalone Fraction arithmetic. Classification completeness is the explicitly
bound primary Chein table, not a claim that these three graphs were enumerated
here. The companion report proves the outside-white and face-count transfers.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite_edge_cusp_tradeoff_20261017'

def sub(A,ids):return [[A[i][j] for j in ids] for i in ids]

def inverse(A):
    n=len(A);M=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        p=next(j for j in range(i,n) if M[j][i]);M[i],M[p]=M[p],M[i]
        d=M[i][i];M[i]=[x/d for x in M[i]]
        for j in range(n):
            if j!=i:
                c=M[j][i];M[j]=[x-c*y for x,y in zip(M[j],M[i])]
    return [r[n:] for r in M]

def determinant(A):
    M=[[F(x) for x in r] for r in A];d=F(1);n=len(M)
    for i in range(n):
        p=next((j for j in range(i,n) if M[j][i]),None)
        if p is None:return F(0)
        if p!=i:M[i],M[p]=M[p],M[i];d=-d
        x=M[i][i];d*=x
        for j in range(i+1,n):
            c=M[j][i]/x
            for k in range(i+1,n):M[j][k]-=c*M[i][k]
    return d

def inertia(A):
    A=[[F(x) for x in r] for r in A];p=m=z=0
    while A:
        n=len(A);i=next((i for i in range(n) if A[i][i]),None)
        if i is None:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if A[i][j]),None)
            if pair is None:z+=n;break
            i,j=pair;ids=[k for k in range(n) if k not in pair];c=A[i][j]
            A=[[A[a][b]-(A[a][i]*A[j][b]+A[a][j]*A[i][b])/c for b in ids] for a in ids]
            p+=1;m+=1
        else:
            d=A[i][i];p+=d>0;m+=d<0;ids=[j for j in range(n) if j!=i]
            A=[[A[a][b]-A[a][i]*A[i][b]/d for b in ids] for a in ids]
    return [p,m,z]

def kernel1(A):
    M=[[F(x) for x in r] for r in A];n=len(M);piv=[];row=0
    for col in range(n):
        j=next((j for j in range(row,n) if M[j][col]),None)
        if j is None:continue
        M[row],M[j]=M[j],M[row];d=M[row][col];M[row]=[x/d for x in M[row]]
        for j in range(n):
            if j!=row:
                c=M[j][col];M[j]=[x-c*y for x,y in zip(M[j],M[row])]
        piv.append(col);row+=1
    free=[i for i in range(n) if i not in piv];assert len(free)==1
    v=[F(0)]*n;v[free[0]]=F(1)
    for i,j in enumerate(piv):v[j]=-M[i][free[0]]
    if all(x<0 for x in v):v=[-x for x in v]
    assert all(x>0 for x in v)
    assert all(sum(x*y for x,y in zip(r,v))==0 for r in A)
    return v

def graph(name):
    path=8 if name=='A' else 9
    edges=[[i,i+1,3] for i in range(path-1)]
    edges+= [[2,8,3],[6,9,3]] if name=='A' else [[6,9,3]]
    if name=='C':edges[0][2]=4
    # Positive diagonal congruence: multiply normalized norm-two root0 by
    # sqrt2 only in C. This makes the entire tree Gram rational.
    G=[[F(2*(i==j)) for j in range(10)] for i in range(10)]
    if name=='C':G[0][0]=F(4)
    for i,j,m in edges:G[i][j]=G[j][i]=F(-2 if m==4 else -1)
    for i,j,m in edges:assert 4*G[i][j]**2/(G[i][i]*G[j][j])=={3:1,4:2}[m]
    assert inertia(G)==[9,1,0]
    deletions=[];aff=[]
    for j in range(10):
        ids=[i for i in range(10) if i!=j];H=sub(G,ids);I=inertia(H)
        assert I in ([9,0,0],[8,0,1])
        rec=dict(deleted=j,inertia=I,determinant=str(determinant(H)))
        if I[-1]:
            aff.append(j);rec['positive_kernel']=[str(x) for x in kernel1(H)]
        deletions.append(rec)
    assert aff=={'A':[0,7,9],'B':[0],'C':[0,8]}[name]
    out=dict(name=name,edges=edges,rational_gram=[[str(x) for x in r] for r in G],
             determinant=str(determinant(G)),inertia=inertia(G),
             positive_diagonal_gauge_sqrt2_at_node0=name=='C',
             deletions=deletions,affine_deletions=aff)
    if name=='B':
        thresholds=[]
        for j in range(1,10):
            ids=[i for i in range(10) if i!=j];H=sub(G,ids)
            t=1/inverse(H)[ids.index(0)][ids.index(0)]
            K=[r[:] for r in H];K[0][0]-=1
            assert inertia(K)==([8,0,1] if j==8 else [9,0,0])
            thresholds.append(dict(deleted=j,threshold=str(t),inertia_at_tau1=inertia(K)))
        assert [x['threshold'] for x in thresholds]==['2','3/2','4/3','5/4','6/5','7/6','9/8','1','10/9']
        out['outside_white_diagonal_thresholds']=thresholds
        out['allowed_face_counts']=[[9,1],[8,2]]
    else:out['allowed_face_counts']=[[10-len(aff),len(aff)]]
    return out

def run():
    diagrams=[graph(x) for x in 'ABC'];rows=[]
    for vertices in range(2,10):
        for f in range(2,vertices+1):
            i=vertices-f;c=F(2)*(i+F(f-2,3));assert c<=14
            rows.append(dict(f=f,i=i,dimension=vertices-1,coefficient=str(c)))
    for f,i in [(7,3),(8,2),(9,1)]:
        c=F(2)*(i+F(f-2,3));assert c<=F(28,3)<14
        rows.append(dict(f=f,i=i,dimension=9,coefficient=str(c)))
    assert [r for r in rows if r['coefficient']=='14']==[dict(f=2,i=7,dimension=8,coefficient='14')]
    sources=[f'scripts/{TAG}.py',f'data/references/{TAG}_chein.pdf',
             'results/fresh_finite_ended_simplex_ports_20261017.md',
             'scripts/fresh_multi_ideal_triangle_scope_20261017.py',
             'results/fresh_multi_ideal_triangle_scope_20261017.json',
             'results/fresh_multi_ideal_triangle_scope_20261017.md',
             'results/fresh_finite_volume_simplex_truncation_20261017.md']
    return dict(status='PASS_EXACT_RANK10_STARS_AND_FINITE_EDGE_TRIANGLE_COEFFICIENT14',
      scope='Uniform finite-edge triangle payment; no ambient dimension bound; Chein completeness is a primary classification input.',
      sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
      classification=dict(url='https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf',
                          location='Appendix, printed p16 / PDF page15, Rang10',
                          independent_full_reclassification=False),
      diagrams=diagrams,allowed_nine_simplex_face_counts=[[7,3],[8,2],[9,1]],
      block_coefficients=rows,uniform_coefficient='14',
      inequality='2(T0+T1) <= 14 E_ff = 7 n v_f - 7 E_fi <= 7 n v_f',
      unresolved_equality_type=dict(f=2,i=7,dimension=8),
      literature=[
       dict(class_scope='simple all-ideal Coxeter',dimension_upper_bound=8,
            url='https://arxiv.org/pdf/math/0502413',location='Theorem1'),
       dict(class_scope='all-ideal right-angled Coxeter',dimension_upper_bound=6,
            url='https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf',location='Corollary7, p1959'),
       dict(class_scope='right-angled finite-volume, preprint claim',finite_vertex_lower_bounds={'7':4,'8':41,'9':1312},
            url='https://arxiv.org/pdf/2603.29925',location='Theorem1.2 and Table2')])

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'],'saved_replay=',args.verify)
