#!/usr/bin/env python3
"""Exact Q(phi) principal-union filters using squareclass vertex rescaling.

No degree-sixteen field arithmetic or numerical eigenvalue tests are needed.
The input union enumeration is a separate completeness obligation.
"""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt5_quad_exhaustion_20261015'
TAG='fresh_sqrt5_quad_filters_20261015'
def sg(x):return (x>0)-(x<0)
class K:
    def __init__(self,a=0,b=0):
        if isinstance(a,K):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,o):o=K(o);return K(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,o):return self+-K(o)
    def __rsub__(self,o):return K(o)+-self
    def __mul__(self,o):
        o=K(o);return K(self.a*o.a+self.b*o.b,self.a*o.b+self.b*o.a+self.b*o.b)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=K(o);norm=o.a*o.a+o.a*o.b-o.b*o.b;assert norm
        return self*K((o.a+o.b)/norm,-o.b/norm)
    def __rtruediv__(self,o):return K(o)/self
    def __eq__(self,o):o=K(o);return self.a==o.a and self.b==o.b
    def __bool__(self):return bool(self.a or self.b)
    def conj(self):return K(self.a+self.b,-self.b)
    def sign(self):
        a,b=2*self.a+self.b,self.b
        if not b:return sg(a)
        if not a or sg(a)==sg(b):return sg(b)
        return sg(a)*sg(a*a-5*b*b)
    def enc(self):return [str(self.a),str(self.b)]
    def integral(self):return self.a.denominator==self.b.denominator==1

PHI=K(0,1)
BITS={3:0,4:1,5:0,6:2,10:4}
RADICANDS=[K(2),K(3),K(2,1)]
SQUARES={2:K(),3:K(1),4:K(2),5:K(1,1),6:K(3),10:K(2,1)}
def dot(a,b):return sum((x*y for x,y in zip(a,b)),K())
def sub(M,ids):return [[M[i][j] for j in ids] for i in ids]
def encode(M):return [[v.enc() for v in row] for row in M]

def coloring(row):
    n=row['nodes'];adj=[[] for _ in range(n)]
    for (i,j),label in zip(combinations(range(n),2),row['labels']):
        assert label!=0
        if label==2:continue
        bit=BITS[label];adj[i].append((j,bit));adj[j].append((i,bit))
    color=[None]*n
    for start in range(n):
        if color[start] is not None:continue
        color[start]=0;queue=deque([start])
        while queue:
            i=queue.popleft()
            for j,bit in adj[i]:
                value=color[i]^bit
                if color[j] is None:color[j]=value;queue.append(j)
                elif color[j]!=value:return dict(passes=False,conflict=[i,j],colors=color)
    return dict(passes=True,colors=color)

def scaled_matrix(row,colors):
    n=row['nodes'];M=[[K() for _ in range(n)] for _ in range(n)]
    for i,mask in enumerate(colors):
        value=K(2)
        for bit,d in enumerate(RADICANDS):
            if mask&(1<<bit):value*=d
        M[i][i]=value;assert value.sign()>0 and value.conj().sign()>0
    for (i,j),label in zip(combinations(range(n),2),row['labels']):
        if label==2:continue
        masks=(colors[i],colors[j],BITS[label]);value=PHI if label==5 else K(1)
        for bit,d in enumerate(RADICANDS):
            exponent=sum(bool(mask&(1<<bit)) for mask in masks);assert exponent in (0,2)
            if exponent:value*=d
        M[i][j]=M[j][i]=-value
        assert 4*value*value/(M[i][i]*M[j][j])==SQUARES[label]
    assert all(x.integral() for row in M for x in row)
    return M

def schur(M,corner):
    ids=corner+[i for i in range(len(M)) if i not in corner];A=sub(M,ids)
    pivots=[]
    for k in range(len(corner)):
        d=A[k][k];assert d.sign()>0
        pivots.append(d)
        for i in range(k+1,len(A)):
            for j in range(k+1,len(A)):A[i][j]-=A[i][k]*A[k][j]/d
    residual=[r[len(corner):] for r in A[len(corner):]]
    return residual,pivots

def inertia(M):
    A=[r[:] for r in M];pos=neg=zero=0
    while A:
        n=len(A);p=next((i for i in range(n) if A[i][i]),None)
        if p is not None:
            A=sub(A,[p]+[i for i in range(n) if i!=p]);d=A[0][0]
            pos+=d.sign()>0;neg+=d.sign()<0
            A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,n)] for i in range(1,n)]
        else:
            p=next(((i,j) for i,j in combinations(range(n),2) if A[i][j]),None)
            if p is None:zero+=n;break
            A=sub(A,list(p)+[i for i in range(n) if i not in p]);b=A[0][1];pos+=1;neg+=1
            A=[[A[i][j]-(A[i][0]*A[1][j]+A[i][1]*A[0][j])/b for j in range(2,n)] for i in range(2,n)]
    return [pos,neg,zero]

def nullspace(M):
    A=[r[:] for r in M];n=len(A);pivots=[];rank=0
    for col in range(n):
        p=next((i for i in range(rank,n) if A[i][col]),None)
        if p is None:continue
        A[rank],A[p]=A[p],A[rank];d=A[rank][col];A[rank]=[x/d for x in A[rank]]
        for i in range(n):
            if i!=rank:
                d=A[i][col];A[i]=[a-d*b for a,b in zip(A[i],A[rank])]
        pivots.append(col);rank+=1
    result=[]
    for free in range(n):
        if free in pivots:continue
        v=[K(i==free) for i in range(n)]
        for i,p in enumerate(pivots):v[p]=-A[i][free]
        assert all(not dot(row,v) for row in M);result.append(v)
    return result

def conditions(row,M=None):
    if M is None:
        colors=coloring(row);assert colors['passes'];M=scaled_matrix(row,colors['colors'])
    corner=row['white']+[0,1];out={}
    for conjugate in (False,True):
        G=[[x.conj() for x in r] for r in M] if conjugate else M
        S,pivots=schur(G,corner);assert len(S)==2
        determinant=S[0][0]*S[1][1]-S[0][1]*S[1][0]
        signs=[S[0][0].sign(),S[1][1].sign(),determinant.sign()]
        if not conjugate:assert signs[:2]==[-1,-1]
        passes=all(s>=0 for s in signs) if conjugate else signs[-1]<=0
        out['conjugate' if conjugate else 'physical']=dict(passes=passes,schur=encode(S),
            determinant=determinant.enc(),signs=signs,corner_pivots=[x.enc() for x in pivots])
    return out

def catalogue_checks():
    d=json.loads((ROOT/'results/fresh_sqrt5_lanner_catalogue_20261015.json').read_text());checked=0
    for block in d['higher_orders']:
        for rec in block['accepted_classes']:
            n=block['order'];row=dict(nodes=n,labels=rec['labels']);c=coloring(row);assert c['passes']
            M=scaled_matrix(row,c['colors']);assert inertia(M)==[n-1,1,0]
            assert inertia([[x.conj() for x in r] for r in M])==[n,0,0]
            for omitted in range(n):assert inertia(sub(M,[i for i in range(n) if i!=omitted]))==[n-1,0,0]
            checked+=1
    # Cross-radical tree tests and every cycle field constraint are exact.
    for a in (3,4,5,6,10):
        for b in (3,4,5,6,10):
            row=dict(nodes=3,labels=[a,2,b]);c=coloring(row);assert c['passes'];scaled_matrix(row,c['colors']);checked+=1
    return checked

def run():
    checks=catalogue_checks();base=json.loads((ROOT/f'results/{BASE}.json').read_text());records=[];counts=Counter()
    for j,cl in enumerate(base['marked_isomorphism_classes']):
        index=cl['representative'];row=base['labelled_union_supports'][index];out=dict(class_index=j,union_index=index)
        if 0 in row['labels']:out['status']='DOTTED_PARAMETER_UNFILTERED'
        else:
            c=coloring(row);out['cyclic_field']=c
            if not c['passes']:out['status']='REJECT_CYCLE_FIELD'
            else:
                M=scaled_matrix(row,c['colors']);out.update(conditions(row,M))
                if not out['physical']['passes']:out['status']='REJECT_TWO_NEGATIVE_DIRECTIONS'
                elif not out['conjugate']['passes']:out['status']='REJECT_CONJUGATE_NOT_PSD'
                else:
                    out['status']='PASS_PRINCIPAL_UNION_NECESSARY_TESTS'
                    out['singular']=out['physical']['signs'][-1]==0
                    out['scaled_gram']=encode(M)
        counts[out['status']]+=1;records.append(out)
        if (j+1)%100==0:print('classes',j+1,'of',len(base['marked_isomorphism_classes']),dict(counts),flush=True)
    files=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',
        'scripts/fresh_sqrt5_lanner_catalogue_20261015.py','results/fresh_sqrt5_lanner_catalogue_20261015.json',
        'results/fresh_sqrt5_lanner_catalogue_20261015_review.json']
    return dict(status='PASS_QSQRT5_PRINCIPAL_UNION_ARITHMETIC_FILTERS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        scaling='Vertex color masks in independent squareclasses2,3,phi+2; positive radical rescaling puts the Gram in Qphi.',
        scaling_sanity_checks=checks,counts=dict(counts),records=records,
        scope='Necessary principal-union obstructions independent of outside roots. Passing cores, dotted distances and arbitrary white attachments remain to be classified. No dimension bound.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--catalogue-check',action='store_true');args=ap.parse_args()
    if args.catalogue_check:print('PASS_QPHI_RADICAL_RESCALING',catalogue_checks())
    else:
        data=run();p=ROOT/f'results/{TAG}.json'
        if args.verify:assert data==json.loads(p.read_text())
        else:p.write_text(json.dumps(data,indent=2)+'\n')
        print(data['status'],data['counts'])
