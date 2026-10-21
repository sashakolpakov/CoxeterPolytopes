#!/usr/bin/env python3
"""Independent complete height100 replay, with reversed positive coordinates.

No producer, Golden, project Schur, or symbolic-algebra imports.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import cmp_to_key
from hashlib import sha256
from math import lcm
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta11_vinberg_20261014_review'
PRIMARY='results/fresh_h9_beta11_vinberg_20261014.json'
BASIS='results/fresh_h9_h8_cap_height_20261014_host_basis.json'
ZEROFILE='results/fresh_h9_beta11_zero_sector_20261014.json'
POOL='results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json'

def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[0],x[1]-y[1]
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def scale(x,k):return x[0]*k,x[1]*k
def sign(x,conjugate=False):
    a,b=x
    if conjugate:a,b=a+b,-b
    u=2*a+b
    if not b:return (u>0)-(u<0)
    if not u:return (b>0)-(b<0)
    if (u>0)==(b>0):return (u>0)-(u<0)
    z=u*u-5*b*b;assert z
    return ((u>0)-(u<0))*((z>0)-(z<0))

class K:
    def __init__(self,a=0,b=0):
        if isinstance(a,K):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,y):
        y=K(y);return K(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,y):return self+-K(y)
    def __rsub__(self,y):return K(y)+-self
    def __mul__(self,y):
        y=K(y);return K(*mul((self.a,self.b),(y.a,y.b)))
    __rmul__=__mul__
    def conjugate(self):return K(self.a+self.b,-self.b)
    def __truediv__(self,y):
        y=K(y);n=y.a*y.a+y.a*y.b-y.b*y.b;assert n
        c=y.conjugate();return self*K(c.a/n,c.b/n)
    def __rtruediv__(self,y):return K(y)/self
    def __eq__(self,y):
        y=K(y);return self.a==y.a and self.b==y.b
    def __bool__(self):return bool(self.a or self.b)
    def __hash__(self):return hash((self.a,self.b))
    def __lt__(self,y):return (self-y).sign()<0
    def sign(self,conjugate=False):return sign((self.a,self.b),conjugate)
    def integral(self):return self.a.denominator==self.b.denominator==1
    def ri(self):assert self.integral();return int(self.a),int(self.b)
    def enc(self):return dict(a=str(self.a),b=str(self.b))

PHI=K(0,1);BETA=3+PHI;DELTA=1+4*PHI
def dec(x):return K(x['a'],x['b'])
def old(x):
    assert all(F(t)==0 for t in x[1])
    a,b=map(F,x[0]);return K(a-b,2*b)
def enc(x):
    if isinstance(x,K):return x.enc()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def inverse(M):
    n=len(M);a=[[F(x)for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(M)]
    for i in range(n):
        j=next(j for j in range(i,n)if a[j][i]);a[i],a[j]=a[j],a[i]
        q=a[i][i];a[i]=[x/q for x in a[i]]
        for j in range(n):
            if j!=i:
                q=a[j][i];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
    return [row[n:]for row in a]
def dot(x,y):return sum((a*b for a,b in zip(x,y)),K())
def mv(M,v):return [dot(row,v)for row in M]
def pair(x,y,G):return dot(x,mv(G,y))
def key(r):return tuple(x.ri()for x in r)
def beta_divisible(x):return (x[0]+8*x[1])%11==0

def run():
    d,basis,zero,pool=map(read,(PRIMARY,BASIS,ZEROFILE,POOL))
    assert F(d['height_bound'])==100
    G=[[old(x)for x in row]for row in basis['gram']]
    assert G==[[dec(x)for x in row]for row in d['form']]
    C=[[G[i+1][j+1].a for j in range(8)]for i in range(8)]
    assert all(not G[i+1][j+1].b for i in range(8)for j in range(8))
    J=inverse(C);assert all(x.denominator==1 and x>=0 for row in J for x in row)
    order=list(range(7,-1,-1))
    # F^{-1} has diagonal1. Prefix minimum energy is the inverse of its
    # corresponding principal block; no eliminated-suffix routine is reused.
    Finv=[[C[i][j]/2 for j in order]+[F(0)]for i in order]+[[F(0)]*8+[F(1)]]
    positive=inverse(Finv)
    assert all(Finv[i][i]==1 for i in range(9))
    assert positive==[[2*J[i][j]for j in order]+[F(0)]for i in order]+[[F(0)]*8+[F(1)]]
    forms=[];minimum_forms=[]
    for n in range(1,10):
        M=inverse([row[:n]for row in Finv[:n]])
        den=lcm(*(x.denominator for row in M for x in row))
        forms.append((den,[[int(x*den)for x in row]for row in M]));minimum_forms.append(M)
    assert all(M[i][i]>0 for i,M in enumerate(minimum_forms))
    # Independent universal coefficient boxes, fixed for T=100.
    assert (K(2)-PHI).sign()>0 and (PHI-1).sign()>0
    assert (K(10)-2*BETA).sign()>0 and (K(16)-4*BETA.conjugate()).sign()>0
    assert (K(400)-400/PHI).sign()>0 and (K(9)-4*PHI).sign()>0
    assert 2020<45*45
    controls=[];control_sets=[]
    for q in (K(2),2*BETA):
        xs=[]
        for a in range(-15,16):
            for b in range(-12,13):
                x=K(a,b);target=2*q+DELTA*x*x
                if x.sign()>0 and (K(100)-DELTA*x*x/(2*q)).sign()>=0 and target.sign(True)>=0:
                    xs.append(x);controls.append((DELTA*x*x/(2*q),q,x))
        xs.sort();control_sets.append((q,xs))
    controls.sort(key=lambda z:z[0])
    assert len(controls)==len(d['levels'])==46
    assert all(a[0]!=b[0]for a,b in zip(controls,controls[1:]))
    for (q,xs),saved in zip(control_sets,d['control_boxes']):
        assert q==dec(saved['norm']) and xs==list(map(dec,saved['values']))
    cone=[list(map(old,r))for r in zero['simple_zero_walls']]
    expected_cone=[]
    for j in range(1,9):expected_cone.append([K(-int(i==j))for i in range(10)])
    expected_cone.append([K(int(i==9))for i in range(10)])
    assert cone==expected_cone
    H=[K(2)]+[K()]*8+[K(1)]
    assert mv(G,H)==[-DELTA]+[K()]*9 and pair(H,H,G)==-2*DELTA
    previous=[];levels=[];total=Counter()
    for index,((height,q,x),saved) in enumerate(zip(controls,d['levels'])):
        assert height==dec(saved['height']) and q==dec(saved['norm']) and x==dec(saved['control'])
        target=2*q+DELTA*x*x;ti=target.ri();xi=x.ri()
        alphabet=[]
        for a in range(-29,30):
            for b in range(-25,26):
                z=(a,b);remaining=sub(ti,mul(z,z))
                if sign(z)>=0 and sign(remaining)>=0 and sign(remaining,True)>=0:alphabet.append(z)
        alphabet.sort(key=cmp_to_key(lambda a,b:sign(sub(a,b))))
        assert set(alphabet)=={dec(v).ri()for v in saved['alphabet']['values']}
        coefficients=[];bounds=[]
        for r in previous:
            L=r[0]-2*r[9]
            coeff=[2*r[i+1]for i in order]+[L]
            assert all(v.sign()>=0 for v in coeff)
            coefficients.append([v.ri()for v in coeff]);bounds.append((DELTA*x*r[0]).ri())
        visited=[0]*10;sp=[0]*10;wp=[0]*10;counts=Counter();roots=[]
        def recurse(prefix,minimum,partials):
            depth=len(prefix);visited[depth]+=1
            if depth:
                num,den=minimum;rem=sub(scale(ti,den),num)
                if sign(rem)<0 or sign(rem,True)<0:sp[depth]+=1;return
                z=prefix[-1]
                partials=[add(p,mul(z,co[depth-1]))for p,co in zip(partials,coefficients)]
                if any(sign(sub(b,p))<0 for b,p in zip(bounds,partials)):wp[depth]+=1;return
            if depth==9:
                num,den=minimum
                if num!=scale(ti,den):counts['norm_mismatch']+=1;return
                weights=[K(*z)for z in reversed(prefix[:8])]
                l=K(*prefix[8]);y=mv(J,weights);t=(x-l)/2
                r=[x,*y,t];assert all(z.integral()for z in r)
                assert pair(r,r,G)==q
                assert all((2*z/q).integral()for z in mv(G,r))
                roots.append(r);return
            den,M=forms[depth];base=(0,0);linear=(0,0)
            for i in range(depth):
                linear=add(linear,scale(prefix[i],2*M[i][depth]))
                for j in range(depth):base=add(base,scale(mul(prefix[i],prefix[j]),M[i][j]))
            for z in alphabet:
                if q==2*BETA and not beta_divisible(z):counts['residue_pruned']+=1;continue
                if depth==8 and any((z[i]-xi[i])%2 for i in range(2)):counts['parity_pruned']+=1;continue
                num=add(add(base,mul(linear,z)),scale(mul(z,z),M[depth][depth]))
                recurse(prefix+[z],(num,den),partials)
        recurse([],None,[(0,0)]*len(previous))
        expected=[list(map(dec,r))for r in saved['roots']]
        assert {key(r)for r in roots}=={key(r)for r in expected}
        assert len(roots)==len(expected)==len({key(r)for r in roots})
        assert all(pair(a,b,G).sign()<=0 for i,a in enumerate(roots)for b in roots[i+1:])
        previous.extend(roots);total.update(counts);total['nodes']+=sum(visited);total['levels']+=1
        levels.append(dict(control=x,norm=q,height=height,alphabet_size=len(alphabet),
            visited_by_depth=visited,schur_pruned_by_depth=sp,wall_pruned_by_depth=wp,
            counts=dict(counts),roots=roots))
        print('level',index+1,'roots',len(roots),'nodes',sum(visited),flush=True)
    allroots=cone+previous;assert len(allroots)==21
    archived=[list(map(dec,r))for r in d['roots']]
    assert {key(r)for r in allroots}=={key(r)for r in archived}
    assert len({key(r)for r in allroots})==21
    for i,a in enumerate(archived):
        for j,b in enumerate(archived):
            assert pair(a,b,G)==dec(d['gram'][i][j])
            assert (2*pair(a,b,G)/pair(b,b,G)).integral()
            if i!=j:assert pair(a,b,G).sign()<=0
    center=list(map(dec,d['strict_center']))
    assert pair(center,center,G).sign()<0 and pair(center,H,G).sign()<0
    assert all(pair(center,n,G).sign()<0 for n in archived)
    for row in d['actual_facet_points']:
        i=row['index'];P=list(map(dec,row['point']));assert pair(P,P,G)==dec(row['norm'])
        assert pair(P,P,G).sign()<0 and pair(P,H,G).sign()<0
        assert all(pair(P,n,G).sign()==(0 if j==i else -1)for j,n in enumerate(archived))
    paths=[PRIMARY,BASIS,ZEROFILE,POOL,f'scripts/{TAG}.py','scripts/fresh_h9_beta11_vinberg_20261014.py']
    return enc(dict(status='PASS_INDEPENDENT_COMPLETE_BOTH_NORM_HEIGHT100_REPLAY',
        sources={p:sha(p)for p in paths},height_bound=100,
        arithmetic='Standalone Fraction field and integer-pair hot loops; no Golden, producer, old Schur or symbolic algebra imports.',
        variable_order=['lambda8','lambda7','lambda6','lambda5','lambda4','lambda3','lambda2','lambda1','l'],
        universal_boxes=dict(control_a=15,control_b=12,coordinate_a=29,coordinate_b=25),
        positive_form=positive,prefix_minimum_forms=minimum_forms,control_sets=[dict(norm=q,values=x)for q,x in control_sets],
        counts=dict(total),levels=levels,wall_count=21,actual_facet_count=21,
        scope='Complete independent positive-height enumeration through100, both norm classes and all46levels. Zero-sector and complete norm-alphabet proofs remain separately certified dependencies. The exact root sets, full441Gramentries, strictcenter and21facetpoints are independently checked. No higher-height or globalvolume assertion.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
