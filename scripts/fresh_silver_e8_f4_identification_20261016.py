#!/usr/bin/env python3
"""Literal integral isometry of the silver E8 and F4-squared Lorentz hosts.

Standalone Fraction-pair arithmetic; no producer or project algebra imports.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_silver_e8_f4_identification_20261016'
E8_DATA = 'results/fresh_silver_e8_h9_continuation_20261016.json'
F4_DATA = 'results/fresh_silver_f4_h9_unit_h30_20261013.json'
CORNER_DATA = 'results/fresh_silver_e8_h9_continuation_20261016_h100_corners.json'
TARGET = [1,0,12,7,4,3,41,6,8]

def k(a=0,b=0): return F(a),F(b)
Z,O=k(),k(1)
def add(x,y): return x[0]+y[0],x[1]+y[1]
def neg(x): return -x[0],-x[1]
def sub(x,y): return add(x,neg(y))
def mul(x,y): return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def inv(x):
    d=x[0]*x[0]-2*x[1]*x[1]
    return x[0]/d,-x[1]/d
def div(x,y): return mul(x,inv(y))
def total(values):
    out=Z
    for x in values: out=add(out,x)
    return out
def sign(x):
    a,b=x
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    return ((a>0)-(a<0))*((a*a>2*b*b)-(a*a<2*b*b))
def dot(a,b):return total(mul(x,y)for x,y in zip(a,b,strict=True))
def tr(a):return tuple(zip(*a))
def mv(a,v):return tuple(dot(r,v)for r in a)
def mm(a,b):return tuple(tuple(dot(r,c)for c in tr(b))for r in a)
def eye(n):return tuple(tuple(O if i==j else Z for j in range(n))for i in range(n))
def integral(x):return all(v.denominator==1 for v in x)
def determinant(a):
    a=[list(r)for r in a];out=O
    for j in range(len(a)):
        t=next((i for i in range(j,len(a))if a[i][j]!=Z),None)
        if t is None:return Z
        if t!=j:a[t],a[j]=a[j],a[t];out=neg(out)
        p=a[j][j];out=mul(out,p)
        for i in range(j+1,len(a)):
            q=div(a[i][j],p)
            a[i]=[sub(x,mul(q,y))for x,y in zip(a[i],a[j])]
    return out
def inverse(a):
    n=len(a);a=[list(r)+list(e)for r,e in zip(a,eye(n))]
    for j in range(n):
        t=next(i for i in range(j,n)if a[i][j]!=Z)
        a[t],a[j]=a[j],a[t];p=a[j][j];a[j]=[div(v,p)for v in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j];a[i]=[sub(x,mul(q,y))for x,y in zip(a[i],a[j])]
    return tuple(tuple(r[n:])for r in a)
def dec(x):return F(x['a']),F(x['b'])
def matrix(x):return tuple(tuple(map(dec,r))for r in x)
def enc(x):
    if isinstance(x,tuple)and len(x)==2 and all(isinstance(v,F)for v in x):
        return dict(a=str(x[0]),b=str(x[1]))
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    if isinstance(x,dict):return {key:enc(v)for key,v in x.items()}
    return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

# Columns are F4-squared-host basis vectors expressed in the E8 host.
P=tuple(tuple(k(*x)for x in r)for r in [
 [(17,12),(0,0),(0,0),(-2,-1),(0,0),(0,0),(0,0),(-4,-3),(0,0),(0,0)],
 [(48,34),(0,0),(1,0),(-5,-4),(0,0),(0,0),(0,0),(-12,-8),(0,0),(0,0)],
 [(96,68),(1,0),(0,0),(-10,-7),(0,0),(0,0),(0,0),(-24,-16),(0,0),(0,0)],
 [(144,102),(0,0),(0,0),(-15,-10),(0,0),(0,0),(0,0),(-36,-24),(0,0),(0,0)],
 [(116,82),(0,0),(0,0),(-12,-8),(0,0),(0,0),(1,0),(-29,-20),(0,0),(0,0)],
 [(88,62),(0,0),(0,0),(-9,-6),(0,0),(1,0),(0,0),(-22,-15),(0,0),(0,0)],
 [(60,42),(0,0),(0,0),(-6,-4),(0,0),(0,0),(0,0),(-15,-10),(0,0),(0,0)],
 [(30,21),(0,0),(0,0),(-3,-2),(0,0),(0,0),(0,0),(-8,-5),(1,0),(0,0)],
 [(72,51),(0,0),(0,0),(-8,-5),(1,0),(0,0),(0,0),(-18,-12),(0,0),(0,0)],
 [(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,0)],
])

def run():
    e8_edges={(i,i+1)for i in range(6)}|{(2,7)}
    E8=tuple(tuple(k(2)if i==j else k(-1)if tuple(sorted((i,j)))in e8_edges else Z
                   for j in range(8))for i in range(8))
    F4=tuple(tuple(k(2)if i==j else k(0,-1)if {i,j}=={1,2} else k(-1)if abs(i-j)==1 else Z
                   for j in range(4))for i in range(4))
    GE=tuple(tuple(k(-1,-1)if i==j==0 else O if i==j==9 else E8[i-1][j-1]
                   if 1<=i<=8 and 1<=j<=8 else Z for j in range(10))for i in range(10))
    GF=tuple(tuple(k(-1,-1)if i==j==0 else O if i==j==9 else F4[(i-1)%4][(j-1)%4]
                   if 1<=i<=8 and 1<=j<=8 and (i-1)//4==(j-1)//4 else Z
                   for j in range(10))for i in range(10))
    e=json.loads((ROOT/E8_DATA).read_text());f=json.loads((ROOT/F4_DATA).read_text())
    c=json.loads((ROOT/CORNER_DATA).read_text())
    assert GE==matrix(e['gram']) and GF==matrix(f['gram'])
    assert determinant(GE)==determinant(GF)==k(-1,-1)
    assert mm(tr(P),mm(GE,P))==GF and determinant(P)==O
    PI=inverse(P);assert mm(P,PI)==mm(PI,P)==eye(10)
    assert all(integral(x)for r in P+PI for x in r)
    H=(O,)+(Z,)*9;image=mv(P,H)
    assert sign(dot(H,mv(GE,image)))<0
    assert dot(image,mv(GE,image))==k(-1,-1)
    assert image[0]==k(17,12) and sign(image[0])>0
    N=matrix(e['all_walls']);V=tuple(map(dec,c['finite_vertices'][51]['point']))
    assert image==tuple(mul(k(17,12),v)for v in V)
    assert [i for i,n in enumerate(N)if dot(image,mv(GE,n))==Z]==sorted(TARGET)
    assert all(sign(dot(image,mv(GE,n)))<=0 for n in N)
    for j,target in enumerate(TARGET,1):
        normal=tuple(k(-1)if i==j else Z for i in range(10))
        assert mv(P,normal)==N[target]
    paths=[f'scripts/{TAG}.py',E8_DATA,F4_DATA,CORNER_DATA,
           'results/fresh_silver_full_chamber_20261015.json',
           'results/fresh_silver_full_chamber_20261015_review.json']
    visited=set()
    def check(p):
        if p in visited:return
        visited.add(p)
        if p.endswith('.json'):
            for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():
                assert sha(q)==h,(p,q);check(q)
    for p in paths:check(p)
    return enc(dict(status='PASS_EXACT_INTEGRAL_E8_F4_SQUARED_LORENTZ_ISOMETRY',
        sources={p:sha(p)for p in paths},field='Q(sqrt2)',integer_ring='Z[sqrt2]',
        source_host='<-1-sqrt2>+F4+F4+<1>',target_host='<-1-sqrt2>+E8+<1>',
        source_gram=GF,target_gram=GE,matrix=P,inverse=PI,determinant=O,
        source_basepoint=H,target_basepoint=image,target_corner_index=51,
        source_normal_indices=list(range(9)),ordered_target_normal_indices=TARGET,
        positive_corner_scale=k(17,12),time_coordinate_mixed=True,
        recursive_manifest_files=len(visited),
        conclusion='These two literal O-lattices are integrally isometric. Their full integral reflection arrangements and groups are conjugate. The E8 host is a duplicate of the previously excluded F4-squared host.',
        scope='An explicit integral determinant-one future Lorentz isometry, not a bounded nonexistence search or an inference from rational equivalence. Frozen prior proofs remain unchanged. No statement about any further lattice is implied.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved literal-host identification replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
