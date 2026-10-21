#!/usr/bin/env python3
"""Producer-free Q(phi) identities and outward integer interval certificates."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from itertools import combinations
from math import isqrt
from pathlib import Path
from scripts import fresh_h4_blocks_independent_verify as K

ROOT=Path(__file__).resolve().parents[1]
PREFIX='data/enumerations/fresh_h9_remote_caps_20261013'
OLD='results/fresh_h9_split_second_frontier_independent_geometry_20261013_input.json'
SCALE=10**80
ZERO=(0,0)
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def ceildiv(a,b):return -((-a)//b)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def rat(x):
    x=Q(x);return x.numerator*SCALE//x.denominator,ceildiv(x.numerator*SCALE,x.denominator)
def mul(a,b):
    p=[x*y for x in a for y in b];return min(p)//SCALE,ceildiv(max(p),SCALE)
def div(a,b):
    assert b[0]>0 or b[1]<0
    return mul(a,(SCALE*SCALE//b[1],ceildiv(SCALE*SCALE,b[0])))
def sqrt(a):
    assert a[0]>0
    return isqrt(a[0]*SCALE),isqrt(a[1]*SCALE)+1
def total(xs):
    ans=ZERO
    for x in xs:ans=add(ans,x)
    return ans
def dot(a,b):return total(mul(x,y)for x,y in zip(a,b,strict=True))
def mv(a,b):return [dot(r,b)for r in a]
def pair(a,g,b):return dot(a,mv(g,b))
def enc(x):return list(map(str,x))
PHI=div(add(rat(1),sqrt(rat(5))),rat(2))
def golden(x):return add(rat(x[0]),mul(rat(x[1]),PHI))
@lru_cache(None)
def radical(field_key):
    f=json.loads(field_key);return sqrt(algebraic(f['base'],f['radicand']))
def algebraic(f,x):
    if f['kind']=='golden':return golden(K.parse(x))
    assert f['kind']=='quadratic'
    return add(algebraic(f['base'],x[0]),mul(algebraic(f['base'],x[1]),radical(json.dumps(f,sort_keys=True))))

def reconstruct(data):
    assert data['source']==OLD and data['source_sha256']==sha(OLD)
    old=read(OLD)
    for p,digest in old['source_sha256'].items():assert sha(p)==digest,(p,'source changed')
    G=K.matrix(old['gram']);I=K.matrix(old['inverse']);h=tuple(map(K.parse,old['center_covector']))
    assert K.mm(G,I)==K.eye(10)
    # Exact diagonal elimination: no numerical signature decision.
    m=[list(r)for r in G];inertia=[]
    for j in range(10):
        z=m[j][j];inertia.append(K.sign(z));assert inertia[-1]
        for a in range(j+1,10):
            for b in range(j+1,10):m[a][b]=K.sub(m[a][b],K.div(K.mul(m[a][j],m[j][b]),z))
    assert inertia.count(1)==9 and inertia.count(-1)==1
    center=tuple(K.neg(x)for x in K.mv(I,h));assert K.sign(K.pair(center,G,center))<0
    B=K.matrix(data['frame']);norms=tuple(map(K.parse,data['frame_norms']))
    assert len(B)==len(norms)==5 and all(len(row)==10 for row in B)
    assert K.sign(norms[0])<0 and all(K.sign(x)>0 for x in norms[1:])
    assert all(K.pair(B[i],G,B[j])==(norms[i]if i==j else K.Z)for i in range(5)for j in range(5))
    source={r['name']:r for r in old['roots']};S=data['orthogonal_old_walls'];assert len(set(S))==5
    seeds=[tuple(map(K.parse,source[n]['coefficients']))for n in S]
    assert all(source[n]['field']=={'kind':'golden'}for n in S)
    assert all(K.pair(a,G,b)==(K.scale(K.O,2)if i==j else K.Z)for i,a in enumerate(seeds)for j,b in enumerate(seeds))
    assert all(K.pair(a,G,b)==K.Z for a in seeds for b in B)
    assert K.sign(K.pair(B[0],G,center))<0
    assert all(K.pair(b,G,center)==K.Z for b in B[1:])
    k=Q(data['boost']);assert k>1
    ct=(k*k+1)/(2*k);st=(k*k-1)/(2*k);assert ct*ct-st*st==1
    signs=data['signs'];assert signs==[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]]
    ds=[K.div(K.scale(K.O,Q(1,2)),K.neg(norms[0])if i==0 else n)for i,n in enumerate(norms)]
    assert list(map(K.enc,ds))==data['radicands']
    # Each scaled frame vector has square +/-1/2. Hadamard signs give
    # exact norm two, six pairings zero and P^2=-1/2 without radical division.
    assert all(Q(1,2)+sum(Q(a*b,2)for a,b in zip(signs[i],signs[j],strict=True))==(2 if i==j else 0)
               for i in range(4)for j in range(4))
    V=[[mul(sqrt(golden(d)),golden(x))for x in row]for d,row in zip(ds,B,strict=True)]
    common=[add(mul(rat(st),x),mul(rat(ct),y))for x,y in zip(V[0],V[1],strict=True)]
    caps=[[add(common[i],total(mul(rat(s),V[j+2][i])for j,s in enumerate(ss)))for i in range(10)]for ss in signs]
    point=[add(mul(rat(ct),x),mul(rat(st),y))for x,y in zip(V[0],V[1],strict=True)]
    old_roots=[[algebraic(r['field'],x)for x in r['coefficients']]for r in old['roots']]
    return G,I,h,old,S,caps,point,old_roots

def check(data):
    G,I,h,old,S,caps,point,old_roots=reconstruct(data);g=[[golden(x)for x in r]for r in G];hi=list(map(golden,h))
    assert data['target_source']=='results/fresh_h9_split_multiple_edge_escape_20261013.json'
    assert sha(data['target_source'])==data['target_source_sha256']
    assert data['target_point']==read(data['target_source'])['point']
    assert all(dot(r,hi)[0]>0 for r in old_roots)
    checks=[]
    for i,r in enumerate(caps):
        assert dot(r,hi)[0]>0
        for other,n in zip(old_roots,(r['name']for r in old['roots']),strict=True):
            if n in S:continue # Zero follows from the exact frame orthogonality.
            p=pair(r,g,other);assert p[1]<-2*SCALE,('not dotted',i,n,p)
            checks.append({'new':i,'old':n,'pairing':enc(p)})
    pp=[]
    for r,n in zip(old_roots,(r['name']for r in old['roots']),strict=True):
        if n in S:continue
        p=pair(point,g,r);assert p[1]<0,('vertex infeasible',n,p);pp.append({'old':n,'pairing':enc(p)})
    # P pairs negatively with the original center; its exact norm is -1/2.
    assert dot(point,hi)[0]>0
    target=tuple(map(K.parse,data['target_point']));assert K.sign(K.pair(target,I,target))>0 and K.sign(K.pair(target,I,h))<0
    assert all(dot(r,list(map(golden,target)))[0]>0 for r in old_roots)
    cuts=[dot(r,list(map(golden,target)))for r in caps];assert all(x[1]<0 for x in cuts)
    source={r['name']:r for r in old['roots']};w=tuple(K.neg(x)for x in K.mv(I,target));q=K.pair(target,I,target)
    seed_products=[K.pair(w,G,tuple(map(K.parse,source['seed'+str(i)]['coefficients'])))for i in range(16)]
    projections={}
    for count in(7,6,5):
        positive=[];tested=0
        for block in(range(8),range(8,16)):
            for ss in combinations(block,count):
                tested+=1;qp=K.sub(q,K.total(K.scale(K.mul(seed_products[i],seed_products[i]),Q(1,2))for i in ss))
                if K.sign(qp)>0:positive.append({'walls':['seed'+str(i)for i in ss],'norm':K.enc(qp)})
        projections[str(count)]={'tested':tested,'spacelike':positive}
    assert not projections['7']['spacelike']and not projections['6']['spacelike']and len(projections['5']['spacelike'])==4
    selected=next(r for r in projections['5']['spacelike']if r['walls']==S)
    assert selected['norm']==data['projected_target_norm']
    out={'status':'PASS_EXACT_REMOTE_FOUR_CAPS','old_facets':79,'new_facets':4,'total_facets':83,
         'new_new_exact_gram':[[2 if i==j else 0 for j in range(4)]for i in range(4)],
         'actual_vertex_gram':'2I9','actual_vertex_norm':'-1/2','new_old_right_pairings':20,
         'new_old_dotted_pairings':len(checks),'dotted_checks':checks,'vertex_nonincident_checks':pp,
         'target_cut_pairings':list(map(enc,cuts)),'interval_scale':str(SCALE),
         'pure_seed_projection_screen':projections,
         'scope':'Exact frame identities and rational interval signs certify four new mutually orthogonal norm-two facets, retaining all79 old facets, and one actual simple finite vertex. A volume conclusion requires a separate certificate.'}
    if 'escape' in data:
        y=tuple(map(K.parse,data['escape']));q=K.pair(y,I,y);future=K.pair(y,I,h)
        assert K.sign(q)>0 and K.sign(future)<0
        vals=[dot(r,list(map(golden,y)))for r in old_roots+caps];assert min(x[0]for x in vals)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(future),'strict_pairings':list(map(enc,vals))}
        out['scope']+=' The saved strict future spacelike escape proves that this83-facet chamber still has infinite volume.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    data=read(PREFIX+'_caps.json');out=check(data);out['source_sha256']=sha(PREFIX+'_caps.json')
    if a.verify:assert read(PREFIX+'_verify.json')==out
    else:(ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['total_facets'],'facets', 'EXACT_ESCAPE'if'remaining_escape'in out else'volume untested')
