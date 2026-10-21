#!/usr/bin/env python3
"""Independent exact whole-H2 geometry of two bridge-face supports.

No project arithmetic imports. Sparse radical multiplication and outward
rational square-root intervals replace the existing recursive field engine.
This does not enumerate or certify any new cap.
"""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from math import isqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_bridge_faces_20261014_review'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
Z=(F(0),)*16;EZ=(Z,Z)
def c(q):return(F(q),)+(F(0),)*15
def ba(a,b):return tuple(x+y for x,y in zip(a,b))
def bs(a,q):return tuple(x*q for x in a)
@lru_cache(maxsize=100000)
def bm(a,b):
    if a==Z or b==Z:return Z
    if not any(a[1:]):return bs(b,a[0])
    if not any(b[1:]):return bs(a,b[0])
    v=[F(0)]*16
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:v[i^j]+=x*y*RAD[i&j]
    return tuple(v)
def ea(a,b):return(ba(a[0],b[0]),ba(a[1],b[1]))
def es(a,q):return(bs(a[0],q),bs(a[1],q))
def et(xs):
    z=EZ
    for x in xs:z=ea(z,x)
    return z
def inv(A):
    n=len(A);B=[list(map(F,r))+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
    for k in range(n):
        p=next(i for i in range(k,n)if B[i][k]);B[k],B[p]=B[p],B[k]
        q=B[k][k];B[k]=[x/q for x in B[k]]
        for i in range(n):
            if i!=k:
                q=B[i][k]
                if q:B[i]=[x-q*y for x,y in zip(B[i],B[k])]
    return[r[n:]for r in B]
def pd(A):
    B=[list(map(F,r))for r in A];piv=[]
    for k in range(len(B)):
        q=B[k][k];assert q>0;piv.append(q)
        for i in range(k+1,len(B)):
            for j in range(k+1,len(B)):B[i][j]-=B[i][k]*B[k][j]/q
    return piv
def inertia(A):
    B=[list(map(F,r))for r in A];pos=neg=zero=0
    while B:
        n=len(B);k=next((i for i in range(n)if B[i][i]),None)
        if k is not None:
            order=[k]+[i for i in range(n)if i!=k];B=[[B[i][j]for j in order]for i in order]
            q=B[0][0];pos+=q>0;neg+=q<0
            B=[[B[i][j]-B[i][0]*B[0][j]/q for j in range(1,n)]for i in range(1,n)]
        else:
            pair=next(((i,j)for i in range(n)for j in range(i+1,n)if B[i][j]),None)
            if pair is None:zero+=n;break
            i,j=pair;order=[i,j]+[k for k in range(n)if k not in pair]
            B=[[B[i][j]for j in order]for i in order];q=B[0][1];pos+=1;neg+=1
            B=[[B[i][j]-(B[i][0]*B[1][j]+B[i][1]*B[0][j])/q for j in range(2,n)]for i in range(2,n)]
    return[pos,neg,zero]
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(tuple,list)):return[enc(y)for y in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    d=json.loads((ROOT/BASE).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    assert d['field']['base_square_root_basis']==RAD
    D=tuple(map(F,d['field']['discriminant']));maxbits=0
    def em(a,b):
        if a==EZ or b==EZ:return EZ
        if a[1]==Z:return(bm(a[0],b[0]),bm(a[0],b[1]))
        if b[1]==Z:return(bm(a[0],b[0]),bm(a[1],b[0]))
        return(ba(bm(a[0],b[0]),bm(D,bm(a[1],b[1]))),ba(bm(a[0],b[1]),bm(a[1],b[0])))
    @lru_cache(None)
    def bounds(bits):
        den=1<<bits;out=[]
        for q in RAD:
            p=isqrt(q*den*den);out.append((F(p,den),F(p+int(p*p!=q*den*den),den)))
        return out
    def bi(a,bits):
        lo=hi=F(0)
        for x,(l,h)in zip(a,bounds(bits)):
            lo+=x*(l if x>=0 else h);hi+=x*(h if x>=0 else l)
        return lo,hi
    @lru_cache(None)
    def dbounds(bits):
        l,h=bi(D,bits);assert l>0;den=1<<bits
        p=isqrt(l.numerator*den*den//l.denominator)
        q=isqrt(h.numerator*den*den//h.denominator)+1
        assert F(p*p,den*den)<=l and F(q*q,den*den)>h
        return F(p,den),F(q,den)
    def sign(a):
        nonlocal maxbits
        if a==EZ:return 0
        bits=192
        while bits<=1536:
            l,h=bi(a[0],bits);bl,bh=bi(a[1],bits);dl,dh=dbounds(bits)
            prod=[bl*dl,bl*dh,bh*dl,bh*dh];l+=min(prod);h+=max(prod)
            if l>0 or h<0:maxbits=max(maxbits,bits);return 1 if l>0 else -1
            bits*=2
        raise AssertionError(('unresolved interval',a))
    assert sign((D,Z))>0
    parse=lambda v:tuple((tuple(map(F,a)),tuple(map(F,b)))for a,b in v)
    N=list(map(parse,d['normals']));H=parse(d['strict_center']);E=parse(d['strict_escape']['point'])
    Q=[list(map(F,r))for r in d['form']];assert len(N)==95 and inertia(Q)==[20,1,0]
    rows=[[(j,q)for j,q in enumerate(r)if q]for r in Q];cache={};dc={}
    def dot(a,b):
        key=a,b
        if key in dc:return dc[key]
        if b not in cache:cache[b]=tuple(et(es(b[j],q)for j,q in row)for row in rows)
        v=et(em(x,y)for x,y in zip(a,cache[b]));dc[key]=v;dc[b,a]=v;return v
    def lin(vs,cs):return tuple(et(em(q,v[k])for v,q in zip(vs,cs))for k in range(21))
    one=(c(1),Z)
    def rational(a):assert a[1]==Z and not any(a[0][1:]);return a[0][0]
    assert sign(dot(H,H))<0 and sign(dot(E,E))>0 and sign(dot(E,H))<0
    assert all(sign(dot(H,n))<0 and sign(dot(E,n))<0 for n in N)
    assert [rational(dot(n,n))for n in N]==list(map(F,d['norms']))
    T=tuple(et(n[k]for n in N)for k in range(21))
    Tpair=[dot(T,n)for n in N];assert sign(dot(T,T))<0 and sign(dot(T,H))<0 and all(sign(x)<0 for x in Tpair)
    records=[]
    for chosen in[61,75]:
        S=list(range(1,11))+[12,13,14,15,17,31,chosen,87]
        roots=[N[i]for i in S];G=[[rational(dot(a,b))for b in roots]for a in roots]
        assert all(G[i][i]==2 for i in range(18))
        assert all(G[i][j]in[-1,0]for i in range(18)for j in range(i))
        piv=pd(G);J=inv(G);assert min(q for r in J for q in r)>=0
        assert all(sum(G[i][k]*J[k][j]for k in range(18))==int(i==j)for i in range(18)for j in range(18))
        def project(v):
            pairs=[dot(n,v)for n in roots]
            coeff=[et(es(x,q)for x,q in zip(pairs,row))for row in J]
            out=lin([v]+roots,[one]+[es(q,-1)for q in coeff])
            assert all(dot(out,n)==EZ for n in roots)
            return out,pairs,coeff
        HP,hp,hc=project(H);EP,ep,ec=project(E);TP,tp,tc=project(T)
        assert all(sign(x)<0 for x in hc+ec)
        assert sign(dot(HP,HP))<0 and sign(dot(HP,H))<0
        assert sign(dot(EP,EP))<0 and sign(dot(EP,HP))<0
        assert sign(dot(TP,TP))<0 and sign(dot(TP,HP))<0
        norms=[];centerpairs=[];escape_pairs=[];projected=[]
        for i,n in enumerate(N):
            np,nn,nc=project(n)
            q=dot(np,np);h=dot(HP,np)
            if i in S:assert np==tuple(EZ for _ in range(21))and q==h==EZ
            else:
                assert sign(q)<0 and sign(h)<0
                assert all(sign(x)<=0 for x in nn)
            norms.append(q);centerpairs.append(h);escape_pairs.append(dot(EP,np));projected.append(np)
        assert TP==tuple(et(v[k]for v in projected)for k in range(21))
        triple=[0,11,28];A=[[dot(projected[i],projected[j])for j in triple]for i in triple]
        det=EZ
        for p in permutations(range(3)):
            v=one
            for i,j in enumerate(p):v=em(v,A[i][j])
            if sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))%2:v=es(v,-1)
            det=ea(det,v)
        assert sign(det)<0
        # Full residual span makes an acute nonzero vector strictly negative
        # on TP=sum projected walls. Thus a separator cutting future EP is
        # neither future nor past causal and must be spacelike.
        rec=dict(support=S,rational_gram=G,inverse=J,PD_pivots=piv,
                 strict_face_point=HP,strict_face_point_norm=dot(HP,HP),
                 projected_escape=EP,projected_escape_norm=dot(EP,EP),escape_future_pairing=dot(EP,HP),
                 projected_positive_wall_sum=TP,wall_sum_norm=dot(TP,TP),wall_sum_future_pairing=dot(TP,HP),
                 all_projected_wall_norms=norms,all_face_center_pairings=centerpairs,
                 all_projected_escape_pairings=escape_pairs,
                 residual_spanning_walls=triple,residual_spanning_gram=A,residual_gram_determinant=det,
                 nonzero_future_timelike_projections=77,
                 conclusion='Every nonsupport wall projection is strictly future timelike. Its inequality is strict at every future hyperbolic point of the entire residual H2 plane. No old wall induces a boundary geodesic on this face.')
        records.append(rec);print('PASS_FACE',chosen,'77 timelike projections',flush=True)
    paths=[BASE,'scripts/fresh_h20_e7_edge_cusp_repair_20261014_exact.py',f'scripts/{TAG}.py']
    return enc(dict(status='PASS_INDEPENDENT_TWO_COMPLETE_H2_BRIDGE_FACES',sources={p:sha(p)for p in paths},
                    field=d['field'],arithmetic='Separate sparse radical multiplication and outward rational intervals; no project arithmetic imports',
                    maximum_sign_enclosure_bits=maxbits,ambient_inertia=[20,1,0],base_norms_checked=95,
                    positive_wall_sum_norm=dot(T,T),positive_wall_sum_all_pairings=Tpair,faces=records,
                    scope='Exact geometry of these two specified PD18 supports and all95 old-wall projections. No numerical cone-LP infeasibility or new cap/angle certificate is claimed. Other PD18 supports are not enumerated. Old95 mutual Coxeter gates are inherited from the hash-bound base; all support, residual, orientation and sum-center calculations are rebuilt.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
