#!/usr/bin/env python3
"""Independent sparse-radical replay of the symmetric bridge-edge closure."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
import fresh_h20_bridge_faces_20261014_review as A
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_bridge_edge_closure_20261014_review'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
PREV='results/fresh_h20_e7_bridge_face_20261014_exact.json'
INPUT='results/fresh_h20_bridge_edge_closure_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def parse(x,level):
    return tuple(map(F,x))if level==0 else tuple(parse(y,level-1)for y in x)
class Field:
    def __init__(self,radicals):
        self.radicals=radicals;self.level=len(radicals);self.zeros=[A.Z];self.bits=0
        for _ in radicals:self.zeros.append((self.zeros[-1],self.zeros[-1]))
    def const(self,q,level=None):
        if level is None:level=self.level
        v=A.c(q)
        for k in range(level):v=(v,self.zeros[k])
        return v
    def lift(self,x,level):
        for k in range(level,self.level):x=(x,self.zeros[k])
        return x
    def add(self,x,y,level=None):
        if level is None:level=self.level
        if not level:return A.ba(x,y)
        return tuple(self.add(a,b,level-1)for a,b in zip(x,y))
    def scale(self,x,q,level=None):
        if level is None:level=self.level
        if not level:return A.bs(x,q)
        return tuple(self.scale(a,q,level-1)for a in x)
    def mul(self,x,y,level=None):
        if level is None:level=self.level
        if x==self.zeros[level]or y==self.zeros[level]:return self.zeros[level]
        if not level:return A.bm(x,y)
        a,b=x;c,d=y;z=self.zeros[level-1]
        if b==z:return(self.mul(a,c,level-1),self.mul(a,d,level-1))
        if d==z:return(self.mul(a,c,level-1),self.mul(b,c,level-1))
        ac=self.mul(a,c,level-1);bd=self.mul(b,d,level-1)
        return(self.add(ac,self.mul(self.radicals[level-1],bd,level-1),level-1),
               self.add(self.mul(a,d,level-1),self.mul(b,c,level-1),level-1))
    def total(self,xs,level=None):
        if level is None:level=self.level
        z=self.zeros[level]
        for x in xs:z=self.add(z,x,level)
        return z
    @lru_cache(None)
    def basebounds(self,bits):
        den=1<<bits;out=[]
        for q in A.RAD:
            n=isqrt(q*den*den);out.append((F(n,den),F(n+int(n*n!=q*den*den),den)))
        return out
    def interval(self,x,level,bits):
        if not level:
            lo=hi=F(0)
            for q,(l,h)in zip(x,self.basebounds(bits)):
                lo+=q*(l if q>=0 else h);hi+=q*(h if q>=0 else l)
            return lo,hi
        l,h=self.interval(x[0],level-1,bits);a,b=self.interval(x[1],level-1,bits)
        c,d=self.sqrtbounds(level,bits);v=[a*c,a*d,b*c,b*d]
        return l+min(v),h+max(v)
    @lru_cache(None)
    def sqrtbounds(self,level,bits):
        l,h=self.interval(self.radicals[level-1],level-1,bits)
        assert l>0
        den=1<<bits;p=isqrt(l.numerator*den*den//l.denominator)
        q=isqrt(h.numerator*den*den//h.denominator)+1
        assert F(p*p,den*den)<=l and F(q*q,den*den)>h
        return F(p,den),F(q,den)
    def sign(self,x,level=None):
        if level is None:level=self.level
        if x==self.zeros[level]:return 0
        for bits in[192,384,768,1536,3072]:
            try:l,h=self.interval(x,level,bits)
            except AssertionError:continue
            if l>0 or h<0:self.bits=max(self.bits,bits);return 1 if l>0 else -1
        raise AssertionError(('unresolved interval',level,x))

def run():
    old,pre,d,face=[json.loads((ROOT/p).read_text())for p in[
        BASE,PREV,INPUT,'results/fresh_h20_bridge_faces_20261014_review.json']]
    for doc in[old,pre,d,face]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert old['field']['base_square_root_basis']==A.RAD
    assert d['field']['base']==pre['field']
    ff=Field([parse(old['field']['discriminant'],0),
              parse(pre['field']['new_positive_sqrt_radical'],1),
              parse(d['field']['new_positive_sqrt_radical'],2)])
    add,mul,scale,total,sign=ff.add,ff.mul,ff.scale,ff.total,ff.sign
    c=ff.const;zero=ff.zeros[3]
    for k,x in enumerate(ff.radicals):assert sign(x,k)>0
    Q=[list(map(F,row))for row in old['form']]
    assert Q==[list(map(F,row))for row in d['form']]
    assert A.inertia(Q)==[20,1,0]
    rows=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={};dc={}
    def vp(raw,level):return tuple(ff.lift(parse(x,level),level)for x in raw)
    def dot(x,y):
        if (x,y)in dc:return dc[x,y]
        if y not in cache:cache[y]=tuple(total(scale(y[j],q)for j,q in row)for row in rows)
        v=total(mul(a,b)for a,b in zip(x,cache[y]));dc[x,y]=v;dc[y,x]=v;return v
    def lin(vs,cs):return tuple(total(mul(q,v[k])for v,q in zip(vs,cs))for k in range(21))
    def rat(x):
        for k in[3,2,1]:assert x[1]==ff.zeros[k-1];x=x[0]
        assert not any(x[1:]);return x[0]
    def eqraw(x,raw,level=3):assert x==ff.lift(parse(raw,level),level)
    N=[vp(v,1)for v in old['normals']]+[vp(pre['new_normal'],2)]
    R=[vp(v,3)for v in d['new_normals']];assert len(R)==2
    C=vp(d['strict_center'],2);E=vp(d['strict_escape']['point'],2)
    assert C==vp(pre['strict_center'],2)and E==vp(pre['strict_escape']['point'],2)
    rec=d['recipe'];S=rec['support'];assert S==pre['recipe']['right_support']
    G=[[rat(dot(N[i],N[j]))for j in S]for i in S];J=A.inv(G);piv=A.pd(G)
    assert min(x for row in J for x in row)>=0
    omega=lin([N[i]for i in S],[c(q)for q in J[S.index(61)]])
    assert omega==vp(rec['dual_weight'],2)
    assert dot(omega,omega)==c(F(2,3))
    assert all(dot(omega,N[i])==c(i==61)for i in S)
    X=vp(rec['X'],2);assert X==vp(pre['actual_edge']['point'],2)
    Av=scale(dot(X,X),-1);eqraw(Av,rec['A'],2);assert sign(Av)>0
    W=vp(rec['W'],2);B=dot(W,W);eqraw(B,rec['W_norm'],2);assert sign(B)>0
    assert dot(W,X)==zero and all(dot(W,N[i])==zero for i in S+[95])
    p=lin([omega,N[95]],[c(-1),c(F(-1,2))]);assert p==vp(rec['particular'],2)
    assert dot(p,p)==c(F(7,6))and dot(p,W)==dot(p,X)==zero
    m=F(rec['boost_integer']);assert m>0 and m.denominator==1
    Delta=mul(Av,add(scale(B,m*m),c(F(-5,6))))
    eqraw(Delta,rec['Delta'],2)
    assert Delta==ff.lift(ff.radicals[2],2)and sign(Delta)>0
    sqrtDelta=(ff.zeros[2],ff.const(1,2))
    norm=scale(mul(Av,Av),2);eqraw(norm,d['new_normal_norm'])
    for k,side in enumerate([1,-1]):
        expected=lin([p,W,X],[Av,scale(Av,side*m),sqrtDelta])
        assert R[k]==expected and dot(R[k],R[k])==norm
    norms=[c(q)for q in map(F,old['norms'])]+[c(2),norm,norm]
    allN=N+R
    assert all(dot(n,n)==q for n,q in zip(allN,norms))
    counts=Counter();newpairs=[]
    for k,r in enumerate(R):
        row=[]
        for i,n in enumerate(allN[:96+k]):
            g=dot(r,n);assert sign(g)<=0;prod=mul(norm,norms[i]);gg=mul(g,g)
            if g==zero:kind='right'
            elif gg==scale(prod,F(1,4)):kind='pi3'
            else:assert sign(add(gg,scale(prod,-1)))>0;kind='dotted'
            counts[kind]+=1;row.append(g)
            if i<96:eqraw(g,d['new_pairings'][k][i])
            else:eqraw(g,d['mutual_pairing'])
        assert [i for i,g in enumerate(row[:96])if g==zero]==[i for i in S if i!=61]
        assert dot(r,N[61])==dot(r,N[95])==scale(Av,-1)
        newpairs.append(row)
    assert dict(counts)==d['new_pair_counts']=={'right':34,'pi3':4,'dotted':155}
    expectedmutual=mul(mul(Av,Av),add(c(2),scale(B,-2*m*m)))
    assert dot(R[0],R[1])==expectedmutual
    assert sign(dot(C,C))<0 and all(sign(dot(C,n))<0 for n in allN)
    oldH=vp(old['strict_center'],1);assert sign(dot(C,oldH))<0
    for k in range(2):eqraw(dot(C,R[k]),d['new_center_pairings'][k])
    # Independently construct facet points without any field division.
    facets=[]
    for r in R:
        fp=lin([C,r],[norm,scale(dot(C,r),-1)])
        assert dot(fp,r)==zero and sign(dot(fp,fp))<0 and sign(dot(fp,C))<0
        assert all(sign(dot(fp,n))<0 for n in allN if n!=r);facets.append(fp)
    endpoints=[]
    for k,side in enumerate([1,-1]):
        P=lin([X,W],[scale(B,m),scale(sqrtDelta,side)])
        row=d['finite_endpoints'][k];assert P==vp(row['point'],3)
        q=dot(P,P);assert q==scale(mul(Av,B),F(-5,6));eqraw(q,row['norm'])
        assert sign(q)<0 and sign(dot(P,C))<0;eqraw(dot(P,C),row['future_pairing'])
        ids=S+[95,96+k];assert row['incident']==ids
        pp=[dot(P,n)for n in allN]
        assert all(sign(g)==(0 if i in ids else -1)for i,g in enumerate(pp))
        assert pp==[parse(g,3)for g in row['all98_pairings']]
        assert dot(P,R[1-k])==scale(mul(mul(Av,B),sqrtDelta),-2*m)
        # Normalize the final cap by positive A via exact cross identities.
        MG=[[rat(dot(allN[i],allN[j]))for j in ids[:-1]]for i in ids[:-1]]
        contact=[F(-1)if i in[61,95]else F(0)for i in ids[:-1]]
        normalized=[row+[x]for row,x in zip(MG,contact)]+[contact+[F(2)]]
        assert normalized==[list(map(F,row))for row in d['finite_endpoints'][k]['normalized_support_gram']]
        A.pd(normalized);assert F(2)-J[S.index(61)][S.index(61)]-F(1,2)==F(5,6)
        endpoints.append(P)
    assert endpoints[0]!=endpoints[1]
    assert d['finite_edge']['shared_support']==S+[95]
    assert d['finite_edge']['endpoint_walls']==[96,97]
    assert sign(dot(E,E))>0 and sign(dot(E,C))<0
    assert all(sign(dot(E,n))<0 for n in allN);eqraw(dot(E,E),d['strict_escape']['norm'],2)
    assert [dot(E,r)for r in R]==[parse(g,3)for g in d['strict_escape']['new_pairings']]
    cusp_pairs=[]
    oldcs=[old['completed_cusp']['point']]+[v['point']for v in old['preserved_old_bounded_cusps']]
    assert len(d['preserved_bounded_cusps'])==len(oldcs)==3
    for raw,row in zip(oldcs,d['preserved_bounded_cusps']):
        nu=tuple(ff.lift(tuple(map(F,x)),0)for x in raw)
        assert nu==vp(row['point'],2)and dot(nu,nu)==zero and sign(dot(nu,C))<0
        pairs=[dot(nu,r)for r in R];assert all(sign(q)<0 for q in pairs)
        assert pairs==[parse(g,3)for g in row['new_pairings']];cusp_pairs.append(pairs)
    print('PASS_193_NEW_PAIRS_AND_TWO_DISTINCT_PD20_ENDPOINTS',flush=True)
    # The other certified entire-plane face survives all three added walls.
    S75=face['faces'][1]['support'];assert S75==[75 if i==61 else i for i in S]
    G75=[[rat(dot(N[i],N[j]))for j in S75]for i in S75];A.pd(G75);J75=A.inv(G75)
    assert J75==[list(map(F,row))for row in face['faces'][1]['inverse']]
    hp=[dot(C,N[i])for i in S75]
    hc=[total(scale(x,q)for x,q in zip(hp,row))for row in J75]
    CP=lin([C]+[N[i]for i in S75],[c(1)]+[scale(q,-1)for q in hc])
    assert sign(dot(CP,CP))<0 and sign(dot(CP,C))<0
    extra=[]
    for i in[95,96,97]:
        pairs=[dot(allN[i],N[j])for j in S75]
        assert all(g==zero for j,g in zip(S75,pairs)if j!=75)
        cost=total(scale(mul(a,b),J75[k][l])for k,a in enumerate(pairs)for l,b in enumerate(pairs))
        q=add(norms[i],scale(cost,-1));h=dot(allN[i],CP)
        assert sign(q)<0 and sign(h)<0
        extra.append(dict(wall=i,projected_norm=q,future_face_center_pairing=h))
    print('PASS_OTHER_WHOLE_H2_FACE_UNCHANGED',flush=True)
    paths=[BASE,PREV,INPUT,'results/fresh_h20_bridge_faces_20261014_review.json',
           f'scripts/{TAG}.py','scripts/fresh_h20_bridge_faces_20261014_review.py',
           'scripts/fresh_h20_bridge_edge_closure_20261014.py']
    return A.enc(dict(status='PASS_INDEPENDENT_98_WALL_EDGE_CLOSURE_REVIEW',sources={p:sha(p)for p in paths},
        arithmetic='Independent sparse Fraction radical multiplication, rational support inversion and outward rational intervals; no producer arithmetic imports',
        maximum_sign_enclosure_bits=ff.bits,new_pairs_checked=193,new_pair_counts=dict(counts),all_norms_checked=98,
        support_PD_pivots=piv,endpoint_final_Schur_pivot=F(5,6),actual_distinct_finite_endpoints=2,
        endpoints_all_inequalities=196,new_facet_points=facets,strict_center_inequalities=98,
        strict_escape_inequalities=98,strict_escape_norm=dot(E,E),preserved_bounded_cusp_new_pairings=cusp_pairs,
        surviving_whole_H2_face=dict(support=S75,strict_face_center=CP,center_norm=dot(CP,CP),additional_projected_walls=extra,
            conclusion='All old95 projections were zero on support or future timelike. Each of the three subsequent wall projections is also strictly future timelike, so this entire H2 remains in the98 chamber.'),
        scope='Exact two-cap/193-pair replay, two actual finite endpoints of one old edge, actual new facets, three bounded cusp signs, and retained strict escape. Old96 pair/corner/link data and the116 old-vertex lower bound are inherited. Additionally proves S75 remains an entire H2. No complete98 vertex census or finite-volume claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
    p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
