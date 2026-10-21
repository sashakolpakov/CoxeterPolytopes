#!/usr/bin/env python3
"""Exact joint cap selection and first endpoints of the forty original edges.

The bounded numerical search proposes recipes and branches only. Every claimed
normal, angle, endpoint and escape is checked using exact ordered radicals.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import argparse,hashlib,json
from fresh_h9_makarov_face_extension_20261013 import K,Z,O,PHI,inv,inertia,mv,pair

ROOT=Path(__file__).resolve().parents[1];P='fresh_h9_makarov_joint_20261013'
SOURCE=ROOT/'results/fresh_h9_makarov_face_extension_20261013.json'
OUT=ROOT/f'results/{P}_exact.json'
def dec(v):return K(*v)
def mul(a,b,ds):
    n=len(a);out=[Z]*n
    for i,x in enumerate(a):
        if not x:continue
        for j,y in enumerate(b):
            if not y:continue
            t=x*y
            for k,d in enumerate(ds):
                if (i&j)>>k&1:t*=d
            out[i^j]+=t
    return tuple(out)
def sign(c,ds):
    if not ds:return c[0].sign()
    h=len(c)//2;a,b=c[:h],c[h:];sa,sb=sign(a,ds[:-1]),sign(b,ds[:-1])
    if not sb:return sa
    if not sa:return sb
    if sa==sb:return sa
    aa=mul(a,a,ds[:-1]);bb=mul(b,b,ds[:-1])
    return sa*sign(tuple(x-ds[-1]*y for x,y in zip(aa,bb)),ds[:-1])
@dataclass(frozen=True)
class X:
    c:tuple
    ds:tuple
    def co(self,o):
        if isinstance(o,X):assert o.ds==self.ds;return o
        if not isinstance(o,K):o=K(o)
        return X((o,)+(Z,)*(len(self.c)-1),self.ds)
    def __add__(self,o):
        o=self.co(o);return X(tuple(a+b for a,b in zip(self.c,o.c)),self.ds)
    __radd__=__add__
    def __neg__(self):return X(tuple(-a for a in self.c),self.ds)
    def __sub__(self,o):return self+-self.co(o)
    def __rsub__(self,o):return self.co(o)+-self
    def __mul__(self,o):return X(mul(self.c,self.co(o).c,self.ds),self.ds)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=self.co(o)
        if not any(o.c[1:]):return X(tuple(a/o.c[0]for a in self.c),self.ds)
        assert len(self.ds)==1
        a,b=o.c;den=a*a-self.ds[0]*b*b;assert den
        return self*X((a/den,-b/den),self.ds)
    def sign(self):return sign(self.c,self.ds)
    def enc(self):return [a.enc()for a in self.c]
    def approx(self):return sum(a.approx()*__import__('math').prod(d.approx()**.5 for k,d in enumerate(self.ds)if i>>k&1)for i,a in enumerate(self.c))
def E(a=Z,b=Z,d=K(2)):return X((a,b),(d,))
def ep(x,y,Q):return sum((x[i]*Q[i][j]*y[j]for i in range(len(x))for j in range(len(y))),x[0]*0)
def cross(v,d,w,e,Q):
    out=[Z]*4
    for i in range(len(v)):
        for j in range(len(w)):
            for a in range(2):
                for b in range(2):out[a+2*b]+=v[i].c[a]*Q[i][j]*w[j].c[b]
    return X(tuple(out),(d,e))
def fsqrt(q):
    if q<0:return None
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    return F(a,b)if a*a==q.numerator and b*b==q.denominator else None
def ksqrt(d):
    n=fsqrt(d.a*d.a-5*d.b*d.b)
    if n is None:return None
    for s in [n,-n]:
        a=fsqrt((d.a+s)/2)
        if a:
            x=K(a,d.b/(2*a))
            if x*x==d:return x if x.sign()>0 else -x
    if not d.b:
        b=fsqrt(d.a/5)
        if b:return K(0,b)
    return None
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def reconstruct(recipe,Q,W):
    if 'old_support'in recipe:
        ids=recipe['old_support']+list(range(7,11));labs=recipe['old_labels']+[2]*4
        for j,m in zip(recipe['transverse_support'],recipe['transverse_labels']):labs[5+j]=m
    else:
        ids=recipe['support']+[recipe['anchor']];lookup=dict(recipe['contacts']);lookup[recipe['anchor']]=recipe['anchor_label'];labs=[lookup.get(i,2)for i in ids]
    assert len(ids)==len(set(ids))==9
    rhs=[]
    for i,m in zip(ids,labs):
        assert m in [2,3,5] and(i!=6 or m==2)
        rhs.append({2:Z,3:-O,5:-PHI}[m])
    M=[mv(Q,W[i])for i in ids];chosen=None
    for free in range(10):
        try:
            B=inv(M+[[K(int(j==free))for j in range(10)]])
            chosen=free;break
        except (StopIteration,AssertionError):pass
    assert chosen is not None
    r0=mv(B,rhs+[Z]);z=[row[-1]for row in B]
    aa=pair(z,z,Q);bb=2*pair(z,r0,Q);cc=pair(r0,r0,Q)-2
    disc=bb*bb-4*aa*cc;assert aa and disc.sign()>0
    a=[r0[i]-bb*z[i]/(2*aa)for i in range(10)];b=[z[i]/(2*aa)for i in range(10)]
    proposals=[[E(x,s*y,disc)for x,y in zip(a,b)]for s in [-1,1]]
    distances=[sum((x.approx()-v)**2 for x,v in zip(p,recipe['normal']))for p in proposals]
    k=min(range(2),key=lambda i:distances[i]);assert distances[k]<1e-10 and distances[1-k]>1e-6
    r=proposals[k];sqrt=ksqrt(disc)
    if sqrt is not None:r=[E(x.c[0]+sqrt*x.c[1])for x in r];disc=K(2)
    assert ep(r,r,Q).sign()>0 and(ep(r,r,Q)-2).sign()==0
    for i,q in zip(ids,rhs):assert(ep(r,[E(x,d=disc)for x in W[i]],Q)-q).sign()==0
    return r,disc,{'constraint_walls':ids,'constraint_labels':labs,'linear_free_coordinate':chosen,'quadratic_sign':[-1,1][k]}

def build():
    src=json.loads(SOURCE.read_text());Q=[[dec(v)for v in row]for row in src['form']];W=[[dec(v)for v in row]for row in src['seed_normals']];center=[dec(v)for v in src['strict_center']]
    assert inertia(Q)==[9,1]
    choice=ROOT/f'results/{P}_cover_simple.json';sel=json.loads(choice.read_text());banks={s:json.loads((ROOT/f'results/{P}{s}.json').read_text())for s in ['', '_profiles']}
    caps=[];recipes=[];fields=[];rows=[]
    for suffix,k in sel['references']:
        recipe=banks[suffix]['caps'][k];r,d,meta=reconstruct(recipe,Q,W);caps.append(r);fields.append(d);recipes.append(recipe)
        rows.append({'source_reference':[suffix,k],'recipe':recipe,'reconstruction':meta,'radicand':d.enc(),'normal':[x.enc()for x in r]})
    # Write the formulas early for a separate arithmetic checker.
    early={'status':'EXACT_NORMALS_PENDING_FULL_GEOMETRY','form':src['form'],'seed_normals':src['seed_normals'],'strict_center':src['strict_center'],'caps':rows}
    OUT.with_name(OUT.stem+'_normals.json').write_text(json.dumps(early,indent=2)+'\n')
    nodes=[[E(x)for x in w]for w in W]+caps;ds=[K(2)]*11+fields
    norms=[pair(w,w,Q)for w in W]+[K(2)]*len(caps)
    cosine={2:Z,3:K(F(1,4)),4:K(F(1,2)),5:(PHI+1)/4,6:K(F(3,4)),10:(PHI+2)/4}
    gates=[];counts={}
    for i in range(len(nodes)):
        for j in range(i):
            p=cross(nodes[i],ds[i],nodes[j],ds[j],Q);assert p.sign()<=0
            q=p*p/(norms[i]*norms[j])
            if(q-1).sign()>0:tag='dotted'
            elif(q-1).sign()==0:tag='parallel'
            else:
                lab=next((m for m,v in cosine.items()if(q-v).sign()==0),None);assert lab is not None,(i,j,p.approx())
                tag='m'+str(lab)
            counts[tag]=counts.get(tag,0)+1;gates.append({'pair':[j,i],'type':tag,'pairing':p.enc()})
    centers=[]
    for r,d in zip(nodes,ds):
        p=ep(r,[E(x,d=d)for x in center],Q);assert p.sign()<0;centers.append(p.enc())
    # Every inequality is a true facet: project the strict timelike center to
    # its mirror. Acuteness makes every other inequality strictly negative.
    assert pair(center,center,Q).sign()<0
    oldvertices=[]
    for row in src['compact_H5_vertices']:
        v=[dec(x)for x in(row['point']if isinstance(row,dict)else row)]
        if len(v)==6:v+= [Z]*4
        oldvertices.append(v)
    assert len(oldvertices)==10
    endpoints=[];old_preservation=[]
    H=[row[6:]for row in Q[6:]];B=inv(H)
    for vi,v in enumerate(oldvertices):
        assert pair(v,v,Q).sign()<0
        assert all(pair(v,w,Q).sign()<=0 for w in W)
        oldactive=[i for i,w in enumerate(W)if not pair(v,w,Q)];assert len(oldactive)==9
        ps=[ep(r,[E(x,d=d)for x in v],Q)for r,d in zip(caps,fields)];assert all(x.sign()<0 for x in ps)
        old_preservation.append([x.enc()for x in ps])
        for j in range(4):
            support=[i for i in oldactive if i!=7+j];dual=[Z]*6+[B[i][j]for i in range(4)]
            D=pair(dual,dual,Q);G=[[pair(W[a],W[b],Q)for b in support]for a in support];assert inertia(G)==[8,0];I=inv(G)
            candidates=[]
            for k,(r,d)in enumerate(zip(caps,fields)):
                a=ep(r,[E(x,d=d)for x in v],Q);b=ep(r,[E(x,d=d)for x in dual],Q)
                if b.sign()>=0:continue
                t=a/b
                if t.sign()>0 and(t*t*D+pair(v,v,Q)).sign()<0:candidates.append((t.approx(),k,t))
            assert candidates,(vi,j)
            _,owner,t=min(candidates);d=fields[owner];point=[E(v[k],d=d)-t*dual[k]for k in range(10)]
            q=ep(point,point,Q);assert q.sign()<0
            assert ep(point,[E(x,d=d)for x in center],Q).sign()<0
            values=[cross(point,d,r,e,Q)for r,e in zip(nodes,ds)];assert all(x.sign()<=0 for x in values)
            zeros=[i for i,x in enumerate(values)if x.sign()==0];assert zeros==sorted(support+[11+owner]),(vi,j,zeros)
            g=[ep(caps[owner],[E(x,d=d)for x in W[i]],Q)for i in support]
            schur=E(K(2),d=d)-sum((g[a]*I[a][b]*g[b]for a in range(8)for b in range(8)),E(d=d));assert schur.sign()>0
            endpoints.append({'old_vertex':vi,'released_H4_wall':7+j,'support':support,'new_owner':11+owner,
                              'radicand':d.enc(),'edge_parameter':t.enc(),'point':[x.enc()for x in point],
                              'norm':q.enc(),'incident':zeros,'schur_complement':schur.enc(),'all_wall_pairings':[x.enc()for x in values]})
    assert len({tuple(row['incident'])for row in endpoints})==40
    escape=[140532998,178861908,130353897,81870848,33432956,12307295,-284558239,-345319008,-230212893,-115106557]
    ev=[K(x)for x in escape];en=pair(ev,ev,Q);et=pair(ev,center,Q);assert en.sign()>0 and et.sign()<0
    eps=[ep(r,[E(x,d=d)for x in ev],Q)for r,d in zip(nodes,ds)];assert all(x.sign()<0 for x in eps)
    out={'status':'EXACT_27_FACET_COXETER_CHAMBER_40_FINITE_OLD_EDGE_ENDPOINTS_WITH_STRICT_ESCAPE',
         'sources':{str(p.relative_to(ROOT)):sha(p)for p in [SOURCE,ROOT/'scripts/fresh_h9_makarov_face_extension_20261013.py',choice,ROOT/f'results/{P}.json',ROOT/f'results/{P}_profiles.json']},
         'base_field':'Q(sqrt(5)); each cap uses one additional positive square root',
         'form':src['form'],'seed_normals':src['seed_normals'],'strict_center':src['strict_center'],'caps':rows,
         'facet_count':27,'pair_counts':counts,'all_pair_gates':gates,'strict_center_pairings':centers,
         'old_vertex_preservation':old_preservation,'new_finite_endpoints':endpoints,
         'escape':{'point':escape,'norm':en.enc(),'center_pairing':et.enc(),'all_wall_pairings':[x.enc()for x in eps]},
         'scope':'All27 actual facets, the10 original finite vertices and40 new first finite endpoints are certified. This is not a full vertex census. Strict spacelike escape proves infinite volume; no compact H9 example.'}
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();result=build()
    if a.verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',result['facet_count'],result['pair_counts'],len(result['new_finite_endpoints']),flush=True)
