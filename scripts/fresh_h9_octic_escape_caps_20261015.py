#!/usr/bin/env python3
"""Exact checks of two norm-two proposals over quadratic extensions of K."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_quartet_polar_second_20261014 as B

A,D,C=B.A,B.D,B.C
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_octic_escape_caps_20261015'
BASE='results/fresh_h9_fourth_escape_extension_20261015.json'
RECIPES='results/fresh_h9_after_fourth_rescore_20261015.json'
DELTA=None

class X:
    __slots__=('a','b')
    def __init__(self,a=0,b=0):
        if isinstance(a,X):self.a,self.b=a.a,a.b
        else:self.a,self.b=A.L(a),A.L(b)
    def __add__(self,y):
        y=X(y);return X(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return X(-self.a,-self.b)
    def __sub__(self,y):return self+-X(y)
    def __rsub__(self,y):return X(y)+-self
    def __mul__(self,y):
        y=X(y)
        if not self.b:return X(self.a*y.a,self.a*y.b)
        if not y.b:return X(self.a*y.a,self.b*y.a)
        return X(self.a*y.a+DELTA*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=X(y)
        if not y.b:return X(self.a/y.a,self.b/y.a)
        q=y.a*y.a-DELTA*y.b*y.b;assert q
        return self*X(y.a/q,-y.b/q)
    def __rtruediv__(self,y):return X(y)/self
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,y):
        y=X(y);return self.a==y.a and self.b==y.b
    def __hash__(self):return hash((self.a,self.b))
    def sign(self,embedding=0):
        assert embedding==0
        a,b=self.a.sign(),self.b.sign()
        if not a:return b
        if not b or a==b:return a
        return a*(self.a*self.a-DELTA*self.b*self.b).sign()
    def enc(self):return [self.a.enc(),self.b.enc()]

Z=X();ONE=X(1)
def enc(x):
    if isinstance(x,X):return x.enc()
    if isinstance(x,(A.L,A.A.Q,F)):return A.enc(x)
    if isinstance(x,(list,tuple)):return [enc(y)for y in x]
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    return x
def dot(v,w):return sum((a*b for a,b in zip(v,w)),Z)
def mv(M,v):return [dot(row,v)for row in M]
def pair(v,w,Q):return dot(v,mv(Q,w))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def setup(d,delta):
    global DELTA
    A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    DELTA=D.dv([delta])[0];assert DELTA.sign()>0
    C.Golden=X;C.ZERO=Z
def embedv(v):return [X(a)for a in D.dv(v)]
def embedm(M):return [embedv(r)for r in M]
def relation(g,q,r,label):
    if not g:return 'right',2,None
    if g.sign()>0:return 'obtuse',None,None
    c=g*g/(q*r)
    if (c-1).sign()>=0:return ('parallel'if c==1 else'dotted'),None,c
    vals={3:X(F(1,4)),4:X(F(1,2)),5:X(A.A.Q(F(3,8),F(1,8))),
          6:X(F(3,4)),10:X(A.A.Q(F(5,8),F(1,8)))}
    if label==3:
        phi=A.L(A.A.Q(F(1,2),F(1,2)))
        assert DELTA==phi*phi*phi*phi*A.L(A.A.Q(10,2))
        vals[20]=X(F(1,2),1/(8*phi*phi))
    for m,v in vals.items():
        if c==v:return 'finite',m,c
    return 'unclassified',None,c

def run():
    global DELTA
    d,recipes=read(BASE),read(RECIPES)
    for doc in [d,recipes]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    C.Golden=A.L;C.ZERO=A.Z
    Q0,N0,H0=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center'])
    active=[0,1,2,3,5,6,8,10,28]
    assert active in d['complete_census']['hyperideal_rays']
    ray,_,_=B.ray(tuple(active),[A.mv(Q0,n)for n in N0])
    epsilon=F(1,1024);common=[x+A.L(epsilon)*h for x,h in zip(ray,H0)]
    common_norm=A.pair(common,common,Q0);common_pairs=[A.pair(common,n,Q0)for n in N0]
    assert common_norm.sign()>0 and A.pair(common,H0,Q0).sign()<0
    assert all(g.sign()<0 for g in common_pairs)
    rows=[]
    for row in recipes['exact_nine_contact_recipes']:
        label=row['support_contact_label'];setup(d,row['delta'])
        Q,N,H,E=embedm(d['form']),embedm(d['normals']),embedv(d['center']),embedv(d['strict_escape']['point'])
        p,v=D.dv(row['particular_p']),D.dv(row['kernel_v'])
        n=[X(a,-b)for a,b in zip(p,v)]
        assert pair(n,n,Q)==2
        norms=[X(a)for a in D.dv(d['norms'])]
        vals=[pair(n,r,Q)for r in N]
        relations=[]
        for i,g in enumerate(vals):
            kind,m,c=relation(g,X(2),norms[i],label)
            relations.append(dict(wall=i,kind=kind,label=m,invariant=c,pairing=g))
        bad=[r for r in relations if r['kind']in ['unclassified','obtuse']]
        assert not bad,enc(bad)
        assert pair(n,H,Q).sign()<0 and pair(n,E,Q).sign()>0
        new_escape_pairing=pair(n,[X(a)for a in common],Q)
        assert new_escape_pairing.sign()<0
        S=row['support'];cov=[mv(Q,r)for r in N]
        V=C._full_rank_kernel([cov[i]for i in S]+[mv(Q,n)])
        assert V is not None
        if pair(V,H,Q).sign()>0:V=[-a for a in V]
        vp=[pair(V,r,Q)for r in N]
        assert pair(V,V,Q).sign()<0 and not pair(V,n,Q)
        assert all(a.sign()==(0 if i in S else -1)for i,a in enumerate(vp))
        NN=N+[n];qq=norms+[X(2)];facets=[]
        for i,r in enumerate(NN):
            t=pair(H,r,Q)/qq[i];w=[a-t*b for a,b in zip(H,r)]
            wv=[pair(w,u,Q)for u in NN]
            assert pair(w,w,Q).sign()<0 and not wv[i]
            assert all(a.sign()<0 for j,a in enumerate(wv)if i!=j)
            facets.append(dict(index=i,point=w,norm=pair(w,w,Q)))
        rows.append(enc(dict(support_contact_label=label,delta=DELTA,normal=n,norm=X(2),
            relations=relations,relation_counts=dict(Counter(r['kind']+str(r['label'])for r in relations)),
            center_pairing=pair(n,H,Q),old_escape_pairing=pair(n,E,Q),
            common_escape_pairing=new_escape_pairing,
            actual_facet_witnesses=facets,
            finite_endpoint=dict(point=V,norm=pair(V,V,Q),support=S,pairings=vp),
            field_encoding='[a,b] means a+b sqrt(delta), each a,b in the inherited quartic encoding')))
        print('exact cap label',label,'PASS',rows[-1]['relation_counts'],flush=True)
    rr=recipes['exact_nine_contact_recipes'];ps=[D.dv(r['particular_p'])for r in rr]
    vs=[D.dv(r['kernel_v'])for r in rr];assert vs[0]==vs[1]
    assert all(not A.pair(p,vs[0],Q0)for p in ps)
    DELTA=D.dv([rr[0]['delta']])[0]*D.dv([rr[1]['delta']])[0]
    mutual=X(A.pair(ps[0],ps[1],Q0),A.pair(vs[0],vs[0],Q0))
    assert DELTA==A.L(A.A.Q(6115,2735))
    assert mutual.a==A.L(A.A.Q(-8,-4))
    assert mutual.b==A.L(A.A.Q(F(-1,2),F(3,10))) and mutual.sign()>0
    paths=[BASE,RECIPES,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py',
           'scripts/fresh_h9_h8_disjoint_walk_20261014_next_frontier.py',
           'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py','src/coxeter_polytopes/verify/chamber.py']
    return dict(status='PASS_TWO_SEPARATE_EXACT_OCTIC_ESCAPE_CAPS',sources={p:sha(p)for p in paths},
        caps=rows,mutual_pair=enc(dict(radical_square=DELTA,pairing=mutual,sign=1,
            conclusion='The two outward normals have positive pairing and cannot be retained together in an acute Coxeter chamber.')),
        common_strict_escape=A.enc(dict(point=common,norm=common_norm,old_wall_pairings=common_pairs,
            center_pairing=A.pair(common,H0,Q0),active_ray=active,epsilon=epsilon,
            encoding='Inherited quartic encoding; both added-wall pairings are recorded in the separate cap rows.')),
        scope='Each cap separately is exactly compatible with all29 retained walls, cuts their saved escape, retains the center and gives30 actualfacets and a finite endpoint. A common strictescape proves both resultingchambers have infinitevolume. The two additions are mutually incompatible. No finite-volume example claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==read(str(out.relative_to(ROOT)))
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
