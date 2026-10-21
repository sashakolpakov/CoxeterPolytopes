#!/usr/bin/env python3
"""Independent SymPy-quartic plus quadratic-extension physical replay.

No producer arithmetic imports. The earlier independent reviewer supplies the
quartic ANP field and rational outward sign bounds; new extension arithmetic is
implemented here. No global cone census is attempted.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_quartet_polar_extension_20261014_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_octic_escape_caps_20261015_review'
BASE='results/fresh_h9_fourth_escape_extension_20261015.json'
PRIMARY='results/fresh_h9_octic_escape_caps_20261015.json'
RECIPES='results/fresh_h9_after_fourth_rescore_20261015.json'
GATE='results/fresh_h9_octic_angle_20261015.json'
K,Z,O=R.K,R.Z,R.O
DELTA=None

def base(a):return a if isinstance(a,K.dtype)else K.convert(a)
class X:
    def __init__(self,a=0,b=0):
        if isinstance(a,X):self.a,self.b=a.a,a.b
        else:self.a,self.b=base(a),base(b)
    def __add__(self,o):
        o=X(o);return X(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return X(-self.a,-self.b)
    def __sub__(self,o):return self+-X(o)
    def __rsub__(self,o):return X(o)+-self
    def __mul__(self,o):
        o=X(o)
        if not self.b:return X(self.a*o.a,self.a*o.b)
        if not o.b:return X(self.a*o.a,self.b*o.a)
        return X(self.a*o.a+DELTA*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=X(o)
        if not o.b:return X(self.a/o.a,self.b/o.a)
        n=o.a*o.a-DELTA*o.b*o.b
        return self*X(o.a/n,-o.b/n)
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,o):
        o=X(o);return self.a==o.a and self.b==o.b
    def sign(self):
        a,b=R.sg(self.a),R.sg(self.b)
        if not a:return b
        if not b or a==b:return a
        return a*R.sg(self.a*self.a-DELTA*self.b*self.b)
def dec(x):return X(R.dec(x[0]),R.dec(x[1]))
def dv(x):return list(map(dec,x))
def emb(v):return [X(x)for x in R.dv(v)]
def dot(v,w):return sum((x*y for x,y in zip(v,w)if x and y),X())
def mv(M,v):return[dot(row,v)for row in M]
def val(x):return [list(map(str,x.a.to_list())),list(map(str,x.b.to_list()))]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def pd(M):
    a=[row[:]for row in M];pivots=[]
    for k in range(len(a)):
        p=a[k][k];assert p.sign()>0;pivots.append(val(p))
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]-=a[i][k]*a[k][j]/p
    return pivots
def relation(g,q,label):
    if not g:return 'right',2,None
    assert g.sign()<0
    c=g*g/(2*q)
    if (c-1).sign()>=0:return ('parallel'if c==1 else'dotted'),None,c
    phi=(O+R.BASIS[1])/2
    vals={3:X(F(1,4)),4:X(F(1,2)),5:X((3*O+R.BASIS[1])/8),
          6:X(F(3,4)),10:X((5*O+R.BASIS[1])/8)}
    if label==3:vals[20]=X(F(1,2),O/(8*phi*phi))
    hits=[m for m,v in vals.items()if c==v];assert len(hits)==1
    return 'finite',hits[0],c

def run():
    global DELTA
    b,d,recipes,gate=map(read,[BASE,PRIMARY,RECIPES,GATE])
    for doc in[b,d,recipes,gate]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert gate['complete_finite_alphabets']['delta3']==[2,3,4,5,6,10,20]
    assert gate['complete_finite_alphabets']['delta5']==[2,3,4,5,6,10]
    BQ,BN,BH,BE=R.dm(b['form']),R.dm(b['normals']),R.dv(b['center']),R.dv(b['strict_escape']['point'])
    bp=lambda v,w:R.C.pair(v,w,BQ)
    es=d['common_strict_escape'];EE=R.dv(es['point'])
    assert bp(EE,EE)==R.dec(es['norm'])and R.sg(bp(EE,EE))>0
    assert bp(EE,BH)==R.dec(es['center_pairing'])and R.sg(bp(EE,BH))<0
    olde=[bp(EE,n)for n in BN];assert olde==R.dv(es['old_wall_pairings'])and all(R.sg(x)<0 for x in olde)
    out=[];base_recipes=[]
    for row,recipe in zip(d['caps'],recipes['exact_nine_contact_recipes']):
        label=row['support_contact_label'];assert label==recipe['support_contact_label']
        DELTA=R.dec(row['delta']);assert DELTA==R.dec(recipe['delta'])and R.sg(DELTA)>0
        pp,vv=R.dv(recipe['particular_p']),R.dv(recipe['kernel_v'])
        assert not bp(pp,vv)and bp(pp,pp)+DELTA*bp(vv,vv)==2*O
        n=[X(a,-v)for a,v in zip(pp,vv)];assert n==dv(row['normal'])
        Q=[[X(x)for x in rr]for rr in BQ];N=[[X(x)for x in rr]for rr in BN]
        H,E,En=[X(x)for x in BH],[X(x)for x in BE],[X(x)for x in EE]
        pair=lambda v,w:dot(v,mv(Q,w))
        assert pair(n,n)==2 and dec(row['norm'])==2
        norms=[pair(x,x)for x in N];assert norms==[X(x)for x in R.dv(b['norms'])]
        g=[pair(n,x)for x in N];counts=Counter();finite=[]
        for i,x in enumerate(g):
            kind,m,c=relation(x,norms[i],label);saved=row['relations'][i]
            assert (kind,m)==(saved['kind'],saved['label'])and x==dec(saved['pairing'])
            assert c==(None if saved['invariant']is None else dec(saved['invariant']))
            counts[kind+str(m)]+=1
            if kind in('right','finite'):finite.append([i,m])
        assert dict(counts)==row['relation_counts']
        for point,field,sgn in [(H,'center_pairing',-1),(E,'old_escape_pairing',1),(En,'common_escape_pairing',-1)]:
            x=pair(n,point);assert x==dec(row[field])and x.sign()==sgn
        assert pair(H,H).sign()<0 and all(pair(H,x).sign()<0 for x in N)
        allN=N+[n];cov=[mv(Q,x)for x in allN]
        assert len(row['actual_facet_witnesses'])==30
        for i,saved in enumerate(row['actual_facet_witnesses']):
            nn=allN[i];P=[a-pair(H,nn)/pair(nn,nn)*z for a,z in zip(H,nn)]
            assert P==dv(saved['point'])and pair(P,P)==dec(saved['norm'])and pair(P,P).sign()<0
            assert pair(P,H).sign()<0
            assert all(dot(P,c).sign()==(0 if j==i else-1)for j,c in enumerate(cov))
        saved=row['finite_endpoint'];V=dv(saved['point']);S=saved['support']
        assert S==recipe['support']and pair(V,V)==dec(saved['norm'])and pair(V,V).sign()<0
        assert pair(V,H).sign()<0 and not pair(V,n)
        vp=[pair(V,x)for x in N];assert vp==dv(saved['pairings'])
        assert all(x.sign()==(0 if i in S else-1)for i,x in enumerate(vp))
        SS=[n]+[N[i]for i in reversed(S)]
        pivots=pd([[pair(x,y)for y in SS]for x in SS])
        out.append(dict(label=label,relation_counts=dict(counts),finite_neighbors=finite,
            endpoint_support=S+[29],reverse_order_PD_pivots=pivots,
            center_pairing=val(pair(n,H)),old_escape_pairing=val(pair(n,E)),common_escape_pairing=val(pair(n,En)),
            actual_facets=30,facet_sign_checks=900,endpoint_sign_checks=30))
        base_recipes.append((pp,vv,DELTA))
        print('independent octic',label,'PASS',dict(counts),flush=True)
    p0,v0,d0=base_recipes[0];p1,v1,d1=base_recipes[1]
    assert not bp(p0,v1)and not bp(v0,p1)
    DELTA=d0*d1;mut=d['mutual_pair'];assert DELTA==R.dec(mut['radical_square'])
    x=X(bp(p0,p1),bp(v0,v1));assert x==dec(mut['pairing'])and x.sign()==1
    paths=[BASE,PRIMARY,RECIPES,GATE,f'scripts/{TAG}.py','scripts/fresh_h9_octic_escape_caps_20261015.py',
        'scripts/fresh_h9_quartet_polar_extension_20261014_review.py',
        'scripts/fresh_h9_h8_all_high_20261014_review.py','scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_TWO_SEPARATE_OCTIC_CAPS_AND_COMMON_ESCAPE',
        sources={p:sha(p)for p in paths},arithmetic='Independent SymPy quartic ANP field/192-bit outward signs plus new quadratic extension, no producer arithmetic imports.',
        caps=out,physical_new_pairs_checked=58,facet_sign_checks=1800,endpoint_sign_checks=60,
        common_old_escape_signs=29,mutual_positive_pairing=val(x),
        common_escape_norm=list(map(str,bp(EE,EE).to_list())),
        scope='Both separate30-wall chambers have actual facets and the specified actual finite endpoint, but a strict common spacelike escape proves each infinite-volume. The two new walls cannot coexist acutely. No full30-wall cone or complete vertex census repeated.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run()
    out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('PASS saved independent octic cap review')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
