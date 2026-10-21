#!/usr/bin/env python3
"""Independent exact saved two-ray certificate checks, no producer imports.

Frozen original cone rays/field geometry remain inputs. Rechecks the full
adjacency domain, sharp cone bounds and every unique anchor solution;
uses completed-square norm evaluation rather than the producer formula.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_two_ray_completion_preflight_20261018'
OLD='fresh_raw23_union_polar_completion_preflight_20261018'
def read(t):return json.loads((ROOT/f'results/{t}.json').read_text())
def bind(d):
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
class Q:
    def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
    def __add__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,o):return self+-o
    def __mul__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return Q(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        if not isinstance(o,Q):o=Q(o)
        n=o.a*o.a-5*o.b*o.b;assert n
        return self*Q(o.a/n,-o.b/n)
    def sign(self):
        if not self.b:return (self.a>0)-(self.a<0)
        if not self.a:return (self.b>0)-(self.b<0)
        if self.a*self.b>0:return (self.a>0)-(self.a<0)
        n=self.a*self.a-5*self.b*self.b
        return ((self.a>0)-(self.a<0))*((n>0)-(n<0))
    def __eq__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return self.a==o.a and self.b==o.b

def mat(x):return [[Q(*y)for y in row]for row in x]
def principal(g,ids):return [[g[i][j]for j in ids]for i in ids]
def positive(g):
    a=[r[:]for r in g]
    for k in range(len(a)):
        if a[k][k].sign()<=0:return False
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]=a[i][j]-a[i][k]*a[k][j]/a[k][k]
    return True

def inverse(g):
    n=len(g);a=[row[:]+[Q(int(i==j))for j in range(n)]for i,row in enumerate(g)]
    for i in range(n):
        k=next(k for k in range(i,n)if a[k][i]!=0);a[i],a[k]=a[k],a[i]
        z=a[i][i];a[i]=[x/z for x in a[i]]
        for k in range(n):
            if k!=i:
                z=a[k][i];a[k]=[x-z*y for x,y in zip(a[k],a[i])]
    return [row[n:]for row in a]


def rational_root(x):
    from math import isqrt
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return F(a,b)if a*a==x.numerator and b*b==x.denominator else None

def root(q):
    if q.sign()<0:return None
    if q==0:return Q()
    if q.b==0:
        d=rational_root(q.a/5)
        if d is not None:return Q(0,d)
    z=rational_root(q.a*q.a-5*q.b*q.b)
    if z is None:return None
    for a in((q.a+z)/2,(q.a-z)/2):
        c=rational_root(a)
        if c is None or not c:continue
        r=Q(c,q.b/(2*c))
        if r*r==q:return r if r.sign()>0 else -r
    return None

class Field:
    """Self-contained ordered multiquadratic arithmetic from the saved-checker design."""
    def __init__(self,radicands):
        self.basis=[];self.powers=[Q(1)]
        for r in radicands:
            assert r.sign()>0
            if any(root(r/s)is not None for s in self.powers):continue
            self.basis.append(r);self.powers +=[s*r for s in self.powers]
        assert len(self.basis)<=3
        self.n=len(self.powers)
    def scalar(self,q):return [q if isinstance(q,Q)else Q(q)]+[Q()for _ in range(self.n-1)]
    def sqrt(self,q):
        if q==0:return self.scalar(0)
        for i,s in enumerate(self.powers):
            z=root(q/s)
            if z is not None:
                a=self.scalar(0);a[i]=z;return a
        raise AssertionError('missing squareclass')
    def mul(self,a,b):
        c=self.scalar(0)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i^j]=c[i^j]+x*y*self.powers[i&j]
        return c
    def sign(self,a):
        def rec(z,k):
            if k==0:return z[0].sign()
            m=1<<(k-1);x=z[:m];y=z[m:];sx,sy=rec(x,k-1),rec(y,k-1)
            if not sx:return sy
            if not sy or sx==sy:return sx
            xx=self.mul(x+[Q()]* (self.n-m),x+[Q()]* (self.n-m))[:m]
            yy=self.mul(y+[Q()]* (self.n-m),y+[Q()]* (self.n-m))[:m]
            norm=[a-self.basis[k-1]*b for a,b in zip(xx,yy)]
            return sx*rec(norm,k-1)
        return rec(a,len(self.basis))

def rank(rows):
    if not rows:return 0
    a=[r[:]for r in rows];k=0
    for j in reversed(range(len(a[0]))):
        z=next((i for i in range(k,len(a))if a[i][j]!=0),None)
        if z is None:continue
        a[k],a[z]=a[z],a[k];p=a[k][j]
        a[k]=[x/p for x in a[k]]
        for i in range(k+1,len(a)):
            z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k

def vector(x):return[Q(*q)for q in x]
def fd(x,f):
    assert len(x)==f.n
    return vector(x)
def add(a,b):return[x+y for x,y in zip(a,b)]
def mulq(a,q):return[x*q for x in a]
def inner(x,g,y):return sum((x[i]*g[i][j]*y[j]for i in range(len(x))for j in range(len(y))),Q())

def main():
    d,old=read(TAG),read(OLD);bind(d);bind(old)
    assert d['status']=='PASS_RAW23_ADJACENT_SPACELIKE_TWO_RAY_PREFLIGHT'
    assert d['relative_interiors_only']and d['boundaries_inherited_from_frozen_polar']
    assert not d['timelike_generator_faces_tested']and not d['higher_dimensional_faces_tested']
    assert d['arbitrary_passive_lifts_relaxed']and not d['unknown_passive_contacts_assumed_zero']
    assert not d['geometric_realization_claimed']and not d['pair_exclusion_claimed']
    assert d['exported_columns']==[]and d['optimizer_invocations']==0
    assert d['candidates_reaching_other_native_or_conjugate_gates']==0 and not d['retained_candidates']
    assert len(d['cases'])==len(old['cases'])==3
    hist=Counter();paircount=facecount=tests=0
    native=(Q(1),Q(2),Q('3/2','1/2'),Q(3),Q('5/2','1/2'))
    for c,o in zip(d['cases'],old['cases']):
        assert c['kinds']==o['kinds']
        cone=o['initial_cone'];g=mat(o['metric']);rows=mat(cone['normal_rows']);sc=vector(o['original_squared_scales'])
        rays={r['index']:r for r in cone['rays']};outside=[r for r in rays.values()if r['norm_sign']>0]
        pairs=list(combinations(outside,2));assert len(pairs)==len(c['adjacency_decisions'])==496
        faces={tuple(z['rays']):z for z in c['faces']};assert len(faces)==len(c['faces'])==112
        seen=set()
        for (x,y),dec in zip(pairs,c['adjacency_decisions']):
            key=(x['index'],y['index']);paircount+=1
            common=sorted(set(x['zero_facets'])&set(y['zero_facets']))
            r=rank([rows[k]for k in common])
            assert dec==dict(rays=list(key),common_walls=common,common_rank=r,adjacent=r==8)
            if r!=8:assert key not in faces;continue
            z=faces[key];seen.add(key);facecount+=1
            xx,yy=vector(x['coordinates']),vector(y['coordinates'])
            h=[[inner(xx,g,xx),inner(xx,g,yy)],[inner(xx,g,yy),inner(yy,g,yy)]]
            det=h[0][0]*h[1][1]-h[0][1]*h[0][1]
            assert h==mat(z['plane_gram'])and h[0][0].sign()>0 and det.sign()>0
            inv=[[h[1][1]/det,-h[0][1]/det],[-h[0][1]/det,h[0][0]/det]]
            contacts=[[-Q(*rr['products'][k])for rr in(x,y)]for k in range(len(rows))]
            assert contacts==mat(z['contact_magnitudes'])and all(q.sign()>=0 for p in contacts for q in p)
            nonzero=[i for i,p in enumerate(contacts)if any(q!=0 for q in p)]
            assert[zr['wall']for zr in z['sharp_cone_bounds']]==nonzero
            eligible=[]
            for i,rec in zip(nonzero,z['sharp_cone_bounds']):
                p=contacts[i];v=[sum((inv[k][l]*p[l]for l in range(2)),Q())for k in range(2)]
                boundary=[2*p[k]*p[k]/(h[k][k]*sc[i])for k in range(2)]
                interior=all(q.sign()>0 for q in v)
                if interior:up=2*(p[0]*v[0]+p[1]*v[1])/sc[i]
                else:up=boundary[0]if(boundary[0]-boundary[1]).sign()>=0 else boundary[1]
                finite=(up-4).sign()<0 or up==4 and not interior
                assert vector(rec['contact'])==p and vector(rec['stationary_coefficients'])==v
                assert vector(rec['boundary_bounds'])==boundary and Q(*rec['upper_square'])==up
                assert rec['upper_attained_in_interior']==interior and rec['strictly_below_four_in_interior']==finite
                if finite:eligible.append(i)
            expected=next((i,j)for i,j in combinations(eligible,2)if contacts[i][0]*contacts[j][1]-contacts[i][1]*contacts[j][0]!=0)
            assert z['anchors']==list(expected);i,j=expected
            assert len(z['candidates'])==25
            for (a,b),rec in zip(product(native,repeat=2),z['candidates']):
                tests+=1;assert vector(rec['anchor_squares'])==[a,b]
                f=Field([a*sc[i],b*sc[j]]);assert f.basis==vector(rec['field_radicands'])
                assert len(rec['coefficients'])==2
                aa,bb=(fd(v,f)for v in rec['coefficients'])
                # The invertible anchor map proves uniqueness of these
                # saved coefficients; no producer solving routine is used.
                assert add(mulq(aa,contacts[i][0]),mulq(bb,contacts[i][1]))==f.sqrt(a*sc[i])
                assert add(mulq(aa,contacts[j][0]),mulq(bb,contacts[j][1]))==f.sqrt(b*sc[j])
                if f.sign(aa)<=0 or f.sign(bb)<=0:
                    assert rec['decision']=='outside_relative_interior'and'projection_norm'not in rec
                else:
                    completed=add(aa,mulq(bb,h[0][1]/h[0][0]))
                    norm=add(mulq(f.mul(completed,completed),h[0][0]),mulq(f.mul(bb,bb),det/h[0][0]))
                    assert norm==fd(rec['projection_norm'],f)and f.sign(add(norm,f.scalar(-2)))>0
                    assert rec['decision']=='physical_norm_exceeds_two'
                hist[rec['decision']]+=1
        assert seen==set(faces)
    assert(paircount,facecount,tests)==(d['complete_pair_domain'],d['adjacent_spacelike_faces'],d['native_anchor_candidates'])==(1488,336,8400)
    assert dict(hist)==d['decision_counts']==dict(outside_relative_interior=4864,physical_norm_exceeds_two=3536)
    print('PASS_STANDALONE_RAW23_TWO_RAY_SAVED',dict(pairs=paircount,faces=facecount,anchor_pairs=tests,decisions=dict(hist)),flush=True)

if __name__=='__main__':main()

