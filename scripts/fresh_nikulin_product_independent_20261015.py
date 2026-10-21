#!/usr/bin/env python3
"""Independent affine constants and scalar review; no producer imports."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_product_independent_20261015'
PRIMARY='results/fresh_nikulin_narrow_20261015_products.json'

def tree(lengths):
    n=sum(lengths)+1;g=2*s.eye(n);k=1
    for length in lengths:
        old=0
        for _ in range(length):
            g[old,k]=g[k,old]=-1;old=k;k+=1
    return g

def exceptional(g,i,j,expected):
    # Pin the last null coefficient to1; invert a proper finite subdiagram.
    J=g[:-1,:-1].inv();d=(-J*g[:-1,-1]).col_join(s.ones(1,1))
    assert g*d==s.zeros(g.rows,1)and all(x>0 for x in d)
    assert all(J[:k,:k].det()>0 for k in range(1,J.rows+1))
    b=s.zeros(g.rows,1);b[i]=1/d[i];b[j]=-1/d[j]
    assert (b.T*d)[0]==0
    z=(J*b[:-1,:]).col_join(s.zeros(1,1))
    assert g*z==b
    e=s.simplify((b.T*z)[0]);w=s.simplify(min(d[i]**2,d[j]**2)*e)
    assert w==expected and w>=s.Rational(2*(g.rows-1),9)
    return dict(rank=g.rows-1,omissions=[i,j],kernel=list(map(str,d)),weighted_energy=str(w))

def classical(family,r):
    e=s.eye(r)
    if family=='C':
        roots=[-s.sqrt(2)*e[0,:]]+[e[i,:]-e[i+1,:]for i in range(r-1)]+[s.sqrt(2)*e[-1,:]]
        choices=[(0,r),(1,r),(0,r-1)];step=s.S.One;scale=s.Rational(1,2)
    else:
        roots=[-e[0,:]-e[1,:]]+[e[i,:]-e[i+1,:]for i in range(r-1)]
        roots += [s.sqrt(2)*e[-1,:]if family=='B'else e[-2,:]+e[-1,:]]
        choices=[(0,r),(1,r),(0,r-1)]if family=='B'else[(0,r-1),(1,r-1),(0,r),(1,r)]
        step=s.Rational(1,2);scale=s.S.One
    a=s.Matrix.vstack(*roots);g=a*a.T
    null=g.nullspace();assert len(null)==1
    d=null[0]
    if d[0]<0:d=-d
    d=s.simplify(d/d[0]*(1/s.sqrt(2)if family=='C'else 1))
    # Difference energies can be computed from contact rows using any
    # independent Euclidean root frame; no affine pseudoinverse is used.
    J=a[1:,:].inv();rows=[]
    for left,right in choices:
        allowed=[]
        for omit in (left,right):
            keep=[i for i in range(r+1)if i!=omit];inv=g.extract(keep,keep).inv()
            assert all(x>=0 for x in inv)
            allowed.append([omit]+[i for k,i in enumerate(keep)if inv[k,k]<=2])
        # The contact-row simplices' extreme differences suffice to check
        # constant coordinates, so arbitrary convex combinations retain them.
        points=[]
        for side in allowed:
            pts=[]
            for i in side:
                b=s.zeros(r+1,1);b[i]=1/d[i]
                pts.append(J*b[1:,:])
            points.append(pts)
        fixed=[]
        for k in range(r):
            differences={s.simplify(x[k]-y[k])for x in points[0]for y in points[1]}
            if len(differences)==1 and next(iter(differences))in(step,-step):fixed.append(k)
        assert len(fixed)>=r-5
        actual_scale=min(d[i]**2 for side in allowed for i in side)
        assert actual_scale>=scale
        lower=actual_scale*len(fixed)*step**2
        assert lower>=s.Rational(r-5,2 if family=='C'else 4)
        rows.append(dict(family=family,rank=r,omissions=[left,right],allowed=allowed,constant_coordinates=fixed,lower_bound=str(lower)))
    return rows

def run():
    old=json.loads((ROOT/PRIMARY).read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert old['status']=='PASS_CONDITIONAL_PRODUCT_CUSP_BOUND_225'
    ff=2*s.eye(5)
    for i,w in enumerate([1,1,s.sqrt(2),1]):ff[i,i+1]=ff[i+1,i]=-w
    gg=s.Matrix([[2,-1,0],[-1,2,-s.sqrt(3)],[0,-s.sqrt(3),2]])
    specs=[('A1',s.Matrix([[2,-2],[-2,2]]),0,1,s.Rational(1,2)),
           ('E6',tree([2,2,2]),2,4,s.Rational(4,3)),
           ('E7',tree([1,3,3]),3,6,s.Integer(2)),
           ('E8',tree([1,2,5]),3,7,s.Integer(2)),
           ('F4',ff,0,4,s.Integer(1)),('G2',gg,0,2,s.Rational(2,3))]
    exc=[dict(type=name,**exceptional(g,i,j,w))for name,g,i,j,w in specs]
    loc=[row for family in('B','C','D')for r in(10,12,17)for row in classical(family,r)]
    # Cycle energy for disjoint arcs after removing radius-two end clusters.
    for r in range(10,201):
        n=r+1;a=n//2-4;b=n-n//2-4
        assert a>0 and b>0 and F(a*b,a+b)>=F(r-8,4)
    # Universal parity tails are identities, rather than these finite samples.
    t=s.symbols('t',nonnegative=True,integer=True)
    for length in(2*t+2,2*t+3):
        a=s.floor(length/2);b=length-a
        assert s.simplify(a*b/length-(length-1)/4)>=0
    rank=10+t;baseline=2*rank/9-2
    surpluses=[s.expand((rank-8)/4-baseline),s.expand((rank-5)/4-baseline),s.expand((rank-5)/2-baseline)]
    assert all(all(c>=0 for c in s.Poly(v,t).all_coeffs())for v in surpluses)
    assert F(2)/F(1,2)==4 and 2*4+1==9
    N=225+t
    bound=s.expand((2*N/9-18)/2-2)
    assert bound==14+t/9
    paths=[PRIMARY,'scripts/fresh_nikulin_narrow_20261015_products.py',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_CONDITIONAL_PRODUCT_CUSP_ARITHMETIC',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        exceptional_reconstructions=exc,localization_checks=loc,large_rank_surpluses=list(map(str,surpluses)),
        final_pairing_lower_bound=str(bound),first_excluded_ambient_dimension=226,
        arithmetic='Separate finite-subdiagram inversion, affine Euclidean coordinates, contact-simplex extreme points and symbolic parity/rank tails. No producer imports.',
        scope='Checks the affine constants and universal scalar estimates. Geometric contact, common-height and framing-wall arguments require the accompanying proof review. The all-finite narrow-facet alternative is not settled here; no unconditional225 claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
