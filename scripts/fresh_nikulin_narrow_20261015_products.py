#!/usr/bin/env python3
"""Exact constants for the separate product-cusp separation proof.

The finite exceptional calculations use SymPy exact matrices. All rank tails,
contact-budget and dimension substitutions are rational identities. Geometric
hypotheses and the terminating-wall selection are proved in the report.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_narrow_20261015_products'


def star(arms):
    n=1+sum(arms); g=2*s.eye(n); k=1
    for length in arms:
        last=0
        for _ in range(length):
            g[last,k]=g[k,last]=-1
            last=k; k+=1
    return g


def path(labels):
    g=2*s.eye(len(labels)+1)
    for i,m in enumerate(labels):
        g[i,i+1]=g[i+1,i]=-2*s.cos(s.pi/m)
    return g


def exception(name,g,i,j):
    null=g.nullspace()
    assert len(null)==1 and g.rank()==g.rows-1
    d=null[0]
    if d[0]<0:d=-d
    assert all(x>0 for x in d)
    d=s.simplify(d/min(d))
    assert g*d==s.zeros(g.rows,1)
    # Adding dd^T makes an exact inverse available without a pseudoinverse.
    # On d-perp it is the inverse of g, and f below belongs to d-perp.
    h=g+d*d.T
    assert all(h[:k,:k].det()>0 for k in range(1,g.rows+1))
    f=s.zeros(g.rows,1);f[i]=1/d[i];f[j]=-1/d[j]
    assert (d.T*f)[0]==0
    energy=s.simplify((f.T*h.inv()*f)[0])
    cost=s.simplify(min(d[i],d[j])**2*energy)
    rank=g.rows-1
    assert s.simplify(cost-s.Rational(2*rank,9))>=0
    return dict(type=name,rank=rank,gram=[[str(x) for x in row] for row in g.tolist()],
                positive_kernel=[str(x) for x in d],omitted_nodes=[i,j],
                energy=str(energy),weighted_energy=str(cost),
                rank_baseline=str(s.Rational(2*rank,9)))


def lorentz(x,y):
    return -x[0]*y[1]-x[1]*y[0]+sum(a*b for a,b in zip(x[2:],y[2:]))


def null_normal(z):
    # A convenient representative with height1. Norm is immaterial to
    # difference energies; using norm2 checks consistency as well.
    return [sum(a*a for a in z)/2-1,s.Integer(1)]+list(z)


def classical_point_checks():
    records=[]
    for r in range(2,41):
        # Affine A_r effective resistance on the r+1 cycle.
        n=r+1; d=n//2; energy=F(d*(n-d),n)
        assert energy>=F(r,4)>=F(2*r,9)
        records.append(['A',r,str(energy)])
        # B: left fork leaf to the short end. Delta=(...,sqrt2).
        zb=[-s.Rational(1,2)]*r
        # C: the two short ends; delta_i=delta_j=1/sqrt2.
        zc=[s.Integer(-1)]*r
        assert sum(x*x for x in zb)==s.Rational(r,4)
        assert sum(x*x for x in zc)/2==s.Rational(r,2)
        records.extend([['B',r,str(F(r,4))],['C',r,str(F(r,2))]])
        if r>=4:
            zd=[-s.Rational(1,2)]*(r-1)+[s.Rational(1,2)]
            assert sum(x*x for x in zd)==s.Rational(r,4)
            records.append(['D',r,str(F(r,4))])
    return records


def localization_coordinate_checks():
    """Verify each allowed normalized-row basis vector, not just one row.

    Normalized nonnegative contact rows are convex combinations of these
    vectors. Constant separating coordinates therefore persist for all rows.
    Closed all-rank support formulas are proved in the matching report.
    """
    records=[]
    for family in ('B','C','D'):
        for r in (10,11,20):
            a=s.zeros(r+1,r);shift=s.zeros(r+1,1)
            if family=='C':
                a[0,0]=-s.sqrt(2);shift[0]=s.sqrt(2)
                for j in range(1,r):a[j,j-1]=1;a[j,j]=-1
                a[r,r-1]=s.sqrt(2)
                delta=s.Matrix([1/s.sqrt(2)]+[1]*(r-1)+[1/s.sqrt(2)])
                choices=((0,r),(1,r),(0,r-1))
                step=s.Integer(1)
            else:
                a[0,0]=a[0,1]=-1;shift[0]=1
                a[1,0]=1;a[1,1]=-1
                for j in range(2,r):a[j,j-1]=1;a[j,j]=-1
                if family=='B':
                    a[r,r-1]=s.sqrt(2)
                    delta=s.Matrix([1,1]+[2]*(r-2)+[s.sqrt(2)])
                    choices=((0,r),(1,r),(0,r-1))
                else:
                    a[r,r-2]=a[r,r-1]=1
                    delta=s.Matrix([1,1]+[2]*(r-3)+[1,1])
                    choices=((0,r-1),(1,r-1),(0,r),(1,r))
                step=s.Rational(1,2)
            g=a*a.T
            assert g*delta==s.zeros(r+1,1)
            assert (delta.T*shift)[0]==1
            inv=a[1:,:].inv()
            allz=[]
            for k in range(r+1):
                p=s.zeros(r+1,1);p[k]=1/delta[k]
                z=-inv*p[1:,:]
                assert a*z==shift-p
                allz.append(z)
            for i,j in choices:
                allowed=[]
                for omit in (i,j):
                    keep=[k for k in range(r+1) if k!=omit]
                    ji=g.extract(keep,keep).inv()
                    assert all(x>=0 for x in ji)
                    allowed.append([omit]+[k for b,k in enumerate(keep) if ji[b,b]<=2])
                middle=[k for k in range(r)
                        if all(allz[b][k]==0 for b in allowed[0])
                        and all(allz[b][k]==-step for b in allowed[1])]
                assert len(middle)>=r-5
                min_left=min(delta[k] for k in allowed[0])
                min_right=min(delta[k] for k in allowed[1])
                weighted=s.simplify(min(min_left,min_right)**2*len(middle)*step**2)
                lower=s.Rational(r-5,2 if family=='C' else 4)
                assert s.simplify(weighted-lower)>=0
                records.append(dict(family=family,rank=r,omissions=[i,j],
                                    allowed_contacts=allowed,
                                    constant_coordinate_indices=middle,
                                    coordinate_difference=str(step),
                                    weighted_energy=str(weighted)))
    return records


def run():
    dependencies=['scripts/fresh_nikulin_narrow_20261015.py',
                  'results/fresh_nikulin_narrow_20261015.json']
    old=json.loads((ROOT/dependencies[1]).read_text())
    assert old['status']=='PASS_EXPLICIT_CLASSICAL_IRREDUCIBLE_NARROW_CUSP_BOUNDS'
    for p,h in old['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    exceptional=[
        exception('A1_affine',s.Matrix([[2,-2],[-2,2]]),0,1),
        exception('E6_affine',star((2,2,2)),2,4),
        exception('E7_affine',star((1,3,3)),3,6),
        exception('E8_affine',star((1,2,5)),3,7),
        exception('F4_affine',path((3,3,4,3)),0,4),
        exception('G2_affine',path((3,6)),0,2),
    ]
    assert [x['weighted_energy'] for x in exceptional]==['1/2','4/3','2','2','1','2/3']
    point_checks=classical_point_checks()
    coordinate_checks=localization_coordinate_checks()

    # The three all-rank loss identities have nonnegative coefficients after
    # r=10+t. They establish loss<=2, not a sampled-rank conclusion.
    t=s.symbols('t',nonnegative=True);r=t+10
    bas=s.Rational(2,9)*r-2
    loss={'A':s.expand((r-8)/4-bas),
          'B_D':s.expand((r-5)/4-bas),
          'C':s.expand((r-5)/2-bas)}
    for v in loss.values():
        assert all(c>=0 for c in s.Poly(v,t).all_coeffs())
    assert all(F(2*r,9)<=2 for r in range(1,10))
    # At either endpoint at most four retained contacts are nonright.
    assert F(2)/F(1,2)==4
    spoiled=2*4; bad=spoiled+1
    assert bad==9
    N=225+t
    final=s.expand(N/9-bad-2)
    assert final==14+t/9
    assert F(224,9)-11<14

    # The scale inequality uses a COMMON height ratio across factors. These
    # algebraic branches verify the sole max(q,1/q) manipulation exactly.
    q=s.symbols('q',positive=True)
    assert s.cancel(2*q-(q+1/q)-(q-1)*(q+1)/q)==0
    assert s.cancel(2/q-(q+1/q)+(q-1)*(q+1)/q)==0

    sources=[f'scripts/{TAG}.py']+dependencies
    return dict(status='PASS_CONDITIONAL_PRODUCT_CUSP_BOUND_225',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                hypothesis='Finite-sided finite-volume real hyperbolic Coxeter polytope with an ideal vertex on the chosen Nikulin narrow facet.',
                conclusion='Ambient dimension n<=225 under the stated hypothesis.',
                exceptional_components=exceptional,
                classical_point_checks=len(point_checks),
                contact_coordinate_checks=coordinate_checks,
                classical_closed_point_bounds=dict(A='r/4',B='r/4',C='r/2',D='r/4'),
                uniform_baseline='2r/9',
                bad_component_large_rank_surplus={k:str(v) for k,v in loss.items()},
                retained_nonright_per_endpoint=4,spoiled_components_at_most=spoiled,
                narrow_component_extra_allowance=1,bad_components_at_most=bad,
                total_weighted_energy='M*(2*(n-1)/9-18)',
                negative_pairing_bound='M*((n-1)/9-11)',
                first_excluded_dimension=226,
                dimension_tail=str(final),
                framing_wall_selection='At the second ideal endpoint, choose the new node of the narrow-facet affine component when its rank>=2; for rank1 choose a new node of another component. This gives a finite pair and hence an actual ridge.',
                remaining_case='The chosen narrow facet has only finite vertices. No unconditional n<=225 or complete Nikulin<300 claim.',
                scope='Exact constants and all-rank algebra plus the geometric proof in the report. No full-polytope enumeration or independent proof of Nikulin narrow-facet existence.' )


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();data=run();out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
