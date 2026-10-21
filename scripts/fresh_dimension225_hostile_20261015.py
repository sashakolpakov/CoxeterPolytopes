#!/usr/bin/env python3
"""Independent exact checks accompanying the dimension-225 hostile proof audit.

No producer arithmetic is imported. The universal geometric proof and its
published classification inputs are reviewed in the matching markdown file.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_dimension225_hostile_20261015'


def encode(x):
    if isinstance(x, S.MatrixBase):
        return [[str(t) for t in row] for row in x.tolist()]
    return str(x)


def affine_checks():
    def star(arms):
        g = 2*S.eye(1+sum(arms)); j = 1
        for length in arms:
            prev = 0
            for _ in range(length):
                g[prev,j] = g[j,prev] = -1
                prev = j; j += 1
        return g
    def path(labels):
        g = 2*S.eye(len(labels)+1)
        for j,m in enumerate(labels):
            g[j,j+1] = g[j+1,j] = -2*S.cos(S.pi/m)
        return g
    cases = [
        ('A1',S.Matrix([[2,-2],[-2,2]]),0,1,S.Rational(1,2)),
        ('E6',star((2,2,2)),2,4,S.Rational(4,3)),
        ('E7',star((1,3,3)),3,6,S.Integer(2)),
        ('E8',star((1,2,5)),3,7,S.Integer(2)),
        ('F4',path((3,3,4,3)),0,4,S.Integer(1)),
        ('G2',path((3,6)),0,2,S.Rational(2,3)),
    ]
    out=[]
    for name,g,i,j,expected in cases:
        d=g.nullspace()[0]
        if d[0]<0: d=-d
        assert all(t>0 for t in d) and g.rank()==g.rows-1
        f=S.zeros(g.rows,1); f[i]=1/d[i];f[j]=-1/d[j]
        assert (d.T*f)[0]==0
        # Ground one potential coordinate, rather than the producer's
        # rank-one regularized pseudoinverse construction.
        keep=[k for k in range(g.rows) if k!=j]
        a=g.extract(keep,keep)
        assert all(a[:k,:k].det()>0 for k in range(1,a.rows+1))
        z=a.inv()*f.extract(keep,[0])
        energy=S.simplify((f.extract(keep,[0]).T*z)[0])
        cost=S.simplify(min(d[i]**2,d[j]**2)*energy)
        assert cost==expected and cost>=S.Rational(2*(g.rows-1),9)
        out.append(dict(type=name,rank=g.rows-1,kernel=encode(d),
                        omissions=[i,j],grounded_energy=str(energy),
                        scale_invariant_cost=str(cost)))
    return out


def tetrahedron():
    # A compact triangular facet whose full ambient star has no Lanner
    # subdiagram. Its opposite ideal vertex is an affine A2 diagram.
    g=S.Matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,-1],[0,-1,-1,2]])
    gi=g.inv(); h=-gi*S.ones(4,1)
    assert g.det()==-3 and (h.T*g*h)[0]==-22
    assert g[:3,:3].det()>0  # plus determinant<0 gives inertia(3,1).
    vertices=[-gi[:,i] for i in range(4)]
    for i,v in enumerate(vertices):
        assert g*v==-S.eye(4)[:,i]
        assert (v.T*g*h)[0]<0
        assert (v.T*g*v)[0] == gi[i,i]
    assert [gi[i,i] for i in range(4)]==[0,-2,S.Rational(-4,3),S.Rational(-4,3)]
    # All proper diagrams PSD; one is affine, so the unique indefinite
    # subdiagram (the whole graph) is not Lanner.
    proper=[]
    for k in range(1,4):
        for ids in combinations(range(4),k):
            a=g.extract(ids,ids)
            minors=[a.extract(ii,ii).det() for kk in range(1,k+1)
                    for ii in combinations(range(k),kk)]
            assert all(t>=0 for t in minors)
            proper.append([list(ids),str(a.det())])
    assert g.extract([1,2,3],[1,2,3])*S.ones(3,1)==S.zeros(3,1)
    # Facet0 has exactly the three finite vertices1,2,3. Its relative
    # interior witness is their sum and is strict against the black walls.
    f=sum(vertices[1:],S.zeros(4,1))
    assert g*f==S.Matrix([0,-1,-1,-1]) and (f.T*g*f)[0]<0
    assert all((vertices[i].T*g*vertices[j])[0]<0 for i,j in combinations(range(1,4),2))
    return dict(gram=encode(g),inverse=encode(gi),determinant='-3',
                strict_center=encode(h),center_square='-22',
                vertex_coordinates=[encode(v) for v in vertices],
                vertex_squares=[str(gi[i,i]) for i in range(4)],
                compact_facet=0,compact_facet_vertices=[1,2,3],
                affine_black_nodes=[1,2,3],proper_subdiagram_determinants=proper,
                lanner_subdiagram_count=0,
                scope='Counterexample to no-affine/Lanner extraction localization, not to the face-weight inequality.')


def scalar_identities():
    s,t,k,l=S.symbols('s t k l',positive=True)
    a,b,c,d=S.symbols('a b c d',real=True)
    x=S.Matrix([(a*a+b*b-2)/(2*s),s,a,b])
    y=S.Matrix([(c*c+d*d-2)/(2*t),t,c,d])
    g=S.Matrix([[0,-1,0,0],[-1,0,0,0],[0,0,1,0],[0,0,0,1]])
    energy=k*k*s*t*(a/(k*s)-c/(k*t))**2+l*l*s*t*(b/(l*s)-d/(l*t))**2
    assert S.cancel(-(x.T*g*y)[0]-(energy/2-s/t-t/s))==0
    assert S.cancel((x.T*g*x)[0]-2)==0
    # All-rank cycle-flow and loss identities; k>=5 in the odd case.
    u=S.symbols('u',nonnegative=True); z=u+5
    even_flow=S.cancel((z-4)/2-(2*z-9)/4)
    odd_flow=S.cancel((z-4)*(z-3)/(2*z-7)-(2*z-8)/4)
    assert even_flow==S.Rational(1,4)
    assert S.cancel(odd_flow-(u+1)/(4*u+6))==0
    r=u+10; baseline=2*r/9-2
    losses=[S.expand((r-8)/4-baseline),S.expand((r-5)/4-baseline),S.expand((r-5)/2-baseline)]
    assert all(all(c>=0 for c in S.Poly(p,u).all_coeffs()) for p in losses)
    assert F(2)/F(1,2)==4 and 2*4+1==9
    assert S.expand((225+u)/9-11)==14+u/9
    # Compact-facet arithmetic is valid only after the repaired face-weight
    # proof. Here ambient n=d+1, so the vertex bound is17(d+1), not17d.
    odd=143+2*u; even=144+2*u
    oddgap=S.expand((odd-1)*(odd-5)/8-17*(odd+1))
    evengap=S.expand(even*(even-6)/8-17*(even+1))
    assert oddgap==u*u/2+36*u+S.Rational(3,2)
    assert evengap==u*u/2+S.Rational(73,2)*u+19
    return dict(independent_scale_identity=True,cycle_even_surplus=str(even_flow),
                cycle_odd_surplus=str(odd_flow),marked_factor_surpluses=[str(p) for p in losses],
                contact_budget=4,marked_factor_budget=9,first_cusp_excluded_dimension=226,
                compact_facet_odd_gap=str(oddgap),compact_facet_even_gap=str(evengap),
                first_excluded_facet_dimension=143,compact_facet_ambient_bound=143)


def distance_counts():
    out=[]
    for r in (36,37,60):
        for family in ('A','D'):
            edges=[(i,i+1) for i in range(r-1)] if family=='A' else [(0,2),(1,2)]+[(i,i+1) for i in range(2,r-1)]
            adj=[[] for _ in range(r)]
            for i,j in edges:adj[i].append(j);adj[j].append(i)
            counts=[0]*18
            for i in range(r):
                ds=[None]*r;ds[i]=0;q=deque([i])
                while q:
                    v=q.popleft()
                    for w in adj[v]:
                        if ds[w] is None:ds[w]=ds[v]+1;q.append(w)
                for j in range(i+1,r):
                    if ds[j]<=17:counts[ds[j]]+=1
            expected=17*r-153+(family=='D')
            assert sum(counts)==expected and sum(counts)<=17*r
            out.append(dict(type=family,rank=r,counts=counts[1:],total=sum(counts)))
    assert all(F(r*(r-1),2)<=17*r for r in range(1,36))
    return out


def run():
    dependencies=[
        'scripts/fresh_nikulin_narrow_20261015.py',
        'results/fresh_nikulin_narrow_20261015.json',
        'scripts/fresh_nikulin_narrow_20261015_products.py',
        'results/fresh_nikulin_narrow_20261015_products.json',
        'scripts/fresh_nikulin_narrow_20261015_compact_facet.py',
        'results/fresh_nikulin_narrow_20261015_compact_facet.json',
    ]
    product=json.loads((ROOT/'results/fresh_nikulin_narrow_20261015_products.json').read_text())
    assert product['status']=='PASS_CONDITIONAL_PRODUCT_CUSP_BOUND_225'
    for p,h in product['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    sources=[f'scripts/{TAG}.py']+dependencies
    return dict(status='PASS_HOSTILE_AUDIT_WITH_COARSE_COMPACT_FACET_REPAIR',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                affine_exception_checks=affine_checks(),
                scalar_identities=scalar_identities(),
                distance_count_checks=distance_counts(),
                localization_counterexample=tetrahedron(),
                published_inputs=[
                    dict(author='Nikulin',location='ICM1986 printed654-655 TheoremA',
                         url='https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf'),
                    dict(author='Prokhorov',location='printed407-408 Propositions1-3',
                         url='https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng'),
                    dict(author='Vinberg',location='printed70-71 Prop1;76-77 Props7-9;81-88 fine-weight proof',
                         url='https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt')],
                rejected_step='Global no-parabolic/Lanner extraction does not follow from compactness of a single face.',
                repaired_step='Distance<=17 weights have face sum>=5-k by minimal-indefinite size<=10 and compact-face geometry.',
                scope='Independent source/proof audit and exact small checks; no independent classification of all Coxeter diagrams, no novelty/current-best claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();dest=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
