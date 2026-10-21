#!/usr/bin/env python3
"""Distance14 obstruction and exact concentrated cusp-link budgets."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb, prod
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance14_unrestricted_dimension156_20261017'
INPUTS=(
 'fresh_high_rank_quasilanner_census_20261017',
 'fresh_simple_cusp_distance_bound_20261017',
 'fresh_cusp_factor_payment_framework_20261017',
 'fresh_cusp_link_shape_curvature_20261017',
 'fresh_nikulin_two_budget_frontier_20261015_dimension163',
)


def digest(path):
    h=sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def gram(diagram):
    g=2*sp.eye(diagram['order'])
    for a,b,m in diagram['edges']:
        assert m in (3,4)
        g[a,b]=g[b,a]=-1 if m==3 else -sp.sqrt(2)
    return g


def literal_obstructions(census,radius):
    rows=radius['black_pair_checks']
    assert len(rows)==384 and max(r['radius']for r in rows)==7
    maximum=[r for r in rows if r['radius']==7]
    assert [(r['diagram'],r['black'],r['farthest'])for r in maximum]==[
        (23,[0,1],[8]),(25,[1,2],[9]),(25,[1,3],[9]),(25,[2,3],[9])]
    diagrams={d['index']:d for d in census['diagrams']}
    g=gram(diagrams[23]);inverses=[]
    for black,expected in[(0,sp.Rational(9,2)),(1,sp.Integer(4))]:
        ids=[i for i in range(10)if i!=black]
        val=sp.simplify(g.extract(ids,ids).inv()[ids.index(8),ids.index(8)])
        assert val==expected and val>2
        inverses.append(dict(deleted_black=black,far_node=8,inverse_diagonal=str(val)))
    # A3 at either end of a proposed cross-edge leaves residual4/3.
    a3=sp.Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]])
    assert 1/a3.inv()[2,2]==sp.Rational(4,3)
    assert sp.Rational(4,3)**2<2
    # At label3, two complete L25 diagrams joined at their far nodes.
    g=sp.diag(gram(diagrams[25]),gram(diagrams[25]))
    g[9,19]=g[19,9]=-1
    centers=[0,10];other=[i for i in range(20)if i not in centers]
    a=g.extract(other,other);x=g.extract(centers,other)
    lower,diagonal=a.LDLdecomposition()
    assert lower*diagonal*lower.T==a and all(diagonal[i,i]>0 for i in range(18))
    schur=g.extract(centers,centers)-x*a.inv()*x.T
    assert schur==sp.Matrix([[-sp.Rational(7,78),-sp.Rational(1,13)],[-sp.Rational(1,13),-sp.Rational(7,78)]])
    assert sorted(schur.eigenvals())==[-sp.Rational(1,6),-sp.Rational(1,78)]
    return dict(radius7_contexts=maximum,L23_inverse_obstructions=inverses,
                bridge_tail_residual='4/3',bridge_label4_squared_margin='2-16/9=2/9>0',
                joined_L25_positive_block_order=18,
                joined_L25_positive_pivots=[str(diagonal[i,i])for i in range(18)],
                joined_L25_schur=[[str(schur[i,j])for j in range(2)]for i in range(2)],
                joined_L25_negative_eigenvalues=['-1/6','-1/78'])


@lru_cache(None)
def hpoly(ranks):
    if not ranks:return (1,)
    p=hpoly(ranks[:-1]);r=ranks[-1];q=[0]*(len(p)+r)
    for i,a in enumerate(p):
        for j in range(r+1):q[i+j]+=a
    return tuple(q)


def independent_absolute_moment(ranks):
    """Truncated bounded-composition sum; no polynomial convolution."""
    d=sum(ranks);s=len(ranks);cut=(d-1)//2;total=0
    for size in range(s+1):
        for subset in combinations(range(s),size):
            shift=sum(ranks[i]+1 for i in subset);m=cut-shift
            if m<0:continue
            term=(d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1)
            total+=(-1)**size*term
    return F(total,prod(r+1 for r in ranks))


def reduced_shapes(d):
    q,a=divmod(d,28)
    yield (28,)*q+((a,)if a else())
    for q in range(d//28+1):
        for a in range(28):
            b=d-28*q-a
            if b>=29:yield (28,)*q+((a,)if a else())+(b,)


def profile(ranks):
    assert all(r>0 for r in ranks)
    d=sum(ranks);h=hpoly(ranks);v=prod(r+1 for r in ranks)
    assert h==h[::-1] and sum(h)==v
    g=sum(a*min(j,d-j)for j,a in enumerate(h))
    curvature=F(5*g-d*v,2*v)
    moment=independent_absolute_moment(ranks)
    assert curvature==F(3*d,4)-F(5,2)*moment
    cost=sum((min(F(14),F(r,2))for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=v,G=g,
                mean_absolute_deviation=str(moment),normalized_curvature=str(curvature),
                normalized_cost_bound=str(cost),margin=str(curvature-cost))


def transfer_checks():
    a,b=sp.symbols('a b',integer=True,positive=True)
    assert sp.expand((a+1)*(b+1)-a*(b+2))==b-a+1
    assert sp.cancel(1/(b+2)-1/(b+1)+1/((b+1)*(b+2)))==0
    probes=[]
    for a in range(1,9):
        for b in range(a,13):
            old=hpoly((a,b));new=hpoly((a-1,b+1));d=a+b
            den0=(a+1)*(b+1);den1=a*(b+2)
            differences=[F(new[j],den1)-F(old[j],den0)for j in range(d+1)]
            assert sum(differences,F())==0 and differences==differences[::-1]
            assert all((z>0 if j<a or j>b else z<0)for j,z in enumerate(differences))
            probes.append([a,b])
    return dict(denominator_difference='b-a+1>0 for1<=a<=b',
                central_mass_difference='-1/((b+1)*(b+2))',
                bounded_sign_probes=probes,
                scope='Exact polynomial identities and the written one-crossing argument prove convex order at every rank; probes only corroborate the distribution convention.')


def run():
    data={name:json.loads((ROOT/f'results/{name}.json').read_text())for name in INPUTS}
    for name,d in data.items():
        assert d['status'].startswith('PASS'),name
        for source,expected in d.get('sources',{}).items():
            assert isinstance(expected,str) and digest(ROOT/source)==expected,(name,source)
    obstructions=literal_obstructions(data[INPUTS[0]],data[INPUTS[1]])
    transfer=transfer_checks()
    records=[];minima=[]
    for d in range(156,163):
        shapes=list(reduced_shapes(d))
        assert len(shapes)==len(set(tuple(sorted(r))for r in shapes))
        group=[profile(r)for r in shapes]
        assert all(F(r['margin'])>0 for r in group)
        minimum=min(F(r['margin'])for r in group)
        minima.append(dict(ambient_dimension=d+1,profiles=len(group),minimum_margin=str(minimum),
                           minimizing_ranks=[r['ranks']for r in group if F(r['margin'])==minimum]))
        records+=group
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>14*n
    assert min(F(r['margin'])for r in records)==F(1476616,20511149)
    boundary=profile((28,)*5+(15,))
    assert F(boundary['margin'])==F(-4172075,82044596)<0
    paths=[Path(__file__)]+[ROOT/f'{folder}/{name}.{extension}'
        for name in INPUTS for folder,extension in [('scripts','py'),('results','json'),('results','md')]]
    paths += [ROOT/'results/fresh_simple_cusp_distance_bound_20261017_review.md',
              ROOT/'results/dimension163_proof_status_20261017.md']
    return dict(status='PASS_DISTANCE14_AND_CONCENTRATED_CUSP_BUDGET_DIMENSION156',
                sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
                dimension_bound=156,unrestricted_finite_volume=True,arithmetic_hypothesis=False,
                simple_cusp_hypothesis=False,factor_rank_or_count_restriction=None,
                inherited_unrestricted_dimension_bound=163,excluded_ambient_dimensions=list(range(157,164)),
                triangle_distance_bound=8,quadrilateral_distance_bound=14,
                finite_vertex_cost='14n',ideal_normalized_cost='sum_i min(14,r_i/2)',
                exact_local_obstructions=obstructions,convex_order_transfer=transfer,
                reduced_profiles=records,dimension_minima=minima,
                minimum_profile_margin='1476616/20511149',distance14_budget_obstruction_at_n156=boundary,
                scope='Written unrestricted geometric payment and convex-order reduction plus exact physical Gram and reduced-profile certificates. The inherited163 theorem covers the upper tail. Not a new example, formal verification or external refereeing; no claim of a155 bound from this payment comparison.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['reduced_profiles']),'concentrated profiles',flush=True)
