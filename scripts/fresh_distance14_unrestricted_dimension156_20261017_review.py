#!/usr/bin/env python3
"""Separate Fraction-only local Gram and cusp-profile review; no producer imports."""
import argparse,json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_distance14_unrestricted_dimension156_20261017'
TAG=BASE+'_review'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
RADIUS='fresh_simple_cusp_distance_bound_20261017'

def digest(path):return sha256(path.read_bytes()).hexdigest()
def matrix(n):return [[F(2)if i==j else F()for j in range(n)]for i in range(n)]
def sub(a,ii,jj=None):return [[a[i][j]for j in(ii if jj is None else jj)]for i in ii]
def inverse(a):
    n=len(a);b=[list(a[i])+[F(i==j)for j in range(n)]for i in range(n)]
    for j in range(n):
        pivot=next(i for i in range(j,n)if b[i][j]);b[j],b[pivot]=b[pivot],b[j]
        q=b[j][j];b[j]=[v/q for v in b[j]]
        for i in range(n):
            if i!=j:
                q=b[i][j];b[i]=[v-q*w for v,w in zip(b[i],b[j])]
    return [r[n:]for r in b]
def dot(a,j,b):return sum((a[i]*j[i][k]*b[k]for i in range(len(a))for k in range(len(b))),F())
def positive_pivots(a):
    b=[list(r)for r in a];out=[]
    while b:
        p=b[0][0];assert p>0;out.append(p)
        b=[[b[i][j]-b[i][0]*b[0][j]/p for j in range(1,len(b))]for i in range(1,len(b))]
    return out

def local_checks(census,radius):
    rows=[]
    for d in census['diagrams']:
        graph=[set()for _ in range(d['order'])]
        for a,b,_ in d['edges']:graph[a].add(b);graph[b].add(a)
        for marks in combinations(d['finite_deletions'],2):
            dist={p:0 for p in marks};queue=deque(marks)
            while queue:
                a=queue.popleft()
                for b in graph[a]:
                    if b not in dist:dist[b]=dist[a]+1;queue.append(b)
            assert len(dist)==d['order'] and max(dist.values())<=7
            rows.append(dict(diagram=d['index'],black=list(marks),radius=max(dist.values()),
                farthest=sorted(i for i in dist if dist[i]==max(dist.values()))))
    assert len(rows)==384
    maximizers=[r for r in rows if r['radius']==7]
    expected=[(23,[0,1],[8]),(25,[1,2],[9]),(25,[1,3],[9]),(25,[2,3],[9])]
    assert [(r['diagram'],r['black'],r['farthest'])for r in maximizers]==expected
    assert [(r['diagram'],r['black'],r['farthest'])for r in radius['radius_maximizers']]==expected
    diagrams={d['index']:d for d in census['diagrams']}
    literal23=[[0,2,3],[1,2,3],[1,9,3],[2,3,3],[3,4,3],[4,5,3],[5,6,3],[6,7,3],[7,8,4]]
    literal25=[[0,1,3],[0,2,3],[0,4,3],[2,3,3],[4,5,3],[5,6,3],[6,7,3],[7,8,3],[8,9,3]]
    assert diagrams[23]['edges']==literal23 and diagrams[25]['edges']==literal25
    # Scale physical node8 by sqrt2: the entire L23 Gram becomes rational.
    g=matrix(10);g[8][8]=4
    for i,j,m in literal23:g[i][j]=g[j][i]=F(-2 if m==4 else -1)
    l23=[]
    for black,expected_value in((0,F(9,2)),(1,F(4))):
        ids=[i for i in range(10)if i!=black];a=sub(g,ids);pivots=positive_pivots(a)
        # Original inverse = D*(scaled inverse)*D; node8 has D88²=2.
        value=2*inverse(a)[ids.index(8)][ids.index(8)]
        assert value==expected_value and value>2
        l23.append(dict(black=black,normalized_inverse_far_diagonal=str(value),
            rational_congruence_positive_pivots=list(map(str,pivots))))
    tail=matrix(3);tail[0][1]=tail[1][0]=tail[1][2]=tail[2][1]=F(-1)
    residual=1/inverse(tail)[2][2];assert residual==F(4,3) and residual**2<2
    g=matrix(20)
    for offset in(0,10):
        for i,j,m in literal25:
            assert m==3;g[offset+i][offset+j]=g[offset+j][offset+i]=F(-1)
    g[9][19]=g[19][9]=F(-1)
    centers=[0,10];ids=[i for i in range(20)if i not in centers]
    a=sub(g,ids);pivots=positive_pivots(a);j=inverse(a);x=sub(g,centers,ids)
    schur=[[g[c][d]-dot(x[i],j,x[k])for k,d in enumerate(centers)]for i,c in enumerate(centers)]
    assert schur==[[F(-7,78),F(-1,13)],[F(-1,13),F(-7,78)]]
    eig=[schur[0][0]+schur[0][1],schur[0][0]-schur[0][1]]
    assert eig==[F(-1,6),F(-1,78)]
    return dict(recomputed_marked_graph_distances=len(rows),radius7_cases=maximizers,
        L23_rational_congruence=l23,A3_tail_residual=str(residual),
        bridge_label4_square_margin=str(2-residual**2),
        joined_L25_positive_pivots=list(map(str,pivots)),joined_L25_schur=[[str(v)for v in r]for r in schur],
        joined_L25_negative_eigenvalues=list(map(str,eig)),
        shared_input_scope='The frozen catalogue supplies completeness and finite deletions; the graph distances and literal new Gram calculations are independently recomputed here.')

def convolution(ranks):
    h=[1]
    for r in ranks:
        q=[];window=0
        for j in range(len(h)+r):
            if j<len(h):window+=h[j]
            if 0<=j-r-1<len(h):window-=h[j-r-1]
            q.append(window)
        h=q
    return h

def inclusion_coefficients(ranks):
    s=len(ranks);shifts=[]
    for mask in range(1<<s):
        ids=[i for i in range(s)if mask>>i&1]
        shifts.append((sum(ranks[i]+1 for i in ids),(-1)**len(ids)))
    return [sum(sign*comb(j-shift+s-1,s-1)for shift,sign in shifts if j>=shift)for j in range(sum(ranks)+1)]

def profile(ranks):
    d=sum(ranks);h=convolution(ranks);other=inclusion_coefficients(ranks)
    assert h==other and h==h[::-1]
    v=prod(r+1 for r in ranks);assert sum(h)==v
    abs_moment=sum((F(a)*abs(F(2*j-d,2))for j,a in enumerate(h)),F())/v
    kappa=F(3*d,4)-F(5,2)*abs_moment
    cost=sum((min(F(14),F(r,2))for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=v,normalized_curvature=str(kappa),
        normalized_cost_bound=str(cost),margin=str(kappa-cost),mean_absolute_deviation=str(abs_moment))

def extremal_profiles():
    rows=[];minima=[];finite=[]
    for d in range(156,163):
        # Enumerate the optional large factor first; all remaining ranks have
        # unique division into28-blocks and one remainder.
        shapes=[]
        for b in[0]+list(range(29,d+1)):
            q,a=divmod(d-b,28)
            shapes.append(tuple(sorted([28]*q+([a]if a else[])+([b]if b else[]))))
        assert len(shapes)==len(set(shapes))==d-27
        group=[profile(r)for r in sorted(shapes)]
        assert all(F(r['margin'])>0 for r in group)
        mm=min(F(r['margin'])for r in group)
        minima.append(dict(ambient_dimension=d+1,profiles=len(group),minimum_margin=str(mm),
            minimizing_ranks=[r['ranks']for r in group if F(r['margin'])==mm]))
        rows.extend(group)
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>14*n;finite.append(dict(dimension=n,finite_budget=14*n,beta=str(beta),strict_margin=str(beta-14*n)))
    assert len(rows)==924 and min(F(r['margin'])for r in rows)==F(1476616,20511149)
    boundary=profile((15,28,28,28,28,28));assert F(boundary['margin'])==F(-4172075,82044596)
    return rows,minima,finite,boundary

def transfer_checks():
    count=0
    # All transfers occurring in the seven target dimensions are within
    # this finite box. The same explicit formulas also give the all-rank proof.
    for a in range(1,82):
        for b in range(a,163-a):
            den0=(a+1)*(b+1);den1=a*(b+2);assert den0-den1==b-a+1>0
            total=F();moment=F()
            for j in range(a+b+1):
                old=min(j+1,a+1,b+1,a+b-j+1)
                new=min(j+1,a,b+2,a+b-j+1)
                diff=F(new,den1)-F(old,den0)
                if j<a or j>b:
                    numerator=min(j+1,a+b-j+1)
                    expected=F(numerator*(b-a+1),den0*den1);assert diff==expected>0
                else:assert diff==F(-1,(b+1)*(b+2))<0
                total+=diff;moment+=j*diff
            assert total==moment==0;count+=1
    return dict(all_pairs_with_total_at_most162=count,
        low_transfer='1<=a<=b<=27, stop at a=0 or b=28',
        high_transfer='28<a<=b, stop at a=28',
        analytic_sign_pattern='Positive outside j=a..b; negative inside; symmetric, mass0 and mean0.',
        convex_order_argument='Choose the chord of an arbitrary convex test function through centered points ±(b−a+1)/2. Its difference from the chord is nonpositive on the negative-mass interval and nonnegative on the positive-mass tails. The chord contributes zero by mass and mean, proving convex order. Independent convolution preserves it.')

def run():
    primary=json.loads((ROOT/'results'/f'{BASE}.json').read_text())
    assert primary['status']=='PASS_DISTANCE14_AND_CONCENTRATED_CUSP_BUDGET_DIMENSION156'
    assert primary['dimension_bound']==156 and primary['unrestricted_finite_volume']
    for name,expected in primary['sources'].items():assert digest(ROOT/name)==expected,name
    census=json.loads((ROOT/'results'/f'{CENSUS}.json').read_text())
    radius=json.loads((ROOT/'results'/f'{RADIUS}.json').read_text())
    local=local_checks(census,radius);transfer=transfer_checks();rows,minima,finite,boundary=extremal_profiles()
    byprofile={tuple(sorted(r['ranks'])):r for r in primary['reduced_profiles']}
    assert len(byprofile)==len(rows)==924
    for row in rows:
        old=byprofile[tuple(row['ranks'])]
        for k in('link_dimension','link_vertices','normalized_curvature','normalized_cost_bound','margin','mean_absolute_deviation'):
            assert row[k]==old[k],(row['ranks'],k,row[k],old[k])
    assert local['joined_L25_schur']==primary['exact_local_obstructions']['joined_L25_schur']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{BASE}.{ext}'for folder,ext in[('scripts','py'),('results','json'),('results','md')]]
    paths +=[ROOT/'results'/f'{s}.{ext}'for s in(CENSUS,RADIUS)for ext in('json','md')]
    paths +=[ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md',ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md',
        ROOT/'results/dimension163_proof_status_20261017.md']
    return dict(status='PASS_INDEPENDENT_DISTANCE14_UNRESTRICTED_DIMENSION156_REVIEW',
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},dimension_bound=156,unrestricted_finite_volume=True,
        excluded_dimensions=list(range(157,164)),inherited_upper_tail_bound=163,local_gram_and_radius_review=local,
        convex_transfer_review=transfer,reduced_profiles=rows,dimension_minima=minima,finite_budgets=finite,
        dimension156_payment_limitation=boundary,
        scope='Separate Fraction-only arithmetic, literal graph distances, probability identities and all924 extremal profiles, with full written geometric/composition review. No producer arithmetic imports, no SymPy, no repeated high-rank classification census. Complete internally reviewed computer-assisted proof; not proof-assistant or external-referee validation.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['reduced_profiles']),'profiles',flush=True)
