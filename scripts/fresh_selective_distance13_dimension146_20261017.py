#!/usr/bin/env python3
"""Selective distance13: isolate the shared-E10 exception, then exact budgets."""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb, prod
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_selective_distance13_dimension146_20261017'
INPUTS=(
 'fresh_high_rank_quasilanner_census_20261017',
 'fresh_simple_cusp_distance_bound_20261017',
 'fresh_cusp_factor_payment_framework_20261017',
 'fresh_cusp_link_shape_curvature_20261017',
 'fresh_distance14_unrestricted_dimension156_20261017',
 'fresh_distance14_unrestricted_dimension156_20261017_review',
)

def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def serial_matrix(m):return [[str(x)for x in m.row(i)]for i in range(m.rows)]

def gram(d):
    g=2*sp.eye(d['order'])
    for a,b,m in d['edges']:
        assert m in(3,4)
        g[a,b]=g[b,a]=-1 if m==3 else -sp.sqrt(2)
    return g

def distances(g,marks):
    ds=[None]*g.rows;q=deque(marks)
    for i in marks:ds[i]=0
    while q:
        i=q.popleft()
        for j in range(g.rows):
            if j!=i and g[i,j]!=0 and ds[j]is None:
                ds[j]=ds[i]+1;q.append(j)
    assert all(d is not None for d in ds)
    return ds

def high_distance_ports(census,old):
    diagrams={d['index']:d for d in census['diagrams']};records=[];count=0
    for row in old['black_pair_checks']:
        g=gram(diagrams[row['diagram']]);ds=distances(g,row['black'])
        assert ds==row['distances'];count+=1
        if max(ds)<6:continue
        ports=[]
        for v,d in enumerate(ds):
            if d<6:continue
            vals=[]
            for b in row['black']:
                ids=[i for i in range(g.rows)if i!=b]
                h=g.extract(ids,ids);diagonal=h.LDLdecomposition()[1]
                assert all(diagonal[i,i]>0 for i in range(h.rows))
                vals.append(sp.simplify(h.inv()[ids.index(v),ids.index(v)]))
            ports.append(dict(node=v,distance=d,inverse_diagonals=list(map(str,vals)),
                              may_contact_external_white=bool(max(vals)<2)))
        records.append(dict(diagram=row['diagram'],black=row['black'],radius=max(ds),ports=ports))
    assert count==384 and len(records)==16
    allowed={(r['diagram'],tuple(r['black']),p['node'])for r in records for p in r['ports']if p['may_contact_external_white']}
    assert allowed=={(25,b,9)for b in[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]}|{(25,(1,2),8)}
    return records

def shared_core(d):
    g=2*sp.eye(19)
    for mp in[list(range(10)),[i+10 if i<9 else 9 for i in range(10)]]:
        for a,b,m in d['edges']:
            assert m==3;g[mp[a],mp[b]]=g[mp[b],mp[a]]=-1
    return g

def local_obstructions(d):
    # Any surviving disjoint far-far join is the literal old distance14 obstruction.
    g=sp.diag(gram(d),gram(d));g[9,19]=g[19,9]=-1
    ids=[i for i in range(20)if i not in(0,10)];a=g.extract(ids,ids)
    ld=a.LDLdecomposition()[1];assert all(ld[i,i]>0 for i in range(18))
    c=g.extract([0,10],ids);s=g.extract([0,10],[0,10])-c*a.inv()*c.T
    assert s==sp.Matrix([[-sp.Rational(7,78),-sp.Rational(1,13)],[-sp.Rational(1,13),-sp.Rational(7,78)]])
    # Far-to-near join: its WHITE T(1,2,5) already fails positive definiteness.
    h=sp.diag(gram(d),gram(d));h[9,18]=h[18,9]=-1
    white=[18,19,17,16,9,8,7,6,5]
    a=h.extract(white,white);ks=a.nullspace();assert len(ks)==1
    k=ks[0];k=k if all(x>0 for x in k)else-k
    assert all(x>0 for x in k)and a*k==sp.zeros(9,1)
    # PSD at the minimal simple cross bond; increasing or adding bonds is non-PD.
    ids=list(range(8));ld=a.extract(ids,ids).LDLdecomposition()[1]
    assert all(ld[i,i]>0 for i in range(8))
    return dict(far_far_schur=serial_matrix(s),far_far_negative_eigenvalues=['-1/6','-1/78'],
                far_near_white_affine_E8_nodes=white,far_near_positive_kernel=list(map(str,k)),
                cross_bond_squared_upper_bound='16/9<2',
                arbitrary_bond_scope='White A3 tails force the far-far bond to be simple. The far-near positive-kernel obstruction covers all magnitudes >=1 and extra nonpositive contacts.')

def isolation_cases(g):
    z=sp.Matrix([-6,-3,-4,-2,-5,-4,-3,-2,-1,0,6,3,4,2,5,4,3,2,1])
    assert g*z==sp.zeros(19,1)and len(g.nullspace())==1
    contexts=[([1,2],2),([1,3],1),([2,3],2)];records=[]
    for (u,p),(v0,q0)in product(contexts,repeat=2):
        v=[i+10 for i in v0];q=q0+10
        omitted=[i for i in u+v if i not in(p,q)]
        ids=[i for i in range(19)if i not in omitted];a=g.extract(ids,ids)
        kernel=a.nullspace();entry={}
        if kernel:
            assert len(kernel)==1;k=kernel[0]
            k=k if all(x>0 for x in k)else-k
            assert all(x>0 for x in k)
            entry['type']='affine_D16';entry['positive_kernel']=list(map(str,k))
            available=[]
        else:
            ld=a.LDLdecomposition()[1];assert all(ld[i,i]>0 for i in range(17))
            inv=a.inv();assert all(x>=0 for x in inv)
            vals=[inv[i,i]for i in range(17)]
            available=[ids[i]for i,x in enumerate(vals)if x<=2]
            entry['type']='A17'if u==[1,2]and v0==[1,2]else'D17'
            assert len(available)==(4 if entry['type']=='A17'else 2)
            entry['inverse_diagonals']=list(map(str,vals))
        allowed=sorted(available+omitted);du=distances(g,u);dv=distances(g,v)
        assert all(z[i]!=0 for i in allowed)
        assert all(du[i]<=1 for i in allowed if z[i]<0)
        assert all(dv[i]<=1 for i in allowed if z[i]>0)
        assert any(z[i]<0 for i in allowed)and any(z[i]>0 for i in allowed)
        # The preferred short endpoints remain at distance14 in this corner.
        adj={i:[j for j in ids if j!=i and g[i,j]!=0]for i in ids}
        dist={p:0};queue=deque([p])
        while queue:
            i=queue.popleft()
            for j in adj[i]:
                if j not in dist:dist[j]=dist[i]+1;queue.append(j)
        assert dist[q]==14
        records.append(dict(left_black=u,right_black=v,preferred_corner=[p,q],
                            corner_nodes=ids,omitted_black=omitted,
                            possible_core_contacts=allowed,
                            negative_radical_contacts=[i for i in allowed if z[i]<0],
                            positive_radical_contacts=[i for i in allowed if z[i]>0],**entry))
    assert Counter(r['type']for r in records)=={'A17':1,'D17':4,'affine_D16':4}
    return dict(radical=list(map(str,z)),cases=records,
                shortcut_bound=4,component_order=17,
                conclusion='Every additional SUPPORTING white is orthogonal to this core if no adjacent marked pair has a white path of length<=13. At an ideal actual corner the varying affine factor is wholly in W+the two sides.')

def known_pyramid(g):
    # Independent geometric existence check in Gram-dual coordinates.
    z=g.nullspace()[0];ids=list(range(1,19));j=g.extract(ids,ids).inv()
    rays=[]
    for i in range(19):
        if z[i]==0:
            p=sp.zeros(19,1);p[i]=-1;rays.append(p)
    for i in range(19):
        for k in range(19):
            if z[i]<0<z[k]:
                p=sp.zeros(19,1);p[i]=-z[k];p[k]=z[i];rays.append(p)
    rows=[];center=sp.zeros(19,1)
    for p in rays:
        assert (z.T*p)[0]==0
        q=p.extract(ids,[0]);norm=(q.T*j*q)[0];assert norm<=0
        center+=p
        rows.append(dict(covector=list(map(str,p)),square=str(norm)))
    c=center.extract(ids,[0]);assert (c.T*j*c)[0]<0 and all(x<0 for x in center)
    assert all((p.extract(ids,[0]).T*j*c)[0]<0 for p in rays)
    counts=Counter('ideal'if r['square']=='0'else'finite'for r in rows)
    assert counts=={'finite':80,'ideal':2}
    # Each inequality is an actual facet: its tight cone rays span a hyperplane.
    for i in range(19):
        assert sp.Matrix.hstack(*[p.extract(ids,[0])for p in rays if p[i]==0]).rank()==17
    return dict(ambient_dimension=17,facets=19,vertices=dict(counts),extreme_rays=rows,
                timelike_strict_center_covector=list(map(str,center)),
                timelike_strict_center_square=str((c.T*j*c)[0]),
                reference='https://arxiv.org/pdf/math/0301133, Theorem2 and Table13',
                novelty='Known finite-volume pyramid; reconstructed here as a counterexample to the unqualified distance13 payment. No new construction dimension.')

def polynomial(ranks):
    h=[1]
    for r in ranks:
        q=[0]*(len(h)+r)
        for i,a in enumerate(h):
            for k in range(r+1):q[i+k]+=a
        h=q
    return h

def independent_moment(ranks):
    d=sum(ranks);s=len(ranks);total=0
    for size in range(s+1):
        for subset in combinations(range(s),size):
            shift=sum(ranks[i]+1 for i in subset);m=(d-1)//2-shift
            if m>=0:total+=(-1)**size*((d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1))
    return F(total,prod(r+1 for r in ranks))

def profiles(d):
    q,a=divmod(d,26);yield (26,)*q+((a,)if a else())
    for q in range(d//26+1):
        for a in range(26):
            b=d-26*q-a
            if b>=27:yield (26,)*q+((a,)if a else())+(b,)

def profile(ranks):
    d=sum(ranks);v=prod(r+1 for r in ranks);h=polynomial(ranks)
    assert h==h[::-1]and sum(h)==v
    G=sum(x*min(i,d-i)for i,x in enumerate(h))
    kappa=F(5*G-d*v,2*v);moment=independent_moment(ranks)
    assert kappa==F(3*d,4)-F(5,2)*moment
    cost=sum((min(F(13),F(r,2))for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=v,G=G,
                normalized_curvature=str(kappa),normalized_cost=str(cost),margin=str(kappa-cost))

def run():
    data={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in INPUTS}
    for n,d in data.items():
        assert d['status'].startswith('PASS'),n
        for path,sha in d.get('sources',{}).items():assert digest(ROOT/path)==sha,(n,path)
    diagrams={d['index']:d for d in data[INPUTS[0]]['diagrams']};d=diagrams[25]
    ports=high_distance_ports(data[INPUTS[0]],data[INPUTS[1]])
    local=local_obstructions(d);g=shared_core(d);isolation=isolation_cases(g)
    # The shared core has17 positive directions,one negative direction,one zero.
    ids=[i for i in range(19)if i not in(0,10)];a=g.extract(ids,ids)
    ld=a.LDLdecomposition()[1];assert all(ld[i,i]>0 for i in range(17))
    c=g.extract([0,10],ids);s=g.extract([0,10],[0,10])-c*a.inv()*c.T
    assert s==sp.ones(2)*(-sp.Rational(1,12))
    known=known_pyramid(g)
    rows=[];mins=[]
    for d in range(146,156):
        rr=list(profiles(d));assert len(rr)==len(set(tuple(sorted(r))for r in rr))
        group=[profile(r)for r in rr];assert all(F(r['margin'])>0 for r in group)
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>13*n
        mins.append(dict(ambient_dimension=n,profiles=len(group),
                         minimum_margin=str(min(F(r['margin'])for r in group)),finite_slack=str(beta-13*n)))
        rows+=group
    assert len(rows)==1255 and min(F(r['margin'])for r in rows)==F(830407,9565938)
    boundary=profile((26,)*5+(15,));assert F(boundary['margin'])==F(-515185,19131876)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{n}.{suffix}'for n in INPUTS
        for folder,suffix in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_SELECTIVE_DISTANCE13_SMALL_COMPONENT_EXCEPTION_DIMENSION146',
                sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
                dimension_bound=146,unrestricted_finite_volume=True,
                arithmetic_hypothesis=False,cusp_restriction=None,
                inherited_dimension_bound=156,excluded_dimensions=list(range(147,157)),
                weight_rule='1 if marked distance<=13 OR both marks lie in one component of order<=17;0 otherwise',
                finite_budget='13n',ideal_budget='sum min(13,r_i/2)',
                radius_contexts=ports,disjoint_obstructions=local,
                shared_core_gram=serial_matrix(g),shared_core_center_schur=serial_matrix(s),
                supporting_white_isolation=isolation,known_sharp_distance14_pyramid=known,
                reduced_profiles=rows,dimension_minima=mins,boundary_at_n146=boundary,
                compact_facet_ambient_bound=111,
                compact_facet_parity_slack_times8=[[12,224,4],[120,228,4]],
                scope='Written all-label geometric classification and core isolation, inherited actual cusp fibers and convex-order reduction,exact scalar certificates. Internal computer-assisted proof; no external referee or proof-assistant verification. The reconstructed H17 pyramid is already known.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['reduced_profiles']),'profiles',flush=True)
