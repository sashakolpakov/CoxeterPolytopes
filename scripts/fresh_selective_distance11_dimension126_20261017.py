#!/usr/bin/env python3
"""Exact distance11 joins/overlaps, small isolated corners and dimension126.

SymPy rational congruences supply exact inertia and nonnegative Z witnesses.
The separately frozen 17/18-root lemmas are bound rather than reimplemented.
No numerical linear algebra or optimizer is used.
"""
import argparse
from collections import Counter,deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
import json
from math import prod
from pathlib import Path
import sympy as sp
import fresh_selective_distance12_dimension136_20261017 as D

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_selective_distance11_dimension126_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
RADIUS='fresh_simple_cusp_distance_bound_20261017'
TRIPLE='fresh_distance11_extension_review_20261017'
MIXED=TRIPLE+'_mixed'
INPUTS=(CENSUS,RADIUS,D.TAG,D.TAG+'_review',TRIPLE,MIXED,
        'fresh_cusp_factor_payment_framework_20261017','fresh_cusp_link_shape_curvature_20261017')
P=D.P


def digest(p):
    h=sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def matrix(g):return [[str(g[i,j])for j in range(g.cols)]for i in range(g.rows)]


def rationalize(g):
    """Positive diagonal sqrt(2) congruence, checked entry by entry."""
    n=g.rows;bits={}
    for root in range(n):
        if root in bits:continue
        bits[root]=0;q=[root]
        for i in q:
            for j in range(n):
                if i==j or not g[i,j]:continue
                assert g[i,j]in(-1,-sp.sqrt(2))
                b=bits[i]^int(g[i,j]==-sp.sqrt(2))
                if j in bits:assert bits[j]==b
                else:bits[j]=b;q.append(j)
    scales=[2**bits[i]for i in range(n)]
    a=sp.zeros(n)
    for i in range(n):
        for j in range(n):
            a[i,j]=sp.simplify(g[i,j]*sp.sqrt(scales[i]*scales[j]))
            assert a[i,j].is_Rational
    return a,scales


def exact_inertia(g):
    """Symmetric 1x1/2x2 Schur congruences, including singular blocks."""
    a=g.copy();ids=list(range(g.rows));pos=neg=zero=0;steps=[]
    while a.rows:
        i=next((i for i in range(a.rows)if a[i,i]),None)
        if i is not None:
            p=a[i,i];pos+=int(bool(p>0));neg+=int(bool(p<0))
            rest=[j for j in range(a.rows)if j!=i]
            steps.append(dict(nodes=[ids[i]],pivot=[[str(p)]]))
            a=a.extract(rest,rest)-a.extract(rest,[i])*a.extract([i],rest)/p
            ids=[ids[j]for j in rest]
        else:
            pair=next(((i,j)for i in range(a.rows)for j in range(i+1,a.rows)if a[i,j]),None)
            if pair is None:
                zero+=a.rows;steps.append(dict(zero_nodes=ids));break
            i,j=pair;rest=[k for k in range(a.rows)if k not in pair]
            block=a.extract(pair,pair);assert block.det()<0
            pos+=1;neg+=1;steps.append(dict(nodes=[ids[i],ids[j]],pivot=matrix(block)))
            a=a.extract(rest,rest)-a.extract(rest,pair)*block.inv()*a.extract(pair,rest)
            ids=[ids[k]for k in rest]
    return dict(signature=[pos,neg,zero],congruence_pivots=steps)


def z_obstruction(a):
    """First nonpositive pivot gives a nonnegative vector of square<=0.

    A positive Schur pivot preserves the Z-matrix sign pattern; its basis
    update adds a nonnegative multiple of the eliminated basis vector.
    Thus the certificate remains an obstruction after extra negative bonds.
    """
    n=a.rows;h=a.copy();basis=sp.eye(n)
    assert a==a.T and all(a[i,j]<=0 for i in range(n)for j in range(n)if i!=j)
    for k in range(n):
        pivot=h[k,k]
        if pivot<=0:
            v=basis.row(k).T
            assert all(x>=0 for x in v)and any(v)
            assert (v.T*a*v)[0]==pivot
            return dict(vector=list(map(str,v)),square=str(pivot),first_nonpositive_pivot=k)
        for i in range(k+1,n):
            coefficient=h[i,k]/pivot;assert coefficient<=0
            basis[i,:]=basis.row(i)-coefficient*basis.row(k)
        rest=list(range(k+1,n))
        if rest:
            updated=h.extract(rest,rest)-h.extract(rest,[k])*h.extract([k],rest)/pivot
            for i,ii in enumerate(rest):
                for j,jj in enumerate(rest):h[ii,jj]=updated[i,j]
    return None


def edge_obstruction(g,black):
    a,scales=rationalize(g);white=[i for i in range(g.rows)if i not in black]
    for mark in black:
        ids=sorted(white+[mark]);ob=z_obstruction(a.extract(ids,ids))
        if ob:return dict(mark=mark,nodes=ids,gauge_squared=scales,**ob)
    return None


def contexts_and_ports(data):
    diagrams={d['index']:d for d in data[CENSUS]['diagrams']}
    grams={i:P.gram(d)for i,d in diagrams.items()}
    contexts=[]
    for i,d in sorted(diagrams.items()):
        for black in combinations(d['finite_deletions'],2):
            contexts.append(dict(diagram=i,black=list(black),distances=P.distances(grams[i],black)))
    old=data[RADIUS]['black_pair_checks'];assert len(contexts)==len(old)==384
    assert max(max(r['distances'])for r in contexts)==7
    assert {(r['diagram'],tuple(r['black']),tuple(r['distances']))for r in contexts}=={
        (r['diagram'],tuple(r['black']),tuple(r['distances']))for r in old}
    @lru_cache(None)
    def deletion_inverse(di,b):
        a,scales=rationalize(grams[di]);ids=[j for j in range(a.rows)if j!=b]
        h=a.extract(ids,ids);D.pd(h)
        return ids,scales,h.inv()
    @lru_cache(None)
    def values(di,black,node):
        out=[]
        for b in black:
            ids,scales,j=deletion_inverse(di,b)
            out.append(scales[node]*j[ids.index(node),ids.index(node)])
        return tuple(out)
    ports=[]
    for r in contexts:
        for node,dist in enumerate(r['distances']):
            if dist>=4:
                vals=values(r['diagram'],tuple(r['black']),node)
                ports.append(dict(diagram=r['diagram'],black=r['black'],node=node,distance=dist,
                    inverse_diagonals=list(map(str,vals)),allowed=bool(max(vals)<2)))
    assert sum(p['allowed']for p in ports)==55
    return diagrams,grams,contexts,ports,values


def disjoint(grams,ports):
    allowed=[p for p in ports if p['allowed']];rows=[];good={};type_cert={}
    for left,right in product(allowed,repeat=2):
        if left['distance']+right['distance']<11:continue
        g,h=grams[left['diagram']],grams[right['diagram']]
        u,mp=D.union(g,h,());u[left['node'],mp[right['node']]]=u[mp[right['node']],left['node']]=-1
        a,scales=rationalize(u);black=left['black']+[mp[i]for i in right['black']]
        white=[i for i in range(u.rows)if i not in black]
        ob=z_obstruction(a.extract(white,white))
        row=dict(left=[left['diagram'],left['black'],left['node']],right=[right['diagram'],right['black'],right['node']])
        if ob:
            row['white_rejection']=dict(nodes=white,gauge_squared=scales,**ob)
        else:
            key=(left['diagram'],tuple(left['black']),right['diagram'],tuple(right['black']))
            good.setdefault(key,[]).append((left['node'],right['node']))
            tails=[]
            for r in(left,right):
                native=grams[r['diagram']];ws=[i for i in range(native.rows)if i not in r['black']];v=r['node']
                choices=[[v,x,y]for x,y in product(ws,repeat=2)if len({v,x,y})==3 and native[v,x]and native[x,y]and not native[v,y]]
                assert choices
                tail=min(choices);tails.append(dict(nodes=tail,native_products=[str(native[tail[0],tail[1]]),str(native[tail[1],tail[2]])]))
            row['white_three_node_tails']=tails
            typ=tuple(sorted(((left['diagram'],left['node']),(right['diagram'],right['node']))))
            if typ not in type_cert:
                (di,i),(dj,j)=typ
                core,cm=D.union(grams[di],grams[dj],());core[i,cm[j]]=core[cm[j],i]=-1
                rat,sigmas=rationalize(core);cert=exact_inertia(rat);assert cert['signature'][1]>=2
                type_cert[typ]=dict(ports=[list(p)for p in typ],gram=matrix(rat),gauge_squared=sigmas,**cert)
            row['inertia_type']=[list(p)for p in typ]
        rows.append(row)
    assert len(rows)==425 and sum('white_rejection'in r for r in rows)==119
    assert len(good)==306 and all(len(es)==1 for es in good.values())and len(type_cert)==13
    # For the minimum six-root two-tail model the positive vector
    # (1,2/3,1/3;1,2/3,1/3) has square8/3-2sqrt2<0 at a4 cross bond.
    assert sp.Rational(8,3)**2<8
    return dict(candidate_joins=425,white_rejections=119,remaining_unique_context_pairs=306,
        rows=rows,underlying_inertia_types=[type_cert[t]for t in sorted(type_cert)],
        bond_forcing=dict(tail_vector=['1','2/3','1/3','1','2/3','1/3'],
            square_at_minimum_nonsimple_bond='8/3-2*sqrt(2)',
            proof='All native tail bonds have magnitude at least1. A cross bond of label>=4 has magnitude>=sqrt2; the nonnegative vector gives a negative square. Thus the unique possible cross bond is simple.'))


def classify_overlap(left,right,match):
    di,dj=left['diagram'],right['diagram']
    if di==dj==25:
        if match==((9,9),):return 'shared_vertex'
        if match==((8,9),(9,8)):return 'shared_edge'
        if match==((7,9),(8,8),(9,7)):return 'shared_three_path'
    if di==22 and dj==25 and match==((8,9),):return 'mixed_22_25'
    if di==25 and dj==22 and match==((9,8),):return 'mixed_25_22'
    raise AssertionError((di,dj,match))


def overlaps(grams,contexts,values):
    high=[r for r in contexts if max(r['distances'])>=5]
    assert max(sum(d>=5 for d in r['distances'])for r in high)==3
    records=[];generated=0
    for left,right in product(high,repeat=2):
        g,h=grams[left['diagram']],grams[right['diagram']]
        cand=[(i,j)for i in range(g.rows)for j in range(h.rows)if left['distances'][i]+right['distances'][j]>=12]
        sizes=min(sum(d>=5 for d in left['distances']),sum(d>=5 for d in right['distances']))
        for size in range(1,sizes+1):
            for match in combinations(cand,size):
                generated+=1;ls=[i for i,j in match];rs=[j for i,j in match]
                if len(set(ls))<size or len(set(rs))<size:continue
                if any(g[i,k]!=h[j,l]for i,j in match for k,l in match):continue
                row=dict(left=[left['diagram'],left['black']],right=[right['diagram'],right['black']],match=[list(p)for p in match])
                rejection=None
                for i,j in match:
                    for side,x,y,rr,other,shared in [('left',i,j,left,h,rs),('right',j,i,right,g,ls)]:
                        ext=[k for k in range(other.rows)if k!=y and other[y,k]and k not in shared]
                        if ext and max(values(rr['diagram'],tuple(rr['black']),x))>=2:
                            other_black=right['black']if side=='left'else left['black']
                            assert all(k not in other_black for k in ext)
                            rejection=dict(side=side,node=x,external_neighbor=ext[0],
                                inverse_diagonals=list(map(str,values(rr['diagram'],tuple(rr['black']),x))))
                            break
                    if rejection:break
                if rejection:row['inverse_rejection']=rejection
                else:
                    u,mp=D.union(g,h,match);black=left['black']+[mp[i]for i in right['black']]
                    ob=edge_obstruction(u,black)
                    if ob:row['edge_pd_rejection']=ob
                    else:row['surviving_type']=classify_overlap(left,right,match)
                records.append(row)
    counts=Counter('inverse_rejection'if'inverse_rejection'in r else'edge_pd_rejection'if'edge_pd_rejection'in r else r['surviving_type']for r in records)
    assert generated==2137 and len(records)==1010
    assert counts=={'inverse_rejection':856,'edge_pd_rejection':40,'shared_vertex':60,'shared_edge':27,
        'shared_three_path':9,'mixed_22_25':9,'mixed_25_22':9}
    return dict(overlap_generated=generated,consistent_overlaps=len(records),overlap_counts=dict(counts),overlap_records=records)


def preferred(black):return 4 if 4 in black else D.preferred(black)


def preferred_corner(g,u,v,p,q):
    omitted=sorted(set(u+v)-{p,q});ids=[i for i in range(g.rows)if i not in omitted]
    h=g.extract(ids,ids);allowed=list(omitted);blocks=[];marked=None
    for part in D.components(h):
        actual=[ids[i]for i in part];a=h.extract(part,part);ks=a.nullspace()
        if ks:
            assert len(ks)==1;k=ks[0]
            if all(x<0 for x in k):k=-k
            assert all(x>0 for x in k);D.pd(a[:-1,:-1])
            blocks.append(dict(nodes=actual,kind='affine',positive_kernel=list(map(str,k))))
        else:
            D.pd(a);j=a.inv();assert all(x>=0 for x in j)
            allowed.extend(actual[i]for i in range(a.rows)if j[i,i]<=2)
            blocks.append(dict(nodes=actual,kind='finite',inverse_diagonals=list(map(str,j.diagonal()))))
        if p in actual:assert q in actual;marked=actual
    assert marked is not None and len(marked)<=17
    # Distances inside this one preferred corner avoid both omitted blacks.
    dist={p:0};todo=deque([p])
    while todo:
        i=todo.popleft()
        for j in ids:
            if i!=j and g[i,j]and j not in dist:dist[j]=dist[i]+1;todo.append(j)
    assert dist[q]>11
    return dict(left_black=u,right_black=v,preferred_corner=[p,q],omitted_black=omitted,
        blocks=blocks,possible_contacts=sorted(allowed),marked_component=marked,marked_distance=dist[q])


def extra_cross_bonds(g,mp,match,u,v,left_order,right_order):
    left=[i for i in range(left_order)if i not in {i for i,j in match}]
    right=[mp[j]for j in range(right_order)if j not in {j for i,j in match}]
    short=0;obstructions=[]
    for x,y in product(left,right):
        if g[x,y]:continue
        h=g.copy();h[x,y]=h[y,x]=-1
        if min(P.distances(h,u)[j]for j in v)<=11:short+=1;continue
        cert=edge_obstruction(h,u+v);assert cert is not None,(match,u,v,x,y)
        obstructions.append(dict(edge=[x,y],**cert))
    return dict(shortcut_rejections=short,edge_pd_z_rejections=obstructions,
        monotonicity='The saved nonnegative vector has square<=0 at the minimum native bond. Larger or additional negative bonds cannot increase its square in the required positive-definite edge support.')


def isolation(grams,table,data):
    records=[];cores={};triple=data[TRIPLE];mixed=data[MIXED]
    assert triple['marked_component_order_bound']==mixed['marked_component_order_bound']==16
    for row in table['overlap_records']:
        typ=row.get('surviving_type')
        if typ is None:continue
        di,u=row['left'];dj,v0=row['right'];match=tuple(map(tuple,row['match']))
        g,mp=D.union(grams[di],grams[dj],match);v=[mp[j]for j in v0]
        assert min(P.distances(g,u)[j]for j in v)>11
        cross=extra_cross_bonds(g,mp,match,u,v,grams[di].rows,grams[dj].rows)
        if typ in('shared_vertex','shared_edge'):
            if typ not in cores:
                ids=[i for i in range(g.rows)if i not in(0,mp[0])];a=g.extract(ids,ids);D.pd(a)
                c=g.extract([0,mp[0]],ids);schur=g.extract([0,mp[0]],[0,mp[0]])-c*a.inv()*c.T
                if typ=='shared_vertex':
                    assert schur==sp.ones(2)*(-sp.Rational(1,12));z=g.nullspace()[0]
                    assert all(z[i]<0 for i in range(9))and z[9]==0 and all(z[i]>0 for i in range(10,19))
                    cores[typ]=dict(gram=matrix(g),schur=matrix(schur),signature=[17,1,1],radical=list(map(str,z)))
                else:
                    assert schur[0,0]<0 and schur.det()<0;j=g.inv()
                    ends=list(range(4));right=[mp[i]for i in ends];block=j.extract(ends,ends)
                    assert all(x>0 for x in block)and min(block.diagonal())>=4
                    assert j.extract(right,right)==block
                    cores[typ]=dict(gram=matrix(g),schur=matrix(schur),signature=[17,1,0],
                        determinant=str(g.det()),end_nodes=[ends,right],end_inverse_block=matrix(block))
            r=preferred_corner(g,u,v,preferred(u),mp[preferred(v0)]);allowed=r['possible_contacts']
            du,dv=P.distances(g,u),P.distances(g,v)
            if typ=='shared_vertex':
                z=sp.Matrix(list(map(sp.Rational,cores[typ]['radical'])))
                assert all(z[i]for i in allowed)
                assert all(du[i]<=3 for i in allowed if z[i]<0)
                assert all(dv[i]<=3 for i in allowed if z[i]>0)
                r['shortcut_bound']=8
            else:
                ends,right=cores[typ]['end_nodes'];assert set(allowed)<=set(ends+right)
                assert all(du[i]<=2 for i in allowed if i in ends)
                assert all(dv[i]<=2 for i in allowed if i in right)
                r['shortcut_bound']=6
            r.update(type=typ,extra_cross_bond_checks=cross);records.append(r)
        else:
            if typ=='shared_three_path':
                target=sp.Matrix(triple['core']['gram']);assert g==target
                matches=[c for c in triple['marked_contexts']if c['left_black']==u and c['right_black']==v]
                assert len(matches)==1;lemma=TRIPLE;mapping=list(range(g.rows))
            else:
                target=sp.Matrix(mixed['core']['gram']);rm=mixed['core']['right_copy_map'];mapping={}
                for i in range(grams[di].rows):mapping[i]=i if di==22 else rm[i]
                for j in range(grams[dj].rows):
                    dest=j if dj==22 else rm[j]
                    if mp[j]in mapping:assert mapping[mp[j]]==dest
                    else:mapping[mp[j]]=dest
                mapping=[mapping[i]for i in range(g.rows)];assert len(set(mapping))==g.rows
                assert g==target.extract(mapping,mapping)
                lu=u if di==22 else v0;rv=[rm[i]for i in(v0 if dj==25 else u)]
                matches=[c for c in mixed['contexts']if c['left_black']==lu and c['right_black']==rv]
                assert len(matches)==1;lemma=MIXED
            records.append(dict(type=typ,left_black=u,right_black=v,local_lemma=lemma,
                local_lemma_context=matches[0],literal_isomorphism=mapping,
                component_order_bound=16,extra_cross_bond_checks=cross))
    assert len(records)==114
    return dict(cores=cores,cases=records,component_order_cap=17,
        type_counts=dict(Counter(r['type']for r in records)))


def profiles(d):
    q,a=divmod(d,22);yield(22,)*q+((a,)if a else())
    for q in range(d//22+1):
        for a in range(22):
            b=d-22*q-a
            if b>=23:yield(22,)*q+((a,)if a else())+(b,)


def budget(rs):
    d=sum(rs);v=prod(r+1 for r in rs);h=P.polynomial(rs)
    assert sum(h)==v and h==h[::-1]
    G=sum(x*min(i,d-i)for i,x in enumerate(h));curvature=F(5*G-d*v,2*v)
    assert curvature==F(3*d,4)-F(5,2)*P.independent_moment(rs)
    cost=sum((min(F(11),F(r,2))for r in rs),F())
    return dict(ranks=list(rs),link_dimension=d,link_vertices=v,G=G,
        normalized_curvature=str(curvature),normalized_cost=str(cost),margin=str(curvature-cost))


def run():
    data={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in INPUTS}
    for n,d in data.items():
        assert d['status'].startswith('PASS'),n
        for path,h in d.get('sources',{}).items():assert digest(ROOT/path)==h,(n,path)
    diagrams,grams,contexts,ports,values=contexts_and_ports(data)
    joins=disjoint(grams,ports);print('DISJOINT_PASS',joins['candidate_joins'],joins['white_rejections'],len(joins['underlying_inertia_types']),flush=True)
    table=overlaps(grams,contexts,values);print('OVERLAP_PASS',table['overlap_counts'],flush=True)
    iso=isolation(grams,table,data);print('ISOLATION_PASS',len(iso['cases']),flush=True)
    rows=[];mins=[]
    for d in range(126,136):
        rr=list(profiles(d));assert len(rr)==len(set(tuple(sorted(rs))for rs in rr))
        group=[budget(rs)for rs in rr];assert all(F(r['margin'])>0 for r in group)
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>11*n
        mins.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=str(min(F(r['margin'])for r in group)),finite_slack=str(beta-11*n)))
        rows.extend(group)
    assert len(rows)==1095 and min(F(r['margin'])for r in rows)==F(40603,559682)
    boundary=budget((22,)*5+(15,));assert F(boundary['margin'])==F(-238535,12872686)
    t=sp.Symbol('t');facet=[]
    for start in(95,96):
        d=start+2*t;beta=(d-1)*(d-5)/8 if start%2 else d*(d-6)/8
        poly=sp.Poly(sp.expand(8*(beta-11*(d+1))),t)
        facet.append([int(poly.nth(k))for k in range(3)])
    assert facet==[[12,192,4],[104,196,4]]
    paths=[Path(__file__)]+[ROOT/f'{folder}/{n}.{ext}'for n in INPUTS for folder,ext in[('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_SELECTIVE_DISTANCE11_SMALL_COMPONENT_EXCEPTION_DIMENSION126',
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},dimension_bound=126,
        inherited_dimension_bound=136,unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        excluded_dimensions=list(range(127,137)),weight_rule='1 if marked distance<=11 OR both marks lie in one component of order<=17;0 otherwise',
        finite_budget='11n',ideal_budget='sum min(11,r_i/2)',
        census=dict(marked_contexts=len(contexts),ports=ports,allowed_ports=sum(r['allowed']for r in ports),**table),
        disjoint_obstructions=joins,isolation=iso,reduced_profiles=rows,dimension_minima=mins,boundary_at_n126=boundary,
        compact_facet_ambient_bound=95,compact_facet_parity_slack_times8=facet,
        scope='Exact finite arithmetic and explicit inherited local lemmas for the selective distance11 argument. Completeness of the published finite/minimal-indefinite classifications, actual corner incidence, all-rank distance budgets and convex-order reduction remains a mathematical proof input. No novelty or new construction claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'profiles',len(d['reduced_profiles']),'overlaps',d['census']['overlap_counts'],flush=True)
