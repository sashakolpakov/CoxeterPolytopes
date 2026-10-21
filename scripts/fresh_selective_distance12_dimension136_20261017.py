#!/usr/bin/env python3
"""Exact distance12 core split, isolated small corners and dimension136."""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import prod
from pathlib import Path
import sympy as sp
import fresh_selective_distance13_dimension146_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_selective_distance12_dimension136_20261017'
INPUTS=('fresh_high_rank_quasilanner_census_20261017',
        'fresh_simple_cusp_distance_bound_20261017',
        'fresh_selective_distance13_dimension146_20261017',
        'fresh_selective_distance13_dimension146_20261017_review',
        'fresh_cusp_factor_payment_framework_20261017',
        'fresh_cusp_link_shape_curvature_20261017')

def digest(p):
    h=sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def matrix(m):return [[str(x)for x in m.row(i)]for i in range(m.rows)]

def pd(g):
    d=g.LDLdecomposition()[1]
    assert all(x>0 for x in d.diagonal())

def components(g):
    left=set(range(g.rows));out=[]
    while left:
        i=min(left);left.remove(i);part={i};q=[i]
        while q:
            i=q.pop()
            for j in sorted(left):
                if g[i,j]!=0:left.remove(j);part.add(j);q.append(j)
        out.append(sorted(part))
    return out

def kernel(g,ids):
    a=g.extract(ids,ids);ks=a.nullspace();assert len(ks)==1
    k=ks[0];k=k if all(x>0 for x in k)else-k
    assert all(x>0 for x in k) and a*k==sp.zeros(len(ids),1)
    pd(a[:-1,:-1])
    return dict(nodes=ids,positive_kernel=list(map(str,k)))

def union(g,h,match):
    mp={j:i for i,j in match};n=g.rows
    for j in range(h.rows):
        if j not in mp:mp[j]=n;n+=1
    a=2*sp.eye(n)
    for i in range(g.rows):
        for j in range(g.rows):a[i,j]=g[i,j]
    for i in range(h.rows):
        for j in range(h.rows):
            if mp[i]<g.rows and mp[j]<g.rows:assert a[mp[i],mp[j]]==h[i,j]
            a[mp[i],mp[j]]=h[i,j]
    return a,[mp[j]for j in range(h.rows)]

def census(data):
    diagrams={d['index']:d for d in data[INPUTS[0]]['diagrams']}
    contexts=data[INPUTS[1]]['black_pair_checks']
    grams={i:P.gram(d)for i,d in diagrams.items()}
    @lru_cache(None)
    def values(di,black,v):
        g=grams[di];out=[]
        for b in black:
            ids=[i for i in range(g.rows)if i!=b];a=g.extract(ids,ids);pd(a)
            out.append(a.inv()[ids.index(v),ids.index(v)])
        return tuple(out)
    ports=[]
    for r in contexts:
        assert P.distances(grams[r['diagram']],r['black'])==r['distances']
        for i,d in enumerate(r['distances']):
            if d>=5:
                vv=values(r['diagram'],tuple(r['black']),i)
                ports.append(dict(diagram=r['diagram'],black=r['black'],node=i,distance=d,
                                  inverse_diagonals=list(map(str,vv)),allowed=bool(max(vv)<2)))
    allowed=[p for p in ports if p['allowed']];assert len(allowed)==18
    joins=[]
    for a,b in product(allowed,repeat=2):
        if a['distance']+b['distance']>=12:
            joins.append(dict(left=[a['diagram'],a['black'],a['node']],
                              right=[b['diagram'],b['black'],b['node']]))
    types=Counter(str(tuple(sorted(((r['left'][0],r['left'][2]),(r['right'][0],r['right'][2])))))for r in joins)
    assert sorted(types.values())==[1,6,6,18,24,60]
    high=[r for r in contexts if max(r['distances'])>=6]
    assert len(high)==16 and all(sum(d>=6 for d in r['distances'])<=2 for r in high)
    overlap=[];generated=0
    for a,b in product(high,repeat=2):
        g=grams[a['diagram']];h=grams[b['diagram']]
        cand=[(i,j)for i in range(g.rows)for j in range(h.rows)if a['distances'][i]+b['distances'][j]>=13]
        for size in(1,2):
            for match in combinations(cand,size):
                generated+=1;lm=[i for i,j in match];rm=[j for i,j in match]
                if len(set(lm))<size or len(set(rm))<size:continue
                if any(g[i,k]!=h[j,l]for i,j in match for k,l in match):continue
                reject=None
                for i,j in match:
                    for swap,x,y,gm,hm,rr,shared in[(False,i,j,g,h,a,rm),(True,j,i,h,g,b,lm)]:
                        ext=[k for k in range(hm.rows)if k!=y and hm[y,k]!=0 and k not in shared]
                        if ext and max(values(rr['diagram'],tuple(rr['black']),x))>=2:
                            reject=dict(side='right'if swap else'left',node=x,external_neighbor=ext[0],
                                        inverse_diagonals=list(map(str,values(rr['diagram'],tuple(rr['black']),x))))
                            break
                    if reject:break
                row=dict(left=[a['diagram'],a['black']],right=[b['diagram'],b['black']],match=[list(p)for p in match])
                if reject:row['inverse_rejection']=reject
                else:
                    assert a['diagram']==b['diagram']==25
                    if match in(((8,9),),((9,8),)):
                        z,mp=union(g,h,match)
                        if match==((9,8),):white=[9,mp[9],mp[7],mp[6],8,7,6,5,4]
                        else:white=[8,9,7,6,mp[8],mp[7],mp[6],mp[5],mp[4]]
                        assert not(set(white)&set(a['black']+[mp[i]for i in b['black']]))
                        row['white_affine_rejection']=kernel(z,white)
                    else:
                        assert match in(((9,9),),((8,9),(9,8)))
                        row['surviving_type']='shared_vertex'if size==1 else'shared_edge'
                overlap.append(row)
    counts=Counter('inverse_rejection'if'inverse_rejection'in r else'white_affine_rejection'if'white_affine_rejection'in r else r['surviving_type']for r in overlap)
    assert generated==208 and len(overlap)==162
    assert counts=={'inverse_rejection':120,'white_affine_rejection':6,'shared_vertex':27,'shared_edge':9}
    return diagrams,contexts,dict(marked_contexts=len(contexts),ports=ports,disjoint_joins=joins,
                                  disjoint_type_counts=dict(types),overlap_generated=generated,
                                  overlap_counts=dict(counts),overlap_records=overlap)

def disjoint(diagrams):
    g=P.gram(diagrams[25]);obstructions=[]
    for ix,v in[(18,7),(22,8),(22,4),(25,9)]:
        h=sp.diag(g,P.gram(diagrams[ix]));h[9,10+v]=h[10+v,9]=-1
        ids=[i for i in range(h.rows)if i not in(0,10)];a=h.extract(ids,ids);pd(a)
        c=h.extract([0,10],ids);s=h.extract([0,10],[0,10])-c*a.inv()*c.T
        assert s[0,0]<0 and s.det()>0
        obstructions.append(dict(right_diagram=ix,right_port=v,positive_nodes=ids,
                                 schur=matrix(s),determinant=str(s.det())))
    h=sp.diag(g,g);h[9,18]=h[18,9]=-1
    far_near=kernel(h,[18,19,17,16,9,8,7,6,5])
    h=sp.diag(g,g);h[8,18]=h[18,8]=-1
    near_near=kernel(h,[8,18,7,9,17,19])
    return dict(simple_leaf_joins=obstructions,far_near_white=far_near,near_near_white=near_near,
                simple_bond_forcing='Two disjoint wholly-white A3 tails have residual4/3; bond squared<16/9<2.',
                additional_bonds='After rejecting the two white-affine types, each marked context pair permits only one cross edge.')

def preferred(u):return 0 if 0 in u else 2 if tuple(u)in((1,2),(2,3))else 1

def corner(g,u,v,p,q):
    omitted=[i for i in u+v if i not in(p,q)];ids=[i for i in range(g.rows)if i not in omitted]
    h=g.extract(ids,ids);available=[];blocks=[];marked=None
    for part in components(h):
        actual=[ids[i]for i in part];a=h.extract(part,part);ks=a.nullspace()
        if ks:
            k=ks[0];k=k if all(x>0 for x in k)else-k
            assert len(ks)==1 and all(x>0 for x in k);pd(a[:-1,:-1])
            blocks.append(dict(nodes=actual,kind='affine',positive_kernel=list(map(str,k))))
        else:
            pd(a);iv=a.inv();assert all(x>=0 for x in iv)
            available.extend(actual[i]for i in range(a.rows)if iv[i,i]<=2)
            blocks.append(dict(nodes=actual,kind='finite',inverse_diagonals=list(map(str,iv.diagonal()))))
        if p in actual:assert q in actual;marked=actual
    assert marked is not None and len(marked)<=17
    # Shortest distance is evaluated only in this corner, not through omitted blacks.
    dist={p:0};queue=deque([p])
    while queue:
        i=queue.popleft()
        for j in ids:
            if i!=j and g[i,j]!=0 and j not in dist:dist[j]=dist[i]+1;queue.append(j)
    assert dist[q]>=13
    return dict(left_black=u,right_black=v,preferred_corner=[p,q],omitted_black=omitted,
                blocks=blocks,possible_contacts=sorted(available+omitted),marked_component=marked,
                marked_distance=dist[q])

def isolation(diagrams,table):
    native=P.gram(diagrams[25]);records=[];core_records={}
    for typ,match in [('shared_vertex',((9,9),)),('shared_edge',((8,9),(9,8)))]:
        g,mp=union(native,native,match);n=g.rows
        ids=[i for i in range(n)if i not in(0,mp[0])];a=g.extract(ids,ids);pd(a)
        c=g.extract([0,mp[0]],ids);s=g.extract([0,mp[0]],[0,mp[0]])-c*a.inv()*c.T
        if typ=='shared_vertex':
            assert s==sp.ones(2)*(-sp.Rational(1,12));z=g.nullspace()[0]
            assert all(z[i]<0 for i in range(9))and z[9]==0 and all(z[i]>0 for i in range(10,19))
            core_records[typ]=dict(gram=matrix(g),schur=matrix(s),radical=list(map(str,z)))
        else:
            assert s[0,0]<0 and s.det()<0;iv=g.inv();ends=[1,2,3];right=[mp[i]for i in ends]
            b=iv.extract(ends,ends);assert b==sp.Matrix([[8,10,5],[10,14,7],[5,7,4]])
            assert iv.extract(right,right)==b
            core_records[typ]=dict(gram=matrix(g),schur=matrix(s),determinant=str(g.det()),
                                   end_nodes=[ends,right],end_inverse_block=matrix(b))
        for row in table['overlap_records']:
            if row.get('surviving_type')!=typ:continue
            u=row['left'][1];v0=row['right'][1];v=[mp[i]for i in v0]
            r=corner(g,u,v,preferred(u),mp[preferred(v0)]);allowed=r['possible_contacts']
            du=P.distances(g,u);dv=P.distances(g,v)
            if typ=='shared_vertex':
                assert all(z[i]!=0 for i in allowed)
                assert all(du[i]<=2 for i in allowed if z[i]<0)
                assert all(dv[i]<=2 for i in allowed if z[i]>0)
                r['shortcut_bound']=6
            else:
                assert set(allowed)<=set(ends+right)
                assert all(du[i]<=1 for i in allowed if i in ends)
                assert all(dv[i]<=1 for i in allowed if i in right)
                r['shortcut_bound']=4
            # No unknown core-cross bond survives the no-short-path hypothesis,
            # except the affine triangle across the shared far vertex.
            left_only=[i for i in range(10)if i not in {a for a,b in match}]
            right_only=[mp[i]for i in range(10)if i not in {b for a,b in match}]
            extra=[]
            for x,y in product(left_only,right_only):
                if g[x,y]!=0:continue
                h=g.copy();h[x,y]=h[y,x]=-1
                if min(P.distances(h,u)[b]for b in v)>12:
                    assert typ=='shared_vertex'and {x,y}=={8,mp[8]}
                    assert x not in u+v and y not in u+v
                    extra.append(kernel(h,[8,9,mp[8]]))
            r['extra_cross_bonds_rejected_by_white_triangle']=extra
            r['type']=typ;records.append(r)
    assert Counter(r['type']for r in records)=={'shared_vertex':27,'shared_edge':9}
    return dict(cores=core_records,cases=records,component_order_cap=17)

def profiles(d):
    q,a=divmod(d,24);yield (24,)*q+((a,)if a else())
    for q in range(d//24+1):
        for a in range(24):
            b=d-24*q-a
            if b>=25:yield (24,)*q+((a,)if a else())+(b,)

def budget(rs):
    d=sum(rs);v=prod(r+1 for r in rs);h=P.polynomial(rs)
    G=sum(x*min(i,d-i)for i,x in enumerate(h));curv=F(5*G-d*v,2*v)
    assert curv==F(3*d,4)-F(5,2)*P.independent_moment(rs)
    cost=sum((min(F(12),F(r,2))for r in rs),F())
    return dict(ranks=list(rs),link_dimension=d,link_vertices=v,G=G,
                normalized_curvature=str(curv),normalized_cost=str(cost),margin=str(curv-cost))

def run():
    data={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in INPUTS}
    for n,d in data.items():
        assert d['status'].startswith('PASS')
        for path,h in d.get('sources',{}).items():assert digest(ROOT/path)==h,(n,path)
    diagrams,contexts,table=census(data);local=disjoint(diagrams);iso=isolation(diagrams,table)
    # All surviving disjoint pairs have a unique non-affine possible cross edge.
    good={((18,7),(25,9)),((22,4),(25,9)),((22,8),(25,9)),((25,9),(25,9))}
    surviving={}
    for r in table['disjoint_joins']:
        a,b=r['left'],r['right'];typ=tuple(sorted(((a[0],a[2]),(b[0],b[2]))))
        if typ in good:
            key=(a[0],tuple(a[1]),b[0],tuple(b[1]));surviving.setdefault(key,[]).append((a[2],b[2]))
            # Each white A3 tail is available and orthogonal to the other before the cross edge.
            for di,bl,v in(a,b):
                gr=P.gram(diagrams[di]);white=[i for i in range(gr.rows)if i not in bl]
                tails=[]
                for x,y in product(white,repeat=2):
                    if len({v,x,y})==3 and gr[v,x]==gr[x,y]==-1 and gr[v,y]==0:tails.append([v,x,y])
                assert tails
    assert all(len(v)==1 for v in surviving.values())
    rows=[];mins=[]
    for d in range(136,146):
        rr=list(profiles(d));assert len(rr)==len(set(tuple(sorted(r))for r in rr))
        group=[budget(r)for r in rr];assert all(F(r['margin'])>0 for r in group)
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>12*n
        mins.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=str(min(F(r['margin'])for r in group)),finite_slack=str(beta-12*n)))
        rows.extend(group)
    assert len(rows)==1175 and min(F(r['margin'])for r in rows)==F(172096,1953125)
    boundary=budget((24,)*5+(15,));assert F(boundary['margin'])==F(-59157,3906250)
    paths=[Path(__file__)]+[ROOT/f'{f}/{n}.{e}'for n in INPUTS for f,e in[('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_SELECTIVE_DISTANCE12_SMALL_COMPONENT_EXCEPTION_DIMENSION136',
                sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
                dimension_bound=136,inherited_dimension_bound=146,unrestricted_finite_volume=True,
                arithmetic_hypothesis=False,cusp_restriction=None,excluded_dimensions=list(range(137,147)),
                census=table,disjoint_obstructions=local,isolation=iso,
                reduced_profiles=rows,dimension_minima=mins,boundary_at_n136=boundary,
                compact_facet_ambient_bound=103,compact_facet_parity_slack_times8=[[12,208,4],[112,212,4]])

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=ROOT/f'results/{TAG}.json';d=run()
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'profiles',len(d['reduced_profiles']),'overlaps',d['census']['overlap_counts'],flush=True)
