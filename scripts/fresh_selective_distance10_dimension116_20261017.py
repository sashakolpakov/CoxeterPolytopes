#!/usr/bin/env python3
"""Exact distance10 joins/overlaps, controlled extensions and dimension116.

SymPy rational congruences supply exact inertia and nonnegative Z witnesses.
A frozen residual-tail lemma supplies the arbitrary-extension rank estimate.
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
TAG='fresh_selective_distance10_dimension116_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
RADIUS='fresh_simple_cusp_distance_bound_20261017'
BASE='fresh_selective_distance11_dimension126_20261017'
EXT='fresh_distance10_extension_20261017'
INPUTS=(CENSUS,RADIUS,BASE,BASE+'_review',EXT,
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
            if dist>=3:
                vals=values(r['diagram'],tuple(r['black']),node)
                ports.append(dict(diagram=r['diagram'],black=r['black'],node=node,distance=dist,
                    inverse_diagonals=list(map(str,vals)),allowed=bool(max(vals)<2)))
    assert sum(p['allowed']for p in ports)==161
    return diagrams,grams,contexts,ports,values


def edge_obstruction(g,black):
    a,scales=rationalize(g);white=[i for i in range(g.rows)if i not in black]
    for mark in black:
        ids=sorted(white+[mark]);ob=z_obstruction(a.extract(ids,ids))
        if ob:return dict(mark=mark,nodes=ids,gauge_squared=scales,**ob)
    return None


def disjoint(grams,ports):
    allowed=[p for p in ports if p['allowed']];rows=[];good={};types={}
    @lru_cache(None)
    def white_diagonal(di,black,node):
        g,sc=rationalize(grams[di]);ids=[i for i in range(g.rows)if i not in black]
        a=g.extract(ids,ids);D.pd(a)
        return sc[node]*a.inv()[ids.index(node),ids.index(node)]
    @lru_cache(None)
    def tail(di,black,node):
        g=grams[di];white=[i for i in range(g.rows)if i not in black]
        choices=[[node,x,y]for x,y in product(white,repeat=2)if len({node,x,y})==3 and g[node,x]and g[x,y]and not g[node,y]]
        assert choices;return min(choices)
    for left,right in product(allowed,repeat=2):
        if left['distance']+right['distance']<10:continue
        di,dj=left['diagram'],right['diagram'];u,v0=left['black'],right['black'];i,j=left['node'],right['node']
        core,mp=D.union(grams[di],grams[dj],());core[i,mp[j]]=core[mp[j],i]=-1
        white=[k for k in range(core.rows)if k not in u+[mp[k]for k in v0]]
        x=white_diagonal(di,tuple(u),i);y=white_diagonal(dj,tuple(v0),j)
        row=dict(left=[di,u,i],right=[dj,v0,j],white_inverse_product=str(x*y))
        if x*y>=1:
            a,sc=rationalize(core);o=z_obstruction(a.extract(white,white));assert o
            row['white_rejection']=dict(nodes=white,gauge_squared=sc,**o)
        else:
            key=(di,tuple(u),dj,tuple(v0));good.setdefault(key,[]).append((i,j))
            row['white_three_node_tails']=[tail(di,tuple(u),i),tail(dj,tuple(v0),j)]
            typ=tuple(sorted(((di,i),(dj,j))));row['inertia_type']=[list(p)for p in typ]
            if typ not in types:
                (a,b),(c,d)=typ;g,gm=D.union(grams[a],grams[c],());g[b,gm[d]]=g[gm[d],b]=-1
                h,sc=rationalize(g);cert=exact_inertia(h);assert cert['signature'][1]==2
                types[typ]=dict(ports=[list(p)for p in typ],gram=matrix(h),gauge_squared=sc,**cert)
        rows.append(row)
    multiple=[]
    for (di,u,dj,v0),edges in sorted(good.items()):
        if len(edges)<2:continue
        base,mp=D.union(grams[di],grams[dj],());black=list(u)+[mp[i]for i in v0]
        for size in range(2,len(edges)+1):
            for chosen in combinations(edges,size):
                g=base.copy()
                for x,y in chosen:g[x,mp[y]]=g[mp[y],x]=-1
                cert=edge_obstruction(g,black);assert cert
                multiple.append(dict(left=[di,list(u)],right=[dj,list(v0)],edges=list(map(list,chosen)),edge_pd_rejection=cert))
    assert len(rows)==1478 and sum('white_rejection'in r for r in rows)==363
    assert Counter(map(len,good.values()))=={1:905,2:96,3:6}
    assert len(types)==35 and len(multiple)==120
    return dict(candidate_joins=1478,white_rejections=363,remaining_joins=1115,
        remaining_context_pair_multiplicities=dict(Counter(map(len,good.values()))),rows=rows,
        underlying_inertia_types=[types[t]for t in sorted(types)],multiple_edge_obstructions=multiple,
        bond_forcing='Every allowed port has a wholly white three-node path. The positive coefficients(1,2/3,1/3) on each side give square<=8/3-2beta, including other nonpositive bonds. Hence every actual cross contact is simple. All nonempty subsets of the remaining at-most-three edges are then exhausted.')


def low_order_disjoint(grams,ports):
    """Complete order5/radius3 supplement, including its all-label proof.

    Radius3 gives white chain4-3-2 with black0/1. Connectivity permits
    three tree attachments; the only additional triangle is forbidden by
    a proper4-node positive-kernel argument. Every bond in a connected
    proper4-node PSD graph has label<=5: a label>=6 produces an affine or
    indefinite connected3-node block, and the fourth-node extension rules
    out the affine equality. Thus labels3..5 here are exhaustive.
    """
    theta={3:sp.Integer(1),4:sp.sqrt(2),5:(1+sp.sqrt(5))/2}
    shapes=[[(0,1),(1,2),(2,3),(3,4)],[(0,1),(0,2),(2,3),(3,4)],[(0,2),(1,2),(2,3),(3,4)]]
    @lru_cache(None)
    def determinant(entries,n):return sp.simplify(sp.Matrix(n,n,entries).det())
    def delta(g,ids):
        h=g.extract(ids,ids);return determinant(tuple(h),len(ids))
    records=[];survivors=[];counts=Counter()
    for si,edges in enumerate(shapes):
        for labels in product((3,4,5),repeat=4):
            g=2*sp.eye(5)
            for(i,j),m in zip(edges,labels):g[i,j]=g[j,i]=-theta[m]
            row=dict(shape=si,edges=[list(e)for e in edges],labels=list(labels));bad=None
            for size in range(3,5):
                for ids in combinations(range(5),size):
                    d=delta(g,ids)
                    if d<0:bad=dict(nodes=list(ids),determinant=str(d));break
                if bad:break
            if bad:row['proper_indefinite']=bad;counts['proper_indefinite']+=1
            elif delta(g,range(5))>=0:
                # All smaller principal minors have passed; nonnegative
                # full determinant makes the entire symmetric Gram PSD.
                row['nonnegative_full_determinant']=str(delta(g,range(5)));counts['not_indefinite']+=1
            else:
                values=[];blackgood=True
                for b in(0,1):
                    ids=[i for i in range(5)if i!=b]
                    if delta(g,ids)<=0:blackgood=False;break
                    h=g.extract(ids,ids);values.append(sp.simplify(h.inv()[ids.index(4),ids.index(4)]))
                if not blackgood:row['black_not_PD']=True;counts['black_not_PD']+=1
                else:
                    row['inverse_diagonals']=list(map(str,values))
                    if max(values)>=2:row['port_rejected']=True;counts['port_rejected']+=1
                    else:
                        counts['survive']+=1;row['gram']=matrix(g);survivors.append(row)
            records.append(row)
    assert counts==dict(proper_indefinite=209,not_indefinite=14,black_not_PD=4,port_rejected=12,survive=4)
    high=[r for r in ports if r['allowed']and r['distance']>=7]
    assert len(high)==3 and all(r['diagram']==25 and r['node']==9 for r in high)
    joins=[]
    for row in survivors:
        g=sp.Matrix(row['gram']);core,mp=D.union(g,grams[25],());core[4,mp[9]]=core[mp[9],4]=-1
        cert=exact_inertia(core);assert cert['signature']==[13,2,0]
        joins.append(dict(low_shape=row['shape'],low_labels=row['labels'],gram=matrix(core),
            high_contexts=[r['black']for r in high],ordered_marked_multiplicity=6,**cert))
    assert sum(r['ordered_marked_multiplicity']for r in joins)==24
    return dict(shapes=[list(map(list,s))for s in shapes],counts=dict(counts),records=records,
        surviving_low_cores=survivors,exact_join_types=joins,ordered_marked_joins=24,
        scope='Complete low-order supplement: orders<=4 have radius<=2; order5 radius3 has exactly the three stated tree supports. Every actual cross bond is the unique4-to9 edge and is simple by the same two-white-tail argument.')


def overlap_type(di,dj,match):
    if di>dj:return overlap_type(dj,di,tuple(sorted((j,i)for i,j in match)))
    if di==dj==25:
        return {((9,9),):'L25_vertex',((8,9),(9,8)):'L25_edge',
            ((7,9),(8,8),(9,7)):'L25_three_path',
            ((6,9),(7,8),(8,7),(9,6)):'L25_four_path'}[match]
    return {(10,25,((5,9),)):'L10_L25_vertex',(18,25,((7,9),)):'L18_L25_vertex',
        (22,25,((4,9),)):'L22_L25_short_vertex',(22,25,((8,9),)):'L22_L25_long_vertex',
        (22,25,((7,9),(8,8))):'L22_L25_edge'}[di,dj,match]


def overlaps(grams,contexts,values):
    high=[r for r in contexts if max(r['distances'])>=4]
    assert all(sum(d>=4 for d in r['distances'])<=4 for r in high)
    records=[];generated=0
    for left,right in product(high,repeat=2):
        g,h=grams[left['diagram']],grams[right['diagram']]
        cand=[(i,j)for i in range(g.rows)for j in range(h.rows)if left['distances'][i]+right['distances'][j]>=11]
        size_bound=min(sum(d>=4 for d in left['distances']),sum(d>=4 for d in right['distances']))
        for size in range(1,size_bound+1):
            for match in combinations(cand,size):
                generated+=1;ls=[i for i,j in match];rs=[j for i,j in match]
                if len(set(ls))<size or len(set(rs))<size:continue
                if any(g[i,k]!=h[j,l]for i,j in match for k,l in match):continue
                row=dict(left=[left['diagram'],left['black']],right=[right['diagram'],right['black']],match=list(map(list,match)))
                bad=None
                for i,j in match:
                    for side,x,y,rr,other,shared,other_black in [
                        ('left',i,j,left,h,rs,right['black']),('right',j,i,right,g,ls,left['black'])]:
                        ext=[k for k in range(other.rows)if k!=y and other[y,k]and k not in shared]
                        if ext and max(values(rr['diagram'],tuple(rr['black']),x))>=2:
                            assert all(k not in other_black for k in ext)
                            bad=dict(side=side,node=x,external_neighbor=ext[0],inverse_diagonals=list(map(str,values(rr['diagram'],tuple(rr['black']),x))))
                            break
                    if bad:break
                if bad:row['inverse_rejection']=bad
                else:
                    core,mp=D.union(g,h,match);black=left['black']+[mp[i]for i in right['black']]
                    ob=edge_obstruction(core,black)
                    if ob:row['edge_pd_rejection']=ob
                    else:row['surviving_type']=overlap_type(left['diagram'],right['diagram'],match)
                records.append(row)
    counts=Counter('inverse_rejection'if'inverse_rejection'in r else'edge_pd_rejection'if'edge_pd_rejection'in r else r['surviving_type']for r in records)
    assert generated==22881 and len(records)==4862
    assert counts==dict(inverse_rejection=4296,edge_pd_rejection=254,L10_L25_vertex=6,
        L18_L25_vertex=6,L22_L25_short_vertex=18,L22_L25_long_vertex=54,L22_L25_edge=18,
        L25_vertex=114,L25_edge=60,L25_three_path=27,L25_four_path=9)
    return dict(overlap_generated=generated,consistent_overlaps=len(records),overlap_counts=dict(counts),overlap_records=records)


def corner(g,sc,u,v,p,q):
    omitted=sorted(set(u+v)-{p,q});ids=[i for i in range(g.rows)if i not in omitted]
    h=g.extract(ids,ids);available=list(omitted);blocks=[];marked=None
    for part in D.components(h):
        actual=[ids[i]for i in part];a=h.extract(part,part);ks=a.nullspace()
        if ks:
            assert len(ks)==1;k=ks[0]
            if all(x<0 for x in k):k=-k
            assert all(x>0 for x in k);D.pd(a[:-1,:-1])
            blocks.append(dict(nodes=actual,kind='affine',positive_kernel=list(map(str,k))))
        else:
            D.pd(a);iv=a.inv();assert all(x>=0 for x in iv)
            vals=[sc[actual[i]]*iv[i,i]for i in range(a.rows)]
            available.extend(actual[i]for i in range(a.rows)if vals[i]<=2)
            blocks.append(dict(nodes=actual,kind='finite',normalized_inverse_diagonals=list(map(str,vals))))
        if p in actual and q in actual:marked=actual
    if marked is None:return None
    return dict(preferred_corner=[p,q],omitted_black=omitted,blocks=blocks,
        possible_contacts=sorted(available),marked_component=marked)


def extra_cross(g,mp,match,u,v,di,dj,grams):
    left=[i for i in range(grams[di].rows)if i not in{i for i,j in match}]
    right=[mp[j]for j in range(grams[dj].rows)if j not in{j for i,j in match}]
    short=0;obs=[]
    for x,y in product(left,right):
        if g[x,y]:continue
        h=g.copy();h[x,y]=h[y,x]=-1
        if min(P.distances(h,u)[j]for j in v)<=10:short+=1;continue
        o=edge_obstruction(h,u+v);assert o
        obs.append(dict(edge=[x,y],**o))
    return dict(shortcut_rejections=short,edge_pd_z_rejections=obs)


def isolation(grams,table):
    cases=[];cores={};short=extra=0
    for row in table['overlap_records']:
        typ=row.get('surviving_type')
        if typ is None:continue
        di,u=row['left'];dj,v0=row['right'];match=tuple(map(tuple,row['match']))
        native,mp=D.union(grams[di],grams[dj],match);v=[mp[i]for i in v0]
        assert min(P.distances(native,u)[j]for j in v)>10
        cross=extra_cross(native,mp,match,u,v,di,dj,grams)
        short+=cross['shortcut_rejections'];extra+=len(cross['edge_pd_z_rejections'])
        g,sc=rationalize(native);key=str((di,dj,match));ks=g.nullspace()
        if key not in cores:
            inert=exact_inertia(g);assert inert['signature'][1]==1
            core=dict(left_diagram=di,right_diagram=dj,match=list(map(list,match)),gram=matrix(g),gauge_squared=sc,**inert)
            if ks:
                assert len(ks)==1;core['radical']=list(map(str,ks[0]))
            else:
                assert set(sc)=={1};core['inverse']=matrix(g.inv())
            cores[key]=core
        du,dv=P.distances(native,u),P.distances(native,v);choices=[]
        for p,q in product(u,v):
            c=corner(g,sc,u,v,p,q)
            if not c:continue
            allowed=c['possible_contacts']
            if ks:
                z=ks[0]
                if any(z[i]==0 for i in allowed):continue
                pos=[i for i in allowed if z[i]>0];neg=[i for i in allowed if z[i]<0]
                paths=[max([du[i]for i in pos],default=0)+max([dv[i]for i in neg],default=0)+2,
                       max([dv[i]for i in pos],default=0)+max([du[i]for i in neg],default=0)+2]
                if min(paths)>10:continue
                c.update(method='radical',shortcut_bound=min(paths),component_order_bound=len(c['marked_component']))
            else:
                iv=g.inv();left=[i for i in allowed if du[i]<=4];right=[i for i in allowed if dv[i]<=4]
                if set(left+right)!=set(allowed)or set(left)&set(right):continue
                if not all(iv[i,j]>0 for block in(left,right)for i in block for j in block):continue
                # A single external root contacting both ends already gives
                # a path of length<=10 through this root.
                rootnodes=[i for i in allowed if iv[i,i]<=2]
                if any(iv[i,i]<=1 for i in rootnodes):continue
                if any(iv[i,i]+2*iv[i,j]+iv[j,j]<=2 for block in(left,right)for i,j in combinations(block,2)):continue
                pairtests=[];valid=True
                for i,j in product(rootnodes,repeat=2):
                    a=iv[i,j];rhoi=2-iv[i,i];rhoj=2-iv[j,j]
                    if not a*a>rhoi*rhoj:valid=False;break
                    shortcut=min(du[i]+dv[j],du[j]+dv[i])+3
                    if a<0 and shortcut>10:valid=False;break
                    pairtests.append(dict(nodes=[i,j],projection_product=str(a),residual_product=str(rhoi*rhoj),
                        reason='positive_product_forbidden'if a>0 else'forced_cross_contact',shortcut_bound=shortcut if a<0 else None))
                if not valid:continue
                extensions=[]
                for i in rootnodes:
                    rho=2-iv[i,i];assert rho in(0,sp.Rational(1,2),sp.Rational(2,3))
                    tail=0 if rho==0 else 1 if rho==sp.Rational(1,2)else 2
                    extensions.append(dict(contact_node=i,contact='-1',projection_square=str(iv[i,i]),
                        residual_square=str(rho),maximum_attached_orthogonal_rank=tail))
                bound=len(c['marked_component'])+max([1+e['maximum_attached_orthogonal_rank']for e in extensions],default=0)
                if bound>17:continue
                c.update(method='lorentz_residual',end_sets=[left,right],single_contact_candidates=extensions,
                    candidate_pair_tests=pairtests,component_order_bound=bound)
            if c['component_order_bound']<=17:choices.append(c)
        assert choices,(di,u,dj,v0,match)
        chosen=min(choices,key=lambda c:(len(c.get('single_contact_candidates',[])),len(c['marked_component']),c['preferred_corner']))
        # Finite/affine decompositions for all four corners have already been
        # checked by corner(), including ones whose marked nodes disconnect.
        cases.append(dict(type=typ,left=[di,u],right=[dj,v0],right_copy_map=mp,match=list(map(list,match)),
            core_key=key,left_black=u,right_black=v,extra_cross_bond_checks=cross,**chosen))
    assert len(cases)==312 and short==21393 and extra==168 and len(cores)==14
    assert max(r['component_order_bound']for r in cases)==17
    return dict(cores=cores,cases=cases,component_order_cap=17,extra_cross_shortcuts=short,extra_cross_pd_rejections=extra,
        residual_tail_lemma='All other supporting whites are core-orthogonal and have PD Gram E. Projection costs c^T E^-1 c<=rho. Two contacted components cost>=1; a connected rank>=3 component costs>=3/4; rank>=2 costs>=2/3. Hence rho0,1/2,2/3 permit at most0,1,2 attached roots. All candidate core contacts are simple, and pair tests force at most one such root.')


def profiles(d):
    q,a=divmod(d,20);yield(20,)*q+((a,)if a else())
    for q in range(d//20+1):
        for a in range(20):
            b=d-20*q-a
            if b>=21:yield(20,)*q+((a,)if a else())+(b,)


def budget(rs):
    d=sum(rs);v=prod(r+1 for r in rs);h=P.polynomial(rs)
    assert sum(h)==v and h==h[::-1]
    G=sum(x*min(i,d-i)for i,x in enumerate(h));curvature=F(5*G-d*v,2*v)
    assert curvature==F(3*d,4)-F(5,2)*P.independent_moment(rs)
    cost=sum((min(F(10),F(r,2))for r in rs),F())
    return dict(ranks=list(rs),link_dimension=d,link_vertices=v,G=G,
        normalized_curvature=str(curvature),normalized_cost=str(cost),margin=str(curvature-cost))


def run():
    data={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in INPUTS}
    for n,d in data.items():
        assert d['status'].startswith('PASS'),n
        for path,h in d.get('sources',{}).items():assert digest(ROOT/path)==h,(n,path)
    diagrams,grams,contexts,ports,values=contexts_and_ports(data);assert len(ports)==527
    joins=disjoint(grams,ports);low=low_order_disjoint(grams,ports);print('DISJOINT_PASS',1478,363,120,'LOW5',24,flush=True)
    table=overlaps(grams,contexts,values);print('OVERLAP_PASS',table['overlap_counts'],flush=True)
    iso=isolation(grams,table);print('ISOLATION_PASS',len(iso['cases']),flush=True)
    rows=[];mins=[]
    for d in range(116,126):
        rr=list(profiles(d));assert len(rr)==len(set(tuple(sorted(rs))for rs in rr))
        group=[budget(rs)for rs in rr];assert all(F(r['margin'])>0 for r in group)
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8);assert beta>10*n
        mins.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=str(min(F(r['margin'])for r in group)),finite_slack=str(beta-10*n)));rows.extend(group)
    assert len(rows)==1015 and min(F(r['margin'])for r in rows)==F(145234,4084101)
    boundary=budget((20,)*5+(15,));assert F(boundary['margin'])==F(-673055,16336404)
    t=sp.Symbol('t');facet=[]
    for start in(87,88):
        d=start+2*t;beta=(d-1)*(d-5)/8 if start%2 else d*(d-6)/8
        poly=sp.Poly(sp.expand(8*(beta-10*(d+1))),t);facet.append([int(poly.nth(k))for k in range(3)])
    assert facet==[[12,176,4],[96,180,4]]
    paths=[Path(__file__),Path(D.__file__)]+[ROOT/f'{folder}/{n}.{ext}'for n in INPUTS for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_SELECTIVE_DISTANCE10_CONTROLLED_EXTENSIONS_DIMENSION116',sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_bound=116,inherited_dimension_bound=126,unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        excluded_dimensions=list(range(117,127)),weight_rule='1 if marked distance<=10 OR both marks lie in one component of order<=17;0 otherwise',
        finite_budget='10n',ideal_budget='sum min(10,r_i/2)',census=dict(marked_contexts=len(contexts),ports=ports,allowed_ports=161,**table),
        disjoint_obstructions=joins,low_order_disjoint=low,isolation=iso,reduced_profiles=rows,dimension_minima=mins,boundary_at_n116=boundary,
        compact_facet_ambient_bound=87,compact_facet_parity_slack_times8=facet,
        scope='Exact finite classification, residual budgets and scalar inequalities with explicit inherited geometry. No construction, external refereeing, formalization or priority claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'1015profiles/4862overlaps/312corners',flush=True)
