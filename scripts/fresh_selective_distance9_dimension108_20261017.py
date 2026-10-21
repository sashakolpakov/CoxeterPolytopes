#!/usr/bin/env python3
"""Exact distance9 joins/overlaps, controlled extensions and dimension108.

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
TAG='fresh_selective_distance9_dimension108_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
RADIUS='fresh_simple_cusp_distance_bound_20261017'
BASE='fresh_selective_distance10_dimension116_20261017'
EXT='fresh_distance9_low_core_20261017'
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
    grams={i:P.gram(d)for i,d in diagrams.items()};contexts=[]
    for i,d in sorted(diagrams.items()):
        for black in combinations(d['finite_deletions'],2):
            contexts.append(dict(diagram=i,black=list(black),distances=P.distances(grams[i],black)))
    assert len(contexts)==384 and max(max(r['distances'])for r in contexts)==7
    @lru_cache(None)
    def deletion_inverse(di,b):
        ids=[j for j in range(grams[di].rows)if j!=b]
        h,scales=rationalize(grams[di].extract(ids,ids));D.pd(h)
        return ids,scales,h.inv()
    @lru_cache(None)
    def values(di,black,node):
        out=[]
        for b in black:
            ids,scales,j=deletion_inverse(di,b);k=ids.index(node)
            out.append(scales[k]*j[k,k])
        return tuple(out)
    ports=[]
    for r in contexts:
        for node,dist in enumerate(r['distances']):
            if dist>=2:
                vals=values(r['diagram'],tuple(r['black']),node)
                ports.append(dict(diagram=r['diagram'],black=r['black'],node=node,distance=dist,
                    inverse_diagonals=list(map(str,vals)),allowed=bool(max(vals)<2)))
    assert len(ports)==1262 and sum(p['allowed']for p in ports)==378
    return diagrams,grams,contexts,ports,values

def edge_obstruction(g,black):
    white=[i for i in range(g.rows)if i not in black]
    for mark in black:
        ids=sorted(white+[mark]);native=g.extract(ids,ids);mode='positive_diagonal_congruence'
        try:a,scales=rationalize(native)
        except AssertionError:
            # An odd label-four cycle prevents a rational gauge. Replacing
            # every nonright bond by1 is a safe Z-matrix upper comparison.
            a=sp.Matrix(native.rows,native.cols,lambda i,j:2 if i==j else -int(bool(native[i,j])))
            scales=[1]*native.rows;mode='unit_bond_relaxation'
        ob=z_obstruction(a)
        if ob:return dict(mark=mark,nodes=ids,gauge_squared=scales,comparison=mode,**ob)
    return None

def disjoint(grams,ports):
    allowed=[p for p in ports if p['allowed']];rows=[];good={};types={}
    @lru_cache(None)
    def white_diagonal(di,black,node):
        ids=[i for i in range(grams[di].rows)if i not in black]
        a,sc=rationalize(grams[di].extract(ids,ids));D.pd(a);k=ids.index(node)
        return sc[k]*a.inv()[k,k]
    for left,right in product(allowed,repeat=2):
        if left['distance']+right['distance']<9:continue
        di,dj=left['diagram'],right['diagram'];u,v0=left['black'],right['black'];i,j=left['node'],right['node']
        core,mp=D.union(grams[di],grams[dj],());core[i,mp[j]]=core[mp[j],i]=-1
        white=[k for k in range(core.rows)if k not in u+[mp[k]for k in v0]]
        x=white_diagonal(di,tuple(u),i);y=white_diagonal(dj,tuple(v0),j)
        row=dict(left=[di,u,i],right=[dj,v0,j],white_inverse_product=str(x*y))
        if x*y>=1:
            a,sc=rationalize(core.extract(white,white));o=z_obstruction(a);assert o
            row['white_rejection']=dict(nodes=white,gauge_squared=sc,**o)
        else:
            # All labels>=4 have squared magnitude>=2. The inverse test
            # or the supporting-white Schur determinant excludes them.
            maxima=[max(map(sp.Rational,left['inverse_diagonals'])),max(map(sp.Rational,right['inverse_diagonals']))]
            assert max(maxima)>=1 or 2*x*y>=1
            row['larger_bond_obstruction']=dict(endpoint_maxima=list(map(str,maxima)),twice_white_inverse_product=str(2*x*y))
            key=(di,tuple(u),dj,tuple(v0));good.setdefault(key,[]).append((i,j))
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
    assert len(rows)==4442 and sum('white_rejection'in r for r in rows)==1103
    assert Counter(map(len,good.values()))=={1:2007,2:348,3:180,4:24}
    assert len(types)==85 and len(multiple)==1332
    return dict(candidate_joins=len(rows),white_rejections=1103,remaining_joins=3339,
        remaining_context_pair_multiplicities=dict(Counter(map(len,good.values()))),rows=rows,
        underlying_inertia_types=[types[t]for t in sorted(types)],multiple_edge_obstructions=multiple,
        bond_forcing='Every label>=4 is excluded by the actual-edge normalized inverse test or the supporting-white Schur product. All nonempty subsets of up tofour remaining simple cross edges are exhausted.')

def overlap_type(di,dj,match):
    return str((di,dj,match))

def overlaps(grams,contexts,values):
    high=[r for r in contexts if max(r['distances'])>=3]
    assert all(sum(d>=3 for d in r['distances'])<=5 for r in high)
    records=[];generated=0
    for left,right in product(high,repeat=2):
        g,h=grams[left['diagram']],grams[right['diagram']]
        cand=[(i,j)for i in range(g.rows)for j in range(h.rows)if left['distances'][i]+right['distances'][j]>=10]
        size_bound=min(sum(d>=3 for d in left['distances']),sum(d>=3 for d in right['distances']))
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
    assert generated==312316 and len(records)==20718
    assert counts['inverse_rejection']==18721 and counts['edge_pd_rejection']==1048
    assert sum(v for k,v in counts.items()if k not in('inverse_rejection','edge_pd_rejection'))==949
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
        if min(P.distances(h,u)[j]for j in v)<=9:short+=1;continue
        o=edge_obstruction(h,u+v);assert o
        obs.append(dict(edge=[x,y],**o))
    return dict(shortcut_rejections=short,edge_pd_z_rejections=obs)


def ordinary_corner(g,sc,iv,ks,c,du,dv):
    allowed=c['possible_contacts'];component=len(c['marked_component'])
    if ks:
        z=ks[0]
        if any(z[i]==0 for i in allowed):return None
        pos=[i for i in allowed if z[i]>0];neg=[i for i in allowed if z[i]<0]
        bounds=[max([du[i]for i in pos],default=0)+max([dv[i]for i in neg],default=0)+2,
                max([dv[i]for i in pos],default=0)+max([du[i]for i in neg],default=0)+2]
        if min(bounds)>9:return None
        return dict(**c,method='radical',shortcut_bound=min(bounds),component_order_bound=component)
    options=[]
    for radius in range(8):
        left=[i for i in allowed if du[i]<=radius];right=[i for i in allowed if dv[i]<=7-radius]
        if set(left+right)==set(allowed)and not(set(left)&set(right))and all(iv[i,j]>0 for part in(left,right)for i in part for j in part):options.append([left,right])
    if not options:return None
    ends=options[0];rootnodes=[i for i in allowed if sc[i]*iv[i,i]<=2]
    if any(sc[i]*iv[i,i]<=1 for i in rootnodes):return None
    if any(sc[i]*iv[i,i]+2*min(sc[i],sc[j])*iv[i,j]+sc[j]*iv[j,j]<=2 for part in ends for i,j in combinations(part,2)):return None
    tests=[]
    for i,j in product(rootnodes,repeat=2):
        a=iv[i,j];rhoi=2-sc[i]*iv[i,i];rhoj=2-sc[j]*iv[j,j];square=sc[i]*sc[j]*a*a
        shortcut=min(du[i]+dv[j],du[j]+dv[i])+3
        if square<=rhoi*rhoj or(a<0 and shortcut>9):return None
        tests.append(dict(nodes=[i,j],gauge_projection_product=str(a),normalized_projection_square=str(square),
            residual_product=str(rhoi*rhoj),reason='positive_product_forbidden'if a>0 else'forced_cross_contact',
            shortcut_bound=shortcut if a<0 else None))
    extensions=[]
    for i in rootnodes:
        rho=2-sc[i]*iv[i,i];assert 0<=rho<1;tail=int(sp.floor(rho/(1-rho)))
        extensions.append(dict(contact_node=i,contact='-1',projection_square=str(sc[i]*iv[i,i]),
            residual_square=str(rho),maximum_attached_orthogonal_rank=tail))
    bound=component+max([1+e['maximum_attached_orthogonal_rank']for e in extensions],default=0)
    if bound>17:return None
    return dict(**c,method='lorentz_residual',end_sets=ends,single_contact_candidates=extensions,
                candidate_pair_tests=tests,component_order_bound=bound)


def relevant_corner(g,sc,iv,c,du,dv):
    assert set(sc)=={1}
    allowed=c['possible_contacts'];ends=None
    for radius in range(8):
        left=[i for i in allowed if du[i]<=radius];right=[i for i in allowed if dv[i]<=7-radius]
        if set(left+right)==set(allowed)and not(set(left)&set(right))and all(iv[i,j]>0 for end in(left,right)for i in end for j in end):ends=[left,right];break
    assert ends
    assert all(iv[i,i]+2*iv[i,j]+iv[j,j]>2 for end in ends for i,j in combinations(end,2))
    candidates=[i for i in allowed if iv[i,i]<=2];assert all(iv[i,i]>=1 for i in candidates)
    relevant=[i for i in candidates if i in c['marked_component']]
    if any(iv[i,i]<=1 for i in relevant):return None
    proofs=[]
    for i in relevant:
        rho=2-iv[i,i];live=[];pair=[]
        for j in candidates:
            residual=2-iv[j,j];q=iv[i,j];distance=min(du[i]+dv[j],du[j]+dv[i])+3
            if q>0 and q*q>rho*residual:pair.append(dict(node=j,kind='positive'));continue
            if q<0 and q*q>rho*residual and distance<=9:pair.append(dict(node=j,kind='short_path'));continue
            if q>0 and rho*residual<(q+1)**2:
                if j in relevant:return None
                live.append(j);pair.append(dict(node=j,kind='orthogonal_by_gap',projection=str(q),residual_product=str(rho*residual)))
            else:return None
        label4=[]
        for j in candidates:
            if iv[j,j]!=1:continue
            q=iv[i,j];distance=min(du[i]+dv[j],du[j]+dv[i])+3
            if not(q>0 or(q<0 and distance<=9)):return None
            label4.append(dict(node=j,reason='positive'if q>0 else'short_path'))
        if len(live)>1:return None
        assert 0<=rho<1;alone_tail=int(sp.floor(rho/(1-rho)));joint=None
        if live:
            j=live[0];residual=2-iv[j,j];q=iv[i,j];delta=rho*residual-q*q
            assert delta>=0 and residual<=1 and q>0
            if delta==0:joint=dict(outside_node=j,determinant='0',kind='positive_residual_relation',tail=0)
            else:
                inverseGG=residual/delta
                if inverseGG<2:return None
                joint=dict(outside_node=j,determinant=str(delta),inverseGG=str(inverseGG),
                    kind='no_tail'if inverseGG>2 else'saturated_single_tail',tail=0 if inverseGG>2 else 1)
        proofs.append(dict(node=i,projection=str(iv[i,i]),rho=str(rho),alone_tail_bound=alone_tail,
            pair_exclusions=pair,label4_exclusions=label4,coexisting_outside=joint))
    bound=len(c['marked_component'])+max([1+max(r['alone_tail_bound'],(r['coexisting_outside']or{}).get('tail',0))for r in proofs],default=0)
    assert bound<=17
    return dict(**c,method='relevant_component',relevant_nodes=relevant,contact_domains=candidates,end_sets=ends,
        extension_cases=proofs,component_order_bound=bound)


def isolation(grams,table):
    cases=[];cores={};short=extra=0
    for row in table['overlap_records']:
        typ=row.get('surviving_type')
        if typ is None:continue
        di,u=row['left'];dj,v0=row['right'];match=tuple(map(tuple,row['match']))
        native,mp=D.union(grams[di],grams[dj],match);v=[mp[i]for i in v0]
        assert min(P.distances(native,u)[j]for j in v)>9
        cross=extra_cross(native,mp,match,u,v,di,dj,grams)
        short+=cross['shortcut_rejections'];extra+=len(cross['edge_pd_z_rejections'])
        g,sc=rationalize(native);key=str((di,dj,match));ks=g.nullspace();iv=None if ks else g.inv()
        if key not in cores:
            inert=exact_inertia(g);assert inert['signature'][1]==1
            core=dict(left_diagram=di,right_diagram=dj,match=list(map(list,match)),gram=matrix(g),gauge_squared=sc,**inert)
            if ks:assert len(ks)==1;core['radical']=list(map(str,ks[0]))
            else:core['inverse']=matrix(iv)
            cores[key]=core
        du,dv=P.distances(native,u),P.distances(native,v);corners=[];choices=[]
        for p,q in product(u,v):
            c=corner(g,sc,u,v,p,q)
            if c:
                corners.append(c);a=ordinary_corner(g,sc,iv,ks,c,du,dv)
                if a:choices.append(a)
        if not choices:
            assert not ks
            for c in corners:
                a=relevant_corner(g,sc,iv,c,du,dv)
                if a:choices.append(a)
        assert choices,(di,u,dj,v0,match)
        if choices[0]['method']=='relevant_component':chosen=min(choices,key=lambda c:(c['component_order_bound'],len(c['extension_cases']),c['preferred_corner']))
        else:chosen=min(choices,key=lambda c:(c['component_order_bound'],len(c['possible_contacts'])))
        cases.append(dict(type=typ,left=[di,u],right=[dj,v0],right_copy_map=mp,match=list(map(list,match)),
            core_key=key,left_black=u,right_black=v,extra_cross_bond_checks=cross,**chosen))
    assert len(cases)==949 and short==60831 and extra==888 and len(cores)==54
    assert sum(r['method']=='relevant_component'for r in cases)==69
    assert max(r['component_order_bound']for r in cases)==17
    return dict(cores=cores,cases=cases,component_order_cap=17,extra_cross_shortcuts=short,extra_cross_pd_rejections=extra,
        residual_tail_lemma='Every connected PD norm2 Coxeter Gram of rankm has inverse diagonal>=m/(m+1). Across all contacted components this gives totalrank<=floor(rho/(1-rho)) forrho<1. When an outside core-contacting root coexists, the exact two-residual Schur tests exclude or seal every tail leading from the preferred marked component.')

def profiles(d):
    q,a=divmod(d,18);yield(18,)*q+((a,)if a else())
    for q in range(d//18+1):
        for a in range(18):
            b=d-18*q-a
            if b>=19:yield(18,)*q+((a,)if a else())+(b,)


def budget(rs):
    d=sum(rs);v=prod(r+1 for r in rs);h=P.polynomial(rs)
    assert sum(h)==v and h==h[::-1]
    G=sum(x*min(i,d-i)for i,x in enumerate(h));curvature=F(5*G-d*v,2*v)
    assert curvature==F(3*d,4)-F(5,2)*P.independent_moment(rs)
    cost=sum((min(F(9),F(r,2))for r in rs),F())
    return dict(ranks=list(rs),link_dimension=d,link_vertices=v,G=G,
        normalized_curvature=str(curvature),normalized_cost=str(cost),margin=str(curvature-cost))


def run():
    data={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in INPUTS}
    for n,d in data.items():
        assert d['status'].startswith('PASS'),n
        for path,h in d.get('sources',{}).items():assert digest(ROOT/path)==h,(n,path)
    diagrams,grams,contexts,ports,values=contexts_and_ports(data)
    joins=disjoint(grams,ports);print('HIGH_DISJOINT_PASS',4442,1103,1332,flush=True)
    table=overlaps(grams,contexts,values);print('OVERLAP_PASS',table['consistent_overlaps'],flush=True)
    iso=isolation(grams,table);print('ISOLATION_PASS',len(iso['cases']),flush=True)
    rows=[];mins=[]
    for d in range(108,116):
        rr=list(profiles(d));assert len(rr)==len(set(tuple(sorted(rs))for rs in rr))
        group=[budget(rs)for rs in rr];assert all(F(r['margin'])>0 for r in group)
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8);assert beta>9*n
        mins.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=str(min(F(r['margin'])for r in group)),finite_slack=str(beta-9*n)));rows.extend(group)
    assert len(rows)==756 and min(F(r['margin'])for r in rows)==F(45873,2476099)
    boundary=budget((18,)*5+(17,));assert F(boundary['margin'])==F(-4659,260642)
    t=sp.Symbol('t');facet=[]
    for start in(79,80):
        d=start+2*t;beta=(d-1)*(d-5)/8 if start%2 else d*(d-6)/8
        poly=sp.Poly(sp.expand(8*(beta-9*(d+1))),t);facet.append([int(poly.nth(k))for k in range(3)])
    assert facet==[[12,160,4],[88,164,4]]
    paths=[Path(__file__),Path(D.__file__)]+[ROOT/f'{folder}/{n}.{ext}'for n in INPUTS for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_SELECTIVE_DISTANCE9_CONTROLLED_EXTENSIONS_DIMENSION108',sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_bound=108,inherited_dimension_bound=116,unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        excluded_dimensions=list(range(109,117)),weight_rule='1 if marked distance<=9 OR both marks lie in one component of order<=17;0 otherwise',
        finite_budget='9n',ideal_budget='sum min(9,r_i/2)',census=dict(marked_contexts=len(contexts),ports=ports,allowed_ports=378,**table),
        disjoint_obstructions=joins,low_core_lemma=EXT,isolation=iso,reduced_profiles=rows,dimension_minima=mins,boundary_at_n108=boundary,
        compact_facet_ambient_bound=79,compact_facet_parity_slack_times8=facet,
        scope='Exact finite classification, controlled relevant-component extensions and scalar inequalities with explicit inherited geometry and complete small-order lemma. No construction, external refereeing, formalization or priority claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'756profiles/20718overlaps/949corners',flush=True)
