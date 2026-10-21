#!/usr/bin/env python3
"""Independent Fraction-only review of selective-distance10 dimension116.

Imports no producer arithmetic, SymPy, NumPy, or optimization package.
"""
import argparse
from collections import deque, Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product, permutations
from math import comb, prod
from functools import lru_cache
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_selective_distance10_dimension116_20261017_review'
BASE='fresh_selective_distance10_dimension116_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'


def sha(p):
    h=sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def matrix(n):return [[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
def extract(a,rows,cols=None):
    cols=rows if cols is None else cols
    return [[a[i][j]for j in cols]for i in rows]
def transpose(a):return list(map(list,zip(*a)))
def multiply(a,b):return [[sum((x*y for x,y in zip(r,c)),F())for c in zip(*b)]for r in a]
def square(a,x):return sum((x[i]*a[i][j]*x[j]for i in range(len(x))for j in range(len(x))),F())
def encode(a):return [[str(x)for x in r]for r in a]


def reduced(a):
    a=[list(map(F,r))for r in a];pivots=[];i=0
    for j in range(len(a[0])):
        k=next((k for k in range(i,len(a))if a[k][j]),None)
        if k is None:continue
        a[i],a[k]=a[k],a[i];x=a[i][j];a[i]=[v/x for v in a[i]]
        for k in range(len(a)):
            if k!=i:
                x=a[k][j];a[k]=[v-x*w for v,w in zip(a[k],a[i])]
        pivots.append(j);i+=1
        if i==len(a):break
    return a,pivots


def inverse(a):
    n=len(a)
    b,p=reduced([r+[F(i==j)for j in range(n)]for i,r in enumerate(a)])
    assert p[:n]==list(range(n))
    return [r[n:]for r in b]


def nullspace(a):
    b,p=reduced(a);n=len(a[0]);free=[i for i in range(n)if i not in p];out=[]
    for j in free:
        x=[F(0)]*n;x[j]=F(1)
        for k,i in enumerate(p):x[i]=-b[k][j]
        out.append(x)
    return out


def psd(a,strict=False):
    a=[r[:]for r in a];p=[]
    while a:
        q=a[0][0];assert q>=0
        if q==0:
            assert not strict and all(v==0 for v in a[0]);p.append(q)
            a=[r[1:]for r in a[1:]]
        else:
            p.append(q)
            a=[[a[i][j]-a[i][0]*a[0][j]/q for j in range(1,len(a))]for i in range(1,len(a))]
    return p


def graph_dist(a,marks,connected=True):
    out={i:0 for i in marks};todo=deque(marks)
    while todo:
        i=todo.popleft()
        for j,v in enumerate(a[i]):
            if j!=i and v and j not in out:out[j]=out[i]+1;todo.append(j)
    if connected:assert len(out)==len(a)
    return [out.get(i)for i in range(len(a))]


def rational_gram(d):
    n=d['order'];adj=[[]for _ in range(n)]
    for i,j,m in d['edges']:
        assert m in(3,4);adj[i].append((j,m==4));adj[j].append((i,m==4))
    bits={}
    for root in range(n):
        if root in bits:continue
        bits[root]=0;todo=deque([root])
        while todo:
            i=todo.popleft()
            for j,flip in adj[i]:
                wanted=bits[i]^flip
                if j in bits:assert bits[j]==wanted
                else:bits[j]=wanted;todo.append(j)
    assert len(bits)==n
    norms=[1<<bits[i]for i in range(n)]
    a=[[F(2*norms[i]if i==j else 0)for j in range(n)]for i in range(n)]
    for i,j,m in d['edges']:
        a[i][j]=a[j][i]=-F(2 if m==4 or bits[i]else 1)
    return a,norms


def bond_graph(d):
    a=matrix(d['order'])
    for i,j,m in d['edges']:a[i][j]=a[j][i]=F(-(m-2))
    return a


def all_contexts(diagrams):
    out=[]
    for d in diagrams:
        g=bond_graph(d)
        for black in combinations(d['finite_deletions'],2):
            out.append(dict(diagram=d['index'],black=tuple(black),distances=graph_dist(g,black)))
    assert len(out)==384
    return out


def overlap_key(row):
    return (row['left'][0],tuple(row['left'][1]),row['right'][0],tuple(row['right'][1]),
            tuple(tuple(x)for x in row['match']))


def unite(a,b,match):
    mapping={j:i for i,j in match};n=len(a)
    for j in range(len(b)):
        if j not in mapping:mapping[j]=n;n+=1
    out=matrix(n)
    for i in range(len(a)):
        for j in range(len(a)):out[i][j]=a[i][j]
    for i in range(len(b)):
        for j in range(len(b)):
            if mapping[i]<len(a)and mapping[j]<len(a):assert out[mapping[i]][mapping[j]]==b[i][j]
            out[mapping[i]][mapping[j]]=b[i][j]
    return out,[mapping[j]for j in range(len(b))]


def positive_kernel(a):
    psd(a);ks=nullspace(a);assert len(ks)==1
    k=ks[0]
    if all(x<0 for x in k):k=[-x for x in k]
    assert all(x>0 for x in k)
    return k



def obstruction(a):
    """First nonpositive Z-matrix pivot, lifted to a nonnegative witness."""
    b=[r.copy()for r in a]
    assert all(a[i][j]<=0 for i in range(len(a))for j in range(i))
    for k in range(len(a)):
        q=b[k][k]
        if q<=0:
            v=[F(0)]*len(a);v[k]=1
            if k:
                inv=inverse(extract(a,list(range(k))))
                v[:k]=[-sum((inv[i][j]*a[j][k]for j in range(k)),F())for i in range(k)]
            assert all(x>=0 for x in v)and square(a,v)==q
            return dict(vector=list(map(str,v)),square=str(q),pivot=k)
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                b[i][j]-=b[i][k]*b[k][j]/q
    return None


def inertia(a):
    a=[r.copy()for r in a];positive=negative=zero=0
    while a:
        k=next((i for i in range(len(a))if a[i][i]),None)
        if k is not None:
            ids=[k]+[i for i in range(len(a))if i!=k];a=extract(a,ids);v=a[0][0]
            positive+=int(v>0);negative+=int(v<0)
            a=[[a[i][j]-a[i][0]*a[0][j]/v for j in range(1,len(a))]for i in range(1,len(a))]
        else:
            pair=next(((i,j)for i in range(len(a))for j in range(i+1,len(a))if a[i][j]),None)
            if pair is None:zero+=len(a);break
            i,j=pair;ids=[i,j]+[k for k in range(len(a))if k not in pair];a=extract(a,ids);v=a[0][1]
            positive+=1;negative+=1
            a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/v for j in range(2,len(a))]for i in range(2,len(a))]
    return [positive,negative,zero]


def joined(d,e,p,q):
    off=d['order'];edges=d['edges']+[[i+off,j+off,m]for i,j,m in e['edges']]+[[p,off+q,3]]
    return rational_gram(dict(order=off+e['order'],edges=edges))


def native_to_rational(g):
    edges=[[i,j,int(2-g[i][j])]for i in range(len(g))for j in range(i+1,len(g))if g[i][j]]
    return rational_gram(dict(order=len(g),edges=edges))


def combined(ds,di,dj,match):
    native,mp=unite(bond_graph(ds[di]),bond_graph(ds[dj]),match)
    g,sc=native_to_rational(native)
    return native,g,sc,mp


def edge_bad(g,black):
    white=[i for i in range(len(g))if i not in black]
    for p in black:
        ids=sorted(white+[p]);o=obstruction(extract(g,ids))
        if o:return dict(mark=p,nodes=ids,**o)
    return None


def classify(diagrams):
    ds={d['index']:d for d in diagrams};graphs={i:bond_graph(d)for i,d in ds.items()};contexts=all_contexts(diagrams)
    assert max(max(r['distances'])for r in contexts)==7
    @lru_cache(None)
    def diagonals(di,black,node):
        a,scale=rational_gram(ds[di]);out=[]
        for b in black:
            ids=[i for i in range(len(a))if i!=b];h=extract(a,ids);psd(h,True)
            j=inverse(h);out.append(scale[node]*j[ids.index(node)][ids.index(node)])
        return tuple(out)
    @lru_cache(None)
    def white_diagonal(di,black,node):
        a,scale=rational_gram(ds[di]);ids=[i for i in range(len(a))if i not in black]
        h=extract(a,ids);psd(h,True);iv=inverse(h)
        return scale[node]*iv[ids.index(node)][ids.index(node)]
    ports=[]
    for r in contexts:
        for node,distance in enumerate(r['distances']):
            if distance<3:continue
            vals=diagonals(r['diagram'],r['black'],node)
            ports.append(dict(diagram=r['diagram'],black=list(r['black']),node=node,distance=distance,
                              inverse_diagonals=list(map(str,vals)),allowed=max(vals)<2))
    allowed=[p for p in ports if p['allowed']];assert len(ports)==527 and len(allowed)==161
    joins=[];good={};types={}
    for a,b in product(allowed,repeat=2):
        if a['distance']+b['distance']<10:continue
        di,dj=a['diagram'],b['diagram'];u,v=a['black'],b['black'];p,q=a['node'],b['node']
        g,sc=joined(ds[di],ds[dj],p,q);off=ds[di]['order'];black=u+[off+x for x in v]
        white=[i for i in range(len(g))if i not in black];bad=obstruction(extract(g,white))
        row=dict(left=[di,u,p],right=[dj,v,q],white=white)
        product_inverse=white_diagonal(di,tuple(u),p)*white_diagonal(dj,tuple(v),q)
        assert bool(bad)==(product_inverse>=1)
        if bad:row['white_obstruction']=bad
        else:
            key=(di,tuple(u),dj,tuple(v));good.setdefault(key,[]).append((p,q))
            tails=[]
            for index,marks,node in ((di,u,p),(dj,v,q)):
                gg=graphs[index];ww=[i for i in range(len(gg))if i not in marks]
                paths=[[node,x,y]for x,y in product(ww,repeat=2)if len({node,x,y})==3 and gg[node][x]and gg[x][y]and not gg[node][y]]
                assert paths;tails.append(min(paths))
            sig=inertia(g);assert sig[1]==2
            row.update(simple_join_inertia=sig,white_tails=tails)
            typ=tuple(sorted(((di,p),(dj,q))))
            if typ not in types:
                (aa,bb),(cc,dd)=typ;h,gauge=joined(ds[aa],ds[cc],bb,dd)
                types[typ]=dict(ports=list(map(list,typ)),inertia=inertia(h),gram=encode(h),gauge=gauge)
        joins.append(row)
    multiple=[]
    for (di,u,dj,v),edges in sorted(good.items()):
        native,g,sc,mp=combined(ds,di,dj,());black=list(u)+[mp[i]for i in v]
        for size in range(2,len(edges)+1):
            for chosen in combinations(edges,size):
                h=[r.copy()for r in native]
                for i,j in chosen:h[i][mp[j]]=h[mp[j]][i]=-1
                a,gauge=native_to_rational(h);ob=edge_bad(a,black);assert ob
                multiple.append(dict(left=[di,list(u)],right=[dj,list(v)],edges=list(map(list,chosen)),gauge=gauge,**ob))
    assert len(joins)==1478 and sum('white_obstruction'in r for r in joins)==363
    assert Counter(map(len,good.values()))=={1:905,2:96,3:6}
    assert len(types)==35 and len(multiple)==120
    print('HIGH_DISJOINT_PASS',len(joins),len(multiple),flush=True)
    high=[r for r in contexts if max(r['distances'])>=4]
    assert all(sum(x>=4 for x in r['distances'])<=4 for r in high)
    records=[];candidate_count=0
    for left,right in product(high,repeat=2):
        ga,gb=graphs[left['diagram']],graphs[right['diagram']]
        choices=[(i,j)for i,d in enumerate(left['distances'])for j,e in enumerate(right['distances'])if d+e>=11]
        size_bound=min(sum(x>=4 for x in left['distances']),sum(x>=4 for x in right['distances']))
        candidate_count+=sum(comb(len(choices),k)for k in range(1,min(size_bound,len(choices))+1))
        ls=sorted({i for i,j in choices});rs=sorted({j for i,j in choices})
        for size in range(1,min(size_bound,len(ls),len(rs))+1):
            for aa,bb in product(combinations(ls,size),combinations(rs,size)):
                for pp in permutations(bb):
                    match=tuple(zip(aa,pp))
                    if any(pair not in choices for pair in match):continue
                    if any(ga[i][k]!=gb[j][l]for i,j in match for k,l in match):continue
                    row=dict(left=[left['diagram'],list(left['black'])],right=[right['diagram'],list(right['black'])],match=[list(x)for x in match])
                    rejects=[]
                    for i,j in match:
                        for side,x,y,rr,other,shared in [('left',i,j,left,gb,set(pp)),('right',j,i,right,ga,set(aa))]:
                            ext=[k for k in range(len(other))if k!=y and other[y][k]and k not in shared]
                            if ext and max(diagonals(rr['diagram'],rr['black'],x))>=2:
                                assert all(k not in(right['black']if side=='left'else left['black'])for k in ext)
                                rejects.append([side,x,ext])
                    if rejects:row.update(kind='inverse_rejection',inverse_rejections=rejects)
                    else:
                        native,g,sc,mp=combined(ds,left['diagram'],right['diagram'],match)
                        u=list(left['black']);v=[mp[i]for i in right['black']];bad=edge_bad(g,u+v)
                        if bad:row.update(kind='edge_rejection',gauge=sc,edge_obstruction=bad)
                        else:row['kind']=overlap_type(left['diagram'],right['diagram'],match)
                    records.append(row)
    counts=Counter(r['kind']for r in records)
    assert candidate_count==22881 and len(records)==4862,(candidate_count,len(records),counts)
    assert counts==dict(inverse_rejection=4296,edge_rejection=254,L10_L25_vertex=6,L18_L25_vertex=6,
        L22_L25_short_vertex=18,L22_L25_long_vertex=54,L22_L25_edge=18,
        L25_vertex=114,L25_edge=60,L25_three_path=27,L25_four_path=9)
    print('OVERLAP_PASS',dict(counts),flush=True)
    return dict(marked_contexts=384,ports=ports,allowed_ports=161,disjoint_joins=joins,
        surviving_disjoint_types=list(types.values()),multiple_edge_obstructions=multiple,
        candidate_matchings=candidate_count,consistent_matchings=len(records),overlap_counts=dict(counts),overlap_records=records)


def overlap_type(di,dj,match):
    if di>dj:return overlap_type(dj,di,tuple(sorted((j,i)for i,j in match)))
    if di==dj==25:
        return {((9,9),):'L25_vertex',((8,9),(9,8)):'L25_edge',
            ((7,9),(8,8),(9,7)):'L25_three_path',((6,9),(7,8),(8,7),(9,6)):'L25_four_path'}[match]
    return {(10,25,((5,9),)):'L10_L25_vertex',(18,25,((7,9),)):'L18_L25_vertex',
        (22,25,((4,9),)):'L22_L25_short_vertex',(22,25,((8,9),)):'L22_L25_long_vertex',
        (22,25,((7,9),(8,8))):'L22_L25_edge'}[di,dj,match]


def pieces(g):
    rest=set(range(len(g)));out=[]
    while rest:
        seen={min(rest)};todo=list(seen);rest-=seen
        for i in todo:
            for j in sorted(rest):
                if g[i][j]:seen.add(j);todo.append(j);rest.remove(j)
        out.append(sorted(seen))
    return out


def inspect_corner(g,sc,black,p,q):
    omitted=sorted(set(black)-{p,q});ids=[i for i in range(len(g))if i not in omitted]
    h=extract(g,ids);psd(h);allowed=omitted.copy();blocks=[];marked=None
    for part in pieces(h):
        actual=[ids[i]for i in part];a=extract(h,part);ks=nullspace(a)
        if ks:blocks.append(dict(nodes=actual,kind='affine',positive_kernel=list(map(str,positive_kernel(a)))))
        else:
            psd(a,True);iv=inverse(a);assert all(x>=0 for row in iv for x in row)
            vals=[sc[actual[i]]*iv[i][i]for i in range(len(a))]
            allowed.extend(actual[i]for i,x in enumerate(vals)if x<=2)
            blocks.append(dict(nodes=actual,kind='finite',normalized_inverse_diagonals=list(map(str,vals))))
        if p in actual and q in actual:marked=actual
    return dict(preferred_corner=[p,q],omitted_black=omitted,blocks=blocks,possible_contacts=sorted(allowed),marked_component=marked)


def isolation_checks(diagrams,data):
    ds={d['index']:d for d in diagrams};out=[];cores={};extra_count=short_count=0
    for row in data['overlap_records']:
        typ=row['kind']
        if typ.endswith('rejection'):continue
        di,u=row['left'];dj,v0=row['right'];match=tuple(map(tuple,row['match']))
        native,g,sc,mp=combined(ds,di,dj,match);v=[mp[i]for i in v0];n=len(g)
        du,dv=graph_dist(native,u),graph_dist(native,v);assert min(du[i]for i in v)>10
        key=str((di,dj,match));ks=nullspace(g);sig=inertia(g);assert sig[1]==1
        if key not in cores:
            cores[key]=dict(gram=encode(g),gauge=sc,signature=sig,radical=list(map(str,ks[0]))if ks else None)
        iv=None if ks else inverse(g)
        if iv:assert set(sc)=={1}
        choices=[];allcorners=[]
        for p,q in product(u,v):
            c=inspect_corner(g,sc,u+v,p,q);allcorners.append(c)
            if c['marked_component']is None:continue
            allowed=c['possible_contacts'];component=len(c['marked_component'])
            if ks:
                assert len(ks)==1;z=ks[0]
                if any(not z[i]for i in allowed):continue
                pos=[i for i in allowed if z[i]>0];neg=[i for i in allowed if z[i]<0]
                bounds=[max([du[i]for i in pos],default=0)+max([dv[i]for i in neg],default=0)+2,
                        max([dv[i]for i in pos],default=0)+max([du[i]for i in neg],default=0)+2]
                if min(bounds)>10:continue
                c.update(method='radical',shortcut_bound=min(bounds),component_order_bound=component)
            else:
                left=[i for i in allowed if du[i]<=4];right=[i for i in allowed if dv[i]<=4]
                if set(left+right)!=set(allowed)or set(left)&set(right):continue
                if not all(iv[i][j]>0 for block in(left,right)for i in block for j in block):continue
                candidates=[i for i in allowed if iv[i][i]<=2]
                if any(iv[i][i]<=1 for i in candidates):continue
                if any(iv[i][i]+2*iv[i][j]+iv[j][j]<=2 for block in(left,right)for i,j in combinations(block,2)):continue
                tests=[];valid=True
                for i,j in product(candidates,repeat=2):
                    a=iv[i][j];rr=(2-iv[i][i])*(2-iv[j][j]);bound=min(du[i]+dv[j],du[j]+dv[i])+3
                    if a*a<=rr or (a<0 and bound>10):valid=False;break
                    tests.append(dict(nodes=[i,j],projection_product=str(a),residual_product=str(rr),
                        reason='positive_product_forbidden'if a>0 else'forced_cross_contact',shortcut_bound=bound if a<0 else None))
                if not valid:continue
                extensions=[]
                for i in candidates:
                    rho=2-iv[i][i];assert rho in(0,F(1,2),F(2,3));tail={F(0):0,F(1,2):1,F(2,3):2}[rho]
                    extensions.append(dict(contact_node=i,contact='-1',projection_square=str(iv[i][i]),residual_square=str(rho),maximum_attached_orthogonal_rank=tail))
                bound=component+max([1+e['maximum_attached_orthogonal_rank']for e in extensions],default=0)
                c.update(method='lorentz_residual',end_sets=[left,right],single_contact_candidates=extensions,
                    candidate_pair_tests=tests,component_order_bound=bound)
            if c['component_order_bound']<=17:choices.append(c)
        assert choices,(di,u,dj,v0,match)
        chosen=min(choices,key=lambda c:(len(c.get('single_contact_candidates',[])),len(c['marked_component']),c['preferred_corner']))
        extras=[];short=0
        left=[i for i in range(ds[di]['order'])if i not in{x for x,y in match}]
        right=[mp[i]for i in range(ds[dj]['order'])if i not in{y for x,y in match}]
        for x,y in product(left,right):
            if native[x][y]:continue
            h=[r.copy()for r in native];h[x][y]=h[y][x]=-1
            if min(graph_dist(h,u)[j]for j in v)<=10:short+=1;continue
            a,gauge=native_to_rational(h);o=edge_bad(a,u+v);assert o
            extras.append(dict(edge=[x,y],gauge=gauge,**o))
        extra_count+=len(extras);short_count+=short
        out.append(dict(type=typ,left=[di,u],right=[dj,v0],match=row['match'],right_copy_map=mp,
            all_corners=allcorners,extra_cross_shortcuts=short,extra_cross_pd_rejections=extras,**chosen))
    assert len(out)==312 and len(cores)==14 and short_count==21393 and extra_count==168
    assert max(r['component_order_bound']for r in out)==17
    print('ISOLATION_PASS',len(out),short_count,extra_count,flush=True)
    return dict(cores=cores,preferred_corners=out,extra_cross_shortcuts=short_count,extra_cross_pd_rejections=extra_count)



def polynomial(ranks):
    p=[1]
    for r in ranks:
        q=[0]*(len(p)+r)
        for i,t in enumerate(p):
            for j in range(r+1):q[i+j]+=t
        p=q
    return p


def profile(ranks):
    d=sum(ranks);h=polynomial(ranks);v=prod(r+1 for r in ranks)
    assert h==h[::-1]and sum(h)==v
    G=sum(t*min(j,d-j)for j,t in enumerate(h));curv=F(5*G-d*v,2*v)
    moment=0;s=len(ranks)
    for k in range(s+1):
        for subset in combinations(ranks,k):
            shift=sum(x+1 for x in subset);m=(d-1)//2-shift
            if m>=0:moment+=(-1)**k*((d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1))
    assert curv==F(3*d,4)-F(5*moment,2*v)
    cost=sum((min(F(10),F(r,2))for r in ranks),F())
    return curv,cost,curv-cost


def budgets():
    rows=[]
    for d in range(116,126):
        shapes=set()
        for k in range(d//20+1):
            for low in range(20):
                high=d-20*k-low
                if high==0 or high>=21:shapes.add(tuple(sorted([20]*k+([low]if low else[])+([high]if high else[]))))
        for ranks in sorted(shapes):
            a,b,c=profile(ranks);assert c>0
            rows.append(dict(ranks=list(ranks),normalized_curvature=str(a),normalized_cost=str(b),margin=str(c)))
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8);assert beta>10*n
    assert len(rows)==1015 and min(F(r['margin'])for r in rows)==F(145234,4084101)
    boundary=profile((20,)*5+(15,))[2];assert boundary==-F(673055,16336404)
    polynomials=[]
    for start,expected in [(87,[12,176,4]),(88,[96,180,4])]:
        vals=[]
        for t in range(4):
            d=start+2*t;beta=F(d*(d-6),8)if d%2==0 else F((d-1)*(d-5),8)
            vals.append(8*(beta-10*(d+1)))
        c2=(vals[2]-2*vals[1]+vals[0])/2;cs=[vals[0],vals[1]-vals[0]-c2,c2]
        assert cs==expected and all(vals[t]==sum(cs[i]*t**i for i in range(3))for t in range(4))
        polynomials.append(list(map(str,cs)))
    return dict(profiles=rows,minimum_margin='145234/4084101',boundary_margin=str(boundary),compact_facet_slacks=polynomials)



class K:
    """Exact a+b*sqrt(5), with order decided by rational square comparison."""
    def __init__(self,a=0,b=0):
        if isinstance(a,K):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,o):o=K(o);return K(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,o):return self+-K(o)
    def __rsub__(self,o):return K(o)+-self
    def __mul__(self,o):o=K(o);return K(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=K(o);norm=o.a*o.a-5*o.b*o.b;assert norm
        return self*K(o.a/norm,-o.b/norm)
    def __rtruediv__(self,o):return K(o)/self
    def __bool__(self):return bool(self.a or self.b)
    def sign(self):
        a,b=self.a,self.b
        if not b:return (a>0)-(a<0)
        if not a:return (b>0)-(b<0)
        if a>0 and b>0:return 1
        if a<0 and b<0:return -1
        c=a*a-5*b*b;assert c
        return ((c>0)-(c<0))*(1 if a>0 else -1)
    def __eq__(self,o):o=K(o);return self.a==o.a and self.b==o.b
    def __lt__(self,o):return (self-o).sign()<0
    def __le__(self,o):return (self-o).sign()<=0
    def __gt__(self,o):return (self-o).sign()>0
    def __ge__(self,o):return (self-o).sign()>=0
    def save(self):return [str(self.a),str(self.b)]


def low_five(diagrams):
    # Matching-polynomial recurrence uses squared bonds only; no irrational
    # matrix entry, symbolic package, or producer determinant is imported.
    shapes=[[(0,1),(1,2),(2,3),(3,4)],[(0,1),(0,2),(2,3),(3,4)],[(0,2),(1,2),(2,3),(3,4)]]
    beta={3:K(1),4:K(2),5:K(F(3,2),F(1,2))};rows=[];survivors=[];counts=Counter()
    assert all(0<x<4 for x in beta.values())
    for shape,edges in enumerate(shapes):
        for labels in product((3,4,5),repeat=4):
            weights={tuple(sorted(e)):beta[m]for e,m in zip(edges,labels)}
            @lru_cache(None)
            def det(nodes):
                if not nodes:return K(1)
                i=nodes[0];rest=nodes[1:];value=2*det(rest)
                for j in rest:
                    w=weights.get(tuple(sorted((i,j))))
                    if w:value-=w*det(tuple(k for k in rest if k!=j))
                return value
            row=dict(shape=shape,labels=list(labels));bad=None
            for size in(3,4):
                for ids in combinations(range(5),size):
                    value=det(ids)
                    if value<0:bad=dict(nodes=list(ids),determinant_coefficients=value.save());break
                if bad:break
            if bad:kind='proper_indefinite';row['witness']=bad
            elif det(tuple(range(5)))>=0:kind='not_indefinite';row['full_determinant_coefficients']=det(tuple(range(5))).save()
            else:
                vals=[]
                for b in(0,1):
                    ids=tuple(i for i in range(5)if i!=b);den=det(ids)
                    if den<=0:break
                    vals.append(det(tuple(i for i in ids if i!=4))/den)
                if len(vals)<2:kind='black_not_PD'
                else:
                    row['inverse_diagonals']=[x.save()for x in vals]
                    if max(vals)>=2:kind='port_rejected'
                    else:kind='survive';survivors.append(row)
            counts[kind]+=1;row['kind']=kind;rows.append(row)
    assert counts==dict(proper_indefinite=209,not_indefinite=14,black_not_PD=4,port_rejected=12,survive=4)
    big=next(d for d in diagrams if d['index']==25);joins=[]
    for r in survivors:
        edges=[(i,j,beta[m])for(i,j),m in zip(shapes[r['shape']],r['labels'])]
        edges += [(i+5,j+5,K(1))for i,j,m in big['edges']]+[(4,14,K(1))]
        assert all(m==3 for i,j,m in big['edges'])
        diagonal={i:K(2)for i in range(15)};adj={i:{}for i in range(15)}
        for i,j,w in edges:adj[i][j]=adj[j][i]=w
        pos=neg=zero=0;steps=[]
        while diagonal:
            i=next(i for i in diagonal if len(adj[i])<=1);q=diagonal[i]
            if not adj[i]:
                pos+=int(q>0);neg+=int(q<0);zero+=int(q==0);steps.append(dict(nodes=[i],pivot=q.save()))
                del adj[i];del diagonal[i]
            else:
                j=next(iter(adj[i]));w=adj[i][j]
                if q:
                    pos+=int(q>0);neg+=int(q<0);diagonal[j]-=w/q
                    steps.append(dict(nodes=[i],pivot=q.save(),neighbor=j,squared_bond=w.save()))
                    del adj[j][i];del adj[i];del diagonal[i]
                else:
                    pos+=1;neg+=1;steps.append(dict(zero_leaf_pair=[i,j],squared_bond=w.save()))
                    for k in list(adj[j]):
                        if k!=i:del adj[k][j]
                    del adj[i];del adj[j];del diagonal[i];del diagonal[j]
        assert [pos,neg,zero]==[13,2,0]
        joins.append(dict(shape=r['shape'],labels=r['labels'],signature=[pos,neg,zero],leaf_schur_steps=steps))
    print('LOW5_PASS',dict(counts),'4 inertia types /24 ordered joins',flush=True)
    return dict(counts=dict(counts),records=rows,survivors=survivors,joined_inertias=joins,ordered_marked_joins=24)


def parse_k(text):
    import ast
    def f(node):
        if isinstance(node,ast.Constant):return K(node.value)
        if isinstance(node,ast.UnaryOp)and isinstance(node.op,ast.USub):return -f(node.operand)
        if isinstance(node,ast.BinOp):
            a,b=f(node.left),f(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
        if isinstance(node,ast.Call)and isinstance(node.func,ast.Name)and node.func.id=='sqrt':
            assert len(node.args)==1 and isinstance(node.args[0],ast.Constant)and node.args[0].value==5
            return K(0,1)
        raise AssertionError(ast.dump(node))
    return f(ast.parse(text,mode='eval').body)


def compare_saved(diagrams,data,iso,bs,low,p):
    assert p['dimension_bound']==116 and p['inherited_dimension_bound']==126
    assert p['unrestricted_finite_volume']is True and p['arithmetic_hypothesis']is False and p['cusp_restriction']is None
    assert p['census']['ports']==data['ports']and p['census']['allowed_ports']==161
    assert p['census']['overlap_generated']==data['candidate_matchings']and p['census']['consistent_overlaps']==data['consistent_matchings']
    ds={d['index']:d for d in diagrams};nchecks=0
    def witness(g,sc,o):
        nonlocal nchecks
        assert o['gauge_squared']==sc;v=list(map(F,o['vector']));ids=o['nodes']
        assert len(v)==len(ids)and all(x>=0 for x in v)and any(v)
        assert square(extract(g,ids),v)==F(o['square'])<=0;nchecks+=1
    for r,s in zip(data['disjoint_joins'],p['disjoint_obstructions']['rows']):
        assert r['left']==s['left']and r['right']==s['right']
        di,u,i=r['left'];dj,v,j=r['right'];g,sc=joined(ds[di],ds[dj],i,j)
        if 'white_obstruction'in r:witness(g,sc,s['white_rejection'])
        else:
            assert 'white_rejection'not in s and r['white_tails']==s['white_three_node_tails']
    assert len(p['disjoint_obstructions']['rows'])==1478
    for s in p['disjoint_obstructions']['underlying_inertia_types']:
        (di,i),(dj,j)=s['ports'];g,sc=joined(ds[di],ds[dj],i,j)
        assert encode(g)==s['gram']and sc==s['gauge_squared']and inertia(g)==s['signature']and s['signature'][1]==2
    saved_multiple={(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['edges']))):r for r in p['disjoint_obstructions']['multiple_edge_obstructions']}
    assert len(saved_multiple)==len(data['multiple_edge_obstructions'])==120
    for r in data['multiple_edge_obstructions']:
        di,u=r['left'];dj,v=r['right'];s=saved_multiple[di,tuple(u),dj,tuple(v),tuple(map(tuple,r['edges']))]
        native,g,sc,mp=combined(ds,di,dj,())
        for i,j in r['edges']:native[i][mp[j]]=native[mp[j]][i]=-1
        g,sc=native_to_rational(native);witness(g,sc,s['edge_pd_rejection'])
    saved={overlap_key(r):r for r in p['census']['overlap_records']};assert len(saved)==4862
    for r in data['overlap_records']:
        s=saved[overlap_key(r)]
        if r['kind']=='inverse_rejection':
            o=s['inverse_rejection'];assert any(o['side']==side and o['node']==node and o['external_neighbor']in ext for side,node,ext in r['inverse_rejections'])
            assert max(map(F,o['inverse_diagonals']))>=2
        elif r['kind']=='edge_rejection':
            di,u=r['left'];dj,v=r['right'];native,g,sc,mp=combined(ds,di,dj,tuple(map(tuple,r['match'])))
            witness(g,sc,s['edge_pd_rejection'])
        else:assert s['surviving_type']==r['kind']
    saved_iso={overlap_key(r):r for r in p['isolation']['cases']};assert len(saved_iso)==312
    for r in iso['preferred_corners']:
        s=saved_iso[overlap_key(r)];di,u=r['left'];dj,v=r['right'];native,g,sc,mp=combined(ds,di,dj,tuple(map(tuple,r['match'])))
        for key in('type','preferred_corner','omitted_black','blocks','possible_contacts','marked_component','method','component_order_bound'):
            assert r[key]==s[key],(overlap_key(r),key)
        if r['method']=='radical':assert r['shortcut_bound']==s['shortcut_bound']
        else:
            for key in('end_sets','single_contact_candidates','candidate_pair_tests'):assert r[key]==s[key],key
        cross=s['extra_cross_bond_checks'];assert cross['shortcut_rejections']==r['extra_cross_shortcuts']
        assert len(cross['edge_pd_z_rejections'])==len(r['extra_cross_pd_rejections'])
        for o in cross['edge_pd_z_rejections']:
            x,y=o['edge'];h=[row.copy()for row in native];h[x][y]=h[y][x]=-1
            gg,gauge=native_to_rational(h);witness(gg,gauge,o)
    lp=p['low_order_disjoint'];assert lp['counts']==low['counts']and lp['ordered_marked_joins']==24
    high=[r for r in data['ports']if r['allowed']and r['distance']>=7]
    assert [(r['diagram'],r['black'],r['node'])for r in high]==[(25,[1,2],9),(25,[1,3],9),(25,[2,3],9)]
    assert all(r['high_contexts']==[q['black']for q in high]and r['ordered_marked_multiplicity']==6 for r in lp['exact_join_types'])
    assert len(lp['records'])==243
    for r,s in zip(low['records'],lp['records']):
        assert (r['shape'],r['labels'])==(s['shape'],s['labels'])
        kind=('proper_indefinite'if'proper_indefinite'in s else'not_indefinite'if'nonnegative_full_determinant'in s else'black_not_PD'if'black_not_PD'in s else'port_rejected'if'port_rejected'in s else'survive')
        assert r['kind']==kind
        if kind=='proper_indefinite':
            assert r['witness']['nodes']==s['proper_indefinite']['nodes']
            assert parse_k(s['proper_indefinite']['determinant']).save()==r['witness']['determinant_coefficients']
        if kind=='not_indefinite':assert parse_k(s['nonnegative_full_determinant']).save()==r['full_determinant_coefficients']
        if 'inverse_diagonals'in r:assert [parse_k(x).save()for x in s['inverse_diagonals']]==r['inverse_diagonals']
    for r,s in zip(low['joined_inertias'],lp['exact_join_types']):assert (r['shape'],r['labels'],r['signature'])==(s['low_shape'],s['low_labels'],s['signature'])
    look={tuple(sorted(r['ranks'])):r for r in p['reduced_profiles']};assert len(look)==1015
    for r in bs['profiles']:
        s=look[tuple(r['ranks'])]
        for key in('normalized_curvature','normalized_cost','margin'):assert F(r[key])==F(s[key])
    assert F(p['boundary_at_n116']['margin'])==F(bs['boundary_margin'])
    assert p['compact_facet_ambient_bound']==87
    assert [list(map(int,row))for row in bs['compact_facet_slacks']]==p['compact_facet_parity_slack_times8']
    assert nchecks==905
    return dict(independently_verified_saved_nonnegative_witnesses=nchecks,
        producer_high_and_low_classification_and_all_extensions_match=True)


def run(unbound=False):
    diagrams=json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']
    data=classify(diagrams);iso=isolation_checks(diagrams,data);bs=budgets();low=low_five(diagrams)
    p=json.loads((ROOT/f'results/{BASE}.json').read_text());comparison=compare_saved(diagrams,data,iso,bs,low,p)
    for path,h in p['sources'].items():assert sha(ROOT/path)==h
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'results/{BASE}.md',
           f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
           'scripts/fresh_distance10_extension_20261017.py','results/fresh_distance10_extension_20261017.json','results/fresh_distance10_extension_20261017.md',
           'results/fresh_selective_distance11_dimension126_20261017.md','results/fresh_cusp_link_shape_curvature_20261017.md']
    if unbound:paths=[p for p in paths if (ROOT/p).exists()]
    return dict(status='PASS_INDEPENDENT_FRACTION_SELECTIVE_DISTANCE10_DIMENSION116_REVIEW'if not unbound else'PASS_UNBOUND_DISTANCE10_REVIEW',
        sources={path:sha(ROOT/path)for path in paths},producer_source_bindings=len(p['sources']),
        classification=data,low_order_disjoint=low,isolation=iso,budgets=bs,saved_comparison=comparison,
        dimension_bound=116,compact_facet_ambient_bound=87,
        scope='Independent Fraction and Q(sqrt5) exact finite reconstruction plus mathematical/source review. The primary finite classifications, face/cusp incidence and previous all-rank curvature/convex-order proofs remain mathematical dependencies. No formalization, external refereeing, priority or new-construction claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--unbound',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run(a.unbound)));path=Path('/private/tmp/distance10_review_unbound.json')if a.unbound else ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'243low/527ports/4862overlaps/312corners/1015profiles',flush=True)
