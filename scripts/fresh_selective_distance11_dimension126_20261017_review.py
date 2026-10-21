#!/usr/bin/env python3
"""Independent Fraction-only review of selective-distance11 dimension126.

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
TAG='fresh_selective_distance11_dimension126_20261017_review'
BASE='fresh_selective_distance11_dimension126_20261017'
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
    bits={0:0};todo=deque([0])
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
    ports=[]
    for r in contexts:
        for node,distance in enumerate(r['distances']):
            if distance<4:continue
            vals=diagonals(r['diagram'],r['black'],node)
            ports.append(dict(diagram=r['diagram'],black=list(r['black']),node=node,distance=distance,
                              inverse_diagonals=list(map(str,vals)),allowed=max(vals)<2))
    allowed=[p for p in ports if p['allowed']];assert len(allowed)==55
    joins=[];unique={};types={}
    for a,b in product(allowed,repeat=2):
        if a['distance']+b['distance']<11:continue
        di,dj=a['diagram'],b['diagram'];u,v=a['black'],b['black'];p,q=a['node'],b['node']
        g,sc=joined(ds[di],ds[dj],p,q);off=ds[di]['order'];black=u+[off+x for x in v]
        white=[i for i in range(len(g))if i not in black];bad=obstruction(extract(g,white))
        row=dict(left=[di,u,p],right=[dj,v,q],white=white)
        if bad:row['white_obstruction']=bad
        else:
            key=(di,tuple(u),dj,tuple(v));unique.setdefault(key,[]).append((p,q))
            # Exact full-white Schur bound on the only possible cross bond.
            values=[]
            for d,bb,nn in [(ds[di],u,p),(ds[dj],v,q)]:
                gg,scale=rational_gram(d);ww=[i for i in range(len(gg))if i not in bb];j=inverse(extract(gg,ww))
                values.append(scale[nn]*j[ww.index(nn)][ww.index(nn)])
            upper=1/(values[0]*values[1]);assert 1<upper<=2
            sig=inertia(g);assert sig[1]>=2
            row.update(cross_bond_square_upper_bound=str(upper),simple_join_inertia=sig)
            typ=tuple(sorted(((di,p),(dj,q))))
            types.setdefault(typ,dict(ports=typ,inertia=sig,gram=encode(g),gauge=sc,count=0))['count']+=1
        joins.append(row)
    assert len(joins)==425 and sum('white_obstruction'in r for r in joins)==119
    assert len(unique)==306 and all(len(v)==1 for v in unique.values())
    high=[r for r in contexts if max(r['distances'])>=5]
    assert all(sum(x>=5 for x in r['distances'])<=3 for r in high)
    records=[];candidate_count=0
    for left,right in product(high,repeat=2):
        ga,gb=graphs[left['diagram']],graphs[right['diagram']]
        choices=[(i,j)for i,d in enumerate(left['distances'])for j,e in enumerate(right['distances'])if d+e>=12]
        size_bound=min(sum(x>=5 for x in left["distances"]),sum(x>=5 for x in right["distances"]))
        candidate_count+=sum(comb(len(choices),k)for k in range(1,min(size_bound,len(choices))+1))
        ls=sorted({i for i,j in choices});rs=sorted({j for i,j in choices})
        for size in range(1,min(len(ls),len(rs))+1):
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
                        g1,s1=rational_gram(ds[left['diagram']]);g2,s2=rational_gram(ds[right['diagram']]);assert set(s1+s2)=={1}
                        g,mp=unite(g1,g2,match);u=list(left['black']);v=[mp[i]for i in right['black']]
                        white=[i for i in range(len(g))if i not in u+v];bad=None
                        for mark in u+v:
                            ids=white+[mark];o=obstruction(extract(g,ids))
                            if o:bad=dict(mark=mark,nodes=ids,**o);break
                        if bad:row.update(kind='edge_rejection',edge_obstruction=bad)
                        else:
                            di,dj=left['diagram'],right['diagram']
                            if(di,dj)==(25,25):
                                kind={((9,9),):'shared_vertex',((8,9),(9,8)):'shared_edge',((7,9),(8,8),(9,7)):'shared_path'}[match]
                            else:
                                assert(di,dj,match)in[(22,25,((8,9),)),(25,22,((9,8),))];kind='mixed'
                            row['kind']=kind
                    records.append(row)
    assert candidate_count==2137 and len(records)==1010,(candidate_count,len(records),Counter(r["kind"]for r in records))
    counts=Counter(r['kind']for r in records)
    assert counts=={'inverse_rejection':856,'edge_rejection':40,'shared_vertex':60,'shared_edge':27,'shared_path':9,'mixed':18}
    return dict(marked_contexts=384,ports=ports,allowed_ports=55,disjoint_joins=joins,
                surviving_disjoint_types=list(types.values()),candidate_matchings=candidate_count,
                consistent_matchings=len(records),overlap_counts=dict(counts),overlap_records=records)



def pieces(g):
    rest=set(range(len(g)));out=[]
    while rest:
        seen={min(rest)};todo=list(seen);rest-=seen
        for i in todo:
            for j in sorted(rest):
                if g[i][j]:seen.add(j);todo.append(j);rest.remove(j)
        out.append(sorted(seen))
    return out


def preferred(black):
    return 4 if 4 in black else 0 if 0 in black else 2 if tuple(black)in((1,2),(2,3))else 1


def edge_bad(g,black):
    white=[i for i in range(len(g))if i not in black]
    for p in black:
        ids=white+[p];o=obstruction(extract(g,ids))
        if o:return dict(mark=p,nodes=ids,**o)
    return None


def isolation_checks(diagrams,data):
    ds={d['index']:d for d in diagrams};out=[];cores={};extra_count=short_count=0
    for r in data['overlap_records']:
        typ=r['kind']
        if typ.endswith('rejection'):continue
        di,u=r['left'];dj,v0=r['right'];match=list(map(tuple,r['match']))
        a,sa=rational_gram(ds[di]);b,sb=rational_gram(ds[dj]);assert set(sa+sb)=={1}
        g,mp=unite(a,b,match);v=[mp[i]for i in v0];n=len(g)
        du,dv=graph_dist(g,u),graph_dist(g,v);assert min(du[i]for i in v)>11
        white=[i for i in range(n)if i not in u+v]
        for p,q in product(u,v):psd(extract(g,white+[p,q]))
        p,q=preferred(u),mp[preferred(v0)];omitted=sorted(set(u+v)-{p,q});ids=[i for i in range(n)if i not in omitted]
        h=extract(g,ids);allowed=omitted.copy();blocks=[];marked=None
        for part in pieces(h):
            actual=[ids[i]for i in part];hh=extract(h,part);ks=nullspace(hh)
            if ks:
                k=positive_kernel(hh);blocks.append(dict(nodes=actual,kind='affine',positive_kernel=list(map(str,k))))
            else:
                psd(hh,True);iv=inverse(hh);assert all(x>=0 for rr in iv for x in rr)
                vals=[iv[i][i]for i in range(len(hh))]
                allowed.extend(actual[i]for i,x in enumerate(vals)if x<=2)
                blocks.append(dict(nodes=actual,kind='finite',inverse_diagonals=list(map(str,vals))))
            if p in actual:assert q in actual;marked=actual
        assert marked is not None and len(marked)<=17
        allowed.sort();distance=graph_dist(h,[ids.index(p)],False)[ids.index(q)];assert distance>11
        sig=inertia(g);radical=nullspace(g);extensions=[]
        if typ in('shared_vertex','mixed'):
            assert len(radical)==1 and sig==[n-2,1,1];z=radical[0]
            if z[0]>0:z=[-x for x in z]
            shared={i for i,j in match};left=set(range(len(a)))-shared;right=set(mp)-shared
            assert all(z[i]<0 for i in left)and all(z[i]>0 for i in right)and all(z[i]==0 for i in shared)
            assert all(z[i]for i in allowed)
            bound_left=max(du[i]for i in allowed if z[i]<0);bound_right=max(dv[i]for i in allowed if z[i]>0)
            bound=bound_left+bound_right+2;assert bound<11
            if typ=='shared_vertex':assert bound_left<=3 and bound_right<=3
            else:assert bound_left<=2 and bound_right<=2 and len(marked)<=16
            core_key=typ if typ=='shared_vertex'else f'mixed_{di}_{dj}'
            cores.setdefault(core_key,dict(gram=encode(g),signature=sig,radical=list(map(str,z))))
        else:
            assert not radical and sig==[n-1,1,0];iv=inverse(g)
            end=list(range(4))if typ=='shared_edge'else[1,2,3];other=[mp[i]for i in end]
            block=extract(iv,end);assert block==extract(iv,other)and all(x>0 for row in block for x in row)
            assert set(allowed)<=set(end+other)
            dl=max(du[i]for i in allowed if i in end);dr=max(dv[i]for i in allowed if i in other);bound=dl+dr+2
            if typ=='shared_edge':assert min(block[i][i]for i in range(len(block)))>=4 and dl<=2 and dr<=2
            else:
                assert block==[[F(7,2),4,2],[4,6,3],[2,3,2]]and dl<=1 and dr<=1
                rr=[]
                for node in (3,mp[3]):
                    root=[-row[node]for row in iv];assert square(g,root)==2
                    assert multiply(g,[[x]for x in root])==[[-F(i==node)]for i in range(n)]
                    rr.append(root);ext=matrix(n+1)
                    for i in range(n):
                        for j in range(n):ext[i][j]=g[i][j]
                    ext[n][node]=ext[node][n]=-1;psd(extract(ext,white+[n]),True)
                    for pp,qq in product(u,v):psd(extract(ext,white+[n,pp,qq]))
                    ec=extract(ext,ids+[n]);ecmarks=None
                    for part in pieces(ec):
                        nodes=[(ids+[n])[i]for i in part]
                        if p in nodes:assert q in nodes;ecmarks=nodes
                    assert ecmarks is not None and len(ecmarks)<=16
                    extensions.append(dict(contact_node=node,coordinates=list(map(str,root)),marked_component=ecmarks))
                assert sum((rr[0][i]*g[i][j]*rr[1][j]for i in range(n)for j in range(n)),F())==-2
            cores.setdefault(typ,dict(gram=encode(g),signature=sig,end_inverse_block=encode(block)))
        assert bound<=8
        extras=[];short=0
        left=[i for i in range(len(a))if i not in {x for x,y in match}]
        right=[mp[j]for j in range(len(b))if j not in {y for x,y in match}]
        for x,y in product(left,right):
            if g[x][y]:continue
            hh=[row.copy()for row in g];hh[x][y]=hh[y][x]=-1
            if min(graph_dist(hh,u)[j]for j in v)<=11:short+=1;continue
            o=edge_bad(hh,u+v);assert o is not None
            extras.append(dict(edge=[x,y],**o))
        extra_count+=len(extras);short_count+=short
        out.append(dict(type=typ,left=[di,u],right=[dj,v0],match=r['match'],left_black=u,right_black=v,
                        preferred_corner=[p,q],blocks=blocks,possible_contacts=allowed,marked_component=marked,
                        marked_distance=distance,shortcut_bound=bound,one_root_extensions=extensions,
                        extra_cross_shortcut_rejections=short,extra_cross_pd_rejections=extras))
    assert len(out)==114 and short_count==8280 and extra_count==45
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
    cost=sum((min(F(11),F(r,2))for r in ranks),F())
    return curv,cost,curv-cost


def budgets():
    rows=[]
    for d in range(126,136):
        shapes=set()
        for k in range(d//22+1):
            for low in range(22):
                high=d-22*k-low
                if high==0 or high>=23:shapes.add(tuple(sorted([22]*k+([low]if low else[])+([high]if high else[]))))
        for ranks in sorted(shapes):
            a,b,c=profile(ranks);assert c>0
            rows.append(dict(ranks=list(ranks),normalized_curvature=str(a),normalized_cost=str(b),margin=str(c)))
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8);assert beta>11*n
    assert len(rows)==1095 and min(F(r['margin'])for r in rows)==F(40603,559682)
    boundary=profile((22,)*5+(15,))[2];assert boundary==-F(238535,12872686)
    polynomials=[]
    for start,expected in [(95,[12,192,4]),(96,[104,196,4])]:
        vals=[]
        for t in range(4):
            d=start+2*t;beta=F(d*(d-6),8)if d%2==0 else F((d-1)*(d-5),8)
            vals.append(8*(beta-11*(d+1)))
        c2=(vals[2]-2*vals[1]+vals[0])/2;cs=[vals[0],vals[1]-vals[0]-c2,c2]
        assert cs==expected and all(vals[t]==sum(cs[i]*t**i for i in range(3))for t in range(4))
        polynomials.append(list(map(str,cs)))
    return dict(profiles=rows,minimum_margin='40603/559682',boundary_margin=str(boundary),compact_facet_slacks=polynomials)



def compare_saved(ds,data,iso,bs,p):
    assert p['dimension_bound']==126 and p['inherited_dimension_bound']==136
    assert p['unrestricted_finite_volume']is True and p['arithmetic_hypothesis']is False and p['cusp_restriction']is None
    assert p['census']['ports']==data['ports']and p['census']['allowed_ports']==55
    assert p['census']['overlap_generated']==data['candidate_matchings']and p['census']['consistent_overlaps']==data['consistent_matchings']
    diagrams={d['index']:d for d in ds};nchecks=0
    def check_witness(g,nodes,o):
        nonlocal nchecks
        v=list(map(F,o['vector']));assert len(v)==len(nodes)and all(x>=0 for x in v)and any(v)
        assert square(extract(g,nodes),v)==F(o['square'])<=0;nchecks+=1
    saved_joins=p['disjoint_obstructions']['rows'];assert len(saved_joins)==425
    for r,s in zip(data['disjoint_joins'],saved_joins):
        assert r['left']==s['left']and r['right']==s['right']
        di,u,i=r['left'];dj,v,j=r['right'];g,scale=joined(diagrams[di],diagrams[dj],i,j)
        if 'white_obstruction'in r:
            o=s['white_rejection'];assert o['nodes']==r['white']and o['gauge_squared']==scale
            check_witness(g,o['nodes'],o)
        else:
            assert 'white_rejection'not in s
            for (di,black,port),tail in zip([r['left'],r['right']],s['white_three_node_tails']):
                nodes=tail['nodes'];gg=bond_graph(diagrams[di]);assert nodes[0]==port and len(set(nodes))==3 and set(nodes).isdisjoint(black)
                assert gg[nodes[0]][nodes[1]]and gg[nodes[1]][nodes[2]]and not gg[nodes[0]][nodes[2]]
            typ=tuple(sorted(((r['left'][0],r['left'][2]),(r['right'][0],r['right'][2]))))
            assert list(map(list,typ))==s['inertia_type']
    assert len(p['disjoint_obstructions']['underlying_inertia_types'])==13
    for r in p['disjoint_obstructions']['underlying_inertia_types']:
        (di,i),(dj,j)=r['ports'];g,scale=joined(diagrams[di],diagrams[dj],i,j)
        assert encode(g)==r['gram']and scale==r['gauge_squared']and inertia(g)==r['signature']and r['signature'][1]==2
    saved={overlap_key(r):r for r in p['census']['overlap_records']};assert len(saved)==1010
    kinds={'shared_three_path':'shared_path','mixed_22_25':'mixed','mixed_25_22':'mixed'}
    for r in data['overlap_records']:
        s=saved[overlap_key(r)]
        if r['kind']=='inverse_rejection':
            o=s['inverse_rejection'];assert any(o['side']==side and o['node']==node and o['external_neighbor']in ext for side,node,ext in r['inverse_rejections'])
            assert max(map(F,o['inverse_diagonals']))>=2
        else:
            di,u=r['left'];dj,v0=r['right'];a,sa=rational_gram(diagrams[di]);b,sb=rational_gram(diagrams[dj]);g,mp=unite(a,b,list(map(tuple,r['match'])))
            if r['kind']=='edge_rejection':
                o=s['edge_pd_rejection'];assert o['gauge_squared']==[1]*len(g);check_witness(g,o['nodes'],o)
            else:assert kinds.get(s['surviving_type'],s['surviving_type'])==r['kind']
    sv={(r['type'],tuple(r['left_black']),tuple(r['right_black'])):r for r in p['isolation']['cases']};assert len(sv)==114
    for r in iso['preferred_corners']:
        di,u=r['left'];dj,v0=r['right'];typ=r['type']
        pt='shared_three_path'if typ=='shared_path'else f'mixed_{di}_{dj}'if typ=='mixed'else typ
        s=sv[pt,tuple(r['left_black']),tuple(r['right_black'])]
        a,sa=rational_gram(diagrams[di]);b,sb=rational_gram(diagrams[dj]);g,mp=unite(a,b,list(map(tuple,r['match'])))
        cross=s['extra_cross_bond_checks'];assert cross['shortcut_rejections']==r['extra_cross_shortcut_rejections']
        assert len(cross['edge_pd_z_rejections'])==len(r['extra_cross_pd_rejections'])
        for o in cross['edge_pd_z_rejections']:
            x,y=o['edge'];h=[row.copy()for row in g];h[x][y]=h[y][x]=-1
            assert o['gauge_squared']==[1]*len(g);check_witness(h,o['nodes'],o)
        if typ in('shared_vertex','shared_edge'):
            for key in ['preferred_corner','blocks','possible_contacts','marked_component','marked_distance']:
                assert r[key]==s[key],(typ,key)
            assert r['shortcut_bound']<=s['shortcut_bound']<=8
        else:
            mapping=s['literal_isomorphism'];target=json.loads((ROOT/f"results/{s['local_lemma']}.json").read_text())
            gg=[list(map(F,row))for row in target['core']['gram']]
            assert g==extract(gg,mapping)
            c=s['local_lemma_context'];corner=c['preferred']if typ=='mixed'else c['preferred']['preferred_corner']
            component=c['marked_component']if typ=='mixed'else c['preferred']['marked_component']
            allowed=c['allowed_contacts']if typ=='mixed'else c['preferred']['possible_contacts']
            assert set(mapping[i]for i in r['preferred_corner'])==set(corner)
            assert set(mapping[i]for i in r['marked_component'])==set(component)
            assert set(mapping[i]for i in r['possible_contacts'])==set(allowed)
            assert s['component_order_bound']==16
    look={tuple(sorted(r['ranks'])):r for r in p['reduced_profiles']};assert len(look)==1095
    for r in bs['profiles']:
        s=look[tuple(r['ranks'])]
        for key in ['normalized_curvature','normalized_cost','margin']:assert F(r[key])==F(s[key])
    assert F(p['boundary_at_n126']['margin'])==F(bs['boundary_margin'])
    assert p['compact_facet_ambient_bound']==95
    assert [list(map(int,row))for row in bs['compact_facet_slacks']]==p['compact_facet_parity_slack_times8']
    assert nchecks==204
    return dict(independently_verified_saved_nonnegative_witnesses=nchecks,producer_classification_and_local_lemma_maps_match=True)


def run():
    ds=json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']
    p=json.loads((ROOT/f'results/{BASE}.json').read_text());data=classify(ds)
    iso=isolation_checks(ds,data);bs=budgets();comparison=compare_saved(ds,data,iso,bs,p)
    for path,h in p['sources'].items():assert sha(ROOT/path)==h
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'results/{BASE}.md',
           f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
           'results/fresh_distance11_extension_review_20261017.md',
           'results/fresh_distance11_extension_review_20261017_mixed.md',
           'results/fresh_selective_distance12_dimension136_20261017.md',
           'results/fresh_cusp_link_shape_curvature_20261017.md']
    return dict(status='PASS_INDEPENDENT_FRACTION_SELECTIVE_DISTANCE11_DIMENSION126_REVIEW',
                sources={path:sha(ROOT/path)for path in paths},producer_source_bindings=len(p['sources']),
                classification=data,isolation=iso,budgets=bs,saved_comparison=comparison,
                dimension_bound=126,compact_facet_ambient_bound=95,
                scope='Independent Fraction-only finite reconstruction and mathematical/source review. Published classifications, actual corner/cusp incidence and prior all-rank curvature/convex-order theorems remain bound mathematical inputs. No formalization, external referee review, novelty or construction claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'384contexts/55ports/425joins/1010overlaps/114corners/1095profiles',flush=True)
