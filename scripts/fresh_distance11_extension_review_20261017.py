#!/usr/bin/env python3
"""Exact supporting-white extensions of the17-root distance11 core.

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
TAG='fresh_distance11_extension_review_20261017'
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


def core(d):
    g=matrix(17);maps=[list(range(10)),list(range(10,17))+[9,8,7]]
    for mp in maps:
        for i,j,m in d['edges']:
            assert m==3;g[mp[i]][mp[j]]=g[mp[j]][mp[i]]=-F(1)
    return g,maps[1]


def blocks(g):
    rest=set(range(len(g)));out=[]
    while rest:
        seen={min(rest)};rest-=seen;todo=list(seen)
        for i in todo:
            for j in sorted(rest):
                if g[i][j]:seen.add(j);todo.append(j);rest.remove(j)
        out.append(sorted(seen))
    return out


def corner(g,u,v,p,q):
    omitted=sorted(set(u+v)-{p,q});ids=[i for i in range(17)if i not in omitted]
    h=extract(g,ids);psd(h)
    parts=blocks(h);allowed=omitted.copy();out=[];marked=None
    for part in parts:
        b=extract(h,part);ks=nullspace(b);actual=[ids[i]for i in part]
        if ks:
            assert len(ks)==1;k=ks[0]
            if all(x<0 for x in k):k=[-x for x in k]
            assert all(x>0 for x in k)
            out.append(dict(nodes=actual,type='affine',positive_kernel=list(map(str,k))))
        else:
            inv=inverse(b);assert all(x>=0 for row in inv for x in row)
            vals=[inv[i][i]for i in range(len(part))]
            allowed.extend(actual[i]for i,x in enumerate(vals)if x<=2)
            out.append(dict(nodes=actual,type='finite',inverse_diagonals=list(map(str,vals))))
        if p in actual:assert q in actual;marked=actual
    assert marked is not None
    return dict(preferred_corner=[p,q],omitted_black=omitted,corner_nodes=ids,
                blocks=out,possible_contacts=sorted(allowed),marked_component=marked,
                marked_distance=graph_dist(h,[ids.index(p)])[ids.index(q)])


def extended(g,p):
    h=matrix(18)
    for i in range(17):
        for j in range(17):h[i][j]=g[i][j]
    h[17][p]=h[p][17]=-F(1)
    return h


def run():
    ds=json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']
    d=next(x for x in ds if x['index']==25);g,mp=core(d);iv=inverse(g)
    ids=[i for i in range(17)if i not in(0,10)];a=extract(g,ids);pivots=psd(a,True)
    x=extract(g,[0,10],ids);cost=multiply(multiply(x,inverse(a)),transpose(x))
    schur=[[g[i][j]-cost[k][l]for l,j in enumerate([0,10])]for k,i in enumerate([0,10])]
    assert schur==[[-F(1,15),-F(1,10)],[-F(1,10),-F(1,15)]]
    det=prod(pivots)*(schur[0][0]*schur[1][1]-schur[0][1]**2)
    assert det==-2
    end=[1,2,3];other=[11,12,13]
    block=[[F(7,2),4,2],[4,6,3],[2,3,2]]
    assert extract(iv,end)==extract(iv,other)==block
    assert all(x>0 for row in block for x in row)
    roots=[]
    for p in [3,13]:
        r=[-row[p]for row in iv]
        pairings=multiply(g,[[x]for x in r]);assert pairings==[[-F(i==p)]for i in range(17)]
        assert square(g,r)==2
        roots.append(dict(contact_node=p,coordinates=list(map(str,r)),square='2',
                          native_pairings=[str(x[0])for x in pairings]))
    assert iv[3][13]==-2
    center=[-sum(row)for row in iv]
    assert square(g,center)==-F(1001,2)
    assert multiply(g,[[x]for x in center])==[[-1]for _ in range(17)]
    assert -center[3]==-center[13]==-13
    candidates=[]
    for u0,v0 in product(combinations(range(7),2),repeat=2):
        u=list(u0);v=[mp[i]for i in v0]
        if min(graph_dist(g,u)[q]for q in v)<=11:continue
        # Check every adjacent marked corner, not only the preferred one.
        allcorners=[]
        for p,q in product(u,v):
            omitted=set(u+v)-{p,q};ids=[i for i in range(17)if i not in omitted]
            h=extract(g,ids);psd(h)
            distance=graph_dist(h,[ids.index(p)],False)[ids.index(q)]
            assert distance is None or distance>=12
            allcorners.append(dict(marked=[p,q],distance=distance))
        pref=lambda pair: 1 if pair==[1,3]else 2
        p,q=pref(u),mp[pref(list(v0))]
        r=corner(g,u,v,p,q)
        assert set(r['possible_contacts'])<=set(end+other)
        du,dv=graph_dist(g,u),graph_dist(g,v)
        assert all((du[i]if i in end else dv[i])<=1 for i in r['possible_contacts'])
        assert len(r['marked_component'])==15 and r['marked_distance']==12
        variants=[]
        for contact in (3,13):
            h=extended(g,contact)
            white=[i for i in range(17)if i not in u+v]+[17]
            psd(extract(h,white),True)
            cornerids=r['corner_nodes']+[17];c=extract(h,cornerids);psd(c)
            marked=None;kind=None
            for part in blocks(c):
                actual=[cornerids[i]for i in part]
                if p in actual:
                    assert q in actual;marked=actual
                    kind='affine'if nullspace(extract(c,part))else'finite'
            assert marked is not None and len(marked)<=16
            variants.append(dict(contact_node=contact,marked_component=marked,order=len(marked),kind=kind))
            for pp,qq in product(u,v):
                deleted=set(u+v)-{pp,qq};cs=[i for i in range(18)if i not in deleted]
                psd(extract(h,cs))
        # Any unspecified direct cross edge already makes a paid path.
        extra=[]
        for i,j in product(range(7),range(10,17)):
            if g[i][j]:continue
            h=[row.copy()for row in g];h[i][j]=h[j][i]=-F(1)
            if min(graph_dist(h,u)[q]for q in v)>11:extra.append([i,j])
        assert not extra
        candidates.append(dict(left_black=u,right_black=v,all_adjacent_corners=allcorners,
                               preferred=r,one_extension_variants=variants))
    assert len(candidates)==9
    assert {tuple(r['left_black'])for r in candidates}=={(1,2),(1,3),(2,3)}
    sizes=Counter((q['order'],q['kind'])for r in candidates for q in r['one_extension_variants'])
    paths=[f'scripts/{TAG}.py',f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
           'results/fresh_cusp_factor_payment_framework_20261017.md',
           'results/fresh_selective_distance12_dimension136_20261017.md',
           'results/fresh_selective_distance12_dimension136_20261017_review.md']
    return dict(status='PASS_EXACT_DISTANCE11_THREE_NODE_OVERLAP_EXTENSION_LEMMA',
                sources={p:sha(ROOT/p)for p in paths},core=dict(order=17,gram=encode(g),right_copy_map=mp,
                    determinant=str(det),schur=encode(schur),end_inverse_block=encode(block),
                    cross_end_inverse_block=encode(extract(iv,end,other))),
                unique_extensions=roots,mutual_product='-2',common_strict_center_square='-1001/2',
                marked_contexts=candidates,one_extension_component_counts=[dict(order=n,kind=k,count=c)for(n,k),c in sorted(sizes.items())],
                marked_component_order_bound=16,
                scope='Conditional exact local lemma for twoL25 cores sharing three terminal nodes in reverse, with no marked white path of length<=11. It does not classify other distance11 overlaps, prove a new dimension bound, or assert finite volume of these local root sets.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=json.loads(json.dumps(run()));target=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(target.read_text())
    else:target.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'9contexts/2uniqueextensions/mutual-2/component<=16',flush=True)

