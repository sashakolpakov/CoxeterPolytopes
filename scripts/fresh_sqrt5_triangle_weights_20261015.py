#!/usr/bin/env python3
"""All Qsqrt5 Lanner triangle markings and a stable spherical weight rule.

The proof for arbitrary white extensions is in the accompanying report.  The
all-rank component formulas, not the finite rank guard, establish the bound.
No producer arithmetic or graph code is imported.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_triangle_weights_20261015'
CAT = 'results/fresh_sqrt5_lanner_catalogue_20261015.json'
CAT_REVIEW = 'results/fresh_sqrt5_lanner_catalogue_20261015_review.json'

def q(a=0, b=0): return F(a), F(b)
def add(x, y): return x[0]+y[0], x[1]+y[1]
def neg(x): return -x[0], -x[1]
def mul(x, y):
    a,b=x; c,d=y
    return a*c+b*d, a*d+b*c+b*d
def div(x, y):
    c,d=y; norm=c*c+c*d-d*d
    assert norm
    return mul(x, ((c+d)/norm, -d/norm))
def sign(x):
    a,b=x; u=2*a+b
    if not b: return (u>0)-(u<0)
    if not u: return (b>0)-(b<0)
    if (u>0)==(b>0): return (u>0)-(u<0)
    d=u*u-5*b*b
    assert d
    return ((u>0)-(u<0))*((d>0)-(d<0))

def graph(nodes, edges):
    g={i:{} for i in nodes}
    for a,b,m in edges:
        assert a!=b and m in (3,4,5,6,10)
        g[a][b]=m; g[b][a]=m
    return g
def restrict(g, nodes):
    nodes=set(nodes)
    return {i:{j:m for j,m in g[i].items() if j in nodes} for i in nodes}
def components(g):
    left=set(g); out=[]
    while left:
        todo=[min(left)]; found=set()
        while todo:
            i=todo.pop()
            if i in found: continue
            found.add(i); todo.extend(g[i])
        left-=found; out.append(sorted(found))
    return out
def distance(g, a, b):
    todo=deque([(a,0)]); seen={a}
    while todo:
        i,d=todo.popleft()
        if i==b: return d
        for j in g[i]:
            if j not in seen: seen.add(j); todo.append((j,d+1))
    return None
def walk(g, start, previous=None):
    out=[start]
    while True:
        nxt=[j for j in g[start] if j!=previous]
        if not nxt: return out
        assert len(nxt)==1
        previous,start=start,nxt[0]; out.append(start)
def describe(g):
    """Intrinsic complete spherical classification, within the field alphabet."""
    r=len(g); edges=[(i,j,m) for i in g for j,m in g[i].items() if i<j]
    assert len(edges)==r-1
    if r==1: return 'A1', {}
    if r==2:
        m=edges[0][2]
        return f'I2({m})', {tuple(sorted((edges[0][0],edges[0][1]))):1} if m==5 else {}
    branches=[i for i in g if len(g[i])==3]
    if not branches:
        ends=[i for i in g if len(g[i])==1]
        assert len(ends)==2 and max(map(len,g.values()))<=2
        path=walk(g,min(ends)); assert len(path)==r
        labels=[g[i][j] for i,j in zip(path,path[1:])]
        if all(m==3 for m in labels): return f'A{r}', {}
        if labels.count(4)==1 and all(m in (3,4) for m in labels):
            k=labels.index(4)
            if r==4 and k==1: return 'F4', {}
            assert k in (0,r-2)
            return f'B{r}', {}
        assert r in (3,4) and labels.count(5)==1 and all(m in (3,5) for m in labels)
        k=labels.index(5); assert k in (0,r-2)
        extras={tuple(sorted((path[k],path[k+1]))):1}
        if r==4: extras[tuple(sorted((path[0],path[-1])))]=1
        return f'H{r}', extras
    assert len(branches)==1 and max(map(len,g.values()))==3
    assert all(m==3 for _,_,m in edges)
    c=branches[0]; lengths=sorted(len(walk(g,j,c)) for j in g[c])
    if lengths[:2]==[1,1]:
        assert lengths==[1,1,r-3]
        return f'D{r}', {}
    assert lengths[:2]==[1,2] and lengths[2] in (2,3,4)
    return f'E{r}', {}
def pd_pivots(g):
    sq={3:q(1),4:q(2),5:q(1,1),6:q(3),10:q(2,1)}
    d={i:q(2) for i in g}; adj={i:set(g[i]) for i in g}; out=[]
    while adj:
        leaves=[i for i in adj if len(adj[i])<=1]
        assert leaves, 'finite corner must be a forest'
        i=min(leaves); p=d[i]; assert sign(p)>0
        out.append([str(z) for z in p])
        if adj[i]:
            j=next(iter(adj[i])); d[j]=add(d[j],neg(div(sq[g[i][j]],p))); adj[j].remove(i)
        del adj[i]
    return out
def weights(g):
    w={pair:int(distance(g,*pair) in (1,2)) for pair in combinations(sorted(g),2)}
    types=[]
    for ids in components(g):
        name,extra=describe(restrict(g,ids)); types.append(name)
        for pair,value in extra.items(): w[pair]+=value
    return w,types
def path(r, first=3):
    return graph(range(r),[(i,i+1,first if i==0 else 3) for i in range(r-1)])
def arms(lengths):
    edges=[]; k=1
    for length in lengths:
        old=0
        for _ in range(length): edges.append((old,k,3)); old=k; k+=1
    return graph(range(k),edges)
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    cat=json.loads((ROOT/CAT).read_text())
    independent=json.loads((ROOT/CAT_REVIEW).read_text())
    assert cat['status'].startswith('PASS') and independent['status'].startswith('PASS')
    assert independent['sources'][CAT]==sha(CAT)
    lanners=[dict(order=3,labels=r['labels']) for r in cat['order_three']['accepted']]
    lanners += [dict(order=group['order'],labels=r['labels']) for group in cat['higher_orders'] for r in group['accepted_classes']]
    assert Counter(r['order'] for r in lanners)=={3:12,4:6,5:4}
    marks=[]; histogram=Counter(); deficiencies=[]
    for index,row in enumerate(lanners):
        n=row['order']; es=list(zip(combinations(range(n),2),row['labels']))
        g=graph(range(n),[(i,j,m) for (i,j),m in es if m!=2])
        for black in combinations(range(n),3):
            corners=[]
            for omitted in black:
                h=restrict(g,set(g)-{omitted}); pair=tuple(x for x in black if x!=omitted)
                w,types=weights(h); d=distance(h,*pair)
                reduced=int(d in (1,2))
                for ids in components(h):
                    name,extra=describe(restrict(h,ids))
                    if name=='H4': reduced+=extra.get(pair,0)
                corner=dict(omitted=omitted,marked=pair,distance=d,base=int(d in (1,2)),
                            weight=w[pair],h4_only_weight=reduced,components=types,positive_pivots=pd_pivots(h))
                corners.append(corner)
            base=sum(c['base'] for c in corners); total=sum(c['weight'] for c in corners)
            assert total>=2
            assert sum(c['h4_only_weight'] for c in corners)>=2
            record=dict(lanner_index=index,order=n,labels=row['labels'],black=black,
                        white=sorted(set(g)-set(black)),corners=corners,base_total=base,total=total)
            marks.append(record); histogram[(n,base,total)]+=1
            if base<2: deficiencies.append(len(marks)-1)
    assert len(marks)==76 and len(deficiencies)==7

    models=[]
    for r in range(1,101):
        models.append((f'A{r}',path(r)))
        if r>=2: models.append((f'B{r}',path(r,4)))
        if r>=4: models.append((f'D{r}',arms([1,1,r-3])))
    models += [('E6',arms([1,2,2])),('E7',arms([1,2,3])),('E8',arms([1,2,4])),
               ('F4',graph(range(4),[(0,1,3),(1,2,4),(2,3,3)])),
               ('I2(5)',path(2,5)),('I2(6)',path(2,6)),('I2(10)',path(2,10)),
               ('H3',path(3,5)),('H4',path(4,5))]
    finite=[]
    for name,g in models:
        rank=len(g); w,types=weights(g); value=sum(w.values()); deficit=2*rank-value
        expected=0 if rank==1 else 2*rank-3
        if name[0] in ('D','E'): expected=2*rank-2
        if name=='I2(5)': expected=2
        if name=='H3': expected=4
        if name=='H4': expected=7
        assert value==expected and deficit>=1
        assert deficit>=2 or name=='H4'
        finite.append(dict(type=name,rank=rank,total=value,deficit=deficit))
        # Exact positivity guards for exceptional finite diagrams and low ranks.
        if rank<=8: pd_pivots(g)
    best=[0]
    for n in range(1,101):
        best.append(max(best[n-r['rank']]+r['total'] for r in finite if r['rank']<=n))
        if n>=5: assert best[n]==2*n-2
    assert best[:6]==[0,0,2,4,7,8]

    # Enumerate all induced embeddings among the possible finite5-components.
    # This independently guards preservation of every correction, including
    # I2(5)->H3->H4. General no-other-extension is the spherical classification.
    five=[('I2(5)',path(2,5)),('H3',path(3,5)),('H4',path(4,5))]
    extensions=[]
    for source,sg in five:
        sw,_=weights(sg)
        for target,tg in five:
            if len(sg)>len(tg): continue
            count=0
            for images in permutations(sorted(tg),len(sg)):
                mp=dict(zip(sorted(sg),images))
                if any(sg[i].get(j,2)!=tg[mp[i]].get(mp[j],2) for i,j in combinations(sg,2)): continue
                count+=1; tw,_=weights(tg)
                for pair,value in sw.items(): assert tw[tuple(sorted(mp[i] for i in pair))]>=value
            extensions.append(dict(source=source,target=target,induced_embeddings=count))
    # Relabelling guard for intrinsic extras; the proof covers all automorphisms.
    for name,g in five:
        original,_=weights(g)
        for p in permutations(sorted(g)):
            mp=dict(zip(sorted(g),p)); h=graph(p,[(mp[i],mp[j],m) for i in g for j,m in g[i].items() if i<j])
            changed,_=weights(h)
            assert {tuple(sorted(mp[i] for i in pair)):value for pair,value in original.items()}==changed
    paths=[f'scripts/{TAG}.py',CAT,'scripts/fresh_sqrt5_lanner_catalogue_20261015.py',CAT_REVIEW,
           'scripts/fresh_sqrt5_lanner_catalogue_20261015_review.py']
    return json.loads(json.dumps(dict(
        status='PASS_COMPLETE_QSQRT5_TRIANGLE_WEIGHT_RULE_AND_ALL_RANK_VERTEX_BOUND',
        sources={p:sha(p) for p in paths},
        weight_rule=dict(base='1 for graph distance1 or2; otherwise0',
                         five_bond='add1 on the unique5-edge in an entire I2(5),H3,H4 component',
                         h4='add1 on the endpoint pair of an entire H4 component',
                         i2_10='no correction'),
        lanner_count_by_order={str(n):sum(r['order']==n for r in lanners) for n in (3,4,5)},
        markings=marks,base_deficient_marking_indices=deficiencies,
        histogram=[dict(order=n,base=b,total=t,count=c) for (n,b,t),c in sorted(histogram.items())],
        finite_component_totals=finite,rank_partition_maxima=best,five_component_extensions=extensions,
        all_rank_proof=dict(path='q2=2r-3 for r>=2',fork='q2=2r-2 for D/E',
                            unique_deficit_below_two='H4, rank4, total7, deficit1',
                            conclusion='rank n>=5: total<=2n-2; single component has deficit>=2, or at least two components each have deficit>=1'),
        arbitrary_white_extension_proof='Spherical graphs are forests, so existing connected paths retain their length. A finite5-component is only I2(5),H3,H4; every5-edge correction persists and any H4 subdiagram is an entire component. Extra connections can only add nonnegative weight to formerly disconnected marked pairs.',
        smaller_triangle_rule='The same triangle conclusion holds with only the two H4 bonuses; I2(5) and H3 bonuses are optional for triangles. This is checked for every marking and is stable because H4 is an entire component.',
        scope='Every triangular2-face of a compact arithmetic Qsqrt5 Coxeter polytope has weight>=2, with arbitrary white extensions and orthogonal components. Vertex total<=2n-2 for n>=5. No quadrilateral coverage or dimension bound claimed.'
    )))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify',action='store_true'); args=parser.parse_args()
    result=run(); output=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert result==json.loads(output.read_text()); print('PASS saved triangle-weight replay')
    else:
        output.write_text(json.dumps(result,indent=2)+'\n'); print(result['status'])
    print('76 markings; 7 q2-deficient; all triangle weights>=2; vertex<=2n-2 for n>=5')
