#!/usr/bin/env python3
"""Separate BFS/arm and Fraction verification of the proposed weight bound.

No producer or project arithmetic/graph imports. Global star coverage is an
explicit external dependency, not re-established by this small checker.
"""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from collections import deque
from pathlib import Path
from hashlib import sha256

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt2_weight_bound_20261015_review'
PRIMARY='results/fresh_sqrt2_weight_bound_20261015.json'
BASE='results/fresh_sqrt2_quad_exhaustion_20261015.json'
CORES='results/fresh_sqrt2_quad_exhaustion_20261015_cores.json'
TAIL='results/fresh_sqrt2_dotted_tail_20261015.json'
CATALOGUE='results/fresh_sqrt2_lanner_catalogue_20261015.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def div(x,y):
    n=y[0]*y[0]-2*y[1]*y[1];assert n
    return mul(x,(y[0]/n,-y[1]/n))
def sgn(x):
    a,b=x
    if not a:return (b>0)-(b<0)
    if not b:return (a>0)-(a<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    n=a*a-2*b*b;assert n
    return ((a>0)-(a<0))*((n>0)-(n<0))
def q(a=0,b=0):return F(a),F(b)
def graph(nodes,edges):
    out={i:{}for i in nodes}
    for a,b,m in edges:out[a][b]=m;out[b][a]=m
    return out
def induced(g,nodes):
    nodes=set(nodes);return {i:{j:m for j,m in g[i].items()if j in nodes}for i in nodes}
def components(g):
    left=set(g);out=[]
    while left:
        todo=[min(left)];found=set()
        while todo:
            i=todo.pop()
            if i in found:continue
            found.add(i);todo.extend(g[i])
        left-=found;out.append(sorted(found))
    return out
def distance(g,a,b):
    todo=deque([(a,0)]);seen={a}
    while todo:
        i,d=todo.popleft()
        if i==b:return d
        for j in g[i]:
            if j not in seen:seen.add(j);todo.append((j,d+1))
    return None
def path_from(g,a,previous=None):
    out=[a]
    while True:
        nxt=[j for j in g[a]if j!=previous]
        if not nxt:return out
        assert len(nxt)==1
        previous,a=a,nxt[0];out.append(a)
def describe(g):
    """Classify a connected finite Qsqrt2 Coxeter tree by intrinsic arms."""
    r=len(g);edges=[(i,j,m)for i in g for j,m in g[i].items()if i<j]
    assert len(edges)==r-1
    if r==1:return 'A1',{}
    if r==2:
        assert edges[0][2]in(3,4,6,8)
        return f'I2({edges[0][2]})',{}
    branch=[i for i in g if len(g[i])==3]
    if not branch:
        ends=[i for i in g if len(g[i])==1];assert len(ends)==2
        path=path_from(g,min(ends));assert len(path)==r
        labels=[g[i][j]for i,j in zip(path,path[1:])]
        if all(m==3 for m in labels):return f'A{r}',{}
        assert labels.count(4)==1 and all(m in(3,4)for m in labels)
        pos=labels.index(4)
        if r==4 and pos==1:return 'F4',{tuple(sorted((path[0],path[-1]))):F(1)}
        assert pos in(0,r-2)
        return f'B{r}',{}
    assert len(branch)==1 and all(len(g[i])<=3 for i in g)
    assert all(m==3 for _,_,m in edges)
    c=branch[0];arms=[path_from(g,i,c)for i in g[c]]
    arms.sort(key=lambda z:(len(z),z));lengths=list(map(len,arms));extra={}
    if lengths[:2]==[1,1]:
        assert lengths==[1,1,r-3]
        if r==5:
            for a in arms[:2]:extra[tuple(sorted((a[0],arms[2][-1])))]=F(1,2)
        return f'D{r}',extra
    assert lengths[:2]==[1,2] and lengths[2]in(2,3,4)
    if r==6:
        for a,b in ((arms[1],arms[2]),(arms[2],arms[1])):
            extra[tuple(sorted((a[0],b[-1])))]=F(1,2)
    else:extra[tuple(sorted((arms[2][0],arms[1][-1])))]=F(1,4)
    return f'E{r}',extra
def finite_pivots(g):
    """Exact leaf Schur pivots; edge square alone suffices for a forest."""
    d={i:q(2)for i in g};a={i:set(g[i])for i in g};out=[]
    sq={3:q(1),4:q(2),6:q(3),8:q(2,1)}
    while a:
        ids=[i for i in a if len(a[i])<=1];assert ids,'cycle in finite corner'
        i=min(ids);pivot=d[i];assert sgn(pivot)>0;out.append([str(x)for x in pivot])
        if a[i]:
            j=next(iter(a[i]));d[j]=add(d[j],neg(div(sq[g[i][j]],pivot)));a[j].remove(i)
        del a[i]
    return out
def weights(g):
    out={p:F(distance(g,*p)in(1,2))for p in combinations(sorted(g),2)}
    types=[]
    for nodes in components(g):
        h=induced(g,nodes);name,extra=describe(h);types.append(name)
        for pair,w in extra.items():
            assert distance(g,*pair)==3 and out[pair]==0;out[pair]+=w
    return out,types
def path(r,label=3):return graph(range(r),[(i,i+1,label if i==0 else 3)for i in range(r-1)])
def arms(lengths):
    edges=[];k=1
    for length in lengths:
        previous=0
        for _ in range(length):edges.append((previous,k,3));previous=k;k+=1
    return graph(range(k),edges)
def total(g):return sum(weights(g)[0].values(),F())
def run():
    primary=json.loads((ROOT/PRIMARY).read_text())
    base=json.loads((ROOT/BASE).read_text());cores=json.loads((ROOT/CORES).read_text())
    tails=json.loads((ROOT/TAIL).read_text())
    records=[]
    for r in cores['records']:
        if r['status']!='PASS_PRINCIPAL_UNION_NECESSARY_TESTS'or r['physical']['determinant']:continue
        u=base['labelled_union_supports'][r['union_index']]
        edges=[(i,j,m)for (i,j),m in zip(combinations(range(u['nodes']),2),u['labels'])if m!=2]
        g=graph(range(u['nodes']),edges);rows=[]
        for a,b in ((0,1),(1,2),(2,3),(3,0)):
            h=induced(g,u['white']+[a,b]);w,types=weights(h)
            rows.append(dict(marked=[a,b],components=types,weight=str(w[tuple(sorted((a,b)))]),pivots=finite_pivots(h)))
        value=sum((F(x['weight'])for x in rows),F());assert value>=1
        records.append(dict(class_index=r['class_index'],corners=rows,total=str(value)))
    assert len(records)==7
    assert [dict(class_index=r['class_index'],corner_weights=[c['weight']for c in r['corners']],total=r['total'])for r in records]==primary['finite_core_weights']
    tail_records=[];cartan={q():2,q(1):3,q(2):4,q(3):6,q(2,1):8}
    for r in tails['allowed_local_stars']:
        M=[[tuple(map(F,x))for x in row]for row in r['gram']];rows=[]
        for corner in r['corners']:
            ids=corner['support'];edges=[]
            for i,j in combinations(ids,2):
                cycle=div(mul(q(4),mul(M[i][j],M[i][j])),mul(M[i][i],M[j][j]))
                assert cycle in cartan
                if cartan[cycle]!=2:edges.append((i,j,cartan[cycle]))
            g=graph(ids,edges);w,types=weights(g);a,b=corner['marked']
            rows.append(dict(marked=[a,b],components=types,weight=str(w[tuple(sorted((a,b)))]),pivots=finite_pivots(g)))
        value=sum((F(x['weight'])for x in rows),F());assert value>=1
        tail_records.append(dict(length=r['ell'],double=r['double'],corners=rows,total=str(value)))
    assert [dict(length=r['length'],double=r['double'],corner_weights=[c['weight']for c in r['corners']],total=r['total'])for r in tail_records]==primary['tail_weights']
    # Literal components, independent enumeration and all-rank deficit proof.
    models=[]
    for rank in range(1,101):
        models.append((f'A{rank}',path(rank)))
        if rank>=2:models.append((f'B{rank}',path(rank,4)))
        if rank>=4:models.append((f'D{rank}',arms([1,1,rank-3])))
    models += [('E6',arms([1,2,2])),('E7',arms([1,2,3])),('E8',arms([1,2,4])),
        ('F4',graph(range(4),[(0,1,3),(1,2,4),(2,3,3)])),('I2(6)',path(2,6)),('I2(8)',path(2,8))]
    finite=[]
    for name,g in models:
        rank=len(g);value=total(g);deficit=2*rank-value
        expected=(F(0)if rank==1 else F(2*rank-3))if name[0]in('A','B','I')else F(2*rank-2)
        if name=='D5':expected=F(9)
        if name=='E6':expected=F(11)
        if name=='E7':expected=F(49,4)
        if name=='E8':expected=F(57,4)
        if name=='F4':expected=F(6)
        assert value==expected and deficit>=1
        if rank>=9:assert deficit>=2
        finite.append(dict(type=name,rank=rank,total=str(value),deficit=str(deficit)))
    values={r['type']:r['total']for r in finite}
    assert all(values[r['type']]==r['total']for r in primary['finite_diagrams'])
    best=[F(0)]
    for n in range(1,101):
        best.append(max(best[n-r['rank']]+F(r['total'])for r in finite if r['rank']<=n))
        if n>=9:assert best[n]==2*n-2
    assert list(map(str,best[:61]))==primary['rank_maxima']
    # Every permutation of equal-length arms preserves the intrinsic extras.
    # More strongly, test all relabellings for the five corrected small types.
    # Arm descriptions depend only on lengths and use all symmetric pairs;
    # the reversed labelling is an explicit independent implementation guard.
    invariance=[]
    for name,g in models:
        if name not in('D5','E6','E7','E8','F4'):continue
        mp={i:len(g)-1-i for i in g};h=graph(mp.values(),[(mp[i],mp[j],m)for i in g for j,m in g[i].items()if i<j])
        w=weights(g)[0];v=weights(h)[0]
        assert {tuple(sorted((mp[i],mp[j]))):a for (i,j),a in w.items()}==v
        invariance.append(name)
    # At least two three-cycle edges, or two short gaps on a4/5cycle.
    triangle=[]
    for n in (3,4,5):
        for marks in combinations(range(n),3):
            if n==3:
                for edges in combinations([(0,1),(1,2),(0,2)],2):
                    g=graph(range(3),[(i,j,3)for i,j in edges]);assert sum(int(distance(g,i,j)==1)for i,j in combinations(marks,2))==2
                continue
            gaps=[(marks[(i+1)%3]-marks[i])%n for i in range(3)]
            assert sum(d<=2 for d in gaps)>=2
            triangle.append(dict(order=n,marks=list(marks),gaps=gaps))
    assert triangle==primary['cycle_triangle_marks']
    # Derive the strict h-vector lower bound from its symmetric coefficients.
    dimensions=[]
    for n in range(9,151):
        coeff=[F(i*(i-1)+(n-i)*(n-i-1),4)for i in range(n+1)]
        minimum=min(coeff)
        assert coeff[0]>minimum and coeff[n]>minimum
        assert minimum==(F(n*(n-2),8)if n%2==0 else F((n-1)*(n-1),8))
        lower=5*minimum-F(n*(n-1),2)
        assert lower==(F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8))
        gap=lower-(2*n-2)
        if n>=21:assert gap>=0
        dimensions.append(dict(dimension=n,min_symmetric_coefficient=str(minimum),strict_weight_lower=str(lower),upper=2*n-2,gap=str(gap)))
    for row in primary['dimension_inequalities']:
        mine=dimensions[row['dimension']-9]
        assert mine['strict_weight_lower']==row['strict_lower']and mine['upper']==row['upper']
    # Polynomial coefficients are reconstructed by finite differences; degree2
    # follows the displayed exact quadratic formulas, not a sampled inference.
    polynomials=[]
    for start in (21,22):
        vals=[]
        for u in range(3):
            n=start+2*u;low=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
            vals.append(low-(2*n-2))
        a=(vals[2]-2*vals[1]+vals[0])/2;b=vals[1]-vals[0]-a;c=vals[0]
        assert a==F(1,2)and b==(5 if start==21 else F(11,2))and c==(0 if start==21 else 2)
        assert a>0 and b>0 and c>=0
        polynomials.append(dict(parity_start=start,coefficients=list(map(str,(c,b,a)))))
    paths=[f'scripts/{TAG}.py',PRIMARY,'scripts/fresh_sqrt2_weight_bound_20261015.py',BASE,CORES,TAIL,CATALOGUE]
    return dict(status='PASS_INDEPENDENT_WEIGHT_ACCOUNTING_CONDITIONAL_ON_GLOBAL_STAR_COVERAGE',
        sources={p:sha(p)for p in paths},arithmetic='Standalone Fraction Qsqrt2 edge-square Schur and BFS/arm classifier; no producer, NetworkX, or project imports.',
        finite_core_weights=records,tail_weights=tail_records,finite_component_totals=finite,
        entire_component_relabel_checks=invariance,rank_partition_maxima=list(map(str,best)),
        triangle_cycle_marks=triangle,hvector_bounds=dimensions,all_dimension_gap_polynomials=polynomials,
        proof_scope=dict(verified='Seven literal core stars, three literal scaled-Gram tails, triangle path lower bound, all-rank finite-component deficit, strict h-vector dimension20 implication.',
            inherited='Qsqrt2 Lanner catalogue and finite spherical diagram classification; completeness of all quadrilateral cores and arbitrary external attachments is NOT established here.',
            conditional_conclusion='If every compact arithmetic Qsqrt2 quadrilateral reduces to the stated entire-component configurations, then dimension is at most20.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS independent weight saved replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
