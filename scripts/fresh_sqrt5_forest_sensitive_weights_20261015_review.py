#!/usr/bin/env python3
"""Independent graph/Fraction replay of the A1-sensitive partial rule.
No project imports; no LP optimizer; sealing and the finite-type theorem are
explicit frozen proof dependencies, not silently recomputed here.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_forest_sensitive_weights_20261015_review'

def graph(nodes,labels):
    g={i:{}for i in nodes}
    for(i,j),m in labels.items():
        if i in g and j in g and m!=2:g[i][j]=m;g[j][i]=m
    return g

def pieces(g):
    left=set(g);out=[]
    while left:
        seen={min(left)};todo=list(seen)
        while todo:
            i=todo.pop()
            for j in g[i]:
                if j not in seen:seen.add(j);todo.append(j)
        left-=seen;out.append({i:{j:m for j,m in g[i].items()if j in seen}for i in seen})
    return out

def code(g,marks=()):
    # Center-rooted tree invariant, different from producer minimum over all roots.
    assert len(pieces(g))==1 and sum(map(len,g.values()))==2*len(g)-2
    marks=set(marks);alive=set(g)
    while len(alive)>2:
        leaves={i for i in alive if sum(j in alive for j in g[i])<=1}
        assert leaves;alive-=leaves
    def rec(i,parent):
        return (i in marks,tuple(sorted((m,rec(j,i))for j,m in g[i].items()if j!=parent)))
    if len(alive)==1:return (1,rec(next(iter(alive)),None))
    a,b=sorted(alive);return(2,g[a][b],tuple(sorted((rec(a,b),rec(b,a)))))

def decode(k):
    g={};marks=[]
    def visit(node):
        i=len(g);g[i]={}
        if node[0]:marks.append(i)
        for m,child in node[1]:
            j=visit(child);g[i][j]=m;g[j][i]=m
        return i
    visit(json.loads(k));return code(g,marks)

def distances(g,a):
    d={a:0};todo=[a]
    for i in todo:
        for j in g[i]:
            if j not in d:d[j]=d[i]+1;todo.append(j)
    return d

def family(g):
    n=len(g);es=[m for i in g for j,m in g[i].items()if i<j]
    if n==1:return 'A1'
    if n==2:return 'I2(%d)'%es[0]
    fork=[i for i in g if len(g[i])==3]
    if fork:
        assert len(fork)==1
        b=fork[0];arms=[]
        for j in g[b]:
            p=b;k=j;r=1
            while len(g[k])==2:
                z=next(t for t in g[k]if t!=p);p,k=k,z;r+=1
            arms.append(r)
        return ('D'if sorted(arms)[:2]==[1,1]else'E')+str(n)
    if 5 in es:return 'H'+str(n)
    if 4 in es:
        a,b=next((i,j)for i in g for j,m in g[i].items()if m==4)
        return 'B'+str(n)if len(g[a])==1 or len(g[b])==1 else'F4'
    return'A'+str(n)

def base(g,a,b):
    d=distances(g,a).get(b)
    if d is None:return F(0)
    return F(int(d<=2)+int(family(g)=='H4'and(d==3 or g[a].get(b)==5)))

def models():
    out=[]
    def add(name,n,edges):out.append((name,graph(range(n),{tuple(sorted((a,b))):m for a,b,m in edges})))
    for r in range(1,30):add('A'+str(r),r,[(i,i+1,3)for i in range(r-1)])
    for r in range(2,30):add('B'+str(r),r,[(i,i+1,4 if i==0 else 3)for i in range(r-1)])
    for r in range(4,30):add('D'+str(r),r,[(i,i+1,3)for i in range(r-2)]+[(1,r-1,3)])
    for r in(6,7,8):add('E'+str(r),r,[(i,i+1,3)for i in range(r-2)]+[(2,r-1,3)])
    add('F4',4,[(0,1,3),(1,2,4),(2,3,3)])
    for r in(3,4):add('H'+str(r),r,[(i,i+1,5 if i==0 else 3)for i in range(r-1)])
    for m in(5,6,10):add('I2(%d)'%m,2,[(0,1,m)])
    assert len(out)==92;return out

def source(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    tag='fresh_sqrt5_forest_sensitive_weights_20261015'
    primary=source('results/'+tag+'.json')
    old=source('results/fresh_sqrt5_sealed_weights_20261015.json')
    unions=source('results/fresh_sqrt5_quad_exhaustion_20261015.json')
    contacts=source('results/fresh_sqrt5_quad_contacts_20261015.json')
    tri=source('results/fresh_sqrt5_triangle_weights_20261015.json')
    # Read literal pair lists, never using the producer canonical-key algorithm.
    bonus={}
    for row in old['component_rows']:
        g=graph(row['nodes'],dict(zip(combinations(row['nodes'],2),row['labels'])))
        for orbit in row['pair_orbits']:
            vals={code(g,p)for p in orbit['pairs']};assert len(vals)==1
            k=vals.pop();v=F(orbit['weight']);assert v>=0
            if k in bonus:assert bonus[k]==v
            bonus[k]=v
            assert decode(orbit['key'])==k
    selected=[r for r in contacts['records']if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and r['singular']and not r['accepted']]
    assert len(selected)==112
    def corners(ci):
        row=unions['labelled_union_supports'][next(r['union_index']for r in selected if r['class_index']==ci)]
        es=dict(zip(combinations(range(row['nodes']),2),row['labels']))
        ans=[]
        for a,b in product((0,2),(1,3)):
            comps=pieces(graph(row['white']+[a,b],es));c=next(g for g in comps if a in g)
            ans.append((a,b,c,comps))
        return ans
    pure={code(c,(a,b))for ci in(535,563)for a,b,c,_ in corners(ci)if b in c}
    half={code(c,(a,b))for ci in(785,860)for a,b,c,_ in corners(ci)if b in c and family(c)=='D7'}
    assert len(pure)==len(half)==2 and not pure&half
    assert pure=={decode(k)for k in primary['pure_d7_orbits']}
    assert half=={decode(k)for k in primary['half_d7_orbits']}
    def weight(g,a,b,flag):
        k=code(g,(a,b));v=bonus.get(k,F(0))
        if family(g)=='D7'and not flag:v=F(1)if k in pure else F(1,2)if k in half else F(0)
        assert v>=0;return base(g,a,b)+v
    facedata=[];minima=[None,None]
    for r in selected:
        ci=r['class_index'];expected=next(v for v in primary['sealed_quad_rows']if v['class_index']==ci)
        sums=[F(0),F(0)];rows=[]
        assert not any(t['status']=='PASS_EXTRA_WHITE_NECESSARY_TESTS'for t in r['trials'])
        for a,b,c,cs in corners(ci):
            flag=any(len(g)==1 for g in cs)
            vals=[weight(c,a,b,flag or outside)if b in c else F(0)for outside in(False,True)]
            ex=next(z for z in expected['corners']if z['marked']==[a,b])
            assert vals==list(map(F,ex['weights']))and flag==ex['core_has_a1']
            assert sorted(map(family,cs))==sorted(ex['core_forest'])
            sums=[x+y for x,y in zip(sums,vals)]
            rows.append(dict(marked=[a,b],full_core_types=sorted(map(family,cs)),has_core_A1=flag,weights=list(map(str,vals))))
        assert sums==list(map(F,expected['sums_without_with_external_a1']))and min(sums)>=1
        for i in range(2):minima[i]=sums[i]if minima[i]is None else min(minima[i],sums[i])
        facedata.append(dict(class_index=ci,sums=list(map(str,sums)),corners=rows))
    componentdata=[]
    for name,g in models():
        ex=next(r for r in primary['component_rows']if r['type']==name);sums=[F(0),F(0)]
        for a,b in combinations(g,2):
            vals=[weight(g,a,b,bool(i))for i in range(2)]
            saved=next(z for z in ex['pair_weights']if z['marked']==[a,b]);assert vals==list(map(F,saved['weights']))
            sums=[x+y for x,y in zip(sums,vals)]
        assert sums==list(map(F,ex['totals_without_with_a1']))
        componentdata.append(dict(type=name,rank=len(g),totals=list(map(str,sums))))
    # Independent unbounded coin-change over component TYPES, with a separate
    # state tracking actual A1 occurrence, instead of reserving one A1 first.
    dps=[]
    for context in(0,1):
        dp={(0,False):(F(0),())}
        for r in componentdata:
            for n in range(r['rank'],30):
                for oldflag in(False,True):
                    oldstate=dp.get((n-r['rank'],oldflag))
                    if oldstate is None:continue
                    state=(n,oldflag or r['type']=='A1');value=(oldstate[0]+F(r['totals'][context]),tuple(sorted(oldstate[1]+(r['type'],))))
                    if state not in dp or value>dp[state]:dp[state]=value
        dps.append(dp)
    no=dps[0][29,False];yes=dps[1][29,True]
    assert no[0]==F(159,2)and yes[0]==84
    assert no[0]==F(primary['rank29_without_a1']['maximum'])and yes[0]==F(primary['rank29_with_a1']['maximum'])
    assert list(no[1])==primary['rank29_without_a1']['forest']and list(yes[1])==sorted(primary['rank29_with_a1']['forest'])
    # Check each saved dynamic-program value, not just its terminal maximum.
    assert primary['no_a1_rank_maxima']==[None if(n,False)not in dps[0]else str(dps[0][n,False][0])for n in range(30)]
    assert primary['with_a1_tail_rank_maxima']==[str(max(v[0]for(k,_),v in dps[1].items()if k==n))for n in range(29)]
    triangles=[]
    for idx,row in enumerate(tri['markings']):
        es=dict(zip(combinations(range(row['order']),2),row['labels']));s=F(0)
        for saved in row['corners']:
            a,b=saved['marked'];cs=pieces(graph(row['white']+[a,b],es));g=next(c for c in cs if a in c)
            v=base(g,a,b)if b in g else F(0);assert v==saved['h4_only_weight'];s+=v
        assert s>=2;triangles.append(str(s))
    paths=[f'scripts/{TAG}.py',f'scripts/{tag}.py',f'results/{tag}.json']
    for stem in('fresh_sqrt5_sealed_weights_20261015','fresh_sqrt5_quad_exhaustion_20261015','fresh_sqrt5_quad_contacts_20261015','fresh_sqrt5_quad_contacts_20261015_review','fresh_sqrt5_triangle_weights_20261015','fresh_sqrt5_triangle_weights_20261015_review'):
        paths.append('results/'+stem+'.json')
    # Direct and inherited manifests remain exact; no outputs are changed.
    for p in paths:
        if p.endswith('.json'):
            for q,h in source(p).get('sources',{}).items():assert sha(q)==h,(p,q)
    return dict(status='PASS_INDEPENDENT_EXACT_FOREST_SENSITIVE_PARTIAL_RULE',sources={p:sha(p)for p in paths},
        arithmetic='Standalone Fraction arithmetic and center-rooted marked-tree invariants; no producer or project imports.',
        sealed_quadrilaterals=112,external_A1_cases=224,minima=list(map(str,minima)),face_rows=facedata,
        finite_component_types=92,component_totals=componentdata,rank29_without_A1=dict(maximum=str(no[0]),forest=list(no[1])),
        rank29_with_A1=dict(maximum=str(yes[0]),forest=list(yes[1])),triangle_markings=76,triangle_sums=triangles,
        inherited='Conditional unpaid first-external-white sealing and complete arithmetic finite-type/triangle-extension premises are frozen independently reviewed dependencies; no new all-contact enumeration here.',
        unresolved='130 finite and24 symbolic dotted unpaid principal classes remain. This is a partial face rule with a valid rank29 vertex budget, not a dimension theorem.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['rank29_without_A1'],out['rank29_with_A1'])
