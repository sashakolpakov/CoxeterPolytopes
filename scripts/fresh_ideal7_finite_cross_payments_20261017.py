#!/usr/bin/env python3
"""Selective order17 cusp payments and finite cross-component coupons.

This package supplies local transfer and all-rank graph/finite cost inputs.
The ambient dimension comparison is a separate certificate.
"""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_selective_distance9_dimension108_20261017_review as R
import fresh_affine_distance8_feature_payment_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal7_finite_cross_payments_20261017'
OLD='fresh_selective_distance9_dimension108_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
MIXED='fresh_ideal_distance7_mixed_replacement_20261017'

def labelled(n,edges):
    g=[[0]*n for _ in range(n)]
    for i,j,m in edges:g[i][j]=g[j][i]=m
    return g

def embeddings(pattern,target):
    """All injective induced labelled graph maps, without graph-library code.

    At each step every unused target vertex consistent with all preceding
    edges AND nonedges is tried. Degree tests are necessary only.
    """
    n,m=len(pattern),len(target);pd=[sum(bool(x)for x in r)for r in pattern];td=[sum(bool(x)for x in r)for r in target]
    result=[];mp={};used=set()
    def visit():
        if len(mp)==n:result.append([mp[i]for i in range(n)]);return
        i=max((i for i in range(n)if i not in mp),key=lambda i:(sum(bool(pattern[i][j])for j in mp),pd[i],-i))
        for j in range(m):
            if j in used or td[j]<pd[i]:continue
            if any(pattern[i][k]!=target[j][v]for k,v in mp.items()):continue
            mp[i]=j;used.add(j);visit();used.remove(j);del mp[i]
    visit()
    assert len({tuple(r)for r in result})==len(result)
    return result

def distance_data(g):
    n=len(g);adj=[{j for j in range(n)if g[i][j]}for i in range(n)];ds=[]
    for i in range(n):
        d={i:0};todo=deque([i])
        while todo:
            v=todo.popleft()
            for j in adj[v]:
                if j not in d:d[j]=d[v]+1;todo.append(j)
        assert len(d)==n;ds.append(d)
    features={i for i in range(n)if len(adj[i])>=3 or any(g[i][j]>=4 for j in adj[i])}
    return adj,ds,features

def basic_pairs(g):
    adj,ds,features=distance_data(g)
    near={i for i in range(len(g))if any(ds[i][j]<=7 for j in features)}
    return {p for p in combinations(range(len(g)),2)if ds[p[0]][p[1]]<=7 or(ds[p[0]][p[1]]in(8,9)and(set(p)&near))}

def dropped_d17(g):
    adj,ds,features=distance_data(g)
    return{p for p in combinations(range(17),2)if ds[p[0]][p[1]]>=13 or(len(adj[p[0]])==len(adj[p[1]])==3)}

def component(g,ids,start):
    out={start};todo=[start]
    while todo:
        i=todo.pop()
        for j in ids:
            if j not in out and g[i][j]:out.add(j);todo.append(j)
    return out

def coupon(native,ids,a,b):
    ca,cb=component(native,ids,a),component(native,ids,b)
    if ca==cb:return F()
    for small,large,x,y in((ca,cb,a,b),(cb,ca,b,a)):
        if(len(small),len(large))not in((1,16),(2,14)):continue
        if any(native[i][j]not in(0,-1)for i,j in combinations(small|large,2)):continue
        if len(small)==2 and not native[x][next(i for i in small if i!=x)]:continue
        adj={i:{j for j in large if j!=i and native[i][j]}for i in large}
        forks=[i for i in large if len(adj[i])==3]
        if len(forks)!=1 or any(len(adj[i])>3 for i in large):continue
        f=forks[0];short=[j for j in adj[f]if len(adj[j])==1]
        # A finite D14 or D16 has precisely these two short arms.
        if len(short)!=2:continue
        if y==f or y in short:return F(1,2)
    return F()

def case_payments(native,scaled,black,left,right):
    white=[i for i in range(len(native))if i not in black];rows=[]
    for a,b in product(left,right):
        c=R.inspect_corner(native,scaled,black,a,b);ids=white+[a,b]
        finite=all(z['kind']=='finite'for z in c['blocks'])
        cp=component(native,ids,a);dist=R.graph_dist([[native[i][j]for j in ids]for i in ids],[ids.index(a)],False)[ids.index(b)]
        base=F(int(finite and b in cp and(dist<=9 or len(cp)<=17)))
        extra=coupon(native,ids,a,b)if finite else F()
        rows.append(dict(marks=[a,b],finite=finite,component_orders=[len(component(native,ids,v))for v in(a,b)],
            distance=dist,base=str(base),cross_coupon=str(extra),weight=str(base+extra),blocks=c['blocks']))
    return rows

def order17():
    old=json.loads((ROOT/f'results/{OLD}.json').read_text());diagrams={r['index']:r for r in json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']}
    patterns=[]
    for index,r in enumerate(old['isolation']['cases']):
        assert r['component_order_bound']<=17
        if r['component_order_bound']<17:continue
        di,u=r['left'];dj,v=r['right'];g,mp=R.unite(R.bond_graph(diagrams[di]),R.bond_graph(diagrams[dj]),tuple(map(tuple,r['match'])))
        ids=r['marked_component'];h=[[int(2-g[i][j])if i!=j and g[i][j]else 0 for j in ids]for i in ids]
        patterns.append((index,r,g,ids,h,[ids.index(x)for x in r['preferred_corner']]))
    assert len(patterns)==80
    targets={kind:labelled(17,A.P.graph(kind,17)[0])for kind in'ABCD'}
    dropped=dropped_d17(targets['D']);assert len(dropped)==9
    records={};lost_cases=set();saved_patterns=[]
    for index,r,g,ids,h,marks in patterns:
        saved_patterns.append(dict(index=index,nodes=ids,labels=h,marks=marks,method=r['method'],order_bound=r['component_order_bound']))
    for kind,target in targets.items():
        basic=basic_pairs(target);extra=set();rows=[];allmaps=0
        for index,r,g,ids,h,marks in patterns:
            maps=embeddings(h,target);allmaps+=len(maps)
            pairs={tuple(sorted(m[i]for i in marks))for m in maps}
            missing=pairs-basic;extra|=missing
            if pairs:rows.append(dict(case=index,embeddings=maps,marked_images=[list(p)for p in sorted(pairs)],extra_images=[list(p)for p in sorted(missing)]))
            if kind=='D'and missing&dropped:lost_cases.add(index)
        paid=(basic|extra)-(dropped if kind=='D'else set())
        records[kind]=dict(basic_pairs=[list(p)for p in sorted(basic)],extra_pairs=[list(p)for p in sorted(extra)],
            dropped_pairs=[list(p)for p in sorted(dropped)]if kind=='D'else[],paid_pairs=[list(p)for p in sorted(paid)],
            cost=len(paid),embedding_count=allmaps,case_embeddings=rows)
    assert{kind:r['cost']for kind,r in records.items()}==dict(A=119,B=111,C=108,D=127)
    assert lost_cases=={383,393,400,647,662,675,813,828,841}
    lost=[]
    for index in sorted(lost_cases):
        r=old['isolation']['cases'][index];di,u=r['left'];dj,v=r['right']
        assert di==dj==25 and r['match']==[[9,9]]and r['method']=='radical'and r['possible_contacts']==[3,13]
        native,gram,sc,mp=R.combined(diagrams,di,dj,((9,9),));right=[mp[i]for i in v];black=u+right
        radical=R.nullspace(gram);assert len(radical)==1 and radical[0][3]*radical[0][13]<0
        assert all(radical[0][i]for i in r['possible_contacts'])
        payments=case_payments(native,sc,black,u,right)
        # Every external supporting white must be orthogonal to the core:
        # otherwise radical balance contacts both3 and13, a length2 payment.
        assert sum(F(p['weight'])for p in payments)>=1
        lost.append(dict(case=index,black=black,possible_contacts=r['possible_contacts'],
            radical=list(map(str,radical[0])),short_path_if_nonorthogonal=[3,'external_white',13],
            corners=payments,finite_total=str(sum(F(p['weight'])for p in payments))))
    return dict(patterns=saved_patterns,targets=records,dropped_case_indices=sorted(lost_cases),replacement_cases=lost,
        embedding_completeness='Explicit injective backtracking tries every target vertex compatible with all previously mapped edges and nonedges; no graph-isomorphism library is imported.')

def affine_costs():
    rows=[]
    for m in range(18,129):
        for kind in'ABCD':
            g=labelled(m,A.P.graph(kind,m)[0]);n=len(basic_pairs(g))
            formula=7*m if kind=='A'else 7*m-{'B':23,'C':28,'D':18}[kind]+min(m-{'B':10,'C':8,'D':12}[kind],{'B':17,'C':18,'D':16}[kind])+min(m-{'B':11,'C':9,'D':13}[kind],{'B':17,'C':18,'D':16}[kind])
            assert n==formula and 2*n<=15*m
            rows.append(dict(kind=kind,order=m,cost=n))
    small=[]
    for kind,m,es in A.small_catalogue():
        if m>16:continue
        n=m*(m-1)//2;assert 2*n<=15*m
        small.append(dict(kind=kind,order=m,cost=n))
    branches=[]
    # In each interval the affine cost is a*m+b. The final slope7 has
    # increasing slack, so its left endpoint controls the infinite tail.
    for kind,s1,s2,b1,b2,b3 in[('C',26,27,-45,-19,8),('B',27,28,-44,-17,11),('D',28,29,-43,-15,14)]:
        for lo,hi,a,b in[(18,s1,9,b1),(s1,s2,8,b2),(s2,None,7,b3)]:
            at=hi if a>F(15,2)else lo;slack=F(15,2)*at-(a*at+b);assert slack>=0
            branches.append(dict(kind=kind,lower_order=lo,upper_order=hi,slope=a,intercept=b,minimum_at=at,minimum_slack=str(slack)))
    return dict(classical_checks=rows,whole_small=small,all_rank_linear_branches=branches,
        all_rank_formulas=dict(A='7m',B='7m-23+min(m-10,17)+min(m-11,17)',C='7m-28+min(m-8,18)+min(m-9,18)',D='7m-18+min(m-12,16)+min(m-13,16)'))

def finite_budgets():
    rows=[]
    for n in range(92,109):
        worst=F(-10**9);witness=None;count=0
        for a in range(n+1):
            for b in range((n-a)//2+1):
                for c in range((n-a-2*b)//14+1):
                    for d in range((n-a-2*b-14*c)//16+1):
                        diff=-9*a-17*b-35*c-24*d+F(3,2)*a*d+3*b*c
                        assert diff<=0;count+=1
                        if diff>worst:worst=diff;witness=[a,b,c,d]
        rows.append(dict(dimension=n,component_count_tuples=count,max_excess_over_9n=str(worst),maximizer=witness))
    return dict(rows=rows,formula='cost<=9n-9a-17b-35c-24d+(3/2)ad+3bc',
        count_meaning='a=#A1,b=#A2,c=#D14,d=#D16; a+2b+14c+16d<=n',
        theorem='For n<=108, d<=6 makes a*((3/2)d-9)<=0. If c<=5 then b*(3c-17)-35c<=0. If c=6 or7, b<=(108-14c)/2 gives the same strict inequality. Thus all displayed dimensions have total finite cost<=9n, independent of the remaining irreducible component ranks.')

def run():
    special=order17();affine=affine_costs();finite=finite_budgets()
    stems=[OLD,OLD+'_review',CENSUS,MIXED,'fresh_affine_distance8_feature_payment_20261017','fresh_cusp_factor_payment_framework_20261017']
    paths=[Path(__file__)]
    for stem in stems:
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            p=ROOT/f'{folder}/{stem}.{ext}'
            if p.is_file():paths.append(p)
    return dict(status='PASS_IDEAL7_SELECTIVE17_AND_FINITE_CROSS_COMPONENT_PAYMENTS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        order17=special,affine_cost=affine,finite_cost=finite,
        scope='Local face-payment transfer and all-rank ideal cost<=min(15m/2,m(m-1)/2), with finite cost<=9n for92<=n<=108. Requires the written all-rank graph-count proof and inherited complete108 classification. No ambient dimension bound is inferred by this source alone.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
