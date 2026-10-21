#!/usr/bin/env python3
"""Exact finite weight LP for the112 sealed Qsqrt5 principal stars.

Optimization proposes rational coefficients; every primal and dual inequality
is verified with Fraction. This does not cover nonsingular or dotted frontiers.
"""
import argparse,json
from collections import Counter,defaultdict,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_sealed_weights_20261015'

def edges(row):return dict(zip(combinations(range(row['nodes']),2),row['labels']))
def pair(i,j):return tuple(sorted((i,j)))
def graph(nodes,g):return {i:{j:g[pair(i,j)]for j in nodes if j!=i and g[pair(i,j)]!=2}for i in nodes}
def component(H,start):
    seen={start};queue=[start]
    for i in queue:
        for j in H[i]:
            if j not in seen:seen.add(j);queue.append(j)
    return sorted(seen)
def distances(H,start):
    out={start:0};queue=[start]
    for i in queue:
        for j in H[i]:
            if j not in out:out[j]=out[i]+1;queue.append(j)
    return out

def tree_key(H,marks=()):
    assert sum(map(len,H.values()))==2*(len(H)-1)
    marks=set(marks)
    def rooted(i,parent):return (int(i in marks),tuple(sorted((m,rooted(j,i))for j,m in H[i].items()if j!=parent)))
    return min(rooted(i,None)for i in H)
def keystr(k):return json.dumps(k,separators=(',',':'))

def typename(H):
    n=len(H);labels=[m for i in H for j,m in H[i].items()if i<j]
    if n==1:return 'A1'
    if n==2:return f'I2({labels[0]})'
    branches=[i for i in H if len(H[i])>=3]
    if not branches:
        if 5 in labels:return f'H{n}'
        if 4 in labels:
            e=next((i,j)for i in H for j,m in H[i].items()if m==4)
            return f'B{n}'if any(len(H[i])==1 for i in e)else'F4'
        return f'A{n}'
    assert len(branches)==1 and len(H[branches[0]])==3
    center=branches[0];arms=[]
    for start in H[center]:
        prev,i=center,start;length=1
        while len(H[i])==2:
            j=next(j for j in H[i]if j!=prev);prev,i=i,j;length+=1
        arms.append(length)
    return f'D{n}'if sorted(arms)[:2]==[1,1]else f'E{n}'

def baseline(H,i,j):
    d=distances(H,i).get(j)
    if d is None:return 0
    value=int(d<=2)
    if typename(H)=='H4' and (d==3 or H.get(i,{}).get(j)==5):value+=1
    return value

def component_record(H):
    ids=sorted(H);n=len(ids);g={(i,j):H[i].get(j,2)for i,j in combinations(ids,2)}
    canonical=keystr(tree_key(H));orbits=defaultdict(list)
    for i,j in combinations(ids,2):
        if distances(H,i)[j]>=3 and not baseline(H,i,j):orbits[keystr(tree_key(H,(i,j)))].append([i,j])
    return canonical,dict(rank=n,type=typename(H),nodes=ids,
        labels=[g[p]for p in combinations(ids,2)],
        base_total=sum(baseline(H,i,j)for i,j in combinations(ids,2)),
        pair_orbits=[dict(key=k,multiplicity=len(v),pairs=v)for k,v in sorted(orbits.items())])

def build():
    union=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    contacts=json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    selected=[r for r in contacts['records']if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and r['singular']and not r['accepted']]
    assert len(selected)==112
    components={};face_rows=[]
    for rec in selected:
        row=union['labelled_union_supports'][rec['union_index']];g=edges(row);corners=[];coef=Counter();base=0
        for a,b in product((0,2),(1,3)):
            H=graph(row['white']+[a,b],g);C=component(H,a)
            if b not in C:
                corners.append(dict(marked=[a,b],connected=False));continue
            H={i:{j:m for j,m in H[i].items()if j in C}for i in C}
            ck,cr=component_record(H)
            if ck in components:
                old=components[ck];assert old['rank']==cr['rank']and old['base_total']==cr['base_total']
                assert [(z['key'],z['multiplicity'])for z in old['pair_orbits']]==[(z['key'],z['multiplicity'])for z in cr['pair_orbits']]
            else:components[ck]=cr
            bvalue=baseline(H,a,b);base+=bvalue
            pk=keystr(tree_key(H,(a,b)))
            if not bvalue:coef[(ck,pk)]+=1
            corners.append(dict(marked=[a,b],connected=True,component=ck,pair_orbit=pk,base=bvalue,type=cr['type']))
        assert base==0
        face_rows.append(dict(class_index=rec['class_index'],union_index=rec['union_index'],corners=corners,
                              coefficients=[dict(component=c,orbit=p,count=v)for(c,p),v in sorted(coef.items())],base=base))
    variables=sorted((c,o['key'])for c,r in components.items()for o in r['pair_orbits'])
    index={p:i for i,p in enumerate(variables)};A=[];b=[];constraints=[]
    for row in face_rows:
        line=[0]*(len(variables)+1)
        for c in row['coefficients']:line[index[c['component'],c['orbit']]]=-c['count']
        A.append(line);b.append(row['base']-1);constraints.append(dict(kind='face',class_index=row['class_index']))
    for ck,rec in sorted(components.items()):
        line=[0]*(len(variables)+1);line[-1]=1
        for o in rec['pair_orbits']:line[index[ck,o['key']]]=o['multiplicity']
        A.append(line);b.append(2*rec['rank']-rec['base_total']);constraints.append(dict(kind='component',component=ck,type=rec['type']))
    return components,variables,face_rows,A,b,constraints

def solve(A,b):
    n=len(A[0]);c=[0]*(n-1)+[-1]
    result=linprog(c,A_ub=np.array(A,dtype=float),b_ub=np.array(b,dtype=float),bounds=[(0,None)]*(n-1)+[(None,None)],method='highs')
    assert result.success,result.message
    x=[F(float(v)).limit_denominator(1000000)for v in result.x]
    y=[F(float(v)).limit_denominator(1000000)for v in result.ineqlin.marginals]
    assert all(z>=0 for z in x[:-1])and all(z<=0 for z in y)
    slack=[F(rhs)-sum(F(a)*z for a,z in zip(row,x))for row,rhs in zip(A,b)]
    assert all(z>=0 for z in slack),('primal rationalization failed',min(slack))
    dual=[sum(y[i]*A[i][j]for i in range(len(A)))for j in range(n)]
    assert all(dual[j]<=c[j]for j in range(n-1))and dual[-1]==-1
    primal=-x[-1];dualvalue=sum(F(rhs)*z for rhs,z in zip(b,y));assert primal==dualvalue
    return x,y,slack

def component_minimum(ck,rec,variables,A,b,faces):
    n=len(variables);c=[0]*n
    for o in rec['pair_orbits']:c[variables.index((ck,o['key']))]=o['multiplicity']
    AA=[row[:-1]for row in A[:len(faces)]];bb=b[:len(faces)]
    result=linprog(c,A_ub=np.array(AA,dtype=float),b_ub=np.array(bb,dtype=float),bounds=(0,None),method='highs')
    assert result.success
    x=[F(float(v)).limit_denominator(1000000)for v in result.x]
    y=[F(float(v)).limit_denominator(1000000)for v in result.ineqlin.marginals]
    assert all(v>=0 for v in x)and all(v<=0 for v in y)
    assert all(sum(F(a)*z for a,z in zip(row,x))<=rhs for row,rhs in zip(AA,bb))
    assert all(sum(y[i]*AA[i][j]for i in range(len(AA)))<=c[j]for j in range(n))
    value=sum(F(a)*z for a,z in zip(c,x));assert value==sum(F(rhs)*z for rhs,z in zip(bb,y))
    return dict(type=rec['type'],base_total=rec['base_total'],minimum_extra=str(value),
                minimum_total=str(value+rec['base_total']),primal=list(map(str,x)),
                dual=[dict(class_index=faces[i]['class_index'],multiplier=str(z))for i,z in enumerate(y)if z])

def run():
    components,variables,faces,A,b,constraints=build();x,y,slack=solve(A,b)
    rows=[]
    for ck,rec in sorted(components.items()):
        vals=[];extra=F(0)
        for o in rec['pair_orbits']:
            weight=x[variables.index((ck,o['key']))];extra+=weight*o['multiplicity']
            vals.append(dict(**o,weight=str(weight)))
        rows.append(dict(key=ck,**{k:v for k,v in rec.items()if k!='pair_orbits'},pair_orbits=vals,
                         extra_total=str(extra),total=str(extra+rec['base_total']),deficit=str(2*rec['rank']-rec['base_total']-extra)))
    dual=[dict(constraint=constraints[i],multiplier=str(z))for i,z in enumerate(y)if z]
    minima=[component_minimum(ck,rec,variables,A,b,faces)for ck,rec in sorted(components.items())if rec['type']in('D7','E8')]
    assert {r['type']:r['minimum_total']for r in minima}=={'D7':'21','E8':'23'}
    files=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_quad_exhaustion_20261015.py','results/fresh_sqrt5_quad_exhaustion_20261015.json',
           'scripts/fresh_sqrt5_quad_contacts_20261015.py','results/fresh_sqrt5_quad_contacts_20261015.json',
           'results/fresh_sqrt5_quad_filters_20261015_review.json']
    return dict(status='PASS_EXACT_SEALED_QSQRT5_WEIGHT_LP',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},
                sealed_faces=len(faces),variables=len(variables),component_types=len(components),
                baseline='distance<=2 plus marked5-edge and endpoint bonuses on entireH4 only',
                common_component_deficit=str(x[-1]),component_rows=rows,face_rows=faces,
                exact_primal=list(map(str,x)),exact_dual=dual,minimum_primal_slack=str(min(slack)),
                component_minima=minima,rank29_forest_barrier=dict(forest=['D7','D7','D7','E8'],rank=29,forced_minimum_weight=86,target_upper=84),
                proof='All integer-coefficient inequalities and strong-duality equality verified with Fraction; optimizer is proposal-only.',
                scope='Component-local weights only; disconnected marked pairs have weight zero; fixed distance-two baseline plus nonnegative within-component long-pair corrections. Only112 sealed singular finite cores. Contacts census sealing is inherited, not independently reproved. Nonsingular extensions and all dotted cases remain open. No global dimension bound.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'faces',out['sealed_faces'],'variables',out['variables'],'components',out['component_types'],'max deficit',out['common_component_deficit'])
    print([(r['type'],r['total'],r['deficit'])for r in out['component_rows']])
