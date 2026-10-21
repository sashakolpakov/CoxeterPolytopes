#!/usr/bin/env python3
"""Optimistic exact free-weight LP:112 quadrilaterals+98 actual triangle stars.

Every connected marked-pair orbit receives an independent nonnegative weight.
This finite relaxation omits other triangle external-white extensions and all other
quadrilateral frontiers. It is a diagnostic, not a dimension proof.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import fresh_sqrt5_sealed_weights_20261015 as T
import fresh_sqrt5_quad_filters_20261015_review as Q
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_free_weights_20261015'

def component_record(H):
    ids=sorted(H);n=len(ids);key=T.keystr(T.tree_key(H));orbits=defaultdict(list)
    for a,b in combinations(ids,2):orbits[T.keystr(T.tree_key(H,(a,b)))].append([a,b])
    return key,dict(type=T.typename(H),rank=n,nodes=ids,labels=[H[a].get(b,2)for a,b in combinations(ids,2)],
                   pair_orbits=[dict(key=k,multiplicity=len(v),pairs=v)for k,v in sorted(orbits.items())])

def actual_triangle(n,g,white):
    row=dict(nodes=n,labels=[g.get(p,2)for p in combinations(range(n),2)])
    M=Q.matrix(row,Q.coloring(row));assert Q.inertia(M)==[n-1,1,0]
    assert Q.inertia([[Q.conj(x)for x in r]for r in M])==[n,0,0]
    S=Q.residual(M,white);assert len(S)==3 and Q.inertia(S)==[2,1,0]
    def dot(a,b):
        out=Q.Z
        for x,y in zip(a,b):out=Q.add(out,Q.mul(x,y))
        return out
    vertices=[]
    for a,b in combinations(range(3),2):
        v,w=S[a],S[b]
        x=[Q.sub(Q.mul(v[1],w[2]),Q.mul(v[2],w[1])),
           Q.sub(Q.mul(v[2],w[0]),Q.mul(v[0],w[2])),
           Q.sub(Q.mul(v[0],w[1]),Q.mul(v[1],w[0]))]
        products=[dot(r,x)for r in S]
        if next(Q.sign(z)for z in products if z!=Q.Z)>0:x=list(map(Q.neg,x));products=list(map(Q.neg,products))
        norm=dot(x,products);assert Q.sign(norm)<0 and all(Q.sign(z)<=0 for z in products)
        assert [i for i,z in enumerate(products)if z==Q.Z]==[a,b]
        vertices.append(x)
    center=[Q.Z]*3
    for x in vertices:center=[Q.add(a,b)for a,b in zip(center,x)]
    products=[dot(r,center)for r in S];assert all(Q.sign(z)<0 for z in products)and Q.sign(dot(center,products))<0
    assert all(Q.sign(dot(x,products))<0 for x in vertices)
    return dict(projected_gram=[[Q.enc(x)for x in r]for r in S],
                finite_corner_norms=[Q.enc(dot(x,[dot(r,x)for r in S]))for x in vertices])

def build():
    union=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    contacts=json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    triangles=json.loads((ROOT/'results/fresh_sqrt5_triangle_weights_20261015.json').read_text())
    selected=[r for r in contacts['records']if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and r['singular']and not r['accepted']]
    assert len(selected)==112 and len(triangles['markings'])==76
    comps={};faces=[]
    def addface(kind,ident,white,cornerpairs,g):
        coef=Counter();corners=[]
        for a,b in cornerpairs:
            H=T.graph(white+[a,b],g);C=T.component(H,a)
            if b not in C:
                corners.append(dict(marked=[a,b],connected=False));continue
            H={i:{j:m for j,m in H[i].items()if j in C}for i in C};ck,cr=component_record(H)
            if ck in comps:
                old=comps[ck];assert old['type']==cr['type']
                assert [(o['key'],o['multiplicity'])for o in old['pair_orbits']]==[(o['key'],o['multiplicity'])for o in cr['pair_orbits']]
            else:comps[ck]=cr
            pk=T.keystr(T.tree_key(H,(a,b)));coef[(ck,pk)]+=1
            corners.append(dict(marked=[a,b],connected=True,component=ck,pair_orbit=pk,type=cr['type']))
        faces.append(dict(kind=kind,index=ident,required=2 if kind=='triangle'else 1,corners=corners,
                          coefficients=[dict(component=c,orbit=p,count=v)for(c,p),v in sorted(coef.items())]))
    for rec in selected:
        row=union['labelled_union_supports'][rec['union_index']]
        addface('quadrilateral',rec['class_index'],row['white'],list(product((0,2),(1,3))),T.edges(row))
    for i,rec in enumerate(triangles['markings']):
        g=dict(zip(combinations(range(rec['order']),2),rec['labels']))
        addface('triangle',i,rec['white'],list(combinations(rec['black'],2)),g)
    known=json.loads((ROOT/'results/fresh_arithmetic_quadratic_bound_20261015.json').read_text());geometry=[]
    for record in known['records']:
        n=int(record['type'][1:]);old={p:2 for p in combinations(range(n+1),2)}
        for i,j in [(i,i+1)for i in range(n-2)]+[(2,n-1)]:old[i,j]=3
        for pi,profile in enumerate(record['profiles']):
            for a,b in profile['triangle_marked_pairs']:
                g=dict(old)
                for i,m in profile['contacts']:g[i,n]=m
                white=[i for i in range(n)if i not in(a,b)]
                ident=f"{record['type']}:profile{pi}:pair{a},{b}"
                geometry.append(dict(index=ident,**actual_triangle(n+1,g,white)))
                addface('triangle',ident,white,list(combinations([a,b,n],2)),g)
    assert len(geometry)==22
    variables=sorted((c,o['key'])for c,r in comps.items()for o in r['pair_orbits']);index={p:i for i,p in enumerate(variables)}
    A=[];b=[];meta=[]
    for f in faces:
        line=[0]*(len(variables)+1)
        for c in f['coefficients']:line[index[c['component'],c['orbit']]]=-c['count']
        A.append(line);b.append(-f['required']);meta.append(dict(kind='face',face_kind=f['kind'],index=f['index']))
    for ck,r in sorted(comps.items()):
        line=[0]*(len(variables)+1);line[-1]=1
        for o in r['pair_orbits']:line[index[ck,o['key']]]=o['multiplicity']
        A.append(line);b.append(2*r['rank']);meta.append(dict(kind='component',component=ck,type=r['type']))
    return comps,variables,faces,A,b,meta,geometry

def run():
    comps,V,faces,A,b,meta,geometry=build();x,y,slack=T.solve(A,b);rows=[]
    for ck,r in sorted(comps.items()):
        vals=[];total=F(0)
        for o in r['pair_orbits']:
            weight=x[V.index((ck,o['key']))];total+=o['multiplicity']*weight;vals.append(dict(**o,weight=str(weight)))
        rows.append(dict(key=ck,**{k:v for k,v in r.items()if k!='pair_orbits'},pair_orbits=vals,total=str(total),deficit=str(2*r['rank']-total)))
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_sealed_weights_20261015.py',
           'results/fresh_sqrt5_quad_exhaustion_20261015.json','results/fresh_sqrt5_quad_contacts_20261015.json',
           'scripts/fresh_sqrt5_triangle_weights_20261015.py','results/fresh_sqrt5_triangle_weights_20261015.json',
           'scripts/fresh_arithmetic_quadratic_bound_20261015.py','results/fresh_arithmetic_quadratic_bound_20261015.json',
           'results/fresh_arithmetic_quadratic_bound_20261015_review.json',
           'scripts/fresh_sqrt5_quad_filters_20261015_review.py','results/fresh_sqrt5_quad_filters_20261015_review.json']
    return dict(status='PASS_EXACT_OPTIMISTIC_FREE_WEIGHT_RELAXATION',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                variables=len(V),faces=len(faces),components=len(comps),common_component_deficit=str(x[-1]),
                component_rows=rows,face_rows=faces,known_triangle_geometry=geometry,exact_primal=list(map(str,x)),minimum_slack=str(min(slack)),
                exact_dual=[dict(constraint=meta[i],multiplier=str(z))for i,z in enumerate(y)if z],
                scope='Optimistic necessary LP only, component-local and zero on disconnected pairs:112 sealed no-short quad cores,76 minimal triangle markings and22 known actual E7/E8 triangles. Omits all short-path-paid quadrilateral unions (whose old q2 payment is invalid after lowering q2), further finite triangle extensions, nonsingular/dotted quad frontiers and absent component types. No weight theorem or dimension bound.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'faces',out['faces'],'variables',out['variables'],'components',out['components'],'max deficit',out['common_component_deficit'])
    print([(r['type'],r['total'])for r in out['component_rows']])
