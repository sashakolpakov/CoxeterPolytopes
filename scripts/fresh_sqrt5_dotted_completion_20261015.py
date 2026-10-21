#!/usr/bin/env python3
"""Exact missing-entry rank completions of the unpaid dotted seed frontier.

This enumerates zero-further-contact completions of the24 seed cores and
their369 inherited first contacts. It does not close arbitrary attachments.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as R
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import networkx as nx
import fresh_sqrt5_singular_extension_frontier_20261015 as F

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_dotted_completion_20261015'
Q=F.Q;K=Q.K
def dot(a,b):return sum((x*y for x,y in zip(a,b)),K())
def enc(M):return [[x.enc()for x in r]for r in M]
def inverse(M):
    n=len(M);A=[r[:]+[K(i==j)for j in range(n)]for i,r in enumerate(M)]
    for j in range(n):
        p=next(i for i in range(j,n)if A[i][j]);A[j],A[p]=A[p],A[j]
        d=A[j][j];A[j]=[x/d for x in A[j]]
        for i in range(n):
            if i!=j:
                c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[j])]
    B=[r[n:]for r in A]
    assert all(dot(M[i],[B[k][j]for k in range(n)])==int(i==j)for i in range(n)for j in range(n))
    return B
def rational_root(x):
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return R(a,b)if a*a==x.numerator and b*b==x.denominator else None
def square_root(x):
    if x.sign()<0:return None
    if not x:return K()
    if not x.b:
        a=rational_root(x.a)
        if a is not None:return K(a)
        b=rational_root(x.a/5)
        if b is not None:return K(0,b)
        return None
    norm=rational_root(x.a*x.a-5*x.b*x.b)
    if norm is None:return None
    for s in(-1,1):
        a=rational_root((x.a+s*norm)/2)
        if a:
            y=K(a,x.b/(2*a))
            if y*y==x:return y if y.sign()>0 else-y
    return None
def integral(x):return (2*x.b).denominator==1 and (x.a-x.b).denominator==1
def cycles(G):
    n=len(G);H=nx.Graph();H.add_nodes_from(range(n));H.add_edges_from((i,j)for i,j in combinations(range(n),2)if G[i][j])
    rows=[];cs=[list(e)for e in H.edges]+list(nx.simple_cycles(H))
    for c in cs:
        z=K(1)
        for i,j in zip(c,c[1:]+c[:1]):z=z*2*G[i][j]/G[i][i]
        rows.append(dict(cycle=c,value=z.enc(),integral=integral(z)))
    return rows

def evaluate(row):
    n=row['nodes'];g=Q.edges(row)
    assert [p for p,m in g.items()if m==0]==[(0,2)]
    M=Q.tree_scaled(n,{p:2 if m==0 else m for p,m in g.items()})
    if M is None:return dict(status='REJECT_FINITE_CYCLE_FIELD')
    cliques=[]
    for omitted in(0,2):
        ids=[i for i in range(n)if i!=omitted];B=Q.restrict(M,ids)
        p=Q.inertia(B);c=Q.inertia([[x.conj()for x in r]for r in B])
        cliques.append(dict(omitted=omitted,physical=p,conjugate=c))
        if p[1]>1:return dict(status='REJECT_PHYSICAL_PRINCIPAL',cliques=cliques)
        if c[1]:return dict(status='REJECT_CONJUGATE_PRINCIPAL',cliques=cliques)
    ids=[i for i in range(n)if i not in(0,2)];B=Q.restrict(M,ids);p=Q.inertia(B)
    if p[2]:return dict(status='SINGULAR_BASE_REQUIRES_RELATION_ANALYSIS',cliques=cliques,base_inertia=p)
    assert p==[n-3,1,0]
    Bi=inverse(B);u=[M[0][i]for i in ids];v=[M[2][i]for i in ids]
    a=M[0][0]-dot(u,[dot(r,u)for r in Bi]);b=M[2][2]-dot(v,[dot(r,v)for r in Bi]);c=dot(u,[dot(r,v)for r in Bi])
    assert a.sign()>=0 and b.sign()>=0
    # Schur residual has off-diagonal -x-c. Corank one requires (x+c)^2=ab.
    H=nx.Graph();H.add_nodes_from(range(n));H.add_edges_from(p for p,m in g.items()if m not in(0,2))
    linked=nx.has_path(H,0,2);out=dict(cliques=cliques,base_inertia=p,schur=dict(a=a.enc(),b=b.enc(),c=c.enc()),finite_path=linked,candidates=[])
    if linked:
        root=square_root(a*b)
        if root is None:return dict(**out,status='NO_FIELD_VALUED_ZERO_CONTACT_COMPLETION')
        choices=[]
        for sign in(-1,1):
            x=-c+sign*root
            if x.sign()<=0 or (x*x-M[0][0]*M[2][2]).sign()<=0:continue
            G=[r[:]for r in M];G[0][2]=G[2][0]=-x
            choices.append((G,dict(dotted_entry=x.enc(),squared_dotted_ratio=(x*x/(M[0][0]*M[2][2])).enc())))
    else:
        assert not c;z=a*b;choices=[]
        if z.conj().sign()>0 and(z-M[0][0]*M[2][2]).sign()>0:
            # Multiply the whole0-component by sqrt(z). The new matrix is inK.
            component=set(nx.node_connected_component(H,0));G=[r[:]for r in M]
            for i in component:
                for j in component:G[i][j]=z*M[i][j]
            G[0][2]=G[2][0]=-z
            choices.append((G,dict(dotted_entry_square=z.enc(),scaled_component=sorted(component),squared_dotted_ratio=(z/(M[0][0]*M[2][2])).enc())))
    for G,record in choices:
        pi=Q.inertia(G);co=Q.inertia([[x.conj()for x in r]for r in G]);assert pi==[n-2,1,1]
        record.update(gram=enc(G),physical=pi,conjugate=co)
        if co[1]:record['status']='REJECT_CONJUGATE_COMPLETION'
        else:
            assert co==[n-1,0,1];cy=cycles(G);record['cycles']=cy
            record['status']='PASS_ARITHMETIC_NECESSARY_COMPLETION'if all(r['integral']for r in cy)else'REJECT_NONINTEGRAL_CYCLIC_PRODUCT'
        out['candidates'].append(record)
    out['status']='ZERO_CONTACT_COMPLETIONS_ENUMERATED'
    return out

def run():
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_singular_extension_frontier_20261015.py','scripts/fresh_sqrt5_quad_contacts_20261015_review.py','results/fresh_sqrt5_quad_exhaustion_20261015.json','results/fresh_sqrt5_quad_contacts_20261015.json']
    Q.manifest_check(paths)
    base=json.loads((ROOT/paths[-2]).read_text())['labelled_union_supports'];old=json.loads((ROOT/paths[-1]).read_text())
    records=[];counts=Counter();cc=Counter()
    for r in old['records']:
        if r['status']!='UNPAID_DOTTED_CONTACT_DOMAIN':continue
        core=base[r['union_index']]
        for ti,c in [(-1,None)]+list(enumerate(r['accepted'])):
            row=core if c is None else F.extension(core,c);z=evaluate(row);counts[z['status']]+=1;cc.update(v['status']for v in z.get('candidates',[]))
            records.append(dict(source_class_index=r['class_index'],union_index=r['union_index'],contact_index=ti,contacts=c,row=row,result=z))
    assert len(records)==393
    return dict(status='PASS_EXACT_DOTTED_ZERO_CONTACT_RANK_COMPLETIONS',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=dict(counts),candidate_counts=dict(cc),records=records,scope='All24 unpaid one-dotted minimal cores and their369 inherited first-white contacts. Principal restrictions apply to any completion. Zero-further-contact rank roots and necessary arithmetic tests are complete for nonsingular common bases. They do not classify arbitrary further attachments or establish a dimension bound.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();data=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],data['candidate_counts'])
