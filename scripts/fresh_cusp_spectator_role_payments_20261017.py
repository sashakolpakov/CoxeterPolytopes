#!/usr/bin/env python3
"""Exact role-specific PSD extensions; no full nonsimple-cusp payment claim."""
import argparse
from collections import deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cusp_spectator_role_payments_20261017'
INPUTS=['results/fresh_simple_cusp_affine_tree_completion_20261017.md',
        'results/fresh_simple_cusp_affine_tree_terminal_audit_20261017.json',
        'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.json',
        'results/fresh_cusp_factor_payment_framework_20261017.json']

def arms(lengths):
    edges=[];leaves=[];z=1
    for length in lengths:
        old=0
        for _ in range(length):edges.append((old,z));old=z;z+=1
        leaves.append(old)
    return z,edges,leaves

def gram(n,edges):
    G=2*s.eye(n)
    for i,j in edges:G[i,j]=G[j,i]=-1
    return G

def distances(n,edges):
    adj=[[]for _ in range(n)]
    for i,j in edges:adj[i].append(j);adj[j].append(i)
    D={}
    for i in range(n):
        ds={i:0};q=deque([i])
        while q:
            x=q.popleft()
            for y in adj[x]:
                if y not in ds:ds[y]=ds[x]+1;q.append(y)
        assert len(ds)==n
        for j in range(i+1,n):D[i,j]=ds[j]
    return D

def affine_extension(G,contact):
    n=G.rows;v=s.Matrix(contact);A=G.row_join(-v).col_join((-v.T).row_join(s.Matrix([[2]])))
    kernel=A.nullspace();assert len(kernel)==1
    k=kernel[0]
    if k[0]<0:k=-k
    assert all(x>0 for x in k)
    for j in range(n+1):
        B=A.minor_submatrix(j,j)
        assert all(B[:i,:i].det()>0 for i in range(1,n+1))
    return dict(gram=[[str(x)for x in row]for row in A.tolist()],
                positive_kernel=list(map(str,k)),all_one_deletions_PD=True)

def run():
    n,e,leaves=arms([1,2,3]);G=gram(n,e);J=G.inv();d=leaves[1]
    assert n==7 and d==3 and J[d,d]==2
    ds=distances(n,e);assert ds[0,d]==2 and ds[leaves[0],d]==3
    v=[0]*7;v[d]=1;e7aff=affine_extension(G,v)
    n8,e8,leaves8=arms([1,1,5]);D8=gram(n8,e8);I=D8.inv();a,b=leaves8[:2]
    assert n8==8 and I.extract([a,b],[a,b])==s.Matrix([[2,s.Rational(3,2)],[s.Rational(3,2),2]])
    assert I[a,a]+2*I[a,b]+I[b,b]==7
    exts=[]
    for z in (a,b):
        v=[0]*8;v[z]=1;exts.append(affine_extension(D8,v))
    assert distances(n8,e8)[a,b]==2
    small=[]
    # For order<=7 every pair costs at most3m. For orders8..13 these
    # literal classical shapes exhaust the small gap before inherited all-rank formulas.
    for m in range(8,14):
        for kind in ('cycle','path','single_fork','double_fork'):
            if kind=='cycle':edges=[(i,(i+1)%m)for i in range(m)];ends=[]
            elif kind=='path':edges=[(i,i+1)for i in range(m-1)];ends=[list(range(6)),list(range(m-6,m))]
            elif kind=='single_fork':
                t=m-2;edges=[(i,i+1)for i in range(t-1)]+[(0,t),(0,t+1)]
                ends=[list(range(4))+[t,t+1],list(range(t-6,t))]
            else:
                t=m-4;edges=[(i,i+1)for i in range(t-1)]+[(0,t),(0,t+1),(t-1,t+2),(t-1,t+3)]
                ends=[list(range(4))+[t,t+1],list(range(t-4,t))+[t+2,t+3]]
            assert all(0<=i<m for E in ends for i in E)
            dist=distances(m,edges);pairs={p for p,d in dist.items()if d<=3}
            for E in ends:pairs.update(combinations(sorted(E),2))
            assert len(pairs)<=3*m
            small.append(dict(shape=kind,order=m,cost=len(pairs),budget=3*m))
    exceptional=[]
    for name,lens in [('affine_E6',[2,2,2]),('affine_E7',[1,3,3]),('affine_E8',[1,2,5])]:
        m,edges,_=arms(lens);dist=distances(m,edges)
        cost=m*(m-1)//2 if m<=7 else sum(d<=3 for d in dist.values())
        assert cost<=3*m
        exceptional.append(dict(type=name,order=m,edges=edges,cost=cost,budget=3*m))
    old=json.loads((ROOT/INPUTS[2]).read_text());assert len(old['affine_graph_checks'])==348
    assert all(row['weight']==3*row['order']for row in old['affine_graph_checks'])
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_EXACT_ROLE_SPECIFIC_PSD_SPECTATOR_PAYMENTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        E7=dict(edges=e,new_cap=d,allowed_old_mark_positions=[0,leaves[0]],inverse_diagonal=str(J[d,d]),extension=e7aff),
        D8=dict(edges=e8,marked_fork_leaves=[a,b],inverse_marked_block=[[str(x)for x in row]for row in I.extract([a,b],[a,b]).tolist()],single_extensions=exts),
        small_classical_cost_checks=small,exceptional_affine_cost_checks=exceptional,
        inherited_all_rank_cost_checks=348,
        ideal_rule='All pairs when affine order<=7; otherwise distance<=3 plus the classical fork-invariant end segments of order<=6.',
        scope='Complete PSD extension classification for the two specified role patterns, and a3m affine-factor cost rule in every rank. Face payment for general nonsimple-cusp quadrilaterals is not proved. No new ambient dimension bound.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['small_classical_cost_checks']),'small classical graphs;',len(d['exceptional_affine_cost_checks']),'exceptional graphs;',flush=True)
