#!/usr/bin/env python3
"""Exact radical isolation of the nine surviving L22/8--L25/9 contexts."""
import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import fresh_distance11_extension_review_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance11_extension_review_20261017_mixed'
F=M.F


def run():
    ds={d['index']:d for d in json.loads((ROOT/f'results/{M.CENSUS}.json').read_text())['diagrams']}
    g=M.matrix(18);mp=list(range(9,18))+[8]
    for d,m in [(ds[22],list(range(9))),(ds[25],mp)]:
        for i,j,label in d['edges']:
            assert label==3;g[m[i]][m[j]]=g[m[j]][m[i]]=-F(1)
    ids=[i for i in range(18)if i not in(0,9)]
    a=M.extract(g,ids);M.psd(a,True)
    b=M.extract(g,[0,9],ids);cost=M.multiply(M.multiply(b,M.inverse(a)),M.transpose(b))
    s=[[g[i][j]-cost[k][l]for l,j in enumerate([0,9])]for k,i in enumerate([0,9])]
    assert s==[[-F(3,20),-F(1,10)],[-F(1,10),-F(1,15)]]
    assert s[0][0]*s[1][1]-s[0][1]**2==0
    z=list(map(F,[-4,-2,-3,-2,-1,-3,-2,-1,0,6,3,4,2,5,4,3,2,1]))
    assert M.nullspace(g)==[z]
    assert M.multiply(g,[[x]for x in z])==[[0]for _ in z]
    cases=[]
    for u,v0 in product([[1,2],[1,3],[2,3]],repeat=2):
        v=[mp[i]for i in v0];du,dv=M.graph_dist(g,u),M.graph_dist(g,v)
        assert min(du[i]for i in v)>11
        corners=[]
        for p,q in product(u,v):
            omit=set(u+v)-{p,q};cs=[i for i in range(18)if i not in omit]
            M.psd(M.extract(g,cs));corners.append([p,q])
        white=[i for i in range(18)if i not in u+v]
        for p in u+v:M.psd(M.extract(g,white+[p]),True)
        pref=lambda pair:1 if pair==[1,3]else 2
        p,q=pref(u),mp[pref(v0)]
        omit=sorted(set(u+v)-{p,q});cs=[i for i in range(18)if i not in omit]
        h=M.extract(g,cs);allowed=omit.copy();parts=[];marked=None
        for part in M.blocks(h):
            b=M.extract(h,part);nodes=[cs[i]for i in part];ks=M.nullspace(b)
            if ks:
                assert len(ks)==1;k=ks[0]
                if all(x<0 for x in k):k=[-x for x in k]
                assert all(x>0 for x in k)
                parts.append(dict(nodes=nodes,kind='affine',positive_kernel=list(map(str,k))))
            else:
                j=M.inverse(b);assert all(x>=0 for row in j for x in row)
                diag=[j[i][i]for i in range(len(part))]
                allowed.extend(nodes[i]for i,d in enumerate(diag)if d<=2)
                parts.append(dict(nodes=nodes,kind='finite',inverse_diagonal=list(map(str,diag))))
            if p in nodes:
                assert q in nodes;marked=nodes;kind='affine'if ks else'finite'
        assert marked is not None and len(marked)<=16
        allowed=sorted(allowed);assert all(z[i]for i in allowed)
        neg=[i for i in allowed if z[i]<0];pos=[i for i in allowed if z[i]>0]
        assert max(du[i]for i in neg)<=2 and max(dv[i]for i in pos)<=1
        extras=[]
        for i,j in product(range(8),range(9,18)):
            a=[row.copy()for row in g];a[i][j]=a[j][i]=-F(1)
            if min(M.graph_dist(a,u)[t]for t in v)>11:extras.append([i,j])
        assert not extras
        cases.append(dict(left_black=u,right_black=v,all_adjacent_corners=corners,
                          preferred=[p,q],blocks=parts,allowed_contacts=allowed,
                          negative_radical_contacts=neg,positive_radical_contacts=pos,
                          marked_component=marked,marked_kind=kind,
                          marked_distance=M.graph_dist(h,[cs.index(p)],False)[cs.index(q)],
                          contact_forced_shortcut_bound=5))
    counts=Counter((len(c['marked_component']),c['marked_kind'])for c in cases)
    paths=[f'scripts/{TAG}.py',f'scripts/{M.TAG}.py',f'results/{M.TAG}.json',
           f'scripts/{M.CENSUS}.py',f'results/{M.CENSUS}.json',
           'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_EXACT_DISTANCE11_MIXED_CORE_RADICAL_ISOLATION',
                sources={p:M.sha(ROOT/p)for p in paths},core=dict(gram=M.encode(g),
                right_copy_map=mp,schur=M.encode(s),signature=[16,1,1],radical=list(map(str,z))),
                contexts=cases,marked_component_order_bound=16,
                component_counts=[dict(order=n,kind=k,count=c)for(n,k),c in sorted(counts.items())],
                scope='Nine specified surviving L22/8--L25/9 marked contexts, up to reversing the two copies. Exact local isolation only; no complete overlap classification, new dimension bound or finite-volume realization is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'9contexts/component<=16/shortcut<=5',flush=True)
