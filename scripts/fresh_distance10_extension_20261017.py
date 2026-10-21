#!/usr/bin/env python3
"""Exact local supporting-white extension bounds for two distance10 cores.

Fraction arithmetic only. No complete distance10 classification or dimension
bound is asserted here; arbitrary indirect whites are handled analytically by
the positive residual-space lemma in the accompanying report.
"""
import argparse
from collections import Counter
from itertools import product
from math import prod
from pathlib import Path
import json
import fresh_distance11_extension_review_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance10_extension_20261017'
F=M.F


def build(di,mp,ds):
    n=max(ds[di]['order'],max(mp)+1);g=M.matrix(n)
    for d,m in [(ds[di],list(range(ds[di]['order']))),(ds[25],mp)]:
        for i,j,label in d['edges']:
            assert label==3;g[m[i]][m[j]]=g[m[j]][m[i]]=-F(1)
    return g


def schur(g,centers):
    ids=[i for i in range(len(g))if i not in centers];a=M.extract(g,ids);piv=M.psd(a,True)
    b=M.extract(g,centers,ids);c=M.multiply(M.multiply(b,M.inverse(a)),M.transpose(b))
    s=[[g[i][j]-c[k][l]for l,j in enumerate(centers)]for k,i in enumerate(centers)]
    det=s[0][0]*s[1][1]-s[0][1]**2;assert det<0
    return s,prod(piv)*det


def corner(g,u,v,p,q):
    omitted=sorted(set(u+v)-{p,q});ids=[i for i in range(len(g))if i not in omitted]
    h=M.extract(g,ids);M.psd(h);allowed=omitted.copy();parts=[];marked=None
    for part in M.blocks(h):
        a=M.extract(h,part);nodes=[ids[i]for i in part];ks=M.nullspace(a)
        if ks:
            assert len(ks)==1;k=ks[0]
            if all(x<0 for x in k):k=[-x for x in k]
            assert all(x>0 for x in k)
            parts.append(dict(nodes=nodes,kind='affine',positive_kernel=list(map(str,k))))
        else:
            inv=M.inverse(a);assert all(x>=0 for row in inv for x in row)
            vals=[inv[i][i]for i in range(len(a))]
            allowed.extend(nodes[i]for i,x in enumerate(vals)if x<=2)
            parts.append(dict(nodes=nodes,kind='finite',inverse_diagonals=list(map(str,vals))))
        if p in nodes:assert q in nodes;marked=nodes
    assert marked is not None
    return dict(preferred=[p,q],blocks=parts,possible_contacts=sorted(allowed),
                marked_component=marked,marked_distance=M.graph_dist(h,[ids.index(p)],False)[ids.index(q)])


def attach(g,node,tail):
    n=len(g);h=M.matrix(n+1+tail)
    for i in range(n):
        for j in range(n):h[i][j]=g[i][j]
    h[n][node]=h[node][n]=-F(1)
    for j in range(tail):h[n+j][n+j+1]=h[n+j+1][n+j]=-F(1)
    return h


def family(di,mp,ds):
    g=build(di,mp,ds);n=len(g);j=M.inverse(g);s,det=schur(g,[0,mp[0]])
    kind='four_node_L25_overlap'if di==25 else'mixed_L22_L25_edge'
    if di==25:
        left=[1,2,3];right=[mp[i]for i in left]
        block=[[F(2),2,1],[2,F(10,3),F(5,3)],[1,F(5,3),F(4,3)]]
        assert M.extract(j,left)==M.extract(j,right)==block and det==-3
        assert s==[[-F(1,18),-F(1,9)],[-F(1,9),-F(1,18)]]
        assert j[1][mp[1]]==-3 and j[1][3]==1 and j[1][mp[3]]==-2 and j[3][mp[3]]==-F(4,3)
        assert F(4,3)-F(2,3)>0 and -F(4,3)+F(2,3)<0
        roots=[]
        for node in [1,mp[1]]:
            root=[-row[node]for row in j];assert M.square(g,root)==2
            assert M.multiply(g,[[x]for x in root])==[[-F(i==node)]for i in range(n)]
            roots.append(dict(contact_node=node,coordinates=list(map(str,root))))
    else:
        left=[1,3,4];right=[mp[i]for i in [1,2,3]]
        block=[[F(7,2),3,F(3,2)],[3,4,2],[F(3,2),2,F(3,2)]]
        assert M.extract(j,left)==block and min(j[i][i]for i in right)==4 and det==-2
        assert all(x>0 for row in M.extract(j,right)for x in row)
        assert s==[[-F(5,36),-F(1,9)],[-F(1,9),-F(1,18)]]
        assert F(3,2)-F(1,2)>0;roots=[]
    assert all(x>0 for row in block for x in row)
    rows=[];variant_counts=Counter()
    for u,v0 in product([[1,2],[1,3],[2,3]],repeat=2):
        v=[mp[i]for i in v0];du,dv=M.graph_dist(g,u),M.graph_dist(g,v)
        assert min(du[i]for i in v)>10
        white=[i for i in range(n)if i not in u+v]
        for p in u+v:M.psd(M.extract(g,white+[p]),True)
        for p,q in product(u,v):M.psd(M.extract(g,white+[p,q]))
        pref=lambda pair:1 if pair==[1,3]else 2
        p,q=pref(u),mp[pref(v0)];c=corner(g,u,v,p,q)
        assert set(c['possible_contacts'])<=set(left+right)
        dl=max(du[i]for i in c['possible_contacts']if i in left);dr=max(dv[i]for i in c['possible_contacts']if i in right)
        assert dl+dr+2<=5
        assert len(c['marked_component'])<=(14 if di==25 else 15)
        cases=[]
        if di==25:
            choices=[('fixed_A',node,0)for node in [1,mp[1]]if node in c['possible_contacts']]
            choices += [('residual_B',node,t)for node in[3,mp[3]]for t in range(3)]
        else:choices=[('residual_B',4,t)for t in range(2)]
        for typ,node,tail in choices:
            h=attach(g,node,tail);w=[i for i in range(len(h))if i not in u+v]
            for mark in u+v:M.psd(M.extract(h,w+[mark]),True)
            for a,b in product(u,v):M.psd(M.extract(h,w+[a,b]))
            ext=corner(h,u,v,p,q);assert len(ext['marked_component'])<=17
            variant_counts[typ,tail,len(ext['marked_component'])]+=1
            cases.append(dict(type=typ,core_contact=node,additional_tail_roots=tail,
                              marked_component=ext['marked_component']))
        common=set(range(ds[di]['order']))&set(mp)
        extra=[]
        for a,b in product(sorted(set(range(ds[di]['order']))-common),sorted(set(mp)-common)):
            h=[row.copy()for row in g];h[a][b]=h[b][a]=-F(1)
            if min(M.graph_dist(h,u)[t]for t in v)>10:extra.append([a,b])
        assert not extra
        rows.append(dict(left_black=u,right_black=v,preferred_corner=c,
                         both_end_shortcut_bound=dl+dr+2,allowed_extension_variants=cases))
    return dict(kind=kind,gram=M.encode(g),right_copy_map=mp,signature=[n-1,1],
                determinant=str(det),schur=M.encode(s),end_nodes=[left,right],
                left_inverse_block=M.encode(block),right_inverse_block=M.encode(M.extract(j,right)),
                cross_inverse_block=M.encode(M.extract(j,left,right)),fixed_projection_roots=roots,
                contexts=rows,marked_component_order_bound=17,
                variant_counts=[dict(type=t,tail=l,component_order=o,count=c)for(t,l,o),c in sorted(variant_counts.items())])


def run():
    ds={d['index']:d for d in json.loads((ROOT/f'results/{M.CENSUS}.json').read_text())['diagrams']}
    a1=[[F(2)]];a2=[[F(2),-F(1)],[-F(1),F(2)]]
    a3=[[F(2),-F(1),0],[-F(1),F(2),-F(1)],[0,-F(1),F(2)]]
    assert M.inverse(a1)[0][0]==F(1,2)and M.inverse(a2)[0][0]==F(2,3)
    assert M.inverse(a3)[0][0]==F(3,4)and M.inverse(a3)[1][1]==1
    # Exact residual Grams for maximal permitted indirect extensions.
    residual_three=[[F(2,3),-F(1),0],[-F(1),F(2),-F(1)],[0,-F(1),F(2)]]
    residual_two=[[F(1,2),-F(1)],[-F(1),F(2)]]
    M.psd(residual_three);M.psd(residual_two)
    assert len(M.nullspace(residual_three))==len(M.nullspace(residual_two))==1
    families=[family(25,list(range(10,16))+[9,8,7,6],ds),family(22,list(range(9,17))+[8,7],ds)]
    assert sum(len(c['allowed_extension_variants'])for f in families for c in f['contexts'])==78
    paths=[f'scripts/{TAG}.py',f'scripts/{M.TAG}.py',f'results/{M.TAG}.json',
           f'scripts/{M.CENSUS}.py',f'results/{M.CENSUS}.json',
           'results/fresh_cusp_factor_payment_framework_20261017.md',
           'results/fresh_selective_distance11_dimension126_20261017.md']
    return dict(status='PASS_EXACT_DISTANCE10_TWO_LOCAL_EXTENSION_COMPONENT_BOUNDS',
                sources={p:M.sha(ROOT/p)for p in paths},families=families,
                residual_budget=dict(A1_inverse='1/2',A2_endpoint_inverse='2/3',
                    minimum_connected_three_node_inverse='3/4',
                    maximum_residual_gram_three=M.encode(residual_three),maximum_residual_gram_two=M.encode(residual_two)),
                contexts=18,extension_variants=78,marked_component_order_bound=17,
                scope='Two specified local overlap families under no marked white path of length<=10. Arbitrarily many indirect supporting whites are bounded by the proved residual-space argument, not by an enumeration cutoff. No complete distance10 census, dimension bound, or finite-volume realization is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'18contexts/78extensions/component<=17',flush=True)
