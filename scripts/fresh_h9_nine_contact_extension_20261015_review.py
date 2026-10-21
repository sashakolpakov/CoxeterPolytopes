#!/usr/bin/env python3
"""Independent complete cone replay of the 30-wall replacement-escape extension."""
import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import fresh_h9_quartet_polar_extension_20261014_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_nine_contact_extension_20261015_review'
PRIMARY='results/fresh_h9_nine_contact_extension_20261015.json'
BASE='results/fresh_h9_fourth_escape_extension_20261015.json'
CAPS='results/fresh_h9_nine_contact_cap_20261015.json'
K,Z,O,C=R.K,R.Z,R.O,R.C

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):return list(map(str,x.to_list()))

def run():
    paths=[PRIMARY,BASE,CAPS,'results/fresh_h9_nine_contact_cap_20261015_review.json',
           'results/fresh_h9_fourth_escape_extension_20261015_review.json']
    docs=[json.loads((ROOT/p).read_text())for p in paths]
    for doc in docs:
        for p,h in doc.get('sources',{}).items():assert sha(p)==h,p
    d,b,caps,_,_=docs
    Q,N,H=R.dm(d['form']),R.dm(d['normals']),R.dv(d['center'])
    assert len(N)==30 and Q==R.dm(b['form'])and H==R.dv(b['center'])
    assert N==R.dm(b['normals'])+[R.dv(row['normal'])for row in caps['caps']]
    pair=lambda v,w:C.pair(v,w,Q)
    cov=[C.mv(Q,n)for n in N]
    norms=[pair(n,n)for n in N]
    assert norms==R.dv(d['norms'])and all(R.sg(q)>0 for q in norms)
    G=[[C.dot(n,row)for row in cov]for n in N]
    assert G==R.dm(d['gram'])
    count=Counter();new={tuple(row['pair']):row for row in d['new_pairs']}
    for i,j in combinations(range(30),2):
        typ,label,u=R.relation(G[i][j],norms[i],norms[j]);assert typ!='forbidden'
        count[typ+(str(label)if label else'')]+=1
        if j>=29:
            row=new.pop((i,j));assert(typ,label)==(row['kind'],row['label'])
            assert u==(None if row['invariant']is None else R.dec(row['invariant']))
    assert not new and dict(count)==d['pair_counts']
    assert R.sg(pair(H,H))<0 and all(R.sg(C.dot(H,row))<0 for row in cov)
    assert len(d['facet_points'])==30
    for i,row in enumerate(d['facet_points']):
        P=[a-C.dot(H,cov[i])/norms[i]*n for a,n in zip(H,N[i])]
        assert P==R.dv(row['point'])and pair(P,P)==R.dec(row['norm'])
        assert R.sg(pair(P,P))<0 and R.sg(pair(P,H))<0
        assert all(R.sg(C.dot(P,cc))==(0 if j==i else-1)for j,cc in enumerate(cov))
    cn,stages,seed=R.cone(N,Q,H)
    counts={name:len(rows)for name,rows in cn.items()}
    assert counts==d['counts']
    assert {key:sorted(map(list,rows))for key,rows in cn.items()}==d['complete_census']
    saved={tuple(row['incident']):row for row in d['finite_vertices']}
    assert set(saved)==set(cn['finite_vertices'])and len(saved)==counts['finite_vertices']
    for ids,v in cn['finite_vertices'].items():
        row=saved[ids];sv=R.dv(row['point'])
        assert len(ids)==9 and R.canon(v)==R.canon(sv)
        assert R.pd([[G[i][j]for j in ids]for i in ids])
        assert pair(sv,sv)==R.dec(row['norm'])and R.sg(pair(sv,sv))<0 and R.sg(pair(sv,H))<0
        vals=[C.dot(sv,cc)for cc in cov]
        assert vals==R.dv(row['pairings'])and all(R.sg(g)==(0 if i in ids else-1)for i,g in enumerate(vals))
    old={tuple(row['incident'])for row in b['finite_vertices']}
    assert old<set(saved)and len(saved)-len(old)==d['new_finite_vertices']
    vertices=d['finite_vertices'];edges=[]
    for i,j in combinations(range(len(vertices)),2):
        support=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
        if len(support)==8:
            assert R.pd([[G[a][b]for b in support]for a in support])
            edges.append(dict(vertices=[i,j],support=support))
    assert edges==d['finite_edges']
    closed={(i,tuple(row['support']))for row in edges for i in row['vertices']}
    germs=[dict(vertex=i,support=[j for j in row['incident']if j!=release],released_wall=release)
           for i,row in enumerate(vertices)for release in row['incident']
           if(i,tuple(j for j in row['incident']if j!=release))not in closed]
    assert germs==d['open_finite_germs']
    assert 9*len(saved)==2*len(edges)+len(germs)
    oldE=R.dv(b['strict_escape']['point']);oldvals=[C.dot(oldE,cc)for cc in cov]
    assert oldvals==R.dv(d['old_escape_pairings'])
    assert all(R.sg(g)<0 for g in oldvals[:29])and all(R.sg(g)>0 for g in oldvals[29:])
    escape=d['strict_escape'];E=R.dv(escape['point']);V=R.dv(escape['hyperideal_ray'])
    eps=K.convert(Fraction(escape['center_epsilon']));active=tuple(escape['active'])
    assert active in cn['hyperideal_rays']and R.canon(V)==R.canon(cn['hyperideal_rays'][active])
    assert E==[v+eps*h for v,h in zip(V,H)]and R.sg(eps)>0
    vals=[C.dot(E,cc)for cc in cov]
    assert vals==R.dv(escape['pairings'])and all(R.sg(g)<0 for g in vals)
    assert pair(E,E)==R.dec(escape['norm'])and R.sg(pair(E,E))>0
    assert pair(E,H)==R.dec(escape['center_pairing'])and R.sg(pair(E,H))<0
    paths += [f'scripts/{TAG}.py','scripts/fresh_h9_nine_contact_extension_20261015.py',
              'scripts/fresh_h9_nine_contact_cap_20261015_review.py',
              'scripts/fresh_h9_quartet_polar_extension_20261014_review.py',
              'scripts/fresh_h9_h8_all_high_20261014_review.py',
              'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_FULL_30_WALL_CONE',sources={p:sha(p)for p in paths},
        arithmetic='Separate SymPy ANP field and192-bit rational outward sign bounds; independent clipping by complete-face ray adjacency. No producer arithmetic or chamber-engine imports.',
        cone_seed=seed,clipping_stages=stages,counts=counts,pair_counts=dict(count),pairs=435,
        actual_facets=30,facet_sign_checks=900,finite_vertices=len(saved),vertex_sign_checks=30*len(saved),
        new_finite_vertices=len(saved)-len(old),finite_edges=len(edges),open_finite_germs=len(germs),
        old_escape_cut=True,strict_escape_norm=enc(pair(E,E)),strict_escape_center_pairing=enc(pair(E,H)),
        strict_escape_pairings=[enc(g)for g in vals],escape_sign_checks=30,
        scope='Complete30-wall extreme-ray census, finite vertices and finite-edge/open-germ incidence. The new wall cuts the previous escape, but the independently verified new strict spacelike escape proves this chamber still has infinite volume. No finite-volume conclusion or general-cap exclusion.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('PASS saved independent30 cone replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'])
