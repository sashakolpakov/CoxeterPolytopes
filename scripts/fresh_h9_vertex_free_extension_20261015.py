#!/usr/bin/env python3
"""Global exact27-wall geometry after two vertex-free edge caps cut oldescape."""
import argparse,json
from collections import Counter
from itertools import combinations
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,C,ROOT=B.A,B.D,B.P,B.C,B.ROOT
TAG='fresh_h9_vertex_free_extension_20261015'
BASE='results/fresh_h9_quartet_polar_second_20261014.json'
CAPS='results/fresh_h9_vertex_free_caps_20261015.json'

def run():
    d,p=B.read(BASE),B.read(CAPS);A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
    N.extend(D.dv(r['normal'])for r in p['caps']);assert len(N)==27
    norms=[A.pair(n,n,Q)for n in N];G=[[A.pair(n,m,Q)for m in N]for n in N]
    counts=Counter();relations=[]
    for i,j in combinations(range(27),2):
        kind,m,c=R.relation(G[i][j],norms[i],norms[j]);assert kind!='forbidden'
        counts[kind+(''+str(m)if m else'')]+=1
        if j>=25:relations.append(dict(pair=[i,j],kind=kind,label=m,invariant=c))
    assert len(relations)==51 and not G[25][26]
    assert A.pair(H,H,Q).sign()<0 and all(A.pair(H,n,Q).sign()<0 for n in N)
    old_escape_pairs=[A.pair(E,n,Q)for n in N];assert all(z.sign()<0 for z in old_escape_pairs[:25])and all(z.sign()>0 for z in old_escape_pairs[25:])
    facets=[]
    for i,n in enumerate(N):
        t=A.pair(H,n,Q)/norms[i];v=[x-t*y for x,y in zip(H,n)];vals=[A.pair(v,m,Q)for m in N]
        assert A.pair(v,v,Q).sign()<0 and A.pair(v,H,Q).sign()<0
        assert not vals[i]and all(z.sign()<0 for j,z in enumerate(vals)if j!=i)
        facets.append(dict(index=i,point=v,norm=A.pair(v,v,Q)))
    census=C.diagnose_chamber_cone(N,Q,H,9);names=('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')
    census_counts={k:len(getattr(census,k))for k in names};print('27-wall cone',census_counts,flush=True)
    funcs=[A.mv(Q,n)for n in N];vertices=[]
    for active in census.finite_vertices:
        v,vals,S=B.ray(active,funcs)
        assert len(active)==9 and A.pair(v,v,Q).sign()<0 and A.pair(v,H,Q).sign()<0
        vertices.append(dict(incident=list(active),point=v,norm=A.pair(v,v,Q),pairings=vals))
    old={tuple(v['incident'])for v in d['finite_vertices']};assert old<={tuple(v['incident'])for v in vertices}
    edges=[]
    for i,j in combinations(range(len(vertices)),2):
        S=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
        if len(S)==8:edges.append(dict(vertices=[i,j],support=S))
    closed={(i,tuple(e['support']))for e in edges for i in e['vertices']}
    germs=[]
    for i,v in enumerate(vertices):
        for release in v['incident']:
            S=tuple(j for j in v['incident']if j!=release)
            if(i,S)not in closed:germs.append(dict(vertex=i,support=list(S),released_wall=release))
    escape=None
    if census.hyperideal_rays:
        for active in census.hyperideal_rays:
            v,vals,S=B.ray(active,funcs)
            if A.pair(v,H,Q).sign()<0:break
        else:raise AssertionError('need general future mixing')
        epsilon=A.F(1)
        while True:
            X=[x+A.L(epsilon)*h for x,h in zip(v,H)]
            if A.pair(X,X,Q).sign()>0:break
            epsilon/=2
        xp=[A.pair(X,n,Q)for n in N]
        assert A.pair(X,H,Q).sign()<0 and all(z.sign()<0 for z in xp)
        escape=dict(point=X,norm=A.pair(X,X,Q),center_pairing=A.pair(X,H,Q),pairings=xp,
                    hyperideal_ray=v,active=list(active),center_epsilon=str(epsilon))
    paths=[BASE,CAPS,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','src/coxeter_polytopes/verify/chamber.py']
    return A.enc(dict(status='PASS_EXACT_VERTEX_FREE_CAP_EXTENSION',sources={q:R.sha(q)for q in paths},
                     field=d['field'],form=Q,normals=N,norms=norms,gram=G,center=H,pair_counts=dict(counts),new_pairs=relations,
                     facet_points=facets,old_escape_pairings=old_escape_pairs,
                     complete_census={k:[list(t)for t in getattr(census,k)]for k in names},counts=census_counts,
                     finite_vertices=vertices,finite_edges=edges,open_finite_germs=germs,
                     new_finite_vertices=len(vertices)-len(old),strict_escape=escape,
                     scope='All27actualfacets and all351 Coxeter pairs checked. Full cone census and finite-edge graph reconstructed. Two newwalls cut the oldescape; a new strictescape, if present, proves currentinfinitevolume. No novelty orfinitevolume claim from localedgeclosure alone.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],len(d['finite_edges']),len(d['open_finite_germs']),'escape',bool(d['strict_escape']))
