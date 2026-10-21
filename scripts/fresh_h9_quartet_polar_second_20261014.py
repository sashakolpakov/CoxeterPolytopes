#!/usr/bin/env python3
"""Add the next two compatible exterior polar caps and inspect the cone."""
import argparse,json
from collections import Counter
from itertools import combinations
import fresh_h9_quartet_polar_20261014 as P
A,D,C,ROOT=P.A,P.D,P.C,P.ROOT
TAG='fresh_h9_quartet_polar_second_20261014'
BASE='results/fresh_h9_quartet_polar_extension_20261014.json'
POLAR=BASE

def read(path):
    d=json.loads((ROOT/path).read_text())
    for p,h in d['sources'].items():assert P.sha(p)==h,p
    return d
def ray(active,functionals):
    S=next(s for s in combinations(active,9)if C._full_rank_kernel([functionals[i]for i in s])is not None)
    v=C._full_rank_kernel([functionals[i]for i in S]);vals=[A.dot(v,n)for n in functionals]
    if any(z.sign()>0 for z in vals):v=[-z for z in v];vals=[-z for z in vals]
    assert all(z.sign()<=0 for z in vals)
    assert tuple(i for i,z in enumerate(vals)if not z)==active
    return v,vals,list(S)

def run():
    d,p=read(BASE),read(POLAR);A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
    chosen=[r for r in p['next_polar_candidates']if r['coxeter_compatible']]
    assert [r['index']for r in chosen]==[82,85]
    assert all(D.dv([r['center_pairing']])[0].sign()<0 for r in chosen)
    N.extend(D.dv(r['normal'])for r in chosen)
    norms=[A.pair(n,n,Q)for n in N]
    G=[[A.pair(n,m,Q)for m in N]for n in N]
    pairs=[];pc=Counter()
    for i,j in combinations(range(len(N)),2):
        kind,label,value=P.relation(G[i][j],norms[i],norms[j])
        assert kind!='forbidden'
        pc[kind+(''+str(label)if label else'')]+=1
        if j>=23:pairs.append(dict(pair=[i,j],kind=kind,label=label,invariant=value))
    assert len(pairs)==47
    assert all(A.pair(H,n,Q).sign()<0 for n in N) and A.pair(H,H,Q).sign()<0
    assert all(A.pair(E,n,Q).sign()<0 for n in N) and A.pair(E,E,Q).sign()>0 and A.pair(E,H,Q).sign()<0
    facets=[]
    for i,n in enumerate(N):
        c=A.pair(H,n,Q)/norms[i];v=[h-c*z for h,z in zip(H,n)];vals=[A.pair(v,m,Q)for m in N]
        assert A.pair(v,v,Q).sign()<0 and A.pair(v,H,Q).sign()<0
        assert not vals[i] and all(z.sign()<0 for j,z in enumerate(vals)if j!=i)
        facets.append(dict(index=i,point=v,norm=A.pair(v,v,Q)))
    census=C.diagnose_chamber_cone(N,Q,H,9)
    names=('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')
    counts={k:len(getattr(census,k))for k in names};print('25-wall cone',counts,flush=True)
    funcs=[A.mv(Q,n)for n in N];vertices=[]
    for active in census.finite_vertices:
        v,vals,support=ray(active,funcs)
        assert len(active)==9 and A.pair(v,v,Q).sign()<0 and A.pair(v,H,Q).sign()<0
        vertices.append(dict(incident=list(active),point=v,norm=A.pair(v,v,Q),pairings=vals))
    old={tuple(v['incident'])for v in d['finite_vertices']}
    assert old<={tuple(v['incident'])for v in vertices}
    edges=[]
    for i,j in combinations(range(len(vertices)),2):
        S=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
        if len(S)==8:edges.append(dict(vertices=[i,j],support=S))
    closed={(v,tuple(e['support']))for e in edges for v in e['vertices']}
    germs=[]
    for i,v in enumerate(vertices):
        for release in v['incident']:
            support=tuple(j for j in v['incident']if j!=release)
            if(i,support)not in closed:germs.append(dict(vertex=i,support=list(support),released_wall=release))
    candidates=[];stats=Counter()
    for i,active in enumerate(census.hyperideal_rays):
        v,vals,support=ray(active,funcs);q=A.pair(v,v,Q);assert q.sign()>0
        bad=[];relations=[]
        for j,g in enumerate(vals):
            typ,m,c=P.relation(g,q,norms[j]);relations.append(dict(kind=typ,label=m,invariant=c))
            if typ=='forbidden':bad.append(j)
        stats['tested']+=1;stats['forbidden'if bad else'compatible']+=1
        if not bad:stats['cuts_escape'if A.pair(v,E,Q).sign()>0 else'misses_escape']+=1
        candidates.append(dict(index=i,active=list(active),independent_support=support,
            normal=v,norm=q,pairings=vals,relations=relations,forbidden_indices=bad,
            coxeter_compatible=not bad,center_pairing=A.pair(v,H,Q),escape_pairing=A.pair(v,E,Q)))
    paths=[BASE,POLAR,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_20261014.py']
    return A.enc(dict(status='PASS_EXACT_SECOND_POLAR_CAP_EXTENSION',sources={p:P.sha(p)for p in paths},
        field=d['field'],source_polar_indices=[r['index']for r in chosen],
        form=Q,normals=N,norms=norms,gram=G,new_pairs=pairs,pair_counts=dict(pc),center=H,
        facet_points=facets,complete_census={k:[list(x)for x in getattr(census,k)]for k in names},
        counts=counts,finite_vertices=vertices,finite_edges=edges,open_finite_germs=germs,
        new_finite_vertex_count=len(vertices)-len(old),
        strict_escape=dict(point=E,norm=A.pair(E,E,Q),pairings=[A.pair(E,n,Q)for n in N]),
        next_polar_counts=dict(stats),next_polar_candidates=candidates,
        scope='All25 facets are actual, and the complete cone is enumerated exactly. The inherited strict escape proves infinite volume. Every polar hyperplane of every exterior vertex of the enlarged cone is tested against the full finite-angle alphabet. Nonpolar caps and later geometric extensions remain outside this census.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['next_polar_counts'],len(d['finite_edges']),len(d['open_finite_germs']))
