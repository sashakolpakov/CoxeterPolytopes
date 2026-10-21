#!/usr/bin/env python3
"""Exact native cyclic-product checks on all strongest-deletion cover terms.

No optimizer and no new physical exclusion. Every term is reconstructed in
its complete finite vertex frame. Pair squares and fundamental cycle
squareclasses certify all cyclic products, not physical realizability.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_restricted_290_unowned_pair_frontier_20261017 as A

N=A.N;R=N.R;M=N.M;K=N.K;V=N.V;W=N.W;ROOT=N.ROOT
TAG='fresh_restricted_290_cyclic_cover_20261017'


def normal(x):return json.loads(json.dumps(x,sort_keys=True))


def chosen_cover(row):
    q=row['all_unowned_deletion'];d=list(map(F,row['demands']))
    owner={i:j for j,g in enumerate(row['corner_orbits'])for i in g['literal_stars']}
    if 'certificate'in q:
        terms=[(q['literal_representatives'][i],F(x))for i,x in q['certificate']['mixture']]
        method='attaining_complete_cone_cover'
    else:
        terms=[(cl,F(x))for cl,x in q['split_cover']['terms']]
        method='verified_explicit_cover'
    projected=Counter()
    for cl,x in terms:
        assert x>0;cl=tuple(i for i in cl if d[owner[i]]>0);assert cl
        projected[cl]+=x
    marginal=[F()for _ in d]
    for cl,x in projected.items():
        for i in cl:marginal[owner[i]]+=x
    assert all(x>=y for x,y in zip(marginal,d))
    mass=sum(projected.values(),F());assert mass==F(q['optimal_mass_upper_bound'])
    return sorted(projected.items()),method,mass,marginal


def cycle_certificate(n,edges):
    """Nonzero edge squares are integral K elements; tree potentials suffice."""
    adj=[[]for _ in range(n)]
    for i,j,s in edges:
        assert i<j and s.sign()>0 and M.integral(s)
        adj[i].append((j,s));adj[j].append((i,s))
    potentials={};parents={};roots=[]
    for root in range(n):
        if root in potentials:continue
        roots.append(root);potentials[root]=K(1);parents[root]=None;todo=[root]
        while todo:
            i=todo.pop()
            for j,s in adj[i]:
                if j not in potentials:
                    potentials[j]=potentials[i]*s;parents[j]=i;todo.append(j)
    checks=[]
    for i,j,s in edges:
        square=potentials[i]*potentials[j]*s
        root=M.ksqrt(square);assert root is not None and root.sign()>0
        assert root*root==square and M.integral(root)
        checks.append([i,j,root.out()])
    return dict(spanning_roots=roots,spanning_parents=[parents[i]for i in range(n)],
        tree_squareclass_potentials=[potentials[i].out()for i in range(n)],
        edge_membership_roots=checks,cycle_rank=len(edges)-n+len(roots))


def forest_check(v,record,p,changed,rows,patches):
    old,g,owner,unowned,demand,groups,forest,applied=A.current(v,p,changed,rows,patches)
    assert forest==record['forest']and normal(groups)==record['corner_orbits']
    pairs=sorted(pair for pp in unowned.values()for pair in pp)
    assert A.O.pairhash(pairs)==record['all_unowned_deletion']['pair_domain_sha256']
    g.remove_edges_from(pairs)
    nodes=A.T.current(v,p,changed)[5]
    terms,method,mass,marginal=chosen_cover(record)
    used=sorted({i for cl,x in terms for i in cl});native={i:N.native(nodes[i])for i in used}
    walls=sorted({q for cv in native.values()for q in cv});wallindex={q:i for i,q in enumerate(walls)}
    target=R.P.forest_graph(forest,dict(R.C.models()));J,active,components=R.projection_data(target)
    assert sorted(target)==list(range(13))
    baseedges=[(i,j,V.squared(V.edge(m)))for i,ns in sorted(target.items())for j,m in sorted(ns.items())if i<j]
    heights=[];contactedges=[]
    for wi,q in enumerate(walls):
        assert all(i in active for i,e in q)
        field=N.pair_field((q,));projection=R.projection(q,q,J,field)
        assert not any(projection[1:]);h=projection[0]-2;assert h.sign()>0;heights.append(h)
        contacts=[]
        for i,e in q:
            sq=V.squared(e);assert sq.sign()>0 and M.integral(sq);contacts.append((i,wi,sq))
        contactedges.append(contacts)
    selfchecks=[]
    for i,cv in sorted(native.items()):
        ids=[wallindex[q]for q in cv];hz=[heights[j]for j in ids];st=nodes[i]
        assert [h.out()for h in hz]==st['squared_timelike_heights']
        edge=tuple(st['native_omitted_product']);field=N.pair_field(cv,(edge,))
        projection=R.projection(cv[0],cv[1],J,field);nativeproduct=field.negative_contact(edge)
        assert field.signs_difference(tuple(x-y for x,y in zip(projection,nativeproduct)),hz[0]*hz[1])==0
        selfchecks.append(dict(star=i,walls=ids,native_omitted_product=list(edge)))
    paircache={};certs=[];counts=Counter()
    for cl,amount in terms:
        assert all(g.has_edge(i,j)for i,j in combinations(cl,2))
        ids=sorted({wallindex[q]for i in cl for q in native[i]});local={w:13+j for j,w in enumerate(ids)}
        edges=list(baseedges)
        for w in ids:edges.extend((i,local[w],s)for i,wi,s in contactedges[w])
        for a,b in combinations(ids,2):
            if (a,b)not in paircache:
                field=N.pair_field((walls[a],walls[b]));projection=R.projection(walls[a],walls[b],J,field)
                h=heights[a]*heights[b];sg=field.signs_difference(projection,h)
                sg1=field.signs_difference(R.qaddk(projection,1),h)
                assert sg<=0 and not(sg<0 and sg1>0)
                code,_=field.arithmetic_wall_gate(tuple((z.a,z.b)for z in projection),h.a,h.b)
                assert code is None,(v,cl,a,b,code)
                rt=field.squareclass_root(h.a,h.b)
                if rt is None:
                    assert not any(projection);square=h
                else:
                    i,r=rt;q=list(projection);q[i]=q[i]-r;sq=field.qm(q,q)
                    assert not any(sq[1:]);square=sq[0]
                assert square.sign()>=0 and M.integral(square)
                paircache[a,b]=dict(walls=[a,b],pair_field_radicals=[field.tau.out(),field.sigma.out()],
                    projection=R.qout(projection),height_product=h.out(),product_sign=sg,
                    squared_product=square.out())
            square=K(*map(F,paircache[a,b]['squared_product']))
            if square:edges.append((local[a],local[b],square))
        edges.sort(key=lambda z:(z[0],z[1]));assert len({(i,j)for i,j,s in edges})==len(edges)
        cert=cycle_certificate(13+len(ids),edges)
        certs.append(dict(literal_stars=list(cl),amount=str(amount),walls=ids,
            nonzero_squared_edges=[[i,j,s.out()]for i,j,s in edges],**cert))
        counts.update(cover_terms=1,external_wall_occurrences=len(ids),cycle_rank=cert['cycle_rank'],
            edge_membership_checks=len(edges),singleton_star_terms=int(len(cl)==1),multi_star_terms=int(len(cl)>1))
    counts.update(distinct_external_walls=len(walls),distinct_cross_wall_pairs=len(paircache),native_star_checks=len(selfchecks))
    oldmu=F(record['all_unowned_deletion']['baseline_mu']);newmu=max(oldmu,mass)
    assert newmu==oldmu+F(record['all_unowned_deletion']['required_mu_increase_upper_bound'])
    return normal(dict(vertex=v,forest=forest,removed_unowned_pairs=len(pairs),pair_domain_sha256=A.O.pairhash(pairs),
        demands=list(map(str,demand)),covered_marginals=list(map(str,marginal)),cover_method=method,
        cover_mass=str(mass),preserved_multiplier=str(newmu),
        native_wall_contacts=walls,squared_timelike_heights=[x.out()for x in heights],
        native_star_checks=selfchecks,cross_wall_checks=[paircache[k]for k in sorted(paircache)],
        cover_terms=certs,counts=dict(counts)))


def run(verify=False):
    start=time.monotonic();p,changed,rows,patches,base=A.prepare();diag=W.read(A.TAG)
    sources=A.bindings()
    for folder,ext in(('scripts','py'),('results','json'),('results','md')):
        path=f'{folder}/{A.TAG}.{ext}';sources[path]=W.digest(path)
    for mod in(R,):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    path=f'scripts/{TAG}.py';sources[path]=W.digest(path)
    saved=W.read(TAG)if verify else None
    if saved:assert saved['sources']==sources
    counts=Counter();records=[];mus={i:F(x)for i,x in diag['baseline_restricted_repair']['vertex_multipliers']}
    for i,row in enumerate(diag['forest_records']):
        r=forest_check(row['vertex'],row,p,changed,rows,patches)
        if saved:assert r==saved['forest_records'][i]
        records.append(r);counts.update(r['counts']);counts['forests']+=1
        mus[r['vertex']]=F(r['preserved_multiplier'])
        if (i+1)%40==0:print('FORESTS',i+1,'/571',dict(counts),flush=True)
    assert counts['forests']==571
    normalizer=sum(mus.values(),F());assert normalizer==F(diag['optimistic_all_unowned']['required_normalizer_upper'])
    num=F(diag['baseline_restricted_repair']['objective_numerator']);ratio=num/normalizer
    assert ratio==F(diag['optimistic_all_unowned']['hypothetical_ratio_lower'])>12
    out=normal(dict(status='PASS_RESTRICTED_290_CYCLIC_PRODUCT_COVER_SURVIVAL',sources=sources,
        baseline_source=A.TAG,forest_records=records,counts=dict(counts),
        hypothetical_restricted_witness=dict(objective_numerator=str(num),normalizer=str(normalizer),bound=str(ratio),
            vertex_multipliers=[[i,str(x)]for i,x in sorted(mus.items())if x]),
        actual_new_incompatibilities=0,actual_exported_columns=[],actual_witness_unchanged=True,
        scope='Every positive-demand projected term of the571 strongest genuinely-unowned-deletion covers. Complete common finite frames, actual native covectors/positive heights, all cross-wall2-cycles and all fundamental cycle squareclasses rechecked. Certifies survival of this hypothetical above12 relaxation witness under cyclic-product necessities in these local frames. Does not prove physical realizability, universal pair-gate sufficiency, global facet-cycle consistency, target feasibility, optimality or a dimension bound. No optimizer.'))
    path=ROOT/f'results/{TAG}.json'
    if saved:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'bound',str(ratio),float(ratio),'seconds',round(time.monotonic()-start,2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();run(args.verify)
