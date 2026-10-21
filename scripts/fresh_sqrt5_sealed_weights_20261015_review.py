#!/usr/bin/env python3
"""Independent exact pair-orbit barrier and fifteen compact face geometries.

Uses the separately independent Fraction Qsqrt5 helper; no producer arithmetic,
tree-key routine, LP optimizer, or numerical geometry is imported.
"""
import argparse,importlib.util,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt5_sealed_weights_20261015';TAG=BASE+'_review'
HELPER='fresh_sqrt5_quad_contacts_20261015_review'
spec=importlib.util.spec_from_file_location('independent_phi',ROOT/f'scripts/{HELPER}.py');Q=importlib.util.module_from_spec(spec);spec.loader.exec_module(Q)
K=Q.K


def graph_isos(left,lg,right,rg):
    if len(left)!=len(right):return
    la=Q.graph(left,lg);ra=Q.graph(right,rg)
    signature=lambda adj,g,i:(len(adj[i]),tuple(sorted(g[Q.pair(i,j)]for j in adj[i])))
    candidates={i:[j for j in right if signature(la,lg,i)==signature(ra,rg,j)]for i in left}
    order=sorted(left,key=lambda i:(len(candidates[i]),-len(la[i]),i));mp={};used=set()
    def rec(k):
        if k==len(order):yield dict(mp);return
        i=order[k]
        for j in candidates[i]:
            if j in used:continue
            if any(lg.get(Q.pair(i,a),2)!=rg.get(Q.pair(j,b),2)for a,b in mp.items()):continue
            mp[i]=j;used.add(j);yield from rec(k+1);used.remove(j);del mp[i]
    yield from rec(0)


def base_weight(nodes,g,a,b):
    adj=Q.graph(nodes,g);d=Q.distance(adj,a,b)
    if d is None:return 0
    weight=int(d<=2)
    es=[g.get((i,j),2)for i,j in combinations(sorted(nodes),2)if g.get((i,j),2)!=2]
    if len(nodes)==4 and sorted(es)==[3,3,5]:
        assert sorted(len(adj[i])for i in nodes)==[1,1,2,2]
        if d==3 or g.get(Q.pair(a,b),2)==5:weight+=1
    return weight


def lp_checks(data,union,contacts):
    comps={};orbits={};variables=[];component_records=[]
    for rec in data['component_rows']:
        ids=rec['nodes'];g=dict(zip(combinations(ids,2),rec['labels']));key=rec['key'];assert Q.finite(ids,g)
        autos=list(graph_isos(ids,g,ids,g));assert autos
        base=sum(base_weight(ids,g,a,b)for a,b in combinations(ids,2));assert base==rec['base_total']
        remaining={p for p in combinations(ids,2)if not base_weight(ids,g,*p)};checked=[]
        for orbit in rec['pair_orbits']:
            pairs={tuple(p)for p in orbit['pairs']};seed=min(pairs)
            actual={Q.pair(mp[seed[0]],mp[seed[1]])for mp in autos}
            assert actual==pairs and len(pairs)==orbit['multiplicity']and pairs<=remaining
            remaining-=pairs;orbits[key,orbit['key']]=pairs;variables.append((key,orbit['key']))
            checked.append(dict(key=orbit['key'],pairs=[list(p)for p in sorted(pairs)],multiplicity=len(pairs)))
        assert not remaining
        comps[key]=(ids,g,rec);component_records.append(dict(key=key,type=rec['type'],base=base,automorphisms=len(autos),orbits=checked))
    variables.sort();assert len(variables)==94 and len(comps)==16
    idx={v:i for i,v in enumerate(variables)}
    selected={r['class_index']:r for r in contacts['records']if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and r['singular']and not r['accepted']}
    assert len(selected)==112 and set(selected)=={r['class_index']for r in data['face_rows']}
    constraints=[];face_records=[]
    for saved in data['face_rows']:
        row=union['labelled_union_supports'][saved['union_index']];g=Q.edges(row);coeff=Counter();base=0;connected=[]
        assert selected[saved['class_index']]['union_index']==saved['union_index']
        for a,b in product((0,2),(1,3)):
            ids=row['white']+[a,b];adj=Q.graph(ids,g);C=next(c for c in Q.components(adj)if a in c)
            if b not in C:continue
            matches=[]
            for key,(ref,rg,rec)in comps.items():
                mp=next(graph_isos(C,g,ref,rg),None)
                if mp is not None:matches.append((key,mp))
            assert len(matches)==1;key,mp=matches[0];target=Q.pair(mp[a],mp[b]);bw=base_weight(C,g,a,b);base+=bw
            keys=[o for(c,o),pairs in orbits.items()if c==key and target in pairs]
            if not bw:assert len(keys)==1;coeff[key,keys[0]]+=1
            connected.append(dict(marked=[a,b],component=key,mapped_pair=list(target),base=bw,orbit=keys[0]if keys else None))
        assert base==saved['base']==0
        assert coeff==Counter({(r['component'],r['orbit']):r['count']for r in saved['coefficients']})
        line=[0]*94
        for var,count in coeff.items():line[idx[var]]=count
        constraints.append(line);face_records.append(dict(class_index=saved['class_index'],connected_corners=connected,coefficients=[dict(component=c,orbit=o,count=v)for(c,o),v in sorted(coeff.items())]))
    primal=list(map(F,data['exact_primal']));assert len(primal)==95 and primal[-1]==-7 and all(v>=0 for v in primal[:-1])
    assert all(sum(a*x for a,x in zip(row,primal))>=1 for row in constraints)
    for key,(_,_,rec)in comps.items():
        extra=sum(F(o['multiplicity'])*primal[idx[key,o['key']]]for o in rec['pair_orbits'])
        assert rec['base_total']+extra<=2*rec['rank']-primal[-1]
        assert extra==F(rec['extra_total'])
    minima=[]
    for minimum in data['component_minima']:
        typ=minimum['type'];key=next(k for k,(_,_,r)in comps.items()if r['type']==typ);needed={};witnesses=[]
        for item in minimum['dual']:
            face=next(r for r in face_records if r['class_index']==item['class_index'])
            assert len(face['connected_corners'])==1 and len(face['coefficients'])==1
            corner=face['connected_corners'][0];coef=face['coefficients'][0]
            assert corner['component']==key and corner['base']==0 and coef['count']==1
            orbit=corner['orbit'];assert orbit not in needed
            multiplicity=len(orbits[key,orbit]);assert F(item['multiplier'])==-multiplicity
            needed[orbit]=multiplicity;witnesses.append(dict(class_index=item['class_index'],orbit=orbit,multiplicity=multiplicity,mapped_pair=corner['mapped_pair']))
        lower=sum(needed.values());assert lower==9 and len(needed)==(9 if typ=='E8'else 6)
        x=list(map(F,minimum['primal']));assert len(x)==94 and all(v>=0 for v in x)
        assert all(sum(a*b for a,b in zip(row,x))>=1 for row in constraints)
        value=sum(len(pairs)*x[idx[c,o]]for(c,o),pairs in orbits.items()if c==key)
        assert value==lower==F(minimum['minimum_extra'])
        total=comps[key][2]['base_total']+lower;assert total==F(minimum['minimum_total'])
        minima.append(dict(type=typ,base=comps[key][2]['base_total'],extra=lower,total=total,faces=witnesses))
    assert {m['type']:m['total']for m in minima}=={'D7':21,'E8':23}
    assert 3*21+23==86>84
    return dict(component_orbits=component_records,face_constraints=face_records,sparse_minima=minima,
                exact_common_optimum=-7,common_optimum_proof='Feasible complete94-variable primal atdelta=-7; sparseE8 lower23 forcesdelta<=16-23=-7.',
                rank29_forced_total=86)


def solve(A,b):
    n=len(A);M=[[K(v)for v in row]+[K(z)]for row,z in zip(A,b)]
    for col in range(n):
        p=next(i for i in range(col,n)if M[i][col]);M[col],M[p]=M[p],M[col];v=M[col][col];M[col]=[a/v for a in M[col]]
        for i in range(n):
            if i!=col:
                v=M[i][col];M[i]=[a-v*b for a,b in zip(M[i],M[col])]
    return [r[-1]for r in M]


def determinant(A):
    n=len(A);value=K()
    for p in permutations(range(n)):
        v=K(1)
        for i,j in enumerate(p):v=v*A[i][j]
        if sum(p[i]>p[j]for i,j in combinations(range(n),2))%2:v=-v
        value=value+v
    return value


def null_line(rows):
    # Rectangular RREF, distinct from the primary cross-product construction.
    A=[[K(z)for z in row]for row in rows];piv=[];rank=0
    for col in range(3):
        p=next((i for i in range(rank,2)if A[i][col]),None)
        if p is None:continue
        A[rank],A[p]=A[p],A[rank];z=A[rank][col];A[rank]=[a/z for a in A[rank]]
        for i in range(2):
            if i!=rank:
                z=A[i][col];A[i]=[a-z*b for a,b in zip(A[i],A[rank])]
        piv.append(col);rank+=1
    assert rank==2;free=next(i for i in range(3)if i not in piv);v=[K(i==free)for i in range(3)]
    for i,col in enumerate(piv):v[col]=-A[i][free]
    return v


def dot(v,w):return sum((a*b for a,b in zip(v,w)),K())


def geometry(row,saved):
    G=Q.tree_scaled(row['nodes'],Q.edges(row));assert G is not None
    phys=Q.inertia(G);conj=Q.inertia([[x.conj()for x in r]for r in G])
    assert phys==saved['full_physical_inertia']and conj==saved['full_conjugate_inertia']
    W=row['white'];H=Q.restrict(G,W);assert Q.inertia(H)==[len(W),0,0]
    projections={j:solve(H,[G[i][j]for i in W])for j in range(4)}
    S=[[G[i][j]-dot([G[i][w]for w in W],projections[j])for j in range(4)]for i in range(4)]
    assert Q.inertia(S)==[2,1,1]
    basis=next(list(t)for t in combinations(range(4),3)if determinant(Q.restrict(S,t)))
    B=Q.restrict(S,basis);cov=[[S[i][j]for j in basis]for i in range(4)];feasible=[];rejected=[]
    for pair in combinations(range(4),2):
        v=null_line([cov[i]for i in pair]);p=[dot(r,v)for r in cov]
        if all(x.sign()>=0 for x in p):v=[-x for x in v];p=[-x for x in p]
        norm=dot(v,[dot(r,v)for r in B]);rec=dict(pair=list(pair),vector=[x.enc()for x in v],pairings=[x.enc()for x in p],norm=norm.enc())
        if all(x.sign()<=0 for x in p):
            assert norm.sign()<0 and [i for i,x in enumerate(p)if not x]==list(pair)
            feasible.append((rec,v))
        else:
            assert any(x.sign()>0 for x in p)and any(x.sign()<0 for x in p);rejected.append(rec)
    assert {tuple(r['pair'])for r,v in feasible}=={(0,1),(0,3),(1,2),(2,3)}
    assert len(rejected)==2
    center=[sum((v[i]for _,v in feasible),K())for i in range(3)];cp=[dot(r,center)for r in cov]
    cn=dot(center,[dot(r,center)for r in B]);assert cn.sign()<0 and all(x.sign()<0 for x in cp)
    future=[dot(v,[dot(r,center)for r in B])for _,v in feasible];assert all(x.sign()<0 for x in future)
    assert {tuple(r['pair'])for r in saved['corners']}=={tuple(r['pair'])for r,v in feasible}
    return dict(physical_inertia=phys,conjugate_inertia=conj,basis=basis,projected_gram=[[x.enc()for x in r]for r in S],
                finite_vertices=[r for r,v in feasible],rejected_opposites=rejected,
                strict_center=dict(vector=[x.enc()for x in center],pairings=[x.enc()for x in cp],norm=cn.enc(),future_vertex_products=[x.enc()for x in future]))


def run():
    data=json.loads((ROOT/f'results/{BASE}.json').read_text());geo=json.loads((ROOT/f'results/{BASE}_geometry.json').read_text())
    union=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    contacts=json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    lp=lp_checks(data,union,contacts);ids={r['class_index']for m in lp['sparse_minima']for r in m['faces']}
    assert len(ids)==15 and ids=={r['class_index']for r in geo['records']}
    geometries=[]
    for saved in geo['records']:
        row=union['labelled_union_supports'][saved['union_index']]
        geometries.append(dict(class_index=saved['class_index'],union_index=saved['union_index'],**geometry(row,saved)))
    assert sum(len(r['finite_vertices'])for r in geometries)==60
    paths=[f'scripts/{TAG}.py',f'scripts/{HELPER}.py',f'results/{HELPER}.json']
    for stem in(BASE,BASE+'_geometry','fresh_sqrt5_quad_exhaustion_20261015','fresh_sqrt5_quad_contacts_20261015'):
        paths += [f'scripts/{stem}.py',f'results/{stem}.json']
    manifests=Q.manifest_check(paths)
    return dict(status='PASS_INDEPENDENT_COMPONENT_LOCAL_WEIGHT_BARRIER_AND_FIFTEEN_ACTUAL_QUADRILATERALS',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},recursive_json_manifests=manifests,
                independent_weight_certificate=lp,independent_geometry=geometries,
                counts=dict(sealed_faces=112,variables=94,component_types=16,barrier_faces=15,pair_intersections=90,finite_vertices=60,rejected_opposites=30),
                scope='Fixed distance-two baseline plus nonnegative component-local corrections, zero on disconnected pairs. No obstruction to changed baselines/cross-component/global rules or global compact extensions. Exact local compact quadrilaterals only.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text());print('PASS saved sealed-weight review replay')
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'])
