#!/usr/bin/env python3
"""Joint E8 witnesses and the spectator-sensitive D7 weight obstruction.

Fresh diagnostics only: no dimension exclusion or compact ambient construction.
Uses the frozen independently reviewed Fraction Qsqrt5 arithmetic read-only.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_quad_filters_20261015_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_arithmetic_weight_obstructions_20261015'
Z=A.Z;ONE=A.ONE;TWO=A.TWO;PHI=(F(1,2),F(1,2))

def dot(u,v):
    ans=Z
    for a,b in zip(u,v):ans=A.add(ans,A.mul(a,b))
    return ans
def mv(M,v):return[dot(row,v)for row in M]
def sub(M,ids):return[[M[i][j]for j in ids]for i in ids]
def enc(v):return[A.enc(x)for x in v]
def encm(M):return[enc(v)for v in M]
def labels(row):
    n=row['nodes'];L=[[2]*n for _ in range(n)]
    for(i,j),m in zip(combinations(range(n),2),row['labels']):L[i][j]=L[j][i]=m
    return L
def components(L,ids):
    todo=set(ids);out=[]
    while todo:
        first=min(todo);todo.remove(first);seen=[first]
        for i in seen:
            for j in sorted(todo):
                if L[i][j]!=2:todo.remove(j);seen.append(j)
        out.append(sorted(seen))
    return out
def template(n):
    L=[[2]*n for _ in range(n)]
    edges=[(i,i+1)for i in range(n-2)]+[(2,7)]if n==8 else[(i,i+1)for i in range(5)]+[(4,6)]
    for i,j in edges:L[i][j]=L[j][i]=3
    return L
def maps(L,ids,T):
    deg={i:sum(L[i][j]!=2 for j in ids if j!=i)for i in ids}
    td=[sum(T[i][j]!=2 for j in range(len(T))if j!=i)for i in range(len(T))]
    order=sorted(ids,key=lambda i:(-deg[i],i));out=[]
    def rec(mp,used):
        if len(mp)==len(ids):out.append(dict(mp));return
        i=order[len(mp)]
        for j in range(len(T)):
            if j in used or deg[i]!=td[j]:continue
            if all(L[i][k]==T[j][v]for k,v in mp.items()):rec({**mp,i:j},used|{j})
    rec({},set())
    return out
def canonical_corner(row,n):
    L=labels(row);out=[]
    for a,b in product((0,2),(1,3)):
        ids=row['white']+[a,b]
        for comp in components(L,ids):
            if a not in comp or b not in comp or len(comp)!=n:continue
            iso=maps(L,comp,template(n))
            if iso:out.append((a,b,comp,sorted(set(ids)-set(comp)),iso))
    assert len(out)==1
    return L,out[0]

def nullvector(rows,n):
    M=[r[:]for r in rows];piv=[];k=0
    for j in range(n):
        i=next((i for i in range(k,len(M))if M[i][j]!=Z),None)
        if i is None:continue
        M[k],M[i]=M[i],M[k];a=M[k][j];M[k]=[A.div(v,a)for v in M[k]]
        for i in range(len(M)):
            if i!=k:
                a=M[i][j];M[i]=[A.sub(v,A.mul(a,w))for v,w in zip(M[i],M[k])]
        piv.append(j);k+=1
        if k==len(M):break
    free=[j for j in range(n)if j not in piv]
    if len(free)!=1:return None
    v=[Z]*n;v[free[0]]=ONE
    for i,j in enumerate(piv):v[j]=A.neg(M[i][free[0]])
    assert all(dot(r,v)==Z for r in rows)
    return v

def coherent_E8(union,ids):
    G=[[TWO if i==j else A.neg(ONE)if template(8)[i][j]==3 else Z for j in range(8)]for i in range(8)]
    assert A.inertia(G)==[8,0,0];J=A.inverse(G)
    b0=[PHI if i==0 else Z for i in range(8)];b6=[PHI if i==6 else Z for i in range(8)]
    d0=A.sub(dot(b0,mv(J,b0)),TWO);d6=A.sub(dot(b6,mv(J,b6)),TWO)
    assert d0==(F(4),F(2))and d6==(F(1),F(1))
    assert d0==A.mul(A.mul(PHI,PHI),d6)
    assert A.sub(dot(b0,mv(J,b6)),A.mul(PHI,d6))==Z
    M=[[TWO if i==j else Z for j in range(10)]for i in range(10)]
    for i,j in product(range(8),repeat=2):M[i][j]=G[i][j]
    M[0][8]=M[8][0]=A.neg(PHI);M[6][9]=M[9][6]=A.neg(PHI)
    assert A.inertia(M)==[8,1,1]
    assert A.inertia([[A.conj(z)for z in row]for row in M])==[9,0,1]
    mappings=[]
    for ci in ids:
        ui=union['marked_isomorphism_classes'][ci]['representative'];row=union['labelled_union_supports'][ui]
        L,(a,b,old,spectator,isomorphisms)=canonical_corner(row,8)
        assert not spectator and len(isomorphisms)==1
        mp=isomorphisms[0];new=[i for i in range(10)if i not in old]
        for z in new:
            contacts={mp[i]:L[z][i]for i in old if L[z][i]!=2}
            assert contacts in({0:5},{6:5})
            mp[z]=8 if contacts=={0:5} else 9
        assert sorted(mp.values())==list(range(10))
        raw=A.matrix(row,A.coloring(row))
        assert all(raw[i][j]==M[mp[i]][mp[j]]for i,j in product(range(10),repeat=2))
        mappings.append(dict(class_index=ci,union_index=ui,vertex_pair=sorted([mp[a],mp[b]]),
            full_bijection=[[i,mp[i]]for i in range(10)],literal_Gram_equalities=100))
    # Explicit rank9 realization in the basis oldE8 + extra wall9.
    basis=list(range(8))+[9];Q=sub(M,basis);IQ=A.inverse(Q)
    normals=[mv(IQ,[M[j][i]for j in basis])for i in range(10)]
    assert all(dot(normals[i],mv(Q,normals[j]))==M[i][j]for i,j in product(range(10),repeat=2))
    H=mv(IQ,[Z]*8+[A.neg(ONE)])
    olddual=mv(J,[ONE]*8)+[Z]
    eps=(F(1,10000),F())
    C=[A.sub(a,A.mul(eps,b))for a,b in zip(H,olddual)]
    cp=[dot(v,mv(Q,C))for v in normals];cn=dot(C,mv(Q,C))
    assert A.sign(cn)<0 and all(A.sign(v)<0 for v in cp)
    facets=[]
    for i,v in enumerate(normals):
        P=[A.sub(a,A.mul(A.div(cp[i],TWO),b))for a,b in zip(C,v)]
        signs=[dot(w,mv(Q,P))for w in normals]
        assert A.sign(dot(P,mv(Q,P)))<0 and signs[i]==Z and all(A.sign(z)<0 for j,z in enumerate(signs)if i!=j)
        facets.append(dict(index=i,point=enc(P),norm=A.enc(dot(P,mv(Q,P))),pairings=enc(signs)))
    rays={}
    cov=[mv(Q,v)for v in normals]
    for S in combinations(range(10),8):
        v=nullvector([cov[i]for i in S],9)
        if v is None:continue
        pairings=[dot(r,v)for r in cov]
        if all(A.sign(z)>=0 for z in pairings):v=list(map(A.neg,v));pairings=list(map(A.neg,pairings))
        if not all(A.sign(z)<=0 for z in pairings):continue
        lead=next(z for z in v if z!=Z);canonical=tuple(A.div(z,lead)for z in v)
        norm=dot(v,mv(Q,v));key=str(canonical)
        rays[key]=dict(support=[i for i,z in enumerate(pairings)if z==Z],point=v,norm=norm)
    counts={str(s):sum(A.sign(r['norm'])==s for r in rays.values())for s in(-1,0,1)}
    assert all(A.sign(dot(r['point'],mv(Q,C)))<0 for r in rays.values()if A.sign(r['norm'])<0)
    candidates=[r for r in rays.values()if A.sign(r['norm'])>0 and A.sign(dot(r['point'],mv(Q,C)))<0]
    assert candidates
    v=candidates[0]['point'];epsilon=F(1,1024)
    while True:
        E=[A.add(a,A.mul((epsilon,F()),b))for a,b in zip(v,C)]
        en=dot(E,mv(Q,E));ep=[dot(w,mv(Q,E))for w in normals]
        if A.sign(en)>0 and all(A.sign(z)<0 for z in ep):break
        epsilon/=2
    return dict(standard_E8_edges=[[i,j]for i,j in combinations(range(8),2)if G[i][j]!=Z],
        Gram=encm(M),physical_inertia=[8,1,1],conjugate_inertia=[9,0,1],new_contacts=[{'node':0,'label':5},{'node':6,'label':5}],
        height_squares=enc([d0,d6]),height_ratio= A.enc(PHI),mutual_pairing=A.enc(Z),marked_core_maps=mappings,
        basis=basis,basis_Gram=encm(Q),normals=[enc(v)for v in normals],strict_center=enc(C),center_norm=A.enc(cn),center_pairings=enc(cp),
        facet_witnesses=facets,complete_8_support_tests=45,extreme_ray_sign_counts=counts,
        strict_escape=dict(point=enc(E),norm=A.enc(en),pairings=enc(ep),ray_support=candidates[0]['support'],epsilon=str(epsilon)),
        scope='All nine E8 forcing marked cores are literally the same ten-wall arithmetic-compatible local chamber. All ten facets are actual, but the explicit strict spacelike escape proves this chamber has infinite volume. No compact global extension is asserted.')

def D7_spectators(union,contacts,ids):
    T=template(7);autos=maps(T,list(range(7)),T);assert len(autos)==2
    recs=[];pure=set();dependent=set()
    for ci in ids:
        ui=union['marked_isomorphism_classes'][ci]['representative'];row=union['labelled_union_supports'][ui]
        L,(a,b,old,spectator,isomorphisms)=canonical_corner(row,7);assert len(isomorphisms)==2
        mp=isomorphisms[0];marked=tuple(sorted((mp[a],mp[b])))
        orbit=sorted({tuple(sorted((p[marked[0]],p[marked[1]])))for p in autos})
        assert all(all(L[i][j]==2 for j in row['white']+[a,b]if j!=i)for i in spectator)
        cr=next(z for z in contacts['records']if z['class_index']==ci)
        assert cr['singular']and not cr['accepted']and cr['status']=='UNPAID_FINITE_CONTACT_DOMAIN'
        assert all(z['status']!='ACCEPT_NECESSARY_FINITE_CONTACT'for z in cr['trials'])
        if spectator:dependent.update(orbit)
        else:pure.update(orbit)
        recs.append(dict(class_index=ci,union_index=ui,connected_corner_pair=[a,b],D7_nodes=old,
            D7_bijection=[[i,mp[i]]for i in sorted(mp)],canonical_pair=list(marked),pair_orbit=[list(p)for p in orbit],
            isolated_white_spectators=spectator,spectator_new_black_contacts=[[i,j,L[i][j]]for i in spectator for j in range(4)if j not in(a,b)and L[i][j]!=2],
            nonempty_contact_profiles=cr['structural_profiles'],accepted_unpaid_profiles=0,
            may_become_H4_paid=any(z['status']=='BECOMES_PAID_BY_ENTIRE_H4'for z in cr['trials'])))
    assert pure=={(0,4),(1,4)}
    assert dependent=={(0,3),(0,5),(0,6),(1,5),(1,6),(2,5),(2,6)}
    assert pure.isdisjoint(dependent)
    return dict(records=recs,pure_forced_long_pairs=[list(p)for p in sorted(pure)],
        requires_entire_A1_when_unpaid=[list(p)for p in sorted(dependent)],
        sparse_barrier_D7_totals=dict(distance_two_base=12,pure_forcing_only=14,including_A1_dependent_forcing=21),
        scope='These are the six sparse dual witnesses, not a classification of every D7 quadrilateral. The entire-A1 conclusion inherits the separately reviewed conditional unpaid sealing theorem. If an attachment creates an H4 payment, the face no longer requires the D7 correction.')

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    union=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    contacts=json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    lp=json.loads((ROOT/'results/fresh_sqrt5_sealed_weights_20261015.json').read_text())
    ids={c['type']:sorted(z['class_index']for z in c['dual'])for c in lp['component_minima']}
    e=coherent_E8(union,ids['E8']);d=D7_spectators(union,contacts,ids['D7'])
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_quad_filters_20261015_review.py']
    for tag in('fresh_sqrt5_quad_exhaustion_20261015','fresh_sqrt5_quad_contacts_20261015','fresh_sqrt5_quad_contacts_20261015_review',
        'fresh_sqrt5_sealed_weights_20261015','fresh_sqrt5_sealed_weights_20261015_geometry','fresh_sqrt5_sealed_weights_20261015_review'):
        paths +=[f'scripts/{tag}.py',f'results/{tag}.json']
    for p in paths:
        if p.endswith('.json'):
            for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert sha(q)==h,(p,q)
    return dict(status='PASS_EXACT_COHERENT_E8_AND_SPECTATOR_SENSITIVE_D7_DIAGNOSTIC',sources={p:sha(p)for p in paths},
        coherent_E8=e,D7_spectator_context=d,
        scope='Fresh arithmetic/local-weight obstruction analysis. No new dimension theorem and no claim that the restricted D7 correction suffices for all faces. Frozen component-only barrier remains valid; its extension to a full-forest-sensitive rule requires preserving spectator context.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['coherent_E8']['extreme_ray_sign_counts'])
